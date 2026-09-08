#!/usr/bin/env python
# coding: utf-8

# # Perform single-cell quality control
# 
# > Note: There will be commented out code for displaying the CytoDataFrames as it is helpful to switch between views during optimization.

# In[1]:


import pathlib
import re
import time

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from IPython.display import display

from cytodataframe import CytoDataFrame
from cosmicqc import find_outliers, label_outliers


# In[2]:


# Set parameters for papermill to use for processing
plate_id = "BR00149544"

# Whether to render per-condition image diagnostics (CytoDataFrame previews of
# outlier crops, pulled from the bandicoot network mount). This is the slow
# part of this notebook and is only useful for interactive review -- the
# final label_outliers export below only depends on the *_thresholds dicts,
# not on these previews. Batch/papermill runs across many plates should set
# this to False.
render_diagnostics = True


# In[3]:


# Directory containing the converted profiles after merging
data_dir = pathlib.Path("./data/merged_profiles/")

# Directory to save qc results from cosmicqc
cleaned_dir = pathlib.Path("./data/qc_results")
cleaned_dir.mkdir(exist_ok=True)

# Combined QC summary CSV (percentage of cells flagged per condition, and
# overall, one row per plate/condition) so outlier plates are easy to spot
# across the whole screen without opening each qc_annotations parquet.
# run_pipeline.sh runs this notebook for one plate at a time, in sequence,
# so appending to this shared file plate-by-plate is safe.
qc_summary_path = cleaned_dir / "qc_summary.csv"


# In[4]:


compartments = ["Nuclei", "Cells"]

shared_metadata_columns = [
    "Image_Metadata_Plate",
    "Image_Metadata_Well",
    "Image_Metadata_Site",
    "Image_Metadata_Row",
    "Image_Metadata_Col",
    "Image_FileName_OrigDNA",
    "Image_FileName_OrigAGP",
    "Image_PathName_OrigDNA",
    "Image_PathName_OrigAGP",
    # per-cell identifiers needed to join qc annotations back onto the
    # full single-cell profile after pycytominer's annotate step
    "Metadata_ImageNumber",
    "Metadata_Cells_Number_Object_Number",
    "Metadata_Nuclei_Number_Object_Number",
]

compartment_location_templates = [
    "Metadata_{compartment}_Location_Center_X",
    "Metadata_{compartment}_Location_Center_Y",
    "{compartment}_AreaShape_BoundingBoxMaximum_X",
    "{compartment}_AreaShape_BoundingBoxMaximum_Y",
    "{compartment}_AreaShape_BoundingBoxMinimum_X",
    "{compartment}_AreaShape_BoundingBoxMinimum_Y",
]

compartment_qc_feature_templates = [
    "{compartment}_Intensity_IntegratedIntensity_DNA",
    "{compartment}_AreaShape_Solidity",
    "{compartment}_Intensity_MassDisplacement_DNA",
]

metadata_columns = shared_metadata_columns + [
    template.format(compartment=compartment)
    for compartment in compartments
    for template in compartment_location_templates
]

needed_columns = metadata_columns + [
    template.format(compartment=compartment)
    for compartment in compartments
    for template in compartment_qc_feature_templates
]


# In[5]:


# Construct the file path for the given plate_id
file_path = data_dir / f"{plate_id}.parquet"

if file_path.exists():
    start_time = time.time()  # Start timer for loading

    # Load the DataFrame with pandas
    plate_df = pd.read_parquet(file_path, engine="pyarrow", columns=needed_columns)

    end_time = time.time()  # End timer for loading
    print(
        f"Loaded plate: {plate_id}, Shape: {plate_df.shape}, Time taken: {end_time - start_time:.2f} seconds"
    )
else:
    print(f"Parquet file for plate {plate_id} not found.")


# In[6]:


# Confirm PathName is a column in input DataFrame
assert "Image_PathName_OrigDNA" in plate_df.columns, "Column 'Image_PathName_OrigDNA' is missing from the input DataFrame"


# In[7]:


correct_parent = "/home/jenna/mnt/bandicoot/PCCMA_data"

for col in plate_df.columns:
    if "PathName" in col and "Illum" not in col:
        plate_df[col] = plate_df[col].apply(
            lambda x: (
                re.sub(r"^.*ALSF_screen_data/", correct_parent + "/", x)
                if isinstance(x, str)
                else x
            )
        )

# Print example image path after fix
print(plate_df["Image_PathName_OrigDNA"].dropna().iloc[0])


# ## Run QC with `nuclei` features

# In[8]:


# Set the compartment of choice to perform QC at the start (will change later)
compartment = "Nuclei"

# create an outline and orig mapping dictionary to map original images to outlines
# note: we turn off formatting here to avoid the key-value pairing definition
# from being reformatted by black, which is normally preferred.
# fmt: off
outline_to_orig_mapping = {
    rf"{compartment}Outlines_{record['Image_Metadata_Plate']}_{record['Image_Metadata_Well']}_{record['Image_Metadata_Site']}.tiff": 
    rf"r{int(record['Image_Metadata_Row']):02d}c{int(record['Image_Metadata_Col']):02d}f{int(record['Image_Metadata_Site']):02d}p(\d{{2}})-ch\d+sk\d+fk\d+fl\d+\.tiff"
    for record in plate_df[
        [
            "Image_Metadata_Plate",
            "Image_Metadata_Well",
            "Image_Metadata_Site",
            "Image_Metadata_Row",
            "Image_Metadata_Col",
        ]
    ].to_dict(orient="records")
}
# fmt: on

next(iter(outline_to_orig_mapping.items()))


# In[9]:


# Find large nuclei outliers for the current plate
nuclei_clustered_thresholds = {
    "Nuclei_Intensity_MassDisplacement_DNA": 0.7,
    "Nuclei_Intensity_IntegratedIntensity_DNA": 1.5,
}

nuclei_clustered_outliers = find_outliers(
    df=plate_df,
    metadata_columns=metadata_columns,
    feature_thresholds=nuclei_clustered_thresholds,
)

if render_diagnostics:
    # MUST SET DATA AS DATAFRAME FOR OUTLINE DIR TO WORK
    nuclei_clustered_outliers_cdf = CytoDataFrame(
        data=pd.DataFrame(nuclei_clustered_outliers),
        data_outline_context_dir=f"/home/jenna/mnt/bandicoot/PCCMA_data/CHP-134_repo1_screen_outputs/SQLite_outputs/{plate_id}",
        segmentation_file_regex=outline_to_orig_mapping,
        display_options={
            "brightness": 1,
            "offset_bounding_box": {
                "x_min": -40,
                "y_min": -40,
                "x_max": 40,
                "y_max": 40,
            },
        },
    )[
        [
            "Image_FileName_OrigDNA",
            "Nuclei_Intensity_MassDisplacement_DNA",
            "Nuclei_Intensity_IntegratedIntensity_DNA",
        ]
    ]

    print(nuclei_clustered_outliers_cdf.shape)
    display(
        nuclei_clustered_outliers_cdf.sort_values(
            by="Nuclei_Intensity_IntegratedIntensity_DNA", ascending=False
        ).head(5)
    )
    # nuclei_clustered_outliers_cdf.sample(n=2, random_state=0)


# In[10]:


# Find low nuclei solidity outliers for the current plate
missegmented_nuclei_thresholds = {
    "Nuclei_AreaShape_Solidity": -2.5,
}

solidity_nuclei_outliers = find_outliers(
    df=plate_df,
    metadata_columns=metadata_columns,
    feature_thresholds=missegmented_nuclei_thresholds,
)

if render_diagnostics:
    # MUST SET DATA AS DATAFRAME FOR OUTLINE DIR TO WORK
    solidity_nuclei_outliers_cdf = CytoDataFrame(
        data=pd.DataFrame(solidity_nuclei_outliers),
        data_outline_context_dir=f"/home/jenna/mnt/bandicoot/PCCMA_data/CHP-134_repo1_screen_outputs/SQLite_outputs/{plate_id}",
        segmentation_file_regex=outline_to_orig_mapping,
        display_options={
            "brightness": 1,
            "offset_bounding_box": {
                "x_min": -40,
                "y_min": -40,
                "x_max": 40,
                "y_max": 40,
            },
        },
    )[
        [
            "Image_FileName_OrigDNA",
            "Nuclei_AreaShape_Solidity",
        ]
    ]

    print(solidity_nuclei_outliers_cdf.shape)
    display(
        solidity_nuclei_outliers_cdf.sort_values(
            by="Nuclei_AreaShape_Solidity", ascending=False
        ).head(5)
    )
    # solidity_nuclei_outliers_cdf.sample(n=5)


# In[11]:


# Find background segmented as nuclei outliers for the current plate
background_segmentation_thresholds = {
    "Nuclei_Intensity_IntegratedIntensity_DNA": -1.8,
}

background_segmentation_outliers = find_outliers(
    df=plate_df,
    metadata_columns=metadata_columns,
    feature_thresholds=background_segmentation_thresholds,
)

if render_diagnostics:
    # MUST SET DATA AS DATAFRAME FOR OUTLINE DIR TO WORK
    background_segmentation_outliers_cdf = CytoDataFrame(
        data=pd.DataFrame(background_segmentation_outliers),
        data_outline_context_dir=f"/home/jenna/mnt/bandicoot/PCCMA_data/CHP-134_repo1_screen_outputs/SQLite_outputs/{plate_id}",
        segmentation_file_regex=outline_to_orig_mapping,
        display_options={
            "brightness": 1,
            "offset_bounding_box": {
                "x_min": -40,
                "y_min": -40,
                "x_max": 40,
                "y_max": 40,
            },
        },
    )[
        [
            "Image_FileName_OrigDNA",
            "Nuclei_Intensity_IntegratedIntensity_DNA",
        ]
    ]

    print(background_segmentation_outliers_cdf.shape)
    display(
        background_segmentation_outliers_cdf.sort_values(
            by="Nuclei_Intensity_IntegratedIntensity_DNA", ascending=False
        ).head(5)
    )
    # background_segmentation_outliers_cdf.sample(n=2, random_state=0)


# ## Run QC with `cell` features

# In[12]:


# Set the compartment of choice to perform QC at the start (will change later)
compartment = "Cells"

# create an outline and orig mapping dictionary to map original images to outlines
# note: we turn off formatting here to avoid the key-value pairing definition
# from being reformatted by black, which is normally preferred.
# fmt: off
outline_to_orig_mapping_cells = {
    rf"{compartment}Outlines_{record['Image_Metadata_Plate']}_{record['Image_Metadata_Well']}_{record['Image_Metadata_Site']}.tiff": 
    rf"r{int(record['Image_Metadata_Row']):02d}c{int(record['Image_Metadata_Col']):02d}f{int(record['Image_Metadata_Site']):02d}p(\d{{2}})-ch\d+sk\d+fk\d+fl\d+\.tiff"
    for record in plate_df[
        [
            "Image_Metadata_Plate",
            "Image_Metadata_Well",
            "Image_Metadata_Site",
            "Image_Metadata_Row",
            "Image_Metadata_Col",
        ]
    ].to_dict(orient="records")
}
# fmt: on

next(iter(outline_to_orig_mapping_cells.items()))


# In[13]:


# Find cell outliers for the current plate
cell_outliers_thresholds = {
    # Set low to attempt to detect all instances of abnormally high int in nuclei for whole cells
    "Cells_Intensity_IntegratedIntensity_DNA": 2,
}

cell_outliers = find_outliers(
    df=plate_df,
    metadata_columns=metadata_columns,
    feature_thresholds=cell_outliers_thresholds,
)

if render_diagnostics:
    # MUST SET DATA AS DATAFRAME FOR OUTLINE DIR TO WORK
    cell_outliers_cdf = CytoDataFrame(
        data=pd.DataFrame(cell_outliers),
        data_outline_context_dir=f"/home/jenna/mnt/bandicoot/PCCMA_data/CHP-134_repo1_screen_outputs/SQLite_outputs/{plate_id}",
        segmentation_file_regex=outline_to_orig_mapping_cells,
        display_options={
            "center_dot": True,
            "brightness": 25,
            "offset_bounding_box": {
                "x_min": -120,
                "y_min": -120,
                "x_max": 120,
                "y_max": 120,
            },
        },
    )[
        [
            "Image_FileName_OrigDNA",
            "Image_FileName_OrigAGP",
            "Cells_Intensity_IntegratedIntensity_DNA",
        ]
    ]

    print(cell_outliers_cdf.shape)
    display(
        cell_outliers_cdf.sort_values(
            by="Cells_Intensity_IntegratedIntensity_DNA", ascending=False
        ).head(5)
    )
    # cell_outliers_cdf.sample(n=5, random_state=0)


# ## Save annotations file for plate with QC flags

# In[14]:


# Use label outliers to export qc.parquet file per plate

# Combine every QC condition defined above into one set of named thresholds
qc_feature_thresholds = {
    "nuclei_clustered": nuclei_clustered_thresholds,
    "missegmented_nuclei": missegmented_nuclei_thresholds,
    "background_segmentation": background_segmentation_thresholds,
    "cell_outliers": cell_outliers_thresholds,
}

# export_as_annotations=True keeps only metadata_columns + the generated
# Metadata_cqc_*_is_outlier flags
qc_annotations = label_outliers(
    df=plate_df,
    feature_thresholds=qc_feature_thresholds,
    export_path=str(cleaned_dir / f"{plate_id}_qc_annotations.parquet"),
    export_as_annotations=True,
    annotation_metadata_columns=metadata_columns,
)

outlier_flags = qc_annotations.filter(like="Metadata_cqc_")

# Percentage of cells flagged as an outlier, per QC condition
qc_condition_summary = outlier_flags.sum().to_frame(name="num_outliers").assign(
    percentage_outliers=lambda x: (x["num_outliers"] / len(qc_annotations)) * 100
)
qc_condition_summary.index.name = "qc_condition"

# Overall: a cell fails the plate if it is flagged by at least one condition
any_outlier = outlier_flags.any(axis=1)
qc_condition_summary.loc["any_condition (overall)"] = [
    any_outlier.sum(),
    (any_outlier.sum() / len(qc_annotations)) * 100,
]

qc_condition_summary = qc_condition_summary.reset_index().assign(plate_id=plate_id)[
    ["plate_id", "qc_condition", "num_outliers", "percentage_outliers"]
]

# Append this plate's rows to the combined screen-wide summary, replacing
# any existing rows for this plate so a rerun does not duplicate them.
# Safe because run_pipeline.sh processes plates sequentially, one at a time.
if qc_summary_path.exists():
    existing_summary = pd.read_csv(qc_summary_path)
    existing_summary = existing_summary[existing_summary["plate_id"] != plate_id]
    qc_condition_summary = pd.concat(
        [existing_summary, qc_condition_summary], ignore_index=True
    )

qc_condition_summary.to_csv(qc_summary_path, index=False)

qc_condition_summary[qc_condition_summary["plate_id"] == plate_id]

