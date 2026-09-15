#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import re
import time

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from cosmicqc import find_outliers, label_outliers
from cytodataframe import CytoDataFrame
from IPython.display import display


# In[2]:


# Set parameters for papermill to use for processing
plate_id = "BR00148920"

# Whether to render per-condition image diagnostics (CytoDataFrame previews of
# outlier crops, pulled from the bandicoot network mount). This is the slow
# part of this notebook and is only useful for interactive review -- the
# final label_outliers export below only depends on the *_thresholds dicts,
# not on these previews. Batch/papermill runs across many plates should set
# this to False.
render_diagnostics = True


# In[3]:


# Directory containing the converted profiles
data_dir = pathlib.Path("/media/18tbdrive2/SK-N-AS_repo1_profiles/converted_profiles")

# Directory to save qc results from cosmicqc
cleaned_dir = pathlib.Path("./data/qc_results")
cleaned_dir.mkdir(exist_ok=True)


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
    "{compartment}_Location_Center_X",
    "{compartment}_Location_Center_Y",
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
file_path = data_dir / f"{plate_id}_converted.parquet"

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
    "Nuclei_Intensity_MassDisplacement_DNA": 1.8,
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
        data_outline_context_dir=f"/home/jenna/mnt/bandicoot/PCCMA_data/SK-N-AS_repo1_screen/SQLite_outputs/{plate_id}",
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
            by="Nuclei_Intensity_MassDisplacement_DNA", ascending=True
        ).head(10).T
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
        data_outline_context_dir=f"/home/jenna/mnt/bandicoot/PCCMA_data/SK-N-AS_repo1_screen/SQLite_outputs/{plate_id}",
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
            by="Nuclei_AreaShape_Solidity", ascending=True
        ).head(10).T
    )
    # solidity_nuclei_outliers_cdf.sample(n=5)


# In[11]:


# Find background segmented as nuclei outliers for the current plate
background_segmentation_thresholds = {
    "Nuclei_Intensity_IntegratedIntensity_DNA": -1.5,
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
        data_outline_context_dir=f"/home/jenna/mnt/bandicoot/PCCMA_data/SK-N-AS_repo1_screen/SQLite_outputs/{plate_id}",
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
        ).head(10).T
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
    "Cells_Intensity_IntegratedIntensity_DNA": 2.5,
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
        data_outline_context_dir=f"/home/jenna/mnt/bandicoot/PCCMA_data/SK-N-AS_repo1_screen/SQLite_outputs/{plate_id}",
        segmentation_file_regex=outline_to_orig_mapping_cells,
        display_options={
            "center_dot": True,
            "brightness": 25,
            "offset_bounding_box": {
                "x_min": -80,
                "y_min": -80,
                "x_max": 80,
                "y_max": 80,
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
        ).head(5).T
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

# Print the number of outliers and percentage of outliers for each QC condition
qc_annotations.filter(like="Metadata_cqc_").sum().to_frame(name="num_outliers").assign(
    percentage_outliers=lambda x: (x["num_outliers"] / len(qc_annotations)) * 100
)

