#!/usr/bin/env python
# coding: utf-8

# # Process bulk profiles
# 
# > NOTE: For bulk profiles, we normalize to the negative controls and use MAD robustize as the "standard".

# ## Import libraries

# In[1]:


import csv
import gc
import os
import pathlib
import pprint
import time
from datetime import datetime, timezone

import pandas as pd
from pycytominer import aggregate, annotate, feature_select, normalize
from pycytominer.cyto_utils import output


# ## Set paths and variables

# In[ ]:


# Directory containing one merged profile parquet per plate
merged_dir = pathlib.Path("./data/merged_profiles")

# Directory containing per-plate QC annotation files from 2.single_cell_qc.ipynb
qc_dir = pathlib.Path("./data/qc_results")

# output path for bulk profiles
output_dir = pathlib.Path("./data/bulk_profiles")
output_dir.mkdir(parents=True, exist_ok=True)

# path for platemap directory
platemap_dir = pathlib.Path("../0.download_data/metadata")

# load in barcode platemap
barcode_platemap = pd.read_csv(platemap_dir / "barcode_platemap.csv")

# plate_id always uses underscores (e.g. "Assay_Plate_1_3"), but a few
# barcodes in this file use a space instead (e.g. "Assay Plate_1_3")
barcode_platemap["Plate Barcode"] = barcode_platemap["Plate Barcode"].str.replace(
    " ", "_", regex=False
)

# extract the plate names from the merged profile file names
plate_names = sorted(file.stem for file in merged_dir.glob("*.parquet"))

# Wells with no compound (`Metadata_Batch_Id` is blank) that still received
# the DMSO vehicle are the negative controls
neg_control_query = 'Metadata_Batch_Id.isna() and Metadata_Solvent == "DMSO"'

# operations to perform for feature selection
# NOTE: drop_na_columns runs before correlation_threshold on purpose. Several
# Costes correlation features are NaN for some cells, and pycytominer's
# correlation_threshold falls back to a much slower NaN-aware pandas .corr()
# (instead of a fast BLAS np.corrcoef) if any NaNs are still present in the
# feature columns it's given. Dropping NaN columns first keeps it on the fast path.
feature_select_ops = [
    "drop_na_columns",
    "blocklist", # default block list uses same standard CP naming convention
    "frequency_threshold",
    "variance_threshold",
    "correlation_threshold",
]

plate_names


# ## Set dictionary with plates to process

# In[3]:


# Create plate info dictionary
plate_info_dictionary = {
    plate_id: {
        "profile_path": str(merged_dir / f"{plate_id}.parquet"),

        # QC annotations are produced per-plate by 2.single_cell_qc.ipynb
        "qc_path": (
            str((qc_dir / f"{plate_id}_qc_annotations.parquet").resolve())
            if (qc_dir / f"{plate_id}_qc_annotations.parquet").exists()
            else None
        ),

        # Find the platemap file based on barcode match
        "platemap_path": (
            str(
                platemap_dir
                / barcode_platemap.loc[
                    barcode_platemap["Plate Barcode"] == plate_id, "File Name"
                ].values[0]
            )
            if plate_id in barcode_platemap["Plate Barcode"].values
            else None
        ),
    }
    for plate_id in plate_names
}

# Display the dictionary to verify the entries
pprint.pprint(plate_info_dictionary, indent=4)


# In[4]:


# Run configuration, both driven by environment variables so run_pipeline.sh
# can control them without papermill:
#   PLATE_ID  - restrict to a single plate (used to run each plate as its own
#               process so memory doesn't accumulate across plates). Leave
#               unset to process all plates, e.g. when running interactively.
#   OVERWRITE - if set (1/true/yes), reprocess and overwrite a plate's output
#               even if it already exists. Leave unset (the default) to skip
#               plates that already have output.
plate_id_filter = os.environ.get("PLATE_ID")
if plate_id_filter:
    if plate_id_filter not in plate_info_dictionary:
        raise ValueError(f"Unknown plate_id in PLATE_ID env var: {plate_id_filter}")
    plate_info_dictionary = {plate_id_filter: plate_info_dictionary[plate_id_filter]}

overwrite = os.environ.get("OVERWRITE", "").strip().lower() in ("1", "true", "yes")

plate_info_dictionary


# ## Process data with pycytominer

# In[ ]:


timing_log_path = output_dir / "timing_log.csv"

for plate_id, info in plate_info_dictionary.items():
    if info["qc_path"] is None:
        print(
            f"Skipping {plate_id}: no QC annotations yet "
            "(run 2.single_cell_qc.ipynb for this plate first)"
        )
        continue
    if info["platemap_path"] is None:
        print(f"Skipping {plate_id}: no platemap found in barcode_platemap.csv")
        continue

    # Output file paths for each file
    output_annotated_file = str(output_dir / f"{plate_id}_bulk_annotated.parquet")
    output_normalized_file = str(output_dir / f"{plate_id}_bulk_normalized.parquet")
    output_feature_select_file = str(
        output_dir / f"{plate_id}_bulk_feature_selected.parquet"
    )

    # Already fully processed, so this plate can safely be skipped on a rerun
    # (unless OVERWRITE is set, in which case it's reprocessed regardless)
    if pathlib.Path(output_feature_select_file).exists() and not overwrite:
        print(f"Skipping {plate_id}: already processed (found {output_feature_select_file})")
        continue

    print(f"Now performing pycytominer pipeline for {plate_id}")
    plate_start_time = time.time()

    # Load single-cell profile, its QC annotations, and the platemap
    single_cell_df = pd.read_parquet(info["profile_path"])
    qc_df = pd.read_parquet(info["qc_path"])
    platemap_df = pd.read_csv(info["platemap_path"]).rename(
        columns={"Well_Position": "Well"}
    )

    # Define the columns from the external QC file that indicate poor-quality segmentations
    cqc_cols = [col for col in qc_df.columns if col.startswith("Metadata_cqc_")]
    join_keys = [
        "Image_Metadata_Well", "Image_Metadata_Site",
        "Metadata_Nuclei_Location_Center_X", "Metadata_Nuclei_Location_Center_Y",
    ]
    assert not qc_df.duplicated(subset=join_keys).any(), f"{plate_id}: QC file has duplicate rows for the same cell"
    assert not single_cell_df.duplicated(subset=join_keys).any(), f"{plate_id}: profile file has duplicate rows for the same cell"

    # Step 1: Annotation (platemap + external QC metadata in one call)
    annotate_start_time = time.time()
    annotated_df = annotate(
        profiles=single_cell_df,
        platemap=platemap_df,
        join_on=["Metadata_Well", "Image_Metadata_Well"],
        external_metadata=qc_df[join_keys + cqc_cols],
        external_join_on=join_keys,
    )
    assert annotated_df[cqc_cols].isna().sum().sum() == 0, f"{plate_id}: some cells have no matching QC annotation"

    is_poor_quality = annotated_df[cqc_cols].any(axis=1)
    print(f"  Dropping {is_poor_quality.sum()} / {len(annotated_df)} poor-quality segmentations")
    annotated_df = annotated_df.loc[~is_poor_quality].drop(columns=cqc_cols).reset_index(drop=True)
    annotate_seconds = time.time() - annotate_start_time

    # Step 2: Aggregation
    # aggregate() only keeps `strata` columns plus the
    # aggregated numeric features, dropping everything else. Rather than
    # re-annotating the aggregated (well-level) data afterward to reattach
    # platemap metadata, include every platemap-derived metadata column that's
    # constant within a well directly in `strata`, so it survives aggregation
    aggregate_start_time = time.time()
    per_cell_cols = set(single_cell_df.columns)
    candidate_metadata_cols = [
        col for col in annotated_df.columns
        if col.startswith("Metadata_")
        and col not in ("Metadata_Plate", "Metadata_Well")
        and col not in per_cell_cols
    ]
    strata_cols = ["Metadata_Plate", "Metadata_Well"] + [
        col for col in candidate_metadata_cols
        if annotated_df.groupby("Metadata_Well")[col].nunique(dropna=False).le(1).all()
    ]
    aggregated_df = aggregate(
        population_df=annotated_df,
        operation="median",
        strata=strata_cols,
    )
    output(
        df=aggregated_df,
        output_filename=output_annotated_file,
        output_type="parquet",
    )
    aggregate_seconds = time.time() - aggregate_start_time

    # Clear memory
    del single_cell_df, qc_df, platemap_df, annotated_df
    gc.collect()

    # Step 3: Normalization (whole-plate) -- feed aggregate()'s own in-memory
    # result directly instead of re-reading output_aggregated_file back off disk.
    normalize_start_time = time.time()
    normalize(
        profiles=aggregated_df,
        method="mad_robustize", # use robustize to avoid influence of outliers in the normalization
        samples=neg_control_query, # normalize to negative controls only, not all wells (which would include compound wells)
        output_file=output_normalized_file,
        output_type="parquet",
    )
    normalize_seconds = time.time() - normalize_start_time

    # Clear memory
    del aggregated_df
    gc.collect()

    # Step 4: Feature selection
    feature_select_start_time = time.time()
    feature_select(
        output_normalized_file,
        operation=feature_select_ops,
        corr_threshold=0.90, # keep the same default value for correlation_threshold to be more strict as to identify the most informative features
        freq_cut=0.05, # keep the same default value for freq_cut
        unique_cut=0.01, # keep the same default value for unique_cut
        na_cutoff=0, # update na_cutoff from default to 0 to remove any columns with any NaN values (best for downstream modeling)
        output_file=output_feature_select_file,
        output_type="parquet",
    )
    feature_select_seconds = time.time() - feature_select_start_time

    total_seconds = time.time() - plate_start_time
    print(
        f"Bulk processing completed for {plate_id}! "
        f"(annotate={annotate_seconds:.1f}s, aggregate={aggregate_seconds:.1f}s, "
        f"normalize={normalize_seconds:.1f}s, feature_select={feature_select_seconds:.1f}s, "
        f"total={total_seconds:.1f}s)"
    )

    # Record timing so per-plate performance can be compared across runs
    write_header = not timing_log_path.exists()
    with open(timing_log_path, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(
                [
                    "plate_id",
                    "annotate_seconds",
                    "aggregate_seconds",
                    "normalize_seconds",
                    "feature_select_seconds",
                    "total_seconds",
                    "timestamp",
                ]
            )
        writer.writerow(
            [
                plate_id,
                f"{annotate_seconds:.1f}",
                f"{aggregate_seconds:.1f}",
                f"{normalize_seconds:.1f}",
                f"{feature_select_seconds:.1f}",
                f"{total_seconds:.1f}",
                datetime.now(timezone.utc).isoformat(timespec="seconds"),
            ]
        )


# ## Spherizing step
# 
# Adapted from Erik Serrano's work in the [fibrosis drug screen repository](https://github.com/WayScience/targeted_fibrosis_drug_screen/)

# In[ ]:


# Unlike the source pipeline this step was adapted from, CHP-134 has no
# batches / replicate plate groups to pool -- it's just one screen. So instead
# of looping per plate map and pooling that plate map's replicate plates,
# feature selection and the sphering fit/transform are all done once across
# every plate's pooled normalized profile, producing a single spherized
# profile for the whole screen (not one per plate).

spherized_output_dir = pathlib.Path("./data/spherized_profiles")
spherized_output_dir.mkdir(parents=True, exist_ok=True)

pooled_feature_select_file = output_dir / "SK-N-AS_repo1_screen_pooled_bulk_feature_selected.parquet"
output_spherized_file = spherized_output_dir / "SK-N-AS_repo1_screen_pooled_bulk_spherized.parquet"


# In[ ]:


if plate_id_filter:
    # Sphering needs every plate's normalized profile pooled together, so it
    # can't run off the single-plate slice of plate_info_dictionary that
    # PLATE_ID leaves in place -- skip until run with PLATE_ID unset.
    print(
        f"Skipping sphering: PLATE_ID={plate_id_filter} is set, so only one "
        "plate's data is loaded. Rerun with PLATE_ID unset once every plate "
        "has been bulk-processed."
    )
elif output_spherized_file.exists() and not overwrite:
    print(f"Skipping sphering: already spherized (found {output_spherized_file})")
else:
    normalized_paths = sorted(output_dir.glob("*_bulk_normalized.parquet"))

    print(f"Pooling {len(normalized_paths)} plates for sphering...")

    # step 1: concat every plate's normalized profile before feature selection
    concat_df = pd.concat(
        [pd.read_parquet(path) for path in normalized_paths],
        ignore_index=True,
    ).reset_index(drop=True)

    # Update any "Assay" plate strings to include underscores to match format
    # of other plates strings in the same file
    concat_df["Metadata_Plate"] = concat_df["Metadata_Plate"].str.replace(
        " ", "_", regex=False
    )

    # step 2a: Apply feature selection across the pooled screen to get a
    # common set of features for sphering.
    print("Feature selecting pooled screen...")
    feature_select_df = feature_select(
        profiles=concat_df,
        operation=feature_select_ops,
        na_cutoff=0, # updated from default to 0 to remove any columns with any NaN values (best for downstream modeling)
        corr_threshold=0.95, # increased from default to 0.95 to be more strict as to identify the most informative features
        freq_cut=0.05, # same as default
        output_file=pooled_feature_select_file,
        output_type="parquet",
    )

    # step 2b: Remove features with too little variation inside the exact
    # control population used to fit spherization.
    print(
        "Feature selecting pooled screen with frequency threshold within "
        "negative controls only..."
    )
    zero_negcon_var_fs_df = feature_select(
        profiles=feature_select_df,
        operation="frequency_threshold",
        freq_cut=0.05, # same as default
        unique_cut=0.01, # same as default
        samples=neg_control_query,
    )

    # step 3: Spherize/whiten the whole pooled screen using the pooled
    # negative controls as the reference population
    print("Sphering pooled screen using pooled negative controls...")
    normalize(
        profiles=zero_negcon_var_fs_df,
        method="spherize",
        samples=neg_control_query,
        spherize_center=True,
        spherize_method="ZCA-cor", # same as default
        spherize_epsilon=1e-6, #same as default
        output_file=output_spherized_file,
        output_type="parquet",
    )

    print(f"Saved feature-selected profiles to {pooled_feature_select_file}")
    print(f"Saved whole-screen spherized profile to {output_spherized_file}")

