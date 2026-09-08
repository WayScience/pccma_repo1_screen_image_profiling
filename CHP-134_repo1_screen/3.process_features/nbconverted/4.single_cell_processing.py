#!/usr/bin/env python
# coding: utf-8

# # Process single cell profiles
# 
# > NOTE: We normalize single-cells to the whole plate unlike bulk profiles, given that the sample size is much larger and less likely as impacted by variation.

# ## Import libraries

# In[ ]:


import csv
import gc
import os
import pathlib
import pprint
import time
from datetime import datetime, timezone

import pandas as pd
from pycytominer import annotate, feature_select, normalize
from pycytominer.cyto_utils import output


# ## Set paths and variables

# In[ ]:


# Directory containing one merged profile parquet per plate
merged_dir = pathlib.Path("./data/merged_profiles")

# Directory containing per-plate QC annotation files from 2.single_cell_qc.ipynb
qc_dir = pathlib.Path("./data/qc_results")

# output path for single-cell profiles
output_dir = pathlib.Path("./data/single_cell_profiles")
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

# In[ ]:


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


# In[ ]:


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


# Set up map for renaming metadata column(s)
column_name_mapping = {
    "Image_Metadata_Site": "Metadata_Site",
}

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

    # Set output paths
    output_annotated_file = str(output_dir / f"{plate_id}_sc_annotated.parquet")
    output_normalized_file = str(output_dir / f"{plate_id}_sc_normalized.parquet")
    output_feature_select_file = str(
        output_dir / f"{plate_id}_sc_feature_selected.parquet"
    )

    # Already fully processed, so this plate can safely be skipped on a rerun
    # (unless OVERWRITE is set, in which case it's reprocessed regardless)
    if pathlib.Path(output_feature_select_file).exists() and not overwrite:
        print(f"Skipping {plate_id}: already processed (found {output_feature_select_file})")
        continue

    print(f"Performing pycytominer pipeline for {plate_id}")
    plate_start_time = time.time()

    # Load in profile, its QC annotations, and the platemap
    profile_df = pd.read_parquet(info["profile_path"])
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

    # annotate()'s external_metadata merge (used below) doesn't validate that
    # join keys are unique the way merge(..., validate="one_to_one") does, so
    # check that explicitly here to keep the same data-integrity guarantee.
    assert not qc_df.duplicated(subset=join_keys).any(), f"{plate_id}: QC file has duplicate rows for the same cell"
    assert not profile_df.duplicated(subset=join_keys).any(), f"{plate_id}: profile file has duplicate rows for the same cell"

    print("Performing annotation for", plate_id, "...")
    # Step 1: Annotation -- merge in both the platemap and the per-cell QC flags
    annotate_start_time = time.time()
    annotated_df = annotate(
        profiles=profile_df,
        platemap=platemap_df,
        join_on=["Metadata_Well", "Image_Metadata_Well"],
        external_metadata=qc_df[join_keys + cqc_cols],
        external_join_on=join_keys,
    )

    # Rename Metadata column(s) using the rename() function
    annotated_df.rename(columns=column_name_mapping, inplace=True)

    # Drop any cell flagged by at least one QC condition (clustered/missegmented
    # nuclei, background segmented as a nucleus, whole-cell intensity outliers,
    # etc.)
    assert (
        annotated_df[cqc_cols].isna().sum().sum() == 0
    ), f"{plate_id}: some cells have no matching QC annotation"
    is_poor_quality = annotated_df[cqc_cols].any(axis=1)
    print(
        f"  Dropping {is_poor_quality.sum()} / {len(annotated_df)} "
        "poor-quality segmentations"
    )
    annotated_df = (
        annotated_df.loc[~is_poor_quality]
        .drop(columns=cqc_cols)
        .reset_index(drop=True)
    )

    # Assert "Metadata_Site" is now present after the rename() call above
    assert "Metadata_Site" in annotated_df.columns, f"{plate_id}: Metadata_Site column missing after rename()"

    # Save the modified annotated DataFrame after dropping poor-quality segmentations, so it can be inspected if needed
    output(
        df=annotated_df,
        output_filename=output_annotated_file,
        output_type="parquet",
    )
    annotate_seconds = time.time() - annotate_start_time

    # Clear memory
    del profile_df, qc_df, platemap_df
    gc.collect()

    print("Performing normalization for", plate_id, "...")
    # Step 2: Normalization
    normalize_start_time = time.time()
    normalize(
        profiles=annotated_df,
        method="standardize", # use standardize method as default
        output_file=output_normalized_file,
        output_type="parquet",
        samples="all" # apply normalization based on all samples
    )
    normalize_seconds = time.time() - normalize_start_time

    # Clear memory
    del annotated_df
    gc.collect()

    print("Performing feature selection for", plate_id, "...")
    # Step 3: Feature selection
    feature_select_start_time = time.time()
    feature_select(
        profiles=output_normalized_file,
        operation=feature_select_ops,
        corr_threshold=0.90, # keep the same default value for correlation_threshold to be more strict as to identify the most informative features
        freq_cut=0.05, # keep the same default value for freq_cut
        unique_cut=0.01, # keep the same default value for unique_cut
        na_cutoff=0, # update na_cutoff from default to 0 to remove any columns with any NaN values (best for downstream modeling)
        output_file=output_feature_select_file,
        output_type="parquet",
    )
    feature_select_seconds = time.time() - feature_select_start_time

    # Clear memory
    gc.collect()

    total_seconds = time.time() - plate_start_time
    print(
        f"Preprocessing features completed for {plate_id}! "
        f"(annotate={annotate_seconds:.1f}s, normalize={normalize_seconds:.1f}s, "
        f"feature_select={feature_select_seconds:.1f}s, total={total_seconds:.1f}s)"
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
                f"{normalize_seconds:.1f}",
                f"{feature_select_seconds:.1f}",
                f"{total_seconds:.1f}",
                datetime.now().isoformat(timespec="seconds"),
            ]
        )

