#!/usr/bin/env python
# coding: utf-8

# # Aggregate bulk profiles
# 
# Aggregates each plate's cleaned single cells to one profile per well (the median of the cells) and writes it to `data/bulk_profiles/<plate>_bulk_aggregated.parquet`. These per-plate, un-normalized profiles are the input of `3b.pooled_bulk_processing.ipynb`, where normalization and feature selection happen once over the whole screen.

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
from pycytominer import aggregate, annotate
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

    # Output file path
    output_aggregated_file = str(output_dir / f"{plate_id}_bulk_aggregated.parquet")

    # Already aggregated, so this plate can safely be skipped on a rerun
    # (unless OVERWRITE is set, in which case it's reprocessed regardless)
    if pathlib.Path(output_aggregated_file).exists() and not overwrite:
        print(f"Skipping {plate_id}: already processed (found {output_aggregated_file})")
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
        output_filename=output_aggregated_file,
        output_type="parquet",
    )
    aggregate_seconds = time.time() - aggregate_start_time

    # Clear memory
    del single_cell_df, qc_df, platemap_df, annotated_df
    gc.collect()
    del aggregated_df
    gc.collect()

    total_seconds = time.time() - plate_start_time
    print(
        f"Bulk aggregation completed for {plate_id}! "
        f"(annotate={annotate_seconds:.1f}s, aggregate={aggregate_seconds:.1f}s, "
        f"total={total_seconds:.1f}s)"
    )

    # Record timing so per-plate performance can be compared across runs
    write_header = not timing_log_path.exists()
    with open(timing_log_path, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["plate_id", "annotate_seconds", "aggregate_seconds", "total_seconds", "timestamp"])
        writer.writerow(
            [
                plate_id,
                f"{annotate_seconds:.1f}",
                f"{aggregate_seconds:.1f}",
                f"{total_seconds:.1f}",
                datetime.now(timezone.utc).isoformat(timespec="seconds"),
            ]
        )

