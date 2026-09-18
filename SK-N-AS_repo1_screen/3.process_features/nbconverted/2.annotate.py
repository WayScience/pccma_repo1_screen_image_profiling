#!/usr/bin/env python
# coding: utf-8

# # Annotate single-cell profiles
# 
# > Shared step used by both the bulk (`3.bulk_processing.ipynb`) and single-cell (`4.single_cell_normalize.ipynb`) branches downstream. Splitting this out means `annotate()` only runs once per plate instead of twice, and its memory-heavy step (holding both the raw profile and its annotated copy at once, each tens of GB) is isolated to its own process rather than sharing one with aggregation or normalization.

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
from pycytominer import annotate
from pycytominer.cyto_utils import output


# ## Set paths and variables

# In[ ]:


# Directory containing the converted single-cell profile parquet per plate,
# produced by 0.convert_cytotable.ipynb. On Alpine (HPC), converted profiles
# live on the PetaLibrary "koala" mount; locally they live on the external
# drive since each plate's profile is tens of GB. Mirrors the same
# Alpine-detection branch used for sqlite_dir in 0.convert_cytotable.ipynb.
alpine_scratch_path = pathlib.Path("/scratch/alpine")

if alpine_scratch_path.exists():
    converted_dir = pathlib.Path("/pl/active/koala/ALSF_screen_data/SK-N-AS_repo1_profiles/converted_profiles")
else:
    converted_dir = pathlib.Path("/media/18tbdrive2/SK-N-AS_repo1_profiles/converted_profiles")

if not converted_dir.exists():
    raise FileNotFoundError(f"The converted profiles path {converted_dir} does not exist.")

# Directory containing per-plate QC annotation files from 1.single_cell_qc.ipynb
qc_dir = pathlib.Path("./data/qc_results")

# Output path for annotated profiles. Shared by both the bulk and single-cell
# downstream steps, so annotate() only needs to run once per plate.
output_dir = pathlib.Path("./data/annotated_profiles")
output_dir.mkdir(parents=True, exist_ok=True)

# path for platemap directory
platemap_dir = pathlib.Path("../0.download_data/metadata")

# load in barcode platemap
barcode_platemap = pd.read_csv(platemap_dir / "barcode_platemap.csv")

# extract the plate names from the converted profile file names
plate_names = sorted(
    file.stem.removesuffix("_converted")
    for file in converted_dir.glob("*_converted.parquet")
)

# Rename map applied after annotate() so both downstream branches see a
# consistent site column name (CytoTable's raw output calls it
# Image_Metadata_Site; pycytominer's own convention is Metadata_Site).
column_name_mapping = {
    "Image_Metadata_Site": "Metadata_Site",
}

plate_names


# ## Set dictionary with plates to process

# In[ ]:


# Create plate info dictionary
plate_info_dictionary = {
    plate_id: {
        "profile_path": str(converted_dir / f"{plate_id}_converted.parquet"),

        # QC annotations are produced per-plate by 1.single_cell_qc.ipynb
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
                    barcode_platemap["Plate_Barcode"] == plate_id, "File_Name"
                ].values[0]
            )
            if plate_id in barcode_platemap["Plate_Barcode"].values
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
#               process so memory doesn't accumulate across plates -- and so
#               the raw profile and its annotated copy, both full single-cell
#               scale, aren't held alongside another plate's leftovers). Leave
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


# ## Annotate profiles

# In[ ]:


timing_log_path = output_dir / "timing_log.csv"

for plate_id, info in plate_info_dictionary.items():
    if info["qc_path"] is None:
        print(
            f"Skipping {plate_id}: no QC annotations yet "
            "(run 1.single_cell_qc.ipynb for this plate first)"
        )
        continue
    if info["platemap_path"] is None:
        print(f"Skipping {plate_id}: no platemap found in barcode_platemap.csv")
        continue

    output_annotated_file = str(output_dir / f"{plate_id}_annotated.parquet")

    # Already annotated, so this plate can safely be skipped on a rerun
    # (unless OVERWRITE is set, in which case it's reprocessed regardless)
    if pathlib.Path(output_annotated_file).exists() and not overwrite:
        print(f"Skipping {plate_id}: already annotated (found {output_annotated_file})")
        continue

    print(f"Annotating {plate_id}")
    plate_start_time = time.time()

    # Load single-cell profile, its QC annotations, and the platemap. The
    # platemap's column names are already clean (no spaces/parentheses --
    # see 0.download_data/convert_xlsx_to_csv.ipynb), so no rename is needed
    # here the way there used to be.
    single_cell_df = pd.read_parquet(info["profile_path"])
    qc_df = pd.read_parquet(info["qc_path"])
    platemap_df = pd.read_csv(info["platemap_path"])

    # A few plates' CytoTable conversion added a "Metadata_" prefix to the
    # Nuclei location columns (e.g. Metadata_Nuclei_Location_Center_X instead
    # of Nuclei_Location_Center_X) that most plates don't have. Normalize
    # those back to the canonical (unprefixed) name -- matching the same
    # fallback 1.single_cell_qc.ipynb applies -- so join_keys below works the
    # same regardless of which naming convention this plate's conversion
    # produced.
    single_cell_df = single_cell_df.rename(
        columns={
            f"Metadata_{col}": col
            for col in ("Nuclei_Location_Center_X", "Nuclei_Location_Center_Y")
            if col not in single_cell_df.columns
            and f"Metadata_{col}" in single_cell_df.columns
        }
    )

    # Define the columns from the external QC file that indicate poor-quality segmentations
    cqc_cols = [col for col in qc_df.columns if col.startswith("Metadata_cqc_")]
    join_keys = [
        "Image_Metadata_Well", "Image_Metadata_Site",
        "Nuclei_Location_Center_X", "Nuclei_Location_Center_Y",
    ]
    assert not qc_df.duplicated(subset=join_keys).any(), f"{plate_id}: QC file has duplicate rows for the same cell"
    assert not single_cell_df.duplicated(subset=join_keys).any(), f"{plate_id}: profile file has duplicate rows for the same cell"

    # Drop poor-quality segmentations BEFORE annotating, not after (for memory reasons)
    qc_flags = single_cell_df[join_keys].merge(
        qc_df[join_keys + cqc_cols], on=join_keys, how="left"
    )
    assert qc_flags[cqc_cols].isna().sum().sum() == 0, f"{plate_id}: some cells have no matching QC annotation"
    is_poor_quality = qc_flags[cqc_cols].any(axis=1)
    print(f"  Dropping {is_poor_quality.sum()} / {len(single_cell_df)} poor-quality segmentations")
    single_cell_df = single_cell_df.loc[~is_poor_quality].reset_index(drop=True)
    del qc_df, qc_flags, is_poor_quality

    # Annotation
    annotated_df = annotate(
        profiles=single_cell_df,
        platemap=platemap_df,
        join_on=["Metadata_Well_Position", "Image_Metadata_Well"],
    )
    annotated_df.rename(columns=column_name_mapping, inplace=True)

    # Assert "Metadata_Site" is now present after the rename() call above
    assert "Metadata_Site" in annotated_df.columns, f"{plate_id}: Metadata_Site column missing after rename()"

    output(
        df=annotated_df,
        output_filename=output_annotated_file,
        output_type="parquet",
    )
    annotate_seconds = time.time() - plate_start_time

    # Clear memory
    del single_cell_df, platemap_df, annotated_df
    gc.collect()

    print(f"Annotation completed for {plate_id}! (annotate={annotate_seconds:.1f}s)")

    # Record timing so per-plate performance can be compared across runs
    write_header = not timing_log_path.exists()
    with open(timing_log_path, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["plate_id", "annotate_seconds", "timestamp"])
        writer.writerow(
            [
                plate_id,
                f"{annotate_seconds:.1f}",
                datetime.now(timezone.utc).isoformat(timespec="seconds"),
            ]
        )

