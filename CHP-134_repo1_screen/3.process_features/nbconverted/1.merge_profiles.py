#!/usr/bin/env python
# coding: utf-8

# # Merge row-batch profiles into one profile per plate
# 
# CHP-134 was processed in 16 row batches per plate (`row_A`...`row_P`), and `0.convert_cytotable.ipynb` converts each plate/row-batch SQLite file into its own parquet file. This notebook concatenates all row-batch parquet files belonging to a plate into a single merged parquet profile per plate, ready for downstream processing (e.g. normalization, feature selection).

# ## Import libraries

# In[1]:


import os
import pathlib
import re
import sys
import traceback

import pandas as pd


# ## Set paths and variables

# In[2]:


# Where to read converted row-batch profiles from and save merged plate profiles to.
# Set to "local" to use this repo's data directory, or "bandicoot" to use the
# bandicoot network mount instead (PCCMA_data/CHP-134_repo1_screen_profiles).
save_location = "local"

if save_location == "bandicoot":
    output_dir = pathlib.Path(
        "~/mnt/bandicoot/PCCMA_data/CHP-134_repo1_screen_profiles"
    ).expanduser()
else:
    output_dir = pathlib.Path("data")

# Directory with the per-plate/row-batch parquet files produced by 0.convert_cytotable.ipynb
converted_dir = output_dir / "converted_profiles"

if not converted_dir.exists():
    raise FileNotFoundError(f"The converted profiles path {converted_dir} does not exist.")

print(f"Reading converted row-batch profiles from: {converted_dir}")

# Directory to save one merged parquet profile per plate
merged_dir = output_dir / "merged_profiles"
merged_dir.mkdir(parents=True, exist_ok=True)

print(f"Saving merged plate profiles to: {merged_dir}")

# Group the converted row-batch files by plate ID (files are named
# "<plate_id>_row_<row>_converted.parquet" by 0.convert_cytotable.ipynb)
row_batch_pattern = re.compile(r"^(?P<plate_id>.+)_row_(?P<row>[A-Za-z0-9]+)_converted\.parquet$")

plate_to_files = {}
for file_path in sorted(converted_dir.glob("*_converted.parquet")):
    match = row_batch_pattern.match(file_path.name)
    if not match:
        print(f"Skipping unrecognized file: {file_path.name}")
        continue
    plate_to_files.setdefault(match.group("plate_id"), []).append(file_path)

plate_to_files = dict(sorted(plate_to_files.items()))

print(
    f"Found {sum(len(files) for files in plate_to_files.values())} row-batch files "
    f"across {len(plate_to_files)} plates to merge."
)


# In[ ]:


# OVERWRITE - if set (1/true/yes), re-merge and overwrite a plate's merged
#             profile even if it already exists. Leave unset (the default) to
#             skip plates that are already merged. Set by run_pipeline.sh so it
#             can control this without papermill.
overwrite = os.environ.get("OVERWRITE", "").strip().lower() in ("1", "true", "yes")


# ## Merge row batches into one profile per plate
# 
# Already-merged plates are skipped by default, so the notebook can safely be re-run if interrupted. Set the `OVERWRITE` env var to re-merge and overwrite plates that already have output.
# 

# In[3]:


failed_plates = []

for i, (plate_id, file_paths) in enumerate(plate_to_files.items(), start=1):
    merged_path = merged_dir / f"{plate_id}.parquet"

    if merged_path.exists() and not overwrite:
        print(f"[{i}/{len(plate_to_files)}] Skipping {plate_id} (already merged)")
        continue

    print(f"[{i}/{len(plate_to_files)}] Merging {plate_id} ({len(file_paths)} row batches)...")

    try:
        # Verify we have all rows A-P for this plate and sort file paths by row order
        file_rows = []
        for fp in file_paths:
            m = row_batch_pattern.match(fp.name)
            if not m:
                raise ValueError(f"Unrecognized file name for plate {plate_id}: {fp.name}")
            file_rows.append((fp, m.group("row").upper()))

        expected_rows = [chr(ord("A") + i) for i in range(16)]  # A..P
        found_rows = [r for _, r in file_rows]
        missing = [r for r in expected_rows if r not in found_rows]
        if missing:
            raise ValueError(f"Missing rows for plate {plate_id}: {missing}")

        order = {r: i for i, r in enumerate(expected_rows)}
        sorted_file_paths = [fp for fp, _ in sorted(file_rows, key=lambda x: order.get(x[1], 999))]

        # Merge the row-batch parquet files for this plate into a single parquet file
        batch_dfs = [pd.read_parquet(fp) for fp in sorted_file_paths]
        merged_df = pd.concat(batch_dfs, ignore_index=True)

        # Confirm no rows were dropped or duplicated during the concatenation
        expected_rows = sum(len(batch_df) for batch_df in batch_dfs)
        assert len(merged_df) == expected_rows, "Row count mismatch after concatenation"

        # Save the merged plate profile to a parquet file
        merged_df.to_parquet(merged_path, index=False)

        print(f"[{i}/{len(plate_to_files)}] Done: {plate_id} -> {merged_df.shape}")

    except Exception:  # noqa: BLE001
        print(f"[{i}/{len(plate_to_files)}] FAILED: {plate_id}")
        traceback.print_exc()
        failed_plates.append(plate_id)
        # Remove any partially written output so a rerun doesn't skip it
        if merged_path.exists():
            merged_path.unlink()
        continue

print(
    f"\nFinished: {len(plate_to_files) - len(failed_plates)}/{len(plate_to_files)} "
    "plates merged successfully."
)
if failed_plates:
    print("Failed plate merges:")
    for plate_id in failed_plates:
        print(f"  {plate_id}")

    # A silent partial merge is worse than a loud failure: exit nonzero so
    # callers (e.g. run_pipeline.sh) can detect an incomplete merge step.
    sys.exit(1)

