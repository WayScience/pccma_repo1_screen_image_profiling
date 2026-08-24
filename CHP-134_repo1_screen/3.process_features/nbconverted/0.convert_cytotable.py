#!/usr/bin/env python
# coding: utf-8

# # Convert SQLite outputs to parquet files with cytotable
# 
# CHP-134 was processed in 16 row batches per plate (`row_A`...`row_P`), so each plate's SQLite output directory contains one subdirectory per row batch, each holding its own SQLite file (all < 2 GB). This notebook runs locally and converts each plate/row-batch SQLite file to parquet sequentially, one at a time, producing 432 parquet files (27 plates x 16 row batches).

# ## Import libraries

# In[1]:


import argparse
import logging
import pathlib
import traceback

import pandas as pd

# cytotable will merge objects from SQLite file into single cells and save as parquet file
from cytotable import convert, presets
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.providers import LocalProvider

# Set the logging level to a higher level to avoid outputting unnecessary errors from config file in convert function
logging.getLogger().setLevel(logging.ERROR)

# check if in a jupyter notebook
try:
    cfg = get_ipython().config
    in_notebook = True
except NameError:
    in_notebook = False


# ## Set paths and variables

# In[2]:


if not in_notebook:
    print("Running as script")
    # set up arg parser
    parser = argparse.ArgumentParser(
        description="CytoTable conversion for CHP-134 REPO1 screen (per plate/row-batch SQLite file)"
    )

    parser.add_argument(
        "--plate_id",
        type=str,
        required=False,
        help="Plate ID (SQLite output subdirectory name) to convert. If omitted, all plates are converted.",
    )
    parser.add_argument(
        "--batch_label",
        type=str,
        required=False,
        help="Row batch label (SQLite output subdirectory name, e.g. row_A) to convert. Requires --plate_id. If omitted, all row batches for the given plate(s) are converted.",
    )

    args = parser.parse_args()
    if args.batch_label and not args.plate_id:
        parser.error("--batch_label requires --plate_id to also be set.")

    plate_id_filter = args.plate_id
    batch_label_filter = args.batch_label
else:
    print("Running in a notebook")
    # Convert every plate/row-batch found in the source directory. Set these
    # to a specific plate_id/batch_label (e.g. "BR00149332" / "row_A") to
    # convert a single SQLite file instead.
    plate_id_filter = None
    batch_label_filter = None


# In[3]:


# preset configurations based on typical CellProfiler outputs
preset = "cellprofiler_sqlite_pycytominer"

# update preset to include both the site metadata, cell counts, and PathName columns
joins = presets.config["cellprofiler_sqlite_pycytominer"]["CONFIG_JOINS"].replace(
    "Image_Metadata_Well,",
    "Image_Metadata_Well, Image_Metadata_Site, Image_Count_Cells, Image_Metadata_Row, Image_Metadata_Col, ",
)

# Add the PathName columns separately
joins = joins.replace(
    "COLUMNS('Image_FileName_.*'),",
    "COLUMNS('Image_FileName_.*'),\n COLUMNS('Image_PathName_.*'),",
)

# type of file output from cytotable (currently only parquet)
dest_datatype = "parquet"

# Set path to directory with SQLite files depending on where this notebook is running.
# On Alpine (HPC), source data lives on the PetaLibrary "koala" mount. Otherwise
# (local), source data lives on the bandicoot network mount. Each plate subdirectory
# contains one folder per row batch (row_A, row_B, ..., row_P), each holding its own
# SQLite file.
alpine_scratch_path = pathlib.Path("/scratch/alpine")

if alpine_scratch_path.exists():
    sqlite_dir = pathlib.Path(
        "/pl/active/koala/ALSF_screen_data/CHP-134_repo1_screen_outputs/SQLite_outputs"
    )
else:
    sqlite_dir = pathlib.Path(
        "~/mnt/bandicoot/PCCMA_data/CHP-134_repo1_screen_outputs/SQLite_outputs"
    ).expanduser()

if not sqlite_dir.exists():
    raise FileNotFoundError(f"The SQLite source path {sqlite_dir} does not exist.")

print(f"Reading SQLite files from: {sqlite_dir}")

# directory for processed data
output_dir = pathlib.Path("data")
output_dir.mkdir(parents=True, exist_ok=True)

# Discover every (plate, row batch) combination present, since each plate was
# processed in 16 row batches (one SQLite file per plate/row batch).
batch_pairs = sorted(
    (plate_path.name, batch_path.name)
    for plate_path in sqlite_dir.iterdir()
    if plate_path.is_dir()
    for batch_path in plate_path.iterdir()
    if batch_path.is_dir()
)

# Optionally restrict to a single plate and/or row batch (e.g. for manual reruns)
if plate_id_filter:
    batch_pairs = [pair for pair in batch_pairs if pair[0] == plate_id_filter]
    if batch_label_filter:
        batch_pairs = [pair for pair in batch_pairs if pair[1] == batch_label_filter]

plate_names = sorted({pair[0] for pair in batch_pairs})

# print the plate names and how many plate/row-batch SQLite files there are (confirmation)
print(
    f"There are {len(batch_pairs)} plate/row-batch SQLite files to convert "
    f"across {len(plate_names)} plates. Below are the plate names:"
)
for name in plate_names:
    print(name)


# ## Convert SQLite to parquet files, then clean up columns
# 
# Runs sequentially (one plate/row-batch at a time) since this is a local run. Already-converted outputs are skipped, so the notebook can safely be re-run if interrupted.
# 
# > Note: We do not run this code cell in notebook, only in script for stability.

# In[ ]:


# Directory with converted profiles
converted_dir = output_dir / "converted_profiles"
converted_dir.mkdir(parents=True, exist_ok=True)

# Define the list of columns to prioritize and prefix
prioritized_columns = [
    "Nuclei_Location_Center_X",
    "Nuclei_Location_Center_Y",
    "Cells_Location_Center_X",
    "Cells_Location_Center_Y",
    "Image_Count_Cells",
]

failed_batches = []

for i, (plate_id, batch_label) in enumerate(batch_pairs, start=1):
    output_path = converted_dir / f"{plate_id}_{batch_label}_converted.parquet"

    if output_path.exists():
        print(f"[{i}/{len(batch_pairs)}] Skipping {plate_id} {batch_label} (already converted)")
        continue

    source_path = sqlite_dir / plate_id / batch_label
    print(f"[{i}/{len(batch_pairs)}] Converting {plate_id} {batch_label}...")

    try:
        convert(
            source_path=str(source_path),
            dest_path=str(output_path),
            dest_datatype=dest_datatype,
            preset=preset,
            joins=joins,
            chunk_size=10000,
            parsl_config=Config(
                executors=[
                    HighThroughputExecutor(
                        label="local_htex",
                        # One CPU per Parsl worker so that chunk-level tasks
                        # don't oversubscribe this workstation's cores. Each
                        # plate/row-batch SQLite file is converted one at a
                        # time, so the full local core count is available.
                        cores_per_worker=1,
                        # This local workstation has 16 cores available.
                        max_workers_per_node=16,
                        # Assign workers to distinct CPU cores to reduce contention.
                        cpu_affinity="block",
                        provider=LocalProvider(
                            init_blocks=1,
                            max_blocks=1,
                        ),
                    )
                ],
                # Keep Parsl logs separate for each plate/row batch.
                run_dir=f"runinfo/{plate_id}/{batch_label}",
                strategy="none",
            ),
        )

        # Load the newly converted profile back in to clean up columns
        converted_df = pd.read_parquet(output_path)

        # If any, drop rows where "Metadata_ImageNumber" is NaN (artifact of cytotable)
        converted_df = converted_df.dropna(subset=["Metadata_ImageNumber"])

        # Rearrange columns and add "Metadata" prefix in one line
        converted_df = converted_df[
            prioritized_columns
            + [col for col in converted_df.columns if col not in prioritized_columns]
        ].rename(
            columns=lambda col: "Metadata_" + col if col in prioritized_columns else col
        )

        # assert that there are column names with PathName in the dataset
        assert any("PathName" in col for col in converted_df.columns)

        # Assert that Metadata_Row and Metadata_Col are present for downstream QC
        assert {"Image_Metadata_Row", "Image_Metadata_Col"}.issubset(
            converted_df.columns
        ), "Missing required Metadata columns: Row and/or Col"

        # Save the processed DataFrame as Parquet in the same path
        converted_df.to_parquet(output_path, index=False)

        print(
            f"[{i}/{len(batch_pairs)}] Done: {plate_id} {batch_label} -> {converted_df.shape}"
        )

    except Exception:
        print(f"[{i}/{len(batch_pairs)}] FAILED: {plate_id} {batch_label}")
        traceback.print_exc()
        failed_batches.append((plate_id, batch_label))
        # Remove any partially written output so a rerun doesn't skip it
        if output_path.exists():
            output_path.unlink()
        continue

print(
    f"\nFinished: {len(batch_pairs) - len(failed_batches)}/{len(batch_pairs)} "
    "plate/row-batch files converted successfully."
)
if failed_batches:
    print("Failed plate/row-batch conversions:")
    for plate_id, batch_label in failed_batches:
        print(f"  {plate_id} {batch_label}")


# **To confirm the number of single cells is correct for each converted profile, please use any database browser software to see if the number of rows in the "Per_Cells" compartment of the corresponding SQLite file matches the number of rows in the parquet file.**
