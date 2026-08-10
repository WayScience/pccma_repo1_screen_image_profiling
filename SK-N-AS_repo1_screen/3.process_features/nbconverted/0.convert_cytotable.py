#!/usr/bin/env python
# coding: utf-8

# # Convert SQLite outputs to parquet files with cytotable
# 
# Notebook will not be executed as we will run on HPC.

# ## Import libraries

# In[ ]:


import argparse
import logging
import os
import pathlib

import pandas as pd

# cytotable will merge objects from SQLite file into single cells and save as parquet file
from cytotable import convert, presets

# used to scope each plate's Parsl run_dir separately (see conversion cell below)
from cytotable.utils import CYTOTABLE_THREAD_EXECUTOR_LABEL
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.executors import ThreadPoolExecutor as ParslThreadPoolExecutor

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
        description="CytoTable conversion for SK-N-AS REPO1 screen"
    )

    parser.add_argument(
        "--plate_id",
        type=str,
        required=True,
        help="Plate ID (SQLite output subdirectory name) to convert",
    )

    args = parser.parse_args()
    plate_id = args.plate_id
else:
    print("Running in a notebook")
    plate_id = "BR00148919"


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
# (local), source data lives on the bandicoot network mount.
alpine_scratch_path = pathlib.Path("/scratch/alpine")

if alpine_scratch_path.exists():
    sqlite_dir = pathlib.Path(
        "/pl/active/koala/ALSF_screen_data/SK-N-AS_repo1_profiles/SQLite_outputs"
    )
else:
    sqlite_dir = pathlib.Path(
        "~/mnt/bandicoot/PCCMA_data/SK-N-AS_repo1_screen/SQLite_outputs"
    ).expanduser()

if not sqlite_dir.exists():
    raise FileNotFoundError(f"The SQLite source path {sqlite_dir} does not exist.")

print(f"Reading SQLite files from: {sqlite_dir}")

# directory for processed data
output_dir = pathlib.Path("data")
output_dir.mkdir(parents=True, exist_ok=True)

plate_names = [
    file_path.stem for file_path in sqlite_dir.iterdir() if file_path.is_dir()
]

# print the plate names and how many plates there are (confirmation)
print(f"There are {len(plate_names)} plates in this dataset. Below are the names:")
for name in plate_names:
    print(name)


# ## Convert SQLite to parquet files

# In[ ]:


file_path = sqlite_dir / plate_id
output_path = pathlib.Path(
    f"{output_dir}/converted_profiles/{plate_id}_converted.parquet"
)

# use SLURM-allocated CPU count, not full node count
n_workers = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
print(f"Using {n_workers} parsl workers (based on allocated CPUs)")

# HTEX only -- see cytomining/CytoTable#75, dispatch stall on BR00148919
parsl_config = Config(
    run_dir=f"runinfo/{plate_id}",
    executors=[
        HighThroughputExecutor(
            label="htex_default_for_cytotable",
            max_workers_per_node=n_workers,
        ),
    ],
)

print("Starting conversion with cytotable for plate:", plate_id)
convert(
    source_path=str(file_path),
    dest_path=str(output_path),
    dest_datatype=dest_datatype,
    preset=preset,
    joins=joins,
    chunk_size=30000,
    parsl_config=parsl_config,
)

print(f"Plate {plate_id} has been converted with cytotable!")


# # Load in converted profiles to update

# In[ ]:


# Directory with converted profiles
converted_dir = pathlib.Path(f"{output_dir}/converted_profiles")

# Define the list of columns to prioritize and prefix
prioritized_columns = [
    "Nuclei_Location_Center_X",
    "Nuclei_Location_Center_Y",
    "Cells_Location_Center_X",
    "Cells_Location_Center_Y",
    "Image_Count_Cells",
]

# Load the DataFrame from the Parquet file
file_path = converted_dir / f"{plate_id}_converted.parquet"
converted_df = pd.read_parquet(file_path)

# If any, drop rows where "Metadata_ImageNumber" is NaN (artifact of cytotable)
converted_df = converted_df.dropna(subset=["Metadata_ImageNumber"])

# Rearrange columns and add "Metadata" prefix in one line
converted_df = converted_df[
    prioritized_columns
    + [col for col in converted_df.columns if col not in prioritized_columns]
].rename(columns=lambda col: "Metadata_" + col if col in prioritized_columns else col)

# assert that there are column names with PathName in the dataset
assert any("PathName" in col for col in converted_df.columns)

# Assert that Metadata_Row and Metadata_Col are present for downstream QC
assert {"Image_Metadata_Row", "Image_Metadata_Col"}.issubset(
    converted_df.columns
), "Missing required Metadata columns: Row and/or Col"

# Save the processed DataFrame as Parquet in the same path
converted_df.to_parquet(file_path, index=False)

# print shape and head of dataset
print(converted_df.shape)
converted_df.head()


# **To confirm the number of single cells is correct above, please use any database browser software to see if the number of rows in the "Per_Cells" compartment matches the number of rows in the data frame.**
