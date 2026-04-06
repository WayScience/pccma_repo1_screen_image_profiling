#!/usr/bin/env python
# coding: utf-8

# # Create LoadData CSVs with the paths to IC functions for analysis
# 
# In this notebook, we create LoadData CSVs that contains paths to each channel per image set and associated illumination correction `npy` files per channel for CellProfiler to process. 

# ## Import libraries

# In[1]:


import argparse
import pathlib
import pandas as pd
import re
import os

import sys

sys.path.append("../../utils")
import loaddata_utils as ld_utils
from bandicoot_utils import bandicoot_check


# ## Set paths

# In[2]:


argparse = argparse.ArgumentParser(
    description="Create LoadData CSV files to run CellProfiler on the cluster"
)
argparse.add_argument("--HPC", action="store_true", help="Type of compute to run on")

# Parse arguments
args = argparse.parse_args(args=sys.argv[1:] if "ipykernel" not in sys.argv[0] else [])
HPC = args.HPC

print(f"HPC: {HPC}")


# In[3]:


# Set the index directory based on whether HPC is used or not
if HPC:
    # Path for index directory to make loaddata csvs though compute cluster (HPC)
    index_directory = pathlib.Path(
        "/scratch/alpine/jtomkinson@xsede.org/ALSF_screen_data/SK-N-AS_repo1_screen"
    )
else:
    # Find root directory of the project
    root_dir = pathlib.Path().resolve()

    image_base_dir = bandicoot_check(
        pathlib.Path(os.path.expanduser("~/mnt/bandicoot")).resolve(), root_dir
    )
    # Path for index directory  to make loaddata csv locally
    index_directory = pathlib.Path(
        f"{image_base_dir}/PCCMA_data/SK-N-AS_repo1_screen/"
    ).resolve(strict=True)

# Set all paths that are common to both HPC and local
config_dir_path = pathlib.Path(
    "../1.illumination_correction/load_data_config/"
).resolve(strict=True)
output_csv_dir = pathlib.Path("./loaddata_csvs/").absolute()
output_csv_dir.mkdir(parents=True, exist_ok=True)
illum_directory = pathlib.Path("../1.illumination_correction/illum_directory/").resolve(
    strict=True
)

# Find all 'Images' folders within the directory
images_folders = list(index_directory.rglob("Images"))
print(f"Found {len(images_folders)} 'Images' folders in the directory.")


# ## Create LoadData CSVs with illumination functions for all data

# In[4]:


config_path = pathlib.Path(config_dir_path / "config.yml").resolve(strict=True)

csv_paths = []

for top_level in index_directory.iterdir():
    if not top_level.is_dir():
        continue

    print(f"\nProcessing: {top_level.name}")

    is_screen = "Screen" in top_level.name
    is_row_o_repeat = "Row O Repeat" in top_level.name

    for subfolder in top_level.iterdir():
        if not subfolder.is_dir():
            continue

        br_id = subfolder.name.split("__")[0]

        match = re.search(r"(BR\d+)", br_id)
        if not match:
            print(f"Skipping {subfolder} (no BR ID found)")
            continue

        plate_id = match.group(1)

        print(f"  Processing plate: {plate_id}")

        # Output paths (function writes to BOTH)
        path_to_output_csv = (
            output_csv_dir / f"{plate_id}_loaddata_original.csv"
        ).absolute()

        path_to_output_with_illum_csv = (
            output_csv_dir / f"{plate_id}_loaddata_with_illum.csv"
        ).absolute()

        illum_output_path = (illum_directory / plate_id).absolute().resolve(strict=True)

        # Run loaddata (function handles writing both outputs)
        ld_utils.create_loaddata_illum_csv(
            index_directory=subfolder / "Images",
            config_path=config_path,
            path_to_output=path_to_output_csv,
            illum_directory=illum_output_path,
            plate_id=plate_id,
            illum_output_path=path_to_output_with_illum_csv,
        )

        # ---- Row O filtering ONLY for Screen plates ----
        df = pd.read_csv(path_to_output_with_illum_csv)

        if is_screen:
            df_filtered = df[~df["Metadata_Well"].astype(str).str.startswith("O")]

            assert len(df_filtered) < len(
                df
            ), f"{plate_id} (Screen) expected Row O wells to be removed."

            assert (
                not df_filtered["Metadata_Well"].astype(str).str.startswith("O").any()
            ), f"{plate_id} (Screen) still contains Row O wells."

            df_filtered.to_csv(path_to_output_with_illum_csv, index=False)

        elif is_row_o_repeat:
            assert (
                df["Metadata_Well"].astype(str).str.startswith("O").any()
            ), f"{plate_id} (Row O Repeat) missing Row O wells."

        print(f"  Created LoadData CSVs for {plate_id}")

        csv_paths.append(path_to_output_with_illum_csv)

# ---- Final sanity check ----
unique_plates = {p.stem.split("_")[0] for p in csv_paths}

print(f"\nTotal plates processed: {len(unique_plates)}")

if len(unique_plates) != 29:
    raise ValueError(f"Expected 29 plates, but got {len(unique_plates)}")

print("All LoadData CSVs with illumination created successfully!")

