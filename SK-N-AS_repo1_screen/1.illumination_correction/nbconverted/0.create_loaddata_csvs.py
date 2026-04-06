#!/usr/bin/env python
# coding: utf-8

# # Create LoadData CSVs to use for illumination correction
# 
# In this notebook, we create a LoadData CSV that contains paths to each channel per image set for CellProfiler to process. 
# We can use this LoadData CSV to run illumination correction (IC) pipeline that saves IC functions in `npy` file format and extract image quality metrics.

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

# Paths for parameters to make loaddata csv
config_dir_path = pathlib.Path("./load_data_config").absolute()
output_csv_dir = pathlib.Path("./loaddata_csvs/")
output_csv_dir.mkdir(parents=True, exist_ok=True)

# Find all 'Images' folders within the directory
images_folders = list(index_directory.rglob("Images"))
print(f"Found {len(images_folders)} 'Images' folders in {index_directory}")


# ## Create LoadData CSVs for all data

# In[4]:


config_path = pathlib.Path(config_dir_path / "config.yml").resolve(strict=True)

csv_paths = []

# Iterate through top-level folders (REPO1 Screen, REPO1 Row O Repeat)
for top_level in index_directory.iterdir():
    if not top_level.is_dir():
        continue

    is_screen = "Screen" in top_level.name
    is_row_o_repeat = "Row O Repeat" in top_level.name

    # Iterate through BR ID folders inside each top-level folder
    for subfolder in top_level.iterdir():
        if not subfolder.is_dir():
            continue

        images_dir = subfolder / "Images"

        if not images_dir.exists():
            print(f"Skipping {subfolder} (no Images folder)")
            continue

        xml_files = list(images_dir.rglob("*Index*.xml"))
        if not xml_files:
            print(f"Skipping {subfolder} (no Index XML found)")
            continue

        base_name = subfolder.name.split("__")[0]

        match = re.search(r"(BR\d+)", base_name)
        if not match:
            print(f"Skipping {subfolder} (no BR ID found)")
            continue

        plate_name = match.group(1)

        plate_id_in_xml = None
        with open(xml_files[0], "r") as f:
            for line in f:
                match_xml = re.search(r"<PlateID>(.*?)</PlateID>", line)
                if match_xml:
                    plate_id_in_xml = match_xml.group(1)
                    break

        if plate_id_in_xml and plate_id_in_xml != plate_name:
            print(
                f"Skipping {subfolder} (PlateID mismatch: {plate_id_in_xml} != {plate_name})"
            )
            continue

        path_to_output_csv = (output_csv_dir / f"{plate_name}_loaddata.csv").absolute()

        ld_utils.create_loaddata_csv(
            index_directory=images_dir,
            config_path=config_path,
            path_to_output=path_to_output_csv,
        )

        # ---- Apply Row O filtering ONLY for Screen plates ----
        df = pd.read_csv(path_to_output_csv)

        if is_screen:
            # Drop wells starting with "O"
            df_filtered = df[~df["Metadata_Well"].astype(str).str.startswith("O")]

            # Assert we actually removed something
            assert len(df_filtered) < len(df), (
                f"{plate_name} (Screen) expected to drop Row O wells, "
                f"but no rows were removed."
            )

            # Assert no "O" wells remain
            assert (
                not df_filtered["Metadata_Well"].astype(str).str.startswith("O").any()
            ), f"{plate_name} (Screen) still contains Row O wells after filtering."

            df_filtered.to_csv(path_to_output_csv, index=False)

        elif is_row_o_repeat:
            # Assert that Row O wells ARE present
            assert df["Metadata_Well"].astype(str).str.startswith("O").any(), (
                f"{plate_name} (Row O Repeat) expected to contain Row O wells, "
                f"but none were found."
            )

        print(f"Created LoadData CSV for {plate_name} at {path_to_output_csv}")
        csv_paths.append(path_to_output_csv)

# ---- Final assertion ----
if len(csv_paths) != 29:
    raise ValueError(f"Expected 29 CSVs, but got {len(csv_paths)}")

print("All CSVs created successfully!")

