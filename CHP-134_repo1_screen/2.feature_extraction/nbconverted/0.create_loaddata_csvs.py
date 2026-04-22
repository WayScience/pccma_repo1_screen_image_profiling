#!/usr/bin/env python
# coding: utf-8

# # Create LoadData CSVs with the paths to IC functions for analysis
# 
# In this notebook, we create LoadData CSVs that contains paths to each channel per image set and associated illumination correction `npy` files per channel for CellProfiler to process. 

# ## Import libraries

# In[ ]:


import argparse
import pathlib
import pandas as pd
import re
import os

import sys

sys.path.append("../../utils")
import loaddata_utils as ld_utils
from bandicoot_utils import bandicoot_check

# check if in a jupyter notebook
try:
    cfg = get_ipython().config
    in_notebook = True
except NameError:
    in_notebook = False


# ## Set paths

# In[ ]:


parser = argparse.ArgumentParser(
    description="Create LoadData CSV files to run CellProfiler on the cluster"
)
parser.add_argument("--HPC", action="store_true", help="Type of compute to run on")

# Parse arguments
args = parser.parse_args(args=sys.argv[1:] if "ipykernel" not in sys.argv[0] else [])
HPC = args.HPC

print(f"HPC: {HPC}")

# Set the index directory based on whether HPC is used or not
if HPC:
    # Path for index directory to make loaddata csvs though compute cluster (HPC)
    index_directory = pathlib.Path(
        "/scratch/alpine/jtomkinson@xsede.org/ALSF_screen_data/CHP-134_repo1_screen"
    )
else:
    print("Running in a notebook/local environment")
    root_dir = pathlib.Path().resolve()

    image_base_dir = bandicoot_check(
        pathlib.Path(os.path.expanduser("~/mnt/bandicoot")).resolve(), root_dir
    )
    index_directory = pathlib.Path(
        f"{image_base_dir}/PCCMA_data/CHP-134_repo1_screen/"
    ).resolve(strict=True)

# Paths for parameters to make loaddata csv
config_dir_path = pathlib.Path(
    "../1.illumination_correction/load_data_config/"
).resolve(strict=True)
output_csv_dir = pathlib.Path("./loaddata_csvs/").absolute()
output_csv_dir.mkdir(parents=True, exist_ok=True)
illum_directory = pathlib.Path("../1.illumination_correction/illum_directory/").resolve(
    strict=True
)

# Recursively find Images folders and print how many plates we are working with
images_folders = []
for current_dir, dirnames, _ in os.walk(index_directory, topdown=True):
    if "Images" in dirnames:
        images_folders.append(pathlib.Path(current_dir) / "Images")
        dirnames.remove("Images")

images_folders = sorted(images_folders)
plate_folders = [images_dir.parent for images_dir in images_folders]
direct_plate_folders = [p for p in plate_folders if p.parent == index_directory]
nested_plate_folders = [p for p in plate_folders if p.parent != index_directory]

print(f"Found {len(images_folders)} Images folders across {len(plate_folders)} plates")
print(f"Direct plate folders under index_directory: {len(direct_plate_folders)}")
print(
    f"Nested plate folders in subdirectories such as reimaged data: {len(nested_plate_folders)}"
)


# ## Create LoadData CSVs with illumination functions for all data

# In[ ]:


# Define the one config path to use
config_path = config_dir_path / "config.yml"

EXPECTED_ROWS = 3456

# Iterate over every discovered plate folder, including nested reimaged-data plates
for subfolder in plate_folders:
    images_dir = subfolder / "Images"

    # Use the exact Index.xml path to avoid scanning a directory full of TIFFs
    xml_file = images_dir / "Index.xml"
    if not xml_file.exists():
        print(f"Skipping {subfolder} (no Index XML found)")
        continue
    print(f"Processing {subfolder} with Index XML: {xml_file.name}")

    # ---- Plate naming logic ----
    base_name = subfolder.name.split("__")[0]

    # Try to extract BR plate ID if present, otherwise use the base name
    match = re.search(r"(BR\d+)", base_name)
    if match:
        plate_name = match.group(1)
        print(f"Using BR plate ID: {plate_name}")
    else:
        plate_name = base_name
        print(f"Using assay plate ID: {plate_name}")

    # ---- Validate PlateID from XML (stream read) ----
    plate_id_in_xml = None
    with open(xml_file, "r") as f:
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

    # Create output paths for the temporary original CSV and final illum CSV
    path_to_output_csv = (
        output_csv_dir / f"{plate_name}_loaddata_original.csv"
    ).absolute()
    path_to_output_with_illum_csv = (
        output_csv_dir / f"{plate_name}_loaddata_with_illum.csv"
    ).absolute()
    illum_output_path = (illum_directory / plate_name).absolute().resolve(strict=True)

    # Run loaddata creation with illumination correction functions
    ld_utils.create_loaddata_illum_csv(
        index_directory=images_dir,
        config_path=config_path,
        path_to_output=path_to_output_csv,
        illum_directory=illum_output_path,
        plate_id=plate_name,
        illum_output_path=path_to_output_with_illum_csv,
    )

    print(
        f"Created LoadData CSV with illumination functions for {plate_name} at "
        f"{path_to_output_with_illum_csv}"
    )

    # ---- Validate final CSV row count ----
    try:
        df = pd.read_csv(path_to_output_with_illum_csv)
        row_count = len(df)

        if row_count != EXPECTED_ROWS:
            print(
                f"WARNING: {plate_name} has {row_count} rows "
                f"(expected {EXPECTED_ROWS})"
            )
        else:
            print(f"{plate_name} row count OK ({row_count})")

    except Exception as e:
        print(f"Error reading CSV for {plate_name}: {e}")

