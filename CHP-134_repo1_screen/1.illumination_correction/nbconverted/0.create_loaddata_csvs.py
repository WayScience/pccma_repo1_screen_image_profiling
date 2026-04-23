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

# check if in a jupyter notebook
try:
    cfg = get_ipython().config
    in_notebook = True
except NameError:
    in_notebook = False


# ## Set paths

# In[2]:


if not in_notebook:
    print("Running as script")
    # set up arg parser
    parser = argparse.ArgumentParser(
        description="Generate LoadData CSV file for CellProfiler from a given input CSV containing image paths and metadata"
    )

    parser.add_argument(
        "--index_directory",
        type=str,
        help="Path to the directory containing image folders",
    )

    args = parser.parse_args()
    index_directory = pathlib.Path(args.index_directory).resolve(strict=True)
else:
    print("Running in a notebook")
    # Set the root directory and image base directory for notebook execution
    root_dir = pathlib.Path().resolve()

    image_base_dir = bandicoot_check(
        pathlib.Path(os.path.expanduser("~/mnt/bandicoot")).resolve(), root_dir
    )
    index_directory = pathlib.Path(
        f"{image_base_dir}/PCCMA_data/CHP-134_repo1_screen/"
    ).resolve(strict=True)

# Paths for parameters to make loaddata csv
config_dir_path = pathlib.Path("./load_data_config").absolute()
output_csv_dir = pathlib.Path("./loaddata_csvs/")
output_csv_dir.mkdir(parents=True, exist_ok=True)

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


# ## Create LoadData CSVs for all data

# In[3]:


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

    # Create output path with a filename-safe version of the plate name
    path_to_output_csv = (output_csv_dir / f"{plate_name}_loaddata.csv").absolute()

    # Run loaddata creation
    ld_utils.create_loaddata_csv(
        index_directory=images_dir,
        config_path=config_path,
        path_to_output=path_to_output_csv,
    )

    print(f"Created LoadData CSV for {plate_name} at {path_to_output_csv}")

    # ---- Validate CSV row count ----
    try:
        df = pd.read_csv(path_to_output_csv)
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

