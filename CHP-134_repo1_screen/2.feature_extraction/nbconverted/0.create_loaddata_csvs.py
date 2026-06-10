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
from typing import List

# check if in a jupyter notebook
try:
    cfg = get_ipython().config
    in_notebook = True
except NameError:
    in_notebook = False


# ## Set helper functions

# In[ ]:


# Set expected values for validation for generating manifest of row batches
EXPECTED_ROWS = 3456
EXPECTED_PLATE_ROWS = list("ABCDEFGHIJKLMNOP")


def plate_name_from_csv(path: pathlib.Path) -> str:
    """Collect name of the plate from the CSV file name.

    Args:
        path (pathlib.Path)): Path to the CSV file.

    Returns:
        str: Name of the plate.
    """
    name = path.stem
    suffix = "_loaddata_with_illum"
    return name[: -len(suffix)] if name.endswith(suffix) else name


def create_row_batch_manifest(csv_paths: List[pathlib.Path], output_path: pathlib.Path) -> pd.DataFrame:
    """Create a manifest CSV file that splits plates into batches per well-row for HPC processing.

    Args:
        csv_paths (List[pathlib.Path]): List of paths to the CSV files.
        output_path (pathlib.Path): Path to the output manifest CSV file.

    Returns:
        pd.DataFrame: DataFrame containing the manifest of row batches.
    """
    records = []

    for csv_path in csv_paths:
        df = pd.read_csv(csv_path)
        df["image_set_number"] = range(1, len(df) + 1)
        df["plate_row"] = (
            df["Metadata_Well"].str.extract(r"^([A-Za-z]+)")[0].str.upper()
        )

        for plate_row in EXPECTED_PLATE_ROWS:
            row_df = df[df["plate_row"] == plate_row]
            if row_df.empty:
                records.append(
                    {
                        "plate": plate_name_from_csv(csv_path),
                        "row": plate_row,
                        "batch_label": f"row_{plate_row}",
                        "loaddata_csv": csv_path.resolve(),
                        "first_image_set": "",
                        "last_image_set": "",
                        "image_set_count": 0,
                        "well_count": 0,
                        "is_contiguous": False,
                        "status": "missing_row",
                        "message": "No image sets found for this plate row",
                    }
                )
                continue

            first_image_set = int(row_df["image_set_number"].min())
            last_image_set = int(row_df["image_set_number"].max())
            image_set_count = len(row_df)
            is_contiguous = image_set_count == last_image_set - first_image_set + 1

            records.append(
                {
                    "plate": plate_name_from_csv(csv_path),
                    "row": plate_row,
                    "batch_label": f"row_{plate_row}",
                    "loaddata_csv": csv_path.resolve(),
                    "first_image_set": first_image_set,
                    "last_image_set": last_image_set,
                    "image_set_count": image_set_count,
                    "well_count": row_df["Metadata_Well"].nunique(),
                    "is_contiguous": is_contiguous,
                    "status": "ready" if is_contiguous else "skip_noncontiguous",
                    "message": (
                        ""
                        if is_contiguous
                        else "Image sets for this row are not contiguous"
                    ),
                }
            )

    manifest_df = pd.DataFrame(records)
    manifest_df.to_csv(output_path, index=False)
    ready = (manifest_df["status"] == "ready").sum()
    skipped = (manifest_df["status"] != "ready").sum()
    print(f"Wrote {output_path} with {ready} ready row batches and {skipped} warnings")

    return manifest_df


# ## Set paths

# In[3]:


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
row_batch_manifest_path = pathlib.Path("./row_batch_manifest.csv").absolute()

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

# In[4]:


# Define the one config path to use
config_path = config_dir_path / "config.yml"
# Initialize list to keep track of generated CSV paths for manifest creation
csv_paths = []


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

    # ---------------------------------------------------------------------------------------------
    # FIX: filesystem-safe name for pe2loaddata ONLY (for non-conventional plate names with spaces)
    # ---------------------------------------------------------------------------------------------
    # Create a filesystem-safe version (without space) to avoid issues
    safe_plate_name = plate_name.replace(" ", "_")

    # Create output paths
    path_to_output_csv = (
        output_csv_dir / f"{safe_plate_name}_loaddata_original.csv"
    ).absolute()

    # LoadData CSV output name will not have spaces to avoid issues,
    # but paths to images will still have spaces if present
    path_to_output_with_illum_csv = (
        output_csv_dir / f"{safe_plate_name}_loaddata_with_illum.csv"
    ).absolute()

    # Must use the original plate name (with spaces if present) for the illumination output path
    # to ensure metadata correctness in the final CSV
    illum_output_path = (illum_directory / plate_name).absolute().resolve(strict=True)

    # Run loaddata creation with illumination correction functions
    ld_utils.create_loaddata_illum_csv(
        index_directory=images_dir,
        config_path=config_path,
        path_to_output=path_to_output_csv,
        illum_directory=illum_output_path,
        plate_id=plate_name,  # keep ORIGINAL (important for metadata correctness)
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
        csv_paths.append(path_to_output_with_illum_csv)

        if row_count != EXPECTED_ROWS:
            print(
                f"WARNING: {plate_name} has {row_count} rows "
                f"(expected {EXPECTED_ROWS})"
            )
        else:
            print(f"{plate_name} row count OK ({row_count})")

    except Exception as e:
        print(f"Error reading CSV for {plate_name}: {e}")

if csv_paths:
    create_row_batch_manifest(csv_paths, row_batch_manifest_path)

