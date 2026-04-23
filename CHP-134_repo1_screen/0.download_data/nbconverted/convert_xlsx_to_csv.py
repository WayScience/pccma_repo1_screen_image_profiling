#!/usr/bin/env python
# coding: utf-8

# # Convert the CHP-134 platemap Excel files to layout-level CSVs

# In[1]:


import pandas as pd
from pathlib import Path


# In[2]:


# Define paths
platemap_file = Path("orig_xlsx_files/CHP-134_REPO1_PlateMaps_20260217.xlsx").resolve(
    strict=True
)
corrected_mapping_file = Path(
    "orig_xlsx_files/PedMap_CHP-134_PlateMapping_corrected.xlsx"
).resolve(strict=True)

# Set output folder for platemaps
metadata_folder = Path("metadata")
metadata_folder.mkdir(exist_ok=True)  # Ensure metadata folder exists


# In[3]:


# Load the source platemap workbook and the corrected plate-layout mapping workbook
platemap_df = pd.read_excel(platemap_file)
corrected_mapping_df = pd.read_excel(corrected_mapping_file)

# Print basic info about the loaded dataframes
print("Platemap DataFrame:")
display(platemap_df.head())
print("\nCorrected Mapping DataFrame:")
display(corrected_mapping_df.head())


# In[4]:


# Basic validation
required_platemap_columns = {"Plate Barcode", "Well Position"}
required_mapping_columns = {"Plate Map Name", "DestinationBarcode"}

# Check for missing columns in both dataframes
missing_platemap_columns = required_platemap_columns - set(platemap_df.columns)
missing_mapping_columns = required_mapping_columns - set(corrected_mapping_df.columns)

# Raise errors if required columns are missing
if missing_platemap_columns:
    raise ValueError(
        f"Missing columns in {platemap_file.name}: {sorted(missing_platemap_columns)}"
    )
if missing_mapping_columns:
    raise ValueError(
        f"Missing columns in {corrected_mapping_file.name}: {sorted(missing_mapping_columns)}"
    )

# Check for duplicate destination barcodes in the corrected mapping dataframe
if corrected_mapping_df["DestinationBarcode"].duplicated().any():
    duplicated_barcodes = corrected_mapping_df.loc[
        corrected_mapping_df["DestinationBarcode"].duplicated(), "DestinationBarcode"
    ].tolist()
    raise ValueError(f"Duplicate destination barcodes found: {duplicated_barcodes}")

# Check for consistency between source plate barcodes and mapped destination barcodes
source_plate_barcodes = set(platemap_df["Plate Barcode"].dropna().unique())
mapped_plate_barcodes = set(
    corrected_mapping_df["DestinationBarcode"].dropna().unique()
)


# In[5]:


# Create one platemap CSV per layout (3 destination plates per layout)
barcode_platemap = []

for plate_map_name, layout_mapping_df in corrected_mapping_df.groupby(
    "Plate Map Name", sort=False
):
    destination_barcodes = layout_mapping_df["DestinationBarcode"].dropna().tolist()

    # The assay plates are missing from the platemap workbook, so use the first real plate
    template_plate_barcode = next(
        (
            barcode
            for barcode in destination_barcodes
            if barcode in source_plate_barcodes
        ),
        None,
    )

    if template_plate_barcode is None:
        raise ValueError(f"No source platemap found for layout {plate_map_name}")

    layout_platemap_df = platemap_df[
        platemap_df["Plate Barcode"] == template_plate_barcode
    ].copy()
    layout_platemap_df["Plate Barcode"] = plate_map_name

    output_file = metadata_folder / f"{plate_map_name}_platemap.csv"
    layout_platemap_df.to_csv(output_file, index=False)
    print(
        f"Saved {output_file} using template plate {template_plate_barcode} "
        f"for {len(destination_barcodes)} destination plates."
    )

    for barcode in destination_barcodes:
        barcode_platemap.append(
            {"Plate Barcode": barcode, "File Name": output_file.name}
        )

unused_source_plates = sorted(source_plate_barcodes - mapped_plate_barcodes)
if unused_source_plates:
    print(f"Source plates not used by the corrected mapping: {unused_source_plates}")


# In[6]:


# Create a DataFrame for barcode and file mapping
barcode_platemap_df = pd.DataFrame(barcode_platemap)
barcode_platemap_file = metadata_folder / "barcode_platemap.csv"
barcode_platemap_df.to_csv(barcode_platemap_file, index=False)
print(f"Saved barcode mapping file: {barcode_platemap_file}")

