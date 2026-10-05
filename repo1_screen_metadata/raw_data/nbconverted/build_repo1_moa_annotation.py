#!/usr/bin/env python
# coding: utf-8

# # REPO1 mechanism of action annotation
# 
# Builds `repo1_moa_metadata_annotation.csv`: one row per compound and annotated mechanism of action, with columns `moa` and `compound_id`.
# 
# The `MOA` column of `REPO1_April2025.xlsx` holds comma-separated mechanisms. Three clean-ups are applied: a mechanism whose name contains a comma is kept whole, an en dash that was garbled (in three different ways) when the spreadsheet was written is restored, and spelling variants that differ only in capitalization are merged into the most common spelling.

# In[1]:


import pandas as pd

plate_map_file = "../REPO1_April2025.xlsx"
output_file = "repo1_moa_metadata_annotation.csv"

# a mechanism name that itself contains a comma, which must not be split
comma_in_name = "nuclear factor erythroid derived, like (NRF2) activator"
placeholder = "nuclear factor erythroid derived|like (NRF2) activator"
# the en dash in "serotonin-norepinephrine" is garbled in the spreadsheet, in three different ways (text decoded with the wrong encoding)
garbled_dashes = ["\u201a\u00c4\u00ec", "\u221a\u00a2\u00ac\u00c4\u00ac\u00ec", "\u201a\u00c4\u00f6\u221a\u00a7\u221a\u00a8"]


# ### Helper functions
# 
# All functions used below are defined here.

# In[2]:


def split_moa(text: str, comma_in_name: str, placeholder: str, garbled_dashes: list[str]) -> list[str]:
    """Split the comma-separated MOA text of one compound into its mechanisms.

    Restores the en dash that was garbled in the spreadsheet and keeps the one mechanism whose name contains a comma whole.

    Args:
        text: The MOA text of a compound, for example "CDK inhibitor, EGFR inhibitor".
        comma_in_name: A mechanism name that itself contains a comma, which must not be split.
        placeholder: A stand-in for that name without the comma, used while splitting.
        garbled_dashes: The garbled spellings of the en dash to replace.

    Returns:
        The list of mechanisms, stripped of surrounding spaces.
    """
    for garbled in garbled_dashes:
        text = text.replace(garbled, "\u2013")
    parts = text.replace(comma_in_name, placeholder).split(",")
    return [p.strip().replace(placeholder, comma_in_name) for p in parts]


# In[3]:


plate_map = pd.read_excel(plate_map_file, sheet_name="REPO2025_PlateMap")

# one row per compound (a compound plated in two wells has identical annotations in both)
compounds = plate_map.dropna(subset=["BROAD_CPD_ID"]).drop_duplicates("BROAD_CPD_ID")

garbled_dashes = ["\u201a\u00c4\u00ec", "\u221a\u00a2\u00ac\u00c4\u00ac\u00ec", "\u201a\u00c4\u00f6\u221a\u00a7\u221a\u00a8"]

pairs = (
    compounds.dropna(subset=["MOA"])
    .assign(moa=lambda x: x["MOA"].map(lambda text: split_moa(text, comma_in_name, placeholder, garbled_dashes)))
    .explode("moa")
    .rename(columns={"BROAD_CPD_ID": "compound_id"})
    .loc[lambda x: x["moa"] != "", ["moa", "compound_id"]]
)
print(f"compounds: {len(compounds)} | with a MOA: {pairs['compound_id'].nunique()} | compound-MOA pairs: {len(pairs)} | unique spellings: {pairs['moa'].nunique()}")


# In[4]:


# merge spellings that differ only in capitalization, keeping the most common one
spellings = pairs["moa"].value_counts()
canonical = spellings.groupby(spellings.index.str.lower()).idxmax()
merged = {s: canonical[s.lower()] for s in spellings.index if s != canonical[s.lower()]}
print(f"{len(merged)} spellings merged:")
print(pd.Series(merged).to_string())

annotation = (
    pairs.assign(moa=pairs["moa"].map(lambda s: canonical[s.lower()]))
    .drop_duplicates()
    .sort_values(["compound_id", "moa"])
    .reset_index(drop=True)
)
annotation.to_csv(output_file, index=False)

print(f"{len(annotation)} rows | compounds: {annotation['compound_id'].nunique()} | unique MOAs: {annotation['moa'].nunique()}")
annotation.head(10)

