#!/usr/bin/env python
# coding: utf-8

# # REPO1 target annotation
# 
# Builds `repo1_target_metadata_annotation.csv`: one row per compound and annotated target gene, with columns `gene_symbol`, `entrez_id` and `compound_id`.
# 
# The `Target` column of `REPO1_April2025.xlsx` holds comma-separated gene symbols. Entrez IDs are looked up with the MyGene.info service (only the gene symbols are sent). Symbols that are no longer current are matched through their aliases; symbols that cannot be matched (for example non-human targets) are kept with an empty `entrez_id`.

# In[ ]:


import pandas as pd
import mygene

plate_map_file = "../REPO1_April2025.xlsx"
output_file = "repo1_target_metadata_annotation.csv"


# In[ ]:


plate_map = pd.read_excel(plate_map_file, sheet_name="REPO2025_PlateMap")

# one row per compound (a compound plated in two wells has identical annotations in both)
compounds = plate_map.dropna(subset=["BROAD_CPD_ID"]).drop_duplicates("BROAD_CPD_ID")

# one row per compound and target gene
pairs = (
    compounds.dropna(subset=["Target"])
    .assign(gene_symbol=lambda x: x["Target"].str.split(","))
    .explode("gene_symbol")
    .assign(gene_symbol=lambda x: x["gene_symbol"].str.strip())
    .rename(columns={"BROAD_CPD_ID": "compound_id"})
    .loc[lambda x: x["gene_symbol"] != "", ["gene_symbol", "compound_id"]]
    .drop_duplicates()
)
symbols = sorted(pairs["gene_symbol"].unique())
print(f"compounds: {len(compounds)} | with a target: {pairs['compound_id'].nunique()} | compound-gene pairs: {len(pairs)} | unique symbols: {len(symbols)}")


# In[ ]:


mg = mygene.MyGeneInfo()

def query(terms, scope):
    hits = mg.querymany(terms, scopes=scope, fields="entrezgene,symbol", species="human", verbose=False)
    hits = pd.DataFrame(hits)
    # drop the terms MyGene.info could not find
    return hits[hits["notfound"].ne(True)] if "notfound" in hits else hits

# first pass: official gene symbols
official = query(symbols, "symbol")
official = official[official["symbol"].str.upper() == official["query"].str.upper()]
official = official.drop_duplicates("query")
mapping = dict(zip(official["query"], official["entrezgene"]))
print(f"matched by official symbol: {len(mapping)} of {len(symbols)}")


# In[ ]:


# second pass: symbols that are not (or no longer) official symbols, matched through aliases when unambiguous
remaining = [s for s in symbols if s not in mapping]
alias = query(remaining, "alias")
counts = alias.groupby("query")["entrezgene"].nunique()
unique_alias = alias[alias["query"].isin(counts[counts == 1].index)].drop_duplicates("query")
# family names that MyGene.info resolves to one member gene are not specific to that gene, so they are left empty
family_names = ["NACHR"]
unique_alias = unique_alias[~unique_alias["query"].isin(family_names)]
alias_mapping = dict(zip(unique_alias["query"], unique_alias["entrezgene"]))
mapping.update(alias_mapping)

ambiguous = sorted(counts[counts > 1].index)
unmatched = [s for s in symbols if s not in mapping]
print(f"matched by alias: {len(alias_mapping)} | ambiguous alias (left empty): {len(ambiguous)} | unmatched (left empty): {len(unmatched)}")
print("alias matches (symbol -> current symbol):")
print(unique_alias.set_index("query")["symbol"].to_string())
print("ambiguous:", ambiguous)
print("unmatched:", unmatched)


# In[ ]:


annotation = pairs.assign(entrez_id=pairs["gene_symbol"].map(mapping).astype("string")).loc[:, ["gene_symbol", "entrez_id", "compound_id"]]
annotation = annotation.sort_values(["compound_id", "gene_symbol"]).reset_index(drop=True)
annotation.to_csv(output_file, index=False)

print(f"{len(annotation)} rows | compounds: {annotation['compound_id'].nunique()} | genes: {annotation['gene_symbol'].nunique()} | rows without entrez_id: {annotation['entrez_id'].isna().sum()}")
annotation.head(10)

