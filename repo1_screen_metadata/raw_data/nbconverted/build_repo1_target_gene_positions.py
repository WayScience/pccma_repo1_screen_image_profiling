#!/usr/bin/env python
# coding: utf-8

# # REPO1 target gene positions
# 
# Builds `repo1_target_gene_positions.csv`: the genomic position (human genome build GRCh38/hg38) of every target gene in `repo1_target_metadata_annotation.csv`, with columns `entrez_id`, `gene_symbol`, `chromosome`, `start` and `end`.
# 
# Positions are looked up with the MyGene.info service (only the Entrez IDs are sent). Genes are kept only when they have a position on a primary chromosome (1-22, X, Y); genes on other sequences (mitochondrial, unplaced or alternate contigs) or with no position are listed and left out.

# In[ ]:


import pandas as pd
import mygene

annotation_file = "repo1_target_metadata_annotation.csv"
output_file = "repo1_target_gene_positions.csv"
primary_chromosomes = [str(i) for i in range(1, 23)] + ["X", "Y"]


# In[ ]:


annotation = pd.read_csv(annotation_file, dtype={"entrez_id": "string"})
genes = annotation.dropna(subset=["entrez_id"]).drop_duplicates("entrez_id")[["entrez_id", "gene_symbol"]]
print(f"genes with an Entrez ID: {len(genes)}")

mg = mygene.MyGeneInfo()
hits = mg.querymany(list(genes["entrez_id"]), scopes="entrezgene", fields="genomic_pos,symbol", species="human", verbose=False)

def primary_position(hit):
    """The position of a hit on a primary chromosome (a gene can have several positions, e.g. on alternate contigs)."""
    positions = hit.get("genomic_pos", [])
    positions = [positions] if isinstance(positions, dict) else positions
    for pos in positions:
        if str(pos["chr"]) in primary_chromosomes:
            return {"chromosome": str(pos["chr"]), "start": pos["start"], "end": pos["end"]}
    return None

rows = []
for hit in hits:
    position = primary_position(hit)
    if position is not None:
        rows.append({"entrez_id": hit["query"], **position})

positions = genes.merge(pd.DataFrame(rows), on="entrez_id", how="left")
missing = positions[positions["chromosome"].isna()]
print(f"genes with a primary chromosome position: {positions['chromosome'].notna().sum()} | without: {len(missing)}")
print("without a position:", ", ".join(missing["gene_symbol"]))


# In[ ]:


positions = positions.dropna(subset=["chromosome"]).astype({"start": int, "end": int})
positions["chromosome"] = pd.Categorical(positions["chromosome"], primary_chromosomes, ordered=True)
positions = positions.sort_values(["chromosome", "start"]).reset_index(drop=True)
positions.to_csv(output_file, index=False)

print(positions.groupby("chromosome", observed=True).size().to_string())
positions.head()

