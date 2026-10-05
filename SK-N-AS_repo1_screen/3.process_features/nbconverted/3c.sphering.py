#!/usr/bin/env python
# coding: utf-8

# # Sphere the pooled bulk profiles, SK-N-AS
# 
# Fits a ZCA-cor whitening transform on the pooled negative-control wells and applies it to every well of the screen, as done before, but on the profiles from `3b.pooled_bulk_processing.ipynb`. Features with too little variation among the controls are removed first, since they cannot be whitened.

# ## Import libraries

# In[ ]:


import pathlib
import sys

import pandas as pd
from pycytominer.cyto_utils.features import infer_cp_features

sys.path.append("../../utils")
import bulk_utils


# ## Set paths and variables

# In[ ]:


screen_name = "SK-N-AS_repo1_screen"

input_file = pathlib.Path("./data/bulk_profiles") / f"{screen_name}_pooled_bulk_feature_selected.parquet"

spherized_dir = pathlib.Path("./data/spherized_profiles")
spherized_dir.mkdir(parents=True, exist_ok=True)
output_file = spherized_dir / f"{screen_name}_pooled_bulk_spherized.parquet"

# Wells with no compound (`Metadata_Batch_Id` is blank) that still received the DMSO vehicle are the negative controls
neg_control_query = 'Metadata_Batch_Id.isna() and Metadata_Solvent == "DMSO"'

sphering_epsilon = 1e-6


# ## Sphering

# In[ ]:


if not input_file.exists():
    raise FileNotFoundError(f"{input_file} not found (run 3b.pooled_bulk_processing.ipynb first)")
pooled = pd.read_parquet(input_file)
features = infer_cp_features(pooled)

# Fewer control wells than features leaves part of the feature space unconstrained by the controls
n_controls = int(bulk_utils.control_mask(pooled, neg_control_query).sum())
print(f"{len(pooled)} wells, {n_controls} control wells, {len(features)} features")

spherized = bulk_utils.spherize(pooled, neg_control_query, sphering_epsilon)
print(f"{len(infer_cp_features(spherized))} features after sphering")

spherized.to_parquet(output_file, index=False)
print(f"saved {output_file}")

