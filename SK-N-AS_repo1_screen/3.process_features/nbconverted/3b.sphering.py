#!/usr/bin/env python
# coding: utf-8

# # Sphere whole-screen bulk profiles
# 
# > Split out from 3.bulk_processing.ipynb into its own notebook/script so it can be invoked directly (e.g. by run_pipeline_hpc_sphering.sh) once every plate has been bulk-processed, instead of needing a SPHERING_ONLY flag to skip the per-plate loop inside that same script.
# 
# Adapted from Erik Serrano's work in the [fibrosis drug screen repository](https://github.com/WayScience/targeted_fibrosis_drug_screen/)
# 
# Unlike the source pipeline this step was adapted from, SK-N-AS has no batches / replicate plate groups to pool -- it's just one screen. So instead of looping per plate map and pooling that plate map's replicate plates, feature selection and the sphering fit/transform are all done once across every plate's pooled normalized profile, producing a single spherized profile for the whole screen (not one per plate).

# ## Import libraries

# In[ ]:


import os
import pathlib

import pandas as pd
from pycytominer import feature_select, normalize


# ## Set paths and variables

# In[2]:


# Directory containing per-plate bulk profiles from 3.bulk_processing.ipynb.
# Bulk processing runs on HPC, so its output isn't synced back into this repo
# checkout by default -- prefer this repo's local data/bulk_profiles dir if
# it's already been populated (e.g. copied down manually), and fall back to
# the bandicoot network mount otherwise.
local_bulk_dir = pathlib.Path("./data/bulk_profiles")
bandicoot_bulk_dir = pathlib.Path(
    "~/mnt/bandicoot/PCCMA_data/SK-N-AS_repo1_profiles/bulk_profiles"
).expanduser()

if any(local_bulk_dir.glob("*_bulk_normalized.parquet")):
    output_dir = local_bulk_dir
else:
    output_dir = bandicoot_bulk_dir

if not output_dir.exists():
    raise FileNotFoundError(f"The bulk profiles path {output_dir} does not exist.")

print(f"Reading bulk profiles from: {output_dir}")

# Wells with no compound (`Metadata_Batch_Id` is blank) that still received
# the DMSO vehicle are the negative controls
neg_control_query = 'Metadata_Batch_Id.isna() and Metadata_Solvent == "DMSO"'

# operations to perform for feature selection
# NOTE: drop_na_columns runs before correlation_threshold on purpose. Several
# Costes correlation features are NaN for some cells, and pycytominer's
# correlation_threshold falls back to a much slower NaN-aware pandas .corr()
# (instead of a fast BLAS np.corrcoef) if any NaNs are still present in the
# feature columns it's given. Dropping NaN columns first keeps it on the fast path.
feature_select_ops = [
    "drop_na_columns",
    "blocklist", # default block list uses same standard CP naming convention
    "frequency_threshold",
    "variance_threshold",
    "correlation_threshold",
]

spherized_output_dir = pathlib.Path("./data/spherized_profiles")
spherized_output_dir.mkdir(parents=True, exist_ok=True)

pooled_feature_select_file = output_dir / "SK-N-AS_repo1_screen_pooled_bulk_feature_selected.parquet"
output_spherized_file = spherized_output_dir / "SK-N-AS_repo1_screen_pooled_bulk_spherized.parquet"


# ## Run configuration

# In[3]:


# Run configuration, driven by an environment variable so run_pipeline.sh (or
# its HPC counterpart, run_pipeline_hpc_sphering.sh) can control it without
# papermill:
#   OVERWRITE - if set (1/true/yes), reprocess and overwrite the spherized
#               profile even if it already exists. Leave unset (the default)
#               to skip if it already exists.
overwrite = os.environ.get("OVERWRITE", "").strip().lower() in ("1", "true", "yes")


# ## Sphering step

# In[4]:


if output_spherized_file.exists() and not overwrite:
    print(f"Skipping sphering: already spherized (found {output_spherized_file})")
else:
    normalized_paths = sorted(output_dir.glob("*_bulk_normalized.parquet"))

    print(f"Pooling {len(normalized_paths)} plates for sphering...")

    # step 1: concat every plate's normalized profile before feature selection
    concat_df = pd.concat(
        [pd.read_parquet(path) for path in normalized_paths],
        ignore_index=True,
    ).reset_index(drop=True)

    # step 2a: Apply feature selection across the pooled screen to get a
    # common set of features for sphering.
    print("Feature selecting pooled screen...")
    feature_select_df = feature_select(
        profiles=concat_df,
        operation=feature_select_ops,
        na_cutoff=0, # updated from default to 0 to remove any columns with any NaN values (best for downstream modeling)
        corr_threshold=0.95, # increased from default to 0.95 to be more strict as to identify the most informative features
        freq_cut=0.05, # same as default
        output_file=pooled_feature_select_file,
        output_type="parquet",
    )

    # step 2b: Remove features with too little variation inside the exact
    # control population used to fit spherization.
    print(
        "Feature selecting pooled screen with frequency threshold within "
        "negative controls only..."
    )
    zero_negcon_var_fs_df = feature_select(
        profiles=feature_select_df,
        operation="frequency_threshold",
        freq_cut=0.05, # same as default
        unique_cut=0.01, # same as default
        samples=neg_control_query,
    )

    # step 3: Spherize/whiten the whole pooled screen using the pooled
    # negative controls as the reference population
    print("Sphering pooled screen using pooled negative controls...")
    normalize(
        profiles=zero_negcon_var_fs_df,
        method="spherize",
        samples=neg_control_query,
        spherize_center=True,
        spherize_method="ZCA-cor", # same as default
        spherize_epsilon=1e-6, #same as default
        output_file=output_spherized_file,
        output_type="parquet",
    )

    print(f"Saved feature-selected profiles to {pooled_feature_select_file}")
    print(f"Saved whole-screen spherized profile to {output_spherized_file}")

