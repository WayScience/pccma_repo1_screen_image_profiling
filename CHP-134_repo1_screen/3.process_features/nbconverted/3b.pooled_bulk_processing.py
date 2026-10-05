#!/usr/bin/env python
# coding: utf-8

# # Pooled bulk processing, CHP-134
# 
# Processes the whole screen's well-level profiles the way the [JUMP profiling recipe](https://github.com/broadinstitute/jump-profiling-recipe) processes its compound screens, using pycytominer wherever it provides the step. The per-plate, un-normalized profiles from `3.bulk_processing.ipynb` are pooled and taken through:
# 
# | Step | What | Tool |
# |---|---|---|
# | 1 | Pool the plates into one table | `utils/bulk_utils.py` |
# | 2 | Drop every feature with a NaN or infinite value in any well | `utils/bulk_utils.py` |
# | 3 | Position correction: subtract each feature's plate-column gradient, estimated from the compound wells | `utils/bulk_utils.py` |
# | 4 | Variant filter: keep a feature only if, on every plate, the control wells have MAD ≠ 0 and \|MAD / median\| > 1e-3 | `utils/bulk_utils.py` |
# | 5 | MAD-normalize each plate to its own control wells | pycytominer `normalize` |
# | 6 | Rank-based inverse normal transform of every feature over all wells (quantile method) | pycytominer `normalize` |
# | 7 | Pooled feature selection: frequency filter, correlation filter (0.9), blocklist, drop NaN | pycytominer `feature_select` |
# 
# The result is written to `data/bulk_profiles/`; `3c.sphering.ipynb` then spheres it. Differences from the JUMP compound recipe: it has no step like our position correction (step 3); the inverse normal transform is pycytominer's quantile method rather than the exact-rank (Blom) transform; and there is no Harmony step.

# ## Import libraries

# In[ ]:


import os
import pathlib
import sys

from pycytominer.cyto_utils.features import infer_cp_features

sys.path.append("../../utils")
import bulk_utils


# ## Set paths and variables

# In[ ]:


screen_name = "CHP-134_repo1_screen"

# Per-plate un-normalized well-level profiles from 3.bulk_processing.ipynb. Use this screen's local
# data/bulk_profiles folder if it has them, and fall back to the bandicoot network mount otherwise
# (BULK_PROFILES_DIR overrides both). The CHP-134 files on bandicoot keep their older *_bulk_annotated name.
local_bulk_dir = pathlib.Path("./data/bulk_profiles")
per_plate_glob = "*_bulk_aggregated.parquet"
bandicoot_bulk_dir = pathlib.Path("~/mnt/bandicoot/PCCMA_data/CHP-134_repo1_screen_outputs/profiles/bulk_profiles").expanduser()
bandicoot_glob = "*_bulk_annotated.parquet"

output_file = local_bulk_dir / f"{screen_name}_pooled_bulk_feature_selected.parquet"

plate_column = "Metadata_Plate"
well_column = "Metadata_Well"

# Empty wells (no compound and no solvent, so no treatment at all) are dropped when pooling, as the JUMP recipe
# drops its untreated wells when it loads the data
untreated_query = 'Metadata_Batch_Id.isna() and Metadata_Solvent.isna()'

# Wells with no compound (`Metadata_Batch_Id` is blank) that still received the DMSO vehicle are the negative controls
neg_control_query = 'Metadata_Batch_Id.isna() and Metadata_Solvent == "DMSO"'

# The plate-map grouping used by the position correction (step 3): plates whose compounds overlap by at
# least this fraction (Jaccard) are treated as the same layout
min_jaccard = 0.9

# The JUMP recipe's values
min_abs_coef_var = 1e-3
correlation_threshold = 0.9


# In[ ]:


if os.environ.get("BULK_PROFILES_DIR"):
    plate_files = sorted(pathlib.Path(os.environ["BULK_PROFILES_DIR"]).expanduser().glob(per_plate_glob))
elif any(local_bulk_dir.glob(per_plate_glob)):
    plate_files = sorted(local_bulk_dir.glob(per_plate_glob))
else:
    plate_files = sorted(bandicoot_bulk_dir.glob(bandicoot_glob))
if not plate_files:
    raise FileNotFoundError("no per-plate profile files found")
print(f"Reading {len(plate_files)} per-plate profiles from {plate_files[0].parent}")


def report(step, profiles, features):
    n_controls = int(bulk_utils.control_mask(profiles, neg_control_query).sum())
    print(f"{step}: {len(profiles)} wells ({n_controls} controls) on {profiles[plate_column].nunique()} plates, {len(features)} features")


# ## Step 1: pool the plates
# 
# One row per well, all plates stacked, with the values as aggregated (a per-well median of the single cells). The empty wells (no compound and no solvent) are dropped here, as the JUMP recipe drops its untreated wells when loading the data; they stay in the per-plate files for diagnostics.

# In[ ]:


pooled = bulk_utils.pool_plate_profiles(plate_files, plate_column)
pooled, n_untreated = bulk_utils.drop_wells(pooled, untreated_query)
print(f"dropped {n_untreated} empty wells (no compound, no solvent)")
assert not pooled.duplicated([plate_column, "Metadata_Well"]).any(), "a plate/well appears more than once"

features = infer_cp_features(pooled)
report("1 pooled", pooled, features)


# ## Step 2: drop features with NaN or infinite values
# 
# A feature is dropped if any well on any plate has a NaN or an infinite value. This includes features that some plates lack altogether (they are NaN on those plates).

# In[ ]:


pooled, features, nonfinite = bulk_utils.drop_nonfinite_features(pooled, features)
print(f"dropped {len(nonfinite)} features with NaN/inf values; examples: {nonfinite[:5]}")
report("2 finite features", pooled, features)


# ## Step 3: position correction
# 
# Every plate has the same left-to-right gradient: DMSO controls differ between plate columns 2 and 24 (where they sit), and compound wells drift smoothly between the two across the columns in between (confirmed on these un-normalized profiles: plate column alone explains up to 9% of a feature's variance on average, over compound wells, and up to ~40% for the worst-affected features). The same axis and about the same size of gap show up in the DMSO controls, which have no compound at all, so this looks like a plate-wide technical effect rather than biology.
# 
# For every feature and plate column, the correction is the median, across the plate-map layouts (groups of plates carrying the same compounds), of that column's mean value among the compound wells of that layout — using the median keeps one layout from dominating the curve. It is subtracted from every well, including the controls, so the per-plate baseline the next step computes is no longer stretched by the two control columns sitting at opposite ends of the gradient. Columns with no compound wells (2 and 24, where only controls sit) reuse the nearest compound column's value.
# 
# In a stand-alone test on these profiles, this correction (versus none): closed the gap between the DMSO column-2 and column-24 medians by about 40%, raised the features passing the variant filter, and raised replicate recovery (the share of a well's 10 nearest neighbours that are the same compound on another plate) from 4.1% to 5.3% (SK-N-AS) and 2.6% to 3.5% (CHP-134).

# In[ ]:


pooled = bulk_utils.correct_column_gradient(pooled, features, plate_column, well_column=well_column, min_jaccard=min_jaccard)
report("3 position corrected", pooled, features)


# ## Step 4: variant filter
# 
# Each plate's median and MAD are computed from its own control wells, on the values as aggregated. A feature is kept only if, on **every** plate, its MAD is not 0 and |MAD / median| is above `min_abs_coef_var`. Features with a MAD of 0 in some plate's controls would otherwise be divided by ~0 in the next step and end up with values around 1e17.

# In[ ]:


is_control = bulk_utils.control_mask(pooled, neg_control_query)
stats = bulk_utils.plate_control_stats(pooled, features, is_control, plate_column)

features_kept, features_dropped = bulk_utils.variant_features(stats, min_abs_coef_var)
print(f"dropped {len(features_dropped)} of {len(features)} features that do not vary in the controls of every plate; examples: {features_dropped[:5]}")

pooled = pooled.drop(columns=features_dropped)
features = features_kept
report("4 variant features", pooled, features)


# ## Step 5: MAD normalization per plate
# 
# pycytominer's `normalize(method="mad_robustize")`, fitted on each plate's own control wells and applied to that plate.

# In[ ]:


pooled = bulk_utils.normalize_by_plate(pooled, features, neg_control_query, plate_column)
assert pooled[features].notna().all().all(), "NaN values after MAD normalization"
report("5 mad normalized", pooled, features)


# ## Step 6: inverse normal transform
# 
# pycytominer's `normalize(method="inverse_normal")` over all wells of the screen: each feature is rank-transformed to a standard normal (interpolating between 1,000 quantile landmarks).

# In[ ]:


pooled = bulk_utils.inverse_normal_transform(pooled, features)
report("6 inverse normal", pooled, features)


# ## Step 7: feature selection
# 
# pycytominer's `feature_select`, on all wells together, in the JUMP recipe's order: frequency filter, correlation filter, blocklist, drop features with NaN.

# In[ ]:


pooled = bulk_utils.select_features(pooled, features, correlation_threshold)
features = infer_cp_features(pooled)
report("7 feature selected", pooled, features)

output_file.parent.mkdir(parents=True, exist_ok=True)
pooled.to_parquet(output_file, index=False)
print(f"saved {output_file}")

