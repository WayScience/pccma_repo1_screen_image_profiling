#!/usr/bin/env python
# coding: utf-8

# # Feature select single-cell profiles

# ## Import libraries

# In[ ]:


import csv
import os
import pathlib
import pprint
import time
from datetime import datetime, timezone

from pycytominer import feature_select


# ## Set paths and variables

# In[ ]:


# Directory containing the converted single-cell profile parquet per plate,
# produced by 0.convert_cytotable.ipynb. On Alpine (HPC), converted profiles
# are meant to live on the PetaLibrary "koala" mount, but fall back to the
# relative data/converted_profiles directory that 0.convert_cytotable.ipynb
# writes to in the repo checkout on scratch, in case a plate hasn't been
# synced to koala yet. Locally they live on the external drive. Only used
# here to enumerate which plates exist in the screen.
alpine_scratch_path = pathlib.Path("/scratch/alpine")

if alpine_scratch_path.exists():
    koala_dir = pathlib.Path("/pl/active/koala/ALSF_screen_data/SK-N-AS_repo1_profiles/converted_profiles")
    scratch_dir = pathlib.Path("data/converted_profiles")
    converted_dir = koala_dir if koala_dir.exists() else scratch_dir
else:
    converted_dir = pathlib.Path("/media/18tbdrive2/SK-N-AS_repo1_profiles/converted_profiles")

if not converted_dir.exists():
    raise FileNotFoundError(f"The converted profiles path {converted_dir} does not exist.")

# Directory containing per-plate normalized profiles from
# 4.single_cell_normalize.ipynb
normalized_dir = pathlib.Path("./data/single_cell_profiles")

# output path for feature-selected single-cell profiles (same directory as
# the normalized profiles they're derived from)
output_dir = normalized_dir

# extract the plate names from the converted profile file names
plate_names = sorted(
    file.stem.removesuffix("_converted")
    for file in converted_dir.glob("*_converted.parquet")
)

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

plate_names


# ## Set dictionary with plates to process

# In[ ]:


# Create plate info dictionary
plate_info_dictionary = {
    plate_id: {
        # Normalized profiles are produced per-plate by 4.single_cell_normalize.ipynb
        "normalized_path": (
            str((normalized_dir / f"{plate_id}_sc_normalized.parquet").resolve())
            if (normalized_dir / f"{plate_id}_sc_normalized.parquet").exists()
            else None
        ),
    }
    for plate_id in plate_names
}

# Display the dictionary to verify the entries
pprint.pprint(plate_info_dictionary, indent=4)


# In[ ]:


# Run configuration, both driven by environment variables so run_pipeline.sh
# can control them without papermill:
#   PLATE_ID  - restrict to a single plate (used to run each plate as its own
#               process so memory doesn't accumulate across plates). Leave
#               unset to process all plates, e.g. when running interactively.
#   OVERWRITE - if set (1/true/yes), reprocess and overwrite a plate's output
#               even if it already exists. Leave unset (the default) to skip
#               plates that already have output.
plate_id_filter = os.environ.get("PLATE_ID")
if plate_id_filter:
    if plate_id_filter not in plate_info_dictionary:
        raise ValueError(f"Unknown plate_id in PLATE_ID env var: {plate_id_filter}")
    plate_info_dictionary = {plate_id_filter: plate_info_dictionary[plate_id_filter]}

overwrite = os.environ.get("OVERWRITE", "").strip().lower() in ("1", "true", "yes")

plate_info_dictionary


# ## Feature select profiles

# In[ ]:


timing_log_path = output_dir / "timing_log_feature_select.csv"

for plate_id, info in plate_info_dictionary.items():
    if info["normalized_path"] is None:
        print(
            f"Skipping {plate_id}: not yet normalized "
            "(run 4.single_cell_normalize.ipynb for this plate first)"
        )
        continue

    output_feature_select_file = str(output_dir / f"{plate_id}_sc_feature_selected.parquet")

    # Already fully processed, so this plate can safely be skipped on a rerun
    # (unless OVERWRITE is set, in which case it's reprocessed regardless)
    if pathlib.Path(output_feature_select_file).exists() and not overwrite:
        print(f"Skipping {plate_id}: already processed (found {output_feature_select_file})")
        continue

    print(f"Feature selecting {plate_id}")
    plate_start_time = time.time()

    feature_select_start_time = time.time()
    feature_select(
        info["normalized_path"],
        operation=feature_select_ops,
        corr_threshold=0.90, # keep the same default value for correlation_threshold to be more strict as to identify the most informative features
        freq_cut=0.05, # keep the same default value for freq_cut
        unique_cut=0.01, # keep the same default value for unique_cut
        na_cutoff=0, # update na_cutoff from default to 0 to remove any columns with any NaN values (best for downstream modeling)
        output_file=output_feature_select_file,
        output_type="parquet",
    )
    feature_select_seconds = time.time() - feature_select_start_time

    total_seconds = time.time() - plate_start_time
    print(
        f"Feature selection completed for {plate_id}! "
        f"(feature_select={feature_select_seconds:.1f}s, total={total_seconds:.1f}s)"
    )

    # Record timing so per-plate performance can be compared across runs
    write_header = not timing_log_path.exists()
    with open(timing_log_path, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(
                ["plate_id", "feature_select_seconds", "total_seconds", "timestamp"]
            )
        writer.writerow(
            [
                plate_id,
                f"{feature_select_seconds:.1f}",
                f"{total_seconds:.1f}",
                datetime.now(timezone.utc).isoformat(timespec="seconds"),
            ]
        )

