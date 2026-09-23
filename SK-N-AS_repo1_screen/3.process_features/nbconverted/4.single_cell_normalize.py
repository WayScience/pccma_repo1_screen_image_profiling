#!/usr/bin/env python
# coding: utf-8

# # Normalize single-cell profiles
# 
# > NOTE: We normalize single-cells to the whole plate unlike bulk profiles, given that the sample size is much larger and less likely as impacted by variation.

# ## Import libraries

# In[ ]:


import csv
import gc
import os
import pathlib
import pprint
import time
from datetime import datetime, timezone

import pandas as pd
from pycytominer import normalize


# ## Set paths and variables

# In[ ]:


# Directory containing the converted single-cell profile parquet per plate,
# produced by 0.convert_cytotable.ipynb. On Alpine (HPC), converted profiles
# are meant to live on the PetaLibrary "koala" mount, but fall back to the
# relative data/converted_profiles directory that 0.convert_cytotable.ipynb
# writes to in the repo checkout on scratch, in case a plate hasn't been
# synced to koala yet. Locally they live on the external drive. Only used
# here to enumerate which plates exist in the screen -- the actual
# single-cell data is read from the already-annotated file below, not from
# here.
alpine_scratch_path = pathlib.Path("/scratch/alpine")

if alpine_scratch_path.exists():
    koala_dir = pathlib.Path("/pl/active/koala/ALSF_screen_data/SK-N-AS_repo1_profiles/converted_profiles")
    scratch_dir = pathlib.Path("data/converted_profiles")
    converted_dir = koala_dir if koala_dir.exists() else scratch_dir
else:
    converted_dir = pathlib.Path("/media/18tbdrive2/SK-N-AS_repo1_profiles/converted_profiles")

if not converted_dir.exists():
    raise FileNotFoundError(f"The converted profiles path {converted_dir} does not exist.")

# Directory containing per-plate annotated profiles from 2.annotate.ipynb
annotated_dir = pathlib.Path("./data/annotated_profiles")

# output path for normalized single-cell profiles
output_dir = pathlib.Path("./data/single_cell_profiles")
output_dir.mkdir(parents=True, exist_ok=True)

# extract the plate names from the converted profile file names
plate_names = sorted(
    file.stem.removesuffix("_converted")
    for file in converted_dir.glob("*_converted.parquet")
)

plate_names


# ## Set dictionary with plates to process

# In[ ]:


# Create plate info dictionary
plate_info_dictionary = {
    plate_id: {
        # Annotated profiles are produced per-plate by 2.annotate.ipynb
        "annotated_path": (
            str((annotated_dir / f"{plate_id}_annotated.parquet").resolve())
            if (annotated_dir / f"{plate_id}_annotated.parquet").exists()
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


# ## Normalize profiles

# In[ ]:


timing_log_path = output_dir / "timing_log_normalize.csv"

for plate_id, info in plate_info_dictionary.items():
    if info["annotated_path"] is None:
        print(
            f"Skipping {plate_id}: not yet annotated "
            "(run 2.annotate.ipynb for this plate first)"
        )
        continue

    output_normalized_file = str(output_dir / f"{plate_id}_sc_normalized.parquet")

    # Already normalized, so this plate can safely be skipped on a rerun
    # (unless OVERWRITE is set, in which case it's reprocessed regardless)
    if pathlib.Path(output_normalized_file).exists() and not overwrite:
        print(f"Skipping {plate_id}: already normalized (found {output_normalized_file})")
        continue

    print(f"Normalizing {plate_id}")
    plate_start_time = time.time()

    # Load the already-annotated single-cell profile (shared with
    # 3.bulk_processing.ipynb, so annotate() only ran once for this plate)
    annotated_df = pd.read_parquet(info["annotated_path"])

    normalize_start_time = time.time()
    normalize(
        profiles=annotated_df,
        method="standardize", # use standardize method as default
        output_file=output_normalized_file,
        output_type="parquet",
        samples="all", # apply normalization based on all samples
    )
    normalize_seconds = time.time() - normalize_start_time

    # Clear memory
    del annotated_df
    gc.collect()

    total_seconds = time.time() - plate_start_time
    print(
        f"Normalization completed for {plate_id}! "
        f"(normalize={normalize_seconds:.1f}s, total={total_seconds:.1f}s)"
    )

    # Record timing so per-plate performance can be compared across runs
    write_header = not timing_log_path.exists()
    with open(timing_log_path, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(
                ["plate_id", "normalize_seconds", "total_seconds", "timestamp"]
            )
        writer.writerow(
            [
                plate_id,
                f"{normalize_seconds:.1f}",
                f"{total_seconds:.1f}",
                datetime.now(timezone.utc).isoformat(timespec="seconds"),
            ]
        )

