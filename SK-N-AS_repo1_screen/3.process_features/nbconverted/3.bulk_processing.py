#!/usr/bin/env python
# coding: utf-8

# # Process bulk profiles
# 
# > NOTE: For bulk profiles, we normalize to the negative controls and use MAD robustize as the "standard".

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
import pyarrow.parquet as pq
from pycytominer import aggregate, feature_select, normalize
from pycytominer.cyto_utils import output


# ## Set paths and variables

# In[ ]:


# Directory containing the converted single-cell profile parquet per plate,
# produced by 0.convert_cytotable.ipynb. On Alpine (HPC), converted profiles
# are meant to live on the PetaLibrary "koala" mount, but fall back to the
# relative data/converted_profiles directory that 0.convert_cytotable.ipynb
# writes to in the repo checkout on scratch, in case a plate hasn't been
# synced to koala yet. Locally they live on the external drive since each
# plate's profile is tens of GB. Only used here for a cheap, metadata-only
# schema read (see per_cell_cols below) -- the actual single-cell data is
# read once, in 2.annotate.ipynb, not here.
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

# output path for bulk profiles
output_dir = pathlib.Path("./data/bulk_profiles")
output_dir.mkdir(parents=True, exist_ok=True)

# extract the plate names from the converted profile file names
plate_names = sorted(
    file.stem.removesuffix("_converted")
    for file in converted_dir.glob("*_converted.parquet")
)

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

plate_names


# ## Set dictionary with plates to process

# In[ ]:


# Create plate info dictionary
plate_info_dictionary = {
    plate_id: {
        # Raw converted profile, only used for a metadata-only schema read
        # (see per_cell_cols below)
        "profile_path": str(converted_dir / f"{plate_id}_converted.parquet"),

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


# ## Process data with pycytominer

# In[ ]:


timing_log_path = output_dir / "timing_log.csv"

for plate_id, info in plate_info_dictionary.items():
    if info["annotated_path"] is None:
        print(
            f"Skipping {plate_id}: not yet annotated "
            "(run 2.annotate.ipynb for this plate first)"
        )
        continue

    # Output file paths for each file
    output_aggregated_file = str(output_dir / f"{plate_id}_bulk_aggregated.parquet")
    output_normalized_file = str(output_dir / f"{plate_id}_bulk_normalized.parquet")
    output_feature_select_file = str(
        output_dir / f"{plate_id}_bulk_feature_selected.parquet"
    )

    # Already fully processed, so this plate can safely be skipped on a rerun
    # (unless OVERWRITE is set, in which case it's reprocessed regardless)
    if pathlib.Path(output_feature_select_file).exists() and not overwrite:
        print(f"Skipping {plate_id}: already processed (found {output_feature_select_file})")
        continue

    print(f"Now performing pycytominer pipeline for {plate_id}")
    plate_start_time = time.time()

    # Load the already-annotated single-cell profile (shared with
    # 4.single_cell_normalize.ipynb, so annotate() only ran once for this
    # plate, not once per branch)
    annotated_df = pd.read_parquet(info["annotated_path"])

    # Step 1: Aggregation
    # aggregate() only keeps `strata` columns plus the aggregated numeric
    # features, dropping everything else. Rather than re-annotating the
    # aggregated (well-level) data afterward to reattach platemap
    # metadata, include every platemap-derived metadata column that's
    # constant within a well directly in `strata`, so it survives
    # aggregation. per_cell_cols is read straight from the raw profile's
    # parquet schema (no data, just column names) so we can tell which
    # Metadata_* columns on annotated_df were added by annotate()'s
    # platemap join, without having to also load the raw single-cell file.
    aggregate_start_time = time.time()
    per_cell_cols = set(pq.ParquetFile(info["profile_path"]).schema_arrow.names)
    candidate_metadata_cols = [
        col for col in annotated_df.columns
        if col.startswith("Metadata_")
        and col not in ("Metadata_Plate", "Metadata_Well")
        and col not in per_cell_cols
    ]
    strata_cols = ["Metadata_Plate", "Metadata_Well"] + [
        col for col in candidate_metadata_cols
        if annotated_df.groupby("Metadata_Well")[col].nunique(dropna=False).le(1).all()
    ]
    aggregated_df = aggregate(
        population_df=annotated_df,
        operation="median",
        strata=strata_cols,
    )
    output(
        df=aggregated_df,
        output_filename=output_aggregated_file,
        output_type="parquet",
    )
    aggregate_seconds = time.time() - aggregate_start_time

    # Clear memory -- annotated_df is full single-cell scale; everything
    # from here on operates on the much smaller well-level aggregated_df
    del annotated_df
    gc.collect()

    # Step 2: Normalization (whole-plate) -- feed aggregate()'s own in-memory
    # result directly instead of re-reading output_aggregated_file back off disk.
    normalize_start_time = time.time()
    normalize(
        profiles=aggregated_df,
        method="mad_robustize", # use robustize to avoid influence of outliers in the normalization
        samples=neg_control_query, # normalize to negative controls only, not all wells (which would include compound wells)
        output_file=output_normalized_file,
        output_type="parquet",
    )
    normalize_seconds = time.time() - normalize_start_time

    # Clear memory
    del aggregated_df
    gc.collect()

    # Step 3: Feature selection
    feature_select_start_time = time.time()
    feature_select(
        output_normalized_file,
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
        f"Bulk processing completed for {plate_id}! "
        f"(aggregate={aggregate_seconds:.1f}s, "
        f"normalize={normalize_seconds:.1f}s, feature_select={feature_select_seconds:.1f}s, "
        f"total={total_seconds:.1f}s)"
    )

    # Record timing so per-plate performance can be compared across runs
    write_header = not timing_log_path.exists()
    with open(timing_log_path, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(
                [
                    "plate_id",
                    "aggregate_seconds",
                    "normalize_seconds",
                    "feature_select_seconds",
                    "total_seconds",
                    "timestamp",
                ]
            )
        writer.writerow(
            [
                plate_id,
                f"{aggregate_seconds:.1f}",
                f"{normalize_seconds:.1f}",
                f"{feature_select_seconds:.1f}",
                f"{total_seconds:.1f}",
                datetime.now(timezone.utc).isoformat(timespec="seconds"),
            ]
        )

