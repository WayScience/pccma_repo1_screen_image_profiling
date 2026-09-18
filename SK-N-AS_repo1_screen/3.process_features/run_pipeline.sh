#!/bin/bash
#
# Run the full feature-processing pipeline, in order:
#   1.single_cell_qc.ipynb              -> one QC-annotations parquet per plate (via papermill,
#                                          which keeps the executed notebook per plate for review)
#   2.annotate.ipynb                    -> one annotated single-cell profile per plate, shared
#                                          input for both branches below (via nbconverted script,
#                                          one plate/process at a time)
#   3.bulk_processing.ipynb             -> aggregated/normalized/feature-selected bulk profiles
#                                          (via nbconverted script, one plate/process at a time)
#   3b.sphering.ipynb                   -> whole-screen sphering step, once every plate has been
#                                          bulk-processed (via nbconverted script, a single process
#                                          -- split out of 3.bulk_processing.ipynb so it can be
#                                          invoked directly instead of via a SPHERING_ONLY flag)
#   4.single_cell_normalize.ipynb       -> normalized single-cell profiles (via nbconverted
#                                          script, one plate/process at a time)
#   5.single_cell_feature_select.ipynb  -> feature-selected single-cell profiles (via nbconverted
#                                          script, one plate/process at a time)
#
# Unlike CHP-134_repo1_screen's run_pipeline.sh, there is no merge_profiles
# step here: 0.convert_cytotable.ipynb (run separately, via
# convert_cytotable_hpc_parent.sh/convert_cytotable_hpc_child.sh) already
# writes one converted parquet per plate directly to the external drive,
# with no row-batch chunks to stitch together first.
#
# Annotation used to be duplicated -- once inside bulk processing, once inside
# single-cell processing -- so it ran (and held both a plate's raw profile and
# its annotated copy in memory at once, each tens of GB) twice per plate. It's
# now its own step (2.annotate.ipynb) that both downstream branches read from,
# so that memory-heavy merge only happens once per plate.
#
# Steps 2, 3, 4, and 5 all run one plate per `python` invocation (PLATE_ID env
# var) rather than looping over all plates inside a single long-lived process.
# Each plate's memory is fully released back to the OS when its process
# exits, which avoids the cross-plate memory buildup that can otherwise lead
# to an OOM kill partway through a batch. This matters even more here than
# for CHP-134, since each SK-N-AS plate (~55GB, 2.17M cells) is roughly 3x the
# size of a CHP-134 plate -- a single plate's raw profile alone has been
# measured at ~97GB in memory, so isolating each stage to its own process is
# load-bearing, not just a nice-to-have.
#
# Each plate's per-stage timing (and a total) is appended to a timing_log.csv
# under that step's output directory (data/annotated_profiles,
# data/bulk_profiles, or data/single_cell_profiles -- the latter split into
# timing_log_normalize.csv and timing_log_feature_select.csv since steps 4
# and 5 run as separate processes).
#
# By default, every step skips a plate that already has output. Set
# OVERWRITE=1 (e.g. `OVERWRITE=1 ./run_pipeline.sh`) to reprocess and
# overwrite every plate's output instead.

set -uo pipefail

overwrite="${OVERWRITE:-}"

# Converted single-cell profiles live on the external drive, one parquet per
# plate, produced by 0.convert_cytotable.ipynb -- not under ./data like
# CHP-134's merged_profiles.
converted_dir="/media/18tbdrive2/SK-N-AS_repo1_profiles/converted_profiles"

# -----------------------------
# Total runtime tracking (printed on exit, success or failure, via trap)
# -----------------------------
pipeline_start_time=$(date +%s)

print_total_time() {
    local elapsed=$(( $(date +%s) - pipeline_start_time ))
    echo "======================================"
    printf 'Total pipeline time: %dh %dm %ds\n' \
        $((elapsed / 3600)) $(((elapsed % 3600) / 60)) $((elapsed % 60))
    echo "======================================"
}
trap print_total_time EXIT

# -----------------------------
# Initialize environment
# -----------------------------
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate pccma_repo1_preprocessing_env

if ! command -v papermill >/dev/null 2>&1; then
    echo "papermill not found on PATH after activating pccma_repo1_preprocessing_env; aborting."
    exit 1
fi

if ! python -c "import pandas, pyarrow, pycytominer" >/dev/null 2>&1; then
    echo "Active Python ($(command -v python)) is missing required packages (pandas/pyarrow/pycytominer); aborting."
    exit 1
fi

if [ ! -d "$converted_dir" ]; then
    echo "Converted profiles directory not found: $converted_dir (external drive mounted?); aborting."
    exit 1
fi

failed_plates=()

# convert notebooks to scripts
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb
if [ $? -ne 0 ]; then
    echo "jupyter nbconvert FAILED; aborting pipeline before QC."
    exit 1
fi

mapfile -t plate_ids < <(find "$converted_dir" -maxdepth 1 -name "*_converted.parquet" -printf "%f\n" | sed 's/_converted\.parquet$//' | sort)

echo "Number of plates found: ${#plate_ids[@]}"

# -----------------------------
# Step 1: per-plate single-cell QC via papermill
# -----------------------------
echo "======================================"
echo "Step 1: 1.single_cell_qc.ipynb (papermill, per plate)"
echo "======================================"

qc_notebook_dir="./executed_notebooks/1.single_cell_qc"
mkdir -p "$qc_notebook_dir"

for plate_id in "${plate_ids[@]}"; do
    qc_output="./data/qc_results/${plate_id}_qc_annotations.parquet"
    executed_notebook="${qc_notebook_dir}/${plate_id}.ipynb"

    if [ -f "$qc_output" ] && [ -z "$overwrite" ]; then
        echo "✅ ${plate_id} already QC'd (found ${qc_output})"
        continue
    fi

    echo ">>> Running QC for ${plate_id}"
    # render_diagnostics=False skips the per-condition CytoDataFrame image
    # previews, which pull crops off the bandicoot network mount and are only
    # useful for interactive review -- the QC annotations export at the
    # bottom of the notebook doesn't depend on them. Leaving this on for a
    # batch run across many plates is what causes each plate to take a very
    # long time (or effectively hang) here.
    papermill 1.single_cell_qc.ipynb "$executed_notebook" \
        -p plate_id "$plate_id" \
        -p render_diagnostics False
    exit_code=$?

    if [ "$exit_code" -ne 0 ]; then
        echo "QC FAILED for plate: ${plate_id} (exit code $exit_code)"
        failed_plates+=("qc:${plate_id}")
    else
        echo "QC done for plate: ${plate_id}"
    fi
done

# -----------------------------
# Step 2: annotation (one plate/process at a time; skips plates without QC,
# and already-annotated plates). Shared input for both the bulk and
# single-cell branches below, so annotate() only runs once per plate.
# -----------------------------
echo "======================================"
echo "Step 2: 2.annotate.ipynb (one plate at a time)"
echo "======================================"

for plate_id in "${plate_ids[@]}"; do
    annotate_output="./data/annotated_profiles/${plate_id}_annotated.parquet"

    if [ -f "$annotate_output" ] && [ -z "$overwrite" ]; then
        echo "✅ ${plate_id} already annotated (found ${annotate_output})"
        continue
    fi

    echo ">>> Running annotation for ${plate_id}"
    PLATE_ID="$plate_id" OVERWRITE="$overwrite" python nbconverted/2.annotate.py
    exit_code=$?

    if [ "$exit_code" -ne 0 ]; then
        echo "Annotation FAILED for plate: ${plate_id} (exit code $exit_code)"
        failed_plates+=("annotate:${plate_id}")
    else
        echo "Annotation done for plate: ${plate_id}"
    fi
done

# -----------------------------
# Step 3: bulk processing (one plate/process at a time; skips plates without
# an annotated profile, and already-processed plates)
# -----------------------------
echo "======================================"
echo "Step 3: 3.bulk_processing.ipynb (one plate at a time)"
echo "======================================"

for plate_id in "${plate_ids[@]}"; do
    bulk_output="./data/bulk_profiles/${plate_id}_bulk_feature_selected.parquet"

    if [ -f "$bulk_output" ] && [ -z "$overwrite" ]; then
        echo "✅ ${plate_id} already bulk-processed (found ${bulk_output})"
        continue
    fi

    echo ">>> Running bulk processing for ${plate_id}"
    PLATE_ID="$plate_id" OVERWRITE="$overwrite" python nbconverted/3.bulk_processing.py
    exit_code=$?

    if [ "$exit_code" -ne 0 ]; then
        echo "Bulk processing FAILED for plate: ${plate_id} (exit code $exit_code)"
        failed_plates+=("bulk_processing:${plate_id}")
    else
        echo "Bulk processing done for plate: ${plate_id}"
    fi
done

# -----------------------------
# Step 3b: sphering (fits the whitening transform once on every plate's
# pooled normalized profile -- there are no batches/replicate plate groups in
# this screen -- then applies it once across the whole pooled screen, writing
# a single spherized profile for the whole screen. Its own notebook/script
# (3b.sphering.ipynb), split out of 3.bulk_processing.ipynb, so it can be
# invoked directly here rather than needing a SPHERING_ONLY flag to skip that
# script's per-plate loop.
# -----------------------------
echo "======================================"
echo "Step 3b: 3b.sphering.ipynb (whole screen)"
echo "======================================"

spherized_file="./data/spherized_profiles/SK-N-AS_repo1_screen_pooled_bulk_spherized.parquet"

if [ -f "$spherized_file" ] && [ -z "$overwrite" ]; then
    echo "✅ Whole-screen spherized profile already exists (${spherized_file})"
else
    echo ">>> Running sphering for the whole screen"
    OVERWRITE="$overwrite" python nbconverted/3b.sphering.py
    exit_code=$?

    if [ "$exit_code" -ne 0 ]; then
        echo "Sphering FAILED (exit code $exit_code)"
        failed_plates+=("sphering:whole_screen")
    else
        echo "Sphering done for whole screen"
    fi
fi

# -----------------------------
# Step 4: single-cell normalization (one plate/process at a time; skips
# plates without an annotated profile, and already-normalized plates)
# -----------------------------
echo "======================================"
echo "Step 4: 4.single_cell_normalize.ipynb (one plate at a time)"
echo "======================================"

for plate_id in "${plate_ids[@]}"; do
    sc_normalized_output="./data/single_cell_profiles/${plate_id}_sc_normalized.parquet"

    if [ -f "$sc_normalized_output" ] && [ -z "$overwrite" ]; then
        echo "✅ ${plate_id} already normalized (found ${sc_normalized_output})"
        continue
    fi

    echo ">>> Running single-cell normalization for ${plate_id}"
    PLATE_ID="$plate_id" OVERWRITE="$overwrite" python nbconverted/4.single_cell_normalize.py
    exit_code=$?

    if [ "$exit_code" -ne 0 ]; then
        echo "Single-cell normalization FAILED for plate: ${plate_id} (exit code $exit_code)"
        failed_plates+=("single_cell_normalize:${plate_id}")
    else
        echo "Single-cell normalization done for plate: ${plate_id}"
    fi
done

# -----------------------------
# Step 5: single-cell feature selection (one plate/process at a time; skips
# plates without a normalized profile, and already-processed plates)
# -----------------------------
echo "======================================"
echo "Step 5: 5.single_cell_feature_select.ipynb (one plate at a time)"
echo "======================================"

for plate_id in "${plate_ids[@]}"; do
    sc_output="./data/single_cell_profiles/${plate_id}_sc_feature_selected.parquet"

    if [ -f "$sc_output" ] && [ -z "$overwrite" ]; then
        echo "✅ ${plate_id} already processed (found ${sc_output})"
        continue
    fi

    echo ">>> Running single-cell feature selection for ${plate_id}"
    PLATE_ID="$plate_id" OVERWRITE="$overwrite" python nbconverted/5.single_cell_feature_select.py
    exit_code=$?

    if [ "$exit_code" -ne 0 ]; then
        echo "Single-cell feature selection FAILED for plate: ${plate_id} (exit code $exit_code)"
        failed_plates+=("single_cell_feature_select:${plate_id}")
    else
        echo "Single-cell feature selection done for plate: ${plate_id}"
    fi
done

conda deactivate

echo "======================================"
echo "Pipeline finished."
if [ ${#failed_plates[@]} -gt 0 ]; then
    echo "Failures:"
    printf '  %s\n' "${failed_plates[@]}"
    exit 1
fi

echo "✅ All steps completed successfully."
