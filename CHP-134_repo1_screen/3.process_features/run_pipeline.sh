#!/bin/bash
#
# Run the full feature-processing pipeline, in order:
#   1.merge_profiles.ipynb       -> one merged parquet per plate (via nbconverted script)
#   2.single_cell_qc.ipynb       -> one QC-annotations parquet per plate (via papermill,
#                                    which keeps the executed notebook per plate for review)
#   3.bulk_processing.ipynb      -> aggregated/annotated/normalized/feature-selected bulk profiles
#                                    (via nbconverted script, one plate/process at a time)
#   4.single_cell_processing.ipynb -> annotated/normalized/feature-selected single-cell profiles
#                                    (via nbconverted script, one plate/process at a time)
#
# Steps 3 and 4 run one plate per `python` invocation (PLATE_ID env var) rather
# than looping over all plates inside a single long-lived process. Each plate's
# memory is fully released back to the OS when its process exits, which avoids
# the cross-plate memory buildup that can otherwise lead to an OOM kill partway
# through a batch.
#
# Each plate's per-stage timing (and a total) is appended to
# data/bulk_profiles/timing_log.csv and data/single_cell_profiles/timing_log.csv.
#
# By default, every step skips a plate that already has output. Set
# OVERWRITE=1 (e.g. `OVERWRITE=1 ./run_pipeline.sh`) to reprocess and
# overwrite every plate's output instead.

set -uo pipefail

overwrite="${OVERWRITE:-}"

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
if [ $? -ne 0 ]; then
    echo "Failed to activate conda env pccma_repo1_preprocessing_env; aborting."
    exit 1
fi

failed_plates=()

# convert notebooks to scripts
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb
if [ $? -ne 0 ]; then
    echo "jupyter nbconvert FAILED; aborting pipeline before QC."
    exit 1
fi

# -----------------------------
# Step 1: merge row-batch profiles into one profile per plate
# -----------------------------
echo "======================================"
echo "Step 1: 1.merge_profiles.ipynb"
echo "======================================"
OVERWRITE="$overwrite" python nbconverted/1.merge_profiles.py
if [ $? -ne 0 ]; then
    echo "1.merge_profiles.py FAILED (incomplete merge); aborting pipeline before QC."
    exit 1
fi

# -----------------------------
# Step 2: per-plate single-cell QC via papermill
# -----------------------------
echo "======================================"
echo "Step 2: 2.single_cell_qc.ipynb (papermill, per plate)"
echo "======================================"

qc_notebook_dir="./executed_notebooks/2.single_cell_qc"
mkdir -p "$qc_notebook_dir"

mapfile -t plate_ids < <(find ./data/merged_profiles -maxdepth 1 -name "*.parquet" -printf "%f\n" | sed 's/\.parquet$//' | sort)

echo "Number of plates found: ${#plate_ids[@]}"

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
    papermill 2.single_cell_qc.ipynb "$executed_notebook" \
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
# Step 3: bulk processing (one plate/process at a time; skips plates without
# QC, and already-processed plates)
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
# a single spherized profile for the whole screen. Run once with PLATE_ID
# unset, rather than per plate like the loop above.)
# -----------------------------
echo "======================================"
echo "Step 3b: 3.bulk_processing.ipynb (sphering, whole screen)"
echo "======================================"

spherized_file="./data/spherized_profiles/repo1_screen_bulk_spherized.parquet"

if [ -f "$spherized_file" ] && [ -z "$overwrite" ]; then
    echo "✅ Whole-screen spherized profile already exists (${spherized_file})"
else
    echo ">>> Running sphering for the whole screen"
    OVERWRITE="$overwrite" python nbconverted/3.bulk_processing.py
    exit_code=$?

    if [ "$exit_code" -ne 0 ]; then
        echo "Sphering FAILED (exit code $exit_code)"
        failed_plates+=("sphering:whole_screen")
    else
        echo "Sphering done for whole screen"
    fi
fi

# -----------------------------
# Step 4: single-cell processing (one plate/process at a time; skips plates
# without QC, and already-processed plates)
# -----------------------------
echo "======================================"
echo "Step 4: 4.single_cell_processing.ipynb (one plate at a time)"
echo "======================================"

for plate_id in "${plate_ids[@]}"; do
    sc_output="./data/single_cell_profiles/${plate_id}_sc_feature_selected.parquet"

    if [ -f "$sc_output" ] && [ -z "$overwrite" ]; then
        echo "✅ ${plate_id} already processed (found ${sc_output})"
        continue
    fi

    echo ">>> Running single-cell processing for ${plate_id}"
    PLATE_ID="$plate_id" OVERWRITE="$overwrite" python nbconverted/4.single_cell_processing.py
    exit_code=$?

    if [ "$exit_code" -ne 0 ]; then
        echo "Single-cell processing FAILED for plate: ${plate_id} (exit code $exit_code)"
        failed_plates+=("single_cell_processing:${plate_id}")
    else
        echo "Single-cell processing done for plate: ${plate_id}"
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
