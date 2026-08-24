#!/bin/bash
#
# Run 0.convert_cytotable.ipynb locally, once per plate, sequentially.
# Each invocation converts all 16 row-batch SQLite files for that plate
# (the row-batch loop lives inside 0.convert_cytotable.py itself).
#
# Run this script from within CHP-134_repo1_screen/3.process_features/.

set -uo pipefail

# -----------------------------
# Initialize environment
# -----------------------------
conda init bash
conda activate pccma_repo1_preprocessing_env

# convert notebooks to scripts
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# -----------------------------
# Discover plates from the SQLite source directory
# -----------------------------
sqlite_dir="$HOME/mnt/bandicoot/PCCMA_data/CHP-134_repo1_screen_outputs/SQLite_outputs"

mapfile -t plate_ids < <(find "$sqlite_dir" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort)

echo "Number of plates found: ${#plate_ids[@]}"
for plate_id in "${plate_ids[@]}"; do
    echo "Found: $plate_id"
done

# expected number of converted parquet files per plate (16 row batches)
expected_batch_count=16
converted_dir="./data/converted_profiles"

failed_plates=()

# -----------------------------
# Loop over all plates sequentially
# -----------------------------
for plate_id in "${plate_ids[@]}"; do
    echo "======================================"
    echo "Processing plate: $plate_id"
    echo "======================================"

    parquet_count=0
    if [ -d "$converted_dir" ]; then
        parquet_count=$(find "$converted_dir" -maxdepth 1 -name "${plate_id}_row_*_converted.parquet" | wc -l)
    fi

    if [ "$parquet_count" -eq "$expected_batch_count" ]; then
        echo "✅ Plate ${plate_id} already fully converted (${parquet_count}/${expected_batch_count} row batches)"
        continue
    fi

    echo ">>> Running CytoTable conversion for ${plate_id} (${parquet_count}/${expected_batch_count} row batches done)"
    python nbconverted/0.convert_cytotable.py --plate_id "$plate_id"
    exit_code=$?

    if [ "$exit_code" -ne 0 ]; then
        echo "CytoTable conversion FAILED for plate: $plate_id (exit code $exit_code)"
        failed_plates+=("$plate_id")
    else
        echo "CytoTable conversion done for plate: $plate_id"
    fi
done

conda deactivate

echo "======================================"
echo "All plates processed."
if [ ${#failed_plates[@]} -gt 0 ]; then
    echo "Failed plates:"
    printf '  %s\n' "${failed_plates[@]}"
    exit 1
fi

echo "✅ All plates converted successfully."
