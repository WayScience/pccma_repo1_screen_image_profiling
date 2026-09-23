#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --account=amc-general
#SBATCH --mem=256G
#SBATCH --time=03:00:00
#SBATCH --output=run_pipeline_child-%j.out

# activate preprocessing environment
module load miniforge
conda init bash
conda activate pccma_repo1_preprocessing_env

# prioritize the env's own libstdc++ over the system one in /lib64, which is
# older and missing symbols (e.g. GLIBCXX_3.4.29) required by libzmq.so.5
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

# plate id passed as first argument
plate_id=$1
overwrite="${OVERWRITE:-}"

failed_steps=()

# Step 1: single-cell QC via papermill (skips if already done)
qc_output="./data/qc_results/${plate_id}_qc_annotations.parquet"
qc_notebook_dir="./executed_notebooks/1.single_cell_qc"
mkdir -p "$qc_notebook_dir"

if [ -f "$qc_output" ] && [ -z "$overwrite" ]; then
    echo "✅ ${plate_id} already QC'd (found ${qc_output})"
else
    papermill 1.single_cell_qc.ipynb "${qc_notebook_dir}/${plate_id}.ipynb" \
        -p plate_id "$plate_id" \
        -p render_diagnostics False
    if [ $? -ne 0 ]; then
        echo "QC FAILED for plate: $plate_id"
        failed_steps+=("qc")
    fi
fi

# Step 2: annotation
PLATE_ID="$plate_id" OVERWRITE="$overwrite" python nbconverted/2.annotate.py
if [ $? -ne 0 ]; then
    echo "Annotation FAILED for plate: $plate_id"
    failed_steps+=("annotate")
fi

# Step 3: bulk processing (per-plate aggregate/normalize/feature-select --
# sphering is not run here; run 3b.sphering.ipynb locally instead)
PLATE_ID="$plate_id" OVERWRITE="$overwrite" python nbconverted/3.bulk_processing.py
if [ $? -ne 0 ]; then
    echo "Bulk processing FAILED for plate: $plate_id"
    failed_steps+=("bulk_processing")
fi

# Step 4: single-cell normalization
PLATE_ID="$plate_id" OVERWRITE="$overwrite" python nbconverted/4.single_cell_normalize.py
if [ $? -ne 0 ]; then
    echo "Single-cell normalization FAILED for plate: $plate_id"
    failed_steps+=("single_cell_normalize")
fi

# Step 5: single-cell feature selection
PLATE_ID="$plate_id" OVERWRITE="$overwrite" python nbconverted/5.single_cell_feature_select.py
if [ $? -ne 0 ]; then
    echo "Single-cell feature selection FAILED for plate: $plate_id"
    failed_steps+=("single_cell_feature_select")
fi

conda deactivate

if [ ${#failed_steps[@]} -gt 0 ]; then
echo "Plate $plate_id FAILED at steps: ${failed_steps[*]}"
exit 1
fi

echo "Plate $plate_id completed all steps successfully."
