#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --account=amc-general
#SBATCH --time=15:00
#SBATCH --output=run_pipeline_parent_test-%j.out

# TEMPORARY test script -- submits a single run_pipeline_hpc_child.sh job
# (instead of looping over all plates like run_pipeline_hpc_parent.sh) so the
# --mem/--time/--partition in run_pipeline_hpc_child.sh can be validated
# end-to-end on one real plate before committing the full 29-plate batch to
# it. Check `seff <child_jobid>` once it finishes to get real wall time /
# peak memory, then right-size run_pipeline_hpc_child.sh's --time/--mem
# accordingly. Safe to delete once that validation run is done.

# activate preprocessing environment
module load miniforge
conda init bash
conda activate pccma_repo1_preprocessing_env

# convert all notebooks to python scripts (if any exist)
jupyter nbconvert --to=script --FilesWriter.build_directory=nbconverted/ *.ipynb

# Converted profiles (one parquet per plate) are meant to live on the
# PetaLibrary "koala" mount on Alpine, but fall back to the relative
# data/converted_profiles directory in this repo checkout on scratch, in
# case a plate hasn't been synced to koala yet. Must match the
# Alpine-detection branch in 1.single_cell_qc.ipynb through
# 5.single_cell_feature_select.ipynb.
koala_dir="/pl/active/koala/ALSF_screen_data/SK-N-AS_repo1_profiles/converted_profiles"
scratch_dir="data/converted_profiles"
if [ -d "$koala_dir" ]; then
    converted_dir="$koala_dir"
else
    converted_dir="$scratch_dir"
fi

# plate id can be passed as an argument (sbatch run_pipeline_hpc_parent_test.sh BR00148919);
# otherwise default to the first plate found (sorted) so this runs with no args too.
if [ -n "$1" ]; then
    plate_id="$1"
else
    plate_id=$(find "$converted_dir" -maxdepth 1 -name "*_converted.parquet" -printf "%f\n" | sed 's/_converted\.parquet$//' | sort | head -n 1)
fi

echo "Test run: submitting a single child job for plate: $plate_id"
sbatch run_pipeline_hpc_child.sh "$plate_id"

conda deactivate

echo "Test run_pipeline job submitted for plate: $plate_id"
