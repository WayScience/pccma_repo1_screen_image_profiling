#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --account=amc-general
#SBATCH --time=5:00
#SBATCH --output=convert_cytotable_parent_test-%j.out

# TEMPORARY test script -- submits a single child job (instead of looping over
# all plates like convert_cytotable_hpc_parent.sh) so the new parsl
# worker-count config in 0.convert_cytotable.ipynb can be validated end-to-end
# on one real plate before committing the full 29-plate batch to it. Check
# `seff <child_jobid>` once it finishes to get real wall time / peak memory,
# then right-size convert_cytotable_hpc_child.sh's --time/--mem accordingly.
# Safe to delete once that validation run is done.

# activate preprocessing environment (includes cytotable)
module load miniforge
conda init bash
conda activate pccma_repo1_preprocessing_env

# convert all notebooks to python scripts (if any exist)
jupyter nbconvert --to=script --FilesWriter.build_directory=nbconverted/ *.ipynb

# SQLite outputs (one subdirectory per plate) live on the PetaLibrary "koala"
# mount on Alpine. This must match the sqlite_dir HPC branch in
# 0.convert_cytotable.ipynb.
sqlite_dir="/pl/active/koala/ALSF_screen_data/SK-N-AS_repo1_profiles/SQLite_outputs"

# plate id can be passed as an argument (sbatch convert_cytotable_hpc_parent_test.sh BR00148919);
# otherwise default to the first plate found (sorted) so this runs with no args too.
if [ -n "$1" ]; then
    plate_id="$1"
else
    plate_id=$(find "$sqlite_dir" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort | head -n 1)
fi

echo "Test run: submitting a single child job for plate: $plate_id"
sbatch convert_cytotable_hpc_child.sh "$plate_id"

conda deactivate

echo "Test CytoTable conversion job submitted for plate: $plate_id"
