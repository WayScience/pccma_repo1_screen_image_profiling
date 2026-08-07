#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --account=amc-general
#SBATCH --time=5:00
#SBATCH --output=convert_cytotable_parent-%j.out

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

mapfile -t plate_ids < <(find "$sqlite_dir" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort)

echo "Number of plates found: ${#plate_ids[@]}"
for plate_id in "${plate_ids[@]}"; do
    echo "Found: $plate_id"
done

# loop over each plate and submit a child job
for plate_id in "${plate_ids[@]}"; do
    # check job count for this user
    number_of_jobs=$(squeue -u "$USER" | wc -l)
    while [ "$number_of_jobs" -gt 990 ]; do
        sleep 1s
        number_of_jobs=$(squeue -u "$USER" | wc -l)
    done
    sbatch convert_cytotable_hpc_child.sh "$plate_id"
done

conda deactivate

echo "All CytoTable conversion jobs submitted!"
