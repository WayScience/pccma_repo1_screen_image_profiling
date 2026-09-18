#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --account=amc-general
#SBATCH --time=15:00
#SBATCH --output=run_pipeline_parent-%j.out

# activate preprocessing environment
module load miniforge
conda init bash
conda activate pccma_repo1_preprocessing_env

# convert all notebooks to python scripts (if any exist)
jupyter nbconvert --to=script --FilesWriter.build_directory=nbconverted/ *.ipynb

# Converted profiles (one parquet per plate) live on the PetaLibrary "koala"
# mount on Alpine. This must match the Alpine-detection branch in
# 1.single_cell_qc.ipynb through 5.single_cell_feature_select.ipynb.
converted_dir="/pl/active/koala/ALSF_screen_data/SK-N-AS_repo1_profiles/converted_profiles"

mapfile -t plate_ids < <(find "$converted_dir" -maxdepth 1 -name "*_converted.parquet" -printf "%f\n" | sed 's/_converted\.parquet$//' | sort)

echo "Number of plates found: ${#plate_ids[@]}"
for plate_id in "${plate_ids[@]}"; do
    echo "Found: $plate_id"
done

# loop over each plate and submit a child job. Runs steps 1-5 only --
# sphering is not run on HPC; run 3b.sphering.ipynb locally once every
# plate's bulk-processed output has been synced back.
for plate_id in "${plate_ids[@]}"; do
    # check job count for this user
    number_of_jobs=$(squeue -u "$USER" | wc -l)
    while [ "$number_of_jobs" -gt 990 ]; do
        sleep 1s
        number_of_jobs=$(squeue -u "$USER" | wc -l)
    done
    sbatch run_pipeline_hpc_child.sh "$plate_id"
done

conda deactivate

echo "All run_pipeline jobs submitted!"
