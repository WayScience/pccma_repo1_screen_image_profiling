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

mapfile -t plate_ids < <(find "$converted_dir" -maxdepth 1 -name "*_converted.parquet" -printf "%f\n" | sed 's/_converted\.parquet$//' | sort)

# optionally skip plates that are already finished, so they don't take up a
# job (each child requests a full node's worth of memory just to print
# "Skipping"). Expects a pace-separated list
exclude_plates="${EXCLUDE_PLATES:-}"
if [ -n "$exclude_plates" ]; then
    filtered_plate_ids=()
    for plate_id in "${plate_ids[@]}"; do
        if [[ " $exclude_plates " == *" $plate_id "* ]]; then
            echo "Excluding: $plate_id"
        else
            filtered_plate_ids+=("$plate_id")
        fi
    done
    plate_ids=("${filtered_plate_ids[@]}")
fi

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
