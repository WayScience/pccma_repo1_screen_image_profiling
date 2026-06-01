#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --account=amc-general
#SBATCH --time=15:00
#SBATCH --output=cp_parent-%j.out

# activate cellprofiler environment
module load miniforge
conda init bash
conda activate pccma_repo1_cp_env

# convert all notebooks to python scripts (if any exist)
jupyter nbconvert --to=script --FilesWriter.build_directory=nbconverted/ *.ipynb

# run the LoadData CSV creation script once before submitting jobs
python nbconverted/0.create_loaddata_csvs.py --HPC

# use row batch ranges created with the LoadData CSVs
row_batch_manifest="./row_batch_manifest.csv"

ready_batches=$(awk -F, 'NR > 1 && $10 == "ready" { count++ } END { print count + 0 }' "$row_batch_manifest")
skipped_batches=$(awk -F, 'NR > 1 && $10 != "ready" { count++ } END { print count + 0 }' "$row_batch_manifest")

echo "Ready row batches: $ready_batches"
echo "Skipped or warning row batches: $skipped_batches"

# loop over each ready row batch and submit child jobs
tail -n +2 "$row_batch_manifest" | while IFS=, read -r plate row batch_label loaddata_file first_image_set last_image_set image_set_count well_count is_contiguous status message; do
    # TEST MODE: only run a specific plate and batch
    # if [ "$plate" != "BR00149332" ] || [ "$batch_label" != "row_A" ]; then
    #     continue
    # fi
    
    if [ "$status" != "ready" ]; then
        echo "Skipping ${plate} ${batch_label}: ${status} ${message}"
        continue
    fi

    # check job count for this user
    number_of_jobs=$(squeue -u "$USER" | wc -l)
    while [ "$number_of_jobs" -gt 990 ]; do
        sleep 1s
        number_of_jobs=$(squeue -u "$USER" | wc -l)
    done
    echo "Submitting ${plate} ${batch_label}: image sets ${first_image_set}-${last_image_set} (${image_set_count} image sets, ${well_count} wells)"
    sbatch cp_analysis_hpc_child.sh "$loaddata_file" "$first_image_set" "$last_image_set" "$batch_label"
done

conda deactivate

echo "All CellProfiler jobs submitted!"
