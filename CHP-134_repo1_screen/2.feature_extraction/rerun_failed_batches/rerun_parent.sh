#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --account=amc-general
#SBATCH --time=0-00:2:00
#SBATCH --output=cp_rerun_parent-%j.out

# Resubmits the row batches listed in rerun_manifest.csv (the 5 batches that
# failed on the first full run) via rerun_child.sh, which requests more
# RAM/time than the original cp_analysis_hpc_child.sh.
# Delete this whole rerun_failed_batches/ folder once the reruns succeed.

# activate cellprofiler environment
module load miniforge
conda init bash
conda activate pccma_repo1_cp_env

# assumes this job's (and thus sbatch's) working directory is
# rerun_failed_batches/, i.e. sbatch was run from inside this folder
rerun_manifest="./rerun_manifest.csv"

# loop over each failed row batch and submit rerun child jobs
tail -n +2 "$rerun_manifest" | while IFS=, read -r plate row batch_label loaddata_file first_image_set last_image_set image_set_count well_count is_contiguous status message; do
    if [ "$status" != "ready" ]; then
        echo "Skipping ${plate} ${batch_label}: ${status} ${message}"
        continue
    fi

    echo "Resubmitting ${plate} ${batch_label}: image sets ${first_image_set}-${last_image_set} (${image_set_count} image sets, ${well_count} wells)"
    sbatch rerun_child.sh "$loaddata_file" "$first_image_set" "$last_image_set" "$batch_label"
done

conda deactivate

echo "All rerun CellProfiler jobs submitted!"
