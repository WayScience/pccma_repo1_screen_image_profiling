#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --account=amc-general
#SBATCH --time=14:00:00
#SBATCH --output=run_CP_rerun_child-%j.out

# Same as ../cp_analysis_hpc_child.sh, but with more RAM and time for the
# row batches that failed on the first full run at 4G / 9:00:00.
# Delete this whole rerun_failed_batches/ folder once the reruns succeed.

# activate cellprofiler environment
module load miniforge
conda init bash
conda activate pccma_repo1_cp_env

# input csv and row batch range passed as arguments
csv=$1
first_image_set=$2
last_image_set=$3
batch_label=$4

# assumes this job's working directory is rerun_failed_batches/, i.e. sbatch
# was run from inside this folder (Slurm jobs inherit the submitter's cwd) --
# nbconverted/ lives one level up, in 2.feature_extraction/
command=(python ../nbconverted/1.cp_analysis_hpc.py --input_csv "$csv")

if [ -n "$first_image_set" ]; then
    command+=(--first_image_set "$first_image_set")
fi
if [ -n "$last_image_set" ]; then
    command+=(--last_image_set "$last_image_set")
fi
if [ -n "$batch_label" ]; then
    command+=(--batch_label "$batch_label")
fi

"${command[@]}"

# deactivate conda environment
conda deactivate

echo "CellProfiler rerun done for directory: $csv ($batch_label)"
