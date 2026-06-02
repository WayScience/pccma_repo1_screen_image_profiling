#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6G
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --account=amc-general
#SBATCH --time=10:00:00
#SBATCH --output=run_CP_child-%j.out

# 1 task at 6 GB RAM to run row-level batches on HPC.

# Each well/fov (image-sets) needs about 2 GB of RAM, it fails at 10 GB after approximately 2,313 image sets.
# (about at well K18 per plate)
# When running at a row-batch level, we are processing 24 wells with 216 image sets (9 FOVs per well)
# which we estimate would need potentially a maximum of 4 GB of RAM 
# (216 image sets x 3.5 MB per image set = 756 MB, plus 2 GB for CellProfiler, plus overhead).
# I am requesting 6 GB of RAM at 10 hours to be safe, in case of underestimation.

# activate cellprofiler environment
module load miniforge
conda init bash
conda activate pccma_repo1_cp_env

# input csv and optional row batch range passed as arguments
csv=$1
first_image_set=$2
last_image_set=$3
batch_label=$4

# run your python analysis script with the input csv
command=(python nbconverted/1.cp_analysis_hpc.py --input_csv "$csv")

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

echo "CellProfiler analysis done for directory: $csv ($batch_label)"
