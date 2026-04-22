#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=22G
#SBATCH --partition=amilan
#SBATCH --qos=long
#SBATCH --account=amc-general
#SBATCH --time=7-00:00:00
#SBATCH --output=run_CP_child-%j.out

# 1 task at 22GB RAM for the core (adjust as needed)

# Each well/fov (image-sets) needs about 2 GB of RAM, it fails at 10 GB after approximately 2,313 image sets.
# (about at well K18 per plate)
# We used all 384 wells of the plate with 9 FOVs each plate, so K18 is approximately 2,313 image sets 
# (384 wells x 6 fovs per well = 2,304 image sets, so K18 is the next one after that).
# The 8 extra GB is for memory leakage and overhead, which we calculated from 8 GB of excess divided by 2,313 image sets, 
# which is approximately 3.5 MB per image set.
# Over the total image sets we expect (384 wells x 9 FOVs = 3,456 image sets), we would need about 12 GB of RAM to run the entire plate without crashing
# plus 2 GB to run CellProfiler.
# I am requesting 22 GB of RAM to be safe, in case of underestimation.

# activate cellprofiler environment
module load miniforge
conda init bash
conda activate pccma_repo1_cp_env

# input csv  passed as first argument
csv=$1

# run your python analysis script with the input csv
python nbconverted/1.cp_analysis_hpc.py --input_csv "$csv"

# deactivate conda environment
conda deactivate

echo "CellProfiler analysis done for directory: $csv"
