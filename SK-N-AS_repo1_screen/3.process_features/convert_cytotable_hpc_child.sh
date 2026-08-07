#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --partition=acpu
#SBATCH --qos=cpu-long
#SBATCH --account=amc-general
#SBATCH --time=1-06:00:00
#SBATCH --output=convert_cytotable_child-%j.out

# NOTE on --mem/--time: estimated (not measured on Alpine) from a prior local run of
# this same join workload (parquet backend, no image export) on plate BR00148919
# (54.6GB SQLite), which took ~18h14m under cytotable's default (largely serial,
# single-slot) parsl config, while sitting at ~3.75GB RSS during that phase.
# The 29 plates here range 28.9-58.3GB; scaling the ~18h14m/~3.75GB figures to the
# largest plate (58.3GB) gives ~19.5h / ~4GB. --time is padded to 30h and --mem to
# 16G for margin.

# activate preprocessing environment (includes cytotable)
module load miniforge
conda init bash
conda activate pccma_repo1_preprocessing_env

# plate id passed as first argument
plate_id=$1

python nbconverted/0.convert_cytotable.py --plate_id "$plate_id"

conda deactivate

echo "CytoTable conversion done for plate: $plate_id"
