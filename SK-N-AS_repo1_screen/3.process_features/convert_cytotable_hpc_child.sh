#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=acpu
#SBATCH --qos=cpu-long
#SBATCH --account=amc-general
#SBATCH --time=1-06:00:00
#SBATCH --output=convert_cytotable_child-%j.out

# NOTE on --cpus-per-task/--mem/--time:
# --time is estimated (not measured on Alpine) from a prior local run of this same
# join workload (parquet backend, no image export) on plate BR00148919 (54.6GB
# SQLite), which took ~18h14m under cytotable's default (largely serial,
# single-slot) parsl config. The 29 plates here range 28.9-58.3GB; scaling to the
# largest plate (58.3GB) gives ~19.5h. --time is padded to 30h for margin.
# This estimate predates the parsl worker-count fix below (0.convert_cytotable.ipynb
# now sizes parsl to --cpus-per-task instead of running effectively single-threaded),
# so the real runtime under 16 workers is not yet known -- --time is intentionally
# left unchanged (not shortened) until a real HPC run confirms a speedup, rather
# than risk a plate getting killed mid-conversion on an unverified assumption.
#
# --mem was originally 16G, based on a measured ~3.75GB RSS for the single-worker
# case above. Under the new multi-worker config, a local benchmark (91k-cell
# synthetic test, 16 workers/threads) measured ~5.7GB RSS -- more workers process
# chunks concurrently instead of one at a time, so memory scales with worker count.
# --mem is padded to 32G (~5.6x the measured parallel figure) for margin, since the
# benchmark was on a much smaller dataset than a real plate and the final
# concat/join step (not chunked) is untested at full scale.

# activate preprocessing environment (includes cytotable)
module load miniforge
conda init bash
conda activate pccma_repo1_preprocessing_env

# prioritize the env's own libstdc++ over the system one in /lib64, which is
# older and missing symbols (e.g. GLIBCXX_3.4.29) required by libzmq.so.5
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

# plate id passed as first argument
plate_id=$1

python nbconverted/0.convert_cytotable.py --plate_id "$plate_id"
exit_code=$?

conda deactivate

if [ "$exit_code" -ne 0 ]; then
    echo "CytoTable conversion FAILED for plate: $plate_id (exit code $exit_code)"
    exit "$exit_code"
fi

echo "CytoTable conversion done for plate: $plate_id"
