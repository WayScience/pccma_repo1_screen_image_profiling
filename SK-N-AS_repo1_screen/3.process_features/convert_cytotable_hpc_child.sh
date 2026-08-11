#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --account=amc-general
#SBATCH --time=20:00:00
#SBATCH --output=convert_cytotable_child-%j.out

# NOTE on --cpus-per-task/--mem/--time:
# --time was originally padded to 30h (cpu-long) based on a local single-worker
# run of this same join workload (parquet backend, no image export) on plate
# BR00148919 (54.6GB SQLite), which took ~18h14m under cytotable's default
# (largely serial, single-slot) parsl config, scaled to the largest of the 29
# plates (58.3GB) at ~19.5h.
# A subsequent Alpine run under the fixed multi-worker config (16 workers, 32G)
# showed the pre-join stage completing in ~1h16m, but the job then failed
# immediately at the joins step with an OOM. This means the bulk of the
# previously-estimated runtime (the old serial single-worker figure) no longer
# reflects actual behavior under parallel workers -- the pre-join stage is fast,
# and the joins step is memory-bound rather than time-bound. --time is reduced
# to 20h and --qos switched to cpu-normal accordingly; this remains a rough
# margin rather than a measured full-pipeline runtime, since a successful
# end-to-end run (through the joins step) hasn't yet completed.
#
# --mem was originally 16G, based on a measured ~3.75GB RSS for the single-worker
# case above. Under the new multi-worker config, a local benchmark (91k-cell
# synthetic test, 16 workers/threads) measured ~5.7GB RSS -- more workers process
# chunks concurrently instead of one at a time, so memory scales with worker count.
# The 32G derived from that benchmark was confirmed insufficient on Alpine (OOM
# at the joins step, per above). --mem is now set to 100G for the 8-worker config
# as a conservative placeholder pending a run that succeeds through the joins step;
# it has not itself been derived from a benchmark at 8 workers.

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
