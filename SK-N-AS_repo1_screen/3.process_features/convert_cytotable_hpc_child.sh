#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --account=amc-general
#SBATCH --time=24:00:00
#SBATCH --output=convert_cytotable_child-%j.out

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
