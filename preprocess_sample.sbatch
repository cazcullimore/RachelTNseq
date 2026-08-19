#!/usr/bin/env bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G

set -eo pipefail

module load miniforge3
conda activate TnSeqPreProcess
set -u
# Don't use BASH_SOURCE for sibling scripts under Slurm.
if (( $# != 3 )); then
    echo "Usage: sbatch preprocess_sample.sbatch READS PREPROCESS_DIR SCRIPT_DIR" >&2
    exit 2
fi
READS=$1
PREPROCESS_DIR=$2
SCRIPT_DIR=$3

if [[ ! -f "$READS" ]]; then
    echo "ERROR: Reads file not found: $READS" >&2
    exit 1
fi
if [[ ! -d "$PREPROCESS_DIR" ]]; then
    echo "ERROR: Preprocessing directory not found: $PREPROCESS_DIR" >&2
    exit 1
fi

echo "Preprocessing Slurm job: ${SLURM_JOB_ID:-unknown}"
echo "script dir: ${SCRIPT_DIR}"
RESULTS_DIR="$PREPROCESS_DIR" bash "${SCRIPT_DIR}/preprocess.sh" "$READS"
