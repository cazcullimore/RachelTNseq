#!/usr/bin/env bash

# ============================================================
# If this script is NOT currently running inside SLURM,
# submit one separate SLURM job for every demultiplexed FASTQ.
# ============================================================
module load miniforge3
conda activate TnSeqPreProcess
if [[ -z "${SLURM_JOB_ID:-}" ]]; then

    if [[ $# -lt 1 ]]; then
        echo "Usage:"
        echo "  ./preprocess.sh sample1.fastq.gz sample2.fastq.gz sample3.fastq.gz"
        echo
        echo "You can also use:"
        echo "  ./preprocess.sh *.fastq.gz"
        exit 1
    fi

    echo "Submitting ${#} input file(s) as separate SLURM jobs..."
    echo

    for INPUT in "$@"; do

        if [[ ! -f "$INPUT" ]]; then
            echo "WARNING: File not found, skipping: $INPUT"
            continue
        fi

        echo "Submitting: $INPUT"

        sbatch \
            --job-name="preprocess_$(basename "$INPUT" .fastq.gz)" \
            --output="preprocess_%j.out" \
            --error="preprocess_%j.err" \
            --time=04:00:00 \
            --cpus-per-task=4 \
            --mem=16G \
            "$0" "$INPUT"
    done

    echo
    echo "All jobs submitted."
    exit 0
fi


# ============================================================
# Everything below here runs INSIDE each individual SLURM job.
# ============================================================

set -euo pipefail

INPUT="$1"

echo "========================================"
echo "SLURM job ID: $SLURM_JOB_ID"
echo "Running on: $(hostname)"
echo "Input FASTQ: $INPUT"
echo "Started: $(date)"
echo "========================================"


# ============================================================
# Load software
# ============================================================

# Change these to match your cluster.
#
# module load fastx-toolkit
# module load cutadapt
#
# OR activate a conda environment:
#
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate sequencing


# ============================================================
# Settings
# ============================================================

# Cutadapt uses a long repeated-base adapter to find the longest possible
# variable-length poly-C tail. C_TAIL_LENGTH is the search length, not a
# requirement that every read contain this many Cs (minimum overlap is 3).
C_TAIL_LENGTH="${C_TAIL_LENGTH:-100}"

# Trim variable-length poly-G tails as well. These commonly occur near the
# 3' end of two-color Illumina reads when the signal becomes weak.
G_TAIL_LENGTH="${G_TAIL_LENGTH:-100}"

# Protocol-specific adapters. At the 5' end, the leading X supplied to
# Cutadapt below makes P5 non-internal: full or partial adapter sequence is
# removed at the read start, but an internal match is not removed. P3 is
# searched for at the 3' end of short-insert reads.
P3_ADAPTER="${P3_ADAPTER:-GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTT}"
P5_ADAPTER="${P5_ADAPTER:-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTTTTTCGAGATGTGTATAAGAGACAG}"

# Demultiplexed input is assumed to have its inline barcode removed already.
# If the barcode is still present at the 5' end, run with TRIM_5P=9 (or the
# appropriate barcode length).
TRIM_5P="${TRIM_5P:-9}"

QUALITY=20
PERCENT=90


# ============================================================
# Check input files
# ============================================================

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: Input file not found:"
    echo "$INPUT"
    exit 1
fi

if [[ ! "$C_TAIL_LENGTH" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: C_TAIL_LENGTH must be a positive integer."
    exit 1
fi

if [[ ! "$G_TAIL_LENGTH" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: G_TAIL_LENGTH must be a positive integer."
    exit 1
fi

if [[ ! "$TRIM_5P" =~ ^[0-9]+$ ]]; then
    echo "ERROR: TRIM_5P must be a non-negative integer."
    exit 1
fi

if [[ ! "$P3_ADAPTER" =~ ^[ACGTacgt]+$ ]]; then
    echo "ERROR: P3_ADAPTER must contain only A, C, G, and T."
    exit 1
fi

if [[ ! "$P5_ADAPTER" =~ ^[ACGTacgt]+$ ]]; then
    echo "ERROR: P5_ADAPTER must contain only A, C, G, and T."
    exit 1
fi


# ============================================================
# Determine input sample name
# ============================================================

BASENAME=$(basename "$INPUT")

SAMPLE="${BASENAME%.fastq.gz}"
SAMPLE="${SAMPLE%.fq.gz}"
SAMPLE="${SAMPLE%.fastq}"
SAMPLE="${SAMPLE%.fq}"

echo "Sample name: $SAMPLE"


# ============================================================
# Create separate directories for this sample
# ============================================================

RESULTS_DIR="${RESULTS_DIR:-results}"
OUTDIR="${RESULTS_DIR}/${SAMPLE}"

TRIMMED_DIR="${OUTDIR}/trimmed"
FINAL_DIR="${OUTDIR}/final"

mkdir -p \
    "$TRIMMED_DIR" \
    "$FINAL_DIR"

echo "Output directory: $OUTDIR"

TRIMMED_FASTQ="${TRIMMED_DIR}/${SAMPLE}_adapter_ctail_trimmed.fastq.gz"
FINAL_FASTQ="${FINAL_DIR}/${SAMPLE}_filtered.fastq.gz"


# ============================================================
# Step 1: Remove adapters and homopolymer tails
# ============================================================

echo
echo "Step 1: Removing adapters and 3' poly-C/poly-G tails from $SAMPLE..."

CUTADAPT_ARGS=(
    -g "X${P5_ADAPTER}"
    -a "${P3_ADAPTER}"
    -a "C{${C_TAIL_LENGTH}}"
    -a "G{${G_TAIL_LENGTH}}"
    --times 3
    --minimum-length 16
    -o "$TRIMMED_FASTQ"
)

if (( TRIM_5P > 0 )); then
    echo "Also removing $TRIM_5P base(s) from the 5' end."
    CUTADAPT_ARGS+=( -u "$TRIM_5P" )
fi

cutadapt "${CUTADAPT_ARGS[@]}" "$INPUT"


# ============================================================
# Step 2: Quality filtering
# ============================================================

echo
echo "Step 2: Quality filtering $SAMPLE..."

FINAL_FASTQ_TMP="${FINAL_FASTQ}.tmp.$$"
trap 'rm -f "$FINAL_FASTQ_TMP"' EXIT

gzip -dc -- "$TRIMMED_FASTQ" \
    | fastq_quality_filter \
        -q "$QUALITY" \
        -p "$PERCENT" \
        -Q33 \
    | gzip -c > "$FINAL_FASTQ_TMP"

mv "$FINAL_FASTQ_TMP" "$FINAL_FASTQ"
trap - EXIT


# ============================================================
# Finished
# ============================================================

echo
echo "========================================"
echo "Finished: $INPUT"
echo "Finished at: $(date)"
echo
echo "Final filtered reads:"
echo "$FINAL_FASTQ"
echo "========================================"
