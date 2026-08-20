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
            --time=08:00:00 \
            --cpus-per-task=24 \
            --mem=64G \
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
# og adapter
# P5_ADAPTER="${P5_ADAPTER:-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTTTTTCGAGATGTGTATAAGAGACAG}" # add a G
# adapter containing complete transposon
P5_ADAPTER="${P5_ADAPTER:-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTTTTTCGAGATGTGTATAAGAGACAGG}" 
# Tn5 transposon: # AGATGTGTATAAGAGACAGG
TRANSPOSON="${TRANSPOSON:-AGATGTGTATAAGAGACAGG}"

# Remove the protocol's 9-base inline barcode from the 5' end by default.
# Override with TRIM_5P=0 to retain it, or set the appropriate barcode length.
TRIM_5P="${TRIM_5P:-0}"

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

if [[ ! "$TRANSPOSON" =~ ^[ACGTacgt]+$ ]]; then
    echo "ERROR: TRANSPOSON must contain only A, C, G, and T."
    exit 1
fi
TRANSPOSON=$(printf '%s' "$TRANSPOSON" | tr '[:lower:]' '[:upper:]')


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
TRANSPOSON_FASTQ="${TRIMMED_DIR}/${SAMPLE}_transposon_trimmed.fastq.gz"
FINAL_FASTQ="${FINAL_DIR}/${SAMPLE}_filtered.fastq.gz"


# ============================================================
# Step 1: Remove adapters and homopolymer tails
# ============================================================

echo
echo "Step 1: Removing adapters and 3' poly-C/poly-G tails from $SAMPLE..."

# Previous adapter and homopolymer-tail trimming command:
    # -g "X${P5_ADAPTER}"  # this contains the transposon will try like this and then with p5 instead of transposon
CUTADAPT_ARGS=(
    -a "${P3_ADAPTER}"
    -a "C{${C_TAIL_LENGTH}}"
    -a "G{${G_TAIL_LENGTH}}"
    --times 3
    --minimum-length 16
    -o "$TRIMMED_FASTQ"
    -j 24
)

if (( TRIM_5P > 0 )); then
    echo "Also removing $TRIM_5P base(s) from the 5' end."
    CUTADAPT_ARGS+=( -u "$TRIM_5P" )
fi

cutadapt "${CUTADAPT_ARGS[@]}" "$INPUT"


# ============================================================
# Step 2: Require and remove the 5' transposon
# ============================================================

echo
echo "Step 2: Requiring and removing the 5' transposon from $SAMPLE..."
echo "Transposon: $TRANSPOSON"
echo "Reads without this 5' transposon will be discarded."

# Keeping this separate makes --discard-untrimmed apply specifically to the
# transposon. Cutadapt prints the total, retained, and discarded read counts
# in the report below.
cutadapt \
    -g "${TRANSPOSON}" \
    --discard-untrimmed \
    -o "$TRANSPOSON_FASTQ" \
    -j 24 \
    --minimum-length 16 \
    "$TRIMMED_FASTQ"
#params for 5' removal
# -g "${P5_ADAPTER}" \
# params for piepline removal
# cp "$TRIMMED_FASTQ" "$TRANSPOSON_FASTQ"
    

# # Trim the protocol's 3' adapter and terminal poly-C/poly-G tails. Require a
# # 20-base overlap for each match to avoid trimming short, natural C/G-rich ends.
# # Retain only adapter-trimmed reads that meet the requested length and
# # ambiguous-base thresholds.
# CUTADAPT_5P_ARGS=()
# if (( TRIM_5P > 0 )); then
#     echo "Also removing $TRIM_5P base(s) from the 5' end."
#     CUTADAPT_5P_ARGS+=( -u "$TRIM_5P" )
# fi

# # cutadapt \
# #     -j 4 \
# #     "${CUTADAPT_5P_ARGS[@]}" \
# #     -a "${P3_ADAPTER}" \
# #     -a "C{100};min_overlap=10" \
# #     -a "G{100};min_overlap=10" \
# #     --output "$TRIMMED_FASTQ" \ls
# #     --error-rate 0.1 \
# #     --times 3 \
# #     --overlap 20 \
# #     --action trim \
# #     --discard-untrimmed \
# #     --minimum-length 16 \
# #     --max-n 0.5 \
# #     "$INPUT"


# ============================================================
# Step 3: Quality filtering
# ============================================================

echo
echo "Step 3: Quality filtering $SAMPLE..."

FINAL_FASTQ_TMP="${FINAL_FASTQ}.tmp.$$"
trap 'rm -f "$FINAL_FASTQ_TMP"' EXIT

gzip -dc -- "$TRANSPOSON_FASTQ" \
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
echo removing intermediate files
rm -r "$TRIMMED_DIR"

echo
echo "========================================"
echo "Finished: $INPUT"
echo "Finished at: $(date)"
echo
echo "Final filtered reads:"
echo "$FINAL_FASTQ"
echo "========================================"