#!/usr/bin/env bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=TnSeqJob

set -euo pipefail

# Environment setup is intentionally left as-is for the cluster owner to adjust.
[[ -f ~/.bashrc ]] && source ~/.bashrc
mamba activate py27_env
module load bowtie2

usage() {
    cat <<'EOF'
Usage:
  sbatch run.sbatch READS.fastq[.gz] REFERENCE.fasta ANNOTATION.gff[.gz] [TRANSPOSON] [OUTPUT_DIR]

Arguments:
  READS.fastq[.gz]  Demultiplexed FASTQ input.
  REFERENCE.fasta   Reference FASTA used to build the Bowtie 2 index.
  ANNOTATION.gff    GFF/GFF3 annotation matching the reference FASTA.
  TRANSPOSON        Optional transposon sequence.
                     Default: AGATGTGTATAAGAGACAGG
  OUTPUT_DIR         Optional run directory.
                     Default: <sample>_tnseq

The job preprocesses the reads, builds the Bowtie 2 reference, converts the
GFF to legacy PTT files, writes a legacy tnSeq.config, and runs the pipeline.
EOF
}

if (( $# < 3 || $# > 5 )); then
    usage >&2
    exit 2
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

absolute_file() {
    local path="$1"
    local directory
    local filename

    directory=$(dirname -- "$path")
    filename=$(basename -- "$path")
    (cd -- "$directory" && printf '%s/%s\n' "$PWD" "$filename")
}

absolute_directory() {
    local path="$1"
    mkdir -p -- "$path"
    (cd -- "$path" && printf '%s\n' "$PWD")
}

READS_ARG="$1"
REFERENCE_ARG="$2"
GFF_ARG="$3"
TRANSPOSON="${4:-AGATGTGTATAAGAGACAGG}"

[[ "$READS_ARG" = /* ]] || READS_ARG="${SUBMIT_DIR}/${READS_ARG}"
[[ "$REFERENCE_ARG" = /* ]] || REFERENCE_ARG="${SUBMIT_DIR}/${REFERENCE_ARG}"
[[ "$GFF_ARG" = /* ]] || GFF_ARG="${SUBMIT_DIR}/${GFF_ARG}"

for input_path in "$READS_ARG" "$REFERENCE_ARG" "$GFF_ARG"; do
    if [[ ! -f "$input_path" ]]; then
        echo "ERROR: Input file not found: $input_path" >&2
        exit 1
    fi
done

if [[ ! "$TRANSPOSON" =~ ^[ACGTacgt]+$ ]]; then
    echo "ERROR: TRANSPOSON must contain only A, C, G, and T." >&2
    exit 1
fi
TRANSPOSON=${TRANSPOSON^^}

READS=$(absolute_file "$READS_ARG")
REFERENCE=$(absolute_file "$REFERENCE_ARG")
GFF=$(absolute_file "$GFF_ARG")

READS_NAME=$(basename -- "$READS")
SAMPLE="${READS_NAME%.fastq.gz}"
SAMPLE="${SAMPLE%.fq.gz}"
SAMPLE="${SAMPLE%.fastq}"
SAMPLE="${SAMPLE%.fq}"

OUTPUT_ARG="${5:-${SUBMIT_DIR}/${SAMPLE}_tnseq}"
[[ "$OUTPUT_ARG" = /* ]] || OUTPUT_ARG="${SUBMIT_DIR}/${OUTPUT_ARG}"
RUN_DIR=$(absolute_directory "$OUTPUT_ARG")

PREPROCESS_DIR="${RUN_DIR}/preprocessed"
PTT_DIR="${RUN_DIR}/ptt"
BOWTIE_DIR="${RUN_DIR}/bowtie_reference"
CONFIG_FILE="${RUN_DIR}/tnSeq.config"
INDEX_PREFIX="${BOWTIE_DIR}/${SAMPLE}"
FINAL_READS="${PREPROCESS_DIR}/${SAMPLE}/final/${SAMPLE}_filtered.fastq.gz"

mkdir -p -- "$PREPROCESS_DIR" "$PTT_DIR" "$BOWTIE_DIR"

echo "Run directory: $RUN_DIR"
echo "Preprocessing reads..."
RESULTS_DIR="$PREPROCESS_DIR" TRANSPOSON="$TRANSPOSON" \
    bash "${SCRIPT_DIR}/preprocess.sh" "$READS"

if [[ ! -f "$FINAL_READS" ]]; then
    echo "ERROR: Expected preprocessed reads were not created: $FINAL_READS" >&2
    exit 1
fi

echo "Building Bowtie 2 reference..."
bowtie2-build --threads "${SLURM_CPUS_PER_TASK:-8}" "$REFERENCE" "$INDEX_PREFIX"

echo "Creating legacy PTT annotation files..."
python3 "${SCRIPT_DIR}/gff_to_ptt.py" \
    --gff "$GFF" \
    --fasta "$REFERENCE" \
    --outdir "$PTT_DIR"

# preprocess.sh removes the transposon along with the protocol adapter, so the
# legacy pipeline must not try to find and remove that sequence a second time.
cat > "$CONFIG_FILE" <<EOF
[input]
Reads           = $FINAL_READS
BowtieReference = $INDEX_PREFIX
Ptt             = $PTT_DIR/*.ptt
Out             = $SAMPLE

[parameters]
Transposon      =
MinBaseOffset   = 0
MaxBaseOffset   = 0
Mismatches      = 0
GeneTrim        = 10
ReadLength      = 16
MinimumHopCount = 16

[options]
Debug                   = False
Normalize               = Total
DeleteIntermediateFiles = True
ReverseComplementReads  = False
IGVNormalize            = True
IGVNegateNegStrand      = False
EOF

echo "Legacy config: $CONFIG_FILE"
echo "Running: python2.7 TnSeq-Pipeline_SKCmodifications.py $CONFIG_FILE"
python2.7 "${SCRIPT_DIR}/TnSeq-Pipeline_SKCmodifications.py" "$CONFIG_FILE"

echo "TnSeq pipeline finished. Outputs are in: ${RUN_DIR}/output_files"
