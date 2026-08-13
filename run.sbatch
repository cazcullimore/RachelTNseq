#!/usr/bin/env bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=TnSeqJob

set -eo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

# Environment setup is intentionally left as-is for the cluster owner to adjust.
module load miniforge3
conda activate Tnseq

set -u

usage() {
    cat <<'EOF'
Usage:
  sbatch run.sbatch READS.fastq[.gz] REFERENCE.fasta[.gz] ANNOTATION.gff[.gz] [TRANSPOSON] [OUTPUT_DIR]

Arguments:
  READS.fastq[.gz]  Demultiplexed FASTQ input.
  REFERENCE.fasta[.gz]
                     Reference FASTA used to build the Bowtie 2 index.
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

create_run_directory() {
    local path="$1"
    local parent
    local name

    parent=$(dirname -- "$path")
    name=$(basename -- "$path")
    parent=$(absolute_directory "$parent")
    path="${parent}/${name}"

    # mkdir without -p is an atomic claim on this run directory. This prevents
    # two Slurm jobs from silently writing into the same set of output files.
    if ! mkdir -- "$path"; then
        echo "ERROR: Output directory already exists or could not be created: $path" >&2
        echo "Choose a different OUTPUT_DIR. Existing runs are never overwritten." >&2
        return 1
    fi

    printf '%s\n' "$path"
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
TRANSPOSON=$(printf '%s' "$TRANSPOSON" | tr '[:lower:]' '[:upper:]')

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
RUN_DIR=$(create_run_directory "$OUTPUT_ARG")

PREPROCESS_DIR="${RUN_DIR}/preprocessed"
PTT_DIR="${RUN_DIR}/ptt"
BOWTIE_DIR="${RUN_DIR}/bowtie_reference"
CONFIG_FILE="${RUN_DIR}/tnSeq.config"
INDEX_PREFIX="${BOWTIE_DIR}/${SAMPLE}"
FINAL_READS="${PREPROCESS_DIR}/${SAMPLE}/final/${SAMPLE}_filtered.fastq.gz"

mkdir -p -- "$PREPROCESS_DIR" "$PTT_DIR" "$BOWTIE_DIR"

# Some Bowtie 2 installations do not transparently read compressed references.
# Expand a gzipped FASTA once, then share the temporary plain-text copy between
# bowtie2-build and gff_to_ptt.py.
REFERENCE_FOR_TOOLS="$REFERENCE"
REFERENCE_TMP=""
cleanup() {
    if [[ -n "$REFERENCE_TMP" ]]; then
        rm -f -- "$REFERENCE_TMP"
    fi
}
trap cleanup EXIT

if [[ "$REFERENCE" == *.gz ]]; then
    REFERENCE_TMP=$(mktemp "${BOWTIE_DIR}/reference_fasta.XXXXXX")
    echo "Decompressing reference FASTA..."
    gzip -dc -- "$REFERENCE" > "$REFERENCE_TMP"
    REFERENCE_FOR_TOOLS="$REFERENCE_TMP"
fi

echo "Run directory: $RUN_DIR"
echo "Preprocessing reads..."
RESULTS_DIR="$PREPROCESS_DIR" TRANSPOSON="$TRANSPOSON" \
    bash "${SCRIPT_DIR}/preprocess.sh" "$READS"

if [[ ! -f "$FINAL_READS" ]]; then
    echo "ERROR: Expected preprocessed reads were not created: $FINAL_READS" >&2
    exit 1
fi

echo "Building Bowtie 2 reference..."
bowtie2-build --threads "${SLURM_CPUS_PER_TASK:-8}" "$REFERENCE_FOR_TOOLS" "$INDEX_PREFIX"

echo "Creating legacy PTT annotation files..."
python3 "${SCRIPT_DIR}/gff_to_ptt.py" \
    --gff "$GFF" \
    --fasta "$REFERENCE_FOR_TOOLS" \
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
