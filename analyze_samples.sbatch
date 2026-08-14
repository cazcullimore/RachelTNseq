#!/usr/bin/env bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

set -eo pipefail

module load miniforge3
conda activate Tnseq
module load bowtie2
set -u

if (( $# != 7 )); then
    echo "Usage: sbatch analyze_samples.sbatch REFERENCE.fasta[.gz] ANNOTATION.gff[.gz] PREPROCESS_DIR ANALYSIS_DIR OUTPUT_NAME EXPECTED_SAMPLES SCRIPT_DIR" >&2
    exit 2
fi

REFERENCE=$1
GFF=$2
PREPROCESS_DIR=$3
ANALYSIS_DIR=$4
OUTPUT_NAME=$5
EXPECTED_SAMPLES=$6
SCRIPT_DIR=$7

for input_path in "$REFERENCE" "$GFF"; do
    if [[ ! -f "$input_path" ]]; then
        echo "ERROR: Input file not found: $input_path" >&2
        exit 1
    fi
done
if [[ ! -d "$PREPROCESS_DIR" ]]; then
    echo "ERROR: Preprocessing directory not found: $PREPROCESS_DIR" >&2
    exit 1
fi
if [[ ! "$EXPECTED_SAMPLES" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: EXPECTED_SAMPLES must be a positive integer." >&2
    exit 1
fi

shopt -s nullglob
FINAL_READS=("${PREPROCESS_DIR}"/*/final/*_filtered.fastq.gz)
if (( ${#FINAL_READS[@]} != EXPECTED_SAMPLES )); then
    echo "ERROR: Expected $EXPECTED_SAMPLES filtered FASTQ files, found ${#FINAL_READS[@]}." >&2
    exit 1
fi

# The submitter reserves the batch root, while this atomic mkdir ensures that
# an analysis is never restarted over an existing result set.
if ! mkdir -- "$ANALYSIS_DIR"; then
    echo "ERROR: Analysis directory already exists or could not be created: $ANALYSIS_DIR" >&2
    exit 1
fi

PTT_DIR="${ANALYSIS_DIR}/ptt"
BOWTIE_DIR="${ANALYSIS_DIR}/bowtie_reference"
CONFIG_FILE="${ANALYSIS_DIR}/tnSeq.config"
INDEX_PREFIX="${BOWTIE_DIR}/combined_reference"
mkdir -- "$PTT_DIR" "$BOWTIE_DIR"

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

echo "Building the shared Bowtie 2 reference..."
bowtie2-build --threads "${SLURM_CPUS_PER_TASK:-8}" "$REFERENCE_FOR_TOOLS" "$INDEX_PREFIX"

echo "Creating shared legacy PTT annotation files..."
python3 "${SCRIPT_DIR}/gff_to_ptt.py" \
    --gff "$GFF" \
    --fasta "$REFERENCE_FOR_TOOLS" \
    --outdir "$PTT_DIR"

# The glob deliberately includes every filtered sample in this batch. The
# legacy pipeline treats these files as conditions in one normalized table.
# trying to make read length longer and not do a preemptive trim
cat > "$CONFIG_FILE" <<EOF
[input]
Reads           = ${PREPROCESS_DIR}/*/final/*_filtered.fastq.gz
BowtieReference = $INDEX_PREFIX
Ptt             = $PTT_DIR/*.ptt
Out             = $OUTPUT_NAME

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

echo "Running one combined analysis for ${#FINAL_READS[@]} samples..."
python2.7 "${SCRIPT_DIR}/TnSeq-Pipeline_SKCmodifications.py" "$CONFIG_FILE"

echo "Combined TnSeq analysis finished: ${ANALYSIS_DIR}/output_files"
