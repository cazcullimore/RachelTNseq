#!/usr/bin/env bash

# Submit parallel preprocessing jobs followed by one combined TnSeq analysis.
# commmand bash combined_TnSeq/submit_samples.sh -o full_run WSHY6G_1_M1_1_reference.fasta M1_Giff.gff *fastq.gz
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  ./submit_samples.sh --output-dir OUTPUT_DIR REFERENCE.fasta[.gz] ANNOTATION.gff[.gz] READS.fastq[.gz]...

Example:
  ./submit_samples.sh --output-dir experiment_1 \
    ref.fasta.gz genes.gff.gz reads/*.fastq.gz

Each sample is preprocessed in a separate Slurm job. One combined analysis job
starts only after every preprocessing job succeeds, preserving multi-sample
normalization in the legacy Python pipeline. OUTPUT_DIR must not already exist.
EOF
}

OUTPUT_ROOT=""
while (( $# > 0 )); do
    case "$1" in
        -o|--output-dir)
            if (( $# < 2 )) || [[ -z "$2" ]]; then
                echo "ERROR: $1 requires a directory path." >&2
                exit 2
            fi
            OUTPUT_ROOT=$2
            shift 2
            ;;
        --output-dir=*)
            OUTPUT_ROOT=${1#*=}
            if [[ -z "$OUTPUT_ROOT" ]]; then
                echo "ERROR: --output-dir requires a directory path." >&2
                exit 2
            fi
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "ERROR: Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
        *)
            break
            ;;
    esac
done

if [[ -z "$OUTPUT_ROOT" ]]; then
    echo "ERROR: --output-dir is required." >&2
    usage >&2
    exit 2
fi

if (( $# < 3 )); then
    usage >&2
    exit 2
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
SUBMIT_DIR=$PWD
REFERENCE=$1
GFF=$2
shift 2

[[ "$REFERENCE" = /* ]] || REFERENCE="${SUBMIT_DIR}/${REFERENCE}"
[[ "$GFF" = /* ]] || GFF="${SUBMIT_DIR}/${GFF}"

for common_input in "$REFERENCE" "$GFF"; do
    if [[ ! -f "$common_input" ]]; then
        echo "ERROR: Input file not found: $common_input" >&2
        exit 1
    fi
done

READ_PATHS=()
SAMPLE_NAMES=()
SAMPLE_COUNT=0

# Validate the complete batch before creating directories or submitting jobs.
for reads_arg in "$@"; do
    [[ "$reads_arg" = /* ]] || reads_arg="${SUBMIT_DIR}/${reads_arg}"
    if [[ ! -f "$reads_arg" ]]; then
        echo "ERROR: Reads file not found: $reads_arg" >&2
        exit 1
    fi

    reads_name=$(basename -- "$reads_arg")
    sample="${reads_name%.fastq.gz}"
    sample="${sample%.fq.gz}"
    sample="${sample%.fastq}"
    sample="${sample%.fq}"

    if [[ "$sample" == "$reads_name" ]]; then
        echo "ERROR: Reads must end in .fastq, .fq, .fastq.gz, or .fq.gz: $reads_arg" >&2
        exit 1
    fi
    if [[ ! "$sample" =~ ^[A-Za-z0-9_-]+$ ]]; then
        echo "ERROR: Sample names may contain only letters, numbers, underscores, and hyphens: $sample" >&2
        echo "This prevents filename collisions inside the unmodified legacy Python pipeline." >&2
        exit 1
    fi
    if (( SAMPLE_COUNT > 0 )); then
        for seen_sample in "${SAMPLE_NAMES[@]}"; do
            if [[ "$sample" == "$seen_sample" ]]; then
                echo "ERROR: More than one input resolves to sample name '$sample'." >&2
                echo "Rename one input so every sample has a unique basename." >&2
                exit 1
            fi
        done
    fi

    READ_PATHS+=("$reads_arg")
    SAMPLE_NAMES+=("$sample")
    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
done

[[ "$OUTPUT_ROOT" = /* ]] || OUTPUT_ROOT="${SUBMIT_DIR}/${OUTPUT_ROOT}"
OUTPUT_PARENT=$(dirname -- "$OUTPUT_ROOT")
OUTPUT_NAME=$(basename -- "$OUTPUT_ROOT")
mkdir -p -- "$OUTPUT_PARENT"
OUTPUT_PARENT=$(cd -- "$OUTPUT_PARENT" && printf '%s\n' "$PWD")
OUTPUT_ROOT="${OUTPUT_PARENT}/${OUTPUT_NAME}"

# Claim the complete batch directory atomically. Separate invocations can never
# share preprocessing files, configs, indexes, or final outputs.
if ! mkdir -- "$OUTPUT_ROOT"; then
    echo "ERROR: OUTPUT_ROOT already exists or could not be created: $OUTPUT_ROOT" >&2
    echo "Choose a new OUTPUT_ROOT; existing batches are never overwritten." >&2
    exit 1
fi

PREPROCESS_DIR="${OUTPUT_ROOT}/preprocessed"
ANALYSIS_DIR="${OUTPUT_ROOT}/analysis"
LOG_DIR="${OUTPUT_ROOT}/slurm_logs"
mkdir -- "$PREPROCESS_DIR" "$LOG_DIR"

PREPROCESS_JOB_IDS=()
for index in "${!READ_PATHS[@]}"; do
    reads=${READ_PATHS[$index]}
    sample=${SAMPLE_NAMES[$index]}
    job_label=${sample//[^[:alnum:]_.-]/_}

    echo "Submitting preprocessing: $sample"
    job_id=$(sbatch --parsable \
        --job-name="preprocess_${job_label}" \
        --output="${LOG_DIR}/${job_label}_preprocess_%j.out" \
        --error="${LOG_DIR}/${job_label}_preprocess_%j.err" \
        "${SCRIPT_DIR}/preprocess_sample.sbatch" \
        "$reads" "$PREPROCESS_DIR" "$SCRIPT_DIR")
    job_id=${job_id%%;*}
    if [[ ! "$job_id" =~ ^[0-9]+$ ]]; then
        echo "ERROR: sbatch returned an invalid preprocessing job ID: $job_id" >&2
        exit 1
    fi
    PREPROCESS_JOB_IDS+=("$job_id")
done

dependency=$(IFS=:; printf '%s' "${PREPROCESS_JOB_IDS[*]}")
analysis_label=${OUTPUT_NAME//[^[:alnum:]_.-]/_}

echo "Submitting combined analysis after preprocessing jobs: $dependency"
analysis_job_id=$(sbatch --parsable \
    --dependency="afterok:${dependency}" \
    --job-name="tnseq_${analysis_label}" \
    --output="${LOG_DIR}/combined_analysis_%j.out" \
    --error="${LOG_DIR}/combined_analysis_%j.err" \
    "${SCRIPT_DIR}/analyze_samples.sbatch" \
    "$REFERENCE" "$GFF" "$PREPROCESS_DIR" "$ANALYSIS_DIR" \
    "$OUTPUT_NAME" "$SAMPLE_COUNT" "$SCRIPT_DIR")
analysis_job_id=${analysis_job_id%%;*}

echo "Submitted $SAMPLE_COUNT preprocessing job(s): $dependency"
echo "Combined analysis job: $analysis_job_id"
echo "Batch output: $OUTPUT_ROOT"
echo "The analysis job runs only if every preprocessing job succeeds."
