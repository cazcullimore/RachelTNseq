#!/usr/bin/env bash

# Create a small, valid FASTQ.GZ subset and run preprocess.sh on it.

set -euo pipefail

if [[ $# -lt 1 || $# -gt 3 ]]; then
    echo "Usage:"
    echo "  ./test_preprocess.sh INPUT.fastq.gz [NUMBER_OF_READS] [OUTPUT.fastq.gz]"
    echo
    echo "Example:"
    echo "  ./test_preprocess.sh sample.fastq.gz 1000"
    exit 1
fi

INPUT="$1"
READ_COUNT="${2:-1000}"

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: Input file not found: $INPUT"
    exit 1
fi

if [[ "$INPUT" != *.fastq.gz && "$INPUT" != *.fq.gz ]]; then
    echo "ERROR: Input must end in .fastq.gz or .fq.gz"
    exit 1
fi

if [[ ! "$READ_COUNT" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: NUMBER_OF_READS must be a positive integer."
    exit 1
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
BASENAME=$(basename "$INPUT")
SAMPLE="${BASENAME%.fastq.gz}"
SAMPLE="${SAMPLE%.fq.gz}"
OUTPUT="${3:-${SCRIPT_DIR}/test_inputs/${SAMPLE}_first_${READ_COUNT}_reads.fastq.gz}"

if [[ "$OUTPUT" != *.fastq.gz && "$OUTPUT" != *.fq.gz ]]; then
    echo "ERROR: Output must end in .fastq.gz or .fq.gz"
    exit 1
fi

if [[ -e "$OUTPUT" ]]; then
    echo "ERROR: Output already exists: $OUTPUT"
    echo "Remove it or provide a different output path."
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT")"

TEMP_OUTPUT=$(mktemp "${OUTPUT}.tmp.XXXXXX")
trap 'rm -f "$TEMP_OUTPUT"' EXIT

echo "Extracting up to $READ_COUNT reads from: $INPUT"

python3 - "$INPUT" "$TEMP_OUTPUT" "$READ_COUNT" <<'PY'
import gzip
import sys

source, destination, requested_text = sys.argv[1:]
requested = int(requested_text)
written = 0

with gzip.open(source, "rb") as input_handle, gzip.open(destination, "wb") as output_handle:
    for record_number in range(1, requested + 1):
        record = [input_handle.readline() for _ in range(4)]

        if not record[0]:
            break
        if any(not line for line in record):
            raise SystemExit(
                f"ERROR: Incomplete FASTQ record {record_number} in {source}"
            )
        if not record[0].startswith(b"@") or not record[2].startswith(b"+"):
            raise SystemExit(
                f"ERROR: Invalid FASTQ structure at record {record_number} in {source}"
            )
        if len(record[1].rstrip(b"\r\n")) != len(record[3].rstrip(b"\r\n")):
            raise SystemExit(
                f"ERROR: Sequence and quality lengths differ at record {record_number}"
            )

        output_handle.writelines(record)
        written += 1

if written == 0:
    raise SystemExit(f"ERROR: No FASTQ reads found in {source}")

print(f"Created a subset containing {written} read(s).")
PY

mv "$TEMP_OUTPUT" "$OUTPUT"
trap - EXIT

echo "Subset: $OUTPUT"
echo "Starting preprocessing..."

"${SCRIPT_DIR}/preprocess.sh" "$OUTPUT"
