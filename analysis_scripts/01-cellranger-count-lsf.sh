#!/usr/bin/env bash
# Example LSF job script for Cell Ranger count.
# Customize project/account, queue, resources, module name, and paths for your HPC site.
#BSUB -P YOUR_LSF_PROJECT
#BSUB -J cellranger_count
#BSUB -q YOUR_LSF_QUEUE
#BSUB -n 20
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -oo logs/cellranger_%J.log
#BSUB -eo logs/cellranger_%J.err
#BSUB -W 10:00

set -euo pipefail

module load "${CELLRANGER_MODULE:-cellranger/8.0.1}"

PROJECT_DIR="${PROJECT_DIR:-$PWD}"
REFERENCE="${REFERENCE:-/path/to/refdata-gex-GRCh38-2024-A}"
SAMPLE_TABLE="${SAMPLE_TABLE:-${PROJECT_DIR}/config/samples.example.tsv}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_DIR}/results/cellranger}"

: "${num:?Set num to the 1-based data-row number before submitting the job (for example, export num=2).}"

[[ -d "$PROJECT_DIR" ]] || { echo "ERROR: PROJECT_DIR not found: $PROJECT_DIR" >&2; exit 1; }
[[ -d "$REFERENCE" ]] || { echo "ERROR: REFERENCE not found: $REFERENCE" >&2; exit 1; }
[[ -f "$SAMPLE_TABLE" ]] || { echo "ERROR: SAMPLE_TABLE not found: $SAMPLE_TABLE" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR" "${PROJECT_DIR}/logs"
cd "$PROJECT_DIR"

line=$(sed -n "${num}p" "$SAMPLE_TABLE")
[[ -n "$line" ]] || { echo "ERROR: Line $num not found in $SAMPLE_TABLE" >&2; exit 1; }

IFS=$'\t' read -r sample id data <<< "$line"
[[ "$sample" != "sample" ]] || { echo "ERROR: num points to the header row; use num >= 2" >&2; exit 1; }
[[ -d "$data" ]] || { echo "ERROR: FASTQ directory not found: $data" >&2; exit 1; }

cellranger count \
  --transcriptome="$REFERENCE" \
  --fastqs="$data" \
  --sample="$sample" \
  --id="$id" \
  --nosecondary \
  --create-bam=false \
  --output-dir="$OUTPUT_DIR"

echo "Cell Ranger count completed for sample '$sample' (ID: $id)."
