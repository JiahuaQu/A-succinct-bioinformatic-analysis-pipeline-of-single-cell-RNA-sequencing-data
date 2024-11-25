#!/bin/bash
#BSUB -P scRNAseq
#BSUB -J count
#BSUB -q superdome
#BSUB -n 20
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=10000]"
#BSUB -oo run_%J.log
#BSUB -eo run_%J.error
#BSUB -W 10:00
#BSUB -B 
#BSUB -N 

# Load the cellranger module
module load cellranger/8.0.1

# Set the working directory
cd /home/jqu/project/scRNAseq_pipeline || { echo "Failed to change directory!"; exit 1; }

# Reference genome path
ref="/research/sharedresources/immunoinformatics/common/jqu/reference_genome/cellranger/human/refdata-gex-GRCh38-2024-A"

# Path to the file listing sample details
sample_text="./sample.txt"

# Validate if sample text file exists
if [[ ! -f ${sample_text} ]]; then
  echo "Error: Sample text file '${sample_text}' not found!"
  exit 1
fi

# Get the line corresponding to the current task ID
line=$(sed -n "${num}p" "${sample_text}")
if [[ -z $line ]]; then
  echo "Error: Line $num not found in '${sample_text}'."
  exit 1
fi

# Parse each field from the line safely
sample=$(echo "$line" | awk '{print $1}')   # Original sample name
id=$(echo "$line" | awk '{print $2}')       # Renamed sample id
data=$(echo "$line" | awk '{print $3}')     # Path to FASTQ directory

# Validate if FASTQ directory exists
if [[ ! -d $data ]]; then
  echo "Error: FASTQ directory '$data' not found for sample '$sample'."
  exit 1
fi

# Create output directory if not exists
output_dir="./count/"
mkdir -p "${output_dir}"

# Run cellranger count
cellranger count \
  --transcriptome="$ref" \
  --fastqs="$data" \
  --sample="$sample" \
  --id="$id" \
  --nosecondary \
  --create-bam false \
  --output-dir="$output_dir" || { echo "Cellranger count failed for sample '$sample'."; exit 1; }

# Log success
echo "Cellranger count completed successfully for sample '$sample' (ID: $id)."
