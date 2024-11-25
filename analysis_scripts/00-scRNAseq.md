```bash
bsub -P scRNAseq -J pipeline -n 1 -q interactive -R "rusage[mem=8001]" -W 8:00 -Is "bash"

cd /home/jqu/project/scRNAseq_pipeline

# Download sequencing data
# https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct
mkdir datasets
cd datasets

nohup wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar &

tar -xvf pbmc_1k_v3_fastqs.tar


### Make table to record all samples and their paths
# 因为是从excel另存为的txt，所以行末会有 ^M，UNIX不能正确识别，所以需要去除
cat -A sample.txt
# pbmc_1k_v3^Isample1^I/home/jqu/project/scRNAseq_pipeline/datasets/pbmc_1k_v3_fastqs^M$

sed -i 's/[[:space:]]*$//' sample.txt
cat -A sample.txt
# pbmc_1k_v3^Isample1^I/home/jqu/project/scRNAseq_pipeline/datasets/pbmc_1k_v3_fastqs$


# Reference genome
# https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads

# Human genome for this time
cd /research/sharedresources/immunoinformatics/common/jqu/reference_genome/cellranger/human

nohup curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz" &

# Verify
md5sum refdata-gex-GRCh38-2024-A.tar.gz
# File size: 11 GB
# md5sum: a7b5b7ceefe10e435719edc1a8b8b2fa
# a7b5b7ceefe10e435719edc1a8b8b2fa  refdata-gex-GRCh38-2024-A.tar.gz
nohup tar -xvzf refdata-gex-GRCh38-2024-A.tar.gz &


# Mouse genome for the future use, not this time
cd /research/sharedresources/immunoinformatics/common/jqu/reference_genome/cellranger/mouse

nohup curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz" &

# Verify
md5sum refdata-gex-GRCm39-2024-A.tar.gz
# File size: 9.7 GB
# md5sum: 37c51137ccaeabd4d151f80dc86ce0b3
nohup tar -xvzf refdata-gex-GRCm39-2024-A.tar.gz &


### https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct
cd /home/jqu/project/scRNAseq_pipeline

mkdir -p ./count/

module avail | grep cellranger

touch 01-cellranger_count.sh

### Do not create the second bash script, but run it in command line directly
# Path to the file listing sample details
sample_text="./sample.txt"

# Get the number of lines in sample.txt
total_lines=$(wc -l < ${sample_text})

# Loop from 2 to total_lines
for num in $(seq 2 ${total_lines})
do
  echo $num
done

for num in $(seq 2 ${total_lines})
do
  # Set the environment variable num and submit the job
  export num
  bsub < 01-cellranger_count.sh
done
```



