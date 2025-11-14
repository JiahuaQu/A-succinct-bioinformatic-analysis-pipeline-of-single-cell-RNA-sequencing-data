```bash

cd /
mkdir datasets/
cd datasets/

nohup wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar &

tar -xvf pbmc_1k_v3_fastqs.tar

cd / 
mkdir -p ref/ 
 
# Download the appropriate reference
nohup curl -o ./ref/refdata-gex-GRCh38-2024-A.tar.gz "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz" & 
 
# Verify 
md5sum ./ref/refdata-gex-GRCh38-2024-A.tar.gz 
# File size: 11 GB 
# md5sum: a7b5b7ceefe10e435719edc1a8b8b2fa 
 
# Decompress 
tar -xzvf ./ref/refdata-gex-GRCh38-2024-A.tar.gz -C ./ref/

# Create output directory if not exists 
output_dir="./count/" 
mkdir -p "${output_dir}" 
 
# Path to the file listing sample details 
sample_text="./sample.txt" 
 
# Get the number of lines in sample.txt 
total_lines=$(wc -l < ${sample_text}) 
 
# Loop from 2 to total_lines, because the first line contans the column names 
for num in $(seq 2 ${total_lines}) 
do 
  # Set the environment variable num and submit the job 
  export num 
  bsub < cellranger_count.sh 
done

```




