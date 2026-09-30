# Download public example data

The workflow can be demonstrated with public 10x Genomics PBMC datasets. Run downloads in a directory of your choice rather than at filesystem root.

```bash
mkdir -p datasets reference results/cellranger logs

# 1k PBMC 3' Gene Expression FASTQs
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar -P datasets/
tar -xf datasets/pbmc_1k_v3_fastqs.tar -C datasets/

# Cell Ranger human reference
curl -L -o reference/refdata-gex-GRCh38-2024-A.tar.gz \
  https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
md5sum reference/refdata-gex-GRCh38-2024-A.tar.gz
# Expected MD5 reported when this example was prepared: a7b5b7ceefe10e435719edc1a8b8b2fa

tar -xzf reference/refdata-gex-GRCh38-2024-A.tar.gz -C reference/

# 10k Human PBMC TotalSeqC 5' GEM-X Multiplex FASTQs
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/8.0.0/10k_Human_PBMC_TotalSeqC_5p_gemx_Multiplex/10k_Human_PBMC_TotalSeqC_5p_gemx_Multiplex_fastqs.tar -P datasets/
md5sum datasets/10k_Human_PBMC_TotalSeqC_5p_gemx_Multiplex_fastqs.tar
# Expected MD5 reported when this example was prepared: bb5c0902ccd2850566f95ccda6a897e4

tar -xf datasets/10k_Human_PBMC_TotalSeqC_5p_gemx_Multiplex_fastqs.tar -C datasets/
```

Update `config/samples.example.tsv` (or make your own sample table) with the local FASTQ paths.

For an LSF system, set the required paths/resources and submit one row at a time, for example:

```bash
export PROJECT_DIR="$PWD"
export REFERENCE="$PWD/reference/refdata-gex-GRCh38-2024-A"
export SAMPLE_TABLE="$PWD/config/samples.example.tsv"
export OUTPUT_DIR="$PWD/results/cellranger"

for num in $(seq 2 "$(wc -l < "$SAMPLE_TABLE")"); do
  export num
  bsub < analysis_scripts/01-cellranger-count-lsf.sh
done
```

The URLs above point to third-party public resources. Availability, terms, and versioning are controlled by their respective providers.
