# A succinct single-cell RNA-seq analysis workflow

A practical, end-to-end example workflow for **10x Genomics single-cell RNA-seq** data, from Cell Ranger quantification through downstream analysis with **Seurat v5** and **Harmony**.

The repository is intended as a concise, reusable analysis reference rather than a turnkey clinical or production pipeline. It demonstrates common analysis steps using public PBMC example data and is designed so that paths, scheduler settings, references, and biological decisions can be adapted to each project.

## Workflow overview

```text
Public or user-provided 10x FASTQs
          |
          v
   Cell Ranger count
          |
          v
 Seurat object creation
          |
          v
       QC/filtering
          |
          v
  normalization + merge
          |
          v
   Harmony integration
          |
          v
 PCA / UMAP / clustering
          |
          v
 marker identification
          |
          v
 module scores / annotation
          |
          v
 MAST differential expression
          |
          v
    volcano visualization
```

## What is included

- `analysis_scripts/00-download-public-example-data.md` — download instructions for public 10x PBMC examples and a Cell Ranger reference.
- `analysis_scripts/01-cellranger-count-lsf.sh` — parameterized example Cell Ranger `count` job for an LSF HPC environment.
- `analysis_scripts/02-seurat-analysis.R` — downstream Seurat/Harmony analysis, including QC, integration, clustering, markers, module scoring, example annotation, MAST differential expression, and volcano plotting.
- `analysis_scripts/sessionInfo.txt` — recorded R/package environment for reproducibility of the example downstream analysis.
- `config/samples.example.tsv` — sanitized example sample manifest.

## Requirements

### Cell Ranger preprocessing

The example job script was written for an **LSF-based HPC environment**. Users must provide their own Cell Ranger installation/module and site-specific scheduler settings.

The example script defaults to `cellranger/8.0.1`, but the exact module name is site dependent. Update the LSF project/account, queue, CPU/memory/runtime requests, reference path, and FASTQ paths for your environment.

### R analysis

The downstream script uses packages including:

- Seurat
- harmony
- clustree
- tidyverse
- MAST
- EnhancedVolcano
- ggplot2 and related plotting/data-manipulation packages

See `analysis_scripts/sessionInfo.txt` for the recorded package versions from one successful analysis environment. Package versions evolve, so users should review compatibility when reproducing the workflow in a newer environment.

## Quick start

### 1. Obtain example data or prepare your own data

See:

```text
analysis_scripts/00-download-public-example-data.md
```

The included examples use publicly available 10x Genomics PBMC datasets. No research participant or patient data are distributed in this repository.

### 2. Configure the sample table

Copy the example manifest and replace the paths with your own:

```bash
cp config/samples.example.tsv config/samples.tsv
```

Format:

```text
sample    id       path
sampleA   sample1  /path/to/fastqs
```

For Cell Ranger, `sample` corresponds to the FASTQ sample prefix(es), `id` is the output ID, and `path` is the FASTQ directory.

### 3. Run Cell Ranger on LSF

Set site/project-specific values before submission:

```bash
export PROJECT_DIR="$PWD"
export REFERENCE="/path/to/refdata-gex-GRCh38-2024-A"
export SAMPLE_TABLE="$PWD/config/samples.tsv"
export OUTPUT_DIR="$PWD/results/cellranger"
export CELLRANGER_MODULE="cellranger/8.0.1"
```

Edit these LSF directives in `analysis_scripts/01-cellranger-count-lsf.sh` for your site:

```bash
#BSUB -P YOUR_LSF_PROJECT
#BSUB -q YOUR_LSF_QUEUE
```

Then submit rows from the sample table, for example:

```bash
for num in $(seq 2 "$(wc -l < "$SAMPLE_TABLE")"); do
  export num
  bsub < analysis_scripts/01-cellranger-count-lsf.sh
done
```

### 4. Perform downstream analysis

Use `analysis_scripts/02-seurat-analysis.R` as a worked analysis template. The script intentionally exposes analysis decisions (QC thresholds, dimensionality choices, clustering resolutions, markers, and annotation logic) so they can be reviewed and adapted rather than treated as universal defaults.

## Analysis scope

The downstream R workflow covers, in broad order:

1. loading Cell Ranger feature-barcode matrices;
2. creating Seurat objects and calculating QC metrics;
3. filtering cells and inspecting QC distributions;
4. normalization and merging datasets;
5. highly variable feature selection and scaling;
6. PCA and Harmony batch integration;
7. neighborhood graph construction, clustering at multiple resolutions, and UMAP;
8. cluster marker identification;
9. gene-set/module scoring;
10. example cell-type annotation;
11. differential expression with MAST; and
12. volcano-plot visualization.

These steps are examples, not biological guarantees. QC thresholds, covariates, integration strategy, cluster resolution, annotation markers, and statistical comparisons should be chosen for the scientific question and dataset at hand.

## Reproducibility and portability

The repository separates site-specific paths from the analysis code. Public examples use placeholder paths, and generated data/results are excluded by `.gitignore`.

The Cell Ranger example is HPC/LSF-oriented. The recorded downstream R session was run in a Windows environment, illustrating that the R/Seurat portion is not tied to the HPC scheduler used for preprocessing.

## Data and privacy

This public repository does **not** contain FASTQ files, patient/research participant data, internal sample identifiers, institutional project/account numbers, credentials, or internal filesystem paths. The example datasets referenced in the documentation are public third-party datasets.

Users are responsible for ensuring that their own input data and generated outputs are handled according to applicable institutional, ethical, privacy, and data-governance requirements.

## Third-party software and data

This repository orchestrates or calls third-party tools and R packages; it does not redistribute Cell Ranger, Seurat, Harmony, MAST, or the public datasets themselves. Users should obtain those resources from their official providers and comply with their respective licenses and terms.

When publishing results, please cite the underlying software and datasets as appropriate in addition to citing this repository.

## Citation

If this workflow contributes to your analysis, please cite this repository and acknowledge the workflow development by **Jia-Hua (George) Qu**. GitHub can also display citation metadata from `CITATION.cff`.

Suggested repository citation before a DOI is assigned:

> Qu, J.-H. (George). *A succinct single-cell RNA-seq analysis workflow with Cell Ranger, Seurat, and Harmony*. Version 1.0.0. GitHub repository. https://github.com/JiahuaQu/A-succinct-bioinformatic-analysis-pipeline-of-single-cell-RNA-sequencing-data

Please also cite Cell Ranger/10x Genomics resources and the relevant R packages according to their authors' recommendations.

## Rights and licensing status

No open-source license is currently granted by this repository. See [`COPYRIGHT.md`](COPYRIGHT.md) for the current conservative rights notice. This can be replaced with an explicit software license later if/when the appropriate ownership and licensing framework is confirmed.

## Version

**v1.0.0** — publication-ready restructuring of the original concise scRNA-seq workflow, with sanitized paths, parameterized HPC configuration, expanded documentation, and reproducibility/citation metadata. The core downstream analysis logic is retained from the original repository.
