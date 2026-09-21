# Software environments

All versions below are the versions used on the MD Anderson seadragon cluster for
the analyses in this repository (RHEL8; conda via miniforge3 24.5.0-0; environment
modules).

## Command-line tools

| Tool | Version | Activation on seadragon | Used in |
|---|---|---|---|
| STAR | 2.7.4a | `conda activate star-2.7.4a` | 03 phenotyping (alignment) |
| Salmon | 1.10.2 | `conda activate salmon-1.10.2` | 02 quantification, 03 phenotyping |
| samtools / htslib (bgzip, tabix) | 1.16.1 | `conda activate samtools-1.16.1` | 03, 05 |
| Trim Galore | 0.6.10 | `conda activate trim-galore-0.6.10` | 02 |
| Subread (featureCounts) | module | `module load subread` | 03 (stability) |
| RegTools | direct binary | — | 03 (splice junctions) |
| leafCutter | venv | `source .../leafcutter/bin/activate` | 03 (splicing clusters) |
| MAJIQ | venv, academic license | `source .../MAJIQ/bin/activate` | 03 (intron retention) |
| gffread | venv | `source .../gffread/bin/activate` | 03 (reference prep) |
| plink / plink2 | module | `module load plink` | 01, 05 |
| bcftools | module | `module load bcftools` | 01 |
| BWA | module | — | 01 (NIGMS WXS) |
| GATK | module | — | 01 (gVCF merge) |
| sratoolkit | module | `module add sratoolkit` | 01, 02 |
| TOPMed imputation server | — | off-cluster web service | 01 |

## Python

- **Phenotyping stages (03)**: Python 3 from the MAJIQ env stack
  (`samtools-1.16.1` conda env + MAJIQ venv): pandas, numpy, scipy, gtfparse,
  matplotlib, PyYAML.
- **Mapping stages (05)**: dedicated python-only conda env `tensorqtl`:
  **tensorqtl 1.0.10**, pandas, numpy, pyarrow, genotypeio, matplotlib.
  Installed by `05_qtl_mapping/21_install_tensorqtl.sh`. The env has **no working
  R/rpy2** — Storey q-values go through the singularity-R bridge
  (`compute_qvalues.R`, `QVALUE_RSCRIPT`).

## R

R **4.3.1** via singularity
(`singularity exec --bind /rsrch5 --bind /rsrch9 .../rstudio_4.3.1.sif Rscript`),
packages from `R_LIBS_USER=/rsrch5/.../R_package_library/ubuntu/4.3.1`.

| Package | Used for |
|---|---|
| data.table, dplyr, tidyr, optparse, parallel | plumbing across stages |
| ggplot2, patchwork, RColorBrewer | QC/diagnostic plots |
| sva (ComBat) | cross-cohort batch correction (03, 05) |
| Rhcpp | HCP latent-factor estimation (05) |
| edgeR, tximport | `catchSalmon` QU correction of isoforms (03) |
| rtracklayer | GTF parsing / tx2gene maps (03) |
| MuSiC, SingleCellExperiment, Matrix | cell-type deconvolution (03) |
| qvalue | Storey q-values (05) |
| GENESIS, GWASTools, SNPRelate | PC-AiR ancestry assignment (01) |
| vcfR | VCF QC (01) |
| rmarkdown, knitr + tidyverse | report rendering (reports/) |

## Reference genomes/annotations

GRCh38 no-alt analysis set (GCA_000001405.15); HPLRv2 long-read placental
transcriptome (SQANTI3-annotated; see
[lr-placenta-transcriptomics](https://github.com/sbresnahan/lr-placenta-transcriptomics));
GENCODE v45 (deduplicated); REDIportal TABLE1_hg38_v3; 1000 Genomes phase 3
(ancestry reference). Paths are configured in `03_phenotyping/config.yml`.
