# Placenta xQTL atlas

Multi-ancestry, multi-modality cis-xQTL mapping of the human placenta on an
LSF cluster (Bhattacharya Lab, MD Anderson; NIH R21).

Citation: Bresnahan ST, Bhattacharya A, et al. A multi-ancestry xQTL atlas of
the human placenta. In preparation.

The pipeline spans raw sequence and array data through cis-xQTL mapping and
summary reporting: per-cohort genotype processing and TOPMed imputation,
cross-cohort ancestry-stratified genotype pooling, RNA-seq quantification,
PANTRY-style phenotyping of 8 RNA modalities (gene expression, isoform
expression, splice junction usage, intron retention, alternative TSS,
alternative polyA, RNA editing, RNA stability), RNA-to-DNA sample linkage,
ancestry-stratified cis-xQTL mapping with
[tensorQTL](https://github.com/broadinstitute/tensorqtl) under GTEx
conventions, and cross-ancestry fine-mapping with SuSHiE.

Phenotypes are quantified against HPLRv2, the long-read placental
transcriptome annotation from
[sbresnahan/lr-placenta-transcriptomics](https://github.com/sbresnahan/lr-placenta-transcriptomics)
([bioRxiv 2025.06.26.661362](https://doi.org/10.1101/2025.06.26.661362)).

## Modules

- **01_genotype_imputation/**: per-cohort genotype retrieval (array + WXS),
  liftover to GRCh38, pre-imputation filtering, TOPMed imputation, Rsq QC,
  PC-AiR ancestry assignment against 1000 Genomes, and cross-cohort pooling
  into per-ancestry plink2 files (`{ANC}_pooled.pgen`).
- **02_rnaseq_quantification/**: SRA retrieval, adapter trimming, and Salmon
  quantification wrappers for cohorts fetched or quantified outside the
  phenotyping driver.
- **03_phenotyping/**: rewrite of [PejLab/Pantry](https://github.com/PejLab/Pantry)
  for LSF. Shared reference prep, per-cohort STAR alignment and per-modality
  quantification/aggregation into unnormalized BEDs, cross-cohort
  splicing/intron-retention coordinate harmonization, and MuSiC cell-type
  deconvolution.
- **04_sample_linkage/**: builds the per-cohort RNA-seq to genotype sample ID
  map (`rnaseq_to_array_id_map.csv`) consumed by the mapping stage.
- **05_qtl_mapping/**: cross-cohort pooling with QN + within-cohort INT, HCP
  latent-factor estimation and per-modality k optimization, genotype PCA,
  sample intersection, outlier exclusion, covariate assembly, tensorQTL cis
  mapping (grouped, ungrouped, combined, stepwise-conditional) with Storey
  q-values, SuSHiE cross-ancestry fine-mapping, and functional enrichment of
  fine-mapped variants.
- **06_colocalization_twas/**: Objective 1.6 — GWAS catalog +
  fetch/harmonize, genome-wide nominal cis stats, SuSiE-coloc at
  fine-mapped loci with ancestry-matched LD, colocBoost multi-trait
  colocalization, isoTWAS weight training (multivariate elastic net across
  isoforms + per-gene expression), FUSION TWAS against harmonized GWAS, and
  cross-trait aggregation with gene-level ACAT combination and cell-type
  annotation.
- **07_gxe_mapping/**: Objective 2.1 — SNP × exposure (GxE) interaction
  mapping. Pooled multi-ancestry tier-1 discovery (tensorQTL-mirroring
  interaction model, adaptive phenotype permutations, Storey q-values),
  ancestry-stratified tier 2 at Aim-1 prioritized loci with IVW
  meta-analysis, SNP × covariate sensitivity refits (config-gated), and
  transmitted/non-transmitted decomposition (status-gated on maternal
  genotypes).
- **reports/**: analysis reports (`report_placenta_xqtl.Rmd`,
  `report_finemap.Rmd`) rendered to self-contained HTML.

## Workflow overview

```
01_genotype_imputation (per cohort, then pooled)
  retrieval -> liftover -> pre-imputation filtering -> TOPMed imputation
  -> Rsq QC -> PC-AiR ancestry assignment -> {ANC}_pooled.pgen
        |
02_rnaseq_quantification        03_phenotyping (per cohort)
  SRA fetch / trim / Salmon  ->   001-003  reference prep (once, shared)
                                  01-06    per-sample processing
                                  10-15    per-modality aggregation
                                  16       index outputs
                                  17       cross-cohort splicing/IR
                                           harmonization (per ancestry)
                                  18       cell-type deconvolution
        |
04_sample_linkage
  rnaseq_to_array_id_map.csv
        |
05_qtl_mapping (per ancestry)
  19/20   pool + QN + within-cohort INT + HCP factors
  21      tensorqtl env setup (once)
  22-24   genotype PCA, sample intersection, outlier exclusion
  25      covariate assembly
  25a/25b HCP k optimization
  26/30   modality harmonization, combined BED
  27/28   tensorQTL mapping (one LSF job per ancestry x modality)
  29      top tables
  31-33   SuSHiE cross-ancestry fine-mapping + aggregation
  34      annotation enrichment of high-PIP variants
        |
06_colocalization_twas (Objective 1.6)
  36-38   coloc env, GWAS catalog fetch + harmonization
  39/40   genome-wide nominal stats, coloc locus tasks
  41-45   SuSiE-coloc + colocBoost (sharded LSF arrays)
  46-48   isoTWAS weight training + FUSION TWAS
  49      aggregation (coloc/TWAS tables, ACAT, cell-type annotation)
        |
07_gxe_mapping (Objective 2.1)
  50/51   pooled multi-ancestry inputs + pooled HCPs
  52/53   tier-1 GxE scan (chromosome arrays) + merge + q-values
  54      tier-2 ancestry-stratified + IVW meta-analysis
  55/56   SNP x covariate sensitivity (gated), T/NT decomposition (gated)
  57      aggregation
        |
reports/
  summary report (Rmd -> HTML)
```

## Configuration

Stages 3 and 5 run from a git clone of this repository on seadragon and share
one `03_phenotyping/config.yml` (shared reference paths plus one block per
cohort), read via `config_get.py`. Data and reference paths are
config-driven; infrastructure paths (conda envs, singularity image, R
library, `/rsrch5`/`/rsrch9` roots) are hardcoded in the scripts.

## Documentation

- 01_genotype_imputation/README.md: imputation, ancestry assignment, pooling
- 02_rnaseq_quantification/README.md: RNA-seq fetch/trim/quant
- 03_phenotyping/README.md: phenotyping walkthrough with commands
- 04_sample_linkage/README.md: sample ID linkage
- 05_qtl_mapping/README.md: QTL mapping and fine-mapping walkthrough with
  commands
- reports/README.md: report contents and re-rendering
- docs/data_availability.md: cohorts, accessions, sample counts, references
- docs/software_environments.md: tool, environment, and package versions
- docs/porting.md: adapting the pipeline to a different cluster
- docs/gxe_power_analysis.html: power calculations for SNP x exposure scans
- MANIFEST.md: file inventory

## Quick start

Each module README documents the full order of operations for its stage.
In brief:

```bash
# Phenotyping (per cohort; see 03_phenotyping/README.md):
cd 03_phenotyping
bsub -env "CONFIG=$PWD/config.yml,SCRIPTS_DIR=$PWD" < 001_refprep_pre.sh   # then 002, 003
bsub -env "CONFIG=$PWD/config.yml,SCRIPTS_DIR=$PWD,COHORT=cohort1" < 00_run_pipeline.sh

# QTL mapping (per ancestry; see 05_qtl_mapping/README.md):
cd ../05_qtl_mapping
bsub -env "CONFIG=$PWD/../03_phenotyping/config.yml,SCRIPTS_DIR=$PWD" < 19_hcp_factors.sh   # then 20-28
```

A standalone, config-driven re-implementation of modules 03 and 05 (scripts
19-30) for general use is at
[sbresnahan/pantry-seadragon](https://github.com/sbresnahan/pantry-seadragon).
This repository is the record of the placenta project as run.

## Data availability

Cohort data are controlled-access or available on request; see
[docs/data_availability.md](docs/data_availability.md) for accession-level
detail (SNUH: SRA PRJNA820329; NIGMS: SRA PRJNA671171; NIEHS_RICHS: dbGaP;
GUSTO: study investigators). No sample-level data are distributed in this
repository.

## Computational environment

All stages run on the MD Anderson seadragon HPC cluster (LSF scheduler,
RHEL8, conda via miniforge3, R 4.3.1 via singularity). Tool and package
versions: [docs/software_environments.md](docs/software_environments.md).
Porting notes (LSF to SLURM, path replacement):
[docs/porting.md](docs/porting.md).

## License

GPL v3 (see [LICENSE](LICENSE)). The bundled MAJIQ scripts require an
academic MAJIQ license, which is not distributed here; obtain one from the
[MAJIQ project](https://majiq.biociphers.org/) and point `majiq_license` in
`03_phenotyping/config.yml` at it.

## Contact

Sean T. Bresnahan — stbresnahan@mdanderson.org
Bhattacharya Lab for Computational Genomics
