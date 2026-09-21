# Multi-ancestry xQTL mapping of the human placenta

This repository contains the complete analysis codebase and documentation for the
multi-ancestry, multi-modal placental xQTL mapping project of the Bhattacharya Lab
(NIH R21; PI: Arjun Bhattacharya):

> **A multi-ancestry xQTL atlas of the human placenta**
> Sean T. Bresnahan, Arjun Bhattacharya, et al.
> *in preparation*

The pipeline spans raw sequence and array data through cis-xQTL mapping and summary
reporting: per-cohort genotype processing and TOPMed imputation, cross-cohort
ancestry-stratified genotype pooling, RNA-seq quantification, PANTRY-style phenotyping
of **8 RNA modalities** (gene expression, isoform expression, splice junction usage,
intron retention, alternative TSS, alternative polyA, RNA editing, RNA stability),
RNA→DNA sample linkage, and ancestry-stratified cis-xQTL mapping with
[tensorQTL](https://github.com/broadinstitute/tensorqtl) under GTEx conventions.

Phenotypes are quantified against **HPLRv2**, the long-read placental transcriptome
annotation from [sbresnahan/lr-placenta-transcriptomics](https://github.com/sbresnahan/lr-placenta-transcriptomics)
([bioRxiv 2025.06.26.661362](https://doi.org/10.1101/2025.06.26.661362)).

---

## The project manifest

📄 [`MANIFEST.md`](MANIFEST.md) is the reproducibility centerpiece: a stage-by-stage
account of the full workflow — raw data → imputation → quantification → phenotyping →
sample linkage → QTL mapping → report — with the scripts, inputs, outputs, software,
and run order for each stage, plus the operational landmines learned across analysis
rounds. **Start here to reproduce or extend the analysis.**

## Repository structure

### I. Genotype imputation and pooling
📁 [`01_genotype_imputation/`](01_genotype_imputation/)
Per-cohort genotype retrieval (array + WXS), liftover to GRCh38, pre-imputation
filtering, TOPMed imputation, Rsq QC, PC-AiR ancestry assignment against 1000 Genomes,
and cross-cohort pooling into per-ancestry plink2 files (`{ANC}_pooled.pgen`).

### II. RNA-seq quantification
📁 [`02_rnaseq_quantification/`](02_rnaseq_quantification/)
SRA retrieval, adapter trimming (Trim Galore), and Salmon quantification wrappers for
the cohorts requiring fetch/quant outside the PANTRY driver.

### III. Transcriptomic phenotyping (PANTRY, seadragon rewrite)
📁 [`03_phenotyping/`](03_phenotyping/)
In-house rewrite of [PejLab/Pantry](https://github.com/PejLab/Pantry) for the MD
Anderson seadragon HPC (LSF): shared reference prep (HPLRv2 GTF normalization,
txrevise event construction, REDIportal edit sites), per-cohort STAR alignment and
per-modality quantification/aggregation into normalized BEDs, cross-cohort splicing /
intron-retention coordinate harmonization, and MuSiC cell-type deconvolution.

### IV. RNA→DNA sample linkage
📁 [`04_sample_linkage/`](04_sample_linkage/)
Builds the verified per-cohort RNA-seq → genotype sample ID map
(`rnaseq_to_array_id_map.csv`) consumed by the QTL mapping stage.

### V. cis-xQTL mapping (tensorQTL, GTEx conventions)
📁 [`05_qtl_mapping/`](05_qtl_mapping/)
HCP latent-factor estimation, cross-cohort ComBat + inverse-normal transformation of
all modalities, cohort-only genotype PCA, genotype×phenotype intersection with a
minor-allele-count floor, PC-outlier exclusion, covariate assembly/optimization,
modality harmonization, and tensorQTL cis mapping (grouped, ungrouped, combined, and
stepwise-conditional modes) with Storey q-values.

### VI. Summary report
📁 [`reports/`](reports/)
Analysis report (`report_placenta_xqtl.Rmd` → self-contained HTML):
mapping yield, calibration diagnostics, power analysis, the FDR-significant
associations, near-miss landscape, and cross-ancestry concordance.

### Supporting documentation
📁 [`docs/`](docs/)

- 📄 `data_availability.md` — cohorts, accessions, ancestry × cohort sample counts, reference data
- 📄 `round_history.md` — design evolution across mapping rounds 1–3 and why conventions changed
- 📄 `software_environments.md` — tool, environment, and package versions
- 📄 `gxe_power_analysis.html` — analytical power calculations for Aim-2 SNP × exposure
  (G×E) scans at the full planned cohort (N = 2,126; GDM binary and pre-pregnancy BMI
  continuous exposures): minimum detectable interaction effects by MAF, ancestry-stratum
  power, and binary-vs-continuous exposure modeling

## Current status

GTEx-conventions mapping is complete for the **EAS** (n ≈ 272) and **EUR**
(n ≈ 136) ancestry strata: 32 ancestry × modality × layer cells, two FDR ≤ 5%
associations (an EAS isoform xQTL at *TSPAN3* and a EUR RNA-editing xQTL at
*NCOA4*), clean calibration, and a power-limited null elsewhere (minimum
detectable r² ≈ 0.10 at EAS n, ≈ 0.19 at EUR n). AFR/AMR/SAS
pooled genotypes exist and are the next mapping targets. See
[`docs/round_history.md`](docs/round_history.md) and the report in
[`reports/`](reports/).

## Data availability

Cohort data are controlled-access or available on request; see
[`docs/data_availability.md`](docs/data_availability.md) for accession-level detail
(SNUH: SRA PRJNA820329; NIGMS: SRA PRJNA671171; NIEHS_RICHS: dbGaP; GUSTO: study
investigators). No sample-level data are distributed in this repository.

## Computational environment

All pipeline stages were developed and run on the MD Anderson seadragon HPC cluster
(LSF scheduler, RHEL8). Scripts contain seadragon-specific paths (`/rsrch5`,
`/rsrch9`), conda environments, and a singularity R image; adapting the pipeline to
another system requires editing these paths (see the porting notes in
[`MANIFEST.md`](MANIFEST.md)). R 4.3.1 package dependencies and tool versions are
listed in [`docs/software_environments.md`](docs/software_environments.md).

## License

GPL v3 (see [`LICENSE`](LICENSE)). Note that the bundled MAJIQ scripts require an
academic MAJIQ license, which is **not** distributed here; obtain one from the
[MAJIQ project](https://majiq.biociphers.org/) and point `majiq_license` in
`03_phenotyping/config.yml` at it.

## Contact

For questions, please contact:

📧 Sean T. Bresnahan — stbresnahan@mdanderson.org

🧬 Bhattacharya Lab for Computational Genomics
