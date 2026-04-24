# Placental Gene-Environment Interaction xQTL Study

A multi-ancestry investigation of how adverse gestational conditions moderate fetal genetic effects on placental isoform regulation and their contribution to childhood metabolic trait heritability.

## Overview

This repository contains analytical code and documentation for a placental expression quantitative trait loci (xQTL) study examining gene-environment (G×E) interactions. We investigate whether gestational diabetes mellitus (GDM) and pre-pregnancy obesity moderate fetal genetic effects on placental gene expression across 2,126 placentas from 9 cohorts spanning 5 ancestry groups.

## Study Design

- **Sample Size**: N = 2,126 placentas
- **Ancestry Groups**: White (N=1,425), East Asian (N=342), Black (N=238), Hispanic/Latino (N=67), South Asian (N=54)
- **Cohorts**: 9 cohorts across multiple countries
- **Exposures**: Gestational diabetes mellitus (GDM, ~20% prevalence), pre-pregnancy BMI
- **Molecular Data**: Bulk RNA-seq quantified against placenta-specific long-read transcriptome assembly
- **RNA Modalities**: Gene expression, isoform expression, splicing, alternative TSS/polyA

## Key Features

- **Multi-ancestry xQTL mapping** across 5 ancestry groups
- **Isoform-resolution analysis** using tissue-specific long-read reference transcriptome
- **Gene-environment interaction testing** for GDM and pre-pregnancy obesity
- **Cell-type deconvolution** to attribute effects to trophoblast, endothelial, and stromal compartments
- **GWAS colocalization** with birth weight and childhood metabolic traits
- **Mediated heritability decomposition** to quantify GxE contributions

## Repository Contents

### Power Analysis
- [power_analysis](https://seantbresnahan.com/Placenta_xQTL_atlas/power_analysis.html) - Comprehensive power calculations for G×E interaction detection
  - Analytical power framework using non-central F-distribution
  - Power curves by minor allele frequency
  - Ancestry-stratified power estimates
  - Minimum detectable effect sizes
  - Binary vs continuous exposure modeling comparisons

## Aims

### Aim 1: Multi-ancestry xQTL Atlas
Build a comprehensive placental xQTL atlas spanning multiple RNA modalities, perform cell-type attribution, and identify placental isoforms underlying GWAS loci through colocalization and isoform-resolution TWAS.

### Aim 2: G×E Interaction Mapping
Map SNP-by-exposure interaction xQTLs for GDM and pre-pregnancy obesity, and estimate their mediated contribution to childhood metabolic trait heritability using MESC and individual-level mediation analysis.

## Data Availability

Data from this study will be made available upon publication through:
- **xQTL summary statistics**: Zenodo
- **Interactive web portal**: [URL to be added]
- **isoTWAS weights**: GitHub release
- **Individual-level data**: dbGaP (controlled access)

## Funding

This work is supported by NIH R21 grant [grant number].

## Citation

If you use data or methods from this study, please cite:

> [Citation to be added upon publication]

## Contact

For questions about this repository or the study:
- **Principal Investigators**: [Names and emails]
- **Issues**: Please use GitHub Issues for code-related questions

## License

Code in this repository is licensed under MIT License. See [LICENSE](LICENSE) for details.

---

**Status**: Active development  
**Last Updated**: April 2026
