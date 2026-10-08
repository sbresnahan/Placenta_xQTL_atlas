# Data availability

## Cohorts

| Cohort | Data types | Access | Accession / route |
|---|---|---|---|
| NIEHS_RICHS | Placenta RNA-seq + array genotypes | Controlled access | dbGaP *(accession number to be added on publication)* |
| GUSTO | Placenta RNA-seq + array genotypes | Available on request | GUSTO study investigators |
| SNUH | Placenta RNA-seq + genotypes | Public | SRA BioProject **PRJNA820329** (102 RNA-seq runs) |
| NIGMS | Placenta RNA-seq (2 sampling sites/subject) + WXS | Public | SRA BioProject **PRJNA671171** (154 RNA-seq + 95 WXS runs) |
| MALI_G3A | Array genotypes only | Controlled access | Pooled but excluded (no matching RNA) |

No sample-level genotype, phenotype, or covariate data are distributed in this
repository. The sample-linkage script (`04_sample_linkage/`) and all mapping
scripts expect the user to have obtained cohort data through the routes above.

## Pooled genotype sample counts (ancestry × cohort)

Counts from the pooled genotype table (`pooled_sample_ancestry.tsv`) after
imputation QC and PC-AiR ancestry assignment:

| Cohort | AFR | AMR | EAS | EUR | SAS |
|---|---|---|---|---|---|
| NIEHS_RICHS | 16 | 6 | 8 | 119 | – |
| GUSTO | – | – | 886 | – | 187 |
| SNUH | – | – | 102 | – | – |
| NIGMS | 19 | 2 | 5 | 27 | 6 |
| MALI_G3A *(excluded)* | 746 | – | – | – | – |

After genotype×phenotype intersection and outlier exclusion, the mapped
strata carry n ≈ 272 (EAS) and n ≈ 136 (EUR) samples. AFR/AMR/SAS pooled
genotypes exist.

## GWAS summary statistics (module 06)

`06_colocalization_twas/gwas_catalog.tsv` lists the seven Objective-1.6
traits with access routes. Most are public downloads fetched by
`37_fetch_gwas.py`; two require manual placement after registration:

| Trait | Source | Access |
|---|---|---|
| Childhood anthropometrics | Richardson 2020 (UKB recalled body size), GWAS Catalog GCST010988 | Public |
| Birth weight / gestational age (fetal) | EGG consortium (SEM-partitioned files are EUR-only) | Public |
| Birth weight / gestational age (trans-ethnic) | JECS | Request form (manual placement) |
| Gestational diabetes / T2D glycemia | ProDiGY via the T2D Knowledge Portal | Registration (manual placement) |

## Cohort metadata extension contract (module 07)

GxE exposures and sensitivity covariates are read from
`placenta_QTL_cohort_metadata.txt`. Current columns: `rnaseq_id`,
`array_id`, `ancestry`, `cohort`, `sex`, `GA`, `ppBMI`. To enable the
planned GDM/OGTT scans and the SNP x covariate sensitivity analysis, add
columns named exactly `gdm` (0/1), `ogtt` (continuous), `maternal_age`,
`parity`, `smoking`, then flip the corresponding `enabled` flags in
`07_gxe_mapping/gxe_config.tsv` / `sensitivity_config.tsv`. The
transmitted/non-transmitted decomposition additionally needs maternal
pgens (`{COHORT}_maternal.pgen`) and a mother-child pairs file
(`child_array_id`, `mother_id`, `cohort`).

## Reference data

| Resource | Version / accession | Use |
|---|---|---|
| GRCh38 no-alt analysis set | GCA_000001405.15 | Alignment, imputation normalization |
| HPLRv2 placental long-read transcriptome | [bioRxiv 2025.06.26.661362](https://doi.org/10.1101/2025.06.26.661362); [github.com/sbresnahan/lr-placenta-transcriptomics](https://github.com/sbresnahan/lr-placenta-transcriptomics) | Phenotyping annotation (SQANTI3-classified) |
| GENCODE | v45 (deduplicated) | Gene-boundary override in GTF normalization |
| REDIportal | TABLE1_hg38_v3 | RNA-editing site catalog |
| TOPMed imputation server | — | Genotype imputation |
| 1000 Genomes Project | phase 3 (plink2 pgen) | PC-AiR ancestry reference |

## Software

Tool and package versions: [`software_environments.md`](software_environments.md).
