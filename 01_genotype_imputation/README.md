# 01_genotype_imputation

Per-cohort genotype processing: retrieval (array + WXS), liftover to GRCh38,
pre-imputation filtering, TOPMed imputation, Rsq QC, PC-AiR ancestry
assignment against 1000 Genomes, and cross-cohort pooling into per-ancestry
plink2 files (`{ANC}_pooled.pgen`) for QTL mapping. All scripts are LSF jobs
for the seadragon cluster. Cohort accessions:
[../docs/data_availability.md](../docs/data_availability.md).

## Prerequisites

- Cohort genotype data obtained through the routes in
  docs/data_availability.md (dbGaP, SRA, or study investigators)
- TOPMed imputation server account (off-cluster web service)
- Software: plink/plink2, bcftools, BWA, GATK, sratoolkit, R with vcfR,
  data.table, GENESIS, GWASTools, SNPRelate (see
  [../docs/software_environments.md](../docs/software_environments.md))

## Steps

Run per cohort, then pooled. Scripts are submitted with `bsub < script.lsf`;
paths at the top of each script are set for the seadragon layout.

1. **Retrieval** — `pull_RICHS_genotypes.lsf` (dbGaP array),
   `pull_Yale_genotypes.lsf` (GUSTO), `retrieve_SNUH.lsf`,
   `pull_NIGMS_exomes.lsf` (SRA WXS).
2. **NIGMS WXS to genotypes** — `loop_NIGMS_preprocessing.sh` submits one
   `NIGMS_preprocessing.lsf` job per sample (SRA -> FASTQ -> alignment ->
   gVCF); then `merge_NIGMS_gvcfs_setup.lsf` (index gVCFs) and
   `merge_NIGMS_gvcfs_chr.lsf` (per-chromosome joint genotyping).
3. **Liftover** — `liftover_RICHS_genotypes.lsf` (array genotypes to
   GRCh38).
4. **Reference prep** — `BWA_ref.lsf`, `prep_dbSNP.lsf`,
   `extract_topmed.lsf`, `extract_G3A.lsf`. The two extract scripts take
   TOPMed download credentials via the environment:

   ```bash
   bsub -env "DIR=<download dir>,PASSWORD=<topmed password>" < extract_topmed.lsf
   ```

5. **Pre-imputation filtering** — `filter_vcf_for_imputation.lsf` and the
   per-cohort variants (`_NIEHS`, `_NIGMS`, `_SNUH`, `_MALI`): SNPs only,
   biallelic only, genotype missingness / MAF / HWE filters.
   `check_vcf_quality.R` runs pre-imputation VCF checks.
6. **TOPMed imputation** — run on the TOPMed imputation server
   (off-cluster); pull results with
   `pull_imputed_{RICHS,NIGMS,MALI,MALI_G3A}_genotypes.lsf`.
7. **Post-imputation QC** — `prep_QC_{GUSTO,MALI,MALI_G3A,NIEHS_RICHS,NIGMS,SNUH}.lsf`
   (Rsq >= 0.4, convert to plink2 pgen). QC summaries:
   `cohort_imputationQC.R` (per cohort), `mega_imputationQC.R` (pooled).
8. **Ancestry assignment** — `prep_pcAiR.lsf` (build the 1000 Genomes
   reference pgen), then `PC_AiR.R` (GENESIS PC-AiR of cohort samples
   against 1KG, writes `assigned_ancestry`). `AFR_check.R` support AFR-stratum checks.
9. **Cross-cohort pooling** — `mega_merge_and_filter.lsf` (+
   `mega_merge_and_filter_check.lsf`): builds `pooled_sample_ancestry.tsv`,
   subsets each cohort's Rsq-passing VCF to ancestry-stratum samples, splits
   multiallelics per cohort (DS correctly partitioned), merges across
   cohorts with `bcftools merge --merge none`, normalizes ref/alt against
   GRCh38, converts to pgen, applies pooled MAF/HWE/missingness filters,
   drops palindromic duplicates.
10. **Housekeeping** — `cleanup.sh`.

## Outputs

- `{ANC}_pooled.{pgen,pvar,psam}` for AFR, AMR, EAS, EUR, SAS
- `pooled_sample_ancestry.tsv` (sample_id, assigned_ancestry, cohort) — the
  genotype ID space for all downstream linkage
- Per-ancestry passing-variant lists and QC TSVs

## Operational notes

- The PC-AiR scripts here perform ancestry assignment for pooling. The
  genotype PCs used as mapping covariates are computed separately in stage
  5 (`22_genotype_pca.sh`).
- `extract_topmed.lsf` / `extract_G3A.lsf` require `DIR` and `PASSWORD`
  via `bsub -env` (TOPMed server download credentials).

Next: [../02_rnaseq_quantification/README.md](../02_rnaseq_quantification/README.md).
