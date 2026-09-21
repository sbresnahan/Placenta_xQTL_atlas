# 01_genotype_imputation — per-cohort imputation, ancestry assignment, pooling

Takes each cohort from raw genotype data (array genotypes or WXS) to imputed,
QC'd, ancestry-assigned genotypes, then pools cohorts into per-ancestry plink2
files (`{ANC}_pooled.pgen`) for QTL mapping. All scripts are LSF jobs for the
seadragon cluster. See [`../MANIFEST.md`](../MANIFEST.md) stage 1 for the full
flow and [`../docs/data_availability.md`](../docs/data_availability.md) for
cohort accessions.

## Flow

1. **Retrieval** — `pull_RICHS_genotypes.lsf`, `pull_Yale_genotypes.lsf`,
   `retrieve_SNUH.lsf`, `pull_NIGMS_exomes.lsf`.
2. **NIGMS WXS → genotypes** — `NIGMS_preprocessing.lsf` (+
   `loop_NIGMS_preprocessing.sh`): SRA → FASTQ → alignment → gVCF;
   `merge_NIGMS_gvcfs_setup.lsf` + `merge_NIGMS_gvcfs_chr.lsf`: per-chromosome
   joint genotyping.
3. **Liftover** — `liftover_RICHS_genotypes.lsf` (→ GRCh38).
4. **Reference prep** — `BWA_ref.lsf`, `prep_dbSNP.lsf`, `extract_topmed.lsf`,
   `extract_G3A.lsf`.
5. **Pre-imputation filtering** — `filter_vcf_for_imputation.lsf` and per-cohort
   variants (`_NIEHS`, `_NIGMS`, `_SNUH`, `_MALI`): SNPs only, biallelic only,
   missingness/MAF/HWE filters.
6. **TOPMed imputation** — off-cluster; results pulled with
   `pull_imputed_{RICHS,NIGMS,MALI,MALI_G3A}_genotypes.lsf`.
7. **Post-imputation QC** — `prep_QC_{GUSTO,MALI,MALI_G3A,NIEHS_RICHS,NIGMS,SNUH}.lsf`
   (+ `prep_GUSTO_redo.lsf`): Rsq ≥ 0.4 → plink2 pgen. QC summaries:
   `cohort_imputationQC.R`, `mega_imputationQC.R`, `check_vcf_quality.R`.
8. **Ancestry assignment** — `prep_pcAiR.lsf` (build 1000 Genomes reference pgen)
   → `PC_AiR.R` (GENESIS PC-AiR vs 1KG → `assigned_ancestry`);
   `AFR_check.R` + `AFR_LD.bed` support AFR-stratum checks.
9. **Pooling** — `mega_merge_and_filter.lsf` (+ `mega_merge_and_filter_check.lsf`):
   per-ancestry subsetting, per-cohort multiallelic splitting, `bcftools merge
   --merge none`, GRCh38 ref/alt normalization, pgen conversion, pooled
   MAF/HWE/missingness filters → `{ANC}_pooled.{pgen,pvar,psam}` +
   `pooled_sample_ancestry.tsv`.
10. **Housekeeping** — `cleanup.sh`.

## Notes

- The PC-AiR scripts here perform **ancestry assignment for pooling** — distinct
  from the retired mapping-stage PC-AiR genotype PCA (see
  [`../docs/round_history.md`](../docs/round_history.md)).
- `extract_topmed.lsf` / `extract_G3A.lsf` take `DIR` and `PASSWORD` via
  `bsub -env` (TOPMed server download credentials).
- JVM crash logs and other run artifacts were intentionally excluded from this
  directory.
