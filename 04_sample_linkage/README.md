# 04_sample_linkage — RNA-seq → genotype sample ID map

📜 `build_RNA_to_DNA_map.R` builds `rnaseq_to_array_id_map.csv`
(columns: `rnaseq_id`, `array_id`, `ancestry`, `cohort`), linking PANTRY BED
sample IDs to pooled-genotype sample IDs across the four cohorts. Base R only;
run interactively. The QTL-mapping intersection (stage 5, script 23) consumes
this map, so it runs **after** phenotyping and genotype pooling and **before**
QTL mapping.

## Verified per-cohort linkage logic

| Cohort | RNA-seq ID | Genotype (pooled) ID | Key |
|---|---|---|---|
| NIEHS_RICHS | SRR Run | SUBJECT_ID_Array (S1, S2, …) | dbGaP key file |
| GUSTO | "J" + covars$ID | `SubjectID_B<SubjectID>` | cohort covariates file |
| SNUH | SRR Run | isolate (OGF###) | SRA metadata |
| NIGMS | SRR Run (placenta only) | WXS Run (SRR) | SRA metadata via submitted_subject_id |

NIGMS subjects have two placenta sampling sites (PLAC-RNA1/-2), so each NIGMS
`array_id` legitimately appears twice; replicates are averaged per individual
during QTL-input harmonization (stage 5).

The script prints per-cohort match diagnostics, guards against duplicate/NA IDs,
and optionally validates coverage against PANTRY `expression.bed.gz` headers.
Edit the path block at the top of the script before running.
