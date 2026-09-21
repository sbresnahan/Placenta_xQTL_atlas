# 03_phenotyping — PANTRY phenotyping (seadragon rewrite)

In-house rewrite of [PejLab/Pantry](https://github.com/PejLab/Pantry) for the MD
Anderson seadragon HPC (LSF/BSUB instead of Snakemake; Salmon instead of
kallisto), targeting the long-read **HPLRv2** transcript annotation. Generates
**8 RNA-modality phenotype BEDs** per cohort — quantile-normalized + rank-based
inverse-normal transformed, bgzipped + tabix-indexed, with
`phenotype_groups.txt` gene groupings — ready for tensorQTL.

Full stage documentation: [`../MANIFEST.md`](../MANIFEST.md) stage 3.

## Architecture

**Tier 1 — shared reference prep (run once).** Three LSF jobs chained by
job-name dependencies (submit back-to-back; see MANIFEST landmine 5):

```
001_refprep_pre.sh      GTF normalization (prep_PANTRY_gtf.R), REDIportal
                        (rediportal_preprocess.py), exonic/intronic GTFs
                        (exonic_intronic_from_gtf.py), MAJIQ gff3
                        (gtf_to_majiq_gff3.py), edit-site→gene map, txrevise prep
002_txrevise_array.sh   txrevise event construction, 100-batch job array
003_refprep_post.sh     merge per-batch GFF3s; build 6 txrevise Salmon indices
```

**Tier 2 — per cohort.** `run_pipeline.py` submits per-sample and aggregation
jobs with LSF `-w done()` dependencies and monitors them (polls `bjobs`, reports
EXIT jobs plus stranded downstream dependents). `00_run_pipeline.sh` wraps the
driver as a single LSF job (`COHORT=cohort1 bsub < 00_run_pipeline.sh`).

```
01_star_align.sh             per-sample STAR alignment
02_salmon_expression.sh      Salmon quant (expression index) + 20 bootstraps
03_salmon_alt_tss_polya.sh   Salmon quant × 6 txrevise indices
04_regtools_junctions.sh     splice-junction extraction
05_featureCounts.sh          exonic + intronic counts (stability)
06_rna_editing_pileup.sh     RNA-editing levels (mpileup at REDIportal sites)
10–15                        per-modality aggregation → normalized BEDs
16_index_outputs.sh          verify + tabix-index all outputs
```

`10_aggregate_expression.sh` runs `qu_correct_salmon.R` (edgeR `catchSalmon`
quantification-uncertainty correction from the 20 bootstraps; Chen et al., NAR
2023, doi:10.1093/nar/gkad1167) — isoform BEDs are built from the corrected
counts, gene-level expression from the original Salmon counts.

**Stage 17 — cross-cohort harmonization (required for splicing + intron
retention).** `17_harmonize_within_ancestry.sh` runs
`harmonize_within_ancestry.py` per ancestry: pools per-cohort leafCutter/MAJIQ
intermediates and harmonizes features by stable genomic coordinates (cohort
cluster/event IDs are not comparable across cohorts). Output: pre-pooled unnorm
BEDs with namespaced sample IDs, consumed by stage-5 script 20 via
`--pre-pooled-dir $HARMONIZE_DIR`. `test_harmonize.py` is the test suite.

**Stage 18 — cell-type deconvolution.** `18_deconvolution.sh` →
`18a_build_reference.sh` (`build_music_reference.R`,
`reference_gene_symbols.txt`), `18b_deconvolve_cohort.sh`
(`build_bulk_tpm_matrix.py` + `run_music_deconvolution.R`),
`18c_pool_outputs.sh` (`pool_deconvolution_outputs.py`);
`check_cell_proportions.R` QC. Cell-type fractions become mapping covariates.

## Configuration

`config.yml` holds shared reference paths + one block per cohort
(`strandedness`, `fastq_dir`, `fastq_map`, `samples_file`); all stage scripts
read it via `config_get.py`. Infrastructure paths (conda envs, singularity
image, R library) are hardcoded in scripts — seadragon-specific.

**Deployment note.** `assemble_bed.py` here is the one **modified** PANTRY
script (reads Salmon `quant.sf` instead of kallisto `abundance.tsv`); it must
replace the bundled original in the PANTRY scripts directory. Bundled PANTRY
scripts used unmodified ship in `txrevise/`, `RNA_editing/`,
`intron_retention/`. The MAJIQ academic license file is **not** distributed —
obtain one from the MAJIQ project and set `majiq_license` in `config.yml`.

## Utilities

`combine_modalities.sh` (concatenate modality BEDs → cross-modality BED),
`cleanup_intermediates.sh` / `cleanup_original_run.sh` (disk management),
`test_subset_fastq.py` / `test_prepare_data.sh` (pilot-run helpers),
`write_laddr_config.py` (optional LaDDR config writer; not part of the current
production flow), `regenerate_refflat.sh`.

## Outputs (per cohort, `output/`)

`expression.bed.gz`, `isoforms.bed.gz`, `alt_TSS.bed.gz`, `alt_polyA.bed.gz`,
`splicing.bed.gz`, `intron_retention.bed.gz`, `RNA_editing.bed.gz`,
`stability.bed.gz` (+ `.tbi`; + `phenotype_groups.txt` for grouped modalities;
+ `.site_to_phenotype.tsv` for RNA editing).
