# 03_phenotyping

Per-cohort RNA phenotyping: a rewrite of
[PejLab/Pantry](https://github.com/PejLab/Pantry) for the seadragon HPC
(LSF instead of Snakemake; Salmon instead of kallisto), targeting the
long-read HPLRv2 transcript annotation. Generates 8 RNA-modality phenotype
BEDs per cohort — expression, isoforms, alt_TSS, alt_polyA, splicing,
intron retention, RNA editing, stability — unnormalized, bgzipped, and
tabix-indexed, with `phenotype_groups.txt` gene groupings. Normalization
(QN + INT) is applied once after ancestry-stratified pooling in stage 5
(scripts 19/20), followed by ComBat (devBrain xQTL schema; Wen et al.,
Science 2024, 384:eadh0829).

A standalone config-driven version of this module for general use is at
[sbresnahan/pantry-seadragon](https://github.com/sbresnahan/pantry-seadragon);
the copy here is the record of the placenta project as run.

## Prerequisites

- `config.yml` with the shared reference keys and one block per cohort
  (`strandedness`, `fastq_dir`, `fastq_map`, `samples_file`); all stage
  scripts read it via `config_get.py`
- MAJIQ academic license (not distributed): obtain from the
  [MAJIQ project](https://majiq.biociphers.org/) and set `majiq_license`
  in `config.yml`
- Software environments per
  [../docs/software_environments.md](../docs/software_environments.md):
  STAR 2.7.4a, Salmon 1.10.2, samtools/htslib 1.16.1, RegTools,
  leafCutter, MAJIQ, Subread featureCounts, gffread, R 4.3.1 (singularity;
  edgeR, tximport, rtracklayer, sva, MuSiC, SingleCellExperiment), Python 3
  (pandas, numpy, gtfparse, scipy, matplotlib, PyYAML)
- Create the log directory once before the first driver run; seadragon LSF
  does not auto-create log directories:

  ```bash
  mkdir -p "$OUTPUT_BASE/logs"
  ```

## Step 1: shared reference prep (run once)

Three LSF jobs chained by job-name dependencies. Submit all three
back-to-back: LSF resolves `done(jobname)` only if a job with that name
exists at submission time, so do not wait for one to finish before
submitting the next.

```bash
cd 03_phenotyping
bsub -env "CONFIG=$PWD/config.yml,SCRIPTS_DIR=$PWD" < 001_refprep_pre.sh
bsub -env "CONFIG=$PWD/config.yml,SCRIPTS_DIR=$PWD" < 002_txrevise_array.sh
bsub -env "CONFIG=$PWD/config.yml,SCRIPTS_DIR=$PWD" < 003_refprep_post.sh
```

- `001_refprep_pre.sh`: HPLRv2 GTF normalization (`prep_PANTRY_gtf.R`),
  REDIportal preprocessing (`rediportal_preprocess.py`), exonic/intronic
  GTFs (`exonic_intronic_from_gtf.py`), MAJIQ gff3
  (`gtf_to_majiq_gff3.py`), edit-site-to-gene map, txrevise prep.
- `002_txrevise_array.sh`: txrevise event construction as a 100-batch job
  array. The array range `#BSUB -J "txrevise[1-100]"` must match
  `txrevise_n_batches` in `config.yml`.
- `003_refprep_post.sh`: merge per-batch GFF3s; build the 6 txrevise
  Salmon indices.

## Step 2: per-cohort processing

`run_pipeline.py` is the driver. For one cohort it stages the cohort
inputs, submits per-sample alignment and quantification jobs with LSF
`-w done()` dependencies, submits per-modality aggregation jobs and the
final indexing job, and monitors by polling `bjobs` (reporting EXIT jobs
plus stranded downstream dependents).

```bash
bsub -env "CONFIG=$PWD/config.yml,SCRIPTS_DIR=$PWD,COHORT=cohort1" < 00_run_pipeline.sh
# or interactively:
python3 run_pipeline.py --config config.yml --cohort cohort1
```

Useful flags: `--dry-run` (print bsub commands without submitting),
`--skip-alignment` (reuse existing BAMs), `--skip-rna-editing`,
`--no-monitor`, `--monitor-all`, `--poll-interval N` (default 60 s).
Extra driver args can be passed through the wrapper via `DRIVER_ARGS`.

Per-sample stages: `01_star_align.sh` (STAR) -> `02_salmon_expression.sh`
(Salmon, expression index, 20 bootstraps), `03_salmon_alt_tss_polya.sh`
(Salmon x 6 txrevise indices), `04_regtools_junctions.sh` (splice
junctions), `05_featureCounts.sh` (exonic + intronic counts),
`06_rna_editing_pileup.sh` (mpileup at REDIportal sites),
`07_majiq_build.sh` (MAJIQ splice graph). Aggregation stages 10-15 write
the per-modality BEDs; `16_index_outputs.sh` verifies and tabix-indexes
them.

`10_aggregate_expression.sh` runs `qu_correct_salmon.R` (edgeR
`catchSalmon` quantification-uncertainty correction from the 20
bootstraps; Chen et al., NAR 2023, doi:10.1093/nar/gkad1167): isoform
BEDs are built from the corrected counts, gene-level expression from the
original Salmon counts.

Run the driver once per cohort. On failure the driver reports the failed
stage/sample and stranded dependents, then exits non-zero; fix the cause
and rerun the same command.

## Step 3: cross-cohort harmonization (splicing + intron retention)

Required before cross-cohort pooling of splicing and intron retention:
leafCutter cluster and MAJIQ intron IDs are cohort-specific, so features
are harmonized by stable genomic coordinates. Run once per ancestry after
all cohorts finish step 2:

```bash
bsub -env "CONFIG=$PWD/config.yml,SCRIPTS_DIR=$PWD,ANCESTRY_MAP=<sample_id/ancestry TSV>,ANCESTRY=EUR" \
  < 17_harmonize_within_ancestry.sh
```

`harmonize_within_ancestry.py` pools per-cohort leafCutter/MAJIQ
intermediates and writes pre-pooled unnormalized BEDs with namespaced
sample IDs (`{cohort}_{rnaseq_id}`). Stage-5 script 20 consumes these via
`--pre-pooled-dir $HARMONIZE_DIR`; 17's output root and 20's
`HARMONIZE_DIR` must point at the same directory (both honor env
overrides). `test_harmonize.py` is the test suite.

## Step 4: cell-type deconvolution

```bash
bsub -env "CONFIG=$PWD/config.yml,SCRIPTS_DIR=$PWD" < 18_deconvolution.sh
```

The driver submits three stages: `18a_build_reference.sh` (MuSiC reference
from the Campbell et al. 2023 placenta scRNA-seq reference, GSE182381;
`build_music_reference.R`, `reference_gene_symbols.txt`),
`18b_deconvolve_cohort.sh` per cohort (`build_bulk_tpm_matrix.py` +
`run_music_deconvolution.R`), and `18c_pool_outputs.sh`
(`pool_deconvolution_outputs.py`), with `check_cell_proportions.R` QC.
Cell-type fractions become mapping covariates (stage 5, script 25).

## Outputs

Per cohort, under `output/`: `expression.bed.gz`, `isoforms.bed.gz`,
`alt_TSS.bed.gz`, `alt_polyA.bed.gz`, `splicing.bed.gz`,
`intron_retention.bed.gz`, `RNA_editing.bed.gz`, `stability.bed.gz`
(+ `.tbi`; + `phenotype_groups.txt` for grouped modalities; +
`.site_to_phenotype.tsv` for RNA editing). All BEDs are unnormalized;
normalization happens in stage 5.

## Operational notes

- `assemble_bed.py` here is the one modified PANTRY script (reads Salmon
  `quant.sf` instead of kallisto `abundance.tsv`); it must replace the
  bundled original in the PANTRY scripts directory on deployment. Bundled
  PANTRY scripts used unmodified ship in `txrevise/`, `RNA_editing/`,
  `intron_retention/`.
- Stage-3 `output/<modality>.bed.gz` files are unnormalized. Consumers
  that expect normalized per-cohort BEDs (e.g. `combine_modalities.sh`)
  must be pointed at stage-5 outputs instead.
- Submit 001/002/003 back-to-back (LSF job-name dependency resolution;
  see step 1).
- The txrevise array range in `002_txrevise_array.sh` must equal
  `txrevise_n_batches` in `config.yml`.
- Create `$OUTPUT_BASE/logs/` before the first driver run.

## Utilities

`combine_modalities.sh` (concatenate modality BEDs into a cross-modality
BED), `regenerate_refflat.sh` (11-field refFlat for Picard),
`test_prepare_data.sh` / `test_subset_fastq.py` (pilot-run helpers),
`write_laddr_config.py` (optional LaDDR config writer; not part of the
production flow).

Next: [../04_sample_linkage/README.md](../04_sample_linkage/README.md).
