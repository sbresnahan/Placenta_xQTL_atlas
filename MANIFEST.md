# Project manifest — raw data to xQTL summary report

This manifest is the stage-by-stage reproducibility record for the multi-ancestry
placental xQTL project. Each stage lists its purpose, scripts, key inputs, key
outputs, software, and run-order dependencies. Stage directories are numbered in
execution order; on the seadragon cluster, stages 3 and 5 run directly from a
git clone of this repository and share one `config.yml`. Every entry-point
wrapper defaults `SCRIPTS_DIR` to its own directory (`${BASH_SOURCE[0]}`), so
no path wiring is needed for direct invocation; an explicit `SCRIPTS_DIR` env
var still overrides (and remains required for `bsub < script` submission —
see landmine 12).

> **Documentation convention.** Descriptions here were written against the
> 2026-09-19 script set (the current, GTEx-conventions versions) — script headers
> and usage blocks are the ground truth. Older archived READMEs reflect earlier
> pipeline versions and may not match the shipped scripts.

**Cluster context.** All stages run on the MD Anderson seadragon HPC: LSF scheduler
(`#BSUB` directives, `bsub`, job-name/job-ID dependencies), RHEL8, conda envs via
miniforge3, R 4.3.1 via a singularity image
(`singularity exec --bind /rsrch5 --bind /rsrch9 .../rstudio_4.3.1.sif Rscript`),
and data roots `/rsrch5` (software/references) and `/rsrch9` (project data).
Key cluster paths (as run):

```
OUTPUT_BASE=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY
QTL_DIR=$OUTPUT_BASE/qtl_inputs          # mapping inputs (BEDs, covariates, pgens)
RESULTS_DIR=$OUTPUT_BASE/qtl_results     # mapping outputs (parquets, top tables)
pooled genotypes: /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/pooled/genotypes/{ANC}_pooled.pgen
RNA→DNA map:      /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/pooled/rnaseq_to_array_id_map.csv
```

---

## Stage 0 — Cohorts and raw data

Four cohorts are carried through the current phase; pooled genotypes span five
ancestry strata (AFR, AMR, EAS, EUR, SAS). MALI_G3A genotypes are pooled but
excluded this phase (no matching RNA).

| Cohort | PANTRY dir | RNA-seq ID type | Genotype ID type | Access |
|---|---|---|---|---|
| NIEHS_RICHS | cohort1 | SRR Run | SUBJECT_ID_Array (S1, S2, …) | dbGaP |
| GUSTO | cohort2 | "J"+ID | SubjectID_BSubjectID | study investigators |
| SNUH | cohort3 | SRR Run | isolate (OGF###) | SRA PRJNA820329 |
| NIGMS | cohort4 | SRR Run (placenta only) | WXS Run (SRR) | SRA PRJNA671171 |

NIGMS subjects have two placenta sampling sites (PLAC-RNA1/-2); replicates are
averaged per individual during harmonization. Accession-level detail and
ancestry × cohort counts: [`docs/data_availability.md`](docs/data_availability.md).

---

## Stage 1 — Genotype imputation and pooling
📁 [`01_genotype_imputation/`](01_genotype_imputation/)

**Purpose.** Take each cohort from raw genotype data (array genotypes or WXS) to
imputed, QC'd, ancestry-assigned genotypes, then pool cohorts into per-ancestry
plink2 files for mapping.

**Flow (per cohort, then pooled):**

1. **Retrieval** — `pull_RICHS_genotypes.lsf` (dbGaP array), `pull_Yale_genotypes.lsf`,
   `retrieve_SNUH.lsf`, `pull_NIGMS_exomes.lsf` (SRA WXS).
2. **NIGMS WXS → genotypes** — `NIGMS_preprocessing.lsf` (+ `loop_NIGMS_preprocessing.sh`;
   SRA → FASTQ → alignment → gVCF), then `merge_NIGMS_gvcfs_setup.lsf` +
   `merge_NIGMS_gvcfs_chr.lsf` (index gVCFs; per-chromosome joint genotyping).
3. **Liftover** — `liftover_RICHS_genotypes.lsf` (array genotypes → GRCh38).
4. **Reference prep** — `BWA_ref.lsf`, `prep_dbSNP.lsf`, `extract_topmed.lsf`,
   `extract_G3A.lsf`.
5. **Pre-imputation filtering** — `filter_vcf_for_imputation*.lsf` (per cohort):
   SNPs only, biallelic only, genotype missingness / MAF / HWE filters.
6. **TOPMed imputation** — run on the TOPMed imputation server (off-cluster);
   results pulled with `pull_imputed_{RICHS,NIGMS,MALI,MALI_G3A}_genotypes.lsf`.
7. **Post-imputation QC** — `prep_QC_*.lsf` per cohort (Rsq ≥ 0.4 → plink2 pgen);
   QC summaries via `cohort_imputationQC.R` (per cohort), `mega_imputationQC.R`
   (pooled), `check_vcf_quality.R` (pre-imputation VCF checks).
8. **Ancestry assignment** — `prep_pcAiR.lsf` (build the 1000 Genomes reference
   pgen) → `PC_AiR.R` (GENESIS PC-AiR of cohort samples against 1KG →
   `assigned_ancestry`; `AFR_check.R` + `AFR_LD.bed` support AFR-stratum checks).
   *This is the current ancestry-assignment method and is distinct from the
   retired mapping-stage PC-AiR genotype PCA (see
   [`docs/round_history.md`](docs/round_history.md)).*
9. **Cross-cohort pooling** — `mega_merge_and_filter.lsf` (+ `mega_merge_and_filter_check.lsf`):
   builds `pooled_sample_ancestry.tsv`, subsets each cohort's Rsq-passing VCF to
   ancestry-stratum samples, splits multiallelics per cohort (DS correctly
   partitioned), merges across cohorts with `bcftools merge --merge none`,
   normalizes ref/alt against GRCh38, converts to pgen, applies pooled MAF/HWE/
   missingness filters, drops palindromic duplicates →
   `{ANC}_pooled.{pgen,pvar,psam}` + per-ancestry passing-variant lists and QC TSVs.
10. **Housekeeping** — `cleanup.sh`.

**Key outputs.** `{ANC}_pooled.{pgen,pvar,psam}` (AFR, AMR, EAS, EUR, SAS);
`pooled_sample_ancestry.tsv` (sample_id, assigned_ancestry, cohort — the genotype
ID space for all downstream linkage).

**Software.** plink/plink2, bcftools, BWA, GATK (gVCF merge), sratoolkit,
R (vcfR, data.table, GENESIS, GWASTools, SNPRelate), TOPMed imputation server.

---

## Stage 2 — RNA-seq quantification
📁 [`02_rnaseq_quantification/`](02_rnaseq_quantification/)

**Purpose.** Fetch and quantify RNA-seq for cohorts whose FASTQs were not already
local, ahead of PANTRY phenotyping.

**Scripts.**

| Script | Purpose |
|---|---|
| `pull_sra_SNUH.lsf` | `prefetch` + `fasterq-dump --split-files` + gzip per accession (LSF transfer queue) |
| `loop_salmon_{NIEHS_RICHS,SNUH}.sh` | Submit one quant job per FASTQ |
| `run_salmon_{NIEHS_RICHS,SNUH}.lsf` | Trim Galore → `salmon quant -l A --validateMappings --seqBias --numBootstraps 50` |

**Notes.** GUSTO FASTQs were already local; NIGMS RNA-seq runs are pulled from SRA
(PRJNA671171). Production phenotyping re-quantifies per sample inside Stage 3
(`02_salmon_expression.sh`, 20 bootstraps, HPLRv2 index from `config.yml`); the
wrappers here document the initial fetch/trim/quant for RICHS and SNUH.

**Software.** sratoolkit, Trim Galore 0.6.10, Salmon 1.10.2.

---

## Stage 3 — Transcriptomic phenotyping (PANTRY, seadragon rewrite)
📁 [`03_phenotyping/`](03_phenotyping/)

**Purpose.** From per-cohort FASTQs, generate **8 RNA-modality phenotype BEDs**
(unnormalized, bgzipped + tabix-indexed, with `phenotype_groups.txt` gene
groupings) ready for cross-cohort pooling. Per-cohort QN+INT was retired in the
2026-09 schema revision: normalization is applied once after ancestry-stratified
pooling in stage 5 (scripts 19/20), modeled on the devBrain xQTL atlas
(Wen et al., Science 2024, 384:eadh0829).

| Module | Modalities | Tool(s) |
|---|---|---|
| align | (shared BAM) | STAR, samtools |
| expression | expression, isoforms | Salmon (+ edgeR `catchSalmon` QU correction for isoforms) |
| alt_TSS_polyA | alt_TSS, alt_polyA | txrevise + Salmon |
| splicing | splicing | RegTools + leafCutter |
| intron_retention | intron_retention | MAJIQ |
| RNA_editing | RNA_editing | samtools mpileup (REDIportal sites) |
| stability | stability | Subread featureCounts (exon/intron ratio) |

**Tier 1 — shared reference prep (run once).** Chained LSF jobs:
`001_refprep_pre.sh` (HPLRv2 GTF normalization via `prep_PANTRY_gtf.R`, REDIportal
preprocessing, exonic/intronic GTFs, MAJIQ gff3, edit-site→gene map, txrevise prep)
→ `002_txrevise_array.sh` (100-batch txrevise event construction) →
`003_refprep_post.sh` (merge per-batch GFF3s; build 6 txrevise Salmon indices).

**Tier 2 — per cohort.** `run_pipeline.py` (driver; submits with LSF `-w done()`
dependencies and monitors; `00_run_pipeline.sh` wraps the driver as one LSF job):
`01_star_align.sh` → per-sample quantification (`02_salmon_expression.sh`,
`03_salmon_alt_tss_polya.sh`, `04_regtools_junctions.sh`, `05_featureCounts.sh`,
`06_rna_editing_pileup.sh`) → per-modality aggregation (`10_aggregate_expression.sh`
with `qu_correct_salmon.R` bootstrap QU correction for isoforms,
`11_aggregate_alt_tss_polya.sh`, `12_aggregate_splicing.sh`,
`13_aggregate_intron_retention.sh`, `14_aggregate_rna_editing.sh`,
`15_aggregate_stability.sh`) → `16_index_outputs.sh`.

**Stage 17 — cross-cohort coordinate harmonization (required for splicing +
intron retention).** `17_harmonize_within_ancestry.sh` runs
`harmonize_within_ancestry.py` per ancestry: pools per-cohort leafCutter/MAJIQ
intermediates, harmonizes features by stable genomic coordinates, and writes
pre-pooled unnorm BEDs with namespaced sample IDs (`{cohort}_{rnaseq_id}`).
Stage-5 script 20 consumes these via `--pre-pooled-dir $HARMONIZE_DIR`; without
this step, splicing and intron retention cannot enter cross-cohort ComBat pooling.
(`test_harmonize.py` is its test suite. Deployment detail: 17's output root and
20's `HARMONIZE_DIR` must point at the same directory — both honor env overrides.)

**Stage 18 — cell-type deconvolution.** `18_deconvolution.sh` orchestrates
`18a_build_reference.sh` (MuSiC reference via `build_music_reference.R`),
`18b_deconvolve_cohort.sh` (`run_music_deconvolution.R` on
`build_bulk_tpm_matrix.py` TPM matrices), `18c_pool_outputs.sh`
(`pool_deconvolution_outputs.py`). Cell-type fractions become mapping covariates
(Stage 5, script 25).

**Configuration.** `config.yml` (shared references + per-cohort blocks) read via
`config_get.py`; the driver passes `--config`/`--scripts-dir` to each stage script.
Bundled PANTRY scripts ship in `txrevise/`, `RNA_editing/`, `intron_retention/`;
`assemble_bed.py` is the one **modified** PANTRY script (reads Salmon `quant.sf`
instead of kallisto) and must replace the bundled original on deployment.

**Key outputs (per cohort, `output/`).** `expression.bed.gz`, `isoforms.bed.gz`,
`alt_TSS.bed.gz`, `alt_polyA.bed.gz`, `splicing.bed.gz`, `intron_retention.bed.gz`,
`RNA_editing.bed.gz`, `stability.bed.gz` (+ `.tbi`, + `phenotype_groups.txt` where
applicable).

**Software.** STAR 2.7.4a, Salmon 1.10.2, samtools/htslib 1.16.1, RegTools,
leafCutter, MAJIQ (academic license — not distributed), Subread featureCounts,
gffread, R 4.3.1 (singularity; edgeR, tximport, rtracklayer, sva, MuSiC,
SingleCellExperiment, data.table, optparse, ggplot2), Python 3 (pandas, numpy,
gtfparse, scipy, matplotlib, PyYAML).

---

## Stage 4 — RNA→DNA sample linkage
📁 [`04_sample_linkage/`](04_sample_linkage/)

**Purpose.** Build `rnaseq_to_array_id_map.csv` (columns: `rnaseq_id`, `array_id`,
`ancestry`, `cohort`) linking PANTRY BED sample IDs to pooled-genotype sample IDs.
Runs before QTL mapping — Stage-5 script 23 intersects on this map.

**Script.** `build_RNA_to_DNA_map.R` (base R, no dependencies; run interactively).
Verified per-cohort linkage logic:

- **NIEHS_RICHS**: key file Run (SRR) → SUBJECT_ID_Array
- **GUSTO**: BED IDs are "J"+covars$ID; pooled ID = `SubjectID_B<SubjectID>`
- **SNUH**: SRA metadata Run → isolate (OGF###)
- **NIGMS**: placenta RNA-Seq runs (excludes decidua; both fetal sampling sites)
  → submitted_subject_id → WXS Run

Includes join diagnostics (match rates, duplicate/NA guards) and optional BED-header
coverage validation. NIGMS array IDs legitimately appear twice (two sampling sites);
replicates are averaged per individual during Stage-5 harmonization.

---

## Stage 5 — cis-xQTL mapping (tensorQTL, GTEx conventions)
📁 [`05_qtl_mapping/`](05_qtl_mapping/)

**Purpose.** From pooled genotypes + modality BEDs, produce ancestry-stratified
cis-xQTL summary statistics. Current baseline = the GTEx-conventions design
(see [`docs/round_history.md`](docs/round_history.md) for design evolution).

**Run order** (per ancestry; scripts 19–20 run once across ancestries):

| # | Script | What it does |
|---|---|---|
| 19 | `19_hcp_factors.sh` (+ `19a_picard_sharded.sh`, `19b_fix_missing_metrics.sh`) | HCP pipeline: Picard QC per sample (`picard_qc.py`) → pool metrics (`fix_missing_metrics.py`) → pool expression within ancestry (`pool_expression_within_ancestry.py`) → TPM > 0.1 in > 25% filter → QN + INT → connectivity-outlier removal (signed bicor network, z < −3; writes `{ANC}_expression_outliers.tsv`) → ComBat (batch=cohort, **last**) → HCP estimation, provisional k=15 (`combat_normalize_hcp.R`, `peer_factors.py`) |
| 20 | `20_combat_modalities.sh` | Pool + QN + INT + ComBat (ComBat **last**) for the 7 non-expression modalities (`pool_modalities_within_ancestry.py`, `combat_normalize_modalities.R`); devBrain filters (isoforms TPM > 0.1 in > 25%; proportion modalities detected in ≥ 40%); no log2/logit pre-transform (no-op under rank-based QN+INT); isoforms exclude the expression-outlier samples from 19 (`--exclude-samples`); splicing/IR pass through stage-17 pre-pooled BEDs. HCP factors from 19 are reused, not re-estimated |
| 21 | `21_install_tensorqtl.sh` | One-time env setup: python-only tensorqtl conda env + R `qvalue` into the singularity R library + Storey-bridge smoke test |
| 22 | `22_genotype_pca.sh` | Cohort-only genotype PCA (GTEx convention): LD-prune (`--indep-pairwise 200 50 0.2`) → `plink2 --pca 20 exact` → `genotype_pca_format.py` → `{ANC}_genotype_pcs.tsv` + scree (`PCA_scree.R`) |
| 23 | `23_prepare_intersection.py` | Genotype×phenotype intersection on the RNA→DNA map; pgen filtered to intersection samples with **MAC ≥ 5** (`--mac 5`; `--mac 0` disables); sample columns renamed rnaseq_id → array_id |
| 24 | `24_outlier_exclusion.py` | Select first 5 genotype PCs; 6-SD outlier exclusion on PC1–5; removes outliers from pgen/BED/HCP/deconvolution/metadata. **Edits intersection files in place — always rerun 23 before 24** |
| 25 | `25_build_covariates.py` | Covariates: PC1–5 + HCP_1–k + sex + GA + cell types (dominant type as compositional reference); near-zero-variance pre-filter; iterative \|r\| > 0.9 pruning (priority sex/GA > PCs > cell types > HCPs); cap ≤ 25. Supports `--hcp-file`/`--hcp-k`/`--out-suffix` for the 25a optimization module. **Note: the maternal cell-fraction covariate was manually removed after script 25 (25 → 24 covariates); rerunning 25 restores it — re-apply the edit or script it** |
| 25a | `25a_optimize_hcp.sh` + `optimize_hcp_chr1.py` | Pre-mapping HCP-count optimization (devBrain §4.2), run after 23/24, before canonical 25: per ancestry, re-estimate HCP at each k ∈ {0,5,10,15,20,25,30}, build per-k covariates, map **chr1 expression only** with tensorQTL (1 Mb window, MAF ≥ 0.01), count eGenes at Storey q ≤ 0.05; k\* = argmax (ties → smaller k). Installs the k\* solution as `{ANC}_hcp_factors_harmonized.tsv` (provisional file backed up to `*.pre25a_backup.tsv`); writes `{ANC}_optimal_hcp.tsv` + `.png` |
| 26 | `26_harmonize_modalities.py` | Harmonize the 7 non-expression modality BEDs to the final array_id sample set (handles stage-17 namespaced IDs for splicing/IR) |
| 30 | `30_combine_modalities.py` | Build the combined cross-modality BED (phenotype IDs namespaced `{modality}__{id}`; cross-modality gene groups; modality sidecar TSV) |
| 27 | `27_run_tensorqtl.sh` + `27_run_tensorqtl.py` | cis mapping per ancestry × modality: grouped (`group_s`, one lead per gene) when a `phenotype_groups.txt` exists, ungrouped otherwise; `--independent` for PANTRY-style stepwise conditional signals; q-values on `pval_beta` via `--qvalue-method storey` (default; R `qvalue` through the `compute_qvalues.R` file bridge, `QVALUE_RSCRIPT` env var) or `bh` (escape hatch) |
| 28 | `28_submit_modalities.sh` | Submission driver: one LSF job per ancestry × modality (`TEST=1` pilot mode; threads `QVALUE_METHOD`/`MAF_THRESHOLD`) |
| 29 | `29_make_top_tables.py` | Rebuild sorted `*_cisqtl_top.tsv` from parquets (no tensorQTL rerun) |

**Diagnostics/utilities.** `hcp_diagnostic.R`, `hcp_diagnostic2.R` (HCP–genotype
correlation diagnostics), `check_hcp_chunks.sh`, `test_hcp.py`,
`test_combat_modalities.py`.

**Key inputs.** Stage-1 `{ANC}_pooled.pgen`; Stage-3 modality BEDs + deconvolution
proportions; Stage-4 `rnaseq_to_array_id_map.csv`; cohort metadata
(`rnaseq_id, array_id, ancestry, cohort, sex, GA, ppBMI`).

**Key outputs.** `$QTL_DIR`: intersected pgens, `{ANC}_genotype_pcs.tsv`,
`{ANC}_covariates.tsv`, harmonized BEDs. `$RESULTS_DIR`:
`{ANC}_{modality}[_ungrouped]_cisqtl.parquet` + `_top.tsv` per layer
(grouped/ungrouped/combined) + `_independent` stepwise outputs.

**Software.** tensorqtl 1.0.10 (python-only conda env: pandas, pyarrow, genotypeio,
matplotlib), plink2, R 4.3.1 singularity (qvalue, sva, Rhcpp, data.table, ggplot2,
patchwork), bgzip/tabix.

---

## Stage 6 — Summary report
📁 [`reports/`](reports/)

**Purpose.** Render the internal QC + results report from mapping outputs.

**Files.** `report_placenta_xqtl.Rmd` (report source),
`xqtl_report_functions.R` (helper module: loaders, QQ/power/concordance helpers),
`gene_map.tsv` (gene_id → symbol/biotype/description), and the rendered
self-contained `report_placenta_xqtl.html`.

**Render.**

```r
rmarkdown::render(
  "reports/report_placenta_xqtl.Rmd",
  params = list(
    data_dir         = "<dir with *_cisqtl_top.tsv>",   # Stage-5 top tables
    qc_dir           = "<dir with QC files>",
    module_path      = "reports/xqtl_report_functions.R",
    gene_map_path    = "reports/gene_map.tsv",
    gene_bodies_path = "<gene_bodies.tsv>",
    fig_dir          = "<fig out>", table_dir = "<table out>"
  )
)
```

---

## Do-not-repeat landmines

1. **rpy2/conda-R is a dead end on seadragon** — two attempts hit
   `GLIBCXX_3.4.30 not found` (libstdc++ loader conflict between the conda env's
   R/libicu and system libstdc++). Storey q-values go through the file-based
   singularity-R bridge (`compute_qvalues.R`), which is statistically identical to
   tensorQTL's `calculate_qvalues` (itself an rpy2 wrapper on `qvalue::qvalue`).
   The `'rfunc' cannot be imported` warning at `import tensorqtl` is cosmetic.
2. **`27_run_tensorqtl.sh` skips bgzip/tabix if `{ANC}_{modality}.bed.gz` exists** —
   stale harmonized BEDs are silently reused unless the `.bed.gz`/`.tbi` files are
   deleted before a rerun.
3. **`24_outlier_exclusion.py` edits the intersection files in place** — always
   rerun 23 before 24.
4. **The covariate files were hand-edited after script 25** (maternal cell-fraction
   covariate removed; 25 → 24 covariates per ancestry). Rerunning 25 restores the
   25-covariate version — re-apply the edit or script it.
5. **Tier-1 LSF name dependencies** — submit `001`/`002`/`003` back-to-back; LSF
   resolves `done(jobname)` only if a job with that name exists at submission time.
6. **seadragon LSF does not auto-create log directories** — create
   `$OUTPUT_BASE/logs/` once before first driver run.
7. **txrevise batch count must match** across `config.yml` (`txrevise_n_batches`)
   and the hardcoded `#BSUB -J "txrevise[1-100]"` array range in `002`.
8. **The modified `assemble_bed.py` must replace the bundled PANTRY original** in
   the PANTRY scripts directory on deployment (Salmon `quant.sf` reader).
9. **`optimize_hcp_chr1.py` (25a) overwrites `{ANC}_hcp_factors_harmonized.tsv`**
   with the k\* solution (the provisional k=15 file is backed up to
   `*.pre25a_backup.tsv`). Always rerun 25 after 25a — and re-apply the
   maternal-fraction removal (landmine 4).
10. **Stage-3 `output/<modality>.bed.gz` files are now UNNORMALIZED** (2026-09
    schema: normalization moved to stage 5). Legacy consumers expecting
    per-cohort normalized BEDs (e.g. `combine_modalities.sh`) must be pointed
    at stage-5 outputs instead.
11. **`picard` must resolve after 19/19a's env stack** — `picard_qc.py` calls
    the literal `picard` executable; if the `picard-2.27.4` activation doesn't
    put it on PATH (wrong env name/path, or a PATH-clobbering venv
    activation), every QC call WARNs `[Errno 2] No such file or directory:
    'picard'` and writes garbage metrics. Both wrappers fail fast with
    diagnostics instead — fix the activation, don't shim around it. Stage 1
    skips existing `{COHORT}_qc_metrics.tsv` — delete tainted QC outputs
    before resubmitting.
12. **`bsub < script` still requires an explicit `SCRIPTS_DIR` export** — LSF
    executes a spool *copy* of the submitted script, so `${BASH_SOURCE[0]}`
    self-location resolves to the spool directory, not the repo. The wrappers
    fail fast with a clear error in that case. Direct invocation
    (`bash script.sh`, interactive) self-locates and needs no SCRIPTS_DIR.
    Cross-directory references are resolved relative to the script:
    `config_get.py` is sought in `SCRIPTS_DIR` then `../03_phenotyping`;
    `picard_qc.py` (regenerate_refflat.sh) in `SCRIPTS_DIR` then
    `../05_qtl_mapping`.

## Porting to another system

- Replace seadragon paths (`/rsrch5`, `/rsrch9`) and infrastructure (conda env
  names/paths, singularity image, `R_LIBS_USER`, module loads) — these are
  intentionally hardcoded in scripts; data/reference paths are config-driven via
  `config.yml`.
- LSF → SLURM: convert `#BSUB` → `#SBATCH`, `bsub` → `sbatch`,
  `-w "done(...)"` → `--dependency=afterok:...`, and the driver's submission calls
  in `run_pipeline.py`.
- Any R ≥ 4.3 with the packages in [`docs/software_environments.md`](docs/software_environments.md)
  works; the singularity image is only needed to reproduce the exact seadragon stack.
