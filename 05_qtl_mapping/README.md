# 05_qtl_mapping

Ancestry-stratified cis-xQTL mapping under GTEx conventions: cross-cohort
pooling with modality-appropriate normalization. Gene and isoform
expression use TMM -> VST, while ratio modalities use pooled QN followed
by per-phenotype rank INT within each cohort (devBrain xQTL schema, Wen
et al., Science 2024, 384:eadh0829). Mapping includes cohort indicator
covariates, HCP latent-factor estimation with expression-specific (25a)
and non-expression/combined (25b) k optimization, cohort-only genotype
PCA, genotype x phenotype intersection with a minor-allele-count floor,
PC-outlier exclusion, covariate assembly, tensorQTL mapping with Storey
q-values, SuSHiE cross-ancestry fine-mapping, and functional enrichment
of high-PIP variants.

These scripts share `config.yml` / `config_get.py` with
`../03_phenotyping/` (one scripts directory on seadragon).

## Prerequisites

- Stage-1 `{ANC}_pooled.pgen`, stage-3 modality BEDs + deconvolution
  proportions, stage-4 `rnaseq_to_array_id_map.csv`, cohort metadata
  (`rnaseq_id, array_id, ancestry, cohort, sex`)
- If rebuilding stage-3 alt-TSS/polyA or stability BEDs manually, first
  follow the per-sample intermediate preflight/recovery instructions in
  `../03_phenotyping/README.md`; scripts 11/15 require complete stage-03
  Salmon and stage-05 featureCounts intermediates.
- Software: tensorqtl 1.0.10 (python-only conda env, installed by
  `21_install_tensorqtl.sh`), sushie >= 0.20 (dedicated conda env, see
  below), plink2, R 4.3.1 singularity image (qvalue, sva, Rhcpp,
  data.table, ggplot2, patchwork), bgzip/tabix, fastVEP 0.4.0 (see
  "Annotation sources" below)
- Environment setup used throughout this walkthrough:

  ```bash
  CONFIG=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/config.yml
  OUTPUT_BASE=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY
  QTL_DIR="$OUTPUT_BASE/qtl_inputs"        # mapping inputs (BEDs, covariates, pgens)
  RESULTS_DIR="$OUTPUT_BASE/qtl_results"   # mapping outputs (parquets, top tables)
  REPO_DIR=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/software/Pantry/phenotyping/scripts
  SCRIPTS_DIR="$REPO_DIR/05_qtl_mapping"
  LOG_DIR=/rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs
  mkdir -p "$LOG_DIR"
  ```

  Pooled genotypes live at
  `/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/pooled/genotypes/{ANC}_pooled.pgen`;
  the RNA-to-DNA map at
  `/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/pooled/rnaseq_to_array_id_map.csv`.

## Run order

Scripts 19-20 run once across ancestries; the rest run per ancestry.

| # | Script(s) | Purpose |
|---|---|---|
| pre | `collapse_replicates.py` | Optional pre-step: collapse technical-replicate runs to one column per individual (see "Technical replicates" below) |
| 19 | `19_hcp_factors.sh`, `19a_picard_sharded.sh` + `picard_persample.py`, `19b_fix_missing_metrics.sh` | HCP pipeline: Picard QC (`picard_qc.py`) -> pool metrics (`fix_missing_metrics.py`) -> pool expression within ancestry (`pool_expression_within_ancestry.py`) -> TPM > 0.1 in >= 25% filter -> TMM -> VST -> connectivity-outlier removal (bicor, z < -3; writes `{ANC}_expression_outliers.tsv`) -> provisional HCP (`normalize_expression_hcp.R`, `peer_factors.py`). 19a is sample-parallel (2026-10-04): a submit driver scans completed legacy chunk (`*.chunk*.qcmetrics.tsv`) and per-sample (`per_sample/<cohort>/<sample>.qcmetrics.tsv`) outputs, then self-submits one array per cohort with one index per REMAINING sample (1 core / 12 GB / 2 h; `DRYRUN=1` preview; `COHORTS`/`MAXCONC`/`QUEUE`/`MEM_GB`/`WALLTIME` overrides); `MERGE=1` keeps chunk rows verbatim, appends new rows, and imputes their missing cells with cohort-level medians (workers run `--no-impute`). 19b is legacy chunk repair only |
| 20 | `20_normalize_modalities.sh` | Pool non-expression modalities (`pool_modalities_within_ancestry.py`, union mode). `isoform_expression` uses pooled counts with TMM -> VST; ratio modalities use pooled QN -> within-cohort INT (`normalize_modalities.R`), with devBrain detection/no-variance/per-cohort guards. Splicing/IR use stage-17 pre-pooled BEDs. Final HCP counts are optimized later by 25b |
| 21 | `21_install_tensorqtl.sh` | One-time setup: tensorqtl conda env + R `qvalue` + Storey-bridge smoke test |
| 22 | `22_genotype_pca.sh` | Cohort-only genotype PCA: LD-prune (`--indep-pairwise 200 50 0.2`) -> `plink2 --pca 20 exact` -> `genotype_pca_format.py` -> `{ANC}_genotype_pcs.tsv` + scree (`PCA_scree.R`) |
| 23 | `23_prepare_intersection.py` | Genotype x phenotype intersection on the RNA-to-DNA map; pgen filtered to intersection samples with MAC >= 5 (`--mac 5`; `--mac 0` disables); sample columns renamed rnaseq_id -> array_id |
| 24 | `24_outlier_exclusion.py` | First-5-PC selection; 6-SD PC outliers removed from pgen/BED/HCP/deconvolution/metadata. Edits intersection files in place — always rerun 23 before 24 |
| 25 | `25_build_covariates.py` | Covariates: PC1-5 + HCP_1-k + sex + cell types (dominant type as compositional reference) + cohort indicators (k-1 dummies from `{ANC}_metadata.tsv`, largest cohort as reference; skipped for single-cohort strata); `--exclude-covariates ct_Maternal` applied before \|r\| > 0.9 pruning; near-zero-variance pre-filter; pruning priority sex/cohort indicators > PCs > cell types > HCPs (pruned HCPs reduce effective k below nominal); no covariate cap. Technical-replicate sample columns (duplicate array_id) are averaged per covariate; discordant-sex replicates are excluded. `--hcp-file`/`--hcp-k`/`--out-suffix` support the 25a/25b optimization modules; the canonical per-modality run writes `{ANC}_covariates_{MOD}.tsv` |
| 25a | `25a_submit_hcp_k_jobs.sh`, `25a_optimize_hcp.sh`, `optimize_hcp_chr1.py`, `finalize_hcp_k_grid.py` | Expression-only HCP-count optimization, run after 23/24 and before canonical 25: one resumable LSF worker per k (production grid 0..100 by 5), chr1 expression mapping per k, k\* maximizing eGenes at Storey q <= 0.05 (ties -> smaller k); installs k\* as `{ANC}_hcp_factors_harmonized.tsv`; supports `LAMBDA1` sensitivity sandboxes and finalization. Every per-k model uses the same pruned fixed covariate set (`EXCLUDE_COVARIATES`, default `ct_Maternal`); HCPs lost to correlation pruning are reported (`n_hcp_used`/`n_hcp_dropped`) |
| 25b | `25b_submit_hcp_k_jobs.sh`, `25b_optimize_hcp_modalities.sh`, `optimize_hcp_modalities.py`, `hcp_from_matrix.R`, `finalize_hcp_k_grid.py` | HCP-count optimization for 8 non-expression modalities plus `combined` (9 arms per ancestry): one resumable LSF worker per k (production grid 0..100 by 5), HCP-only estimation from each harmonized BED, chr1-subset mapping (genome-wide when < 300 chr1 phenotypes), k\* = argmax eGenes at Storey q <= 0.05 (ties -> smaller k); installs `{ANC}_{MOD}_hcp_factors_optimized.tsv`; supports the same chosen `LAMBDA1` as 25a |
| 26 | `26_harmonize_modalities.py` | Harmonize the 7 non-expression modality BEDs to the final array_id sample set (handles stage-17 namespaced IDs for splicing/IR) |
| 30 | `30_combine_modalities.py` | Combined cross-modality BED (`{modality}__{id}` namespacing; cross-modality gene groups; modality sidecar TSV) |
| 27 | `27_run_tensorqtl.sh` + `27_run_tensorqtl.py` | tensorQTL `cis.map_cis` per ancestry x modality (grouped, `group_s`, when `phenotype_groups.txt` exists); `--independent` stepwise conditional mode; `--independent-only` + `--chunk-index`/`--chunk-size` run one ~100-gene chunk of the stepwise scan from the saved map_cis parquet (full `cis_df` is passed so the significance threshold matches the monolithic run; bit-identical with `--seed`); Storey q-values via the `compute_qvalues.R` file bridge (`QVALUE_RSCRIPT`), `--qvalue-method bh` as fallback. Covariates default to `{ANC}_covariates_{MOD}.tsv` (expression from 25a; non-expression/combined from 25b), falling back to `{ANC}_covariates.tsv` with a warning; `COVARIATES_FILE` accepts `{ANC}` and `{MOD}` placeholders |
| 27b/c | `27b_independent_chunk.sh`, `27c_merge_independent.py` | Chunked stepwise layer: 27b runs one array task (`--chunk-index $LSB_JOBINDEX`) -> `independent_chunks/{ANC}_{label}/chunk_KKKK.parquet` (+ `.json` parameter sidecar); 27c verifies completeness/parameters and merges chunks into the same `{ANC}_{label}_cisqtl_independent.parquet`/`_top.tsv` the monolithic mode writes |
| 28 | `28_submit_modalities.sh` | Submission driver: one LSF job per ancestry x modality (`TEST=1` pilot; threads `QVALUE_METHOD`/`MAF_THRESHOLD`). `INDEPENDENT=1` no longer runs the stepwise scan in-job; it hands off to 28b after the map_cis submissions |
| 28b | `28b_submit_independent.sh` | Chunked independent-scan submitter: sizes a `inqtl_*` job array from the number of FDR-significant rows (`CHUNK_SIZE`, default 100 genes), submits it (CPU, or `GPU=1` -> `GPU_QUEUE`/`GPU_OPTS`), and submits a dependent `inqtlm_*` merge job. Chains on running/newly-submitted map_cis jobs when the parquet is missing; `CHUNK_INDICES="3,17"` re-runs failed chunks; skip-if-done/running guards |
| 29 | `29_make_top_tables.py` | Rebuild sorted `*_cisqtl_top.tsv` from parquets (no tensorQTL rerun) |
| 29b-f | `29b_expression_diagnostics.py`, `29c_choi_comparison.py`, `29d_cohort_heterogeneity.py`, `29e_plot_diagnostics.R`, `29f_run_diagnostics.sh` | Validation gate between mapping and fine-mapping: expression sample PCs + per-cohort residual-variance table + k=15-vs-45 lead stability (29b); Choi 2024 gene-top / exact-lead / significant-set retention vs the external SNUH study (29c); per-cohort scans + Cochran Q / I2 + genotype x cohort interaction for pooled-significant pairs (29d); figures + baseline-vs-rerun validation summary (29e); LSF driver (29f) |
| 31 | `31_sushie_finemap.py` | Cross-ancestry fine-mapping: `prepare-loci` builds per-modality locus lists (union of q <= 0.05 grouped-layer lead phenotypes across ancestries; tested windows from the mapping parquets; L = min(10, max(5, n_independent + 2)) from the stepwise layer); `run` fine-maps one shard of loci jointly across ancestries with SuSHiE (individual-level mode, in-sample LD from the intersected pgens; purity 0.5; phenotypes missing from an ancestry's BED drop that ancestry for the locus; per-locus `.ancestries`/`.done` markers + per-shard diagnostics) |
| 32 | `32_submit_sushie.sh` + `32a_run_sushie_shard.sh` | LSF driver: prepare loci per modality -> shard (default 50 loci/shard) -> one array job per shard in the `sushie` conda env (`TEST=1` pilot; skip-if-done/running guards; `FORCE_LOCI`/`FORCE_RUN`; `ANCESTRIES`/`MODALITIES`/`SHARD_SIZE`/`QUEUE`/`WALLTIME`/`THREADS` env overrides) |
| 33 | `33_aggregate_finemap.py` | Aggregate per-locus SuSHiE outputs -> `finemap/aggregated/`: `finemap_pips.tsv.gz` (per-variant PIPs, CS membership, per-ancestry effect weights, locus diagnostics), `finemap_credible_sets.tsv.gz` (per-CS summaries incl. cross-ancestry rho), `finemap_locus_summary.tsv` (convergence, ELBO, n CS, max PIP per locus) |
| 34 | `34_pip_annotation_enrichment.py` | High-PIP (>= 0.9) variant annotation enrichment: `--make-fastvep-input` -> fastVEP consequences (splice/LoF/missense/synonymous/UTR/intron/regulatory/flanking/intergenic); intersects ENCODE SCREEN cCRE classes + placenta OCR BED; Fisher exact (primary), log10(distance)-adjusted logistic (sensitivity), PIP-weighted enrichment; BH-FDR per test family |

Diagnostics/utilities: `hcp_diagnostic.R`, `hcp_diagnostic2.R`,
`check_hcp_chunks.sh`, `test_hcp.py`, `test_normalize_modalities.py`,
`test_build_covariates.py`, `test_collapse_replicates.py`,
`test_optimize_hcp.py`, `test_optimize_hcp_modalities.py`,
`test_sushie_finemap.py` (end-to-end fine-mapping fixture: locus prep,
SuSHiE causal-variant recovery, aggregation, enrichment, report render).

## Technical replicates: collapse policy (optional pre-step)

58 individuals (all NIGMS — paired placental-quadrant samples) have two
sequenced runs. To collapse replicates to one column per individual before
stages 17/19/20:

```bash
COLLAPSE_DIR="$OUTPUT_BASE/replicate_collapsed"
python3 "$SCRIPTS_DIR/collapse_replicates.py" \
  --metadata /path/to/placenta_QTL_cohort_metadata.tsv \
  --ancestry-map "$ANCESTRY_MAP" \
  --config "$CONFIG" \
  --out-dir "$COLLAPSE_DIR" \
  --ref-anno "$NORMALIZED_GTF" \
  --dry-run        # first pass: reports only, no staging tree
```

Inspect `$COLLAPSE_DIR/reports/replicate_collapses.tsv` (action per pair x
modality) and `replicate_concordance.tsv` (Pearson/Spearman per pair), then
rerun without `--dry-run` to write the staging tree. Point the downstream
stages at the collapsed inputs:

```bash
export COLLAPSE_DIR="$OUTPUT_BASE/replicate_collapsed"                     # 17, 19, 20
export ANCESTRY_MAP="$COLLAPSE_DIR/reports/ancestry_map_collapsed.tsv"     # non-primary runs removed
export COLLAPSE_MAP=/path/to/placenta_QTL_cohort_metadata.tsv              # 18b deconvolution
```

Policy (defaults):

- Pairs are discovered from the pooled metadata (rnaseq_id -> array_id);
  second runs are often absent from the ancestry map.
- Primary run = the run present in the ancestry map, else lexicographically
  first. Non-primary runs are dropped from every cohort file (one global
  representative per individual, so pooling can never duplicate a person).
- Count-exact collapse per modality: raw counts are summed across runs and
  every ratio is recomputed from the summed counts (expression/isoform
  BEDs; isoform + alt_TSS/alt_polyA within-gene ratios with the full
  transcript denominator; LeafCutter numers -> within-cluster ratios;
  stability exon/intron ratio with the >= 10 count floor applied after
  summing; RNA editing (n+0.5)/(d+0.5) per site with the original row-mean
  imputation).
- Concordance gate: a pair with Spearman < 0.9 (or < 100 pairwise-complete
  features) on the cohort unnorm BED falls back to keep-primary and is
  flagged in `replicate_concordance.tsv` for review.
- Cross-protocol pairs (runs from > 1 cohort): default keep-primary.
  `--cross-protocol average` allows count-exact collapse only for the
  source-file ratio modalities (isoforms, alt_TSS/alt_polyA, stability,
  RNA editing) — never expression/isoform_expression BEDs, splicing, or IR.
- Intron retention is always keep-primary: MAJIQ PSI has no per-run counts
  to sum, so the primary run's PSI is kept and the pair is flagged
  `keep_primary_no_counts` in the report.
- Splicing and IR are staged as intermediates (numers counts / PSI tsv);
  stage 17 re-derives the pooled features from them. Untouched cohorts are
  symlinked into the staging tree unchanged.

Without collapsing, the downstream guards hold: 23 averages duplicate
sample columns after the rnaseq_id -> array_id rename, and 25 averages
duplicate covariate columns and excludes samples with discordant sex
across runs.

The ancestry map and cohort manifests may already drop second runs in
places — verify where before assuming a pair reaches a given stage.

## Step 1: HCP smoke test (5 min)

Run one HCP estimation interactively before launching the optimization
grid. Expected: a 5-row factor table whose sample columns match the BED
header count.

```bash
export R_LIBS_USER=/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1
singularity exec --bind /rsrch5 --bind /rsrch9 \
  /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif Rscript \
  "$SCRIPTS_DIR/hcp_from_matrix.R" \
  --bed "$QTL_DIR/EAS_expression_harmonized.bed" \
  --qc-metrics "$OUTPUT_BASE/hcp/all_qc_metrics.tsv" \
  --metadata "$QTL_DIR/EAS_metadata.tsv" \
  --k 5 --output /tmp/EAS_expression_hcp_smoke.tsv
```

Verify dimensions (factors = 5; samples must equal the BED's sample
count):

```bash
awk 'NR==1 {print "smoke samples:", NF-1} END {print "smoke factors:", NR}' /tmp/EAS_expression_hcp_smoke.tsv
head -1 "$QTL_DIR/EAS_expression_harmonized.bed" | awk '{print "BED samples:", NF-4}'
```

## Step 2: HCP optimization grids (25a expression; 25b non-expression + combined)

Expression is optimized separately by 25a. Non-expression modalities and
the combined phenotype are optimized by 25b. Both submitters are
resumable and use one LSF worker per `k`; completed grid points are
skipped. Use the same selected `LAMBDA1` for 25a and 25b.

Expression (25a):

```bash
export ANCESTRIES="EAS EUR"
export K_GRID="0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100"
bash 25a_submit_hcp_k_jobs.sh
```

Continue after both expression finalizers finish and
`hcp_optimization/{EAS,EUR}_optimal_hcp.tsv` plus
`{EAS,EUR}_hcp_factors_harmonized.tsv` exist.

Then harmonize all non-expression modality BEDs with script 26, build the
combined phenotype with script 30, and submit 25b:

```bash
export MODALITIES="isoform_expression alt_polyA alt_TSS intron_retention isoforms RNA_editing splicing stability combined"
bash 25b_submit_hcp_k_jobs.sh
```

A fresh production 25b grid is at most 378 workers (2 ancestries x 9
arms x 21 k values), with one finalizer per ancestry x arm. Completion
checks after all 18 finalizers:

```bash
ls "$QTL_DIR"/hcp_optimization_modalities/*_optimal_hcp.tsv | wc -l   # 18
ls "$QTL_DIR"/*_hcp_factors_optimized.tsv | wc -l                    # 18
```

`LAMBDA1` defaults to 0.5. If a non-default lambda wins the expression
sensitivity check, pass that same value through the full 25a and 25b
grids; non-default runs use suffixed sandbox trees until finalized.

## Step 3: canonical per-modality covariates (25)

Builds `{ANC}_covariates_{MOD}.tsv` for all 10 groups per ancestry.
Expression uses the winning 25a HCP set; the other eight modalities plus
combined use their 25b winners. Then copy the expression set to the
canonical `{ANC}_covariates.tsv` (kept for the report archive and other
downstream consumers). Light enough for a login node; wrap in bsub if
preferred (~1-2 min per modality).

```bash
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate tensorqtl

MODS="isoform_expression alt_polyA alt_TSS intron_retention isoforms RNA_editing splicing stability combined"

for ANC in EAS EUR; do
  # Expression winner installed by 25a.
  python3 "$SCRIPTS_DIR/25_build_covariates.py" \
    --qtl-dir "$QTL_DIR" \
    --pcair-dir "$OUTPUT_BASE/genotype_pcs" \
    --ancestries "$ANC" \
    --out-suffix "_expression" \
    --exclude-covariates ct_Maternal \
    --cor-threshold 0.9 --min-sd 1e-8

  # 25b winners: 8 non-expression modalities + combined.
  for MOD in $MODS; do
    KSTAR=$(awk '$NF=="True" {print $3}' "$QTL_DIR/hcp_optimization_modalities/${ANC}_${MOD}_optimal_hcp.tsv")
    echo "== $ANC $MOD (k*=$KSTAR) =="
    python3 "$SCRIPTS_DIR/25_build_covariates.py" \
      --qtl-dir "$QTL_DIR" \
      --pcair-dir "$OUTPUT_BASE/genotype_pcs" \
      --ancestries "$ANC" \
      --hcp-file "$QTL_DIR/${ANC}_${MOD}_hcp_factors_optimized.tsv" \
      --hcp-k "$KSTAR" \
      --out-suffix "_${MOD}" \
      --exclude-covariates ct_Maternal \
      --cor-threshold 0.9 --min-sd 1e-8
  done

  cp "$QTL_DIR/${ANC}_covariates_expression.tsv" "$QTL_DIR/${ANC}_covariates.tsv"
  cp "$QTL_DIR/${ANC}_covariate_pruning_expression.tsv" "$QTL_DIR/${ANC}_covariate_pruning.tsv"
  cp "$QTL_DIR/${ANC}_covariate_correlation_expression.png" "$QTL_DIR/${ANC}_covariate_correlation.png"
done
```

Completion check (expect 20 per-modality files; the canonical
`{ANC}_covariates.tsv` copies are not matched by this glob):

```bash
ls "$QTL_DIR"/*_covariates_*.tsv | wc -l
head -3 "$QTL_DIR/EAS_covariates_splicing.tsv" | cut -f1-3
```

With the replicate-collapse workflow, the collapsed ancestry map should
already yield one RNA run per individual by this point; duplicate-column
handling remains a defensive fallback, not the expected production path.

## Step 4: xQTL mapping (three submission passes)

`28_submit_modalities.sh` skips combos with existing results or running
jobs, so all three passes can be rerun freely. Each mapping job picks up
`{ANC}_covariates_{MOD}.tsv` automatically. Expected wall time: ~1 h per
map_cis job.

With `INDEPENDENT=1`, the stepwise (conditionally independent) scan no
longer runs inside the mapping job (that took >24 h per ancestry x
modality). Instead, 28 hands off to `28b_submit_independent.sh`, which
submits the scan as a chunked LSF job array (~`CHUNK_SIZE` genes per
chunk, default 100) plus a dependent merge job; chunk arrays chain on the
map_cis jobs automatically when the map_cis parquet is not yet on disk.
Merged outputs are the same `{ANC}_{MOD}[_ungrouped]_cisqtl_independent.*`
files as before, statistically identical to the monolithic scan (the full
map_cis table is passed to every chunk so tensorQTL's internal
significance threshold is unchanged; bit-identical with a fixed `SEED`).

```bash
cd "$SCRIPTS_DIR"

# Pass 1: grouped + chunked stepwise (independent) layer, all 9 modalities
MAF_THRESHOLD=0.01 INDEPENDENT=1 \
  MODALITIES="expression isoforms isoform_expression splicing intron_retention alt_TSS alt_polyA RNA_editing stability" \
  bash 28_submit_modalities.sh

# Pass 2: ungrouped layer (per-phenotype lead variants) + chunked stepwise
# (independent) scan, 8 non-expression modalities
MAF_THRESHOLD=0.01 GROUPED=0 INDEPENDENT=1 bash 28_submit_modalities.sh

# Pass 3: combined cross-modality arm (long queue for the map_cis job)
MAF_THRESHOLD=0.01 INDEPENDENT=1 MODALITIES=combined \
  QUEUE=long WALLTIME=48:00 bash 28_submit_modalities.sh
```

The independent layer can also be (re)submitted on its own — e.g. to reuse
the map_cis parquets of a killed monolithic run, or to use GPU nodes
(tensorQTL auto-uses CUDA; typically 10-50x faster per chunk, so a larger
`CHUNK_SIZE` is worthwhile):

```bash
# Chunked independent scans only, from existing map_cis parquets
MAF_THRESHOLD=0.01 ANCESTRIES=EAS MODALITIES="combined splicing" \
  bash 28b_submit_independent.sh

# GPU chunks on seadragon's gpu queue
MAF_THRESHOLD=0.01 GPU=1 GPU_QUEUE=gpu CHUNK_SIZE=500 \
  bash 28b_submit_independent.sh

# Retry failed chunks (the merge job's log lists the missing indices),
# then re-merge automatically
CHUNK_INDICES="3,17" MAF_THRESHOLD=0.01 MODALITIES=combined \
  bash 28b_submit_independent.sh
```

Per-chunk outputs land in `qtl_results/independent_chunks/{ANC}_{label}/`
(with `.json` parameter sidecars; the merge refuses to mix chunks computed
with different `CHUNK_SIZE`/`INDEPENDENT_FDR`). If a map_cis job chained
by 28b fails, the planner job never runs — fix the failure and rerun the
same 28b command.

Completion check (expect 36 primary + 36 independent parquets: per
ancestry, 9 grouped + 8 ungrouped + 1 combined primary, and 9 grouped +
8 ungrouped + 1 combined independent):

```bash
ls "$RESULTS_DIR"/*_cisqtl.parquet | wc -l
ls "$RESULTS_DIR"/*_cisqtl_independent.parquet | wc -l
```

If a job fails, resubmit the same pass — finished combos are skipped.
Check that no run fell back to shared covariates (should print nothing):

```bash
grep -l "falling back to shared covariates" "$LOG_DIR"/tensorqtl.*.out 2>/dev/null
```

To regenerate results from scratch, delete `$RESULTS_DIR` (genotypes,
harmonized BEDs, and metadata under `qtl_inputs` are reused as-is). The
existing `{ANC}_{MOD}.bed.gz` + `.tbi` files in `qtl_inputs` remain valid
and are reused by `27_run_tensorqtl.sh`:

```bash
( set -e
  case "$RESULTS_DIR" in
    /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/qtl_results) ;;
    *) echo "REFUSING: unexpected RESULTS_DIR=$RESULTS_DIR"; exit 1 ;;
  esac
  echo "Deleting $RESULTS_DIR ($(ls "$RESULTS_DIR" | wc -l) entries)"
  rm -rf "$RESULTS_DIR"
  mkdir -p "$RESULTS_DIR"
  echo "Done." )
```

## Step 5: top tables + report extras + report-input archive

```bash
python3 "$SCRIPTS_DIR/29_make_top_tables.py" --results-dir "$RESULTS_DIR"

# Targeted extractions for the DevBrain-style report panels (cross-ancestry
# and cross-modality lead-pair lookups, showcase-locus scans, lead-variant
# annotations, gnomAD constraint, PC-AiR PC staging). Single node, < 1 h.
# The mapping parquets are lead-only (cis.map_cis), so these panels need
# targeted re-computation rather than file copies.
# Compute nodes have no internet access: pre-download the gnomAD constraint
# file (~30 MB) on a login node into the script's cache path first:
#   wget -q https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
#     -O "$RESULTS_DIR/report_extras/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
python3 "$SCRIPTS_DIR/35_extract_report_extras.py" \
  --results-dir "$RESULTS_DIR" --qtl-dir "$QTL_DIR" \
  --ccre $CRE_DIR/GRCh38-cCREs.CA.bed:CA $CRE_DIR/GRCh38-cCREs.CA-CTCF.bed:CA_CTCF \
         $CRE_DIR/GRCh38-cCREs.CA-H3K4me3.bed:CA_H3K4me3 $CRE_DIR/GRCh38-cCREs.CA-TF.bed:CA_TF \
         $CRE_DIR/GRCh38-cCREs.CTCF-bound.bed:CTCF_bound $CRE_DIR/GRCh38-cCREs.dELS.bed:dELS \
         $CRE_DIR/GRCh38-cCREs.pELS.bed:pELS $CRE_DIR/GRCh38-cCREs.PLS.bed:PLS \
         $CRE_DIR/GRCh38-cCREs.TF.bed:TF \
  --placenta-ocr $ENCODE_DIR/placenta_dnase_merged_bothstrands.bed \
  --pcair-dir "$OUTPUT_BASE/pcair_pcs"   # staging dir of symlinked
  # *_pcair_pcs.tsv — PC_AiR.R writes per-cohort under
  # $PCAIR_SRC/{COHORT}/genotypes/imputed/qc/; symlink them into one dir:
  #   mkdir -p "$OUTPUT_BASE/pcair_pcs"
  #   for c in NIEHS_RICHS GUSTO SNUH NIGMS; do
  #     ln -sf "/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/$c/genotypes/imputed/qc/${c}_pcair_pcs.tsv" \
  #       "$OUTPUT_BASE/pcair_pcs/"
  #   done
# then annotate the emitted VCF with fastVEP and re-run section 4:
#   fastvep annotate -i "$RESULTS_DIR/report_extras/lead_variants_fastvep_input.vcf" \
#     -o "$RESULTS_DIR/report_extras/lead_variants_fastvep.txt" --output-format tab \
#     --gff3 "$REF_DIR/Homo_sapiens.GRCh38.115.gff3" --fasta "$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
#   python3 "$SCRIPTS_DIR/35_extract_report_extras.py" --results-dir "$RESULTS_DIR" \
#     --qtl-dir "$QTL_DIR" --only annotate --fastvep-output "$RESULTS_DIR/report_extras/lead_variants_fastvep.txt" ...

# QC summaries for the report's processing-QC section (before/after
# phenotype violins + PCAs per modality, genotype stage counts, Picard
# metrics). Single node, < 30 min. All inputs optional (graceful skips);
# see reports/runbook_report_extras_seadragon.md step 5 for details.
python3 "$SCRIPTS_DIR/35_extract_report_extras.py" \
  --results-dir "$RESULTS_DIR" --qtl-dir "$QTL_DIR" --only qc \
  --pooled-bed-dir "$OUTPUT_BASE/normalized_modalities/pooled" \
  --ancestry-map "/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/pooled/pooled_sample_ancestry_RNAseq.tsv" \
  --geno-ancestry-map "/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/pooled/genotypes/pooled_sample_ancestry.tsv" \
  --cohort-qc-glob "/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/*/genotypes/imputed/qc" \
  --geno-pooled-dir "/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/pooled/genotypes" \
  --picard-glob "$OUTPUT_BASE/hcp/qc_metrics/*_qc_metrics.tsv"

cd "$REPO_DIR"
# Optional: GENO_QC_DIR stages module-01 genotype-QC tables (symlink the
# pooled/genotypes/report pooled_* tables and the per-cohort
# genotypes/imputed/qc/report/*_rsq_pass_per_chr.tsv into one dir first);
# POOLED_ANCESTRY_TSV stages the genotype sample->ancestry map.
CONFIG="$CONFIG" GENO_QC_DIR="$OUTPUT_BASE/geno_qc_stage" \
  bash ../reports/make_report_archive.sh
```

The archive (`placenta_xqtl_report_inputs_<date>.tar.gz`, written to the
repo root) includes `data/qc/hcp_optimization_modalities/` with the
per-modality k\* tables and curves alongside `data/qc/hcp_optimization/`.
Review the MANIFEST.txt it prints: no core inputs may be missing.

## Step 6: cross-ancestry fine-mapping (SuSHiE)

SuSHiE joint multi-ancestry fine-mapping of every FDR <= 5% lead
phenotype from the **phenotype-level discovery layer** — the ungrouped
scan for the 8 non-expression modalities, and the primary (per-gene)
scan for expression — using in-sample LD from the intersected pgens
(individual-level mode — LD and genotypes are guaranteed consistent).
cis windows are reconstructed from the phenotype BEDs (ungrouped
parquets do not carry window bounds), and contig names are stripped of
any `chr` prefix for SuSHiE `--chrom`. Each locus is fine-mapped jointly
across all ancestries in `ANCESTRIES`, regardless of which ancestry
reached significance. Outputs: 95% credible
sets, per-variant PIPs, per-ancestry effect weights, and cross-ancestry
effect-size correlations (rho) per credible set.

### 6.0 One-time setup

SuSHiE requires python > 3.11 and runs in its own conda env:

```bash
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda create -n sushie python=3.11
conda activate sushie
pip install sushie pandas pyarrow fastparquet
sushie finemap --help   # smoke test
```

Annotation sources for the enrichment step (6.4):

- **fastVEP** (variant consequences): Rust reimplementation of Ensembl
  VEP; no multi-GB cache download. On seadragon:
  `conda create -n fastvep -c conda-forge rust c-compiler`, then
  `cargo install --path crates/fastvep-cli --root "$CONDA_PREFIX"` from a
  clone of https://github.com/Huang-lab/fastVEP at tag v0.4.0.
  Requires two reference files (both bare `1`-style contigs):
  Ensembl release-115 GFF3 (`Homo_sapiens.GRCh38.115.gff3`, ~50 MB from
  ftp.ensembl.org) and the GRCh38 primary-assembly FASTA with `.fai`
  (already at `$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa`).
  fastVEP v0.4.0 does not annotate Ensembl regulatory-build features, so
  the `consequence_regulatory_region` class is expected to be empty.
  Regulatory signal is captured by the cCRE/OCR intersections below.
- **ENCODE SCREEN cCREs**: download GRCh38 "cCREs by class" BEDs from
  https://screen.wenglab.org/downloads — one file per class: PLS, pELS,
  dELS, CTCF-bound, CA-TF (e.g. `GRCh38-cCREs.PLS.bed.gz`).
- **Placenta open chromatin**: query the ENCODE portal (type=Experiment,
  assay = DNase-seq or ATAC-seq, biosample = placenta, assembly GRCh38,
  file type = bed narrowPeak). Prefer replicated / IDR-thresholded peak
  calls per the ENCODE4 standards. No single canonical accession — pick
  the experiment(s) that best match the cohort's gestational-age range
  and record the accession(s) used.

### 6.1 Pilot (measure throughput before the full run)

```bash
cd "$SCRIPTS_DIR"
TEST=1 bash 32_submit_sushie.sh
```

Builds per-modality locus lists (`$RESULTS_DIR/finemap/loci/{MOD}_loci.tsv`:
union of q <= 0.05 lead phenotypes across ancestries from the
phenotype-level discovery layer, tested windows reconstructed from the
phenotype BEDs, L = min(10, max(5, n_independent + 2))), shards them
(50 loci/shard), and submits ONE pilot shard. Check `$LOG_DIR/sushie_*.out`
for per-locus wall time to size `SHARD_SIZE` / `WALLTIME` for the full run
(fixture rate ~ 17 s per 300-variant locus).

### 6.2 Full submission

```bash
bash 32_submit_sushie.sh          # all 9 modalities, EAS+EUR
```

Skip-if-done / skip-if-running guards make reruns safe. `FORCE_RUN=1`
re-runs loci with existing `.done` markers; `FORCE_LOCI=1` rebuilds the
locus lists. Per-locus outputs land in
`$RESULTS_DIR/finemap/{MOD}/{phenotype_id}/` (SuSHiE
`.sushie.weights.tsv`, `.sushie.cs.tsv`, `.sushie.corr.tsv`, `.log`, plus
`.ancestries` and `.done` markers); per-shard diagnostics in
`$RESULTS_DIR/finemap/{MOD}/logs/`.

### 6.3 Aggregate

```bash
python3 "$SCRIPTS_DIR/33_aggregate_finemap.py" \
  --finemap-dir "$RESULTS_DIR/finemap" --ancestries EAS EUR
```

Writes to `$RESULTS_DIR/finemap/aggregated/`:
`finemap_pips.tsv.gz` (per-variant PIPs + CS membership + per-ancestry
effect weights), `finemap_credible_sets.tsv.gz` (per-CS summaries incl.
cross-ancestry rho), `finemap_locus_summary.tsv` (per-locus diagnostics:
convergence, ELBO, n CS, max PIP).

### 6.4 Annotation enrichment

Install fastVEP

```bash
conda create -n fastvep -c conda-forge -y rust c-compiler
conda activate fastvep
cargo install --path crates/fastvep-cli --root "$CONDA_PREFIX"
conda install pip
pip install pandas scipy
```

```bash
AGG="$RESULTS_DIR/finemap/aggregated"
GENOME_DIR=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/GenomicReferences/genome
REF_DIR=/rsrch9/home/epi/bhattacharya_lab/data/fastVEP/references
mkdir -p $REF_DIR
cd $REF_DIR

wget https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz && gunzip Homo_sapiens.GRCh38.115.gff3.gz

mkdir -p /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/qtl_results/finemap_annotate
cd /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/qtl_results/finemap_annotate

# 1. fastVEP input VCF (one row per unique fine-mapped variant)
python3 "$SCRIPTS_DIR/34_pip_annotation_enrichment.py" \
  --pips "$AGG/finemap_pips.tsv.gz" --make-fastvep-input --out fastvep_input.vcf

# 2. Run fastVEP (command printed by the previous step; seconds-to-minutes,
#    single node, <1 GB RAM — no batch job needed)
fastvep annotate -i fastvep_input.vcf -o fastvep_output.txt \
  --output-format tab \
  --gff3 "$REF_DIR/Homo_sapiens.GRCh38.115.gff3" \
  --fasta "$GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# 3. Enrichment (high-PIP >= 0.9 vs all fine-mapped variants as background)
CRE_DIR=/rsrch9/home/epi/bhattacharya_lab/data/functional_annotations/placenta/cCREs
ENCODE_DIR=/rsrch9/home/epi/bhattacharya_lab/data/functional_annotations/placenta/ENCODE

python3 "$SCRIPTS_DIR/34_pip_annotation_enrichment.py" \
  --pips "$AGG/finemap_pips.tsv.gz" \
  --fastvep fastvep_output.txt \
  --ccre $CRE_DIR/GRCh38-cCREs.CA.bed:CA \
         $CRE_DIR/GRCh38-cCREs.CA-CTCF.bed:CA_CTCF \
         $CRE_DIR/GRCh38-cCREs.CA-H3K4me3.bed:CA_H3K4me3 \
         $CRE_DIR/GRCh38-cCREs.CA-TF.bed:CA_TF \
         $CRE_DIR/GRCh38-cCREs.CTCF-bound.bed:CTCF_bound \
         $CRE_DIR/GRCh38-cCREs.dELS.bed:dELS \
         $CRE_DIR/GRCh38-cCREs.pELS.bed:pELS \
         $CRE_DIR/GRCh38-cCREs.PLS.bed:PLS \
         $CRE_DIR/GRCh38-cCREs.TF.bed:TF \
  --placenta-ocr $ENCODE_DIR/placenta_dnase_merged_bothstrands.bed \
  --out "$AGG/finemap_enrichment.tsv"
```

Per annotation class: Fisher exact test (primary), logistic regression
adjusting for log10(distance to phenotype start) (sensitivity), and a
PIP-weighted enrichment; BH-FDR within each test family.

### 6.5 Diagnostic report

```bash
cd /rsrch5/home/epi/stbresnahan/bhattacharya_lab/software/Pantry/phenotyping/scripts
singularity exec --bind /rsrch5 --bind /rsrch9 \
  /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif \
  Rscript -e '.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths())); rmarkdown::render(
    "reports/report_finemap.Rmd",
    params=list(
      agg_dir="'$RESULTS_DIR'/finemap/aggregated",
      enrichment_path="'$RESULTS_DIR'/finemap/aggregated/finemap_enrichment.tsv",
      fig_dir="reports/fig_finemap",
      table_dir="reports/tables_finemap"
    )
  )'
```

Sections: overview, max-PIP-per-locus and PIP distributions, credible-set
sizes and CSs per locus, cross-ancestry rho, per-ancestry effect-weight
concordance, top-locus PIP tracks, enrichment forest plot, run
diagnostics.

### 6.6 Fallback

Loci with `status=failed` or non-converged SuSHiE runs (see
`finemap_locus_summary.tsv`) are candidates for single-ancestry FINEMAP,
or restricting to the strongest single signal. Not implemented here —
flag and handle case-by-case.

## Operational notes

- **Storey q-values use the singularity-R bridge.** rpy2/conda-R fails on
  seadragon with a `GLIBCXX_3.4.30 not found` libstdc++ loader conflict;
  `compute_qvalues.R` (file bridge, `QVALUE_RSCRIPT`) is statistically
  identical to tensorQTL's `calculate_qvalues` (itself an rpy2 wrapper on
  `qvalue::qvalue`). The `'rfunc' cannot be imported` warning at
  `import tensorqtl` is cosmetic.
- **`27_run_tensorqtl.sh` skips bgzip/tabix if `{ANC}_{modality}.bed.gz`
  exists** — stale harmonized BEDs are silently reused unless the
  `.bed.gz`/`.tbi` files are deleted before a rerun.
- **`24_outlier_exclusion.py` edits the intersection files in place** —
  always rerun 23 before 24.
- **Covariate exclusions are applied in code** — `ct_Maternal` via
  `25_build_covariates.py --exclude-covariates` (and the
  `EXCLUDE_COVARIATES` env var in 25a) before correlation pruning; there
  is no covariate cap.
- **25a installs the expression HCP winner and 25b installs the
  non-expression/combined winners.** `optimize_hcp_chr1.py` / the 25a
  finalizer installs `{ANC}_hcp_factors_harmonized.tsv`; 25b installs
  `{ANC}_{MOD}_hcp_factors_optimized.tsv` for the eight non-expression
  modalities plus `combined`. Always rebuild canonical covariates after
  HCP optimization so each `{ANC}_covariates_{MOD}.tsv` matches its
  installed k\*.
- **`picard` must resolve after 19/19a's env stack** — `picard_qc.py`
  calls the literal `picard` executable; if the `picard-2.27.4`
  activation does not put it on PATH (wrong env name/path, or a
  PATH-clobbering venv activation), every QC call WARNs `[Errno 2] No
  such file or directory: 'picard'` and writes garbage metrics. Both
  wrappers fail fast with diagnostics — fix the activation. Stage 19
  skips existing `{COHORT}_qc_metrics.tsv`, so delete tainted QC outputs
  before resubmitting.
- **Picard QC extracts the full MultiQC-style field set** (all rows of
  all five collectors; primary row = `CATEGORY=PAIR` else `UNPAIRED`
  else first, unsuffixed; other rows suffixed e.g. `.FIRST_OF_PAIR`).
  The 24 legacy column names are unchanged, but `AlignMetrics.*` values
  now come from the library-level PAIR row (the legacy whitelist read
  the FIRST_OF_PAIR row for paired-end data). Raw Picard outputs persist
  under `hcp/qc_metrics/raw/{COHORT}/`; `picard_qc.py --parse-only`
  (or `PARSE_ONLY=1` in 19a) re-extracts metrics from cached raw output
  without re-running Picard.
- **Pooling takes the feature union within each ancestry stratum**
  (`--pool-mode union`, default, in `pool_expression_within_ancestry.py`
  and `pool_modalities_within_ancestry.py`): per-cohort BEDs are
  reindexed to the union feature set and NaN-filled where a cohort lacks
  a feature, so per-cohort prefilters no longer decide the pooled
  feature set — the NaN-aware pooled detection filters at the
  normalization stage do. `--pool-mode intersection` restores the legacy
  behavior. Pooling summary TSVs report
  intersection/union/gained/NaN-cell counts per stratum.
- **`bsub < script` requires an explicit `SCRIPTS_DIR` export** — LSF
  executes a spool copy of the submitted script, so
  `${BASH_SOURCE[0]}` self-location resolves to the spool directory. The
  wrappers fail fast with a clear error in that case. Direct invocation
  (`bash script.sh`, interactive) self-locates and needs no SCRIPTS_DIR.
  Cross-directory references resolve relative to the script:
  `config_get.py` is sought in `SCRIPTS_DIR`, then `../03_phenotyping`
  (shell wrappers resolve it once as `CONFIG_GET`; the python callers
  `picard_qc.py`, `pool_expression_within_ancestry.py`,
  `pool_modalities_within_ancestry.py`, and `23_prepare_intersection.py`
  apply the same chain, with a cwd/PATH fallback last); `picard_qc.py`
  (regenerate_refflat.sh) in `SCRIPTS_DIR` then `../05_qtl_mapping`.
- **The design is ancestry-agnostic**: when AFR/AMR/SAS inputs land, add
  the labels to the ANC loops in steps 2-4 (and `ANCESTRIES` in 28 and
  32) with no code changes. SuSHiE then fine-maps jointly across all
  listed ancestries.

## Inputs / outputs

**Inputs** (from earlier stages): `{ANC}_pooled.pgen` (01), modality BEDs
+ deconvolution proportions (03), `rnaseq_to_array_id_map.csv` (04),
cohort metadata (`rnaseq_id, array_id, ancestry, cohort, sex, GA, ppBMI`).

**Outputs**: `$QTL_DIR` (intersected pgens, `{ANC}_genotype_pcs.tsv`,
`{ANC}_covariates.tsv`, harmonized BEDs) and `$RESULTS_DIR`
(`{ANC}_{modality}[_ungrouped]_cisqtl.parquet` + `_top.tsv` per layer —
grouped / ungrouped / combined — plus `_independent` stepwise outputs;
`finemap/` per-locus SuSHiE outputs + `aggregated/` PIP, credible-set,
locus-summary, and enrichment tables).

Next: [../reports/README.md](../reports/README.md).

---

*Schema-change notes: on 2026-10-02, ComBat batch correction was dropped
in favor of cohort indicator covariates. On 2026-10-03, count modalities
(gene and isoform expression) moved to pooled TMM -> VST, while ratio
modalities retained pooled QN -> within-cohort rank INT. The last commit
using the ComBat schema is `ea2b344ab9cd7cdf116b742d97f072d3f7597186`
(main, 2026-10-02, "Remove stale files").*
