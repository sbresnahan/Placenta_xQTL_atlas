# Runbook: per-modality HCP optimization + xQTL remapping

End-to-end procedure for seadragon: optimize HCP hidden-covariate count k
separately for each modality and the combined cross-modality arm (module 25b),
rebuild covariates per modality, rerun all mapping layers with the matched
covariate sets, and regenerate the report-input archive.

Each tensorQTL run uses the covariate set optimized for exactly the matrix it
maps: 10 HCP sets per ancestry (9 modalities + combined). HCPs are estimated
from the final harmonized BEDs (`{ANC}_{MOD}_harmonized.bed`, already
QN+INT+ComBat'd), so phenotype matrices are unchanged — only HCPs, covariates,
and mapping outputs are regenerated. k* per group maximizes cis-eGenes
(Storey q <= 0.05) on a chr1 subset, or genome-wide for modalities with
< 300 chr1 phenotypes (RNA_editing). HCP estimation for the combined arm uses
a deterministic 40k-phenotype subsample to bound cost (estimation only;
mapping always uses the full matrix).

Expected wall time: ~14 h for the longest optimization job (combined arm),
~1 day for mapping. Everything is resumable: rerun any failed piece with
`SKIP_EXISTING=1` (25b) or by resubmitting (28 skips finished combos).

## 0. Setup

```bash
CONFIG=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/config.yml
OUTPUT_BASE=/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY
QTL_DIR="$OUTPUT_BASE/qtl_inputs"
RESULTS_DIR="$OUTPUT_BASE/qtl_results"
REPO_DIR=/rsrch5/home/epi/stbresnahan/bhattacharya_lab/software/Pantry/phenotyping/scripts
SCRIPTS_DIR="$REPO_DIR/05_qtl_mapping"
LOG_DIR=/rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs
mkdir -p "$LOG_DIR"
```

Deploy the current repo (pick one):

```bash
cd "$REPO_DIR" && git pull
# or, from the delivered tarball:
# tar xzf /path/to/Placenta_xQTL_atlas_patched.tar.gz -C "$REPO_DIR" --strip-components=1
```

Delete the existing mapping results (they are regenerated below; nothing
else under qtl_inputs is touched — genotypes, harmonized BEDs, and metadata
are reused as-is):

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

The existing `{ANC}_{MOD}.bed.gz` + `.tbi` files in qtl_inputs are still valid
(the harmonized BEDs are unchanged), so 27_run_tensorqtl.sh will reuse them.

## Technical replicates: collapse policy (optional pre-step for future runs)

58 individuals each have two sequenced runs (all in the NIGMS cohort — paired
placental-quadrant samples). The EAS/EUR results in the current report were
mapped WITHOUT collapsing: both runs entered ComBat and the QTL inputs, and
the duplicate individual columns were averaged downstream (scripts 23/25/26).
That is why run-level and individual-level counts differ in the report (e.g.,
EUR 145 individuals at ComBat vs 161 runs at QTL input). For AFR/AMR/SAS — and
any regeneration of EAS/EUR — collapse replicates to one column per individual
BEFORE stages 17/19/20:

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
re-run without `--dry-run` to write the staging tree. Point the downstream
stages at the collapsed inputs:

```bash
export COLLAPSE_DIR="$OUTPUT_BASE/replicate_collapsed"                          # 17, 19, 20
export ANCESTRY_MAP="$COLLAPSE_DIR/reports/ancestry_map_collapsed.tsv"          # non-primary runs removed
export COLLAPSE_MAP=/path/to/placenta_QTL_cohort_metadata.tsv                   # 18b deconvolution
```

Policy summary (defaults):

- Pairs are discovered from the pooled metadata (rnaseq_id -> array_id), NOT
  the ancestry map — second runs are often absent from the ancestry map.
- Primary run = the run present in the ancestry map, else lexicographically
  first. Non-primary runs are dropped from every cohort file (one global
  representative per individual, so pooling can never duplicate a person).
- Count-exact collapse per modality: raw counts are SUMMED across runs and
  every ratio is recomputed from the summed counts (expression/isoform BEDs;
  isoform + alt_TSS/alt_polyA within-gene ratios with the full transcript
  denominator; LeafCutter numers -> within-cluster ratios; stability
  exon/intron ratio with the >=10 count floor applied AFTER summing; RNA
  editing (n+0.5)/(d+0.5) per site with the original row-mean imputation).
- Concordance gate: a pair with Spearman < 0.9 (or < 100 pairwise-complete
  features) on the cohort unnorm BED falls back to keep-primary and is
  flagged in `replicate_concordance.tsv` for review.
- Cross-protocol pairs (runs from > 1 cohort): default keep-primary.
  `--cross-protocol average` allows count-exact collapse only for the
  source-file ratio modalities (isoforms, alt_TSS/alt_polyA, stability,
  RNA editing) — never expression/isoform_expression BEDs, splicing, or IR.
- CAVEAT — intron_retention is ALWAYS keep-primary: MAJIQ PSI has no
  per-run counts to sum, so the primary run's PSI is kept and the pair is
  flagged `keep_primary_no_counts` in the report.
- Splicing and IR are staged as intermediates (numers counts / PSI tsv);
  stage 17 re-derives the pooled features from them. Untouched cohorts are
  symlinked into the staging tree unchanged.

Even without collapsing, the downstream guards hold: 23 averages duplicate
sample columns after the rnaseq_id -> array_id rename, and 25 averages
duplicate covariate columns and EXCLUDES samples with discordant sex across
runs (previously warn-only).

Seadragon-side check: the ancestry map and cohort manifests may already drop
second runs in places — verify where before assuming a pair reaches a given
stage.

## 1. Smoke test: hcp_from_matrix.R (5 min)

Run one HCP estimation interactively before launching the grid. Expected:
a 5-row factor table whose sample columns match the BED header count.

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

Verify dimensions (factors = 5; samples must equal the BED's sample count):

```bash
awk 'NR==1 {print "smoke samples:", NF-1} END {print "smoke factors:", NR}' /tmp/EAS_expression_hcp_smoke.tsv
head -1 "$QTL_DIR/EAS_expression_harmonized.bed" | awk '{print "BED samples:", NF-4}'
```

## 2. HCP optimization grid (20 LSF jobs)

One job per ancestry x modality group. isoform_expression and combined go to
the long queue; the rest fit in medium.

```bash
for ANC in EAS EUR; do
  for MOD in expression isoforms isoform_expression splicing intron_retention alt_TSS alt_polyA RNA_editing stability combined; do
    Q=medium; W=24:00
    case "$MOD" in isoform_expression|combined) Q=long; W=48:00 ;; esac
    bsub -J "hcpopt_${ANC}_${MOD}" -q "$Q" -n 4 -M 32 -R "rusage[mem=32]" -W "$W" \
      -o "$LOG_DIR/hcpopt_${ANC}_${MOD}.%J.out" -e "$LOG_DIR/hcpopt_${ANC}_${MOD}.%J.err" \
      -env "CONFIG=$CONFIG,SCRIPTS_DIR=$SCRIPTS_DIR,ANCESTRIES=$ANC,MODALITIES=$MOD" \
      < "$SCRIPTS_DIR/25b_optimize_hcp_modalities.sh"
  done
done
```

Monitor:

```bash
bjobs -w | grep hcpopt
```

If any job fails (check the `.err` log), fix the cause and resubmit that pair
with SKIP_EXISTING=1 appended to the -env list — completed k grid points are
reused:

```bash
# example: resubmit EAS splicing
bsub -J hcpopt_EAS_splicing -q medium -n 4 -M 32 -R "rusage[mem=32]" -W 24:00 \
  -o "$LOG_DIR/hcpopt_EAS_splicing.%J.out" -e "$LOG_DIR/hcpopt_EAS_splicing.%J.err" \
  -env "CONFIG=$CONFIG,SCRIPTS_DIR=$SCRIPTS_DIR,ANCESTRIES=EAS,MODALITIES=splicing,SKIP_EXISTING=1" \
  < "$SCRIPTS_DIR/25b_optimize_hcp_modalities.sh"
```

Completion checks (both counts must be 20):

```bash
ls "$QTL_DIR"/hcp_optimization_modalities/*_optimal_hcp.tsv | wc -l
ls "$QTL_DIR"/*_hcp_factors_optimized.tsv | wc -l
```

Chosen k* per ancestry x modality:

```bash
for f in "$QTL_DIR"/hcp_optimization_modalities/*_optimal_hcp.tsv; do
  awk -v f="$(basename "$f" _optimal_hcp.tsv)" '$NF=="True" {print f, "k*="$3, "eGenes="$4, "scope="$6}' "$f"
done
```

## 3. Canonical per-modality covariates (~1-2 min per modality)

Builds `{ANC}_covariates_{MOD}.tsv` for all 10 groups per ancestry from the
installed k* HCP sets, then copies the expression set to the canonical
`{ANC}_covariates.tsv` (kept for the report archive and other downstream
consumers). Light enough for a login node; wrap in bsub if preferred.

```bash
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate tensorqtl

for ANC in EAS EUR; do
  for MOD in expression isoforms isoform_expression splicing intron_retention alt_TSS alt_polyA RNA_editing stability combined; do
    KSTAR=$(awk '$NF=="True" {print $3}' "$QTL_DIR/hcp_optimization_modalities/${ANC}_${MOD}_optimal_hcp.tsv")
    echo "== $ANC $MOD (k*=$KSTAR) =="
    python3 "$SCRIPTS_DIR/25_build_covariates.py" \
      --qtl-dir "$QTL_DIR" \
      --pcair-dir "$OUTPUT_BASE/genotype_pcs" \
      --ancestries "$ANC" \
      --hcp-file "$QTL_DIR/${ANC}_${MOD}_hcp_factors_optimized.tsv" \
      --hcp-k "$KSTAR" \
      --out-suffix "_${MOD}" \
      --exclude-covariates ct_Maternal
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

Expected log output: each 25 run prints a NOTE line listing the
technical-replicate samples it averages (4 EAS samples: SRR13696945,
SRR13696964, SRR13696996, SRR13697029 — the same individuals processed in two
cohort batches). A `WARN: discordant sex across technical replicates` line
would instead indicate a metadata inconsistency — investigate before
proceeding if it appears.

## 4. xQTL mapping (three submission passes)

28_submit_modalities.sh skips combos with existing results or running jobs,
so all three passes can be rerun freely. Each mapping job picks up
`{ANC}_covariates_{MOD}.tsv` automatically.

```bash
cd "$SCRIPTS_DIR"

# Pass 1: grouped + stepwise (independent) layer, all 9 modalities
MAF_THRESHOLD=0.01 INDEPENDENT=1 \
  MODALITIES="expression isoforms isoform_expression splicing intron_retention alt_TSS alt_polyA RNA_editing stability" \
  bash 28_submit_modalities.sh

# Pass 2: ungrouped layer (per-phenotype lead variants), 8 non-expression modalities
MAF_THRESHOLD=0.01 GROUPED=0 bash 28_submit_modalities.sh

# Pass 3: combined cross-modality arm (long queue)
MAF_THRESHOLD=0.01 INDEPENDENT=1 MODALITIES=combined \
  QUEUE=long WALLTIME=48:00 bash 28_submit_modalities.sh
```

Completion check (expect 36 primary + 20 independent parquets: per ancestry,
9 grouped + 8 ungrouped + 1 combined primary, 9 grouped-independent +
1 combined-independent):

```bash
ls "$RESULTS_DIR"/*_cisqtl.parquet | wc -l
ls "$RESULTS_DIR"/*_cisqtl_independent.parquet | wc -l
```

If a job fails, resubmit the same pass — finished combos are skipped. Check
that no run fell back to shared covariates (should print nothing):

```bash
grep -l "falling back to shared covariates" "$LOG_DIR"/tensorqtl.*.out 2>/dev/null
```

## 5. Top tables + report-input archive

```bash
python3 "$SCRIPTS_DIR/29_make_top_tables.py" --results-dir "$RESULTS_DIR"

cd "$REPO_DIR"
CONFIG="$CONFIG" bash reports/make_report_archive.sh
```

The archive (`placenta_xqtl_report_inputs_<date>.tar.gz`, written to the
repo root) now includes `data/qc/hcp_optimization_modalities/` with the
per-modality k* tables and curves alongside the existing
`data/qc/hcp_optimization/`. Review the MANIFEST.txt it prints: no core
inputs may be missing.

## 6. Handoff

Upload `placenta_xqtl_report_inputs_<date>.tar.gz` back to Biomni for the
report revision (per-modality k* table/figure, updated covariate QC, and
refreshed results sections).

## 7. Cross-ancestry fine-mapping (Objective 1.5)

SuSHiE joint multi-ancestry fine-mapping of every grouped-layer lead
phenotype at FDR ≤ 5%, using in-sample LD from the intersected pgens
(individual-level mode — LD and genotypes are guaranteed consistent).
Each locus is fine-mapped jointly across all ancestries in `ANCESTRIES`,
regardless of which ancestry reached significance. Outputs: 95% credible
sets, per-variant PIPs, per-ancestry effect weights, and cross-ancestry
effect-size correlations (rho) per credible set.

### 7.0 One-time setup

*Note for future cleanup: sushie requires python > 3.11; it requires
its own conda env. All sushie-related scripts have been modified to
use `conda activate sushie` instead.*

```bash
# SuSHiE into the tensorqtl env (same env scripts 27/28 use)
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate tensorqtl
pip install sushie
sushie finemap --help   # smoke test
```

Annotation sources for the enrichment step (7.4):

- **VEP** (variant consequences): assumed available on seadragon
  (`vep --help` to confirm; needs a GRCh38 offline cache).
- **ENCODE SCREEN cCREs**: download GRCh38 "cCREs by class" BEDs from
  https://screen.wenglab.org/downloads — one file per class: PLS, pELS,
  dELS, CTCF-bound, CA-TF (e.g. `GRCh38-cCREs.PLS.bed.gz`).
- **Placenta open chromatin**: query the ENCODE portal
  (type=Experiment, assay = DNase-seq or ATAC-seq, biosample = placenta,
  assembly GRCh38, file type = bed narrowPeak). Prefer replicated /
  IDR-thresholded peak calls per the ENCODE4 standards. No single
  canonical accession — pick the experiment(s) that best match the
  cohort's gestational-age range and record the accession(s) used.

### 7.1 Pilot (measure throughput before the full run)

```bash
cd "$SCRIPTS_DIR"
TEST=1 bash 32_submit_sushie.sh
```

Builds per-modality locus lists (`$RESULTS_DIR/finemap/loci/{MOD}_loci.tsv`:
union of q ≤ 0.05 lead phenotypes across ancestries, exact tested windows
from the mapping parquets, L = min(10, max(5, n_independent + 2))), shards
them (50 loci/shard), and submits ONE pilot shard. Check
`$LOG_DIR/sushie_*.out` for per-locus wall time to size `SHARD_SIZE` /
`WALLTIME` for the full run (fixture rate ≈ 17 s per 300-variant locus).

### 7.2 Full submission

```bash
bash 32_submit_sushie.sh          # all 9 modalities, EAS+EUR
```

Skip-if-done / skip-if-running guards make reruns safe. `FORCE_RUN=1`
re-runs loci with existing `.done` markers; `FORCE_LOCI=1` rebuilds the
locus lists. Per-locus outputs land in
`$RESULTS_DIR/finemap/{MOD}/{phenotype_id}/` (SuSHiE `.sushie.weights.tsv`,
`.sushie.cs.tsv`, `.sushie.corr.tsv`, `.log`, plus `.ancestries` and
`.done` markers); per-shard diagnostics in `$RESULTS_DIR/finemap/{MOD}/logs/`.

### 7.3 Aggregate

```bash
python3 "$SCRIPTS_DIR/33_aggregate_finemap.py" \
  --finemap-dir "$RESULTS_DIR/finemap" --ancestries EAS EUR
```

Writes to `$RESULTS_DIR/finemap/aggregated/`:
`finemap_pips.tsv.gz` (per-variant PIPs + CS membership + per-ancestry
effect weights), `finemap_credible_sets.tsv.gz` (per-CS summaries incl.
cross-ancestry rho), `finemap_locus_summary.tsv` (per-locus diagnostics:
convergence, ELBO, n CS, max PIP).

### 7.4 Annotation enrichment

```bash
AGG="$RESULTS_DIR/finemap/aggregated"

# 1. VEP input (one row per unique fine-mapped variant)
python3 "$SCRIPTS_DIR/34_pip_annotation_enrichment.py" \
  --pips "$AGG/finemap_pips.tsv.gz" --make-vep-input --out vep_input.tsv

# 2. Run VEP (command printed by the previous step)
vep -i vep_input.tsv --cache --offline --assembly GRCh38 \
  --output_file vep_output.txt --force_overwrite

# 3. Enrichment (high-PIP ≥ 0.9 vs all fine-mapped variants as background)
python3 "$SCRIPTS_DIR/34_pip_annotation_enrichment.py" \
  --pips "$AGG/finemap_pips.tsv.gz" \
  --vep vep_output.txt \
  --ccre GRCh38-cCREs.PLS.bed.gz:PLS GRCh38-cCREs.pELS.bed.gz:pELS \
         GRCh38-cCREs.dELS.bed.gz:dELS GRCh38-cCREs.CTCF-bound.bed.gz:CTCF_bound \
         GRCh38-cCREs.CA-TF.bed.gz:CA_TF \
  --placenta-ocr placenta_ocr.bed.gz \
  --out "$AGG/finemap_enrichment.tsv"
```

Per annotation class: Fisher exact test (primary), logistic regression
adjusting for log10(distance to phenotype start) (sensitivity), and a
PIP-weighted enrichment; BH-FDR within each test family.

### 7.5 Diagnostic report

```bash
singularity exec --bind /rsrch5 --bind /rsrch9 \
  /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif \
  Rscript -e 'rmarkdown::render("reports/report_finemap.Rmd",
    params=list(agg_dir="'$RESULTS_DIR'/finemap/aggregated",
                enrichment_path="'$RESULTS_DIR'/finemap/aggregated/finemap_enrichment.tsv",
                fig_dir="reports/fig_finemap", table_dir="reports/tables_finemap"))'
```

Sections: overview, max-PIP-per-locus and PIP distributions, credible-set
sizes and CSs per locus, cross-ancestry rho, per-ancestry effect-weight
concordance, top-locus PIP tracks, enrichment forest plot, run diagnostics.

### 7.6 Fallback

Loci with `status=failed` or non-converged SuSHiE runs (see
`finemap_locus_summary.tsv`) are candidates for the documented fallback:
single-ancestry FINEMAP, or restricting to the strongest single signal.
Not implemented here — flag and handle case-by-case.

## Notes

- k grid: 0 5 10 15 20 25 30 for every group (same as the expression-only
  module, for comparability across modalities).
- HCP parameters: lambda1 = 0.5, lambda2 = lambda3 = 1, QC |r| > 0.9 pruning —
  identical to script 19.
- `ct_Maternal` is excluded before correlation pruning in every covariate
  build; no covariate cap.
- The design is ancestry-agnostic: when AFR/AMR/SAS inputs land, add the
  labels to the ANC loop in steps 2-4 (and ANCESTRIES in 28 and 32) with
  no code changes. SuSHiE then fine-maps jointly across all listed
  ancestries.
- Per-k staging lives under `$QTL_DIR/hcp_optimization_modalities/{ANC}/{MOD}/`
  (logs, per-k HCPs, covariates, mapping parquets) — keep it until the
  results are signed off; it is the audit trail for each k* choice.
