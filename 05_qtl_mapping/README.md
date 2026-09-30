# 05_qtl_mapping

Ancestry-stratified cis-xQTL mapping under GTEx conventions: cross-cohort
pooling with QN + INT + ComBat (normalization first, ComBat last; devBrain
xQTL schema, Wen et al., Science 2024, 384:eadh0829), HCP latent-factor
estimation with per-modality k optimization, cohort-only genotype PCA,
genotype x phenotype intersection with a minor-allele-count floor,
PC-outlier exclusion, covariate assembly, tensorQTL mapping with Storey
q-values, SuSHiE cross-ancestry fine-mapping, and functional enrichment of
high-PIP variants.

These scripts share `config.yml` / `config_get.py` with
`../03_phenotyping/` (one scripts directory on seadragon).

## Prerequisites

- Stage-1 `{ANC}_pooled.pgen`, stage-3 modality BEDs + deconvolution
  proportions, stage-4 `rnaseq_to_array_id_map.csv`, cohort metadata
  (`rnaseq_id, array_id, ancestry, cohort, sex, GA, ppBMI`)
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
| 19 | `19_hcp_factors.sh`, `19a_picard_sharded.sh`, `19b_fix_missing_metrics.sh` | HCP pipeline: Picard QC (`picard_qc.py`) -> pool metrics (`fix_missing_metrics.py`) -> pool expression within ancestry (`pool_expression_within_ancestry.py`) -> TPM > 0.1 in > 25% filter -> QN + INT -> connectivity-outlier removal (bicor, z < -3; writes `{ANC}_expression_outliers.tsv`) -> ComBat last -> HCP, provisional k=15 (`combat_normalize_hcp.R`, `peer_factors.py`) |
| 20 | `20_combat_modalities.sh` | Pool + QN + INT + ComBat (last), 7 non-expression modalities (`pool_modalities_within_ancestry.py`, `combat_normalize_modalities.R`); devBrain filters (isoforms TPM > 0.1/>25%; others detected >= 40%); isoforms exclude expression outliers via `--exclude-samples`; splicing/IR via stage-17 pre-pooled BEDs. HCP factors from 19 are reused |
| 21 | `21_install_tensorqtl.sh` | One-time setup: tensorqtl conda env + R `qvalue` + Storey-bridge smoke test |
| 22 | `22_genotype_pca.sh` | Cohort-only genotype PCA: LD-prune (`--indep-pairwise 200 50 0.2`) -> `plink2 --pca 20 exact` -> `genotype_pca_format.py` -> `{ANC}_genotype_pcs.tsv` + scree (`PCA_scree.R`) |
| 23 | `23_prepare_intersection.py` | Genotype x phenotype intersection on the RNA-to-DNA map; pgen filtered to intersection samples with MAC >= 5 (`--mac 5`; `--mac 0` disables); sample columns renamed rnaseq_id -> array_id |
| 24 | `24_outlier_exclusion.py` | First-5-PC selection; 6-SD PC outliers removed from pgen/BED/HCP/deconvolution/metadata. Edits intersection files in place — always rerun 23 before 24 |
| 25 | `25_build_covariates.py` | Covariates: PC1-5 + HCP_1-k + sex + GA + cell types (dominant type as compositional reference); `--exclude-covariates ct_Maternal` applied before \|r\| > 0.9 pruning; near-zero-variance pre-filter; pruning priority sex/GA > PCs > cell types > HCPs (pruned HCPs reduce effective k below nominal); no covariate cap. Technical-replicate sample columns (duplicate array_id) are averaged per covariate; discordant-sex replicates are excluded. `--hcp-file`/`--hcp-k`/`--out-suffix` support the 25a/25b optimization modules; the canonical per-modality run writes `{ANC}_covariates_{MOD}.tsv` |
| 25a | `25a_optimize_hcp.sh`, `optimize_hcp_chr1.py` | Expression-only HCP-count optimization (devBrain section 4.2), run after 23/24 and before canonical 25: re-estimate HCP per k in {0,5,10,15,20,25,30}, map chr1 expression per k (1 Mb window, MAF >= 0.01), pick k\* maximizing eGenes at Storey q <= 0.05 (ties -> smaller k); installs k\* as `{ANC}_hcp_factors_harmonized.tsv` (provisional file backed up to `*.pre25a_backup.tsv`); writes `{ANC}_optimal_hcp.tsv` + `.png`. Every per-k model uses the same pruned fixed covariate set (`EXCLUDE_COVARIATES`, default `ct_Maternal`); HCPs lost to correlation pruning are reported (`n_hcp_used`/`n_hcp_dropped`) |
| 25b | `25b_optimize_hcp_modalities.sh`, `optimize_hcp_modalities.py`, `hcp_from_matrix.R` | Per-modality HCP-count optimization: per ancestry x modality group (9 modalities + combined), HCP-only estimation from that group's harmonized BED per k (deterministic 40k-phenotype subsample caps cost for the combined arm), chr1-subset mapping (genome-wide when < 300 chr1 phenotypes), k\* = argmax eGenes at Storey q <= 0.05 (ties -> smaller k); installs `{ANC}_{MOD}_hcp_factors_optimized.tsv`, writes `{ANC}_{MOD}_optimal_hcp.tsv`/`.png` to `qtl_inputs/hcp_optimization_modalities/`. One LSF job per ancestry x modality; full procedure below |
| 26 | `26_harmonize_modalities.py` | Harmonize the 7 non-expression modality BEDs to the final array_id sample set (handles stage-17 namespaced IDs for splicing/IR) |
| 30 | `30_combine_modalities.py` | Combined cross-modality BED (`{modality}__{id}` namespacing; cross-modality gene groups; modality sidecar TSV) |
| 27 | `27_run_tensorqtl.sh` + `27_run_tensorqtl.py` | tensorQTL `cis.map_cis` per ancestry x modality (grouped, `group_s`, when `phenotype_groups.txt` exists); `--independent` stepwise conditional mode; Storey q-values via the `compute_qvalues.R` file bridge (`QVALUE_RSCRIPT`), `--qvalue-method bh` as fallback. Covariates default to `{ANC}_covariates_{MOD}.tsv` (25b), falling back to `{ANC}_covariates.tsv` with a warning; `COVARIATES_FILE` accepts `{ANC}` and `{MOD}` placeholders |
| 28 | `28_submit_modalities.sh` | Submission driver: one LSF job per ancestry x modality (`TEST=1` pilot; threads `QVALUE_METHOD`/`MAF_THRESHOLD`) |
| 29 | `29_make_top_tables.py` | Rebuild sorted `*_cisqtl_top.tsv` from parquets (no tensorQTL rerun) |
| 31 | `31_sushie_finemap.py` | Cross-ancestry fine-mapping: `prepare-loci` builds per-modality locus lists (union of q <= 0.05 grouped-layer lead phenotypes across ancestries; tested windows from the mapping parquets; L = min(10, max(5, n_independent + 2)) from the stepwise layer); `run` fine-maps one shard of loci jointly across ancestries with SuSHiE (individual-level mode, in-sample LD from the intersected pgens; purity 0.5; phenotypes missing from an ancestry's BED drop that ancestry for the locus; per-locus `.ancestries`/`.done` markers + per-shard diagnostics) |
| 32 | `32_submit_sushie.sh` + `32a_run_sushie_shard.sh` | LSF driver: prepare loci per modality -> shard (default 50 loci/shard) -> one array job per shard in the `sushie` conda env (`TEST=1` pilot; skip-if-done/running guards; `FORCE_LOCI`/`FORCE_RUN`; `ANCESTRIES`/`MODALITIES`/`SHARD_SIZE`/`QUEUE`/`WALLTIME`/`THREADS` env overrides) |
| 33 | `33_aggregate_finemap.py` | Aggregate per-locus SuSHiE outputs -> `finemap/aggregated/`: `finemap_pips.tsv.gz` (per-variant PIPs, CS membership, per-ancestry effect weights, locus diagnostics), `finemap_credible_sets.tsv.gz` (per-CS summaries incl. cross-ancestry rho), `finemap_locus_summary.tsv` (convergence, ELBO, n CS, max PIP per locus) |
| 34 | `34_pip_annotation_enrichment.py` | High-PIP (>= 0.9) variant annotation enrichment: `--make-fastvep-input` -> fastVEP consequences (splice/LoF/missense/synonymous/UTR/intron/regulatory/flanking/intergenic); intersects ENCODE SCREEN cCRE classes + placenta OCR BED; Fisher exact (primary), log10(distance)-adjusted logistic (sensitivity), PIP-weighted enrichment; BH-FDR per test family |

Diagnostics/utilities: `hcp_diagnostic.R`, `hcp_diagnostic2.R`,
`check_hcp_chunks.sh`, `test_hcp.py`, `test_combat_modalities.py`,
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

## Step 2: HCP optimization grid (25b; 20 LSF jobs)

One job per ancestry x modality group. isoform_expression and combined go
to the long queue; the rest fit in medium. Expected wall time: ~14 h for
the longest job (combined arm).

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

Monitor with `bjobs -w | grep hcpopt`. If a job fails, fix the cause and
resubmit that pair with `SKIP_EXISTING=1` appended to the `-env` list —
completed k grid points are reused.

Completion checks (both counts must be 20):

```bash
ls "$QTL_DIR"/hcp_optimization_modalities/*_optimal_hcp.tsv | wc -l
ls "$QTL_DIR"/*_hcp_factors_optimized.tsv | wc -l
```

Chosen k\* per ancestry x modality:

```bash
for f in "$QTL_DIR"/hcp_optimization_modalities/*_optimal_hcp.tsv; do
  awk -v f="$(basename "$f" _optimal_hcp.tsv)" '$NF=="True" {print f, "k*="$3, "eGenes="$4, "scope="$6}' "$f"
done
```

HCP parameters match script 19: lambda1 = 0.5, lambda2 = lambda3 = 1,
QC |r| > 0.9 pruning; k grid {0,5,10,15,20,25,30} for every group. Per-k
staging lives under `$QTL_DIR/hcp_optimization_modalities/{ANC}/{MOD}/`
(logs, per-k HCPs, covariates, mapping parquets) — keep it until results
are signed off; it is the audit trail for each k\* choice.

## Step 3: canonical per-modality covariates (25)

Builds `{ANC}_covariates_{MOD}.tsv` for all 10 groups per ancestry from
the installed k\* HCP sets, then copies the expression set to the
canonical `{ANC}_covariates.tsv` (kept for the report archive and other
downstream consumers). Light enough for a login node; wrap in bsub if
preferred (~1-2 min per modality).

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
SRR13696964, SRR13696996, SRR13697029 — the same individuals processed in
two cohort batches). A `WARN: discordant sex across technical replicates`
line indicates a metadata inconsistency — investigate before proceeding
if it appears.

## Step 4: xQTL mapping (three submission passes)

`28_submit_modalities.sh` skips combos with existing results or running
jobs, so all three passes can be rerun freely. Each mapping job picks up
`{ANC}_covariates_{MOD}.tsv` automatically. Expected wall time: ~1 day.

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

Completion check (expect 36 primary + 20 independent parquets: per
ancestry, 9 grouped + 8 ungrouped + 1 combined primary, 9
grouped-independent + 1 combined-independent):

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

## Step 5: top tables + report-input archive

```bash
python3 "$SCRIPTS_DIR/29_make_top_tables.py" --results-dir "$RESULTS_DIR"

cd "$REPO_DIR"
CONFIG="$CONFIG" bash ../reports/make_report_archive.sh
```

The archive (`placenta_xqtl_report_inputs_<date>.tar.gz`, written to the
repo root) includes `data/qc/hcp_optimization_modalities/` with the
per-modality k\* tables and curves alongside `data/qc/hcp_optimization/`.
Review the MANIFEST.txt it prints: no core inputs may be missing.

## Step 6: cross-ancestry fine-mapping (SuSHiE)

SuSHiE joint multi-ancestry fine-mapping of every grouped-layer lead
phenotype at FDR <= 5%, using in-sample LD from the intersected pgens
(individual-level mode — LD and genotypes are guaranteed consistent).
Each locus is fine-mapped jointly across all ancestries in `ANCESTRIES`,
regardless of which ancestry reached significance. Outputs: 95% credible
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
union of q <= 0.05 lead phenotypes across ancestries, tested windows from
the mapping parquets, L = min(10, max(5, n_independent + 2))), shards them
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
- **`optimize_hcp_chr1.py` (25a) overwrites
  `{ANC}_hcp_factors_harmonized.tsv`** with the k\* solution (the
  provisional k=15 file is backed up to `*.pre25a_backup.tsv`). Always
  rerun 25 after 25a. The per-modality module (25b) instead installs
  `{ANC}_{MOD}_hcp_factors_optimized.tsv` per ancestry x modality group;
  rerun the canonical 25 per modality (with
  `--hcp-file {ANC}_{MOD}_hcp_factors_optimized.tsv --out-suffix _{MOD}`)
  so every `{ANC}_covariates_{MOD}.tsv` matches its installed k\*.
- **`picard` must resolve after 19/19a's env stack** — `picard_qc.py`
  calls the literal `picard` executable; if the `picard-2.27.4`
  activation does not put it on PATH (wrong env name/path, or a
  PATH-clobbering venv activation), every QC call WARNs `[Errno 2] No
  such file or directory: 'picard'` and writes garbage metrics. Both
  wrappers fail fast with diagnostics — fix the activation. Stage 19
  skips existing `{COHORT}_qc_metrics.tsv`, so delete tainted QC outputs
  before resubmitting.
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
