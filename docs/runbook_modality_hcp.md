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

## Notes

- k grid: 0 5 10 15 20 25 30 for every group (same as the expression-only
  module, for comparability across modalities).
- HCP parameters: lambda1 = 0.5, lambda2 = lambda3 = 1, QC |r| > 0.9 pruning —
  identical to script 19.
- `ct_Maternal` is excluded before correlation pruning in every covariate
  build; no covariate cap.
- The design is ancestry-agnostic: when AFR/AMR/SAS inputs land, add the
  labels to the ANC loop in steps 2-4 (and ANCESTRIES in 28) with no code
  changes.
- Per-k staging lives under `$QTL_DIR/hcp_optimization_modalities/{ANC}/{MOD}/`
  (logs, per-k HCPs, covariates, mapping parquets) — keep it until the
  results are signed off; it is the audit trail for each k* choice.
