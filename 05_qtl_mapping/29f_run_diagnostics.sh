#!/bin/bash
# =============================================================================
# 29f_run_diagnostics.sh — LSF driver for the 29b-29e validation-gate
#                          diagnostics (within-cohort-INT schema)
# =============================================================================
# Chains, per ancestry:
#   29b expression_diagnostics.py  (sample PCs, per-cohort residual
#                                   variance, HCP k15-vs-k45 exact pairs)
#   29c choi_comparison.py         (Choi 2024 comparison + retention/catalog;
#                                   only for CHOI_ANCESTRY, default EAS)
#   29d cohort_heterogeneity.py    (per-cohort scans, Cochran Q / I2, GxC)
#   29e plot_diagnostics.R         (figures + validation_summary.tsv)
#
# Run AFTER 28/29 (tensorQTL results + top tables) and BEFORE 32 (SuSHiE).
# Review $DIAG_DIR/<ANC>/validation_summary.tsv before fine-mapping.
#
# Usage (submitter — one LSF job per ancestry):
#   CONFIG=config.yml bash 29f_run_diagnostics.sh
#   TEST=1 CONFIG=config.yml bash 29f_run_diagnostics.sh   # first ancestry only
#
# Required env vars:
#   CONFIG        — path to config.yml (provides OUTPUT_BASE)
# Optional env vars:
#   SCRIPTS_DIR   — default: this script's own directory
#   ANCESTRIES    — default "EAS EUR"
#   QTL_DIR       — default ${OUTPUT_BASE}/qtl_inputs
#   RESULTS_DIR   — default ${OUTPUT_BASE}/qtl_results
#   DIAG_DIR      — default ${RESULTS_DIR}/diagnostics
#   HCP_OPT_DIR   — 25a staging dir (default ${QTL_DIR}/hcp_optimization)
#   K_LOW/K_HIGH  — HCP comparison grid points (default 15 / 45)
#   COHORT_A/COHORT_B — comparison cohorts (default: two largest per ancestry)
#   CHOI_ANCESTRY — ancestry with a Choi external reference (default EAS;
#                   empty string disables 29c everywhere)
#   CHOI_SUMSTATS / CHOI_SIGNIFICANT — Choi full summary stats (gzip TSV) and
#                   significant-eGene list (required when CHOI_ANCESTRY set)
#   BASELINE_DIR  — archived pre-fix diagnostics dir; 29e reads
#                   <BASELINE_DIR>/<ANC>/diagnostic_summary_stats.tsv and
#                   <BASELINE_DIR>/baseline_metrics.tsv when present
#   QUEUE         — default medium; WALLTIME — default 12:00
#   N_THREADS     — default 8
# =============================================================================

set -eo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG env var required (path to config.yml)}"
SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
CONFIG_GET="${SCRIPTS_DIR}/config_get.py"
[ -f "$CONFIG_GET" ] || CONFIG_GET="${SCRIPTS_DIR}/../03_phenotyping/config_get.py"
if [ ! -f "$CONFIG_GET" ]; then
    echo "ERROR: config_get.py not found (set SCRIPTS_DIR)" >&2
    exit 1
fi
eval "$(python3 "$CONFIG_GET" "$CONFIG")"

ANCESTRIES="${ANCESTRIES:-EAS EUR}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
DIAG_DIR="${DIAG_DIR:-${RESULTS_DIR}/diagnostics}"
HCP_OPT_DIR="${HCP_OPT_DIR:-${QTL_DIR}/hcp_optimization}"
K_LOW="${K_LOW:-15}"
K_HIGH="${K_HIGH:-45}"
CHOI_ANCESTRY="${CHOI_ANCESTRY-EAS}"   # note: '-' not ':-' so '' disables
QUEUE="${QUEUE:-medium}"
WALLTIME="${WALLTIME:-12:00}"
N_THREADS="${N_THREADS:-8}"
LOG_DIR="${LOG_DIR:-$(dirname "$OUTPUT_BASE")/PANTRY/logs}"
mkdir -p "$LOG_DIR" "$DIAG_DIR"

# ---------------------------------------------------------------------------
# Submitter mode: one job per ancestry
# ---------------------------------------------------------------------------
if [ -z "${DIAG_ANCESTRY:-}" ]; then
    n_sub=0
    for ANC in $ANCESTRIES; do
        job="diag29_${ANC}"
        if bjobs -noheader -J "$job" 2>/dev/null \
            | awk '$3 ~ /^(PEND|RUN|PSUSP|USUSP|SSUSP|WAIT)$/ {found=1} END {exit !found}'; then
            echo "  skip ${ANC}: job '${job}' already running/pending"
            continue
        fi
        if [ -f "${DIAG_DIR}/${ANC}/validation_summary.tsv" ] && \
           [ -z "${FORCE:-}" ]; then
            echo "  skip ${ANC}: validation_summary.tsv exists (FORCE=1 to redo)"
            continue
        fi
        # Export the ancestry-specific worker flag and all resolved paths/options so
        # `bsub -env all` passes clean values to the worker.  Resources remain
        # explicitly defined on the bsub command as required by the cluster.
        export DIAG_ANCESTRY="$ANC"
        export CONFIG SCRIPTS_DIR QTL_DIR RESULTS_DIR DIAG_DIR HCP_OPT_DIR
        export K_LOW K_HIGH CHOI_ANCESTRY QUEUE WALLTIME N_THREADS LOG_DIR
        export CHOI_SUMSTATS CHOI_SIGNIFICANT BASELINE_DIR COHORT_A COHORT_B

        bsub -J "$job" -q "$QUEUE" -n "$N_THREADS" -W "$WALLTIME" \
            -M 32 -R "rusage[mem=32]" \
            -o "${LOG_DIR}/diag29_${ANC}.%J.out" \
            -e "${LOG_DIR}/diag29_${ANC}.%J.err" \
            -env all \
            < "${BASH_SOURCE[0]}"

        unset DIAG_ANCESTRY
        n_sub=$((n_sub + 1))
        [ "${TEST:-0}" = "1" ] && { echo "TEST=1: submitted one job"; break; }
    done
    echo "Submitted $n_sub diagnostic job(s)."
    exit 0
fi

# ---------------------------------------------------------------------------
# Worker mode: run 29b -> 29c -> 29d -> 29e for $DIAG_ANCESTRY
# ---------------------------------------------------------------------------
ANC="$DIAG_ANCESTRY"
OUT="${DIAG_DIR}/${ANC}"
mkdir -p "$OUT"

eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
SING_R="singularity exec --bind /rsrch5 --bind /rsrch9 /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif Rscript"

echo "==================================================================="
echo "29b-29e diagnostics: ${ANC}"
echo "  QTL_DIR:     $QTL_DIR"
echo "  RESULTS_DIR: $RESULTS_DIR"
echo "  OUT:         $OUT"
echo "==================================================================="

# ---- 29b: expression diagnostics (tensorqtl env for the k-comparison) ----
conda activate tensorqtl
python3 "${SCRIPTS_DIR}/29b_expression_diagnostics.py" \
    --qtl-dir "$QTL_DIR" --ancestry "$ANC" --output-dir "$OUT" \
    --hcp-opt-dir "$HCP_OPT_DIR" --k-low "$K_LOW" --k-high "$K_HIGH"

# ---- 29c: Choi comparison (reference ancestry only) ----
CHOI_PREFIX="${OUT}/Choi_vs_pooled_${ANC}"
EXTRA_PAIRS_ARG=""
if [ "$ANC" = "$CHOI_ANCESTRY" ]; then
    [ -n "${CHOI_SUMSTATS:-}" ] && [ -n "${CHOI_SIGNIFICANT:-}" ] || {
        echo "ERROR: CHOI_ANCESTRY=${ANC} but CHOI_SUMSTATS/CHOI_SIGNIFICANT unset" >&2
        exit 1
    }
    python3 "${SCRIPTS_DIR}/29c_choi_comparison.py" --ancestry "$ANC" \
        --choi "$CHOI_SUMSTATS" --choi-significant "$CHOI_SIGNIFICANT" \
        --top-table "${RESULTS_DIR}/${ANC}_expression_cisqtl_top.tsv" \
        --out-prefix "$CHOI_PREFIX"
    EXTRA_PAIRS_ARG="--extra-pairs ${CHOI_PREFIX}.choi_defined_pairs.tsv"
else
    echo "  29c skipped: no external reference for ${ANC} (CHOI_ANCESTRY='${CHOI_ANCESTRY}')"
fi

# ---- 29d: cohort heterogeneity ----
COHORT_ARGS=""
[ -n "${COHORT_A:-}" ] && [ -n "${COHORT_B:-}" ] && \
    COHORT_ARGS="--cohorts ${COHORT_A},${COHORT_B}"
python3 "${SCRIPTS_DIR}/29d_cohort_heterogeneity.py" \
    --qtl-dir "$QTL_DIR" --results-dir "$RESULTS_DIR" --ancestry "$ANC" \
    --output-dir "$OUT" $COHORT_ARGS $EXTRA_PAIRS_ARG

# ---- 29e: figures + validation summary (R) ----
R_ARGS=(--ancestry "$ANC" --outdir "$OUT"
        --pcs "${OUT}/${ANC}_expression_sample_PCs.tsv"
        --variance "${OUT}/${ANC}_gene_residual_variance_by_cohort.tsv.gz"
        --hcp "${OUT}/${ANC}_expression_hcp${K_LOW}_vs_hcp${K_HIGH}_exact_pairs.tsv.gz"
        --heterogeneity "${OUT}/${ANC}_expression_cohort_heterogeneity.tsv.gz"
        --top-table "${RESULTS_DIR}/${ANC}_expression_cisqtl_top.tsv")
[ -n "$COHORT_ARGS" ] && R_ARGS+=(--cohort-a "$COHORT_A" --cohort-b "$COHORT_B")
if [ "$ANC" = "$CHOI_ANCESTRY" ]; then
    R_ARGS+=(--gene-top "${CHOI_PREFIX}.gene_top_comparison.tsv.gz"
             --exact-pairs "${CHOI_PREFIX}.exact_${ANC}_lead_comparison.tsv.gz"
             --retention "${CHOI_PREFIX}.Choi_significant_retention.tsv.gz"
             --catalog "${CHOI_PREFIX}.common_tested_catalog.tsv.gz"
             --extra-pairs "${OUT}/${ANC}_expression_extra_pairs_withincohort.tsv.gz")
fi
if [ -n "${BASELINE_DIR:-}" ]; then
    [ -f "${BASELINE_DIR}/${ANC}/diagnostic_summary_stats.tsv" ] && \
        R_ARGS+=(--baseline-summary "${BASELINE_DIR}/${ANC}/diagnostic_summary_stats.tsv")
    [ -f "${BASELINE_DIR}/baseline_metrics.tsv" ] && \
        R_ARGS+=(--baseline-metrics "${BASELINE_DIR}/baseline_metrics.tsv")
fi
$SING_R "${SCRIPTS_DIR}/29e_plot_diagnostics.R" "${R_ARGS[@]}"

echo "==================================================================="
echo "29b-29e done for ${ANC}. Review ${OUT}/validation_summary.tsv"
echo "before submitting fine-mapping (32_submit_sushie.sh)."
echo "==================================================================="
