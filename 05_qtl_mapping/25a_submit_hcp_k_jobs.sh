#!/bin/bash
# =============================================================================
# 25a_submit_hcp_k_jobs.sh — expression HCP optimization, one LSF job per k
# =============================================================================
# Submission driver only; run from a login shell after Steps 23/24.
#
# Reuses completed k points from either:
#   * the legacy serial staging tree: qtl_inputs/hcp_optimization/{ANC}/; or
#   * prior isolated per-k jobs:      qtl_inputs/hcp_optimization_per_k/{ANC}/kK/.
#
# A legacy k is considered complete only if its HCP file is non-empty and its
# tensorQTL parquet has intact PAR1 header/footer magic. New isolated jobs are
# complete only after the current 25a wrapper has written both its one-row
# optimal_hcp TSV and sandbox-installed HCP file.
#
# Required env: CONFIG, SCRIPTS_DIR, OUTPUT_BASE, ANCESTRY_MAP, LOGS_DIR
# Optional: QTL_DIR, ANCESTRIES, K_GRID, FDR, EXCLUDE_COVARIATES, LAMBDA1, FINALIZE
#   LAMBDA1 — HCP prior strength passed to optimize_hcp_chr1.py (default 0.5).
#   When LAMBDA1 != 0.5 the per-k sandboxes live in a lambda-suffixed tree
#   (hcp_optimization_per_k_lam<LAMBDA1>/), the lambda=0.5 legacy tree is
#   never consulted, job names carry a _lam<LAMBDA1> suffix, and the
#   canonical finalizer is SKIPPED unless FINALIZE=1 is set explicitly
#   (sensitivity scans read per-k parquets directly; set FINALIZE=1 when a
#   non-default lambda becomes the production choice).
# =============================================================================

set -eo pipefail

: "${CONFIG:?ERROR: CONFIG required}"
: "${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR required}"
: "${OUTPUT_BASE:?ERROR: OUTPUT_BASE required}"
: "${ANCESTRY_MAP:?ERROR: ANCESTRY_MAP required}"
: "${LOGS_DIR:?ERROR: LOGS_DIR required}"

REAL_QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
K_GRID="${K_GRID:-0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100}"
FDR="${FDR:-0.05}"
LAMBDA1="${LAMBDA1:-0.5}"
FINALIZE="${FINALIZE:-0}"
mkdir -p "$LOGS_DIR"

# Lambda-suffixed sandbox tree + job-name suffix for non-default lambda1.
SANDBOX_SUFFIX=""
JOB_SUFFIX=""
if [ "$LAMBDA1" != "0.5" ]; then
  SANDBOX_SUFFIX="_lam${LAMBDA1}"
  JOB_SUFFIX="_lam${LAMBDA1}"
fi
echo "lambda1: $LAMBDA1 (sandbox tree: hcp_optimization_per_k${SANDBOX_SUFFIX})"

valid_parquet() {
  local f="$1" head tail
  [ -s "$f" ] || return 1
  head=$(head -c 4 "$f" 2>/dev/null || true)
  tail=$(tail -c 4 "$f" 2>/dev/null || true)
  [ "$head" = "PAR1" ] && [ "$tail" = "PAR1" ]
}

tsv_has_k() {
  local f="$1" wanted="$2"
  [ -s "$f" ] || return 1
  awk -F '\t' -v wanted="$wanted" '
    NR==1 { for (i=1; i<=NF; i++) if ($i=="k") kc=i; next }
    kc && $kc==wanted { found=1 }
    END { exit(found ? 0 : 1) }
  ' "$f"
}

summary_has_grid() {
  local f="$1" k
  [ -s "$f" ] || return 1
  for k in $K_GRID; do
    tsv_has_k "$f" "$k" || return 1
  done
}

legacy_done() {
  # Legacy serial tree holds lambda1=0.5 results only — never reuse it for
  # non-default-lambda runs.
  [ -z "$SANDBOX_SUFFIX" ] || return 1
  local anc="$1" k="$2"
  local stage="$REAL_QTL_DIR/hcp_optimization/$anc"
  [ -s "$stage/${anc}_hcp_k${k}_harmonized.tsv" ] || return 1
  valid_parquet "$stage/results_k${k}/${anc}_expression_cisqtl.parquet"
}

sandbox_root() {
  echo "$REAL_QTL_DIR/hcp_optimization_per_k${SANDBOX_SUFFIX}/$1/k$2"
}

sandbox_done() {
  local anc="$1" k="$2" root
  root=$(sandbox_root "$anc" "$k")
  [ -s "$root/qtl_inputs/${anc}_hcp_factors_harmonized.tsv" ] || return 1
  tsv_has_k "$root/work/${anc}_optimal_hcp.tsv" "$k"
}

active_job_ids() {
  local name="$1"
  bjobs -J "$name" 2>/dev/null | awk 'NR>1 && $1 ~ /^[0-9]+$/ {print $1}'
}

prepare_sandbox() {
  local anc="$1" k="$2" root qdir f
  root=$(sandbox_root "$anc" "$k")
  qdir="$root/qtl_inputs"
  mkdir -p "$qdir" "$root/work"

  for f in \
    "${anc}_metadata.tsv" \
    "${anc}_expression_harmonized.bed" \
    "${anc}_qtl.pgen" "${anc}_qtl.pvar" "${anc}_qtl.psam" \
    "${anc}_selected_pcs.txt" \
    "${anc}_deconvolution_harmonized.tsv"; do
    [ -e "$REAL_QTL_DIR/$f" ] || { echo "ERROR: missing $REAL_QTL_DIR/$f" >&2; exit 1; }
    ln -sfn "$REAL_QTL_DIR/$f" "$qdir/$f"
  done
}

parse_job_id() {
  sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' | head -1
}

for ANC in $ANCESTRIES; do
  canonical_summary="$REAL_QTL_DIR/hcp_optimization/${ANC}_optimal_hcp.tsv"
  canonical_hcp="$REAL_QTL_DIR/${ANC}_hcp_factors_harmonized.tsv"

  # Canonical outputs belong to the production lambda; sensitivity runs
  # (non-default LAMBDA1) never consult them.
  if [ -z "$SANDBOX_SUFFIX" ] && summary_has_grid "$canonical_summary" && [ -s "$canonical_hcp" ]; then
    echo "[$ANC] already finalized for the full k grid; nothing to submit."
    continue
  fi

  deps=()
  echo "[$ANC] scanning k grid: $K_GRID"

  for K in $K_GRID; do
    if sandbox_done "$ANC" "$K"; then
      echo "  SKIP k=$K (completed isolated per-k job)"
      continue
    fi
    if legacy_done "$ANC" "$K"; then
      echo "  SKIP k=$K (completed legacy serial grid point)"
      continue
    fi

    JOBNAME="hcp25a_${ANC}_k${K}${JOB_SUFFIX}"
    mapfile -t ACTIVE < <(active_job_ids "$JOBNAME")
    if [ "${#ACTIVE[@]}" -gt 0 ]; then
      echo "  ACTIVE k=$K (${JOBNAME}: ${ACTIVE[*]}); not resubmitting"
      deps+=("${ACTIVE[@]}")
      continue
    fi

    prepare_sandbox "$ANC" "$K"
    ROOT=$(sandbox_root "$ANC" "$K")
    JOB_QTL="$ROOT/qtl_inputs"
    JOB_WORK="$ROOT/work"

    OUT=$(
      export ANCESTRIES="$ANC"
      export K_GRID="$K"
      export SKIP_EXISTING=1
      export QTL_DIR="$JOB_QTL"
      export WORK_DIR="$JOB_WORK"
      export FDR
      export EXTRA_ARGS="--lambda1 ${LAMBDA1}"
      bsub -J "$JOBNAME" \
        -q long -n 4 -M 32 -R "rusage[mem=32]" -W 48:00 \
        -o "$LOGS_DIR/${JOBNAME}.%J.out" \
        -e "$LOGS_DIR/${JOBNAME}.%J.err" \
        -env all < "$SCRIPTS_DIR/25a_optimize_hcp.sh"
    )
    echo "  $OUT"
    JID=$(printf '%s\n' "$OUT" | parse_job_id)
    [ -n "$JID" ] || { echo "ERROR: could not parse LSF job id for $JOBNAME" >&2; exit 1; }
    deps+=("$JID")
  done

  # Sensitivity runs (non-default LAMBDA1) skip the canonical finalizer
  # unless FINALIZE=1 is set explicitly (i.e. this lambda is being promoted
  # to production).
  if [ -n "$SANDBOX_SUFFIX" ] && [ "$FINALIZE" != "1" ]; then
    echo "  Sensitivity mode (lambda1=$LAMBDA1): per-k jobs submitted;"
    echo "  finalizer skipped. Compare per-k parquets under"
    echo "    $REAL_QTL_DIR/hcp_optimization_per_k${SANDBOX_SUFFIX}/${ANC}/"
    echo "  against the production tree, then re-run with FINALIZE=1 if this"
    echo "  lambda is adopted."
    continue
  fi

  FINAL_NAME="hcp25a_${ANC}_finalize${JOB_SUFFIX}"
  mapfile -t ACTIVE_FINAL < <(active_job_ids "$FINAL_NAME")
  if [ "${#ACTIVE_FINAL[@]}" -gt 0 ]; then
    echo "  Finalizer already active (${FINAL_NAME}: ${ACTIVE_FINAL[*]}); not resubmitting"
    continue
  fi

  DEP_EXPR=""
  for jid in "${deps[@]}"; do
    [ -z "$DEP_EXPR" ] || DEP_EXPR+="&&"
    DEP_EXPR+="ended(${jid})"
  done

  BSUB=(bsub -J "$FINAL_NAME" -q medium -n 2 -M 8 -R "rusage[mem=8]" -W 4:00
        -o "$LOGS_DIR/${FINAL_NAME}.%J.out" -e "$LOGS_DIR/${FINAL_NAME}.%J.err")
  [ -z "$DEP_EXPR" ] || BSUB+=(-w "$DEP_EXPR")

  OUT=$(
    export REAL_QTL_DIR="$REAL_QTL_DIR"
    export HCP_FINALIZE_MODE=expression
    export ANCESTRY="$ANC"
    export K_GRID FDR
    export SANDBOX_SUFFIX
    "${BSUB[@]}" -env all < "$SCRIPTS_DIR/25_hcp_k_finalize.sh"
  )
  echo "  $OUT"
done
