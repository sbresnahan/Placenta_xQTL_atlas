#!/bin/bash
# =============================================================================
# 25b_submit_hcp_k_jobs.sh — modality HCP optimization, one LSF job per k
# =============================================================================
# Submission driver only; run after 26_harmonize_modalities.py.
#
# Reuses completed k points from either:
#   * legacy serial staging: qtl_inputs/hcp_optimization_modalities/{ANC}/{MOD}/; or
#   * isolated per-k jobs:  qtl_inputs/hcp_optimization_modalities_per_k/{ANC}/{MOD}/kK/.
#
# Required env: CONFIG, SCRIPTS_DIR, OUTPUT_BASE, LOGS_DIR
# Optional: QTL_DIR, ANCESTRIES, MODALITIES, K_GRID, FDR, CHR1_MIN,
#           MAX_HCP_PHENOTYPES, EXCLUDE_COVARIATES, LAMBDA1, FINALIZE
#   LAMBDA1 — HCP prior strength passed to optimize_hcp_modalities.py
#   (default 0.5). When LAMBDA1 != 0.5 the per-k sandboxes live in a
#   lambda-suffixed tree (hcp_optimization_modalities_per_k_lam<LAMBDA1>/),
#   the lambda=0.5 legacy tree is never consulted, job names carry a
#   _lam<LAMBDA1> suffix, and the canonical finalizer is SKIPPED unless
#   FINALIZE=1 is set explicitly.
# =============================================================================

set -euo pipefail

: "${CONFIG:?ERROR: CONFIG required}"
: "${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR required}"
: "${OUTPUT_BASE:?ERROR: OUTPUT_BASE required}"
: "${LOGS_DIR:?ERROR: LOGS_DIR required}"

REAL_QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
MODALITIES="${MODALITIES:-isoform_expression alt_polyA alt_TSS intron_retention isoforms RNA_editing splicing stability}"
K_GRID="${K_GRID:-0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100}"
FDR="${FDR:-0.05}"
CHR1_MIN="${CHR1_MIN:-300}"
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
echo "lambda1: $LAMBDA1 (sandbox tree: hcp_optimization_modalities_per_k${SANDBOX_SUFFIX})"

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
  local anc="$1" mod="$2" k="$3"
  local stage="$REAL_QTL_DIR/hcp_optimization_modalities/$anc/$mod"
  [ -s "$stage/${anc}_${mod}_hcp_k${k}.tsv" ] || return 1
  valid_parquet "$stage/results_k${k}/${anc}_${mod}_cisqtl.parquet"
}

sandbox_root() {
  echo "$REAL_QTL_DIR/hcp_optimization_modalities_per_k${SANDBOX_SUFFIX}/$1/$2/k$3"
}

sandbox_done() {
  local anc="$1" mod="$2" k="$3" root
  root=$(sandbox_root "$anc" "$mod" "$k")
  [ -s "$root/qtl_inputs/${anc}_${mod}_hcp_factors_optimized.tsv" ] || return 1
  tsv_has_k "$root/work/${anc}_${mod}_optimal_hcp.tsv" "$k"
}

active_job_ids() {
  local name="$1"
  bjobs -J "$name" 2>/dev/null | awk 'NR>1 && $1 ~ /^[0-9]+$/ {print $1}'
}

prepare_sandbox() {
  local anc="$1" mod="$2" k="$3" root qdir f groups
  root=$(sandbox_root "$anc" "$mod" "$k")
  qdir="$root/qtl_inputs"
  mkdir -p "$qdir" "$root/work"

  for f in \
    "${anc}_metadata.tsv" \
    "${anc}_${mod}_harmonized.bed" \
    "${anc}_qtl.pgen" "${anc}_qtl.pvar" "${anc}_qtl.psam" \
    "${anc}_selected_pcs.txt" \
    "${anc}_deconvolution_harmonized.tsv"; do
    [ -e "$REAL_QTL_DIR/$f" ] || { echo "ERROR: missing $REAL_QTL_DIR/$f" >&2; exit 1; }
    ln -sfn "$REAL_QTL_DIR/$f" "$qdir/$f"
  done

  groups="$REAL_QTL_DIR/${anc}_${mod}.phenotype_groups.txt"
  if [ -e "$groups" ]; then
    ln -sfn "$groups" "$qdir/${anc}_${mod}.phenotype_groups.txt"
  fi
}

parse_job_id() {
  sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' | head -1
}

for ANC in $ANCESTRIES; do
  for MOD in $MODALITIES; do
    canonical_summary="$REAL_QTL_DIR/hcp_optimization_modalities/${ANC}_${MOD}_optimal_hcp.tsv"
    canonical_hcp="$REAL_QTL_DIR/${ANC}_${MOD}_hcp_factors_optimized.tsv"

    # Canonical outputs belong to the production lambda; sensitivity runs
    # (non-default LAMBDA1) never consult them.
    if [ -z "$SANDBOX_SUFFIX" ] && summary_has_grid "$canonical_summary" && [ -s "$canonical_hcp" ]; then
      echo "[$ANC $MOD] already finalized for the full k grid; nothing to submit."
      continue
    fi

    # The old runbook used this exact monolithic job name. Refuse to launch
    # per-k workers while one of those is still active, because it writes into
    # the legacy staging tree and can overwrite the canonical HCP on completion.
    mapfile -t OLD_ACTIVE < <(active_job_ids "hcpopt_${ANC}_${MOD}")
    if [ "${#OLD_ACTIVE[@]}" -gt 0 ]; then
      echo "ERROR: legacy monolithic job hcpopt_${ANC}_${MOD} is still active: ${OLD_ACTIVE[*]}" >&2
      echo "Cancel/finish it before launching isolated per-k jobs for this pair." >&2
      exit 1
    fi

    deps=()
    echo "[$ANC $MOD] scanning k grid: $K_GRID"

    for K in $K_GRID; do
      if sandbox_done "$ANC" "$MOD" "$K"; then
        echo "  SKIP k=$K (completed isolated per-k job)"
        continue
      fi
      if legacy_done "$ANC" "$MOD" "$K"; then
        echo "  SKIP k=$K (completed legacy serial grid point)"
        continue
      fi

      JOBNAME="hcp25b_${ANC}_${MOD}_k${K}${JOB_SUFFIX}"
      mapfile -t ACTIVE < <(active_job_ids "$JOBNAME")
      if [ "${#ACTIVE[@]}" -gt 0 ]; then
        echo "  ACTIVE k=$K (${JOBNAME}: ${ACTIVE[*]}); not resubmitting"
        deps+=("${ACTIVE[@]}")
        continue
      fi

      prepare_sandbox "$ANC" "$MOD" "$K"
      ROOT=$(sandbox_root "$ANC" "$MOD" "$K")
      JOB_QTL="$ROOT/qtl_inputs"
      JOB_WORK="$ROOT/work"

      Q=medium
      W=24:00
      if [ "$MOD" = "isoform_expression" ]; then
        Q=long
        W=48:00
      fi

      OUT=$(
        export ANCESTRIES="$ANC"
        export MODALITIES="$MOD"
        export K_GRID="$K"
        export SKIP_EXISTING=1
        export QTL_DIR="$JOB_QTL"
        export WORK_DIR="$JOB_WORK"
        export FDR CHR1_MIN
        export EXTRA_ARGS="--lambda1 ${LAMBDA1}"
        bsub -J "$JOBNAME" \
          -q "$Q" -n 4 -M 32 -R "rusage[mem=32]" -W "$W" \
          -o "$LOGS_DIR/${JOBNAME}.%J.out" \
          -e "$LOGS_DIR/${JOBNAME}.%J.err" \
          -env all < "$SCRIPTS_DIR/25b_optimize_hcp_modalities.sh"
      )
      echo "  $OUT"
      JID=$(printf '%s\n' "$OUT" | parse_job_id)
      [ -n "$JID" ] || { echo "ERROR: could not parse LSF job id for $JOBNAME" >&2; exit 1; }
      deps+=("$JID")
    done

    # Sensitivity runs (non-default LAMBDA1) skip the canonical finalizer
    # unless FINALIZE=1 is set explicitly (lambda promoted to production).
    if [ -n "$SANDBOX_SUFFIX" ] && [ "$FINALIZE" != "1" ]; then
      echo "  Sensitivity mode (lambda1=$LAMBDA1): per-k jobs submitted;"
      echo "  finalizer skipped for $ANC $MOD."
      continue
    fi

    FINAL_NAME="hcp25b_${ANC}_${MOD}_finalize${JOB_SUFFIX}"
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
      export HCP_FINALIZE_MODE=modality
      export ANCESTRY="$ANC"
      export MODALITY="$MOD"
      export K_GRID FDR CHR1_MIN
      export SANDBOX_SUFFIX
      "${BSUB[@]}" -env all < "$SCRIPTS_DIR/25_hcp_k_finalize.sh"
    )
    echo "  $OUT"
  done
done
