#!/bin/bash
# =============================================================================
# 47_submit_isotwas.sh — Train isoTWAS/TWAS weights (one array per weight set)
# =============================================================================
# Objective 1.6: multivariate elastic-net isoform weights (isoTWAS) plus
# per-gene expression weights, trained per ancestry (EAS, EUR) and pooled
# (ancestry-stacked with ancestry-indicator residualization).
#
# The pooled weight set stacks per-ancestry dosage blocks in R (no merged
# pgen required).
# Step 1 counts genes per weight set. Step 2 submits one LSF job array per
# weight set; each element runs
# 47a_run_isotwas_shard.sh -> 46_isotwas_train.R on a round-robin shard of
# the gene universe (expression BED phenotype_ids; intersection across
# ancestries for the pooled set).
#
# Prerequisites:
#   - 36_install_coloc_env.sh completed (glmnet in the R library)
#   - module 05 inputs: {QTL_DIR}/{ANC}_qtl.pgen, {ANC}_expression.bed.gz,
#     {ANC}_isoform_expression.bed.gz, {ANC}_covariates_{MOD}.tsv
#
# Usage:
#   TEST=1 bash 47_submit_isotwas.sh   # ONE pilot shard (EAS, 25 genes)
#   bash 47_submit_isotwas.sh          # submit all weight sets
#
# Optional overrides (environment):
#   WEIGHT_SETS  — default "EAS EUR pooled"
#   SHARD_SIZE   — genes per array element, default 25 (measure with TEST=1)
#   QUEUE        — default medium
#   WALLTIME     — default 04:00
#   THREADS      — default 2
#   R2_MIN       — retention gate on CV R^2, default 0.01
#   FORCE_RUN=1  — retrain genes with existing .wgt.RDat (worker deletes them)
#   SKIP_COLLAPSE_CHECK=1 — skip the collapsed-replicate input audit
#
# TEST=1 submits only array index 1 while keeping the true N_SHARDS in the
# worker environment, so the pilot trains ~SHARD_SIZE genes (NOT the whole
# gene universe).
# =============================================================================
set -euo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"
REPO_ROOT="$(cd "${SCRIPTS_DIR}/.." && pwd)"
OUTPUT_BASE="${OUTPUT_BASE:-$(dirname "$CONFIG")}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
LOG_DIR="${LOG_DIR:-${OUTPUT_BASE}/logs}"
ISOTWAS_DIR="${ISOTWAS_DIR:-${RESULTS_DIR}/isotwas}"
WEIGHT_SETS="${WEIGHT_SETS:-EAS EUR pooled}"
SHARD_SIZE="${SHARD_SIZE:-25}"
QUEUE="${QUEUE:-medium}"
WALLTIME="${WALLTIME:-04:00}"
THREADS="${THREADS:-2}"
R2_MIN="${R2_MIN:-0.01}"
TEST="${TEST:-0}"
FORCE_RUN="${FORCE_RUN:-0}"
RSCRIPT="${RSCRIPT:-${REPO_ROOT}/bin/Rscript_sif}"
PLINK2="${PLINK2:-plink2}"

mkdir -p "$ISOTWAS_DIR" "$LOG_DIR"

echo "=== 47_submit_isotwas.sh ==="
echo "  WEIGHT_SETS: $WEIGHT_SETS"
echo "  ISOTWAS_DIR: $ISOTWAS_DIR"
echo "  SHARD_SIZE:  $SHARD_SIZE   QUEUE: $QUEUE   WALLTIME: $WALLTIME"

# --- Collapsed-replicate input audit ------------------------------------------
# Weight training consumes sample-level qtl_inputs (BEDs, covariates, pgen).
# Verify those sample sets match the module-05 collapsed-replicate contract
# (the same ancestry_map_collapsed.tsv authority module 07 enforces) before
# training weights on them. Stdlib-only script: no conda env required.
if [ "${SKIP_COLLAPSE_CHECK:-0}" != "1" ]; then
    python3 "${SCRIPTS_DIR}/check_collapsed_inputs.py" \
        --qtl-dir "$QTL_DIR" \
        --output-base "$OUTPUT_BASE" \
        --ancestries EAS EUR
fi

# --- Step 1: gene counts per weight set ---------------------------------------
gene_count() {
    local WS="$1"
    if [ "$WS" = "pooled" ]; then
        python3 - "$QTL_DIR" $WEIGHT_SETS <<'PYEOF'
import sys, subprocess
qtl_dir = sys.argv[1]
ancs = [a for a in sys.argv[2:] if a != "pooled"]
sets = []
for anc in ancs:
    out = subprocess.run(
        f"zcat {qtl_dir}/{anc}_expression.bed.gz | cut -f4 | tail -n +2",
        shell=True, capture_output=True, text=True)
    sets.append(set(out.stdout.split()))
print(len(set.intersection(*sets)))
PYEOF
    else
        echo $(( $(zcat "${QTL_DIR}/${WS}_expression.bed.gz" | cut -f4 | tail -n +2 | wc -l) ))
    fi
}

# --- Step 2: submit one array per weight set ----------------------------------
N=0
for WS in $WEIGHT_SETS; do
    N_GENES=$(gene_count "$WS")
    if [ "$N_GENES" -le 0 ]; then
        echo "  SKIP ${WS}: 0 genes"
        continue
    fi
    N_SHARDS=$(( (N_GENES + SHARD_SIZE - 1) / SHARD_SIZE ))
    # TEST=1: submit ONLY array index 1 but keep the true N_SHARDS in the
    # worker environment (round-robin assignment; N_SHARDS=1 would give the
    # pilot shard EVERY gene).
    ARRAY_SPEC="1-${N_SHARDS}"
    if [ "$TEST" = "1" ]; then ARRAY_SPEC="1"; fi

    DONE_SHARDS=0
    for DIAG in "${ISOTWAS_DIR}/diagnostics/isotwas_${WS}".shard-*.diagnostics.tsv; do
        [ -f "$DIAG" ] && DONE_SHARDS=$((DONE_SHARDS + 1))
    done
    if [ "$DONE_SHARDS" -ge "$N_SHARDS" ] && [ "$FORCE_RUN" != "1" ]; then
        echo "  SKIP ${WS}: all ${N_SHARDS} shard diagnostics exist"
        continue
    fi
    if bjobs -J "isotwas_${WS}" 2>/dev/null | grep -q "isotwas_${WS}"; then
        echo "  SKIP ${WS}: array already running/pending"
        continue
    fi

    ENV_STR="CONFIG=${CONFIG},SCRIPTS_DIR=${SCRIPTS_DIR},REPO_ROOT=${REPO_ROOT},QTL_DIR=${QTL_DIR},ISOTWAS_DIR=${ISOTWAS_DIR},WEIGHT_SET=${WS},N_SHARDS=${N_SHARDS},R2_MIN=${R2_MIN},FORCE_RUN=${FORCE_RUN},RSCRIPT=${RSCRIPT}"
    bsub -J "isotwas_${WS}[${ARRAY_SPEC}]" -q "$QUEUE" -n "$THREADS" -W "$WALLTIME" \
         -o "${LOG_DIR}/isotwas_${WS}_%J_%I.out" \
         -e "${LOG_DIR}/isotwas_${WS}_%J_%I.err" \
         -env "$ENV_STR" \
         < "${SCRIPTS_DIR}/47a_run_isotwas_shard.sh"
    echo "  ${WS}: ${N_GENES} genes -> array of ${N_SHARDS} shards (submitted: ${ARRAY_SPEC})"
    N=$((N + 1))
    if [ "$TEST" = "1" ]; then
        echo ""
        echo "TEST=1: submitted one pilot shard (isotwas_${WS}[1] of ${N_SHARDS}, ~SHARD_SIZE genes)."
        echo "Check per-gene wall times before submitting the rest:"
        echo "  tail ${LOG_DIR}/isotwas_${WS}_<jobid>_1.out"
        exit 0
    fi
done
echo ""
echo "Submitted isoTWAS training arrays for $N weight sets."
