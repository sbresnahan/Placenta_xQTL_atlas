#!/bin/bash
# =============================================================================
# 45_submit_colocboost.sh — Submit colocBoost region array jobs (per ancestry)
# =============================================================================
# Objective 1.6: multi-trait colocBoost over merged fine-mapping regions.
#
# Step 1 builds region/outcome manifests (43_prepare_colocboost.py) unless
# present. Step 2 submits one LSF job array per ancestry; each element runs
# 45a_run_colocboost_shard.sh -> 44_colocboost.R on a round-robin shard of
# the region list. Per-region .done sentinels make reruns incremental.
#
# Prerequisites:
#   - 36_install_coloc_env.sh completed (colocboost in the R library)
#   - module 05 fine-mapping loci: {RESULTS_DIR}/finemap/loci/{MOD}.loci.tsv
#   - 39 completed: merged nominal TSVs per ancestry x modality
#   - 37/38 completed (or manual GWAS placed)
#   - 1KG pgen + sample->superpopulation map
#
# Usage:
#   TEST=1 bash 45_submit_colocboost.sh   # build manifests + ONE pilot shard
#   bash 45_submit_colocboost.sh          # submit all arrays
#
# Optional overrides (environment):
#   ANCESTRIES   — default "EAS EUR"
#   MODALITIES   — default all 9 per-modality layers
#   SHARD_SIZE   — regions per array element, default 10 (colocBoost regions
#                  are heavier than pairwise coloc tasks; measure with TEST=1)
#   QUEUE        — default medium
#   WALLTIME     — default 06:00
#   THREADS      — default 2
#   MERGE_GAP    — locus-merge gap in bp, default 100000
#   FORCE_PREP=1 — rebuild manifests even if present
#   FORCE_RUN=1  — re-run regions with existing .done markers
#   PYENV      — conda env providing pandas for step 43 (default tensorqtl;
#                PYENV=none skips activation)
#
# TEST=1 submits only array index 1 while keeping the true N_SHARDS in the
# worker environment, so the pilot processes ~SHARD_SIZE regions (NOT every
# region in the ancestry).
# =============================================================================
set -eo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"
REPO_ROOT="$(cd "${SCRIPTS_DIR}/.." && pwd)"

# --- Python environment ------------------------------------------------------
# Step 1 runs 43_prepare_colocboost.py, which requires pandas (absent from
# the login-node system python3). Activate the pipeline conda env.
PYENV="${PYENV:-tensorqtl}"
if [ "$PYENV" != "none" ]; then
    source /etc/profile.d/modules.sh
    eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
    conda activate "$PYENV"
fi
if ! python3 -c "import pandas" 2>/dev/null; then
    echo "ERROR: python3 cannot import pandas (needed by 43_prepare_colocboost.py)." >&2
    echo "  conda activate tensorqtl   (module-05 script 21), or PYENV=<env>" >&2
    exit 1
fi
OUTPUT_BASE="${OUTPUT_BASE:-$(dirname "$CONFIG")}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
LOG_DIR="${LOG_DIR:-${OUTPUT_BASE}/logs}"
COLOC_DIR="${COLOC_DIR:-${RESULTS_DIR}/coloc}"
CB_DIR="${CB_DIR:-${COLOC_DIR}/colocboost}"
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
MODALITIES="${MODALITIES:-expression isoforms isoform_expression splicing intron_retention alt_TSS alt_polyA RNA_editing stability}"
SHARD_SIZE="${SHARD_SIZE:-10}"
QUEUE="${QUEUE:-medium}"
WALLTIME="${WALLTIME:-06:00}"
THREADS="${THREADS:-2}"
MERGE_GAP="${MERGE_GAP:-100000}"
TEST="${TEST:-0}"
FORCE_PREP="${FORCE_PREP:-0}"
FORCE_RUN="${FORCE_RUN:-0}"
RSCRIPT="${RSCRIPT:-${REPO_ROOT}/bin/Rscript_sif}"

GWAS_DIR="${GWAS_DIR:-${OUTPUT_BASE}/gwas}"
KG_PGEN="${KG_PGEN:-/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/1kGP/1kGP_hg38}"
KG_SAMPLE_MAP="${KG_SAMPLE_MAP:-${GWAS_DIR}/1kg_sample_superpop.tsv}"

mkdir -p "$CB_DIR" "$LOG_DIR"

echo "=== 45_submit_colocboost.sh ==="
echo "  ANCESTRIES:  $ANCESTRIES"
echo "  CB_DIR:      $CB_DIR"
echo "  SHARD_SIZE:  $SHARD_SIZE   QUEUE: $QUEUE   WALLTIME: $WALLTIME"

# --- Step 1: region/outcome manifests ---------------------------------------
if [ -f "${CB_DIR}/manifest.done" ] && [ "$FORCE_PREP" != "1" ]; then
    echo "  manifests exist (manifest.done); skipping prep"
else
    python3 "${SCRIPTS_DIR}/43_prepare_colocboost.py" \
        --results-dir "$RESULTS_DIR" \
        --qtl-dir "$QTL_DIR" \
        --gwas-dir "$GWAS_DIR" \
        --catalog "${SCRIPTS_DIR}/gwas_catalog.tsv" \
        --ancestries $ANCESTRIES \
        --modalities $MODALITIES \
        --merge-gap "$MERGE_GAP" \
        --out-dir "$CB_DIR"
    touch "${CB_DIR}/manifest.done"
fi

# --- Step 2: submit one array per ancestry -----------------------------------
N=0
for ANC in $ANCESTRIES; do
    REGIONS="${CB_DIR}/${ANC}.regions.tsv"
    if [ ! -s "$REGIONS" ]; then
        echo "  SKIP ${ANC}: no regions file"
        continue
    fi
    N_REGIONS=$(($(wc -l < "$REGIONS") - 1))
    if [ "$N_REGIONS" -le 0 ]; then
        echo "  SKIP ${ANC}: 0 regions"
        continue
    fi
    N_SHARDS=$(( (N_REGIONS + SHARD_SIZE - 1) / SHARD_SIZE ))
    # TEST=1: submit ONLY array index 1 but keep the true N_SHARDS in the
    # worker environment (round-robin assignment; N_SHARDS=1 would give the
    # pilot shard EVERY region).
    ARRAY_SPEC="1-${N_SHARDS}"
    if [ "$TEST" = "1" ]; then ARRAY_SPEC="1"; fi

    DONE_SHARDS=0
    for DIAG in "${CB_DIR}/diagnostics/${ANC}".shard-*.diagnostics.tsv; do
        [ -f "$DIAG" ] && DONE_SHARDS=$((DONE_SHARDS + 1))
    done
    if [ "$DONE_SHARDS" -ge "$N_SHARDS" ] && [ "$FORCE_RUN" != "1" ]; then
        echo "  SKIP ${ANC}: all ${N_SHARDS} shard diagnostics exist"
        continue
    fi
    if bjobs -J "colocboost_${ANC}" 2>/dev/null | grep -q "colocboost_${ANC}"; then
        echo "  SKIP ${ANC}: array already running/pending"
        continue
    fi

    ENV_STR="CONFIG=${CONFIG},SCRIPTS_DIR=${SCRIPTS_DIR},REPO_ROOT=${REPO_ROOT},CB_DIR=${CB_DIR},REGIONS=${REGIONS},OUTCOMES=${CB_DIR}/${ANC}.outcomes.tsv,N_SHARDS=${N_SHARDS},LD_XQTL_PGEN=${QTL_DIR}/${ANC}_qtl,LD_GWAS_PGEN=${KG_PGEN},LD_GWAS_KEEP=${COLOC_DIR}/loci/${ANC}.1kg.keep,FORCE_RUN=${FORCE_RUN},RSCRIPT=${RSCRIPT}"
    bsub -J "colocboost_${ANC}[${ARRAY_SPEC}]" -q "$QUEUE" -n "$THREADS" -W "$WALLTIME" \
         -o "${LOG_DIR}/colocboost_${ANC}_%J_%I.out" \
         -e "${LOG_DIR}/colocboost_${ANC}_%J_%I.err" \
         -env "$ENV_STR" \
         < "${SCRIPTS_DIR}/45a_run_colocboost_shard.sh"
    echo "  ${ANC}: ${N_REGIONS} regions -> array of ${N_SHARDS} shards (submitted: ${ARRAY_SPEC})"
    N=$((N + 1))
    if [ "$TEST" = "1" ]; then
        echo ""
        echo "TEST=1: submitted one pilot shard (colocboost_${ANC}[1] of ${N_SHARDS})."
        echo "The pilot processes ~SHARD_SIZE regions. Check per-region wall"
        echo "times before submitting the rest:"
        echo "  tail ${LOG_DIR}/colocboost_${ANC}_<jobid>_1.out"
        exit 0
    fi
done
echo ""
echo "Submitted colocBoost arrays for $N ancestries."
