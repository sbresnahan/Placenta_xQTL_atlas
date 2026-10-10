#!/bin/bash
# =============================================================================
# 42_submit_coloc.sh — Submit SuSiE-coloc array jobs (one array per modality)
# =============================================================================
# Objective 1.6: pairwise SuSiE-coloc of every fine-mapped xQTL locus against
# every ancestry-matched GWAS in the catalog.
#
# Step 1 builds the per-modality task lists (40_prepare_coloc_loci.py) unless
# they already exist. Step 2 submits one LSF job array per modality; each
# array element runs 42a_run_coloc_shard.sh, which prepares its round-robin
# task shard on the host and then runs the pure-R 41_susie_coloc.R worker.
# Per-task .done sentinels make reruns
# incremental, so it is safe to resubmit at any time.
#
# Prerequisites:
#   - 36_install_coloc_env.sh completed (coloc/susieR in the R library)
#   - 37/38 completed (or manual GWAS placed) for every cataloged trait
#   - 39 completed: {RESULTS_DIR}/nominal/{ANC}/{ANC}_{MOD}.nominal.tsv.gz+.tbi
#   - 40 completed (or run here): {RESULTS_DIR}/coloc/loci/{MOD}.tasks.tsv
#   - compute-node modules 'samtools' (tabix) and 'plink' (plink2); these are
#     used by host-side preprocessing, never invoked from R/Singularity
#
# Usage:
#   TEST=1 bash 42_submit_coloc.sh   # build task lists + submit ONE pilot shard
#   bash 42_submit_coloc.sh          # submit all arrays
#
# Optional overrides (environment):
#   MODALITIES   — default all 9 per-modality layers (combined excluded)
#   SHARD_SIZE   — tasks per array element, default 25 (measure with TEST=1)
#   QUEUE        — default medium
#   WALLTIME     — default 04:00
#   THREADS      — default 2
#   MEM_GB       — LSF memory request per shard, default 8 GB
#   PLINK_MEMORY_MB — memory cap passed to plink2, default 2048 MiB
#   MIN_VARIANTS — default 50
#   PP_H4        — colocalization call threshold, default 0.7
#   FORCE_TASKS=1 — rebuild task lists even if present
#   FORCE_RUN=1   — re-run tasks with existing .done markers (worker deletes
#                   them first via FORCE_RUN passthrough)
#   PYENV       — conda env providing pandas for step 40 (default tensorqtl;
#                 PYENV=none skips activation)
#
# TEST=1 submits only array index 1 while keeping the true N_SHARDS in the
# worker environment, so the pilot processes ~SHARD_SIZE tasks (NOT the whole
# modality).
# =============================================================================
set -eo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"
REPO_ROOT="$(cd "${SCRIPTS_DIR}/.." && pwd)"

# --- Python environment ------------------------------------------------------
# Step 1 runs 40_prepare_coloc_loci.py, which requires pandas (absent from
# the login-node system python3). Activate the pipeline conda env.
PYENV="${PYENV:-tensorqtl}"
if [ "$PYENV" != "none" ]; then
    source /etc/profile.d/modules.sh
    eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
    conda activate "$PYENV"
fi
if ! python3 -c "import pandas" 2>/dev/null; then
    echo "ERROR: python3 cannot import pandas (needed by 40_prepare_coloc_loci.py)." >&2
    echo "  conda activate tensorqtl   (module-05 script 21), or PYENV=<env>" >&2
    exit 1
fi
OUTPUT_BASE="${OUTPUT_BASE:-$(dirname "$CONFIG")}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
LOG_DIR="${LOG_DIR:-${OUTPUT_BASE}/logs}"
COLOC_DIR="${COLOC_DIR:-${RESULTS_DIR}/coloc}"
MODALITIES="${MODALITIES:-expression isoforms isoform_expression splicing intron_retention alt_TSS alt_polyA RNA_editing stability}"
SHARD_SIZE="${SHARD_SIZE:-25}"
QUEUE="${QUEUE:-medium}"
WALLTIME="${WALLTIME:-04:00}"
THREADS="${THREADS:-2}"
MEM_GB="${MEM_GB:-8}"
PLINK_MEMORY_MB="${PLINK_MEMORY_MB:-2048}"
MIN_VARIANTS="${MIN_VARIANTS:-50}"
PP_H4="${PP_H4:-0.7}"
TEST="${TEST:-0}"
FORCE_TASKS="${FORCE_TASKS:-0}"
FORCE_RUN="${FORCE_RUN:-0}"
RSCRIPT="${RSCRIPT:-${REPO_ROOT}/bin/Rscript_sif}"

# GWAS + 1KG reference locations (38_harmonize_gwas.py output; 1KG pgen)
GWAS_DIR="${GWAS_DIR:-${OUTPUT_BASE}/gwas}"
KG_PGEN="${KG_PGEN:-/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/1kGP/1kGP_hg38}"
KG_SAMPLE_MAP="${KG_SAMPLE_MAP:-${GWAS_DIR}/1kg_sample_superpop.tsv}"

mkdir -p "$COLOC_DIR/loci" "$LOG_DIR"

echo "=== 42_submit_coloc.sh ==="
echo "  MODALITIES:  $MODALITIES"
echo "  COLOC_DIR:   $COLOC_DIR"
echo "  GWAS_DIR:    $GWAS_DIR"
echo "  KG_PGEN:     $KG_PGEN"
echo "  SHARD_SIZE:  $SHARD_SIZE   QUEUE: $QUEUE   WALLTIME: $WALLTIME"
echo "  MEMORY:      ${MEM_GB}G LSF; ${PLINK_MEMORY_MB} MiB plink2 cap"

# --- Step 1: per-modality task lists ----------------------------------------
for MOD in $MODALITIES; do
    TASKS="${COLOC_DIR}/loci/${MOD}.tasks.tsv"
    if [ -f "$TASKS" ] && [ "$FORCE_TASKS" != "1" ]; then
        echo "  task list exists: $TASKS ($(($(wc -l < "$TASKS") - 1)) tasks)"
        continue
    fi
    echo "  building task list: $MOD"
    python3 "${SCRIPTS_DIR}/40_prepare_coloc_loci.py" \
        --results-dir "$RESULTS_DIR" \
        --qtl-dir "$QTL_DIR" \
        --gwas-dir "$GWAS_DIR" \
        --catalog "${SCRIPTS_DIR}/gwas_catalog.tsv" \
        --kg-pgen "$KG_PGEN" \
        --kg-sample-map "$KG_SAMPLE_MAP" \
        --modalities "$MOD" \
        --out-dir "$COLOC_DIR"
done

# --- Step 2: submit one array per modality -----------------------------------
N=0
for MOD in $MODALITIES; do
    TASKS="${COLOC_DIR}/loci/${MOD}.tasks.tsv"
    if [ ! -s "$TASKS" ]; then
        echo "  SKIP ${MOD}: no task list (missing loci or GWAS inputs)"
        continue
    fi
    N_TASKS=$(($(wc -l < "$TASKS") - 1))
    if [ "$N_TASKS" -le 0 ]; then
        echo "  SKIP ${MOD}: 0 tasks"
        continue
    fi
    N_SHARDS=$(( (N_TASKS + SHARD_SIZE - 1) / SHARD_SIZE ))
    # TEST=1: submit ONLY array index 1 but keep the true N_SHARDS in the
    # worker environment. The host preparer assigns tasks round-robin by
    # N_SHARDS, so N_SHARDS=1 would assign EVERY task to the pilot shard.
    ARRAY_SPEC="1-${N_SHARDS}"
    if [ "$TEST" = "1" ]; then ARRAY_SPEC="1"; fi

    # skip if every shard already has diagnostics (all tasks attempted)
    DONE_SHARDS=0
    for DIAG in "${COLOC_DIR}/diagnostics/${MOD}".shard-*.diagnostics.tsv; do
        [ -f "$DIAG" ] && DONE_SHARDS=$((DONE_SHARDS + 1))
    done
    if [ "$DONE_SHARDS" -ge "$N_SHARDS" ] && [ "$FORCE_RUN" != "1" ]; then
        echo "  SKIP ${MOD}: all ${N_SHARDS} shard diagnostics exist"
        continue
    fi
    if bjobs -J "coloc_${MOD}" 2>/dev/null | grep -q "coloc_${MOD}"; then
        echo "  SKIP ${MOD}: array already running/pending"
        continue
    fi

    ENV_STR="CONFIG=${CONFIG},SCRIPTS_DIR=${SCRIPTS_DIR},REPO_ROOT=${REPO_ROOT},RESULTS_DIR=${RESULTS_DIR},COLOC_DIR=${COLOC_DIR},TASKS=${TASKS},N_SHARDS=${N_SHARDS},MIN_VARIANTS=${MIN_VARIANTS},PP_H4=${PP_H4},FORCE_RUN=${FORCE_RUN},KEEP_PREP=${KEEP_PREP:-0},PLINK_MEMORY_MB=${PLINK_MEMORY_MB},RSCRIPT=${RSCRIPT}"
    bsub -J "coloc_${MOD}[${ARRAY_SPEC}]" -q "$QUEUE" -n "$THREADS" -W "$WALLTIME" \
         -M "${MEM_GB}G" -R "rusage[mem=${MEM_GB}G]" \
         -o "${LOG_DIR}/coloc_${MOD}_%J_%I.out" \
         -e "${LOG_DIR}/coloc_${MOD}_%J_%I.err" \
         -env "$ENV_STR" \
         < "${SCRIPTS_DIR}/42a_run_coloc_shard.sh"
    echo "  ${MOD}: ${N_TASKS} tasks -> array of ${N_SHARDS} shards (submitted: ${ARRAY_SPEC})"
    N=$((N + 1))
    if [ "$TEST" = "1" ]; then
        echo ""
        echo "TEST=1: submitted one pilot shard (coloc_${MOD}[1] of ${N_SHARDS})."
        echo "The pilot processes ~SHARD_SIZE tasks. Check per-task wall times"
        echo "before submitting the rest:"
        echo "  tail ${LOG_DIR}/coloc_${MOD}_<jobid>_1.out"
        exit 0
    fi
done
echo ""
echo "Submitted arrays for $N modalities. Aggregate when done:"
echo "  python3 ${SCRIPTS_DIR}/49_aggregate_coloc_twas.py \\"
echo "      --results-dir ${RESULTS_DIR} --qtl-dir ${QTL_DIR}"
