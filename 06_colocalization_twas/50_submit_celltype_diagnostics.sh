#!/bin/bash
# =============================================================================
# 50_submit_celltype_diagnostics.sh — LSF-array Module-06 GxC diagnostics
# =============================================================================
# Architecture:
#   1) build deterministic phenotype-shard plans locally
#   2) one small preparation job caches PLINK dosages + selected BED phenotypes
#   3) ONE LSF JOB ARRAY PER MODALITY; array elements are phenotype shards
#      containing both ancestries' shards for that modality
#   4) one modality aggregation job after each array; this is where downsample
#      BH/FDR is computed across ALL phenotype shards for that modality
#   5) one finalization job after all modality aggregators creates the bundle
#
# This replaces the old single 12-core job that looped serially over modalities.
# It also uses the v3 batched FWL downsampling engine, avoiding millions of
# independent Python/lstsq calls.
#
# Typical production run:
#   DOWNSAMPLE_REPS=100 SHARD_FEATURES=100 CPUS_PER_SHARD=2 \
#     ARRAY_CONCURRENCY=8 bash 50_submit_celltype_diagnostics.sh
#
# Main overrides:
#   SHARD_FEATURES=100       unique xQTL features per array element
#   CPUS_PER_SHARD=2         local workers inside each phenotype shard
#   ARRAY_CONCURRENCY=8      max simultaneously running elements PER modality
#   SHARD_MEM_GB=8           memory per array element
#   SHARD_WALLTIME=04:00     walltime per array element
#   DOWNSAMPLE_REPS=100
#   PERMUTATIONS=100
#   HC3=1
#   FORCE=1                  rebuild existing cache/shards/aggregates
# =============================================================================
set -eo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"
OUTPUT_BASE="${OUTPUT_BASE:-$(dirname "$CONFIG")}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
LOG_DIR="${LOG_DIR:-${OUTPUT_BASE}/logs}"
OUT_DIR="${OUT_DIR:-${RESULTS_DIR}/coloc/aggregated/celltype_diagnostics}"
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
PYENV="${PYENV:-tensorqtl}"
QUEUE="${QUEUE:-medium}"

SHARD_FEATURES="${SHARD_FEATURES:-100}"
CPUS_PER_SHARD="${CPUS_PER_SHARD:-2}"
ARRAY_CONCURRENCY="${ARRAY_CONCURRENCY:-8}"
SHARD_MEM_GB="${SHARD_MEM_GB:-8}"
SHARD_WALLTIME="${SHARD_WALLTIME:-04:00}"
PREP_MEM_GB="${PREP_MEM_GB:-16}"
PREP_WALLTIME="${PREP_WALLTIME:-02:00}"
AGG_MEM_GB="${AGG_MEM_GB:-12}"
AGG_WALLTIME="${AGG_WALLTIME:-02:00}"
FINAL_MEM_GB="${FINAL_MEM_GB:-12}"
FINAL_WALLTIME="${FINAL_WALLTIME:-02:00}"

PERMUTATIONS="${PERMUTATIONS:-100}"
PERMUTE_TOP="${PERMUTE_TOP:-20}"
PERMUTE_RANDOM="${PERMUTE_RANDOM:-40}"
HC3="${HC3:-1}"
DOWNSAMPLE_REPS="${DOWNSAMPLE_REPS:-100}"
DOWNSAMPLE_N="${DOWNSAMPLE_N:-}"
DOWNSAMPLE_ANCESTRY="${DOWNSAMPLE_ANCESTRY:-EAS}"
DOWNSAMPLE_REFERENCE="${DOWNSAMPLE_REFERENCE:-EUR}"
MIN_N="${MIN_N:-20}"
SEED="${SEED:-20261010}"
FORCE="${FORCE:-0}"
TEST="${TEST:-0}"

mkdir -p "$LOG_DIR" "$OUT_DIR/work"
BUNDLE="${OUT_DIR}/celltype_diagnostics_bundle.zip"
if [ -s "$BUNDLE" ] && [ "$FORCE" != "1" ]; then
    echo "SKIP: diagnostic bundle already exists: $BUNDLE"
    exit 0
fi

if [ -n "${THREADS:-}" ]; then
    echo "NOTE: THREADS is no longer the main parallelism control."
    echo "      Using CPUS_PER_SHARD=${CPUS_PER_SHARD}; ARRAY_CONCURRENCY=${ARRAY_CONCURRENCY}."
fi

if [ "$TEST" = "1" ]; then
    ANCESTRIES="EAS"
    SHARD_FEATURES=10
    CPUS_PER_SHARD=2
    ARRAY_CONCURRENCY=2
    PERMUTATIONS=0
    HC3=0
    DOWNSAMPLE_REPS=2
    DOWNSAMPLE_N=40
    echo "TEST=1: EAS only; 10 features/shard; 2 downsample reps to N=40"
fi

INTER="${RESULTS_DIR}/coloc/aggregated/xqtl_celltype_interactions.tsv.gz"
ANN="${RESULTS_DIR}/coloc/aggregated/xqtl_celltype_annotation.tsv.gz"
for f in "$INTER" "$ANN"; do
    [ -s "$f" ] || { echo "ERROR: required stage-49 output missing: $f" >&2; exit 1; }
done

FORCE_FLAG=()
[ "$FORCE" = "1" ] && FORCE_FLAG=(--force)
HC3_FLAG=()
[ "$HC3" = "0" ] && HC3_FLAG=(--skip-hc3)
DOWNSAMPLE_N_FLAG=()
[ -n "$DOWNSAMPLE_N" ] && DOWNSAMPLE_N_FLAG=(--downsample-n "$DOWNSAMPLE_N")

# The plan is cheap to create and is needed before LSF arrays can be sized.
python3 "${SCRIPTS_DIR}/50a_prepare_celltype_diagnostics.py" \
    --results-dir "$RESULTS_DIR" --qtl-dir "$QTL_DIR" --out-dir "$OUT_DIR" \
    --ancestries $ANCESTRIES --shard-features "$SHARD_FEATURES" --plan-only \
    "${FORCE_FLAG[@]}"

PLAN_ALL="${OUT_DIR}/work/plan_all.tsv"
[ -s "$PLAN_ALL" ] || { echo "ERROR: plan was not generated: $PLAN_ALL" >&2; exit 1; }

# Export once; LSF inherits the submission environment.  Avoid -env entirely,
# since this cluster rejects some empty-value -env syntax.
export CONFIG SCRIPTS_DIR OUTPUT_BASE QTL_DIR RESULTS_DIR LOG_DIR OUT_DIR ANCESTRIES PYENV
export PERMUTATIONS PERMUTE_TOP PERMUTE_RANDOM HC3 DOWNSAMPLE_REPS DOWNSAMPLE_N
export SHARD_FEATURES ARRAY_CONCURRENCY
export DOWNSAMPLE_ANCESTRY DOWNSAMPLE_REFERENCE MIN_N SEED FORCE CPUS_PER_SHARD

parse_jobid() {
    sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' | head -n1
}

PREP_DEP=""
if [ -s "${OUT_DIR}/work/.prepare_done" ] && [ "$FORCE" != "1" ]; then
    echo "PREP: reusing existing diagnostic input cache"
else
    echo "Submitting preparation/cache job..."
    PREP_SUBMIT=$(bsub -J "ctd_prepare" -q "$QUEUE" -n 1 -W "$PREP_WALLTIME" \
        -M "${PREP_MEM_GB}G" -R "rusage[mem=${PREP_MEM_GB}G]" \
        -o "${LOG_DIR}/ctd_prepare.%J.out" -e "${LOG_DIR}/ctd_prepare.%J.err" <<'EOF'
#!/bin/bash
set -eo pipefail
source /etc/profile.d/modules.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ "$PYENV" != "none" ]; then
    eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
    conda activate "$PYENV"
fi
module load plink
ARGS=(--results-dir "$RESULTS_DIR" --qtl-dir "$QTL_DIR" --out-dir "$OUT_DIR" --ancestries $ANCESTRIES --shard-features "$SHARD_FEATURES")
[ "$FORCE" = "1" ] && ARGS+=(--force)
python3 "${SCRIPTS_DIR}/50a_prepare_celltype_diagnostics.py" "${ARGS[@]}"
EOF
)
    echo "$PREP_SUBMIT"
    PREP_ID=$(printf '%s\n' "$PREP_SUBMIT" | parse_jobid)
    [ -n "$PREP_ID" ] || { echo "ERROR: could not parse prepare job ID" >&2; exit 1; }
    PREP_DEP="done(${PREP_ID})"
fi

AGG_IDS=()
shopt -s nullglob
for PLAN_PATH in "${OUT_DIR}/work/plans/"*.tsv; do
    MODALITY=$(basename "$PLAN_PATH" .tsv)
    NTASKS=$(( $(wc -l < "$PLAN_PATH") - 1 ))
    [ "$NTASKS" -gt 0 ] || continue
    MOD_DONE="${OUT_DIR}/work/modality/${MODALITY}/.done"
    if [ -s "$MOD_DONE" ] && [ "$FORCE" != "1" ]; then
        echo "MODALITY ${MODALITY}: aggregate already complete; skipping array"
        continue
    fi
    export PLAN_PATH MODALITY
    JOB="ctd_${MODALITY}"
    DEP_ARGS=()
    [ -n "$PREP_DEP" ] && DEP_ARGS=(-w "$PREP_DEP")
    echo "Submitting ${JOB}[1-${NTASKS}]%${ARRAY_CONCURRENCY}: ${SHARD_FEATURES} features/shard, ${CPUS_PER_SHARD} CPU/shard"
    ARRAY_SUBMIT=$(bsub -J "${JOB}[1-${NTASKS}]%${ARRAY_CONCURRENCY}" -q "$QUEUE" \
        -n "$CPUS_PER_SHARD" -W "$SHARD_WALLTIME" \
        -M "${SHARD_MEM_GB}G" -R "rusage[mem=${SHARD_MEM_GB}G]" \
        "${DEP_ARGS[@]}" \
        -o "${LOG_DIR}/${JOB}.%J.%I.out" -e "${LOG_DIR}/${JOB}.%J.%I.err" <<'EOF'
#!/bin/bash
set -eo pipefail
source /etc/profile.d/modules.sh
# Small linear-algebra batches: one BLAS thread per Python worker is much faster
# than CPUS_PER_SHARD workers each creating their own BLAS thread pool.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 MALLOC_ARENA_MAX=2
if [ "$PYENV" != "none" ]; then
    eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
    conda activate "$PYENV"
fi
ARGS=(
  --results-dir "$RESULTS_DIR" --qtl-dir "$QTL_DIR" --out-dir "$OUT_DIR"
  --plan "$PLAN_PATH" --array-index "$LSB_JOBINDEX" --workers "$CPUS_PER_SHARD"
  --min-n "$MIN_N" --permutations "$PERMUTATIONS" --permute-top "$PERMUTE_TOP"
  --permute-random "$PERMUTE_RANDOM" --downsample-reps "$DOWNSAMPLE_REPS"
  --downsample-ancestry "$DOWNSAMPLE_ANCESTRY"
  --downsample-reference-ancestry "$DOWNSAMPLE_REFERENCE" --seed "$SEED"
)
[ -n "${DOWNSAMPLE_N:-}" ] && ARGS+=(--downsample-n "$DOWNSAMPLE_N")
[ "$HC3" = "0" ] && ARGS+=(--skip-hc3)
[ "$FORCE" = "1" ] && ARGS+=(--force)
python3 "${SCRIPTS_DIR}/50b_run_celltype_diag_shard.py" "${ARGS[@]}"
EOF
)
    echo "$ARRAY_SUBMIT"
    ARRAY_ID=$(printf '%s\n' "$ARRAY_SUBMIT" | parse_jobid)
    [ -n "$ARRAY_ID" ] || { echo "ERROR: could not parse array job ID for $MODALITY" >&2; exit 1; }

    AGG_JOB="ctd_${MODALITY}_agg"
    AGG_SUBMIT=$(bsub -J "$AGG_JOB" -q "$QUEUE" -n 1 -W "$AGG_WALLTIME" \
        -M "${AGG_MEM_GB}G" -R "rusage[mem=${AGG_MEM_GB}G]" \
        -w "done(${ARRAY_ID})" \
        -o "${LOG_DIR}/${AGG_JOB}.%J.out" -e "${LOG_DIR}/${AGG_JOB}.%J.err" <<'EOF'
#!/bin/bash
set -eo pipefail
source /etc/profile.d/modules.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ "$PYENV" != "none" ]; then
    eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
    conda activate "$PYENV"
fi
ARGS=(--out-dir "$OUT_DIR" --plan "$PLAN_PATH" --modality "$MODALITY"
      --downsample-ancestry "$DOWNSAMPLE_ANCESTRY"
      --downsample-reference-ancestry "$DOWNSAMPLE_REFERENCE")
[ "$FORCE" = "1" ] && ARGS+=(--force)
python3 "${SCRIPTS_DIR}/50c_aggregate_celltype_diag_modality.py" "${ARGS[@]}"
EOF
)
    echo "$AGG_SUBMIT"
    AGG_ID=$(printf '%s\n' "$AGG_SUBMIT" | parse_jobid)
    [ -n "$AGG_ID" ] || { echo "ERROR: could not parse aggregate job ID for $MODALITY" >&2; exit 1; }
    AGG_IDS+=("$AGG_ID")
done
shopt -u nullglob

FINAL_DEP=""
if [ "${#AGG_IDS[@]}" -gt 0 ]; then
    for jid in "${AGG_IDS[@]}"; do
        if [ -z "$FINAL_DEP" ]; then FINAL_DEP="done(${jid})"; else FINAL_DEP="${FINAL_DEP} && done(${jid})"; fi
    done
fi
FINAL_DEP_ARGS=()
[ -n "$FINAL_DEP" ] && FINAL_DEP_ARGS=(-w "$FINAL_DEP")
export PLAN_PATH="$PLAN_ALL"
FINAL_SUBMIT=$(bsub -J "ctd_finalize" -q "$QUEUE" -n 1 -W "$FINAL_WALLTIME" \
    -M "${FINAL_MEM_GB}G" -R "rusage[mem=${FINAL_MEM_GB}G]" \
    "${FINAL_DEP_ARGS[@]}" \
    -o "${LOG_DIR}/ctd_finalize.%J.out" -e "${LOG_DIR}/ctd_finalize.%J.err" <<'EOF'
#!/bin/bash
set -eo pipefail
source /etc/profile.d/modules.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ "$PYENV" != "none" ]; then
    eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
    conda activate "$PYENV"
fi
ARGS=(
  --results-dir "$RESULTS_DIR" --qtl-dir "$QTL_DIR" --out-dir "$OUT_DIR"
  --plan "$PLAN_PATH" --workers "$CPUS_PER_SHARD" --min-n "$MIN_N"
  --permutations "$PERMUTATIONS" --permute-top "$PERMUTE_TOP"
  --permute-random "$PERMUTE_RANDOM" --downsample-reps "$DOWNSAMPLE_REPS"
  --downsample-ancestry "$DOWNSAMPLE_ANCESTRY"
  --downsample-reference-ancestry "$DOWNSAMPLE_REFERENCE" --seed "$SEED"
)
[ -n "${DOWNSAMPLE_N:-}" ] && ARGS+=(--downsample-n "$DOWNSAMPLE_N")
[ "$HC3" = "0" ] && ARGS+=(--skip-hc3)
[ "$FORCE" = "1" ] && ARGS+=(--force)
python3 "${SCRIPTS_DIR}/50d_finalize_celltype_diagnostics.py" "${ARGS[@]}"
EOF
)
echo "$FINAL_SUBMIT"
echo
echo "Submitted sharded diagnostics."
echo "  shard size:             $SHARD_FEATURES features"
echo "  CPUs per shard:         $CPUS_PER_SHARD"
echo "  concurrency / modality: $ARRAY_CONCURRENCY"
echo "  downsample reps:        $DOWNSAMPLE_REPS"
echo "  final bundle:           $BUNDLE"
echo "Monitor with: bjobs -w | grep ctd_"
