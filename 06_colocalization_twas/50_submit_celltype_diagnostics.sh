#!/bin/bash
# =============================================================================
# 50_submit_celltype_diagnostics.sh — submit Module-06 GxC QC diagnostics
# =============================================================================
# Runs ancestry-specific GxC power/calibration diagnostics for the stage-49
# cell-type attribution results.  EUR N=136 is treated as the known full EUR
# cohort size; the centerpiece is repeated EAS -> EUR-sized downsampling.
#
# Prerequisites:
#   - stage 49 completed, including:
#       $RESULTS_DIR/coloc/aggregated/xqtl_celltype_interactions.tsv.gz
#       $RESULTS_DIR/coloc/aggregated/xqtl_celltype_annotation.tsv.gz
#   - qtl_inputs from module 05
#   - tensorqtl conda env (pandas/numpy/scipy)
#   - `module load plink` provides the `plink2` command
#
# Usage:
#   TEST=1 bash 50_submit_celltype_diagnostics.sh
#   bash 50_submit_celltype_diagnostics.sh
#
# Main overrides:
#   THREADS=8          worker processes, default 8
#   MEM_GB=32          LSF memory request, default 32 GB
#   QUEUE=medium       default medium
#   WALLTIME=08:00     default 8 hours
#   PERMUTATIONS=100   residual permutations per selected test; 0 disables
#   PERMUTE_TOP=20     top nominal tests per ancestry x modality
#   PERMUTE_RANDOM=40  random tests per ancestry x modality
#   HC3=1              1=run robust-SE refits; 0=skip
#   DOWNSAMPLE_REPS=100 repeated EAS -> EUR-sized exact refits; 0 disables
#   DOWNSAMPLE_N=       optional target N; default is EUR_qtl.psam count (136)
#   DOWNSAMPLE_ANCESTRY=EAS
#   DOWNSAMPLE_REFERENCE=EUR
#   FORCE=1            rebuild even if the final bundle exists
#   ANCESTRIES="EAS EUR"
#   OUT_DIR=...        default RESULTS_DIR/coloc/aggregated/celltype_diagnostics
#
# The final file to upload for review is:
#   $OUT_DIR/celltype_diagnostics_bundle.zip
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
THREADS="${THREADS:-8}"
MEM_GB="${MEM_GB:-32}"
QUEUE="${QUEUE:-medium}"
WALLTIME="${WALLTIME:-08:00}"
PERMUTATIONS="${PERMUTATIONS:-100}"
PERMUTE_TOP="${PERMUTE_TOP:-20}"
PERMUTE_RANDOM="${PERMUTE_RANDOM:-40}"
HC3="${HC3:-1}"
DOWNSAMPLE_REPS="${DOWNSAMPLE_REPS:-100}"
DOWNSAMPLE_N="${DOWNSAMPLE_N:-}"
DOWNSAMPLE_ANCESTRY="${DOWNSAMPLE_ANCESTRY:-EAS}"
DOWNSAMPLE_REFERENCE="${DOWNSAMPLE_REFERENCE:-EUR}"
FORCE="${FORCE:-0}"
TEST="${TEST:-0}"
PYENV="${PYENV:-tensorqtl}"

mkdir -p "$LOG_DIR" "$OUT_DIR"
BUNDLE="${OUT_DIR}/celltype_diagnostics_bundle.zip"
if [ -s "$BUNDLE" ] && [ "$FORCE" != "1" ]; then
    echo "SKIP: diagnostic bundle already exists: $BUNDLE"
    exit 0
fi

INTER="${RESULTS_DIR}/coloc/aggregated/xqtl_celltype_interactions.tsv.gz"
ANN="${RESULTS_DIR}/coloc/aggregated/xqtl_celltype_annotation.tsv.gz"
for f in "$INTER" "$ANN"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: required stage-49 output missing: $f" >&2
        exit 1
    fi
done

if [ "$TEST" = "1" ]; then
    # Fast smoke test: EAS only, no permutations, no HC3.
    ANCESTRIES="EAS"
    PERMUTATIONS=0
    HC3=0
    DOWNSAMPLE_REPS=2
    DOWNSAMPLE_N=40
    THREADS=2
    WALLTIME="04:00"
    echo "TEST=1: EAS only, HC3/permutations off, 2 downsample replicates to N=40"
fi

JOB="celltype_diag"
if bjobs -J "$JOB" 2>/dev/null | grep -q "$JOB"; then
    echo "SKIP: $JOB already running/pending"
    exit 0
fi

# LSF's -env parser rejects empty NAME= entries on this cluster.  Keep the
# environment list non-empty-valued and append optional DOWNSAMPLE_N only when
# explicitly set.  Boolean command-line flags are reconstructed in the worker
# from HC3/FORCE rather than passed as possibly-empty strings.
ENV_STR="CONFIG=${CONFIG},SCRIPTS_DIR=${SCRIPTS_DIR},QTL_DIR=${QTL_DIR},RESULTS_DIR=${RESULTS_DIR},OUT_DIR=${OUT_DIR},ANCESTRIES=${ANCESTRIES},THREADS=${THREADS},PERMUTATIONS=${PERMUTATIONS},PERMUTE_TOP=${PERMUTE_TOP},PERMUTE_RANDOM=${PERMUTE_RANDOM},DOWNSAMPLE_REPS=${DOWNSAMPLE_REPS},DOWNSAMPLE_ANCESTRY=${DOWNSAMPLE_ANCESTRY},DOWNSAMPLE_REFERENCE=${DOWNSAMPLE_REFERENCE},HC3=${HC3},FORCE=${FORCE},PYENV=${PYENV}"
if [ -n "$DOWNSAMPLE_N" ]; then
    ENV_STR="${ENV_STR},DOWNSAMPLE_N=${DOWNSAMPLE_N}"
fi

bsub -J "$JOB" -q "$QUEUE" -n "$THREADS" -W "$WALLTIME" \
     -M "${MEM_GB}G" -R "rusage[mem=${MEM_GB}G]" \
     -o "${LOG_DIR}/${JOB}.%J.out" -e "${LOG_DIR}/${JOB}.%J.err" \
     -env "$ENV_STR" <<'EOF'
#!/bin/bash
set -eo pipefail
source /etc/profile.d/modules.sh
if [ "$PYENV" != "none" ]; then
    eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
    conda activate "$PYENV"
fi
module load plink
command -v plink2 >/dev/null || { echo "ERROR: plink2 not found after module load plink" >&2; exit 1; }
python3 - <<'PY'
import numpy, pandas, scipy
print("python diagnostics deps OK", numpy.__version__, pandas.__version__, scipy.__version__)
PY
DOWNSAMPLE_N_FLAG=()
[ -n "${DOWNSAMPLE_N:-}" ] && DOWNSAMPLE_N_FLAG=(--downsample-n "$DOWNSAMPLE_N")
HC3_FLAG=()
[ "$HC3" = "0" ] && HC3_FLAG=(--skip-hc3)
FORCE_FLAG=()
[ "$FORCE" = "1" ] && FORCE_FLAG=(--force)
python3 "${SCRIPTS_DIR}/50_celltype_diagnostics.py" \
    --results-dir "$RESULTS_DIR" \
    --qtl-dir "$QTL_DIR" \
    --out-dir "$OUT_DIR" \
    --ancestries $ANCESTRIES \
    --workers "$THREADS" \
    --permutations "$PERMUTATIONS" \
    --permute-top "$PERMUTE_TOP" \
    --permute-random "$PERMUTE_RANDOM" \
    --downsample-reps "$DOWNSAMPLE_REPS" \
    --downsample-ancestry "$DOWNSAMPLE_ANCESTRY" \
    --downsample-reference-ancestry "$DOWNSAMPLE_REFERENCE" \
    "${DOWNSAMPLE_N_FLAG[@]}" \
    "${HC3_FLAG[@]}" "${FORCE_FLAG[@]}"
EOF

echo "Submitted $JOB"
echo "  workers:      $THREADS"
echo "  HC3:          $HC3"
echo "  permutations: $PERMUTATIONS/test; top=$PERMUTE_TOP random=$PERMUTE_RANDOM per stratum"
echo "  downsampling: $DOWNSAMPLE_ANCESTRY -> $DOWNSAMPLE_REFERENCE-sized; reps=$DOWNSAMPLE_REPS target_N=${DOWNSAMPLE_N:-auto}"
echo "  output:       $OUT_DIR"
echo "When complete, upload:"
echo "  $BUNDLE"
