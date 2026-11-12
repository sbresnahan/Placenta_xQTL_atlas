#!/bin/bash
# =============================================================================
# 39c_submit_cs_nominal.sh — credible-set variant x gene-feature nominal stats
# =============================================================================
# For every variant in ANY SuSHiE credible set, pull the stage-39 nominal
# slope/SE/p against every PANTRY feature (all modalities) of the
# fine-mapped phenotype's gene, per ancestry, and package an archive.
# See 39b_cs_nominal_lookup.py for the output schema.
#
# Jobs: csnom_index -> csnom_x[1-N] (one per ancestry x modality, tabix
# lookups) -> csnom_pkg (README/MANIFEST + tar.gz; runs when the array ends
# and refuses to package if any task failed unless ALLOW_MISSING=1).
#
# Prerequisites: 33_aggregate_finemap.py inputs ($RESULTS_DIR/finemap/{MOD}/
#   *.sushie.cs.tsv) and 39_run_nominal.sh merged stores
#   ($RESULTS_DIR/nominal/{ANC}/{ANC}_{MOD}.nominal.tsv.gz + .tbi).
#
# Usage:
#   bash 39c_submit_cs_nominal.sh
#   SCOPE=cis bash 39c_submit_cs_nominal.sh        # every feature tested vs
#                                                  # the variant, not just the
#                                                  # CS gene (much larger)
#   LEAD_ONLY=1 bash 39c_submit_cs_nominal.sh      # CS lead variants only
#
# Paths default to the production PANTRY tree; override with OUTPUT_BASE,
#   RESULTS_DIR, QTL_DIR, LOG_DIR; CSNOM_SCRIPTS_DIR (default: this
#   script's dir — an inherited SCRIPTS_DIR is ignored).
# Optional overrides: ANCESTRIES, MODALITIES, OUT_DIR, SCOPE (gene|cis),
#   CS_SOURCE (cs|pips), MIN_PIP, LEAD_ONLY=1, QUEUE, WALLTIME, MEM,
#   FORCE=1, ALLOW_MISSING=1, SKIP_INDEX=1 (reuse an existing index).
# =============================================================================
set -eo pipefail

# Paths: defaults are the production PANTRY tree (05_qtl_mapping/README.md).
# CONFIG is NOT used to derive paths — a CONFIG pointing elsewhere silently
# redirected every input on the first run.
# Always this script's own directory: a SCRIPTS_DIR exported by other
# projects' shell setup must not leak in (override: CSNOM_SCRIPTS_DIR).
SCRIPTS_DIR="${CSNOM_SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
OUTPUT_BASE="${OUTPUT_BASE:-/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
LOG_DIR="${LOG_DIR:-${OUTPUT_BASE}/logs}"
OUT_DIR="${OUT_DIR:-${RESULTS_DIR}/cs_nominal_lookup}"
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
MODALITIES="${MODALITIES:-expression isoforms isoform_expression splicing intron_retention alt_TSS alt_polyA RNA_editing stability}"
SCOPE="${SCOPE:-gene}"
CS_SOURCE="${CS_SOURCE:-cs}"
MIN_PIP="${MIN_PIP:-0}"
LEAD_ONLY="${LEAD_ONLY:-0}"
# seadragon runlimits: short <3h, medium 3-24h
QUEUE="${QUEUE:-medium}"
WALLTIME="${WALLTIME:-06:00}"
MEM="${MEM:-32G}"
FORCE="${FORCE:-0}"
ALLOW_MISSING="${ALLOW_MISSING:-0}"

echo "Resolved paths:"
echo "  SCRIPTS_DIR = $SCRIPTS_DIR"
echo "  RESULTS_DIR = $RESULTS_DIR"
echo "  QTL_DIR     = $QTL_DIR"
echo "  OUT_DIR     = $OUT_DIR"
echo "  LOG_DIR     = $LOG_DIR"

# --- preflight: abort before submitting anything ---------------------------
fail=0
[ -f "${SCRIPTS_DIR}/39b_cs_nominal_lookup.py" ] || {
    echo "ERROR: ${SCRIPTS_DIR}/39b_cs_nominal_lookup.py not found"; fail=1; }
[ -d "$RESULTS_DIR" ] || { echo "ERROR: RESULTS_DIR not found: $RESULTS_DIR"; fail=1; }
if [ "${SKIP_INDEX:-0}" != "1" ]; then
    if [ "$CS_SOURCE" = "pips" ]; then
        [ -f "${RESULTS_DIR}/finemap/aggregated/finemap_pips.tsv.gz" ] || {
            echo "ERROR: ${RESULTS_DIR}/finemap/aggregated/finemap_pips.tsv.gz not found"; fail=1; }
    else
        n_cs=$(find "${RESULTS_DIR}/finemap/" -mindepth 2 -maxdepth 2 \
                    -name '*.sushie.cs.tsv' 2>/dev/null | head -1 | wc -l || true)
        [ "$n_cs" -gt 0 ] || {
            echo "ERROR: no ${RESULTS_DIR}/finemap/{MOD}/*.sushie.cs.tsv files"; fail=1; }
    fi
fi
n_missing=0
for ANC in $ANCESTRIES; do
  for MOD in $MODALITIES; do
    F="${RESULTS_DIR}/nominal/${ANC}/${ANC}_${MOD}.nominal.tsv.gz"
    if [ ! -f "$F" ] || [ ! -f "${F}.tbi" ]; then
        echo "  missing nominal store or index: $F(.tbi)"
        n_missing=$((n_missing + 1))
    fi
    [ -f "${QTL_DIR}/${ANC}_${MOD}.bed.gz" ] || \
        echo "  WARNING: ${QTL_DIR}/${ANC}_${MOD}.bed.gz not found (no NA rows for untested features)"
  done
done
if [ "$n_missing" -gt 0 ] && [ "$ALLOW_MISSING" != "1" ]; then
    echo "ERROR: $n_missing nominal stores missing (run 39_run_nominal.sh, or"
    echo "       restrict ANCESTRIES/MODALITIES, or ALLOW_MISSING=1)"
    fail=1
fi
if bjobs -w 2>/dev/null | grep -q "csnom_"; then
    echo "ERROR: csnom_* jobs already pending/running (bjobs -w | grep csnom_)"
    fail=1
fi
[ "$fail" = "0" ] || exit 1

mkdir -p "$LOG_DIR" "$OUT_DIR"

# ancestry x modality task list (array index -> line); with ALLOW_MISSING=1
# tasks whose store is missing are left out
TASKS="${OUT_DIR}/tasks.txt"
: > "$TASKS"
for ANC in $ANCESTRIES; do
  for MOD in $MODALITIES; do
    F="${RESULTS_DIR}/nominal/${ANC}/${ANC}_${MOD}.nominal.tsv.gz"
    [ -f "$F" ] && [ -f "${F}.tbi" ] && echo "$ANC $MOD" >> "$TASKS"
  done
done
N=$(wc -l < "$TASKS")

ENV_STR="SCRIPTS_DIR=${SCRIPTS_DIR},QTL_DIR=${QTL_DIR},RESULTS_DIR=${RESULTS_DIR},OUT_DIR=${OUT_DIR},SCOPE=${SCOPE},CS_SOURCE=${CS_SOURCE},MIN_PIP=${MIN_PIP},LEAD_ONLY=${LEAD_ONLY},FORCE=${FORCE},ALLOW_MISSING=${ALLOW_MISSING},TASKS=${TASKS},ANCESTRIES=${ANCESTRIES// /:},MODALITIES=${MODALITIES// /:}"

DEP=""
if [ "${SKIP_INDEX:-0}" != "1" ]; then
    bsub -J csnom_index -q "$QUEUE" -n 4 -W "$WALLTIME" -M 8G -R "rusage[mem=8G]" \
         -o "${LOG_DIR}/csnom_index.%J.out" -e "${LOG_DIR}/csnom_index.%J.err" \
         -env "$ENV_STR" <<'EOF'
#!/bin/bash
set -eo pipefail
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate tensorqtl
LEAD_FLAG=""; [ "$LEAD_ONLY" = "1" ] && LEAD_FLAG="--lead-only"
python3 -u "${SCRIPTS_DIR}/39b_cs_nominal_lookup.py" index \
    --out-dir "$OUT_DIR" --finemap-dir "${RESULTS_DIR}/finemap" \
    --cs-source "$CS_SOURCE" --min-pip "$MIN_PIP" $LEAD_FLAG --threads 16 \
    --modalities ${MODALITIES//:/ }
EOF
    DEP='-w done(csnom_index)'
    echo "  submitted csnom_index"
fi

bsub -J "csnom_x[1-${N}]" -q "$QUEUE" -n 2 -W "$WALLTIME" -M "$MEM" -R "rusage[mem=${MEM}]" \
     $DEP \
     -o "${LOG_DIR}/csnom_x.%J.%I.out" -e "${LOG_DIR}/csnom_x.%J.%I.err" \
     -env "$ENV_STR" <<'EOF'
#!/bin/bash
set -eo pipefail
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate tensorqtl
module load samtools
read ANC MOD < <(sed -n "${LSB_JOBINDEX}p" "$TASKS")
FORCE_FLAG=""; [ "$FORCE" = "1" ] && FORCE_FLAG="--force"
python3 -u "${SCRIPTS_DIR}/39b_cs_nominal_lookup.py" extract \
    --out-dir "$OUT_DIR" --results-dir "$RESULTS_DIR" --qtl-dir "$QTL_DIR" \
    --ancestry "$ANC" --modality "$MOD" --scope "$SCOPE" \
    --backend tabix $FORCE_FLAG
EOF
echo "  submitted csnom_x[1-${N}]"

bsub -J csnom_pkg -q "$QUEUE" -n 1 -W "$WALLTIME" -M 8G -R "rusage[mem=8G]" \
     -w "ended(csnom_x)" \
     -o "${LOG_DIR}/csnom_pkg.%J.out" -e "${LOG_DIR}/csnom_pkg.%J.err" \
     -env "$ENV_STR" <<'EOF'
#!/bin/bash
set -eo pipefail
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate tensorqtl
MISS_FLAG=""; [ "$ALLOW_MISSING" = "1" ] && MISS_FLAG="--allow-missing"
python3 -u "${SCRIPTS_DIR}/39b_cs_nominal_lookup.py" package \
    --out-dir "$OUT_DIR" --ancestries ${ANCESTRIES//:/ } \
    --modalities ${MODALITIES//:/ } $MISS_FLAG
EOF
echo "  submitted csnom_pkg (archive -> $(dirname "$OUT_DIR")/$(basename "$OUT_DIR")_<date>.tar.gz)"
