#!/bin/bash
# =============================================================================
# 40b_submit_1kg_ld_reference.sh — Build ancestry/chromosome 1KG LD references
# =============================================================================
# One-time preprocessing for stages 42, 45, and 48. The full 70M-variant 1KG PGEN
# is subset ONCE by superpopulation and chromosome into compact references:
#
#   $KG_LD_REF_DIR/EAS/chr1.{pgen,pvar,psam} ... chr22
#   $KG_LD_REF_DIR/EUR/chr1.{pgen,pvar,psam} ... chr22
#
# Runtime coloc/colocBoost jobs consume these compact references directly and
# never run --keep against the full 1KG PGEN.
#
# This submitter creates/refreshes private ancestry keep files with the
# stdlib-only 40b_prepare_1kg_reference.py helper, then submits one LSF element per
# ancestry x autosome. Re-running is safe: completed chromosomes are skipped.
#
# Usage:
#   TEST=1 bash 40b_submit_1kg_ld_reference.sh   # one EAS chr1 pilot
#   bash 40b_submit_1kg_ld_reference.sh          # all EAS/EUR autosomes
#
# Optional overrides:
#   ANCESTRIES       — default "EAS EUR"
#   KG_PGEN           — full 1KG PGEN prefix
#   KG_SAMPLE_MAP     — sample -> superpopulation map
#   KG_LD_REF_DIR     — default $COLOC_DIR/ld_reference/1kg
#   QUEUE             — default medium
#   WALLTIME          — default 04:00
#   MEM_GB            — default 16
#   PLINK_MEMORY_MB   — default 12000
#   MAXCONC           — max simultaneous full-PGEN readers, default 4
#   FORCE_REF=1       — rebuild existing chromosome references
#   TEST=1            — submit only array index 1
# =============================================================================
set -eo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"
OUTPUT_BASE="${OUTPUT_BASE:-$(dirname "$CONFIG")}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
GWAS_DIR="${GWAS_DIR:-${OUTPUT_BASE}/gwas}"
COLOC_DIR="${COLOC_DIR:-${RESULTS_DIR}/coloc}"
LOG_DIR="${LOG_DIR:-${OUTPUT_BASE}/logs}"
KG_PGEN="${KG_PGEN:-/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/1kGP/1kGP_hg38}"
KG_SAMPLE_MAP="${KG_SAMPLE_MAP:-${GWAS_DIR}/1kg_sample_superpop.tsv}"
KG_LD_REF_DIR="${KG_LD_REF_DIR:-${COLOC_DIR}/ld_reference/1kg}"
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
QUEUE="${QUEUE:-medium}"
WALLTIME="${WALLTIME:-04:00}"
MEM_GB="${MEM_GB:-16}"
PLINK_MEMORY_MB="${PLINK_MEMORY_MB:-12000}"
MAXCONC="${MAXCONC:-4}"
FORCE_REF="${FORCE_REF:-0}"
TEST="${TEST:-0}"

mkdir -p "$KG_LD_REF_DIR" "$LOG_DIR"

# Build private one-time keep files under the reference directory. Runtime
# coloc/colocBoost/FUSION jobs never consume these files directly.
python3 "${SCRIPTS_DIR}/40b_prepare_1kg_reference.py" \
    --sample-map "$KG_SAMPLE_MAP" \
    --kg-pgen "$KG_PGEN" \
    --ancestries $ANCESTRIES \
    --out-dir "$KG_LD_REF_DIR"

# Pass ancestries as a colon-separated scalar through LSF -env.
ANCESTRIES_CSV="$(printf '%s\n' $ANCESTRIES | paste -sd: -)"
N_ANC=$(echo "$ANCESTRIES" | wc -w)
N_JOBS=$((N_ANC * 22))
ARRAY_RANGE="1-${N_JOBS}"
ARRAY_LIMIT="%${MAXCONC}"
if [ "$TEST" = "1" ]; then
    ARRAY_RANGE="1"
    ARRAY_LIMIT=""
fi

# LSF array concurrency syntax places the limiter *after* the closing bracket:
#   name[1-44]%4
# not name[1-44%4], which LSF rejects as "Bad job name".
ARRAY_DISPLAY="${ARRAY_RANGE}${ARRAY_LIMIT}"

echo "=== 40b_submit_1kg_ld_reference.sh ==="
echo "  ANCESTRIES:      $ANCESTRIES"
echo "  KG_PGEN:         $KG_PGEN"
echo "  KG_LD_REF_DIR:   $KG_LD_REF_DIR"
echo "  ARRAY:           $ARRAY_DISPLAY (${N_JOBS} total; max ${MAXCONC} concurrent)"
echo "  QUEUE/WALLTIME:  $QUEUE / $WALLTIME"
echo "  MEMORY:          ${MEM_GB}G LSF; ${PLINK_MEMORY_MB} MiB plink2 cap"

if bjobs -J "kg_ldref" 2>/dev/null | grep -q "kg_ldref"; then
    echo "ERROR: kg_ldref array is already running/pending" >&2
    exit 1
fi

ENV_STR="SCRIPTS_DIR=${SCRIPTS_DIR},KG_PGEN=${KG_PGEN},KG_LD_REF_DIR=${KG_LD_REF_DIR},ANCESTRIES_CSV=${ANCESTRIES_CSV},PLINK_MEMORY_MB=${PLINK_MEMORY_MB},FORCE_REF=${FORCE_REF}"
bsub -J "kg_ldref[${ARRAY_RANGE}]${ARRAY_LIMIT}" -q "$QUEUE" -n 1 -W "$WALLTIME" \
     -M "${MEM_GB}G" -R "rusage[mem=${MEM_GB}G]" \
     -o "${LOG_DIR}/kg_ldref_%J_%I.out" \
     -e "${LOG_DIR}/kg_ldref_%J_%I.err" \
     -env "$ENV_STR" \
     < "${SCRIPTS_DIR}/40c_run_1kg_ld_reference.sh"

if [ "$TEST" = "1" ]; then
    echo "TEST=1: submitted only $(echo "$ANCESTRIES" | awk '{print $1}') chr1 (array index 1)."
    echo "After it finishes, inspect: ${KG_LD_REF_DIR}/$(echo "$ANCESTRIES" | awk '{print $1}')/chr1.done"
else
    echo "Submitted ${N_JOBS} ancestry/chromosome reference jobs."
    echo "Rerun 42_submit_coloc.sh only after every *.done marker exists."
fi
