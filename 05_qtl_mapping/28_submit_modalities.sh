#!/bin/bash
# =============================================================================
# 28_submit_modalities.sh — Submit tensorQTL jobs for all ancestries × modalities
# =============================================================================
# Submits one LSF job per ancestry × modality. Each job runs
# 27_run_tensorqtl.sh with MODALITY set.
#
# Prerequisites:
#   - Scripts 23/24/25 completed (qtl_inputs has metadata, covariates, pgen)
#   - Script 26 completed (qtl_inputs has {ANC}_{modality}_harmonized.bed
#     and {ANC}_{modality}.phenotype_groups.txt)
#   - For MODALITIES including "combined": script 30 completed
#     ({ANC}_combined_harmonized.bed + {ANC}_combined.phenotype_groups.txt)
#
# Usage:
#   TEST=1 bash 28_submit_modalities.sh   # submit ONE pilot job (first combo)
#   bash 28_submit_modalities.sh          # submit all remaining combos
#
# Modes (combine freely):
#   GROUPED=0    — ungrouped: per-phenotype lead variants (transcript-level
#                  driver resolution: which isoform/junction/site carries
#                  each signal). Outputs {ANC}_{MOD}_ungrouped_cisqtl.*
#                  Default: GROUPED=1 (one lead variant per gene).
#   INDEPENDENT=1 — after map_cis, run cis.map_independent stepwise
#                  regression (PANTRY-style conditionally independent
#                  xQTLs). Outputs {ANC}_{MOD}[_ungrouped]_cisqtl_independent.*
#                  Default: 0.
#
# Examples (GTEx-conventions round, MAF 0.01 + MAC>=5 carrier floor baked
# into the qtl pgen by 23/24):
#   MAF_THRESHOLD=0.01 GROUPED=0 bash 28_submit_modalities.sh
#       # 14 ungrouped per-modality jobs (transcript-level driver layer)
#   MAF_THRESHOLD=0.01 INDEPENDENT=1 \
#       MODALITIES="expression alt_polyA alt_TSS intron_retention isoforms RNA_editing splicing stability" \
#       bash 28_submit_modalities.sh
#       # 16 grouped+stepwise jobs (PANTRY separate-modality layer)
#   MAF_THRESHOLD=0.01 INDEPENDENT=1 MODALITIES="combined" \
#       QUEUE=long WALLTIME=48:00 bash 28_submit_modalities.sh
#       # 2 combined cross-modality jobs (PANTRY cross-modality layer)
#
# Combos with existing results or a running/pending job of the same name
# are skipped, so it is safe to rerun this script at any time.
#
# Optional overrides (environment):
#   ANCESTRIES    — default "EAS EUR"
#   MODALITIES    — default all 7 non-expression modalities
#   QUEUE         — default medium
#   WALLTIME      — default 12:00 (combined runs are ~5x larger than the
#                   biggest per-modality run; use QUEUE=long WALLTIME=48:00)
#   CIS_WINDOW    — default 1000000
#   MAF_THRESHOLD — REQUIRED (no default): set explicitly, e.g. 0.01.
#                   (Deliberately required so the threshold is never
#                   silently inherited. Note the qtl pgen already carries
#                   a MAC>=5 carrier floor from 23/24, so MAF 0.01 no
#                   longer admits 2-4-carrier variants.)
#   COVARIATES_FILE — covariates TSV path; the literal placeholder {ANC} is
#                   replaced per ancestry by the wrapper (e.g.
#                   .../qtl_inputs/{ANC}_covariates_optimized.tsv).
#                   Default: {ANC}_covariates.tsv in qtl_inputs.
#
# NOTE: the -env string deliberately contains NO inner quotes — LSF
# preserves them literally in the variable values, which mangles paths.
# All values passed here are space-free (one ancestry per job).
# =============================================================================

set -eo pipefail

CONFIG="${CONFIG:-/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/config.yml}"
# Default SCRIPTS_DIR to this script's own directory (repo-clone layout);
# passed through to 27_run_tensorqtl.sh via bsub -env.
SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
OUTPUT_BASE="${OUTPUT_BASE:-/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY}"
LOG_DIR="/rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs"

ANCESTRIES="${ANCESTRIES:-EAS EUR}"
MODALITIES="${MODALITIES:-alt_polyA alt_TSS intron_retention isoforms RNA_editing splicing stability}"
QUEUE="${QUEUE:-medium}"
WALLTIME="${WALLTIME:-12:00}"
CIS_WINDOW="${CIS_WINDOW:-1000000}"
TEST="${TEST:-0}"
GROUPED="${GROUPED:-1}"
INDEPENDENT="${INDEPENDENT:-0}"
COVARIATES_FILE="${COVARIATES_FILE:-}"
QVALUE_METHOD="${QVALUE_METHOD:-storey}"
MAF_THRESHOLD="${MAF_THRESHOLD:?ERROR: set MAF_THRESHOLD explicitly (e.g. MAF_THRESHOLD=0.01) — it is never defaulted, to avoid silently running the wrong threshold}"

# ---- Mode-dependent suffixes ----
# Output base name: {ANC}_{MOD}{BASE_SUFFIX}_cisqtl{IND_SUFFIX}.parquet
# Job name:         eqtl_{ANC}_{MOD}{JOB_SUFFIX}
BASE_SUFFIX=""
MODE_DESC="grouped (one lead variant per gene)"
if [ "$GROUPED" != "1" ]; then
    NO_GROUPS=1
    BASE_SUFFIX="_ungrouped"
    MODE_DESC="UNGROUPED (per-phenotype lead variants)"
else
    NO_GROUPS=0
fi
IND_SUFFIX=""
JOB_SUFFIX="$BASE_SUFFIX"
if [ "$INDEPENDENT" = "1" ]; then
    IND_SUFFIX="_independent"
    JOB_SUFFIX="${JOB_SUFFIX}_ind"
    MODE_DESC="$MODE_DESC + stepwise (map_independent)"
fi

mkdir -p "$LOG_DIR"

echo "Submitting tensorQTL jobs:"
echo "  CONFIG:      $CONFIG"
echo "  SCRIPTS_DIR: $SCRIPTS_DIR"
echo "  OUTPUT_BASE: $OUTPUT_BASE"
echo "  Ancestries:  $ANCESTRIES"
echo "  Modalities:  $MODALITIES"
echo "  Mode:        $MODE_DESC"
echo "  MAF:         $MAF_THRESHOLD"
echo "  Queue: $QUEUE, walltime: $WALLTIME, TEST: $TEST"
echo ""

# Sanity checks before submitting anything
if [ ! -f "${SCRIPTS_DIR}/27_run_tensorqtl.sh" ]; then
    echo "ERROR: ${SCRIPTS_DIR}/27_run_tensorqtl.sh not found"
    exit 1
fi
if [ ! -f "${SCRIPTS_DIR}/27_run_tensorqtl.py" ]; then
    echo "ERROR: ${SCRIPTS_DIR}/27_run_tensorqtl.py not found"
    exit 1
fi

N=0
for ANC in $ANCESTRIES; do
    for MOD in $MODALITIES; do
        # Skip if harmonized BED is missing
        if [ ! -f "${OUTPUT_BASE}/qtl_inputs/${ANC}_${MOD}_harmonized.bed" ]; then
            if [ "$MOD" = "combined" ]; then
                echo "  SKIP ${ANC}/${MOD}: combined BED not found (run 30_combine_modalities.py first)"
            else
                echo "  SKIP ${ANC}/${MOD}: harmonized BED not found (run 26 first)"
            fi
            continue
        fi
        # Skip if results already exist for this mode (check the FINAL product)
        if [ -f "${OUTPUT_BASE}/qtl_results/${ANC}_${MOD}${BASE_SUFFIX}_cisqtl${IND_SUFFIX}.parquet" ]; then
            echo "  SKIP ${ANC}/${MOD}${JOB_SUFFIX}: results already exist"
            continue
        fi
        # Skip if a job with the same name is already running/pending
        if bjobs -J "eqtl_${ANC}_${MOD}${JOB_SUFFIX}" 2>/dev/null | grep -q "eqtl_${ANC}_${MOD}${JOB_SUFFIX}"; then
            echo "  SKIP ${ANC}/${MOD}${JOB_SUFFIX}: job already running/pending"
            continue
        fi
        ENV_STR="CONFIG=${CONFIG},SCRIPTS_DIR=${SCRIPTS_DIR},ANCESTRIES=${ANC},MODALITY=${MOD},NO_GROUPS=${NO_GROUPS},INDEPENDENT=${INDEPENDENT},MAF_THRESHOLD=${MAF_THRESHOLD},CIS_WINDOW=${CIS_WINDOW},QVALUE_METHOD=${QVALUE_METHOD}"
        if [ -n "$COVARIATES_FILE" ]; then
            ENV_STR="${ENV_STR},COVARIATES_FILE=${COVARIATES_FILE}"
        fi
        bsub -J "eqtl_${ANC}_${MOD}${JOB_SUFFIX}" -q "$QUEUE" -n 8 -W "$WALLTIME" \
             -o "${LOG_DIR}/eqtl_${ANC}_${MOD}${JOB_SUFFIX}.%J.out" \
             -e "${LOG_DIR}/eqtl_${ANC}_${MOD}${JOB_SUFFIX}.%J.err" \
             -env "$ENV_STR" \
             < "${SCRIPTS_DIR}/27_run_tensorqtl.sh"
        N=$((N + 1))
        if [ "$TEST" = "1" ]; then
            echo ""
            echo "TEST=1: submitted one pilot job (eqtl_${ANC}_${MOD}${JOB_SUFFIX})."
            echo "Check its log before submitting the rest:"
            echo "  tail ${LOG_DIR}/eqtl_${ANC}_${MOD}${JOB_SUFFIX}.<jobid>.out"
            if [ "$INDEPENDENT" = "1" ]; then
                echo "Look for 'Running cis.map_cis...' AND 'Running cis.map_independent' —"
                echo "then rerun without TEST=1."
            else
                echo "Look for 'Running cis.map_cis...' — then rerun without TEST=1."
            fi
            exit 0
        fi
    done
done

echo ""
echo "Submitted $N jobs ($MODE_DESC). Monitor with: bjobs"
echo "Logs: ${LOG_DIR}/eqtl_<ANC>_<MOD>${JOB_SUFFIX}.<jobid>.out"
