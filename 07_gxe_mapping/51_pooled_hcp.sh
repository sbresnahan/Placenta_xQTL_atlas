#!/bin/bash
# =============================================================================
# 51_pooled_hcp.sh — ancestry-specific HCP estimation (one ancestry x modality)
# =============================================================================
# Historical filename retained for compatibility. Primary Module-07 discovery
# is ancestry-stratified; this worker estimates HCP factors separately within
# ANCESTRY for MODALITY, then finalizes inputs/{ANC}_covariates_{MOD}.tsv.
#
# Required env: CONFIG, SCRIPTS_DIR, ANCESTRY. MODALITY may be provided
# directly or resolved from LSB_JOBINDEX against MODALITIES.
# =============================================================================
set -eo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"
ANCESTRY="${ANCESTRY:?ERROR: ANCESTRY env var required}"
REPO_ROOT="$(cd "${SCRIPTS_DIR}/.." && pwd)"
OUTPUT_BASE="${OUTPUT_BASE:-$(dirname "$CONFIG")}"
QTL_DIR="${QTL_DIR:-${OUTPUT_BASE}/qtl_inputs}"
RESULTS_DIR="${RESULTS_DIR:-${OUTPUT_BASE}/qtl_results}"
GXE_DIR="${GXE_DIR:-${RESULTS_DIR}/gxe}"
INPUTS="${GXE_DIR}/inputs"
MODALITIES="${MODALITIES:-expression isoforms isoform_expression splicing intron_retention alt_TSS alt_polyA RNA_editing stability}"
COHORT_DUMMY_MODALITIES="${COHORT_DUMMY_MODALITIES:-expression isoform_expression}"
HCP_K="${HCP_K:-}"
MAX_PHENOTYPES="${MAX_PHENOTYPES:-40000}"
RSCRIPT="${RSCRIPT:-${REPO_ROOT}/bin/Rscript_sif}"
R_PACKAGE_LIB="${R_PACKAGE_LIB:-/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1}"
QC_METRICS="${QC_METRICS:-${OUTPUT_BASE}/hcp/all_qc_metrics.tsv}"
MODULE05="${MODULE05:-${REPO_ROOT}/05_qtl_mapping}"
CONDA_EXE="${CONDA_EXE:-/risapps/rhel8/miniforge3/24.5.0-0/bin/conda}"
CONDA_ENV="${CONDA_ENV:-tensorqtl}"

source /etc/profile.d/modules.sh
[ -x "$CONDA_EXE" ] || { echo "ERROR: conda executable not found/executable: $CONDA_EXE"; exit 1; }
eval "$("$CONDA_EXE" shell.bash hook)"
conda activate "$CONDA_ENV"
command -v python3 >/dev/null || { echo "ERROR: python3 not found after conda activation"; exit 1; }

if [ -z "${MODALITY:-}" ]; then
    MODALITY="$(echo $MODALITIES | awk -v i="${LSB_JOBINDEX:?set MODALITY or run under LSF}" '{print $i}')"
fi
echo "=== 51_pooled_hcp.sh: ${ANCESTRY} ${MODALITY} ==="

BED_GZ="${INPUTS}/${ANCESTRY}_${MODALITY}.bed.gz"
META="${INPUTS}/${ANCESTRY}_metadata.tsv"
HCP_OUT="${INPUTS}/${ANCESTRY}_hcp_${MODALITY}.tsv"
FINAL_COV="${INPUTS}/${ANCESTRY}_covariates_${MODALITY}.tsv"
BASE_COV="${INPUTS}/${ANCESTRY}_covariates_base.tsv"
for f in "$BED_GZ" "$META" "$BASE_COV"; do
    [ -f "$f" ] || { echo "ERROR: missing $f — run Stage 0 first"; exit 1; }
done
if [ -f "$FINAL_COV" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "  ${FINAL_COV} exists — skipping (FORCE=1 to redo)"; exit 0
fi

# Match Module-05's ancestry/modality optimized HCP count when available.
if [ -z "$HCP_K" ]; then
    HCP_K=0
    for COV in "${QTL_DIR}/${ANCESTRY}_covariates_${MODALITY}.tsv" "${QTL_DIR}/${ANCESTRY}_covariates.tsv"; do
        if [ -f "$COV" ]; then
            HCP_K=$(cut -f1 "$COV" | grep -c '^HCP_' || true)
            break
        fi
    done
    [ "$HCP_K" -eq 0 ] && HCP_K=20
fi
echo "  ${ANCESTRY} HCP k = ${HCP_K}"

if [ "$HCP_K" -ge 1 ]; then
    SCRATCH_BED="${TMPDIR:-/tmp}/${ANCESTRY}_${MODALITY}.$$.bed"
    zcat "$BED_GZ" > "$SCRATCH_BED"
    EXTRA=()
    case " $COHORT_DUMMY_MODALITIES " in
        *" ${MODALITY} "*) EXTRA+=(--cohort-dummies);;
    esac

    HCP_R_SRC="${MODULE05}/hcp_from_matrix.R"
    [ -f "$HCP_R_SRC" ] || { echo "ERROR: missing $HCP_R_SRC"; exit 1; }
    SCRATCH_R="${TMPDIR:-/tmp}/hcp_from_matrix.module07.$$.R"
    R_PACKAGE_LIB_R=${R_PACKAGE_LIB//\\/\\\\}
    R_PACKAGE_LIB_R=${R_PACKAGE_LIB_R//\"/\\\"}
    {
        printf '%s\n' \
            ".module07_r_lib <- \"${R_PACKAGE_LIB_R}\"" \
            'if (!dir.exists(.module07_r_lib)) stop("Required Module-07 R package library does not exist: ", .module07_r_lib)' \
            '.module07_base_libPaths <- base::.libPaths' \
            '.module07_base_libPaths(c(.module07_r_lib, .module07_base_libPaths()))' \
            '.libPaths <- function(new) {' \
            '    if (missing(new)) return(.module07_base_libPaths())' \
            '    .module07_base_libPaths(c(.module07_r_lib, new))' \
            '}' \
            'if (normalizePath(.libPaths()[1], mustWork=TRUE) != normalizePath(.module07_r_lib, mustWork=TRUE)) stop("Failed to prepend required Module-07 R package library")'
        cat "$HCP_R_SRC"
    } > "$SCRATCH_R"

    cleanup_hcp_scratch() { rm -f "$SCRATCH_BED" "$SCRATCH_R"; }
    trap cleanup_hcp_scratch EXIT

    "$RSCRIPT" "$SCRATCH_R" \
        --bed "$SCRATCH_BED" \
        --qc-metrics "$QC_METRICS" \
        --metadata "$META" \
        --k "$HCP_K" \
        --max-phenotypes "$MAX_PHENOTYPES" \
        "${EXTRA[@]}" \
        --output "$HCP_OUT"

    cleanup_hcp_scratch
    trap - EXIT
else
    python3 - "$BED_GZ" "$HCP_OUT" <<'PY'
import gzip, sys
bed, out = sys.argv[1:]
with gzip.open(bed, 'rt') as fh:
    samples = fh.readline().rstrip('\n').split('\t')[4:]
with open(out, 'w') as fh:
    fh.write('covariate\t' + '\t'.join(samples) + '\n')
PY
fi

python3 "${SCRIPTS_DIR}/50_build_gxe_inputs.py" \
    --qtl-dir "$QTL_DIR" --results-dir "$RESULTS_DIR" --gxe-dir "$GXE_DIR" \
    --finalize-covariates --ancestry "$ANCESTRY" --modality "$MODALITY" \
    --hcp-file "$HCP_OUT"
echo "=== done: ${ANCESTRY} ${MODALITY} ==="
