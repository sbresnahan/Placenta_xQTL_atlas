#!/bin/bash
# =============================================================================
# 51_pooled_hcp.sh — pooled multi-ancestry HCP estimation (one modality)
# =============================================================================
# Objective 2.1 tier-1 covariates: HCP latent factors re-estimated on the
# POOLED (within-ancestry z-scored) modality BED, using the module-05
# estimator (05_qtl_mapping/hcp_from_matrix.R — same Rhcpp settings, pooled
# Picard QC metrics as the known-covariate matrix, cohort dummies for the
# VST-scale modalities). k defaults to the max HCP count across the
# per-ancestry optimized covariate tables ({QTL_DIR}/{ANC}_covariates_{MOD}.tsv),
# falling back to the shared {ANC}_covariates.tsv, then to 20.
#
# After estimation, 50_build_gxe_inputs.py --finalize-covariates appends the
# HCP rows to the pooled base covariates (HCPs lowest pruning priority) and
# writes inputs/pooled_covariates_{MOD}.tsv.
#
# Worker contract: runs ONE modality. MODALITY env var, or $LSB_JOBINDEX
# (1-based) indexing the MODALITIES list when submitted as an array by
# 53_submit_gxe.sh.
#
# Usage:
#   MODALITY=expression bash 51_pooled_hcp.sh     # standalone
#   (53_submit_gxe.sh submits the array form)
#
# Optional env: CONFIG, SCRIPTS_DIR, OUTPUT_BASE, QTL_DIR, RESULTS_DIR,
#   GXE_DIR, MODALITIES, HCP_K (override), COHORT_DUMMY_MODALITIES
#   (default "expression isoform_expression"), RSCRIPT, R_PACKAGE_LIB,
#   MAX_PHENOTYPES, CONDA_EXE, CONDA_ENV.
# =============================================================================
set -eo pipefail

CONFIG="${CONFIG:?ERROR: CONFIG env var required}"
SCRIPTS_DIR="${SCRIPTS_DIR:?ERROR: SCRIPTS_DIR env var required}"
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

# LSF workers must initialize their own software environment. Do not rely on
# the interactive submit shell having tensorqtl activated.
source /etc/profile.d/modules.sh
[ -x "$CONDA_EXE" ] || { echo "ERROR: conda executable not found/executable: $CONDA_EXE"; exit 1; }
eval "$("$CONDA_EXE" shell.bash hook)"
conda activate "$CONDA_ENV"
command -v python3 >/dev/null || { echo "ERROR: python3 not found after conda activation"; exit 1; }

# ---- resolve modality (env var or LSF array index) ----
if [ -z "${MODALITY:-}" ]; then
    MODALITY="$(echo $MODALITIES | awk -v i="${LSB_JOBINDEX:?set MODALITY or run under LSF}" '{print $i}')"
fi
echo "=== 51_pooled_hcp.sh: ${MODALITY} ==="

BED_GZ="${INPUTS}/pooled_${MODALITY}.bed.gz"
META="${INPUTS}/pooled_metadata.tsv"
HCP_OUT="${INPUTS}/pooled_hcp_${MODALITY}.tsv"
FINAL_COV="${INPUTS}/pooled_covariates_${MODALITY}.tsv"
for f in "$BED_GZ" "$META" "${INPUTS}/pooled_covariates_base.tsv"; do
    [ -f "$f" ] || { echo "ERROR: missing $f — run 50_build_gxe_inputs.py first"; exit 1; }
done
if [ -f "$FINAL_COV" ] && [ "${FORCE:-0}" != "1" ]; then
    echo "  ${FINAL_COV} exists — skipping (FORCE=1 to redo)"; exit 0
fi

# ---- k: max HCP rows across per-ancestry optimized covariates ----
if [ -z "$HCP_K" ]; then
    K=0
    for ANC in ${ANCESTRIES:-EAS EUR}; do
        for COV in "${QTL_DIR}/${ANC}_covariates_${MODALITY}.tsv" "${QTL_DIR}/${ANC}_covariates.tsv"; do
            if [ -f "$COV" ]; then
                N=$(cut -f1 "$COV" | grep -c '^HCP_' || true)
                [ "$N" -gt "$K" ] && K=$N
                break
            fi
        done
    done
    [ "$K" -eq 0 ] && K=20
    HCP_K=$K
fi
echo "  pooled HCP k = ${HCP_K}"

# ---- HCP estimation (hcp_from_matrix.R needs a plain-text BED) ----
if [ "$HCP_K" -ge 1 ]; then
    SCRATCH_BED="${TMPDIR:-/tmp}/pooled_${MODALITY}.$$.bed"
    zcat "$BED_GZ" > "$SCRATCH_BED"
    EXTRA=()
    case " $COHORT_DUMMY_MODALITIES " in
        *" ${MODALITY} "*) EXTRA+=(--cohort-dummies);;
    esac
    # Module-07 contract: the required package library must be set from
    # inside R. Do not rely on R_LIBS_USER/R_LIBS_SITE exported by the shell
    # or container wrapper. Run a temporary copy of the Module-05 script with
    # an R preamble that keeps the required library first even if the sourced
    # script later calls .libPaths() itself.
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
    # k=0: header-only HCP file (finalize appends nothing)
    (printf 'covariate'; zcat "$BED_GZ" | head -1 | cut -f5- | tr '\t' '\n' | sed 's/^/\t/'; echo) \
        | tr -d '\n' | sed 's/\t/\t/g' > "$HCP_OUT"
    echo "" >> "$HCP_OUT"
fi

# ---- finalize per-modality covariates ----
python3 "${SCRIPTS_DIR}/50_build_gxe_inputs.py" \
    --qtl-dir "$QTL_DIR" --results-dir "$RESULTS_DIR" \
    --gxe-dir "$GXE_DIR" \
    --finalize-covariates --modality "$MODALITY" --hcp-file "$HCP_OUT"
echo "=== done: ${MODALITY} ==="
