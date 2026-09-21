#!/bin/bash
# =============================================================================
# cleanup_original_run.sh — remove/move original-run (MAF 0.01) QTL files
# =============================================================================
# Prepares qtl_inputs/ and qtl_results/ for the PANTRY-style MAF 0.05 rerun.
#
#   DRY-RUN by default — prints every action without touching anything:
#       bash cleanup_original_run.sh
#   Execute for real:
#       bash cleanup_original_run.sh --yes
#
# Actions:
#   DELETE (superseded or regenerable):
#     - qtl_results/{ANC}_{mod}_cisqtl.parquet + _top.tsv   (first-run grouped
#       per-modality results, 7 modalities x 2 ancestries = 28 files; replaced
#       by the grouped+stepwise rerun)
#     - qtl_results/{ANC}_expression_cisqtl_top.tsv         (buggy head(1000)
#       top tables; the parquets are MOVED, not deleted)
#     - qtl_inputs/{ANC}_{mod}.bed.gz + .tbi and
#       qtl_inputs/{ANC}_expression.bed.gz + .tbi           (regenerable
#       bgzip/tabix intermediates; 27_run_tensorqtl.sh rebuilds them from the
#       harmonized BEDs)
#   MOVE to qtl_results/archive_maf01/ (stale MAF 0.01 results: the
#   expression parquets would otherwise be overwritten by the 0.05 rerun,
#   and the old ungrouped parquets would block the 0.05 ungrouped rerun
#   via the driver's skip-logic. Delete the archive manually once the
#   0.05 results are validated):
#     - qtl_results/{ANC}_expression_cisqtl.parquet
#     - qtl_results/{ANC}_{mod}_ungrouped_cisqtl.parquet + _top.tsv
#   KEEP (needed for the rerun/report — never touched):
#     - *_harmonized.bed, *.phenotype_groups.txt, *.phenotype_modalities.tsv
#     - genotypes (*.pgen/*.pvar/*.psam), covariates, metadata,
#       *_replicate_merges.tsv
#
# Optional env:
#   OUTPUT_BASE — default /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY
#   ANCESTRIES  — default "EAS EUR"
#   MODALITIES  — default all 7 non-expression modalities
# =============================================================================

set -eo pipefail

OUTPUT_BASE="${OUTPUT_BASE:-/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY}"
ANCESTRIES="${ANCESTRIES:-EAS EUR}"
MODALITIES="${MODALITIES:-alt_polyA alt_TSS intron_retention isoforms RNA_editing splicing stability}"
QTL_DIR="${OUTPUT_BASE}/qtl_inputs"
RESULTS_DIR="${OUTPUT_BASE}/qtl_results"
ARCHIVE_DIR="${RESULTS_DIR}/archive_maf01"

EXECUTE=0
if [ "${1:-}" = "--yes" ]; then
    EXECUTE=1
fi

if [ "$EXECUTE" = "1" ]; then
    echo "MODE: EXECUTE (--yes given) — files WILL be deleted/moved"
else
    echo "MODE: DRY-RUN — no files will be touched (rerun with --yes to execute)"
fi
echo "  OUTPUT_BASE: $OUTPUT_BASE"
echo "  Archive dir: $ARCHIVE_DIR"
echo ""

N_DEL=0
N_DEL_MISSING=0
N_MOVE=0
N_MOVE_MISSING=0

delete_file() {
    if [ -f "$1" ]; then
        echo "  DELETE  $1"
        if [ "$EXECUTE" = "1" ]; then rm -f "$1"; fi
        N_DEL=$((N_DEL + 1))
    else
        N_DEL_MISSING=$((N_DEL_MISSING + 1))
    fi
}

move_file() {
    if [ -f "$1" ]; then
        echo "  MOVE    $1  ->  $ARCHIVE_DIR/"
        if [ "$EXECUTE" = "1" ]; then mv "$1" "$ARCHIVE_DIR/"; fi
        N_MOVE=$((N_MOVE + 1))
    else
        N_MOVE_MISSING=$((N_MOVE_MISSING + 1))
    fi
}

if [ "$EXECUTE" = "1" ]; then
    mkdir -p "$ARCHIVE_DIR"
fi

for ANC in $ANCESTRIES; do
    echo "------------------------------------------------------------"
    echo "Ancestry: $ANC"
    echo "------------------------------------------------------------"

    # ---- DELETE: first-run grouped per-modality results ----
    for MOD in $MODALITIES; do
        delete_file "${RESULTS_DIR}/${ANC}_${MOD}_cisqtl.parquet"
        delete_file "${RESULTS_DIR}/${ANC}_${MOD}_cisqtl_top.tsv"
    done

    # ---- DELETE: buggy expression top tables (head(1000) bug) ----
    delete_file "${RESULTS_DIR}/${ANC}_expression_cisqtl_top.tsv"

    # ---- DELETE: regenerable bgzip/tabix intermediates ----
    delete_file "${QTL_DIR}/${ANC}_expression.bed.gz"
    delete_file "${QTL_DIR}/${ANC}_expression.bed.gz.tbi"
    for MOD in $MODALITIES; do
        delete_file "${QTL_DIR}/${ANC}_${MOD}.bed.gz"
        delete_file "${QTL_DIR}/${ANC}_${MOD}.bed.gz.tbi"
    done

    # ---- MOVE: stale MAF 0.01 results that block the 0.05 rerun ----
    move_file "${RESULTS_DIR}/${ANC}_expression_cisqtl.parquet"
    for MOD in $MODALITIES; do
        move_file "${RESULTS_DIR}/${ANC}_${MOD}_ungrouped_cisqtl.parquet"
        move_file "${RESULTS_DIR}/${ANC}_${MOD}_ungrouped_cisqtl_top.tsv"
    done
done

echo ""
echo "============================================================"
echo "Summary"
echo "============================================================"
echo "  To delete: $N_DEL files ($N_DEL_MISSING already absent)"
echo "  To move:   $N_MOVE files ($N_MOVE_MISSING already absent)"
echo ""
echo "  KEPT (never touched):"
echo "    *_harmonized.bed, *.phenotype_groups.txt, *.phenotype_modalities.tsv"
echo "    genotypes (*.pgen/*.pvar/*.psam), covariates, metadata,"
echo "    *_replicate_merges.tsv"
echo ""
if [ "$EXECUTE" = "0" ]; then
    echo "DRY-RUN only. Inspect the list above, then execute with:"
    echo "  bash cleanup_original_run.sh --yes"
else
    echo "Done. Stale MAF 0.01 results are in: $ARCHIVE_DIR"
    echo "Delete the archive manually once the 0.05 results are validated:"
    echo "  rm -rf $ARCHIVE_DIR"
fi
