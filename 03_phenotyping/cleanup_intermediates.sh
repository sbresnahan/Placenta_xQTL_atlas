#!/bin/bash
# =============================================================================
# cleanup_intermediates.sh — Delete per-cohort intermediate files that are no
# longer needed after the phenotyping pipeline (stages 01–16) completes.
#
# Keeps intermediates still required by downstream cross-cohort steps:
#   - intermediate/star_out/*.Aligned.sortedByCoord.out.bam  (Picard QC / HCP)
#   - intermediate/expression/<sample>/quant.sf               (deconv + HCP bias)
#   - intermediate/splicing/leafcutter_perind_numers.counts.gz (harmonize splicing)
#   - intermediate/intron_retention/retained_intron_psi.tsv.gz (harmonize IR)
#   - output/unnorm/<modality>.bed                            (cross-cohort pooling)
#   - output/  (all final bgzipped+tabix BEDs + phenotype_groups)
#
# Usage:
#   bash cleanup_intermediates.sh <DIRECTORY> <COHORT> [--dry-run]
#
# Example:
#   bash cleanup_intermediates.sh \
#       /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY cohort1
#
#   bash cleanup_intermediates.sh \
#       /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY cohort1 --dry-run
#
# Arguments:
#   DIRECTORY — output_base from config.yml (parent of cohort directories)
#   COHORT    — cohort name (subdirectory under DIRECTORY)
#   --dry-run — list what would be deleted with sizes, but delete nothing
# =============================================================================

set -euo pipefail

# ---- Args ----
if [ $# -lt 2 ]; then
    echo "Usage: bash $0 <DIRECTORY> <COHORT> [--dry-run]" >&2
    echo "  DIRECTORY — output_base (parent of cohort dirs)" >&2
    echo "  COHORT    — cohort name" >&2
    echo "  --dry-run — preview only, no deletion" >&2
    exit 1
fi

DIRECTORY="$1"
COHORT="$2"
DRY_RUN=false
if [ "${3:-}" = "--dry-run" ]; then
    DRY_RUN=true
fi

COHORT_DIR="${DIRECTORY}/${COHORT}"
INTERM_DIR="${COHORT_DIR}/intermediate"

# ---- Validate ----
if [ ! -d "$COHORT_DIR" ]; then
    echo "ERROR: cohort directory not found: $COHORT_DIR" >&2
    exit 1
fi
if [ ! -d "$INTERM_DIR" ]; then
    echo "ERROR: intermediate directory not found: $INTERM_DIR" >&2
    exit 1
fi

# ---- Helpers ----
total_freed=0

delete_path() {
    # Delete a file or directory if it exists. Reports size. Respects DRY_RUN.
    local path="$1"
    local label="$2"

    if [ ! -e "$path" ]; then
        return 0
    fi

    local size
    if [ -d "$path" ]; then
        size=$(du -sh "$path" 2>/dev/null | cut -f1)
    else
        size=$(du -sh "$path" 2>/dev/null | cut -f1)
    fi

    if $DRY_RUN; then
        echo "  [DRY-RUN] would delete: $path ($size) — $label"
    else
        rm -rf "$path"
        echo "  deleted: $path ($size) — $label"
    fi
}

delete_glob() {
    # Delete files matching a glob pattern within a directory.
    # Args: dir, pattern, label
    local dir="$1"
    local pattern="$2"
    local label="$3"

    if [ ! -d "$dir" ]; then
        return 0
    fi

    local count
    count=$(find "$dir" -maxdepth 1 -name "$pattern" 2>/dev/null | wc -l)
    if [ "$count" -eq 0 ]; then
        return 0
    fi

    local size
    size=$(find "$dir" -maxdepth 1 -name "$pattern" -print0 2>/dev/null \
           | xargs -0 du -ch 2>/dev/null | tail -1 | cut -f1)

    if $DRY_RUN; then
        echo "  [DRY-RUN] would delete: $dir/$pattern ($count files, $size) — $label"
    else
        find "$dir" -maxdepth 1 -name "$pattern" -delete 2>/dev/null || true
        echo "  deleted: $dir/$pattern ($count files, $size) — $label"
    fi
}

# ---- Banner ----
echo "================================================================"
echo "cleanup_intermediates.sh"
echo "  Cohort dir: $COHORT_DIR"
echo "  Dry run:    $DRY_RUN"
echo "================================================================"
echo ""

# ---- SAFE TO DELETE (consumed only within per-cohort pipeline) ----
echo "--- Deleting intermediates safe to remove ---"
echo ""

# 1. Shrunk BAMs (intermediate/bam/) — consumed by RegTools, featureCounts, MAJIQ
delete_path "${INTERM_DIR}/bam" "shrunk BAMs (RegTools/featureCounts/MAJIQ done)"

# 2. Per-sample splice junctions (intermediate/splicing/*.junc)
#    KEEP leafcutter_perind_numers.counts.gz (needed by harmonize)
delete_glob "${INTERM_DIR}/splicing" "*.junc" "per-sample junctions (leafCutter clustering done)"
delete_glob "${INTERM_DIR}/splicing" "juncfiles.txt" "juncfile list"

# 3. Per-sample featureCounts (intermediate/stability/) — consumed by step 15
delete_path "${INTERM_DIR}/stability" "per-sample exon/intron counts (stability aggregation done)"

# 4. Per-sample RNA editing levels (intermediate/RNA_editing/)
#    The whole directory: edit_levels/, edit_site_matrix.tsv, phenotype_matrix.tsv
#    All consumed within step 14. Pooling reads output/unnorm/RNA_editing.bed instead.
delete_path "${INTERM_DIR}/RNA_editing" "RNA editing intermediates (aggregation done; unnorm BED kept)"

# 5. Per-sample alt TSS/polyA Salmon quant (intermediate/alt_TSS_polyA/)
#    All 6 txrevise indices x all samples. Consumed by step 11.
delete_path "${INTERM_DIR}/alt_TSS_polyA" "alt TSS/polyA Salmon quant (aggregation done)"

# 6. QU-corrected quant.sf (intermediate/expression_qu/) — consumed within step 10
delete_path "${INTERM_DIR}/expression_qu" "QU-corrected quant.sf (isoform aggregation done)"

# 7. Salmon bootstraps (intermediate/expression/<sample>/aux_info/)
#    Consumed by edgeR::catchSalmon in step 10. KEEP quant.sf (deconv + HCP).
echo "  Salmon bootstraps (aux_info/ per sample):"
if [ -d "${INTERM_DIR}/expression" ]; then
    bootstrap_count=0
    for sample_dir in "${INTERM_DIR}/expression"/*/; do
        aux="${sample_dir}aux_info"
        if [ -d "$aux" ]; then
            bootstrap_count=$((bootstrap_count + 1))
            if $DRY_RUN; then
                size=$(du -sh "$aux" 2>/dev/null | cut -f1)
                echo "    [DRY-RUN] would delete: $aux ($size)"
            else
                rm -rf "$aux"
            fi
        fi
    done
    if [ "$bootstrap_count" -gt 0 ]; then
        if ! $DRY_RUN; then
            echo "    deleted aux_info/ from $bootstrap_count sample dirs"
        fi
    else
        echo "    (none found)"
    fi
else
    echo "    (expression dir not found)"
fi
echo ""

# 8. MAJIQ IR intermediates (intermediate/intron_retention/)
#    DELETE: build/, psicov/, majiq.psi.tsv, experiments.tsv
#    KEEP:   retained_intron_psi.tsv.gz (needed by harmonize)
delete_path "${INTERM_DIR}/intron_retention/build" "MAJIQ build intermediates"
delete_path "${INTERM_DIR}/intron_retention/psicov" "MAJIQ psi-coverage batches"
delete_path "${INTERM_DIR}/intron_retention/majiq.psi.tsv" "MAJIQ psi TSV"
delete_path "${INTERM_DIR}/intron_retention/experiments.tsv" "MAJIQ experiments list"

# 9. STAR log/tab files (intermediate/star_out/)
#    KEEP: *.Aligned.sortedByCoord.out.bam (Picard QC needs full BAM with SEQ/QUAL)
delete_glob "${INTERM_DIR}/star_out" "*.Log.final.out" "STAR final logs"
delete_glob "${INTERM_DIR}/star_out" "*.Log.out" "STAR logs"
delete_glob "${INTERM_DIR}/star_out" "*.Log.progress.out" "STAR progress logs"
delete_glob "${INTERM_DIR}/star_out" "*.SJ.out.tab" "STAR junction tabs"
delete_glob "${INTERM_DIR}/star_out" "*_STARpass1*" "STAR pass1 logs"

# 10. LSF logs (cohort_dir/logs/) — not consumed by anything
delete_path "${COHORT_DIR}/logs" "LSF job logs"

# ---- KEPT (printed for clarity) ----
echo ""
echo "--- Kept (required by downstream cross-cohort steps) ---"
echo ""

keep_check() {
    if [ -e "$1" ]; then
        local size
        size=$(du -sh "$1" 2>/dev/null | cut -f1)
        echo "  KEEP: $1 ($size) — $2"
    fi
}

keep_check "${INTERM_DIR}/star_out" \
    "full STAR BAMs (Picard QC / HCP — delete after 19_hcp_factors.sh)"
keep_check "${INTERM_DIR}/expression" \
    "Salmon quant.sf (deconvolution + HCP bias — delete after both done)"
keep_check "${INTERM_DIR}/splicing/leafcutter_perind_numers.counts.gz" \
    "leafCutter counts (harmonize splicing)"
keep_check "${INTERM_DIR}/intron_retention/retained_intron_psi.tsv.gz" \
    "MAJIQ IR PSI (harmonize intron retention)"
keep_check "${COHORT_DIR}/output/unnorm" \
    "unnorm BEDs (cross-cohort pooling — all modalities)"
keep_check "${COHORT_DIR}/output" \
    "final bgzipped+tabix BEDs + phenotype_groups"

# ---- Summary ----
echo ""
echo "================================================================"
if $DRY_RUN; then
    echo "DRY RUN complete — nothing was deleted."
    echo "Re-run without --dry-run to actually delete."
else
    echo "Cleanup complete for cohort: $COHORT"
fi
echo "================================================================"
