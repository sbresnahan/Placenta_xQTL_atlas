#!/bin/bash
# =============================================================================
# 16_index_outputs.sh — Final check: ensure all .bed.gz outputs are tabix-indexed
# =============================================================================
# The aggregation scripts (10-15) already bgzip + tabix each output. This
# script is a safety net: it verifies all expected outputs exist and are
# indexed, and re-indexes any that are missing their .tbi.
#
# Usage: 16_index_outputs.sh --config <config.yml> --scripts-dir <dir> <cohort>
# =============================================================================

#BSUB -q short
#BSUB -n 2
#BSUB -M 8
#BSUB -R "rusage[mem=8]"
#BSUB -W 3:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/index.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/index.%J.err

set -eo pipefail

# ---- Config loading (Tier 2: --config and --scripts-dir as args 1-2) ----
CONFIG="${1:?Usage: --config <config.yml> --scripts-dir <dir> <cohort>}"
SCRIPTS_DIR="${2:?}"
shift 2
COHORT="${1:?Usage: <cohort>}"

# ---- Global init ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load all config values as shell vars
eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}" --cohort "${COHORT}")"

# ---- Paths from config ----
OUTPUT_BASE="${OUTPUT_BASE}"
OUTPUT_DIR="${OUTPUT_BASE}/${COHORT}/output"

echo "[$(date)] Verifying + indexing outputs (cohort $COHORT)"

conda activate samtools-1.16.1

# Expected .bed.gz outputs (stability has no phenotype_groups)
EXPECTED_BEDS=(
    "expression.bed.gz"
    "isoforms.bed.gz"
    "alt_TSS.bed.gz"
    "alt_polyA.bed.gz"
    "splicing.bed.gz"
    "intron_retention.bed.gz"
    "RNA_editing.bed.gz"
    "stability.bed.gz"
)

EXPECTED_GROUPS=(
    "isoforms.phenotype_groups.txt"
    "alt_TSS.phenotype_groups.txt"
    "alt_polyA.phenotype_groups.txt"
    "splicing.phenotype_groups.txt"
    "intron_retention.phenotype_groups.txt"
    "RNA_editing.phenotype_groups.txt"
)

ALL_OK=true

for bed in "${EXPECTED_BEDS[@]}"; do
    bedpath="${OUTPUT_DIR}/${bed}"
    if [ ! -f "$bedpath" ]; then
        echo "  MISSING: $bedpath"
        ALL_OK=false
        continue
    fi
    if [ ! -f "${bedpath}.tbi" ]; then
        echo "  Re-indexing: $bed (missing .tbi)"
        tabix -p bed "$bedpath"
    fi
    lines=$(zcat < "$bedpath" | wc -l)
    echo "  OK: $bed ($lines lines)"
done

echo ""
echo "Phenotype group files:"
for grp in "${EXPECTED_GROUPS[@]}"; do
    grppath="${OUTPUT_DIR}/${grp}"
    if [ ! -f "$grppath" ]; then
        echo "  MISSING: $grp"
        ALL_OK=false
    else
        lines=$(wc -l < "$grppath")
        echo "  OK: $grp ($lines lines)"
    fi
done

# RNA editing also has site_to_phenotype.tsv
if [ -f "${OUTPUT_DIR}/RNA_editing.site_to_phenotype.tsv" ]; then
    echo "  OK: RNA_editing.site_to_phenotype.tsv"
else
    echo "  MISSING: RNA_editing.site_to_phenotype.tsv"
    ALL_OK=false
fi

conda deactivate 2>/dev/null || true

if [ "$ALL_OK" = true ]; then
    echo ""
    echo "[$(date)] All outputs verified for cohort $COHORT"
else
    echo ""
    echo "[$(date)] WARNING: Some outputs missing for cohort $COHORT (see above)"
fi
