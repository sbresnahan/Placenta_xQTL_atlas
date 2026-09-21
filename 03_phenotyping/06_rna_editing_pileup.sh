#!/bin/bash
# =============================================================================
# 06_rna_editing_pileup.sh — Per-sample RNA editing level quantification
# =============================================================================
# Runs PANTRY's query_editing_level.py which uses samtools mpileup at the
# pre-defined edit sites. Uses the STAR BAM (full, with sequences) because
# mpileup needs base qualities — the shrunk BAM has SEQ/QUAL stripped.
#
# Usage: 06_rna_editing_pileup.sh --config <config.yml> --scripts-dir <dir> <cohort> <sample_id>
# Outputs:
#   <cohort_dir>/intermediate/RNA_editing/edit_levels/<sample>.rnaeditlevel.tsv.gz
# =============================================================================

#BSUB -q medium
#BSUB -n 4
#BSUB -M 16
#BSUB -R "rusage[mem=16]"
#BSUB -W 8:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/rnaedit.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/rnaedit.%J.err

set -eo pipefail

# ---- Config loading (Tier 2: --config and --scripts-dir as args 1-2) ----
CONFIG="${1:?Usage: --config <config.yml> --scripts-dir <dir> <cohort> <sample>}"
SCRIPTS_DIR="${2:?}"
shift 2
COHORT="${1:?Usage: <cohort> <sample_id>}"
SAMPLE="${2:?Usage: <cohort> <sample_id>}"

# ---- Global init ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load all config values as shell vars
eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}" --cohort "${COHORT}")"

# ---- Paths from config ----
OUTPUT_BASE="${OUTPUT_BASE}"
REFERENCE_DIR="${REFERENCE_DIR}"
REF_GENOME="${REF_GENOME}"
PANTRY_SCRIPTS="${PANTRY_SCRIPTS}"
COHORT_DIR="${OUTPUT_BASE}/${COHORT}"
INTERM_DIR="${COHORT_DIR}/intermediate"

# Use the STAR output BAM (has SEQ/QUAL), NOT the shrunk BAM
STAR_BAM="${INTERM_DIR}/star_out/${SAMPLE}.Aligned.sortedByCoord.out.bam"
EDIT_SITES="${REFERENCE_DIR}/rediportal_hg38.bed"
EDIT_LEVELS_DIR="${INTERM_DIR}/RNA_editing/edit_levels"
OUT="${EDIT_LEVELS_DIR}/${SAMPLE}.rnaeditlevel.tsv.gz"
mkdir -p "$EDIT_LEVELS_DIR"

if [ ! -f "$STAR_BAM" ]; then
    echo "ERROR: STAR BAM not found for $SAMPLE: $STAR_BAM" >&2
    exit 1
fi
if [ ! -f "$EDIT_SITES" ]; then
    echo "ERROR: Edit sites BED not found: $EDIT_SITES" >&2
    echo "  Run Tier 1 reference prep (001/002/003) with rediportal_input set." >&2
    exit 1
fi

echo "[$(date)] RNA editing pileup for $SAMPLE (cohort $COHORT)"

# query_editing_level.py calls samtools mpileup internally, so we need
# samtools on PATH. It also needs Python (pandas not required here, but
# the MAJIQ env provides a clean Python).
conda activate samtools-1.16.1
source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate

python3 "$PANTRY_SCRIPTS/RNA_editing/query_editing_level.py" \
    --edit_sites "$EDIT_SITES" \
    --ref "$REF_GENOME" \
    --bam "$STAR_BAM" \
    --output "$OUT"

conda deactivate 2>/dev/null || true

echo "[$(date)] Done: $SAMPLE"
echo "  Output: $OUT"
