#!/bin/bash
# =============================================================================
# 17_harmonize_within_ancestry.sh — Cross-cohort coordinate harmonization
# for splicing & intron retention within an ancestry stratum.
# =============================================================================
# Pools per-cohort splicing (leafCutter) and intron retention (MAJIQ)
# intermediates, harmonizes features by stable genomic coordinates, and writes
# pooled unnorm/ BED files ready for downstream ComBat + normalization.
#
# No modifications to scripts 12/13 or assemble_bed.py. Reuses assemble_bed.py
# on pooled intermediates.
#
# Required env vars (passed via LSF -env or shell):
#   CONFIG        — path to config.yml
#   SCRIPTS_DIR   — directory containing config_get.py, harmonize_within_ancestry.py,
#                   and the stage scripts
#   ANCESTRY_MAP  — TSV mapping sample_id -> ancestry (columns: sample_id, ancestry)
#   ANCESTRY      — ancestry group to process (e.g. EUR, EAS, AFR, HIS, SAS)
#
# Optional env vars:
#   MODALITY      — splicing | intron-retention | both  (default: both)
#   OUTPUT_BASE_OVERRIDE — override output base from config (rarely needed)
#
# Usage:
#   bsub -env "CONFIG=\"config.yml\",SCRIPTS_DIR=\"/path/to/scripts\",\
# ANCESTRY_MAP=\"/path/to/samples_ancestry.tsv\",ANCESTRY=\"EUR\",MODALITY=\"both\"" \
#        < 17_harmonize_within_ancestry.sh
#
# Outputs (per ancestry, per modality):
#   {output_base}/{ancestry}/splicing/
#     intermediate/leafcutter_harmonized_perind_numers.counts.gz
#     intermediate/harmonization_map.tsv
#     unnorm/splicing.bed
#     unnorm/splicing.phenotype_groups.txt
#   {output_base}/{ancestry}/intron_retention/
#     intermediate/retained_intron_psi_harmonized.tsv.gz
#     intermediate/harmonization_map.tsv
#     unnorm/intron_retention.bed
#     unnorm/intron_retention.phenotype_groups.txt
#   {output_base}/{ancestry}/sample_namespace_map.tsv
# =============================================================================

#BSUB -q medium
#BSUB -n 8
#BSUB -M 32
#BSUB -R "rusage[mem=32]"
#BSUB -W 4:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/harmonize.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/harmonize.%J.err

set -eo pipefail

# ---- Required env vars ----
CONFIG="${CONFIG:?Usage: set CONFIG, SCRIPTS_DIR, ANCESTRY_MAP, ANCESTRY env vars}"
SCRIPTS_DIR="${SCRIPTS_DIR:?}"
ANCESTRY_MAP="${ANCESTRY_MAP:?}"
ANCESTRY="${ANCESTRY:?}"
MODALITY="${MODALITY:-both}"

# ---- Global init (same as scripts 10-15) ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load config values as shell vars (need OUTPUT_BASE for output paths)
eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}")"

OUTPUT_BASE="${OUTPUT_BASE}"
ANCESTRY_OUT="${OUTPUT_BASE}/${ANCESTRY}"

echo "[$(date)] Harmonizing within ancestry=${ANCESTRY}, modality=${MODALITY}"
echo "  CONFIG=${CONFIG}"
echo "  ANCESTRY_MAP=${ANCESTRY_MAP}"
echo "  OUTPUT_BASE=${OUTPUT_BASE}"
echo "  ANCESTRY_OUT=${ANCESTRY_OUT}"

# ---- Activate Python env (pandas, numpy, scipy, gtfparse) ----
# Same env stack as the aggregation scripts: samtools-1.16.1 + MAJIQ venv
conda activate samtools-1.16.1
source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate

HARM_SCRIPT="${SCRIPTS_DIR}/harmonize_within_ancestry.py"

# ---- Run harmonization ----
if [ "${MODALITY}" = "splicing" ] || [ "${MODALITY}" = "both" ]; then
    echo "[$(date)] === Splicing harmonization ==="
    python3 "${HARM_SCRIPT}" splicing \
        --config "${CONFIG}" \
        --ancestry-map "${ANCESTRY_MAP}" \
        --ancestry "${ANCESTRY}" \
        --scripts-dir "${SCRIPTS_DIR}" \
        --output-dir "${ANCESTRY_OUT}/splicing"
fi

if [ "${MODALITY}" = "intron-retention" ] || [ "${MODALITY}" = "both" ]; then
    echo "[$(date)] === Intron retention harmonization ==="
    python3 "${HARM_SCRIPT}" intron-retention \
        --config "${CONFIG}" \
        --ancestry-map "${ANCESTRY_MAP}" \
        --ancestry "${ANCESTRY}" \
        --scripts-dir "${SCRIPTS_DIR}" \
        --output-dir "${ANCESTRY_OUT}/intron_retention"
fi

conda deactivate 2>/dev/null || true

# ---- bgzip + tabix the unnorm BED outputs ----
conda activate samtools-1.16.1

for mod in splicing intron_retention; do
    if [ "${MODALITY}" = "${mod}" ] || [ "${MODALITY}" = "both" ]; then
        BED="${ANCESTRY_OUT}/${mod}/unnorm/${mod}.bed"
        if [ -f "${BED}" ]; then
            echo "[$(date)] bgzip + tabix: ${BED}"
            bgzip "${BED}"
            tabix -p bed "${BED}.gz"
        else
            echo "[WARN] ${BED} not found, skipping index" >&2
        fi
    fi
done

conda deactivate 2>/dev/null || true

echo "[$(date)] Done: harmonization for ancestry=${ANCESTRY}, modality=${MODALITY}"
echo "Outputs under: ${ANCESTRY_OUT}/"
