#!/bin/bash
# =============================================================================
# 05_featureCounts.sh — Per-sample exon/intron counts (for stability)
# =============================================================================
# Runs featureCounts twice: once on exonic.gtf (constitutive exons) and once
# on intronic.gtf (introns). The exon/intron read ratio = mRNA stability.
#
# Strandedness comes from the cohort config: -s 0 (unstranded), 1 (stranded),
# 2 (reverse stranded). The driver passes the -s value as $3 (after
# --config/--scripts-dir shift).
#
# Usage: 05_featureCounts.sh --config <config.yml> --scripts-dir <dir> <cohort> <sample_id> <feature_type> <strandedness>
#   feature_type:  exonic or intronic
#   strandedness:  0, 1, or 2
# Outputs:
#   <cohort_dir>/intermediate/stability/<sample>.<feature_type>.counts.txt
# =============================================================================

#BSUB -q medium
#BSUB -n 8
#BSUB -M 16
#BSUB -R "rusage[mem=16]"
#BSUB -W 4:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/featcounts.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/featcounts.%J.err

set -eo pipefail

# ---- Config loading (Tier 2: --config and --scripts-dir as args 1-2) ----
CONFIG="${1:?Usage: --config <config.yml> --scripts-dir <dir> <cohort> <sample> <feature_type> <strandedness>}"
SCRIPTS_DIR="${2:?}"
shift 2
COHORT="${1:?Usage: <cohort> <sample> <feature_type> <strandedness>}"
SAMPLE="${2:?Usage: <cohort> <sample> <feature_type> <strandedness>}"
FEATURE_TYPE="${3:?Usage: feature_type must be exonic or intronic}"
# NOTE: named FC_STRAND, not STRANDEDNESS, to avoid being clobbered by the
# config_get.py eval below, which exports the cohort's STRANDEDNESS config
# field (the raw word, e.g. "unstranded") as a shell var.
FC_STRAND="${4:?Usage: strandedness must be 0, 1, or 2}"

# ---- Global init ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# Load all config values as shell vars
eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}" --cohort "${COHORT}")"

# ---- Paths from config ----
OUTPUT_BASE="${OUTPUT_BASE}"
REFERENCE_DIR="${REFERENCE_DIR}"
COHORT_DIR="${OUTPUT_BASE}/${COHORT}"
INTERM_DIR="${COHORT_DIR}/intermediate"
FASTQ_MAP="${COHORT_DIR}/fastq_map.txt"

BAM="${INTERM_DIR}/bam/${SAMPLE}.bam"
GTF="${REFERENCE_DIR}/${FEATURE_TYPE}.gtf"
STAB_DIR="${INTERM_DIR}/stability"
OUT="${STAB_DIR}/${SAMPLE}.${FEATURE_TYPE}.counts.txt"
mkdir -p "$STAB_DIR"

if [ ! -f "$BAM" ]; then
    echo "ERROR: BAM not found for $SAMPLE: $BAM" >&2
    exit 1
fi
if [ ! -f "$GTF" ]; then
    echo "ERROR: GTF not found: $GTF" >&2
    echo "  Run Tier 1 reference prep (001/002/003) first." >&2
    exit 1
fi

# Determine feature ID and overlap fraction per PANTRY:
#   exonic  -> -t exon  --fracOverlap 1  (read must be fully within exon)
#   intronic -> -t intron --fracOverlap 0 (read need only overlap intron by 1 base)
if [ "$FEATURE_TYPE" = "exonic" ]; then
    FEATURE_ID="exon"
    FRAC_OVERLAP=1
elif [ "$FEATURE_TYPE" = "intronic" ]; then
    FEATURE_ID="intron"
    FRAC_OVERLAP=0
else
    echo "ERROR: feature_type must be 'exonic' or 'intronic', got: $FEATURE_TYPE" >&2
    exit 1
fi

# Determine paired-end flag from fastq_map
MAP_LINE=$(grep -P "\t${SAMPLE}\s*$" "$FASTQ_MAP" | head -1)
NFIELDS=$(echo "$MAP_LINE" | awk -F'\t' '{print NF}')
if [ "$NFIELDS" -eq 3 ]; then
    PAIRED_FLAG="-p"
else
    PAIRED_FLAG=""
fi

echo "[$(date)] featureCounts $SAMPLE $FEATURE_TYPE (cohort $COHORT, -s $FC_STRAND)"

module load subread

featureCounts \
    "$BAM" \
    $PAIRED_FLAG \
    -a "$GTF" \
    -t "$FEATURE_ID" \
    -s "$FC_STRAND" \
    --fracOverlap "$FRAC_OVERLAP" \
    -T 8 \
    -o "$OUT"

echo "[$(date)] Done: $SAMPLE $FEATURE_TYPE"
echo "  Output: $OUT"
tail -3 "$OUT"