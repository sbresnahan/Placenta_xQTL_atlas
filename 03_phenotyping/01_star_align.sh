#!/bin/bash
# =============================================================================
# 01_star_align.sh — Per-sample STAR alignment + shrink BAM + index
# =============================================================================
# Usage: 01_star_align.sh --config <config.yml> --scripts-dir <dir> <cohort> <sample_id>
# The driver passes --config and --scripts-dir as the first two args.
#
# Reads FASTQ paths from the cohort's fastq_map.txt (tab-delimited:
# R1_path [R2_path] sample_id, paths relative to fastq_dir).
# Outputs:
#   <cohort_dir>/intermediate/star_out/<sample>.Aligned.sortedByCoord.out.bam
#   <cohort_dir>/intermediate/bam/<sample>.bam          (shrunk: SEQ/QUAL removed)
#   <cohort_dir>/intermediate/bam/<sample>.bam.bai
# =============================================================================

#BSUB -q medium
#BSUB -n 16
#BSUB -M 60
#BSUB -R "rusage[mem=60]"
#BSUB -W 8:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/align.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/align.%J.err

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

# Load all config values as shell vars (OUTPUT_BASE, STAR_INDEX, etc.)
eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}" --cohort "${COHORT}")"

# ---- Paths from config ----
OUTPUT_BASE="${OUTPUT_BASE}"
STAR_INDEX="${STAR_INDEX}"
COHORT_DIR="${OUTPUT_BASE}/${COHORT}"
INTERM_DIR="${COHORT_DIR}/intermediate"
FASTQ_MAP="${COHORT_DIR}/fastq_map.txt"
FASTQ_DIR="${COHORT_DIR}/fastq"

STAR_OUT_DIR="${INTERM_DIR}/star_out"
BAM_DIR="${INTERM_DIR}/bam"
mkdir -p "$STAR_OUT_DIR" "$BAM_DIR"

echo "[$(date)] Aligning $SAMPLE (cohort $COHORT)"

# ---- Parse FASTQ paths for this sample from fastq_map ----
# fastq_map format: R1_path \t R2_path \t sample_id  (paired)
#              or:   R1_path \t sample_id            (single)
# Paths are relative to FASTQ_DIR.
MAP_LINE=$(grep -P "\t${SAMPLE}\s*$" "$FASTQ_MAP" | head -1)
if [ -z "$MAP_LINE" ]; then
    echo "ERROR: Sample $SAMPLE not found in $FASTQ_MAP" >&2
    exit 1
fi

NFIELDS=$(echo "$MAP_LINE" | awk -F'\t' '{print NF}')
if [ "$NFIELDS" -eq 3 ]; then
    # Paired-end
    R1_REL=$(echo "$MAP_LINE" | awk -F'\t' '{print $1}')
    R2_REL=$(echo "$MAP_LINE" | awk -F'\t' '{print $2}')
    R1="${FASTQ_DIR}/${R1_REL}"
    R2="${FASTQ_DIR}/${R2_REL}"
    READ_FILES_IN="--readFilesIn ${R1} ${R2}"
    echo "  Paired-end: $R1 , $R2"
elif [ "$NFIELDS" -eq 2 ]; then
    # Single-end
    R1_REL=$(echo "$MAP_LINE" | awk -F'\t' '{print $1}')
    R1="${FASTQ_DIR}/${R1_REL}"
    READ_FILES_IN="--readFilesIn ${R1}"
    echo "  Single-end: $R1"
else
    echo "ERROR: Unexpected field count in fastq_map for $SAMPLE" >&2
    exit 1
fi

# ---- Determine readFilesCommand based on file extension ----
READ_FILES_CMD=""
if [[ "$R1" == *.gz ]]; then
    READ_FILES_CMD="--readFilesCommand zcat"
    echo "  Input: gzipped (using zcat)"
else
    echo "  Input: uncompressed (no readFilesCommand)"
fi

# ---- STAR alignment ----
conda activate star-2.7.4a

STAR_PREFIX="${STAR_OUT_DIR}/${SAMPLE}."
STAR \
    --runMode alignReads \
    --genomeDir "$STAR_INDEX" \
    $READ_FILES_IN \
    $READ_FILES_CMD \
    --twopassMode Basic \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "$STAR_PREFIX" \
    --runThreadN 16

conda deactivate 2>/dev/null || true

STAR_BAM="${STAR_OUT_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam"
if [ ! -f "$STAR_BAM" ]; then
    echo "ERROR: STAR output BAM not found: $STAR_BAM" >&2
    exit 1
fi

# ---- Shrink BAM (remove SEQ and QUAL fields to reduce size) ----
conda activate samtools-1.16.1

SHRUNK_BAM="${BAM_DIR}/${SAMPLE}.bam"
samtools view -h "$STAR_BAM" \
    | awk -v OFS="\t" '{if (substr($0, 1, 1) != "@") {$10="*"; $11="*"}; print}' \
    | samtools view -h -b \
    > "$SHRUNK_BAM"

# ---- Index the shrunk BAM ----
samtools index -@ 15 "$SHRUNK_BAM"

conda deactivate 2>/dev/null || true

echo "[$(date)] Done: $SAMPLE"
echo "  BAM: $SHRUNK_BAM"
echo "  BAI: ${SHRUNK_BAM}.bai"
ls -lh "$SHRUNK_BAM" "${SHRUNK_BAM}.bai"
