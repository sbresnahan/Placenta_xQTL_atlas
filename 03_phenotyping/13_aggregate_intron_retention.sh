#!/bin/bash
# =============================================================================
# 13_aggregate_intron_retention.sh — MAJIQ build + psi + IR extraction + assemble
# =============================================================================
# This is the most complex aggregation: MAJIQ build runs at the COHORT level
# (needs all BAMs in one experiments.tsv), then psi-coverage in batches,
# then psi, then extract retained introns, then assemble + normalize.
#
# Usage: 13_aggregate_intron_retention.sh --config <config.yml> --scripts-dir <dir> <cohort>
# Outputs:
#   <cohort_dir>/intermediate/intron_retention/splicegraph.zarr/
#   <cohort_dir>/intermediate/intron_retention/build/*.sj
#   <cohort_dir>/intermediate/intron_retention/psicov/batch_*.psicov.zarr/
#   <cohort_dir>/intermediate/intron_retention/majiq.psi.tsv
#   <cohort_dir>/intermediate/intron_retention/retained_intron_psi.tsv.gz
#   <cohort_dir>/output/unnorm/intron_retention.bed
#   <cohort_dir>/output/intron_retention.bed.gz (+ .tbi)
#   <cohort_dir>/output/intron_retention.phenotype_groups.txt
# =============================================================================

#BSUB -q long
#BSUB -n 16
#BSUB -M 64
#BSUB -R "rusage[mem=64]"
#BSUB -W 25:00
#BSUB -o /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/agg_ir.%J.out
#BSUB -e /rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/PANTRY/logs/agg_ir.%J.err

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
REFERENCE_DIR="${REFERENCE_DIR}"
# Use the NORMALIZED GTF (has gene_name, gene_biotype, gene features required
# by assemble_bed.py), not the raw SQANTI3 GTF.
REF_ANNO="${NORMALIZED_GTF}"
PANTRY_SCRIPTS="${PANTRY_SCRIPTS}"
MAJIQ_LICENSE="${MAJIQ_LICENSE}"
COHORT_DIR="${OUTPUT_BASE}/${COHORT}"
INTERM_DIR="${COHORT_DIR}/intermediate"
OUTPUT_DIR="${COHORT_DIR}/output"
UNNORM_DIR="${OUTPUT_DIR}/unnorm"
SAMPLES_FILE="${COHORT_DIR}/samples.txt"

# MAJIQ parameters (from config)
PSICOV_BATCH_SIZE="${INTRON_RETENTION_PSICOV_BATCH_SIZE}"
MIN_EXPERIMENTS="${INTRON_RETENTION_MIN_EXPERIMENTS}"
MAJIQ_THREADS=16

MAJIQ_GFF3="${REFERENCE_DIR}/majiq/annotation.gff3"
MAJIQ_DIR="${INTERM_DIR}/intron_retention"
MAJIQ_BUILD_DIR="${MAJIQ_DIR}/build"
MAJIQ_PSICOV_DIR="${MAJIQ_DIR}/psicov"
BAM_DIR="${INTERM_DIR}/bam"

mkdir -p "$MAJIQ_BUILD_DIR" "$MAJIQ_PSICOV_DIR" "$UNNORM_DIR" "$OUTPUT_DIR"

if [ ! -f "$MAJIQ_GFF3" ]; then
    echo "ERROR: MAJIQ GFF3 not found: $MAJIQ_GFF3" >&2
    echo "  Run Tier 1 reference prep (001/002/003) first." >&2
    exit 1
fi

echo "[$(date)] Aggregating intron retention (cohort $COHORT)"

# ---- Load MAJIQ env (samtools first, then MAJIQ) ----
# MAJIQ's rna_majiq Python package (Python 3.12) is compiled against HTSlib
# and dynamically links libhts.so.3 at runtime. The samtools-1.16.1 conda env
# provides libhts.so.3 in its lib/ directory, but conda activate does NOT add
# that directory to LD_LIBRARY_PATH on seadragon (it only puts the LSF
# scheduler's lib dir there). Without the explicit export below, MAJIQ fails
# with: ImportError: libhts.so.3: cannot open shared object file.
#
# Fix: after conda activate, explicitly prepend $CONDA_PREFIX/lib to
# LD_LIBRARY_PATH. $CONDA_PREFIX is set by conda to the active env's root
# (e.g. /risapps/rhel8/miniforge3/24.5.0-0/envs/samtools-1.16.1). MAJIQ's own
# venv activate script does not clobber LD_LIBRARY_PATH, so the path persists.
conda activate samtools-1.16.1
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
# MAJIQ activate path (infra — hardcoded)
source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate

# ---- Write experiments.tsv (group \t bam_path) ----
EXPERIMENTS_TSV="${MAJIQ_DIR}/experiments.tsv"
echo -e "group\tpath" > "$EXPERIMENTS_TSV"
while IFS= read -r sample; do
    [ -z "$sample" ] && continue
    echo -e "pantry\t${BAM_DIR}/${sample}.bam" >> "$EXPERIMENTS_TSV"
done < "$SAMPLES_FILE"
echo "  experiments.tsv: $(($(wc -l < "$EXPERIMENTS_TSV") - 1)) samples"

# ---- MAJIQ build (all BAMs at once) ----
echo "[$(date)] MAJIQ build"
rm -rf "${MAJIQ_BUILD_DIR}/splicegraph.zarr"
rm -f "${MAJIQ_BUILD_DIR}"/*.sj

majiq build \
    --license "$MAJIQ_LICENSE" \
    --all-introns \
    --min-experiments "$MIN_EXPERIMENTS" \
    --nthreads "$MAJIQ_THREADS" \
    --overwrite \
    "$MAJIQ_GFF3" \
    "$EXPERIMENTS_TSV" \
    "$MAJIQ_BUILD_DIR"

# ---- Determine batches ----
N_SAMPLES=$(($(wc -l < "$SAMPLES_FILE")))
N_BATCHES=$(( (N_SAMPLES + PSICOV_BATCH_SIZE - 1) / PSICOV_BATCH_SIZE ))
echo "  $N_SAMPLES samples -> $N_BATCHES batches (batch size $PSICOV_BATCH_SIZE)"

# ---- MAJIQ psi-coverage (batched) ----
echo "[$(date)] MAJIQ psi-coverage ($N_BATCHES batches)"
SAMPLES=($(cat "$SAMPLES_FILE"))
for b in $(seq 0 $((N_BATCHES - 1))); do
    start=$((b * PSICOV_BATCH_SIZE))
    batch_samples=("${SAMPLES[@]:start:PSICOV_BATCH_SIZE}")
    prefixes=$(IFS=' '; echo "${batch_samples[*]}")
    sj_files=""
    for s in "${batch_samples[@]}"; do
        sj_files="${sj_files} ${MAJIQ_BUILD_DIR}/${s}.sj"
    done

    psicov_out="${MAJIQ_PSICOV_DIR}/batch_${b}.psicov.zarr"
    rm -rf "$psicov_out"
    echo "  Batch $b: ${#batch_samples[@]} samples"
    majiq psi-coverage \
     --license "$MAJIQ_LICENSE" \
     --nthreads "$MAJIQ_THREADS" \
     --prefixes $prefixes \
     --overwrite \
     "${MAJIQ_BUILD_DIR}/splicegraph.zarr" \
     "$psicov_out" \
     $sj_files
done

# ---- MAJIQ psi (across all batches) ----
echo "[$(date)] MAJIQ psi"
PSICOV_ARGS=""
for b in $(seq 0 $((N_BATCHES - 1))); do
    PSICOV_ARGS="${PSICOV_ARGS} ${MAJIQ_PSICOV_DIR}/batch_${b}.psicov.zarr"
done

majiq psi \
    --license "$MAJIQ_LICENSE" \
    --splicegraph "${MAJIQ_BUILD_DIR}/splicegraph.zarr" \
    --output-tsv "${MAJIQ_DIR}/majiq.psi.tsv" \
    --quantiles 0.025 0.5 0.975 \
    --nthreads "$MAJIQ_THREADS" \
    --overwrite \
    $PSICOV_ARGS

# ---- Extract retained introns ----
echo "[$(date)] Extracting retained intron PSI"
python3 "${PANTRY_SCRIPTS}/intron_retention/extract_ir_psi.py" \
    --input "${MAJIQ_DIR}/majiq.psi.tsv" \
    --output "${MAJIQ_DIR}/retained_intron_psi.tsv.gz"

# ---- Assemble BED ----
python3 "${PANTRY_SCRIPTS}/assemble_bed.py" intron-retention \
    --input "${MAJIQ_DIR}/retained_intron_psi.tsv.gz" \
    --ref-anno "$REF_ANNO" \
    --output "${UNNORM_DIR}/intron_retention.bed"

# ---- Canonical BED = unnorm (2026-09 schema: normalization moved to stage 5) ----
cp "${UNNORM_DIR}/intron_retention.bed" "${OUTPUT_DIR}/intron_retention.bed"

conda deactivate 2>/dev/null || true

# ---- bgzip + tabix + phenotype groups ----
conda activate samtools-1.16.1
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

bgzip "${OUTPUT_DIR}/intron_retention.bed"
tabix -p bed "${OUTPUT_DIR}/intron_retention.bed.gz"

zcat < "${OUTPUT_DIR}/intron_retention.bed.gz" \
    | tail -n +2 \
    | cut -f4 \
    | awk '{ g=$1; sub(/__.*$/, "", g); print $1 "\t" g }' \
    > "${OUTPUT_DIR}/intron_retention.phenotype_groups.txt"

conda deactivate 2>/dev/null || true

echo "[$(date)] Done: intron retention"
ls -lh "${OUTPUT_DIR}/intron_retention.bed.gz" "${OUTPUT_DIR}/intron_retention.phenotype_groups.txt"
