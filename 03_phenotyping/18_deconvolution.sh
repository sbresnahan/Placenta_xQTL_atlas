#!/bin/bash
# =============================================================================
# 18_deconvolution.sh — LSF driver for cell-type deconvolution (Objective 1.4)
# =============================================================================
# Submits three stages of LSF jobs for placental cell-type deconvolution using
# MuSiC with the Campbell et al. 2023 scRNA-seq reference (GSE182381):
#
#   Stage 1 (single job):    Build MuSiC reference SCE from GEO GSE182381
#   Stage 2 (per-cohort):    Build bulk TPM matrix + run MuSiC deconvolution
#   Stage 3 (single job):    Pool per-cohort outputs + QC report
#
# Dependencies: Stage 2 jobs depend on Stage 1 (done(refbuild_jobid)).
#               Stage 3 depends on all Stage 2 jobs (done(j1) && done(j2) && ...).
#
# This driver is submitted via stdin (bsub < 18_deconvolution.sh), so the
# #BSUB directives below ARE parsed by LSF. The stage scripts (18a/18b/18c)
# are submitted via the argument form (bsub ... script.sh), so their #BSUB
# directives are NOT parsed — resources are passed on the bsub command line.
#
# Usage:
#   bsub < 18_deconvolution.sh
#
# Or with explicit env (if your site restricts default propagation):
#   bsub -env "CONFIG=\"/path/to/config.yml\",SCRIPTS_DIR=\"/path/to/scripts\"" < 18_deconvolution.sh
#
# Required environment:
#   CONFIG       Path to config.yml (default: seadragon PANTRY config)
#   SCRIPTS_DIR  Directory containing the deconvolution scripts
#                (build_music_reference.R, build_bulk_tpm_matrix.py,
#                 run_music_deconvolution.R, pool_deconvolution_outputs.py,
#                 18a_build_reference.sh, 18b_deconvolve_cohort.sh,
#                 18c_pool_outputs.sh, config_get.py)
#
# Optional environment:
#   N_PSEUDO           Number of pseudo-subjects for reference (default: 5)
#   MATERNAL_THRESHOLD  Maternal fraction flag threshold (default: 0.10)
#   ANCESTRY_MAP       Optional TSV: sample_id, ancestry (for QC stratification)
#   QUEUES             LSF queue for stage jobs (default: medium)
#
# Wrapper resources: long queue, 2 cores, 8 GB, 4h. The wrapper itself is
# lightweight (bsub submissions only) and exits after submitting all jobs.
#
# Create log dir once before first use:
#   mkdir -p /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/logs
#
# PREREQUISITE: MuSiC, TOAST, and SingleCellExperiment must be installed in
#   /rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1
#   (the seadragon R singularity image provides base R 4.3.1 but NOT these
#   Bioconductor packages). Install with:
#     SING_R="singularity exec --bind /rsrch5 --bind /rsrch9 \
#       /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif Rscript"
#     export R_LIBS_USER="/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"
#     $SING_R -e 'BiocManager::install(c("MuSiC","TOAST","SingleCellExperiment"))'
# =============================================================================

#BSUB -J deconv_run
#BSUB -q long
#BSUB -n 2
#BSUB -M 8
#BSUB -W 25:00
#BSUB -o /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/logs/deconvolution.%J.out
#BSUB -e /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/logs/deconvolution.%J.err

set -eo pipefail

# --- Defaults ---
CONFIG="${CONFIG:-/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/config.yml}"
# Default SCRIPTS_DIR to this script's own directory (repo-clone layout).
SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
N_PSEUDO="${N_PSEUDO:-5}"
MATERNAL_THRESHOLD="${MATERNAL_THRESHOLD:-0.10}"
ANCESTRY_MAP="${ANCESTRY_MAP:-}"
QUEUES="${QUEUES:-medium}"

LOG_DIR="/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/logs"
mkdir -p "${LOG_DIR}"

# --- Validate ---
if [ ! -f "${CONFIG}" ]; then
    echo "ERROR: config.yml not found at ${CONFIG}" >&2
    exit 1
fi

for script in 18a_build_reference.sh 18b_deconvolve_cohort.sh 18c_pool_outputs.sh \
              build_music_reference.R build_bulk_tpm_matrix.py \
              run_music_deconvolution.R pool_deconvolution_outputs.py \
              config_get.py; do
    if [ ! -f "${SCRIPTS_DIR}/${script}" ]; then
        echo "ERROR: ${script} not found in ${SCRIPTS_DIR}" >&2
        exit 1
    fi
done

# --- Global init (seadragon environment) ---
# The driver needs Python with PyYAML to read config.yml for the cohort list.
# config_get.py uses a built-in parser (no PyYAML), but we use PyYAML here
# for robust YAML parsing of the cohorts section.
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate samtools-1.16.1
source /rsrch5/home/epi/bhattacharya_lab/software/MAJIQ/bin/activate

# Load top-level config values (OUTPUT_BASE, etc.)
eval "$(python3 "${SCRIPTS_DIR}/config_get.py" "${CONFIG}")"

OUTPUT_BASE="${OUTPUT_BASE}"
DECONV_DIR="${OUTPUT_BASE}/deconvolution"
REF_RDS="${DECONV_DIR}/placenta_music_reference.rds"
mkdir -p "${DECONV_DIR}"

echo "============================================================"
echo "Cell-type deconvolution pipeline (Objective 1.4)"
echo "  Config:        ${CONFIG}"
echo "  Scripts dir:   ${SCRIPTS_DIR}"
echo "  Output dir:    ${DECONV_DIR}"
echo "  Pseudo-subj:   ${N_PSEUDO}"
echo "  Mat threshold: ${MATERNAL_THRESHOLD}"
echo "  Ancestry map:  ${ANCESTRY_MAP:-<none>}"
echo "  Queue:         ${QUEUES}"
echo "  Started:       $(date)"
echo "  LSF job:       ${LSB_JOBID:-<not under LSF>}"
echo "============================================================"

# --- Read cohort list from config.yml ---
COHORTS=$(python3 -c "
import yaml
with open('${CONFIG}') as f:
    cfg = yaml.safe_load(f)
cohorts = list(cfg.get('cohorts', {}).keys())
print(' '.join(cohorts))
")
echo "Cohorts: ${COHORTS}"

# =============================================================================
# Helper: parse job ID from bsub output
# =============================================================================
parse_job_id() {
    echo "$1" | sed -n 's/.*<\([0-9]*\)>.*/\1/p'
}

# =============================================================================
# Helper: build LSF done() dependency expression from a list of job IDs
# =============================================================================
build_dep_expr() {
    # Args: space-separated list of job IDs
    # Returns: done(j1) && done(j2) && ...  (or empty if no IDs)
    local ids="$1"
    local expr=""
    for id in ${ids}; do
        if [ -z "${expr}" ]; then
            expr="done(${id})"
        else
            expr="${expr} && done(${id})"
        fi
    done
    echo "${expr}"
}

# =============================================================================
# Stage 1: Build MuSiC reference (single job)
# =============================================================================
echo ""
echo "--- Submitting Stage 1: Reference build ---"

if [ -f "${REF_RDS}" ]; then
    echo "  Reference already exists: ${REF_RDS}"
    echo "  Skipping Stage 1. (Delete file to rebuild.)"
    REFBUILD_JOB_ID=""
else
    REFBUILD_CMD=$(bsub -J "refbuild" \
         -q "${QUEUES}" \
         -n 4 \
         -M 16 \
         -R "rusage[mem=16]" \
         -W 4:00 \
         -o "${LOG_DIR}/deconv_refbuild.%J.out" \
         -e "${LOG_DIR}/deconv_refbuild.%J.err" \
         "${SCRIPTS_DIR}/18a_build_reference.sh" \
             --config "${CONFIG}" \
             --scripts-dir "${SCRIPTS_DIR}" \
             --output "${REF_RDS}" \
             --n-pseudo-subjects "${N_PSEUDO}" \
             --seed 42 2>&1)

    REFBUILD_JOB_ID=$(parse_job_id "${REFBUILD_CMD}")

    if [ -z "${REFBUILD_JOB_ID}" ]; then
        echo "ERROR: Failed to submit refbuild job" >&2
        echo "  bsub output: ${REFBUILD_CMD}" >&2
        exit 1
    fi
    echo "  Submitted refbuild -> job ${REFBUILD_JOB_ID}"
fi

# =============================================================================
# Stage 2: Per-cohort deconvolution (depends on refbuild)
# =============================================================================
echo ""
echo "--- Submitting Stage 2: Per-cohort deconvolution ---"

# Build dependency flag for Stage 2 jobs (all depend on refbuild)
STAGE2_DEP_IDS="${REFBUILD_JOB_ID}"
STAGE2_DEP_EXPR=$(build_dep_expr "${STAGE2_DEP_IDS}")

DECONV_JOB_IDS=""
SUBMITTED_COHORTS=""

for cohort in ${COHORTS}; do
    echo "  Cohort ${cohort}:"

    # Build bsub command — conditionally include -w if there's a dependency
    if [ -n "${STAGE2_DEP_EXPR}" ]; then
        DECONV_CMD=$(bsub -w "${STAGE2_DEP_EXPR}" \
             -J "deconv_${cohort}" \
             -q "${QUEUES}" \
             -n 4 \
             -M 16 \
             -R "rusage[mem=16]" \
             -W 4:00 \
             -o "${LOG_DIR}/deconv_${cohort}.%J.out" \
             -e "${LOG_DIR}/deconv_${cohort}.%J.err" \
             "${SCRIPTS_DIR}/18b_deconvolve_cohort.sh" \
                 --config "${CONFIG}" \
                 --scripts-dir "${SCRIPTS_DIR}" \
                 --cohort "${cohort}" \
                 --reference "${REF_RDS}" \
                 --maternal-threshold "${MATERNAL_THRESHOLD}" 2>&1)
    else
        DECONV_CMD=$(bsub \
             -J "deconv_${cohort}" \
             -q "${QUEUES}" \
             -n 4 \
             -M 16 \
             -R "rusage[mem=16]" \
             -W 4:00 \
             -o "${LOG_DIR}/deconv_${cohort}.%J.out" \
             -e "${LOG_DIR}/deconv_${cohort}.%J.err" \
             "${SCRIPTS_DIR}/18b_deconvolve_cohort.sh" \
                 --config "${CONFIG}" \
                 --scripts-dir "${SCRIPTS_DIR}" \
                 --cohort "${cohort}" \
                 --reference "${REF_RDS}" \
                 --maternal-threshold "${MATERNAL_THRESHOLD}" 2>&1)
    fi

    JOB_ID=$(parse_job_id "${DECONV_CMD}")

    if [ -z "${JOB_ID}" ]; then
        echo "    ERROR: Failed to submit deconv_${cohort}" >&2
        echo "    bsub output: ${DECONV_CMD}" >&2
        # Continue with other cohorts rather than aborting
    else
        echo "    Submitted deconv_${cohort} -> job ${JOB_ID}"
        DECONV_JOB_IDS="${DECONV_JOB_IDS} ${JOB_ID}"
        SUBMITTED_COHORTS="${SUBMITTED_COHORTS} ${cohort}"
    fi
done

# =============================================================================
# Stage 3: Pool outputs + QC report (depends on all Stage 2 jobs)
# =============================================================================
echo ""
echo "--- Submitting Stage 3: Pool outputs + QC report ---"

if [ -z "${DECONV_JOB_IDS}" ]; then
    echo "  WARNING: No Stage 2 jobs submitted. Skipping Stage 3." >&2
    POOL_JOB_ID=""
else
    # Build dependency: done(job1) && done(job2) && ...
    POOL_DEP_EXPR=$(build_dep_expr "${DECONV_JOB_IDS}")

    # Build pool script args
    POOL_SCRIPT_ARGS=(
        --config "${CONFIG}"
        --scripts-dir "${SCRIPTS_DIR}"
    )
    if [ -n "${ANCESTRY_MAP}" ]; then
        POOL_SCRIPT_ARGS+=(--ancestry-map "${ANCESTRY_MAP}")
    fi

    POOL_CMD=$(bsub -w "${POOL_DEP_EXPR}" \
         -J "deconv_pool" \
         -q "short" \
         -n 2 \
         -M 8 \
         -R "rusage[mem=8]" \
         -W 3:00 \
         -o "${LOG_DIR}/deconv_pool.%J.out" \
         -e "${LOG_DIR}/deconv_pool.%J.err" \
         "${SCRIPTS_DIR}/18c_pool_outputs.sh" \
             "${POOL_SCRIPT_ARGS[@]}" 2>&1)

    POOL_JOB_ID=$(parse_job_id "${POOL_CMD}")

    if [ -z "${POOL_JOB_ID}" ]; then
        echo "  ERROR: Failed to submit pool job" >&2
        echo "  bsub output: ${POOL_CMD}" >&2
    else
        echo "  Submitted deconv_pool -> job ${POOL_JOB_ID}"
    fi
fi

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "============================================================"
echo "All jobs submitted."
echo "  Stage 1: refbuild"
if [ -n "${REFBUILD_JOB_ID}" ]; then
    echo "    job ${REFBUILD_JOB_ID}"
else
    echo "    (skipped — reference already exists)"
fi
echo "  Stage 2: ${SUBMITTED_COHORTS}"
for jid in ${DECONV_JOB_IDS}; do
    echo "    job ${jid}"
done
echo "  Stage 3: deconv_pool"
if [ -n "${POOL_JOB_ID}" ]; then
    echo "    job ${POOL_JOB_ID}"
else
    echo "    (skipped)"
fi
echo ""
echo "  Results will appear in: ${DECONV_DIR}"
echo "  Monitor: bjobs"
echo "  Ended: $(date)"
echo "============================================================"
