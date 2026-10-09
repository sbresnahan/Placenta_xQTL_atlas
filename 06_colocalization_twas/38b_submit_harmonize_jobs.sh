#!/usr/bin/env bash
# 38b_submit_harmonize_jobs.sh — submit one LSF job per pending GWAS trait
# for 38_harmonize_gwas.py.
#
# 38 skips traits whose {trait_id}.sumstats.tsv.gz (+ .tbi) already exists,
# and this script does the same check before submitting, so it is safe to
# re-run after failures: only missing traits are submitted. Each job writes
# its own per-trait QC file ($GWAS_DIR/qc/{trait_id}.qc.tsv) and refreshes
# the aggregated harmonization_qc.tsv when it finishes.
#
# Usage:
#   bash 38b_submit_harmonize_jobs.sh                 # all pending traits
#   bash 38b_submit_harmonize_jobs.sh <trait_id> ...  # only these traits
#
# Env overrides: GWAS_DIR PGEN_DIR CHAIN QUEUE MEM_GB WALL CONDA_ENV
#
# Environment: jobs run in the tensorqtl conda env (python3 + pandas) with
# the samtools module loaded (bgzip/tabix). The submitter itself only runs
# bsub and a tiny catalog parse, so it is fine on a login node.
set -euo pipefail

SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
OUTPUT_BASE="${OUTPUT_BASE:-/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY}"
GWAS_DIR="${GWAS_DIR:-${OUTPUT_BASE}/gwas}"
PGEN_DIR="${PGEN_DIR:-/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/pooled/genotypes}"
CHAIN="${CHAIN:-/home/stbresnahan/bhattacharya_lab/data/GenomicReferences/liftover/hg19ToHg38.over.chain}"
CONDA_BIN="${CONDA_BIN:-/risapps/rhel8/miniforge3/24.5.0-0/bin/conda}"
CONDA_ENV="${CONDA_ENV:-tensorqtl}"
QUEUE="${QUEUE:-medium}"
MEM_GB="${MEM_GB:-24}"   # peak ~10-12 GB for the largest (JECS) files
WALL="${WALL:-4:00}"

LOG_DIR="${GWAS_DIR}/logs"
mkdir -p "${LOG_DIR}"

if [[ $# -gt 0 ]]; then
    TRAITS=("$@")
else
    mapfile -t TRAITS < <(python3 - "${SCRIPTS_DIR}/gwas_catalog.tsv" <<'PY'
import csv, sys
with open(sys.argv[1]) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        print(row["trait_id"])
PY
)
fi

n_sub=0; n_skip=0
for tid in "${TRAITS[@]}"; do
    out="${GWAS_DIR}/${tid}.sumstats.tsv.gz"
    if [[ -s "${out}" && -s "${out}.tbi" ]]; then
        echo "[submit] ${tid}: output exists, skipping"
        n_skip=$((n_skip+1))
        continue
    fi
    JOB_CMD=$(cat <<EOF
eval "\$(${CONDA_BIN} shell.bash hook)"
conda activate ${CONDA_ENV}
module load samtools
python3 "${SCRIPTS_DIR}/38_harmonize_gwas.py" \
    --catalog "${SCRIPTS_DIR}/gwas_catalog.tsv" \
    --raw-dir "${GWAS_DIR}/raw" \
    --out-dir "${GWAS_DIR}" \
    --chain "${CHAIN}" \
    --pgen-dir "${PGEN_DIR}" \
    --ancestries EAS EUR \
    --traits "${tid}"
EOF
)
    bsub -q "${QUEUE}" -n 1 -M "${MEM_GB}" -R "rusage[mem=${MEM_GB}]" \
         -W "${WALL}" -J "harm_${tid}" \
         -o "${LOG_DIR}/38_${tid}.out" -e "${LOG_DIR}/38_${tid}.err" \
         bash -c "${JOB_CMD}"
    n_sub=$((n_sub+1))
done
echo "[submit] ${n_sub} job(s) submitted, ${n_skip} trait(s) already done"
echo "[submit] monitor with: bjobs -w | grep harm_"
echo "[submit] logs: ${LOG_DIR}/38_<trait_id>.{out,err}"
