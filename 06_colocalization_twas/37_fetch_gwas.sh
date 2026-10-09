#!/bin/bash
# =============================================================================
# 37_fetch_gwas.sh — fetch Table-2 GWAS summary statistics (LOGIN NODE ONLY)
# =============================================================================
# Objective 1.6: downloads the GWAS summary statistics catalogued in
# gwas_catalog.tsv into $OUTPUT_BASE/gwas/raw/. Compute nodes have no
# internet access — run this on a login node, not via bsub.
#
#   bash 37_fetch_gwas.sh                 # all catalog traits
#   TRAITS="egg_bw_fetal_2019" bash 37_fetch_gwas.sh   # subset
#
# JECS (request form) and ProDiGY (T2D Knowledge Portal registration) are
# manual-placement traits: the script prints instructions and verifies files
# you place at $OUTPUT_BASE/gwas/raw/{trait_id}.txt.gz.
#
# Optional overrides: OUTPUT_BASE, GWAS_DIR (default $OUTPUT_BASE/gwas),
#   FORCE=1 to re-download, PYENV (conda env providing pandas; default
#   tensorqtl from module-05 script 21; PYENV=none skips activation).
# =============================================================================
set -eo pipefail

SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "$0")" && pwd)}"
OUTPUT_BASE="${OUTPUT_BASE:-/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY}"
GWAS_DIR="${GWAS_DIR:-${OUTPUT_BASE}/gwas}"
TRAITS="${TRAITS:-}"
FORCE="${FORCE:-0}"

# --- Python environment ------------------------------------------------------
# 37_fetch_gwas.py requires pandas, which the login-node system python3 does
# not provide. Activate the pipeline conda env (module-05 script 21 creates
# 'tensorqtl': python 3.10 + pandas/numpy/scipy). Override with PYENV=<env>;
# PYENV=none skips activation (e.g. when run inside another suitable env).
PYENV="${PYENV:-tensorqtl}"
if [ "$PYENV" != "none" ]; then
    source /etc/profile.d/modules.sh
    eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
    conda activate "$PYENV"
fi
if ! python3 -c "import pandas" 2>/dev/null; then
    echo "ERROR: python3 cannot import pandas." >&2
    echo "  Activate an env that provides pandas, e.g.:" >&2
    echo "    conda activate tensorqtl   (module-05 script 21)" >&2
    echo "  or rerun with PYENV=<env_name>." >&2
    exit 1
fi

mkdir -p "${GWAS_DIR}/raw"

ARGS=(--catalog "${SCRIPTS_DIR}/gwas_catalog.tsv" --out-dir "${GWAS_DIR}/raw")
if [ -n "$TRAITS" ]; then
    # shellcheck disable=SC2206
    ARGS+=(--traits $TRAITS)
fi
if [ "$FORCE" = "1" ]; then
    ARGS+=(--force)
fi

echo "=== 37_fetch_gwas.sh (login node) ==="
python3 "${SCRIPTS_DIR}/37_fetch_gwas.py" "${ARGS[@]}"

echo ""
echo "Next: harmonize to GRCh38 + pooled-pgen variant IDs"
echo "(needs the pandas env above plus bgzip/tabix, e.g. 'module load samtools'):"
echo "  python3 ${SCRIPTS_DIR}/38_harmonize_gwas.py \\"
echo "    --catalog ${SCRIPTS_DIR}/gwas_catalog.tsv \\"
echo "    --raw-dir ${GWAS_DIR}/raw --out-dir ${GWAS_DIR} \\"
echo "    --chain /home/stbresnahan/bhattacharya_lab/data/GenomicReferences/liftover/hg19ToHg38.over.chain \\"
echo "    --pgen-dir /rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/pooled/genotypes"
echo "NOTE: --out-dir is ${GWAS_DIR} itself (not ${GWAS_DIR}/harmonized);"
echo "  steps 40/42/43/45/48 expect ${GWAS_DIR}/{trait_id}.sumstats.tsv.gz."
