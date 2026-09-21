#!/bin/bash
# =============================================================================
# 21_install_tensorqtl.sh — One-time setup: conda env for tensorQTL
#                           + R qvalue for Storey q-values (singularity R)
# =============================================================================
# Run interactively or as a short LSF job:
#   bash 21_install_tensorqtl.sh
# Or:
#   bsub -q short -n 1 -M 8 -R "rusage[mem=8]" -W 0:30 bash 21_install_tensorqtl.sh
#
# Two independent pieces:
#   1. The tensorqtl conda env (python-only — UNCHANGED from the working
#      original; no conda R, no rpy2. The conda-R route was abandoned after
#      hitting libstdc++/GLIBCXX conflicts with the system /lib64 on
#      seadragon.)
#   2. The R 'qvalue' package (Storey q-values, GTEx convention), installed
#      into the singularity R library (R_LIBS_USER) — the same R container
#      scripts 19/20 already use. 27_run_tensorqtl.py calls it through
#      compute_qvalues.R (file-based bridge, no rpy2).
#
# NOTE: if you previously ran the rpy2 version of this script, your conda
# env contains r-base/rpy2/bioconductor-qvalue that are now UNUSED. They are
# harmless — the 'rfunc cannot be imported' warning at 'import tensorqtl'
# is cosmetic and does not affect mapping. (Optional cleanup:
# conda remove -n tensorqtl r-base rpy2 bioconductor-qvalue.)
# =============================================================================

set -euo pipefail

SCRIPT_HOME="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- Conda init ----
source /etc/profile.d/modules.sh
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

# ---- Create env ----
ENV_NAME="tensorqtl"

if conda env list | grep -q "^${ENV_NAME} "; then
    echo "Conda env '${ENV_NAME}' already exists. Activating."
    conda activate "${ENV_NAME}"
else
    echo "Creating conda env '${ENV_NAME}'..."
    conda create -n "${ENV_NAME}" python=3.10 -y
    conda activate "${ENV_NAME}"
fi

# ---- Install packages ----
echo "Installing tensorQTL and dependencies..."
pip install tensorqtl pandas numpy scipy pyarrow matplotlib

# pandas_plink — required by tensorqtl.genotypeio (imported eagerly at
# 'import tensorqtl'; not always pulled in as a dependency)
pip install pandas_plink

# plinkio — needed by tensorqtl for pgen I/O
pip install plinkio || {
    echo "pip install plinkio failed; trying from GitHub..."
    pip install git+https://github.com/broadinstitute/plinkio.git
}

# ---- Verify python side ----
echo ""
echo "Verifying installation..."
python - <<'EOF'
import importlib.metadata as im
import tensorqtl
print(f'  tensorqtl: {tensorqtl.__version__}')
for pkg in ['pandas', 'numpy', 'scipy', 'pyarrow', 'matplotlib']:
    try:
        print(f'  {pkg}: {im.version(pkg)}')
    except im.PackageNotFoundError:
        print(f'  {pkg}: NOT FOUND')
for mod in ['plinkio', 'pandas_plink']:
    try:
        __import__(mod); print(f'  {mod}: OK')
    except ImportError:
        print(f'  {mod}: NOT FOUND')
EOF

# ---- Install R qvalue into the singularity R library ----
echo ""
echo "Installing R qvalue package into the singularity R library..."
export R_LIBS_USER="${R_LIBS_USER:-/rsrch5/home/epi/stbresnahan/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1}"
SING_R="singularity exec --bind /rsrch5 --bind /rsrch9 /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif Rscript"

$SING_R -e '
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
if (!requireNamespace("qvalue", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("qvalue", lib = Sys.getenv("R_LIBS_USER"), update = FALSE, ask = FALSE)
}
stopifnot(requireNamespace("qvalue", quietly = TRUE))
cat("  qvalue:", as.character(packageVersion("qvalue")), "installed in", Sys.getenv("R_LIBS_USER"), "\n")
'

# ---- Smoke test: the Storey bridge end-to-end through singularity R ----
echo ""
echo "Smoke test: compute_qvalues.R bridge (Storey, singularity R qvalue)..."
TMPD=$(mktemp -d)
python -c "
import numpy as np, pandas as pd
rng = np.random.default_rng(1)
p = np.concatenate([rng.uniform(0, 1, 900), rng.beta(0.5, 50, 100)])
pd.DataFrame({'pval_beta': p}).to_csv('${TMPD}/pvals.tsv', sep='\t', index=False)
"
if $SING_R "${SCRIPT_HOME}/compute_qvalues.R" "${TMPD}/pvals.tsv" "${TMPD}/qvals.tsv"; then
    python -c "
import pandas as pd
df = pd.read_csv('${TMPD}/qvals.tsv', sep='\t')
assert len(df) == 1000 and df['qval'].between(0, 1).all()
print('  Storey bridge: OK (1000 q-values, all in [0,1])')
"
else
    echo "ERROR: compute_qvalues.R bridge FAILED — qvalue is not usable in the singularity R."
    exit 1
fi
rm -rf "$TMPD"

echo ""
echo "Done. Activate with: conda activate ${ENV_NAME}"
