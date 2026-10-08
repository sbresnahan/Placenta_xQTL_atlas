#!/bin/bash
# =============================================================================
# 36_install_coloc_env.sh — one-time environment setup for module 06
# =============================================================================
# Objective 1.6 (GWAS colocalization and isoTWAS) needs R packages that are not
# in the lab singularity library (coloc >= 5.2 for coloc.susie, susieR,
# colocboost, plink2R) plus the FUSION TWAS scripts. This script installs them
# into the same R 4.3.1 library used across the pipeline (via bin/Rscript_sif)
# and clones FUSION.
#
# Run ONCE on a login node (compute nodes have no internet access):
#
#   bash 36_install_coloc_env.sh
#
# Optional overrides (environment):
#   FUSION_DIR  — where to clone fusion_twas
#                 (default: /rsrch5/home/epi/stbresnahan/bhattacharya_lab/software/fusion_twas)
#   SKIP_R=1    — skip R package installation (FUSION clone only)
#
# Notes:
#   - colocboost is GitHub-only (StatFunGen/colocboost; Cao et al. 2025,
#     medRxiv 2025.04.17.25326042) and is installed with remotes::install_github.
#   - plink2R (gabraham/plink2R) is required by FUSION.assoc_test.R to read
#     LD reference panels in plink format.
#   - susieR >= 0.12.35 and coloc >= 5.2.1 are required for coloc.susie.
# =============================================================================
set -euo pipefail

SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "$0")" && pwd)}"
REPO_ROOT="$(dirname "$SCRIPTS_DIR")"
RSCRIPT="${RSCRIPT:-${REPO_ROOT}/bin/Rscript_sif}"
FUSION_DIR="${FUSION_DIR:-/rsrch5/home/epi/stbresnahan/bhattacharya_lab/software/fusion_twas}"
SKIP_R="${SKIP_R:-0}"

echo "=== 36_install_coloc_env.sh ==="
echo "  RSCRIPT:    $RSCRIPT"
echo "  FUSION_DIR: $FUSION_DIR"

if [ "$SKIP_R" != "1" ]; then
    echo ""
    echo "== Step 1: R packages (CRAN) =="
    "$RSCRIPT" - <<'EOF'
R_LIB <- "/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"
if (!dir.exists(R_LIB)) stop("R package library does not exist: ", R_LIB)
if (file.access(R_LIB, 2L) != 0L) stop("R package library is not writable: ", R_LIB)
.libPaths(c(R_LIB, .libPaths()))
message("R package library: ", .libPaths()[1])

options(repos = c(CRAN = "https://cloud.r-project.org"))
need_cran <- c("coloc", "susieR", "glmnet", "optparse", "data.table",
               "R.utils", "remotes", "matrixStats", "irlba")
have <- rownames(installed.packages(lib.loc = R_LIB))
for (p in need_cran[!need_cran %in% have]) {
  message("installing ", p, " -> ", R_LIB)
  install.packages(p, lib = R_LIB)
}
# Version gates: coloc.susie requires coloc >= 5.2; susieR >= 0.12.35
stopifnot(packageVersion("coloc") >= "5.2.0")
stopifnot(packageVersion("susieR") >= "0.12.0")
message("CRAN packages OK: coloc ", as.character(packageVersion("coloc")),
        ", susieR ", as.character(packageVersion("susieR")),
        ", glmnet ", as.character(packageVersion("glmnet")))
EOF

    echo ""
    echo "== Step 2: GitHub packages (colocboost, plink2R) =="
    "$RSCRIPT" - <<'EOF'
R_LIB <- "/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"
if (!dir.exists(R_LIB)) stop("R package library does not exist: ", R_LIB)
if (file.access(R_LIB, 2L) != 0L) stop("R package library is not writable: ", R_LIB)
.libPaths(c(R_LIB, .libPaths()))
message("R package library: ", .libPaths()[1])

suppressPackageStartupMessages(library(remotes))
have <- rownames(installed.packages(lib.loc = R_LIB))
if (!"colocboost" %in% have) {
  message("installing StatFunGen/colocboost -> ", R_LIB)
  remotes::install_github("StatFunGen/colocboost", lib = R_LIB,
                          upgrade = "never")
}
if (!"plink2R" %in% have) {
  message("installing gabraham/plink2R (FUSION LD-reference reader) -> ", R_LIB)
  remotes::install_github("gabraham/plink2R", subdir = "plink2R",
                          lib = R_LIB, upgrade = "never")
}
stopifnot("colocboost" %in% rownames(installed.packages(lib.loc = R_LIB)))
stopifnot("plink2R" %in% rownames(installed.packages(lib.loc = R_LIB)))
message("GitHub packages OK")
EOF
fi

echo ""
echo "== Step 3: FUSION TWAS scripts =="
if [ -d "$FUSION_DIR" ] && [ -f "$FUSION_DIR/FUSION.assoc_test.R" ]; then
    echo "  FUSION already present: $FUSION_DIR"
else
    git clone https://github.com/gusevlab/fusion_twas "$FUSION_DIR"
fi

echo ""
echo "== Step 4: smoke tests =="
"$RSCRIPT" - <<'EOF'
R_LIB <- "/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"
if (!dir.exists(R_LIB)) stop("R package library does not exist: ", R_LIB)
.libPaths(c(R_LIB, .libPaths()))

suppressPackageStartupMessages({
  library(coloc); library(susieR); library(glmnet); library(colocboost)
})
# coloc.susie smoke test on the package's built-in simulated data
data(coloc_test_data, package = "coloc")
f <- coloc::coloc.susie(coloc_test_data$D1, coloc_test_data$D2)
stopifnot(is.finite(f$summary["PP.H4.abf"]))
message("coloc.susie smoke test OK (PP.H4 = ",
        round(f$summary["PP.H4.abf"], 3), ")")
EOF

echo ""
echo "Done. Record installed versions in docs/software_environments.md:"
"$RSCRIPT" -e '.libPaths(c("/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1", .libPaths())); cat(sprintf("coloc %s | susieR %s | glmnet %s | colocboost %s\n",
  packageVersion("coloc"), packageVersion("susieR"),
  packageVersion("glmnet"), packageVersion("colocboost")))'
echo "FUSION: $(cd "$FUSION_DIR" && git rev-parse --short HEAD)"
