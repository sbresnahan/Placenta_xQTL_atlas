#!/usr/bin/env Rscript
# =============================================================================
# 48b_fusion_assoc.R — wrapper around FUSION.assoc_test.R
# =============================================================================
# FUSION.assoc_test.R (gusevlab/fusion_twas clone, not part of this repo) loads
# plink2R/optparse at startup but does NOT set .libPaths(), and bash-level
# R_LIBS_* / R_LIBS_USER are not reliably propagated to R on this cluster.
# This wrapper points .libPaths() at the lab R 4.3.1 package library INSIDE R
# (same contract as 41_susie_coloc.R / 44_colocboost.R / 46_isotwas_train.R)
# and then sources FUSION.assoc_test.R.
#
# FUSION.assoc_test.R parses its flags with optparse::parse_args(), which
# reads commandArgs(trailingOnly = TRUE); arguments passed after this wrapper
# on the Rscript command line are therefore seen by FUSION unchanged:
#
#   "$RSCRIPT" 48b_fusion_assoc.R --sumstats X --weights Y ... --out Z
#
# FUSION_DIR must be in the environment (48a_run_fusion.sh exports it; it is
# propagated by 48_fusion_twas.sh via bsub -env).
# =============================================================================

R_LIB <- "/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"
if (!dir.exists(R_LIB)) stop("R package library does not exist: ", R_LIB)
.libPaths(c(R_LIB, .libPaths()))

fusion_dir <- Sys.getenv("FUSION_DIR", unset = "")
if (!nzchar(fusion_dir)) stop("FUSION_DIR is not set in the environment")
fusion_assoc <- file.path(fusion_dir, "FUSION.assoc_test.R")
if (!file.exists(fusion_assoc)) stop("FUSION.assoc_test.R not found: ", fusion_assoc)

# Evaluate in the global environment so FUSION's top-level code (option
# parsing included) runs exactly as if it had been invoked directly.
source(fusion_assoc, local = FALSE)
