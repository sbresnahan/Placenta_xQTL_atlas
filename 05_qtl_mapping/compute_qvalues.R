#!/usr/bin/env Rscript
# =============================================================================
# compute_qvalues.R — Storey q-values via the R qvalue package (GTEx convention)
# =============================================================================
# File-based bridge called by 27_run_tensorqtl.py. This replaces tensorQTL's
# calculate_qvalues() rpy2 path, which is fragile on seadragon (conda-R /
# libstdc++ conflicts). Statistically identical: tensorQTL's calculate_qvalues
# is itself a wrapper around qvalue::qvalue with default lambda.
#
# Usage: Rscript compute_qvalues.R <in.tsv> <out.tsv>
#   in.tsv : single column 'pval_beta' (permutation-calibrated p-values)
#   out.tsv: same rows + 'qval' column. Prints pi0 to stdout.
# Exit status 1 (and no output file) on any failure — the caller hard-errors.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript compute_qvalues.R <in.tsv> <out.tsv>")
}
in_file <- args[1]
out_file <- args[2]

ok <- tryCatch({
  suppressPackageStartupMessages(library(qvalue))

  df <- read.delim(in_file)
  if (!("pval_beta" %in% colnames(df))) {
    stop("input TSV must have a 'pval_beta' column")
  }
  p <- as.numeric(df$pval_beta)
  p[!is.finite(p)] <- 1          # NA/NaN -> 1 (conservative; matches BH fallback)
  p <- pmin(pmax(p, 0), 1)       # clamp to [0,1]
  if (length(unique(p)) < 2) {
    stop("fewer than 2 unique p-values — qvalue cannot estimate pi0")
  }

  q <- qvalue(p)                 # default lambda: automatic smoother (GTEx)
  df$qval <- q$qvalues
  cat(sprintf("  R qvalue: pi0 = %.4f (estimated signal fraction 1-pi0 = %.4f)\n",
              q$pi0, 1 - q$pi0))
  write.table(df, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  TRUE
}, error = function(e) {
  cat(sprintf("ERROR in compute_qvalues.R: %s\n", conditionMessage(e)), file = stderr())
  FALSE
})

if (!ok) quit(save = "no", status = 1)
