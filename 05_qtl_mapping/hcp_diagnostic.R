#!/usr/bin/env Rscript
#
# hcp_diagnostic.R — Diagnose why Rhcpp::hcp produces zero-variance factors
#
# Run interactively (not via LSF):
#   singularity exec --bind /rsrch5 --bind /rsrch9 \
#     /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif Rscript hcp_diagnostic.R \
#     --expression /rsrch9/.../EAS_combat_int_expression.bed \
#     --qc-metrics /rsrch9/.../all_qc_metrics.tsv \
#     --ancestry-map /rsrch9/.../pooled_sample_ancestry_RNAseq.tsv \
#     --ancestry EAS
#

.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--expression", type = "character"),
  make_option("--qc-metrics", type = "character"),
  make_option("--ancestry-map", type = "character"),
  make_option("--ancestry", type = "character", default = "EAS"),
  make_option("--qc-cor-threshold", type = "double", default = 0.9)
)
opt <- parse_args(OptionParser(option_list = option_list))

.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))

cat("\n========================================\n")
cat("HCP Diagnostic\n")
cat("========================================\n")

# ---- Load data (same as combat_normalize_hcp.R) ----
cat("\n[1] Loading expression BED...\n")
bed_hdr <- strsplit(readLines(opt[["expression"]], n = 1, warn = FALSE), "\t", fixed = TRUE)[[1]]
bed_cc <- rep("numeric", length(bed_hdr))
bed_cc[1:4] <- c("character", "integer", "integer", "character")
expr_bed <- read.delim(opt[["expression"]], sep = "\t", check.names = FALSE, colClasses = bed_cc)
meta_cols <- c("#chr", "start", "end", "phenotype_id")
sample_cols <- setdiff(colnames(expr_bed), meta_cols)
cat(sprintf("  Expression: %d genes x %d samples\n", nrow(expr_bed), length(sample_cols)))

cat("\n[2] Loading QC metrics...\n")
qc_all <- read.delim(opt[["qc-metrics"]], sep = "\t", row.names = 1, check.names = FALSE)
cat(sprintf("  QC metrics: %d samples x %d metrics\n", nrow(qc_all), ncol(qc_all)))

cat("\n[3] Loading ancestry map...\n")
anc_map <- read.delim(opt[["ancestry-map"]], sep = "\t", stringsAsFactors = FALSE)
colnames(anc_map) <- tolower(colnames(anc_map))
stratum_samples <- anc_map[anc_map$assigned_ancestry == opt[["ancestry"]], ]
expr_samples <- intersect(sample_cols, stratum_samples$sample_id)
cat(sprintf("  Ancestry %s: %d samples in expression\n", opt[["ancestry"]], length(expr_samples)))

# ---- Build matrices ----
expr_data <- t(as.matrix(expr_bed[, expr_samples, drop = FALSE]))
rownames(expr_data) <- expr_samples
qc_subset <- qc_all[expr_samples, , drop = FALSE]

# Impute QC NAs
for (col in colnames(qc_subset)) {
  if (any(is.na(qc_subset[[col]]))) {
    med <- median(qc_subset[[col]], na.rm = TRUE)
    if (is.na(med)) med <- 0
    qc_subset[[col]][is.na(qc_subset[[col]])] <- med
  }
}

# Remove zero-variance QC
zero_var_qc <- sapply(qc_subset, function(x) sd(as.numeric(x)) == 0)
if (any(zero_var_qc)) {
  cat(sprintf("  Removing %d zero-variance QC metrics: %s\n",
              sum(zero_var_qc), paste(colnames(qc_subset)[zero_var_qc], collapse = ", ")))
  qc_subset <- qc_subset[, !zero_var_qc, drop = FALSE]
}

# Correlation pruning
qc_cor_threshold <- opt[["qc-cor-threshold"]]
if (ncol(qc_subset) > 1) {
  qc_cor_mat <- cor(as.matrix(qc_subset), use = "pairwise.complete.obs")
  diag(qc_cor_mat) <- 0
  drop_cols <- character(0)
  while (any(abs(qc_cor_mat) > qc_cor_threshold, na.rm = TRUE)) {
    high_cor_counts <- colSums(abs(qc_cor_mat) > qc_cor_threshold, na.rm = TRUE)
    worst <- names(which.max(high_cor_counts))
    drop_cols <- c(drop_cols, worst)
    qc_cor_mat <- qc_cor_mat[rownames(qc_cor_mat) != worst,
                              colnames(qc_cor_mat) != worst, drop = FALSE]
    if (ncol(qc_cor_mat) <= 1) break
  }
  if (length(drop_cols) > 0) {
    cat(sprintf("  Removing %d highly correlated QC metrics (|r| > %.2f): %s\n",
                length(drop_cols), qc_cor_threshold, paste(drop_cols, collapse = ", ")))
    qc_subset <- qc_subset[, !(colnames(qc_subset) %in% drop_cols), drop = FALSE]
  }
}
cat(sprintf("  QC metrics after pruning: %d\n", ncol(qc_subset)))

# ---- Standardize ----
standardize_for_hcp <- function(mat) {
  mat <- as.matrix(mat)
  mat <- mat - matrix(colMeans(mat), nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE)
  ss <- sqrt(colSums(mat^2))
  ss[ss == 0] <- 1
  sweep(mat, 2, ss, "/")
}

qc_std <- standardize_for_hcp(qc_subset)
expr_std <- standardize_for_hcp(expr_data)

# ---- Diagnostics ----
cat("\n[4] Expression matrix properties:\n")
cat(sprintf("  Dimensions: %d samples x %d genes\n", nrow(expr_std), ncol(expr_std)))
cat(sprintf("  ||Y||^2 = %.1f (should equal n_genes = %d)\n", sum(expr_std^2), ncol(expr_std)))
cat(sprintf("  Column SS range: [%.4f, %.4f] (should all be 1.0)\n",
            min(colSums(expr_std^2)), max(colSums(expr_std^2))))
cat(sprintf("  Any NaN/Inf: %s\n", any(!is.finite(expr_std))))

# SVD of expression
cat("\n[5] SVD of expression matrix:\n")
sv <- svd(expr_std, nv = 0)
cat(sprintf("  Top 20 singular values:\n"))
for (i in 1:min(20, length(sv$d))) {
  cat(sprintf("    d[%d] = %.4f\n", i, sv$d[i]))
}
cat(sprintf("  Rank (d > 1e-10): %d\n", sum(sv$d > 1e-10)))
cat(sprintf("  Condition number (d[1]/d[min]): %.1f\n",
            sv$d[1] / sv$d[min(length(sv$d), sum(sv$d > 1e-10))]))

cat("\n[6] QC metrics (Z) properties:\n")
cat(sprintf("  Dimensions: %d samples x %d metrics\n", nrow(qc_std), ncol(qc_std)))
cat(sprintf("  ||Z||^2 = %.1f (should equal n_metrics = %d)\n", sum(qc_std^2), ncol(qc_std)))
ZtZ <- t(qc_std) %*% qc_std
cat(sprintf("  Z'Z condition number: %.1f\n", rcond(ZtZ)))
cat(sprintf("  Z'Z eigenvalues:\n"))
eig_ZtZ <- eigen(ZtZ, only.values = TRUE)$values
for (i in 1:min(20, length(eig_ZtZ))) {
  cat(sprintf("    lambda[%d] = %.6f\n", i, eig_ZtZ[i]))
}

# ---- Check Rhcpp source ----
cat("\n[7] Rhcpp package info:\n")
cat(sprintf("  Rhcpp version: %s\n", as.character(packageVersion("Rhcpp"))))
cat(sprintf("  Source files:\n"))
pkg_path <- system.file(package = "Rhcpp")
cat(sprintf("  Package path: %s\n", pkg_path))
R_files <- list.files(pkg_path, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
for (f in R_files) {
  cat(sprintf("    %s\n", f))
}

# Try to see the hcp function body
cat("\n[8] Rhcpp::hcp function:\n")
tryCatch({
  cat("  Methods:\n")
  print(methods(hcp))
  cat("\n  Function body (first 50 lines):\n")
  body_text <- deparse(body(hcp))
  cat(paste(head(body_text, 50), collapse = "\n"))
  cat("\n")
}, error = function(e) {
  cat(sprintf("  Could not inspect hcp body: %s\n", e$message))
})

# Check for Rcpp source
cpp_files <- list.files(file.path(pkg_path, "src"), pattern = "\\.cpp$",
                         recursive = TRUE, full.names = TRUE)
if (length(cpp_files) > 0) {
  cat(sprintf("\n  C++ source files:\n"))
  for (f in cpp_files) {
    cat(sprintf("    %s\n", f))
    # Look for initialization code
    lines <- readLines(f, warn = FALSE)
    init_lines <- grep("init|W\\s*=|X\\s*=|svd|SVD|random|start", lines,
                       ignore.case = TRUE, value = TRUE)
    if (length(init_lines) > 0) {
      cat(sprintf("    Lines matching init/W=/svd/random:\n"))
      for (l in head(init_lines, 20)) {
        cat(sprintf("      %s\n", l))
      }
    }
  }
}

# ---- Try HCP with multiple parameter combinations ----
cat("\n[9] Testing HCP with different parameters:\n")
suppressPackageStartupMessages(library(Rhcpp))

param_grid <- list(
  list(k = 15, lambda1 = 20, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 5, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 1, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 0.5, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 0.1, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 0.01, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 0, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 0, lambda2 = 1, lambda3 = 0.1),
  list(k = 15, lambda1 = 0, lambda2 = 0.1, lambda3 = 0.1),
  list(k = 15, lambda1 = 0, lambda2 = 0.01, lambda3 = 0.01),
  list(k = 5, lambda1 = 20, lambda2 = 1, lambda3 = 1),
  list(k = 5, lambda1 = 1, lambda2 = 1, lambda3 = 1),
  list(k = 5, lambda1 = 0, lambda2 = 1, lambda3 = 0.1),
  list(k = 1, lambda1 = 20, lambda2 = 1, lambda3 = 1),
  list(k = 1, lambda1 = 1, lambda2 = 1, lambda3 = 1),
  list(k = 1, lambda1 = 0, lambda2 = 1, lambda3 = 0.1)
)

cat(sprintf("  %-6s %-8s %-8s %-8s %-8s %-10s %-15s %s\n",
            "k", "lambda1", "lambda2", "lambda3", "iter", "status", "zero_var_factors", "W_var_range"))
cat(paste(rep("-", 90), collapse = ""), "\n")

for (params in param_grid) {
  k <- params$k; l1 <- params$lambda1; l2 <- params$lambda2; l3 <- params$lambda3
  result <- tryCatch({
    fit <- hcp(Z = qc_std, Y = expr_std, k = k,
               lambda1 = l1, lambda2 = l2, lambda3 = l3,
               iter = 100, stand = FALSE, log = FALSE, verbose = FALSE)
    W <- fit$W
    wvars <- apply(W, 2, var)
    n_zero <- sum(wvars < 1e-10)
    iter_used <- attr(fit, "niter")
    if (is.null(iter_used)) iter_used <- NA
    sprintf("  %-6d %-8.2f %-8.2f %-8.2f %-8s %-10s %-15d [%.4e, %.4e]\n",
            k, l1, l2, l3, as.character(iter_used), "OK", n_zero,
            min(wvars), max(wvars))
  }, error = function(e) {
    sprintf("  %-6d %-8.2f %-8.2f %-8.2f %-8s %-10s %-15s %s\n",
            k, l1, l2, l3, "-", "ERROR", "-", substr(e$message, 1, 40))
  })
  cat(result)
}

cat("\n[10] Done.\n")
