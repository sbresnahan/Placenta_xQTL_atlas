#!/usr/bin/env Rscript
#
# hcp_diagnostic2.R — Test HCP with orthogonalized Z (PCA of QC metrics)
#
# Orthogonalizing Z via PCA makes Z'Z diagonal (condition number 1.0),
# eliminating the numerical singularity while preserving the same column
# space (same information). Also tests finer lambda1 grid around the
# transition point (0.6–0.95).
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
cat("HCP Diagnostic 2 — Orthogonalized Z\n")
cat("========================================\n")

# ---- Load data (same as diagnostic 1) ----
cat("\n[1] Loading data...\n")
bed_hdr <- strsplit(readLines(opt[["expression"]], n = 1, warn = FALSE), "\t", fixed = TRUE)[[1]]
bed_cc <- rep("numeric", length(bed_hdr))
bed_cc[1:4] <- c("character", "integer", "integer", "character")
expr_bed <- read.delim(opt[["expression"]], sep = "\t", check.names = FALSE, colClasses = bed_cc)
meta_cols <- c("#chr", "start", "end", "phenotype_id")
sample_cols <- setdiff(colnames(expr_bed), meta_cols)

qc_all <- read.delim(opt[["qc-metrics"]], sep = "\t", row.names = 1, check.names = FALSE)
anc_map <- read.delim(opt[["ancestry-map"]], sep = "\t", stringsAsFactors = FALSE)
colnames(anc_map) <- tolower(colnames(anc_map))
stratum_samples <- anc_map[anc_map$assigned_ancestry == opt[["ancestry"]], ]
expr_samples <- intersect(sample_cols, stratum_samples$sample_id)

expr_data <- t(as.matrix(expr_bed[, expr_samples, drop = FALSE]))
qc_subset <- qc_all[expr_samples, , drop = FALSE]

for (col in colnames(qc_subset)) {
  if (any(is.na(qc_subset[[col]]))) {
    med <- median(qc_subset[[col]], na.rm = TRUE)
    if (is.na(med)) med <- 0
    qc_subset[[col]][is.na(qc_subset[[col]])] <- med
  }
}

zero_var_qc <- sapply(qc_subset, function(x) sd(as.numeric(x)) == 0)
if (any(zero_var_qc)) qc_subset <- qc_subset[, !zero_var_qc, drop = FALSE]

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
  if (length(drop_cols) > 0)
    qc_subset <- qc_subset[, !(colnames(qc_subset) %in% drop_cols), drop = FALSE]
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

cat(sprintf("  Original Z'Z condition number: %.1f\n", 1/rcond(t(qc_std) %*% qc_std)))

# ---- Orthogonalize Z via PCA ----
cat("\n[2] Orthogonalizing Z via PCA...\n")
# PCA: Z_orth = Z * V, where V are eigenvectors of Z'Z
# Then Z_orth'Z_orth is diagonal (eigenvalues)
ZtZ <- t(qc_std) %*% qc_std
eig <- eigen(ZtZ, symmetric = TRUE)
# Keep all components (preserve full column space)
V <- eig$vectors
qc_orth <- qc_std %*% V  # n x d, orthogonal columns

# Re-standardize to unit SS (PCA changes the scale)
qc_orth_std <- standardize_for_hcp(as.data.frame(qc_orth))

# Verify
ZtZ_orth <- t(qc_orth_std) %*% qc_orth_std
cat(sprintf("  Orthogonalized Z'Z condition number: %.6f\n", 1/rcond(ZtZ_orth)))
cat(sprintf("  Z'Z diagonal? (off-diagonal max: %.2e)\n",
            max(abs(ZtZ_orth - diag(diag(ZtZ_orth))))))
cat(sprintf("  Z dimensions: %d x %d (unchanged)\n", nrow(qc_orth_std), ncol(qc_orth_std)))

# ---- Also find Rhcpp source ----
cat("\n[3] Looking for Rhcpp source files...\n")
pkg_path <- system.file(package = "Rhcpp")
cat(sprintf("  Package path: %s\n", pkg_path))
# Check for compiled code
src_dir <- file.path(pkg_path, "src")
if (dir.exists(src_dir)) {
  src_files <- list.files(src_dir, recursive = TRUE, full.names = TRUE)
  cat(sprintf("  Source files in src/:\n"))
  for (f in src_files) cat(sprintf("    %s (%d bytes)\n", f, file.info(f)$size))
  # Look for hcp-related code
  for (f in src_files) {
    if (grepl("\\.(cpp|c|h|hpp)$", f)) {
      lines <- readLines(f, warn = FALSE)
      hcp_lines <- grep("hcp|inv|solve|singular", lines, ignore.case = TRUE, value = TRUE)
      if (length(hcp_lines) > 0) {
        cat(sprintf("\n  --- %s (hcp/inv/solve lines) ---\n", basename(f)))
        for (l in head(hcp_lines, 30)) cat(sprintf("    %s\n", l))
      }
    }
  }
} else {
  cat("  No src/ directory found.\n")
  # Check if it's a source package
  desc <- readLines(file.path(pkg_path, "DESCRIPTION"), warn = FALSE)
  cat(sprintf("  DESCRIPTION:\n    %s\n", paste(desc, collapse = "\n    ")))
}

# ---- Parameter sweep with BOTH original and orthogonalized Z ----
cat("\n[4] Parameter sweep (original Z vs orthogonalized Z):\n")
suppressPackageStartupMessages(library(Rhcpp))

param_grid <- list(
  list(k = 15, lambda1 = 20, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 5, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 2, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 1, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 0.95, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 0.9, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 0.8, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 0.7, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 0.6, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 0.5, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 0.3, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 0.1, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 0, lambda2 = 1, lambda3 = 1),
  list(k = 15, lambda1 = 0, lambda2 = 1, lambda3 = 0.1),
  list(k = 15, lambda1 = 0, lambda2 = 0.1, lambda3 = 0.1),
  list(k = 5, lambda1 = 1, lambda2 = 1, lambda3 = 1),
  list(k = 5, lambda1 = 0.5, lambda2 = 1, lambda3 = 1),
  list(k = 5, lambda1 = 0, lambda2 = 1, lambda3 = 0.1)
)

cat(sprintf("  %-6s %-8s %-8s %-8s | %-26s | %-26s\n",
            "k", "lambda1", "lambda2", "lambda3",
            "Original Z", "Orthogonalized Z"))
cat(paste(rep("-", 100), collapse = ""), "\n")

for (params in param_grid) {
  k <- params$k; l1 <- params$lambda1; l2 <- params$lambda2; l3 <- params$lambda3

  # Original Z
  res_orig <- tryCatch({
    fit <- hcp(Z = qc_std, Y = expr_std, k = k,
               lambda1 = l1, lambda2 = l2, lambda3 = l3,
               iter = 100, stand = FALSE, log = FALSE, verbose = FALSE)
    W <- fit$W
    wvars <- apply(W, 2, var)
    n_zero <- sum(wvars < 1e-10)
    sprintf("OK zero=%dd var=[%.1e,%.1e]", n_zero, min(wvars), max(wvars))
  }, error = function(e) {
    sprintf("ERROR: %s", substr(e$message, 1, 22))
  })

  # Orthogonalized Z
  res_orth <- tryCatch({
    fit <- hcp(Z = qc_orth_std, Y = expr_std, k = k,
               lambda1 = l1, lambda2 = l2, lambda3 = l3,
               iter = 100, stand = FALSE, log = FALSE, verbose = FALSE)
    W <- fit$W
    wvars <- apply(W, 2, var)
    n_zero <- sum(wvars < 1e-10)
    sprintf("OK zero=%dd var=[%.1e,%.1e]", n_zero, min(wvars), max(wvars))
  }, error = function(e) {
    sprintf("ERROR: %s", substr(e$message, 1, 22))
  })

  cat(sprintf("  %-6d %-8.2f %-8.2f %-8.2f | %-26s | %-26s\n",
              k, l1, l2, l3, res_orig, res_orth))
}

cat("\n[5] Done.\n")
