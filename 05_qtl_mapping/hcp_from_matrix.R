#!/usr/bin/env Rscript
#
# hcp_from_matrix.R — HCP-only estimation from an already-normalized BED
#
# Companion to combat_normalize_hcp.R for the per-modality HCP optimization
# module (25b). The input BED is a FINAL harmonized mapping matrix
# ({ANC}_{MOD}_harmonized.bed: QN + INT + ComBat already applied, array_id
# sample columns), so this script runs ONLY the HCP step of the
# combat_normalize_hcp.R schema — no QN/INT, no outlier removal, no ComBat.
# Phenotype values are therefore byte-identical to the matrix tensorQTL maps,
# which isolates the HCP effect when comparing against the previous round.
#
# Steps (mirrors combat_normalize_hcp.R lines ~398-520):
#   1. Load BED -> samples x phenotypes matrix
#   2. Optional deterministic phenotype subsample (--max-phenotypes; bounds
#      HCP cost for the ~150k-phenotype combined arm. HCP factors live in
#      SAMPLE space, so a fixed random phenotype subset is statistically
#      safe — same convention as PEER/SVD practice on large matrices.)
#   3. Load pooled Picard QC metrics (rnaseq_id rows) -> rename to array_id
#      via --metadata, subset to BED samples
#   4. QC prep: median-impute NAs, drop zero-variance metrics, iteratively
#      drop |r| > --qc-cor-threshold metrics (avoids singular Z'Z)
#   5. Standardize both matrices (center + unit sum of squares)
#   6. Rhcpp::hcp(Z = QC, Y = phenotype matrix, k, lambda1-3)
#   7. Write hidden covariates W as a tensorQTL covariate table
#      (covariate x sample, array_id columns — no harmonization needed)
#
# Usage:
#   Rscript hcp_from_matrix.R \
#       --bed EAS_splicing_harmonized.bed \
#       --qc-metrics all_qc_metrics.tsv \
#       --metadata EAS_metadata.tsv \
#       --k 15 \
#       --output EAS_splicing_hcp_factors.tsv
#
# Dependencies: Rhcpp (GitHub: mvaniterson/Rhcpp), optparse (CRAN)
#

.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
})

# ---- CLI ----
option_list <- list(
  make_option("--bed", type = "character",
              help = "Harmonized phenotype BED (QN+INT+ComBat'd, array_id columns)"),
  make_option("--qc-metrics", type = "character",
              help = "Pooled PicardTools QC metrics TSV (rnaseq_id rows)"),
  make_option("--metadata", type = "character",
              help = "{ANC}_metadata.tsv with rnaseq_id and array_id columns"),
  make_option("--k", type = "integer",
              help = "Number of HCP hidden factors"),
  make_option("--lambda1", type = "double", default = 0.5,
              help = "HCP prior strength (default: 0.5, as in 19_hcp_factors.sh)"),
  make_option("--lambda2", type = "double", default = 1,
              help = "HCP effect regularization (default: 1)"),
  make_option("--lambda3", type = "double", default = 1,
              help = "HCP coefficient regularization (default: 1)"),
  make_option("--qc-cor-threshold", type = "double", default = 0.9,
              help = "Drop QC metrics with |correlation| above this threshold (default: 0.9)"),
  make_option("--max-phenotypes", type = "integer", default = 40000,
              help = "Deterministic phenotype subsample cap for HCP estimation; 0 disables (default: 40000)"),
  make_option("--seed", type = "integer", default = 1,
              help = "Seed for the phenotype subsample (default: 1)"),
  make_option("--output", type = "character",
              help = "Output HCP factors TSV path (covariate x sample)")
)

opt <- parse_args(OptionParser(option_list = option_list))

opt_bed           <- opt[["bed"]]
opt_qc_metrics    <- opt[["qc-metrics"]]
opt_metadata      <- opt[["metadata"]]
opt_k             <- opt[["k"]]
opt_lambda1       <- opt[["lambda1"]]
opt_lambda2       <- opt[["lambda2"]]
opt_lambda3       <- opt[["lambda3"]]
opt_qc_cor_threshold <- opt[["qc-cor-threshold"]]
opt_max_phenotypes <- opt[["max-phenotypes"]]
opt_seed          <- opt[["seed"]]
opt_output        <- opt[["output"]]

if (is.null(opt_bed) || is.null(opt_qc_metrics) || is.null(opt_metadata) ||
    is.null(opt_k) || is.null(opt_output)) {
  stop("--bed, --qc-metrics, --metadata, --k, and --output are required")
}
if (opt_k < 1) {
  stop("--k must be >= 1 (the k=0 case is handled by the driver, which writes a header-only file)")
}

cat(sprintf("[%s] HCP-only estimation from normalized BED\n", date()))
cat(sprintf("  BED: %s\n", opt_bed))
cat(sprintf("  k = %d, lambda1 = %.2f, lambda2 = %.1f, lambda3 = %.1f\n",
            opt_k, opt_lambda1, opt_lambda2, opt_lambda3))

# ---- Helper ----
# Standardize: center to mean 0, scale to unit sum of squares.
# Required by Rhcpp (named differently to avoid conflict with Rhcpp::standardize).
standardize_for_hcp <- function(mat) {
  mat <- as.matrix(mat)
  mat <- mat - matrix(colMeans(mat), nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE)
  ss <- sqrt(colSums(mat^2))
  ss[ss == 0] <- 1
  sweep(mat, 2, ss, "/")
}

# ---- Load BED ----
# check.names=FALSE is critical — without it, R converts "#chr" to "X.chr".
cat(sprintf("[%s] Loading BED: %s\n", date(), opt_bed))
bed_hdr <- strsplit(readLines(opt_bed, n = 1, warn = FALSE), "\t", fixed = TRUE)[[1]]
bed_cc <- rep("numeric", length(bed_hdr))
bed_cc[1:4] <- c("character", "integer", "integer", "character")
bed <- read.delim(opt_bed, sep = "\t", check.names = FALSE, colClasses = bed_cc)
meta_cols <- c("#chr", "start", "end", "phenotype_id")
sample_cols <- setdiff(colnames(bed), meta_cols)
cat(sprintf("  BED: %d phenotypes x %d samples\n", nrow(bed), length(sample_cols)))

# ---- Optional deterministic phenotype subsample (HCP cost cap) ----
if (opt_max_phenotypes > 0 && nrow(bed) > opt_max_phenotypes) {
  set.seed(opt_seed)
  keep_idx <- sort(sample(nrow(bed), opt_max_phenotypes))
  cat(sprintf("  Subsampling phenotypes for HCP estimation: %d -> %d (seed %d)\n",
              nrow(bed), opt_max_phenotypes, opt_seed))
  bed <- bed[keep_idx, , drop = FALSE]
}

# Phenotype matrix: samples x phenotypes
expr_data <- t(as.matrix(bed[, sample_cols, drop = FALSE]))
rownames(expr_data) <- sample_cols
colnames(expr_data) <- bed$phenotype_id
rm(bed)
invisible(gc(verbose = FALSE))

# Drop zero-variance phenotypes (singular HCP guard; rare post-INT)
pheno_var <- apply(expr_data, 2, var, na.rm = TRUE)
n_zero_var <- sum(!is.finite(pheno_var) | pheno_var < 1e-10)
if (n_zero_var > 0) {
  cat(sprintf("  Dropping %d zero-variance phenotypes\n", n_zero_var))
  expr_data <- expr_data[, is.finite(pheno_var) & pheno_var >= 1e-10, drop = FALSE]
}
cat(sprintf("  Phenotype matrix: %d samples x %d phenotypes\n",
            nrow(expr_data), ncol(expr_data)))

# ---- Load QC metrics and map to array_id space ----
cat(sprintf("[%s] Loading QC metrics: %s\n", date(), opt_qc_metrics))
qc_all <- read.delim(opt_qc_metrics, sep = "\t", row.names = 1, check.names = FALSE)
cat(sprintf("  QC metrics: %d samples x %d metrics\n", nrow(qc_all), ncol(qc_all)))

meta <- read.delim(opt_metadata, sep = "\t", stringsAsFactors = FALSE)
if (!all(c("rnaseq_id", "array_id") %in% colnames(meta))) {
  stop("--metadata must have rnaseq_id and array_id columns; got: ",
       paste(colnames(meta), collapse = ", "))
}
id_map <- setNames(meta$array_id, meta$rnaseq_id)
mapped <- id_map[rownames(qc_all)]
rownames(qc_all) <- ifelse(is.na(mapped), rownames(qc_all), as.character(mapped))

# Subset QC to BED samples (array_id space); median-impute any missing
expr_samples <- sample_cols
missing_qc <- setdiff(expr_samples, rownames(qc_all))
if (length(missing_qc) > 0) {
  cat(sprintf("  WARN: %d samples missing from QC metrics, imputing with column median\n",
              length(missing_qc)))
  for (s in missing_qc) {
    qc_all[s, ] <- NA
  }
}
qc_subset <- qc_all[expr_samples, , drop = FALSE]

for (col in colnames(qc_subset)) {
  if (any(is.na(qc_subset[[col]]))) {
    med <- median(qc_subset[[col]], na.rm = TRUE)
    if (is.na(med)) med <- 0
    qc_subset[[col]][is.na(qc_subset[[col]])] <- med
  }
}

# Remove zero-variance QC columns
zero_var_qc <- sapply(qc_subset, function(x) sd(as.numeric(x)) == 0)
if (any(zero_var_qc)) {
  cat(sprintf("  Removing %d zero-variance QC metrics: %s\n",
              sum(zero_var_qc), paste(colnames(qc_subset)[zero_var_qc], collapse = ", ")))
  qc_subset <- qc_subset[, !zero_var_qc, drop = FALSE]
}

# Remove highly correlated QC columns (|r| > threshold) to avoid singular
# Z'Z in HCP. Iteratively drops the column with the most remaining
# high-correlation partners.
if (ncol(qc_subset) > 1) {
  qc_cor_mat <- cor(as.matrix(qc_subset), use = "pairwise.complete.obs")
  diag(qc_cor_mat) <- 0
  drop_cols <- character(0)
  while (any(abs(qc_cor_mat) > opt_qc_cor_threshold, na.rm = TRUE)) {
    high_cor_counts <- colSums(abs(qc_cor_mat) > opt_qc_cor_threshold, na.rm = TRUE)
    worst <- names(which.max(high_cor_counts))
    drop_cols <- c(drop_cols, worst)
    qc_cor_mat <- qc_cor_mat[rownames(qc_cor_mat) != worst,
                              colnames(qc_cor_mat) != worst, drop = FALSE]
    if (ncol(qc_cor_mat) <= 1) break
  }
  if (length(drop_cols) > 0) {
    cat(sprintf("  Removing %d highly correlated QC metrics (|r| > %.2f): %s\n",
                length(drop_cols), opt_qc_cor_threshold,
                paste(drop_cols, collapse = ", ")))
    qc_subset <- qc_subset[, !(colnames(qc_subset) %in% drop_cols), drop = FALSE]
  }
}
cat(sprintf("  QC metrics after pruning: %d\n", ncol(qc_subset)))

# ---- Standardize + HCP ----
qc_std <- standardize_for_hcp(qc_subset)
expr_std <- standardize_for_hcp(expr_data)
cat(sprintf("  QC standardized: %d x %d\n", nrow(qc_std), ncol(qc_std)))
cat(sprintf("  Phenotypes standardized: %d x %d\n", nrow(expr_std), ncol(expr_std)))

cat(sprintf("[%s] Running HCP (k=%d)\n", date(), opt_k))
suppressPackageStartupMessages(library(Rhcpp))

# Rhcpp::hcp: Z = known covariates (n_samples x d_metrics), Y = phenotype
# matrix (n_samples x g_phenotypes). Data already standardized -> stand=FALSE.
# W (n_samples x k) = hidden covariates = QTL covariates.
hcp_result <- tryCatch({
  hcp(Z = qc_std, Y = expr_std, k = opt_k,
      lambda1 = opt_lambda1, lambda2 = opt_lambda2, lambda3 = opt_lambda3,
      iter = 100, stand = FALSE, log = FALSE, fast = FALSE, verbose = TRUE)
}, error = function(e) {
  cat(sprintf("  HCP failed: %s\n", e$message))
  stop(e)
})

if (is.list(hcp_result) && "W" %in% names(hcp_result)) {
  hcp_factors <- hcp_result$W
} else {
  stop("HCP result does not contain 'W' (hidden covariates). Got: ",
       paste(names(hcp_result), collapse = ", "))
}

cat(sprintf("  HCP factors: %d samples x %d factors\n",
            nrow(hcp_factors), ncol(hcp_factors)))
colnames(hcp_factors) <- sprintf("HCP_%d", seq_len(ncol(hcp_factors)))
rownames(hcp_factors) <- expr_samples

factor_vars <- apply(hcp_factors, 2, var)
zero_var_factors <- sum(factor_vars < 1e-10)
if (zero_var_factors > 0) {
  cat(sprintf("  WARN: %d of %d HCP factors have ~zero variance (degenerate). ",
              zero_var_factors, ncol(hcp_factors)))
  cat("Consider reducing k or lowering lambda1.\n")
}

# ---- Write tensorQTL covariate table (covariate x sample) ----
hcp_t <- t(hcp_factors)
hcp_out <- data.frame(
  covariate = rownames(hcp_t),
  hcp_t,
  check.names = FALSE
)
colnames(hcp_out)[1] <- "covariate"

write.table(hcp_out, opt_output, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("[%s] Written: %s (%d factors x %d samples)\n",
            date(), opt_output, nrow(hcp_out), ncol(hcp_out) - 1))
