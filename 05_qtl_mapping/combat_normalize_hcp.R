#!/usr/bin/env Rscript
#
# combat_normalize_hcp.R — ComBat + INT + HCP estimation per ancestry stratum
#
# Batch-correct (ComBat, batch=cohort), quantile-normalize + rank-based INT,
# then estimate HCP latent factors using PicardTools QC metrics as priors.
# HCP factors are output as a tensorQTL covariate table.
#
# Normalization flow (Brain xQTL / Wen et al. 2024 approach):
#   1. Load pooled unnorm expression BED
#   2. log2(TPM + 1)
#   3. ComBat (batch = cohort, par.prior = TRUE)
#   4. Quantile normalization + rank-based INT
#   5. Standardize (center + unit SS) for HCP
#   6. HCP: hcp(Z = QC metrics, Y = INT'd expression, k, lambda1, lambda2, lambda3)
#   7. Extract hidden covariates W (n_samples x k) from HCP result
#
# Usage:
#   Rscript combat_normalize_hcp.R \
#       --expression EUR_pooled_expression.bed \
#       --qc-metrics all_qc_metrics.tsv \
#       --ancestry-map pooled_sample_ancestry_RNAseq.tsv \
#       --ancestry EUR \
#       --k 15 \
#       --output-dir hcp_output/
#
# Dependencies:
#   sva     (Bioconductor) — ComBat
#   Rhcpp   (GitHub: mvaniterson/Rhcpp) — HCP
#   optparse (CRAN) — CLI parsing
#
# Install on seadragon:
#   Rscript -e 'BiocManager::install("sva")'
#   Rscript -e 'remotes::install_github("mvaniterson/Rhcpp")'
#

.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
})

# ---- CLI ----
# NOTE: optparse keeps hyphens in option names, so we must use bracket
# notation opt[["output-dir"]] instead of opt$output-dir (which R parses as
# subtraction: opt$output minus dir).
option_list <- list(
  make_option("--expression", type = "character",
              help = "Pooled unnorm expression BED file"),
  make_option("--qc-metrics", type = "character",
              help = "PicardTools QC metrics TSV (from picard_qc.py)"),
  make_option("--ancestry-map", type = "character",
              help = "Ancestry map TSV: sample_id, assigned_ancestry, cohort"),
  make_option("--ancestry", type = "character",
              help = "Ancestry stratum to process (e.g. EUR, EAS, AFR, HIS)"),
  make_option("--k", type = "integer", default = 15,
              help = "Number of HCP hidden factors (default: 15)"),
  make_option("--lambda1", type = "double", default = 20,
              help = "HCP prior strength (default: 20, from Mostafavi et al.)"),
  make_option("--lambda2", type = "double", default = 1,
              help = "HCP effect regularization (default: 1)"),
  make_option("--lambda3", type = "double", default = 1,
              help = "HCP coefficient regularization (default: 1)"),
  make_option("--output-dir", type = "character", default = ".",
              help = "Output directory"),
  make_option("--nonparametric", action = "store_true", default = FALSE,
              help = "Use non-parametric ComBat (if parametric fails to converge)"),
  make_option("--skip-combat", action = "store_true", default = FALSE,
              help = "Skip ComBat (for single-cohort ancestry strata)"),
  make_option("--min-samples", type = "integer", default = 10,
              help = "Minimum samples required to run a stratum; smaller strata are skipped gracefully (default: 10)"),
  make_option("--skip-hcp", action = "store_true", default = FALSE,
              help = "Skip HCP factor estimation (ComBat+INT expression BED is still written). Use when HCP fails and PEER is run separately."),
  make_option("--qc-cor-threshold", type = "double", default = 0.9,
              help = "Drop QC metrics with |correlation| above this threshold to avoid singular Z'Z in HCP (default: 0.9)")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Use bracket notation for hyphenated option names
opt_expression    <- opt[["expression"]]
opt_qc_metrics    <- opt[["qc-metrics"]]
opt_ancestry_map  <- opt[["ancestry-map"]]
opt_ancestry      <- opt[["ancestry"]]
opt_k             <- opt[["k"]]
opt_lambda1       <- opt[["lambda1"]]
opt_lambda2       <- opt[["lambda2"]]
opt_lambda3       <- opt[["lambda3"]]
opt_output_dir    <- opt[["output-dir"]]
opt_nonparametric <- opt[["nonparametric"]]
opt_skip_combat   <- opt[["skip-combat"]]
opt_min_samples   <- opt[["min-samples"]]
opt_qc_cor_threshold <- opt[["qc-cor-threshold"]]
opt_skip_hcp      <- opt[["skip-hcp"]]

if (is.null(opt_expression) || is.null(opt_qc_metrics) ||
    is.null(opt_ancestry_map) || is.null(opt_ancestry)) {
  stop("--expression, --qc-metrics, --ancestry-map, and --ancestry are required")
}

dir.create(opt_output_dir, showWarnings = FALSE, recursive = TRUE)
ancestry <- opt_ancestry
cat(sprintf("[%s] HCP estimation for ancestry: %s\n", date(), ancestry))

# ---- Helper functions ----

# Quantile normalization: normalize samples (rows) to the same distribution.
# Input mat is samples x genes. Transposes to genes x samples, normalizes
# columns (samples), transposes back. Matches PANTRY normalize_phenotypes.py.
quantile_normalize_rows <- function(mat) {
  # Impute NAs with column (sample) medians so that apply(,2,sort) returns
  # a proper matrix. Without this, sort() drops NAs and columns with
  # different NA counts produce vectors of inconsistent length, causing
  # rowMeans() to fail.
  if (any(is.na(mat))) {
    for (j in seq_len(ncol(mat))) {
      na_idx <- is.na(mat[, j])
      if (any(na_idx)) {
        med <- median(mat[!na_idx, j], na.rm = TRUE)
        if (is.na(med)) med <- 0
        mat[na_idx, j] <- med
      }
    }
  }
  t_mat <- t(mat)  # genes x samples
  sorted_mat <- apply(t_mat, 2, sort)
  mean_dist <- rowMeans(sorted_mat)
  ranked_mat <- apply(t_mat, 2, rank, ties.method = "average")
  result <- matrix(0, nrow = nrow(t_mat), ncol = ncol(t_mat))
  for (j in seq_len(ncol(t_mat))) {
    result[ranked_mat[, j], j] <- mean_dist
  }
  dimnames(result) <- dimnames(t_mat)
  t(result)  # back to samples x genes
}

# Rank-based inverse normal transform per gene (column).
# mat: samples x genes. Matches PANTRY inverse_normal_transform (per gene/row
# in the genes x samples orientation = per column here).
inverse_normal_transform <- function(mat) {
  result <- apply(mat, 2, function(x) {
    r <- rank(x, ties.method = "average", na.last = "keep")
    qnorm(r / (sum(!is.na(x)) + 1))
  })
  dimnames(result) <- dimnames(mat)
  result
}

# Standardize: center to mean 0, scale to unit sum of squares.
# Required by Rhcpp (named differently to avoid conflict with Rhcpp::standardize).
standardize_for_hcp <- function(mat) {
  mat <- as.matrix(mat)
  mat <- mat - matrix(colMeans(mat), nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE)
  ss <- sqrt(colSums(mat^2))
  ss[ss == 0] <- 1
  sweep(mat, 2, ss, "/")
}

# ---- Load data ----
# NOTE: check.names=FALSE is critical — without it, R converts "#chr" to
# "X.chr" because # is not a valid identifier character, breaking all
# downstream column references.
cat(sprintf("[%s] Loading expression BED: %s\n", date(), opt_expression))
# Size colClasses to the actual header (a fixed-length vector longer than the
# column count triggers a spurious "cols != length(data)" warning).
bed_hdr <- strsplit(readLines(opt_expression, n = 1, warn = FALSE), "\t", fixed = TRUE)[[1]]
bed_cc <- rep("numeric", length(bed_hdr))
bed_cc[1:4] <- c("character", "integer", "integer", "character")
expr_bed <- read.delim(opt_expression, sep = "\t", check.names = FALSE,
                       colClasses = bed_cc)
meta_cols <- c("#chr", "start", "end", "phenotype_id")
sample_cols <- setdiff(colnames(expr_bed), meta_cols)
cat(sprintf("  Expression: %d genes x %d samples\n", nrow(expr_bed), length(sample_cols)))

cat(sprintf("[%s] Loading QC metrics: %s\n", date(), opt_qc_metrics))
qc_all <- read.delim(opt_qc_metrics, sep = "\t", row.names = 1, check.names = FALSE)
cat(sprintf("  QC metrics: %d samples x %d metrics\n", nrow(qc_all), ncol(qc_all)))

cat(sprintf("[%s] Loading ancestry map: %s\n", date(), opt_ancestry_map))
anc_map <- read.delim(opt_ancestry_map, sep = "\t", stringsAsFactors = FALSE)
colnames(anc_map) <- tolower(colnames(anc_map))

# ---- Subset to this ancestry stratum ----
stratum_samples <- anc_map[anc_map$assigned_ancestry == ancestry, ]
cat(sprintf("  Ancestry %s: %d samples in map\n", ancestry, nrow(stratum_samples)))

# Verify sample overlap
expr_samples <- intersect(sample_cols, stratum_samples$sample_id)
cat(sprintf("  Samples in expression & ancestry map: %d\n", length(expr_samples)))

if (length(expr_samples) < opt_min_samples) {
  cat(sprintf("  SKIP: too few samples (%d) for ancestry %s (need at least %d). No outputs written for this stratum.\n",
              length(expr_samples), ancestry, opt_min_samples))
  quit(save = "no", status = 0)
}

# Report orphan samples
orphan_expr <- setdiff(sample_cols, stratum_samples$sample_id)
if (length(orphan_expr) > 0) {
  cat(sprintf("  WARN: %d samples in expression but not in ancestry map (excluded)\n",
              length(orphan_expr)))
}
orphan_map <- setdiff(stratum_samples$sample_id, sample_cols)
if (length(orphan_map) > 0) {
  cat(sprintf("  WARN: %d samples in ancestry map but not in expression (excluded)\n",
              length(orphan_map)))
}

# Build expression matrix: samples x genes
expr_samples <- intersect(sample_cols, expr_samples)
expr_data <- t(as.matrix(expr_bed[, expr_samples, drop = FALSE]))
rownames(expr_data) <- expr_samples
colnames(expr_data) <- expr_bed$phenotype_id
cat(sprintf("  Expression matrix: %d samples x %d genes\n",
            nrow(expr_data), ncol(expr_data)))

# ---- Filter genes: remove those with >50% zeros or no variance ----
zero_frac <- colMeans(expr_data == 0, na.rm = TRUE)
no_var <- apply(expr_data, 2, function(x) length(unique(x[!is.na(x)])) <= 1)
keep_genes <- (zero_frac <= 0.5) & !no_var
cat(sprintf("  Gene filtering: %d -> %d (removed %d with >50%% zeros, %d no-variance)\n",
            ncol(expr_data), sum(keep_genes),
            sum(zero_frac > 0.5), sum(no_var & zero_frac <= 0.5)))
expr_data <- expr_data[, keep_genes, drop = FALSE]

# ---- log2(TPM + 1) ----
cat(sprintf("[%s] log2(TPM+1) transform\n", date()))
expr_log <- log2(expr_data + 1)

# ---- ComBat (batch = cohort) ----
cohort_labels <- stratum_samples$cohort[match(expr_samples, stratum_samples$sample_id)]
names(cohort_labels) <- expr_samples

batch_table <- table(cohort_labels)
cat(sprintf("[%s] Batch (cohort) structure:\n", date()))
for (b in names(batch_table)) {
  cat(sprintf("    %s: %d samples\n", b, batch_table[b]))
}

if (opt_skip_combat || length(unique(cohort_labels)) < 2) {
  cat(sprintf("[%s] Skipping ComBat (single cohort or --skip-combat)\n", date()))
  expr_combat <- expr_log
} else {
  # Merge single-sample batches into "other"
  single_batches <- names(batch_table[batch_table < 2])
  if (length(single_batches) > 0) {
    cat(sprintf("  WARN: cohorts with <2 samples: %s\n",
                paste(single_batches, collapse = ", ")))
    cat("  These will be merged into a single 'other' batch for ComBat\n")
    cohort_labels[cohort_labels %in% single_batches] <- "other"
  }

  cat(sprintf("[%s] Running ComBat (batch = cohort, par.prior = %s)\n",
              date(), ifelse(opt_nonparametric, "FALSE", "TRUE")))
  suppressPackageStartupMessages(library(sva))

  mod <- model.matrix(~1, data = data.frame(row.names = expr_samples))

  expr_combat <- tryCatch({
    ComBat(dat = t(expr_log),  # ComBat expects genes x samples
           batch = as.factor(cohort_labels),
           mod = mod,
           par.prior = !opt_nonparametric,
           prior.plots = FALSE)
  }, error = function(e) {
    cat(sprintf("  ComBat (parametric) failed: %s\n", e$message))
    cat("  Retrying with non-parametric prior...\n")
    ComBat(dat = t(expr_log),
           batch = as.factor(cohort_labels),
           mod = mod,
           par.prior = FALSE,
           prior.plots = FALSE)
  })
  expr_combat <- t(expr_combat)  # back to samples x genes
  cat(sprintf("  ComBat done: %d samples x %d genes\n",
              nrow(expr_combat), ncol(expr_combat)))
}

# ---- Quantile normalization + rank-based INT ----
cat(sprintf("[%s] Quantile normalization + rank-based INT\n", date()))
expr_qn <- quantile_normalize_rows(expr_combat)
expr_int <- inverse_normal_transform(expr_qn)
cat(sprintf("  INT done: mean=%.4f, sd=%.4f (should be ~0, ~1)\n",
            mean(expr_int, na.rm = TRUE), sd(expr_int, na.rm = TRUE)))

# ---- Write ComBat + INT expression BED ----
# Match gene metadata from the original BED by phenotype_id.
gene_idx <- match(colnames(expr_int), expr_bed$phenotype_id)
expr_out <- data.frame(
  chr = expr_bed[["#chr"]][gene_idx],
  start = expr_bed$start[gene_idx],
  end = expr_bed$end[gene_idx],
  phenotype_id = colnames(expr_int),
  t(expr_int),
  check.names = FALSE
)
colnames(expr_out)[1] <- "#chr"
expr_out_path <- file.path(opt_output_dir, sprintf("%s_combat_int_expression.bed", ancestry))
write.table(expr_out, expr_out_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Written: %s (%d genes x %d samples)\n",
            expr_out_path, nrow(expr_out), ncol(expr_out) - 4))

# ---- Skip HCP if requested (--skip-hcp) ----
if (opt_skip_hcp) {
  cat(sprintf("\n[%s] --skip-hcp: skipping HCP estimation. ComBat+INT expression BED written.\n", date()))
  cat(sprintf("  Run peer_factors.py separately to extract latent factors.\n"))
  quit(save = "no", status = 0)
}

# ---- Prepare QC metrics for HCP ----
cat(sprintf("[%s] Preparing QC metrics for HCP\n", date()))
qc_subset <- qc_all[expr_samples, , drop = FALSE]
cat(sprintf("  QC metrics: %d samples x %d metrics\n",
            nrow(qc_subset), ncol(qc_subset)))

# Impute missing QC samples with column median
missing_qc <- setdiff(expr_samples, rownames(qc_all))
if (length(missing_qc) > 0) {
  cat(sprintf("  WARN: %d samples missing from QC metrics, imputing with column median\n",
              length(missing_qc)))
  for (s in missing_qc) {
    qc_subset[s, ] <- NA
  }
  qc_subset <- qc_subset[expr_samples, , drop = FALSE]
}

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
# Z'Z in HCP. Iteratively drops the column that has the most remaining
# high-correlation partners, preserving the more "unique" metrics.
qc_cor_threshold <- opt_qc_cor_threshold
if (ncol(qc_subset) > 1) {
  qc_cor_mat <- cor(as.matrix(qc_subset), use = "pairwise.complete.obs")
  diag(qc_cor_mat) <- 0
  drop_cols <- character(0)
  while (any(abs(qc_cor_mat) > qc_cor_threshold, na.rm = TRUE)) {
    # Find the column with the most high-correlation partners
    high_cor_counts <- colSums(abs(qc_cor_mat) > qc_cor_threshold, na.rm = TRUE)
    worst <- names(which.max(high_cor_counts))
    drop_cols <- c(drop_cols, worst)
    qc_cor_mat <- qc_cor_mat[rownames(qc_cor_mat) != worst,
                              colnames(qc_cor_mat) != worst, drop = FALSE]
    if (ncol(qc_cor_mat) <= 1) break
  }
  if (length(drop_cols) > 0) {
    cat(sprintf("  Removing %d highly correlated QC metrics (|r| > %.2f): %s\n",
                length(drop_cols), qc_cor_threshold,
                paste(drop_cols, collapse = ", ")))
    qc_subset <- qc_subset[, !(colnames(qc_subset) %in% drop_cols), drop = FALSE]
  }
}
cat(sprintf("  QC metrics after pruning: %d\n", ncol(qc_subset)))

# Standardize for HCP (center + unit SS)
qc_std <- standardize_for_hcp(qc_subset)
expr_std <- standardize_for_hcp(expr_int)

cat(sprintf("  QC standardized: %d x %d\n", nrow(qc_std), ncol(qc_std)))
cat(sprintf("  Expression standardized: %d x %d\n", nrow(expr_std), ncol(expr_std)))

# ---- HCP estimation ----
cat(sprintf("[%s] Running HCP (k=%d, lambda1=%.1f, lambda2=%.1f, lambda3=%.1f)\n",
            date(), opt_k, opt_lambda1, opt_lambda2, opt_lambda3))
suppressPackageStartupMessages(library(Rhcpp))

# Rhcpp::hcp signature: hcp(Z, Y, k, lambda1, lambda2, lambda3, ...)
#   Z = known covariates matrix (n_samples x d_metrics) — the PRIORS
#   Y = expression matrix (n_samples x g_genes)
# Data is already standardized, so stand=FALSE and log=FALSE.
#
# Return value is a list with:
#   Res = residual expression (n_samples x n_genes) — cleaned data
#   W   = hidden covariates (n_samples x k) — THIS is what we want as QTL covariates
#   B   = effects of hidden covariates (k x n_genes)
#   Y   = input expression (echoed back)
#   Z   = input known covariates (echoed back)
hcp_result <- tryCatch({
  hcp(Z = qc_std, Y = expr_std, k = opt_k,
      lambda1 = opt_lambda1, lambda2 = opt_lambda2, lambda3 = opt_lambda3,
      iter = 100, stand = FALSE, log = FALSE, fast = FALSE, verbose = TRUE)
}, error = function(e) {
  cat(sprintf("  HCP failed: %s\n", e$message))
  stop(e)
})

# Extract hidden covariates W (n_samples x k)
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

# Warn about zero-variance factors (indicates degenerate HCP solution)
factor_vars <- apply(hcp_factors, 2, var)
zero_var_factors <- sum(factor_vars < 1e-10)
if (zero_var_factors > 0) {
  cat(sprintf("  WARN: %d of %d HCP factors have ~zero variance (degenerate). ",
              zero_var_factors, ncol(hcp_factors)))
  cat("Consider reducing k or lowering lambda1.\n")
}

# ---- Write HCP factors as tensorQTL covariate table ----
# tensorQTL expects: rows = covariates, columns = samples, tab-delimited,
# with a header row of sample IDs.
hcp_t <- t(hcp_factors)
hcp_out <- data.frame(
  covariate = rownames(hcp_t),
  hcp_t,
  check.names = FALSE
)
colnames(hcp_out)[1] <- "covariate"

hcp_out_path <- file.path(opt_output_dir, sprintf("%s_hcp_factors.tsv", ancestry))
write.table(hcp_out, hcp_out_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Written: %s (%d factors x %d samples)\n",
            hcp_out_path, nrow(hcp_out), ncol(hcp_out) - 1))

# ---- Diagnostics ----
cat(sprintf("[%s] Generating diagnostics\n", date()))
diag_path <- file.path(opt_output_dir, sprintf("%s_hcp_diagnostics.pdf", ancestry))
pdf(diag_path, width = 12, height = 10)
par(mfrow = c(2, 2))

# 1. Scree plot: variance explained by each HCP factor
# Handle edge case where some factors may have zero/NaN variance (e.g. with
# strong priors on synthetic data).
n_factors <- ncol(hcp_factors)
hcp_var <- apply(hcp_factors, 2, var, na.rm = TRUE)
hcp_var[!is.finite(hcp_var)] <- 0
total_var <- sum(hcp_var)
if (total_var > 0) {
  hcp_var_pct <- hcp_var / total_var * 100
} else {
  hcp_var_pct <- rep(0, n_factors)
}
plot(seq_len(n_factors), hcp_var_pct, type = "b", pch = 19,
     xlab = "HCP factor", ylab = "Variance explained (%)",
     main = "HCP factor variance (scree)",
     ylim = c(0, max(hcp_var_pct, 1)))  # ensure finite ylim

# 2. Correlation heatmap: HCP factors vs QC metrics
cor_mat <- cor(hcp_factors, as.matrix(qc_subset), use = "pairwise.complete.obs")
cor_mat[!is.finite(cor_mat)] <- 0
if (ncol(cor_mat) > 20) {
  max_abs <- apply(abs(cor_mat), 2, max)
  top_idx <- order(max_abs, decreasing = TRUE)[1:20]
  cor_mat_sub <- cor_mat[, top_idx, drop = FALSE]
} else {
  cor_mat_sub <- cor_mat
}
suppressPackageStartupMessages(library(RColorBrewer))
if (ncol(cor_mat_sub) > 0) {
  image(seq_len(ncol(cor_mat_sub)), seq_len(nrow(cor_mat_sub)),
        t(cor_mat_sub)[ncol(cor_mat_sub):1, nrow(cor_mat_sub):1],
        col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
        xlab = "QC metric", ylab = "HCP factor",
        main = "HCP factor vs QC metric correlation",
        axes = FALSE)
  axis(1, seq_len(ncol(cor_mat_sub)), colnames(cor_mat_sub), las = 2, cex.axis = 0.6)
  axis(2, seq_len(nrow(cor_mat_sub)), rev(rownames(cor_mat_sub)), las = 1, cex.axis = 0.7)
}

# 3. PCA of expression before/after ComBat (colored by cohort)
pca_before <- prcomp(expr_log, scale. = FALSE, center = TRUE)
pca_after <- prcomp(expr_int, scale. = FALSE, center = TRUE)
cohort_colors <- as.numeric(as.factor(cohort_labels[rownames(expr_log)]))

plot(pca_before$x[, 1:2], col = cohort_colors, pch = 19, cex = 0.8,
     xlab = sprintf("PC1 (%.1f%%)", summary(pca_before)$importance[2, 1] * 100),
     ylab = sprintf("PC2 (%.1f%%)", summary(pca_before)$importance[2, 2] * 100),
     main = "PCA before ComBat (log2 TPM)")
legend("topright", legend = levels(as.factor(cohort_labels)),
       col = seq_along(levels(as.factor(cohort_labels))), pch = 19, cex = 0.6)

plot(pca_after$x[, 1:2], col = cohort_colors, pch = 19, cex = 0.8,
     xlab = sprintf("PC1 (%.1f%%)", summary(pca_after)$importance[2, 1] * 100),
     ylab = sprintf("PC2 (%.1f%%)", summary(pca_after)$importance[2, 2] * 100),
     main = "PCA after ComBat + INT")
legend("topright", legend = levels(as.factor(cohort_labels)),
       col = seq_along(levels(as.factor(cohort_labels))), pch = 19, cex = 0.6)

dev.off()
cat(sprintf("  Written: %s\n", diag_path))

cat(sprintf("[%s] Done: ancestry %s\n", date(), ancestry))
cat(sprintf("  HCP factors: %s\n", hcp_out_path))
cat(sprintf("  Expression: %s\n", expr_out_path))
cat(sprintf("  Diagnostics: %s\n", diag_path))
