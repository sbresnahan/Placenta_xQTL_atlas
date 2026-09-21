#!/usr/bin/env Rscript
#
# combat_normalize_modalities.R — ComBat + INT for non-expression RNA modalities
#
# Generalizes combat_normalize_hcp.R to all RNA phenotype modalities, WITHOUT
# HCP estimation (HCP factors are estimated once on gene-level expression and
# reused as covariates for every modality).
#
# Per ancestry stratum, per modality:
#   1. Load pooled unnorm BED (samples x features)
#   2. Modality-specific pre-transform:
#        - Proportion modalities (bounded [0,1]): logit with epsilon clipping
#             logit(p) = log((p + eps) / (1 - p + eps)),  eps = 1e-4
#          Maps [0,1] -> R, handles exact 0/1 (common for PSI) without +/-Inf.
#        - Unbounded-nonnegative modalities: log2(x + 1)
#   3. ComBat (batch = cohort, par.prior = TRUE; fallback non-parametric)
#   4. Quantile normalization (across samples) + rank-based INT (per feature)
#   5. Write ComBat+INT BED + phenotype_groups.txt + diagnostics PDF
#
# Normalization order (ComBat -> INT) follows the Brain xQTL approach
# (Wen et al. 2024) used in the expression/HCP pipeline.
#
# Usage:
#   Rscript combat_normalize_modalities.R \
#       --input EUR_alt_TSS_pooled.bed \
#       --ancestry-map pooled_sample_ancestry_RNAseq.tsv \
#       --ancestry EUR \
#       --modality alt_TSS \
#       --output-dir combat_output/
#
# For pre-pooled modalities (splicing, intron_retention), pass the cohort-label
# sidecar produced by pool_modalities_within_ancestry.py:
#       --cohort-labels EUR_splicing_cohort_labels.tsv
#
# Dependencies:
#   sva       (Bioconductor) — ComBat
#   optparse  (CRAN)          — CLI parsing
#

.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
})

# ---- Modality -> transform type map (single source of truth) ----
# "logit"  : bounded [0,1] proportions  (alt_TSS, alt_polyA, splicing,
#             intron_retention, RNA_editing)
# "log2"   : unbounded non-negative     (isoforms, stability; expression
#             handled by HCP script)
# NOTE (GTEx-conventions round): isoforms is isoform EXPRESSION (Salmon
# quantification), not a ratio — it is unbounded and must NOT be
# logit-transformed. Fixed from "logit" to "log2".
MODALITY_TRANSFORMS <- list(
  isoforms         = "log2",
  alt_TSS          = "logit",
  alt_polyA        = "logit",
  splicing         = "logit",
  intron_retention = "logit",
  RNA_editing      = "logit",
  stability        = "log2"
)

# ---- CLI ----
option_list <- list(
  make_option("--input", type = "character",
              help = "Pooled unnorm BED file for one modality"),
  make_option("--ancestry-map", type = "character", default = NA,
              help = "Ancestry map TSV: sample_id, assigned_ancestry, cohort"),
  make_option("--cohort-labels", type = "character", default = NA,
              help = paste("Cohort-label sidecar TSV (sample_id, cohort) for",
                           "pre-pooled modalities (splicing/intron_retention)")),
  make_option("--ancestry", type = "character",
              help = "Ancestry stratum to process (e.g. EUR, EAS, AFR, HIS)"),
  make_option("--modality", type = "character",
              help = paste("Modality name; one of:",
                           paste(names(MODALITY_TRANSFORMS), collapse = ", "))),
  make_option("--eps", type = "double", default = 1e-4,
              help = "Epsilon for logit clipping (default 1e-4)"),
  make_option("--output-dir", type = "character", default = ".",
              help = "Output directory"),
  make_option("--nonparametric", action = "store_true", default = FALSE,
              help = "Use non-parametric ComBat (if parametric fails to converge)"),
  make_option("--skip-combat", action = "store_true", default = FALSE,
              help = "Skip ComBat (for single-cohort strata; INT only)")
)

opt <- parse_args(OptionParser(option_list = option_list))

opt_input          <- opt[["input"]]
opt_ancestry_map   <- opt[["ancestry-map"]]
opt_cohort_labels  <- opt[["cohort-labels"]]
opt_ancestry       <- opt[["ancestry"]]
opt_modality       <- opt[["modality"]]
opt_eps            <- opt[["eps"]]
opt_output_dir     <- opt[["output-dir"]]
opt_nonparametric  <- opt[["nonparametric"]]
opt_skip_combat    <- opt[["skip-combat"]]

if (is.null(opt_input) || is.null(opt_ancestry) || is.null(opt_modality)) {
  stop("--input, --ancestry, and --modality are required")
}
if (!(opt_modality %in% names(MODALITY_TRANSFORMS))) {
  stop("Unknown modality '", opt_modality, "'. Valid: ",
       paste(names(MODALITY_TRANSFORMS), collapse = ", "))
}
transform_type <- MODALITY_TRANSFORMS[[opt_modality]]

dir.create(opt_output_dir, showWarnings = FALSE, recursive = TRUE)
ancestry <- opt_ancestry
cat(sprintf("[%s] ComBat + INT for %s / %s (transform: %s)\n",
            date(), ancestry, opt_modality, transform_type))

# ---- Helper functions (match combat_normalize_hcp.R) ----

quantile_normalize_rows <- function(mat) {
  # Impute NAs with column (sample) medians so that apply(,2,sort) returns
  # a proper matrix. Without this, sort() drops NAs and columns with
  # different NA counts produce vectors of inconsistent length, causing
  # rowMeans() to fail with "'x' must be an array of at least two dimensions".
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
  t_mat <- t(mat)
  sorted_mat <- apply(t_mat, 2, sort)
  mean_dist <- rowMeans(sorted_mat)
  ranked_mat <- apply(t_mat, 2, rank, ties.method = "average")
  result <- matrix(0, nrow = nrow(t_mat), ncol = ncol(t_mat))
  for (j in seq_len(ncol(t_mat))) {
    result[ranked_mat[, j], j] <- mean_dist
  }
  dimnames(result) <- dimnames(t_mat)
  t(result)
}

inverse_normal_transform <- function(mat) {
  result <- apply(mat, 2, function(x) {
    r <- rank(x, ties.method = "average", na.last = "keep")
    qnorm(r / (sum(!is.na(x)) + 1))
  })
  dimnames(result) <- dimnames(mat)
  result
}

# ---- Pre-transforms ----

logit_transform <- function(mat, eps) {
  # Clip to [eps, 1-eps] to avoid +/-Inf at exact 0 or 1 (common for PSI).
  mat <- as.matrix(mat)
  mat[!is.finite(mat)] <- NA
  p <- pmin(pmax(mat, eps), 1 - eps)
  log(p / (1 - p))
}

log2_transform <- function(mat) {
  mat <- as.matrix(mat)
  mat[!is.finite(mat)] <- NA
  mat[mat < 0] <- 0
  log2(mat + 1)
}

# ---- Load data ----
cat(sprintf("[%s] Loading BED: %s\n", date(), opt_input))
# Read header first to size colClasses exactly (avoids the "cols != length"
# warning from over-specifying numeric classes, and prevents phantom columns).
hdr <- readLines(opt_input, n = 1)
ncols <- length(strsplit(hdr, "\t")[[1]])
bed <- read.delim(opt_input, sep = "\t", check.names = FALSE,
                  colClasses = c("character", "integer", "integer",
                                 "character", rep("numeric", ncols - 4)))
meta_cols <- c("#chr", "start", "end", "phenotype_id")
sample_cols <- setdiff(colnames(bed), meta_cols)
cat(sprintf("  %d features x %d samples\n", nrow(bed), length(sample_cols)))

# ---- Resolve cohort labels (batch variable for ComBat) ----
# Two sources:
#   (a) cohort-label sidecar (pre-pooled modalities, namespaced IDs)
#   (b) ancestry map join on sample_id (intersection modalities, original IDs)
if (!is.na(opt_cohort_labels) && file.exists(opt_cohort_labels)) {
  cat(sprintf("[%s] Loading cohort labels: %s\n", date(), opt_cohort_labels))
  cl <- read.delim(opt_cohort_labels, sep = "\t", stringsAsFactors = FALSE)
  cohort_labels <- setNames(cl$cohort, cl$sample_id)
  cohort_labels <- cohort_labels[sample_cols]
} else if (!is.na(opt_ancestry_map) && file.exists(opt_ancestry_map)) {
  cat(sprintf("[%s] Loading ancestry map: %s\n", date(), opt_ancestry_map))
  anc_map <- read.delim(opt_ancestry_map, sep = "\t", stringsAsFactors = FALSE)
  colnames(anc_map) <- tolower(colnames(anc_map))
  stratum <- anc_map[anc_map$assigned_ancestry == ancestry, ]
  cohort_labels <- setNames(stratum$cohort, stratum$sample_id)
  cohort_labels <- cohort_labels[sample_cols]
  missing_labels <- sample_cols[is.na(cohort_labels)]
  if (length(missing_labels) > 0) {
    cat(sprintf("  WARN: %d samples missing cohort label from ancestry map; ",
                length(missing_labels)))
    # Try namespaced fallback (cohort prefix) for any still missing
    for (s in missing_labels) {
      if (grepl("_", s)) {
        cohort_labels[s] <- sub("^[^_]+_.*", "\\1", s)
      } else {
        cohort_labels[s] <- "unknown"
      }
    }
    cat("filled via namespaced prefix or 'unknown'\n")
  }
} else {
  stop("Need either --cohort-labels or --ancestry-map to resolve batch labels")
}

cohort_labels <- as.character(cohort_labels)
names(cohort_labels) <- sample_cols
batch_table <- table(cohort_labels)
cat(sprintf("[%s] Batch (cohort) structure:\n", date()))
for (b in names(batch_table)) {
  cat(sprintf("    %s: %d samples\n", b, batch_table[b]))
}

# ---- Build data matrix: samples x features ----
data_mat <- t(as.matrix(bed[, sample_cols, drop = FALSE]))
rownames(data_mat) <- sample_cols
colnames(data_mat) <- bed$phenotype_id
cat(sprintf("  Data matrix: %d samples x %d features\n",
            nrow(data_mat), ncol(data_mat)))

# ---- Filter features: >50% zeros, no variance, or excessive NAs ----
zero_frac <- colMeans(data_mat == 0, na.rm = TRUE)
no_var <- apply(data_mat, 2, function(x) length(unique(x[!is.na(x)])) <= 1)
# Drop features where any batch has <2 non-NA values (ComBat can't estimate
# batch effects for them → singular design matrix)
na_per_batch <- sapply(split(seq_len(nrow(data_mat)), cohort_labels), function(idx) {
  colSums(!is.na(data_mat[idx, , drop = FALSE]))
})
min_batch_obs <- apply(na_per_batch, 1, min)
na_excessive <- min_batch_obs < 2
n_na_excessive <- sum(na_excessive & !(zero_frac > 0.5 | no_var))
keep <- (zero_frac <= 0.5) & !no_var & !na_excessive
cat(sprintf("  Feature filtering: %d -> %d (removed %d >50%% zeros, %d no-variance, %d insufficient obs per batch)\n",
            ncol(data_mat), sum(keep),
            sum(zero_frac > 0.5), sum(no_var & zero_frac <= 0.5 & !na_excessive),
            n_na_excessive))
data_mat <- data_mat[, keep, drop = FALSE]

# ---- Impute remaining NAs with batch-specific medians ----
if (any(is.na(data_mat))) {
  n_na <- sum(is.na(data_mat))
  for (b in unique(cohort_labels)) {
    b_idx <- which(cohort_labels == b)
    for (f in seq_len(ncol(data_mat))) {
      na_idx <- which(is.na(data_mat[b_idx, f]))
      if (length(na_idx) > 0) {
        med <- median(data_mat[b_idx[-na_idx], f], na.rm = TRUE)
        if (is.na(med)) med <- median(data_mat[, f], na.rm = TRUE)
        if (is.na(med)) med <- 0
        data_mat[b_idx[na_idx], f] <- med
      }
    }
  }
  cat(sprintf("  Imputed %d NA values with batch-specific medians\n", n_na))
}

# ---- Pre-transform ----
if (transform_type == "logit") {
  cat(sprintf("[%s] logit transform (eps=%g)\n", date(), opt_eps))
  data_tx <- logit_transform(data_mat, opt_eps)
} else {
  cat(sprintf("[%s] log2(x+1) transform\n", date()))
  data_tx <- log2_transform(data_mat)
}

# ---- ComBat (batch = cohort) ----
if (opt_skip_combat || length(unique(cohort_labels)) < 2) {
  cat(sprintf("[%s] Skipping ComBat (single cohort or --skip-combat)\n", date()))
  data_combat <- data_tx
} else {
  # Merge single-sample batches into "other"
  single_batches <- names(batch_table[batch_table < 2])
  if (length(single_batches) > 0) {
    cat(sprintf("  WARN: cohorts with <2 samples: %s -> merged into 'other'\n",
                paste(single_batches, collapse = ", ")))
    cohort_labels[cohort_labels %in% single_batches] <- "other"
  }

  cat(sprintf("[%s] Running ComBat (batch = cohort, par.prior = %s)\n",
              date(), ifelse(opt_nonparametric, "FALSE", "TRUE")))
  suppressPackageStartupMessages(library(sva))

  mod <- model.matrix(~1, data = data.frame(row.names = sample_cols))

  data_combat <- tryCatch({
    ComBat(dat = t(data_tx),  # ComBat expects features x samples
           batch = as.factor(cohort_labels),
           mod = mod,
           par.prior = !opt_nonparametric,
           prior.plots = FALSE)
  }, error = function(e) {
    cat(sprintf("  ComBat (parametric) failed: %s\n", e$message))
    cat("  Retrying with non-parametric prior...\n")
    tryCatch({
      ComBat(dat = t(data_tx),
             batch = as.factor(cohort_labels),
             mod = mod,
             par.prior = FALSE,
             prior.plots = FALSE)
    }, error = function(e2) {
      cat(sprintf("  ComBat (non-parametric) also failed: %s\n", e2$message))
      cat("  WARNING: Skipping ComBat for this stratum; using transformed data without batch correction\n")
      t(data_tx)  # return features x samples, same as ComBat output
    })
  })
  data_combat <- t(data_combat)  # back to samples x features
  cat(sprintf("  ComBat done: %d samples x %d features\n",
              nrow(data_combat), ncol(data_combat)))
}

# ---- Quantile normalization + rank-based INT ----
cat(sprintf("[%s] Quantile normalization + rank-based INT\n", date()))
data_qn <- quantile_normalize_rows(data_combat)
data_int <- inverse_normal_transform(data_qn)
cat(sprintf("  INT done: mean=%.4f, sd=%.4f (should be ~0, ~1)\n",
            mean(data_int, na.rm = TRUE), sd(data_int, na.rm = TRUE)))

# ---- Write ComBat + INT BED ----
feat_idx <- match(colnames(data_int), bed$phenotype_id)
out_bed <- data.frame(
  chr = bed[["#chr"]][feat_idx],
  start = bed$start[feat_idx],
  end = bed$end[feat_idx],
  phenotype_id = colnames(data_int),
  t(data_int),
  check.names = FALSE
)
colnames(out_bed)[1] <- "#chr"
bed_path <- file.path(opt_output_dir,
                      sprintf("%s_%s_combat_int.bed", ancestry, opt_modality))
write.table(out_bed, bed_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Written: %s (%d features x %d samples)\n",
            bed_path, nrow(out_bed), ncol(out_bed) - 4))

# ---- Phenotype groups (gene grouping for tensorQTL) ----
# phenotype_id format: {gene_id}__{rest} -> group by gene_id
groups_path <- file.path(opt_output_dir,
                         sprintf("%s_%s.phenotype_groups.txt", ancestry, opt_modality))
gene_ids <- sub("__.*$", "", out_bed$phenotype_id)
groups_df <- data.frame(phenotype_id = out_bed$phenotype_id,
                        gene_id = gene_ids)
write.table(groups_df, groups_path, sep = "\t", row.names = FALSE,
            quote = FALSE, col.names = FALSE)
cat(sprintf("  Written: %s\n", groups_path))

# ---- Diagnostics ----
cat(sprintf("[%s] Generating diagnostics\n", date()))
diag_path <- file.path(opt_output_dir,
                       sprintf("%s_%s_combat_diagnostics.pdf", ancestry, opt_modality))
pdf(diag_path, width = 12, height = 5)
par(mfrow = c(1, 2))

cohort_colors <- as.numeric(as.factor(cohort_labels[rownames(data_tx)]))

# PCA before ComBat — impute NAs with column medians for prcomp
pca_before_data <- data_tx
if (any(is.na(pca_before_data))) {
  for (j in seq_len(ncol(pca_before_data))) {
    na_idx <- is.na(pca_before_data[, j])
    if (any(na_idx)) {
      med <- median(pca_before_data[!na_idx, j], na.rm = TRUE)
      if (is.na(med)) med <- 0
      pca_before_data[na_idx, j] <- med
    }
  }
}
tryCatch({
  pca_before <- prcomp(pca_before_data, scale. = FALSE, center = TRUE)
  plot(pca_before$x[, 1:2], col = cohort_colors, pch = 19, cex = 0.8,
       xlab = sprintf("PC1 (%.1f%%)", summary(pca_before)$importance[2, 1] * 100),
       ylab = sprintf("PC2 (%.1f%%)", summary(pca_before)$importance[2, 2] * 100),
       main = sprintf("PCA before ComBat (%s)", opt_modality))
  legend("topright", legend = levels(as.factor(cohort_labels)),
         col = seq_along(levels(as.factor(cohort_labels))), pch = 19, cex = 0.6)
}, error = function(e) {
  cat(sprintf("  WARN: PCA before ComBat failed: %s\n", e$message))
  plot.new(); text(0.5, 0.5, "PCA before ComBat failed")
})

tryCatch({
  pca_after <- prcomp(data_int, scale. = FALSE, center = TRUE)
  plot(pca_after$x[, 1:2], col = cohort_colors, pch = 19, cex = 0.8,
       xlab = sprintf("PC1 (%.1f%%)", summary(pca_after)$importance[2, 1] * 100),
       ylab = sprintf("PC2 (%.1f%%)", summary(pca_after)$importance[2, 2] * 100),
       main = sprintf("PCA after ComBat + INT (%s)", opt_modality))
  legend("topright", legend = levels(as.factor(cohort_labels)),
         col = seq_along(levels(as.factor(cohort_labels))), pch = 19, cex = 0.6)
}, error = function(e) {
  cat(sprintf("  WARN: PCA after ComBat failed: %s\n", e$message))
  plot.new(); text(0.5, 0.5, "PCA after ComBat failed")
})

dev.off()
cat(sprintf("  Written: %s\n", diag_path))

cat(sprintf("[%s] Done: %s / %s\n", date(), ancestry, opt_modality))
cat(sprintf("  BED:     %s\n", bed_path))
cat(sprintf("  Groups:  %s\n", groups_path))
cat(sprintf("  Diag:    %s\n", diag_path))
