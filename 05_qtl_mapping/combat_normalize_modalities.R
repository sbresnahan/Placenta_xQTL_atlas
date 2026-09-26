#!/usr/bin/env Rscript
#
# combat_normalize_modalities.R — QN + INT + ComBat for non-expression RNA modalities
#
# Generalizes combat_normalize_hcp.R to all RNA phenotype modalities, WITHOUT
# HCP estimation (HCP factors are estimated once on gene-level expression and
# reused as covariates for every modality).
#
# SCHEMA (2026-09 revision, modeled on the devBrain xQTL atlas — Wen et al.,
# Science 2024, 384:eadh0829): normalization FIRST, batch correction LAST.
# Per ancestry stratum, per modality:
#   1. Load pooled unnorm BED (samples x features)
#   2. Optional sample exclusion (--exclude-samples): expression-outlier list
#      from script 19, applied to isoforms only (devBrain §3.3 excludes
#      connectivity outliers from expression/isoforms but not splicing)
#   3. Feature filter (devBrain §3.3/§3.4):
#        - isoforms: TPM > 0.1 in > 25% of stratum samples
#        - proportion/ratio modalities (alt_TSS, alt_polyA, splicing,
#          intron_retention, RNA_editing) and stability: detected (non-NA)
#          in >= 40% of stratum samples
#        - no-variance removal; per-batch >= 2 non-NA guard for ComBat
#      (Zeros are real PSI values — the old >50%-zeros rule is removed.)
#   4. Batch-median NA imputation
#   5. Quantile normalization (across samples) + rank-based INT (per feature)
#   6. ComBat (batch = cohort, par.prior = TRUE; fallback non-parametric)
#   7. Write ComBat'd BED + phenotype_groups.txt + diagnostics PDF
#
# No log2/logit pre-transform: rank-based QN + INT is invariant to monotone
# transforms, so they were no-ops on the final scale.
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
# For isoforms and isoform_expression, pass the expression-outlier list from
# script 19:
#       --exclude-samples EUR_expression_outliers.tsv
#
# Dependencies:
#   sva       (Bioconductor) — ComBat
#   optparse  (CRAN)          — CLI parsing
#

.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
})

# ---- Valid modalities (single source of truth) ----
# expression is handled by combat_normalize_hcp.R (script 19).
# isoforms = within-gene usage ratios; isoform_expression = transcript-level
# TPM abundance (the pre-ratio matrix from assemble_bed.py).
MODALITIES <- c("isoforms", "isoform_expression", "alt_TSS", "alt_polyA",
                "splicing", "intron_retention", "RNA_editing", "stability")

# Modalities filtered on TPM detection (devBrain §3.3); all others use the
# >=40% detection filter (devBrain §3.4). For isoforms (usage ratios) the
# tpm-min threshold acts as a usage > 0.1 filter; for isoform_expression the
# matrix holds actual TPMs, so the filter is dimensionally the devBrain filter.
TPM_MODALITIES <- c("isoforms", "isoform_expression")

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
                           paste(MODALITIES, collapse = ", "))),
  make_option("--exclude-samples", type = "character", default = NA,
              help = paste("Optional TSV of sample IDs to exclude before",
                           "filtering (one column 'sample_id'; e.g. expression",
                           "outliers from script 19, applied to isoforms)")),
  make_option("--tpm-min", type = "double", default = 0.1,
              help = "isoforms filter: minimum TPM (default: 0.1, devBrain §3.3)"),
  make_option("--tpm-min-prop", type = "double", default = 0.25,
              help = "isoforms filter: required fraction of samples above --tpm-min (default: 0.25)"),
  make_option("--min-detect-prop", type = "double", default = 0.4,
              help = "Proportion modalities filter: required fraction of samples with detected (non-NA) values (default: 0.4, devBrain §3.4)"),
  make_option("--output-dir", type = "character", default = ".",
              help = "Output directory"),
  make_option("--nonparametric", action = "store_true", default = FALSE,
              help = "Use non-parametric ComBat (if parametric fails to converge)"),
  make_option("--skip-combat", action = "store_true", default = FALSE,
              help = "Skip ComBat (for single-cohort strata; QN + INT only)")
)

opt <- parse_args(OptionParser(option_list = option_list))

opt_input           <- opt[["input"]]
opt_ancestry_map    <- opt[["ancestry-map"]]
opt_cohort_labels   <- opt[["cohort-labels"]]
opt_ancestry        <- opt[["ancestry"]]
opt_modality        <- opt[["modality"]]
opt_exclude_samples <- opt[["exclude-samples"]]
opt_tpm_min         <- opt[["tpm-min"]]
opt_tpm_min_prop    <- opt[["tpm-min-prop"]]
opt_min_detect_prop <- opt[["min-detect-prop"]]
opt_output_dir      <- opt[["output-dir"]]
opt_nonparametric   <- opt[["nonparametric"]]
opt_skip_combat     <- opt[["skip-combat"]]

if (is.null(opt_input) || is.null(opt_ancestry) || is.null(opt_modality)) {
  stop("--input, --ancestry, and --modality are required")
}
if (!(opt_modality %in% MODALITIES)) {
  stop("Unknown modality '", opt_modality, "'. Valid: ",
       paste(MODALITIES, collapse = ", "))
}

dir.create(opt_output_dir, showWarnings = FALSE, recursive = TRUE)
ancestry <- opt_ancestry
cat(sprintf("[%s] QN + INT + ComBat for %s / %s\n",
            date(), ancestry, opt_modality))

# ---- Helper functions (match combat_normalize_hcp.R) ----

quantile_normalize_rows <- function(mat) {
  # Impute NAs with column (feature) medians so that apply(,2,sort) below
  # returns a proper matrix. Without this, sort() drops NAs and columns with
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
  result <- matrix(0, nrow = nrow(t_mat), ncol = ncol(t_mat))
  # Assign by ORDER (preprocessCore / PANTRY normalize_phenotypes.py
  # convention): the smallest value in a sample gets the smallest mean
  # quantile. Tied values get the mean of the quantiles they span.
  # (The previous implementation assigned result[ranked, j] <- mean_dist,
  # which applies the INVERSE permutation and scrambles values across
  # features — fixed 2026-09-25.)
  for (j in seq_len(ncol(t_mat))) {
    o <- order(t_mat[, j])
    result[o, j] <- mean_dist
    xs <- t_mat[o, j]
    tie_runs <- rle(xs)
    if (any(tie_runs$lengths > 1)) {
      ends <- cumsum(tie_runs$lengths)
      starts <- ends - tie_runs$lengths + 1
      for (g in which(tie_runs$lengths > 1)) {
        result[o[starts[g]:ends[g]], j] <- mean(mean_dist[starts[g]:ends[g]])
      }
    }
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

# ---- Optional sample exclusion (expression outliers -> isoforms) ----
if (!is.na(opt_exclude_samples) && file.exists(opt_exclude_samples)) {
  excl <- read.delim(opt_exclude_samples, sep = "\t", stringsAsFactors = FALSE)
  excl_ids <- excl[[1]]
  drop <- intersect(sample_cols, excl_ids)
  if (length(drop) > 0) {
    cat(sprintf("  Excluding %d samples listed in %s: %s\n",
                length(drop), basename(opt_exclude_samples),
                paste(drop, collapse = ", ")))
    sample_cols <- setdiff(sample_cols, excl_ids)
  } else {
    cat(sprintf("  Exclude list %s: no overlapping samples\n",
                basename(opt_exclude_samples)))
  }
}

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

# Merge single-sample batches into "other" BEFORE filtering: the per-batch
# >= 2 non-NA guard below must reflect the batch structure ComBat will
# actually see. Batches with <2 samples after merging are excluded from the
# guard (ComBat cannot estimate their effects; the fallback chain in the
# ComBat step handles them).
batch_table <- table(cohort_labels)
single_batches <- names(batch_table[batch_table < 2])
if (length(single_batches) > 0 && length(unique(cohort_labels)) > 1) {
  cat(sprintf("  WARN: cohorts with <2 samples: %s -> merged into 'other'\n",
              paste(single_batches, collapse = ", ")))
  cohort_labels[cohort_labels %in% single_batches] <- "other"
}
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

# ---- Filter features (devBrain §3.3/§3.4) ----
# isoforms: TPM > --tpm-min in > --tpm-min-prop of samples.
# All other modalities: detected (non-NA) in >= --min-detect-prop of samples.
# Zeros are real values for proportion modalities (PSI = 0), so the old
# >50%-zeros rule is removed. No-variance and per-batch >= 2 non-NA guards
# are retained.
no_var <- apply(data_mat, 2, function(x) length(unique(x[!is.na(x)])) <= 1)
# Drop features where any ComBat-estimable batch (>= 2 samples) has <2 non-NA
# values (ComBat can't estimate batch effects for them → singular design)
estimable_batches <- names(batch_table[batch_table >= 2])
if (length(estimable_batches) > 0) {
  na_per_batch <- sapply(estimable_batches, function(b) {
    colSums(!is.na(data_mat[cohort_labels == b, , drop = FALSE]))
  })
  min_batch_obs <- apply(na_per_batch, 1, min)
  na_excessive <- min_batch_obs < 2
} else {
  na_excessive <- rep(FALSE, ncol(data_mat))
}

if (opt_modality %in% TPM_MODALITIES) {
  detect_frac <- colMeans(data_mat > opt_tpm_min, na.rm = TRUE)
  low_detect <- detect_frac <= opt_tpm_min_prop
  filter_desc <- sprintf("TPM <= %g in >= %d%% of samples",
                         opt_tpm_min, as.integer(opt_tpm_min_prop * 100))
} else {
  detect_frac <- colMeans(!is.na(data_mat))
  low_detect <- detect_frac < opt_min_detect_prop
  filter_desc <- sprintf("detected in < %d%% of samples",
                         as.integer(opt_min_detect_prop * 100))
}
keep <- !low_detect & !no_var & !na_excessive
cat(sprintf("  Feature filtering: %d -> %d (removed %d %s, %d no-variance, %d insufficient obs per batch)\n",
            ncol(data_mat), sum(keep),
            sum(low_detect), filter_desc,
            sum(no_var & !low_detect & !na_excessive),
            sum(na_excessive & !low_detect & !no_var)))
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

# ---- Quantile normalization + rank-based INT ----
cat(sprintf("[%s] Quantile normalization + rank-based INT\n", date()))
data_qn <- quantile_normalize_rows(data_mat)
data_int <- inverse_normal_transform(data_qn)
cat(sprintf("  INT done: mean=%.4f, sd=%.4f (should be ~0, ~1)\n",
            mean(data_int, na.rm = TRUE), sd(data_int, na.rm = TRUE)))

# ---- ComBat (batch = cohort), on the INT'd values (ComBat LAST) ----
if (opt_skip_combat || length(unique(cohort_labels)) < 2) {
  cat(sprintf("[%s] Skipping ComBat (single cohort or --skip-combat)\n", date()))
  data_combat <- data_int
} else {
  cat(sprintf("[%s] Running ComBat (batch = cohort, par.prior = %s)\n",
              date(), ifelse(opt_nonparametric, "FALSE", "TRUE")))
  suppressPackageStartupMessages(library(sva))

  mod <- model.matrix(~1, data = data.frame(row.names = sample_cols))

  data_combat <- tryCatch({
    ComBat(dat = t(data_int),  # ComBat expects features x samples
           batch = as.factor(cohort_labels),
           mod = mod,
           par.prior = !opt_nonparametric,
           prior.plots = FALSE)
  }, error = function(e) {
    cat(sprintf("  ComBat (parametric) failed: %s\n", e$message))
    cat("  Retrying with non-parametric prior...\n")
    tryCatch({
      ComBat(dat = t(data_int),
             batch = as.factor(cohort_labels),
             mod = mod,
             par.prior = FALSE,
             prior.plots = FALSE)
    }, error = function(e2) {
      cat(sprintf("  ComBat (non-parametric) also failed: %s\n", e2$message))
      cat("  WARNING: Skipping ComBat for this stratum; using INT data without batch correction\n")
      t(data_int)  # return features x samples, same as ComBat output
    })
  })
  data_combat <- t(data_combat)  # back to samples x features
  cat(sprintf("  ComBat done: %d samples x %d features\n",
              nrow(data_combat), ncol(data_combat)))
}

# ---- Write ComBat BED ----
feat_idx <- match(colnames(data_combat), bed$phenotype_id)
out_bed <- data.frame(
  chr = bed[["#chr"]][feat_idx],
  start = bed$start[feat_idx],
  end = bed$end[feat_idx],
  phenotype_id = colnames(data_combat),
  t(data_combat),
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

cohort_colors <- as.numeric(as.factor(cohort_labels[rownames(data_int)]))

tryCatch({
  pca_before <- prcomp(data_int, scale. = FALSE, center = TRUE)
  plot(pca_before$x[, 1:2], col = cohort_colors, pch = 19, cex = 0.8,
       xlab = sprintf("PC1 (%.1f%%)", summary(pca_before)$importance[2, 1] * 100),
       ylab = sprintf("PC2 (%.1f%%)", summary(pca_before)$importance[2, 2] * 100),
       main = sprintf("PCA after QN + INT, pre-ComBat (%s)", opt_modality))
  legend("topright", legend = levels(as.factor(cohort_labels)),
         col = seq_along(levels(as.factor(cohort_labels))), pch = 19, cex = 0.6)
}, error = function(e) {
  cat(sprintf("  WARN: PCA before ComBat failed: %s\n", e$message))
  plot.new(); text(0.5, 0.5, "PCA before ComBat failed")
})

tryCatch({
  pca_after <- prcomp(data_combat, scale. = FALSE, center = TRUE)
  plot(pca_after$x[, 1:2], col = cohort_colors, pch = 19, cex = 0.8,
       xlab = sprintf("PC1 (%.1f%%)", summary(pca_after)$importance[2, 1] * 100),
       ylab = sprintf("PC2 (%.1f%%)", summary(pca_after)$importance[2, 2] * 100),
       main = sprintf("PCA after ComBat, final (%s)", opt_modality))
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
