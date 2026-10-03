#!/usr/bin/env Rscript
# Lab R package library (edgeR/DESeq2/Rhcpp/tximport/...) -- must precede any library() call.
.libPaths(c("/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1", .libPaths()))
#
# normalize_modalities.R — QN + within-cohort INT for non-expression RNA modalities
# (formerly combat_normalize_modalities.R); TMM -> VST for isoform_expression
#
# Generalizes normalize_expression_hcp.R to all RNA phenotype modalities,
# WITHOUT HCP estimation (HCP factors are estimated once on gene-level
# expression and reused as covariates for every modality).
#
# SCHEMA: pooled quantile normalization, then per-feature rank-based INT
# separately within each cohort. ComBat is DROPPED (schema change, Oct 2026) —
# see normalize_expression_hcp.R header for the rationale; the last commit
# using the ComBat schema is ea2b344ab9cd7cdf116b742d97f072d3f7597186.
#
# EXCEPTION — isoform_expression with --counts-input (PsychENCODE/isoTWAS
# schema): the TPM matrix drives the devBrain detection filter EXACTLY as
# below (feature-set parity with the INT path), then the filtered feature set
# is read off the pooled QU-corrected counts BED and normalized by edgeR TMM
# -> DESeq2 VST (terminal; NO QN, NO INT, NO NA imputation). Output is
# {ANC}_isoform_expression_vst.bed. Cohort scale heterogeneity is NOT
# equilibrated by construction in this mode — cohort enters via dummies in
# the HCP known-covariate matrix (hcp_from_matrix.R --cohort-dummies) and the
# tensorQTL covariate model; the 29b-29f gate re-verifies cohort variance.
# Without --counts-input, isoform_expression falls back to the INT path
# (the --mode int-style fallback for the diagnostic gate).
#
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
#        - no-variance removal; per-cohort >= 2 non-NA guard (retained as a
#          conservative filter; it was originally a ComBat design guard)
#      (Zeros are real PSI values — the old >50%-zeros rule is removed.)
#   4. Cohort-median NA imputation (INT path only)
#   5. Quantile normalization (across all stratum samples) (INT path only)
#   6. Rank-based INT per feature SEPARATELY WITHIN EACH COHORT (every cohort,
#      regardless of size), then concatenate back to one pooled matrix
#      (INT path only) — OR TMM -> VST on counts (VST path)
#   7. Write normalized BED + phenotype_groups.txt + diagnostics PDF
#
# No log2/logit pre-transform: rank-based QN + INT is invariant to monotone
# transforms, so they were no-ops on the final scale.
#
# Usage:
#   Rscript normalize_modalities.R \
#       --input EUR_alt_TSS_pooled.bed \
#       --ancestry-map pooled_sample_ancestry_RNAseq.tsv \
#       --ancestry EUR \
#       --modality alt_TSS \
#       --output-dir normalized_output/
#
# For pre-pooled modalities (splicing, intron_retention), pass the cohort-label
# sidecar produced by pool_modalities_within_ancestry.py:
#       --cohort-labels EUR_splicing_cohort_labels.tsv
#
# For isoforms and isoform_expression, pass the expression-outlier list from
# script 19:
#       --exclude-samples EUR_expression_outliers.tsv
#
# For isoform_expression TMM -> VST, additionally pass the pooled counts BED:
#       --counts-input EUR_isoform_expression_counts_pooled.bed
#
# Dependencies:
#   optparse  (CRAN)          — CLI parsing
#   edgeR, DESeq2 (Bioconductor) — VST path only (isoform_expression +
#                                  --counts-input)
#

.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
})

# ---- Valid modalities (single source of truth) ----
# expression is handled by normalize_expression_hcp.R (script 19).
# isoforms = within-gene usage ratios; isoform_expression = transcript-level
# TPM abundance (the pre-ratio matrix from assemble_bed.py).
MODALITIES <- c("isoforms", "isoform_expression", "alt_TSS", "alt_polyA",
                "splicing", "intron_retention", "RNA_editing", "stability")

# Modalities filtered on TPM detection (devBrain §3.3); all others use the
# >=40% detection filter (devBrain §3.4). For isoforms (usage ratios) the
# tpm-min threshold acts as a usage > 0.1 filter; for isoform_expression the
# matrix holds length-normalized TPMs (recomputed from QU-corrected counts x
# Salmon effective lengths by assemble_bed.py --units tpm_from_counts), so
# the filter is dimensionally the devBrain filter.
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
  make_option("--counts-input", type = "character", default = NA,
              help = paste("Pooled counts BED (isoform_expression only): triggers",
                           "the TMM -> VST path instead of QN + within-cohort INT")),
  make_option("--tpm-min", type = "double", default = 0.1,
              help = "isoforms filter: minimum TPM (default: 0.1, devBrain §3.3)"),
  make_option("--tpm-min-prop", type = "double", default = 0.25,
              help = "isoforms filter: required fraction of samples above --tpm-min (default: 0.25)"),
  make_option("--min-detect-prop", type = "double", default = 0.4,
              help = "Proportion modalities filter: required fraction of samples with detected (non-NA) values (default: 0.4, devBrain §3.4)"),
  make_option("--output-dir", type = "character", default = ".",
              help = "Output directory"),
  make_option("--nonparametric", action = "store_true", default = FALSE,
              help = "DEPRECATED no-op: ComBat was removed from this script (within-cohort INT schema). Accepted for backward compatibility."),
  make_option("--skip-combat", action = "store_true", default = FALSE,
              help = "DEPRECATED no-op: ComBat was removed from this script (within-cohort INT schema). Accepted for backward compatibility.")
)

opt <- parse_args(OptionParser(option_list = option_list))

opt_input           <- opt[["input"]]
opt_ancestry_map    <- opt[["ancestry-map"]]
opt_cohort_labels   <- opt[["cohort-labels"]]
opt_ancestry        <- opt[["ancestry"]]
opt_modality        <- opt[["modality"]]
opt_exclude_samples <- opt[["exclude-samples"]]
opt_counts_input    <- opt[["counts-input"]]
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

# VST path: isoform_expression + --counts-input (PsychENCODE/isoTWAS schema).
use_vst <- !is.na(opt_counts_input)
if (use_vst && opt_modality != "isoform_expression") {
  stop("--counts-input is only valid for --modality isoform_expression ",
       "(got '", opt_modality, "'). All other modalities use QN + within-cohort INT.")
}
if (use_vst) {
  for (pkg in c("edgeR", "DESeq2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("--counts-input (VST path) requires R package '%s' (not found in .libPaths(): %s)",
                   pkg, paste(.libPaths(), collapse = ", ")))
    }
  }
}

dir.create(opt_output_dir, showWarnings = FALSE, recursive = TRUE)
ancestry <- opt_ancestry
cat(sprintf("[%s] %s for %s / %s\n",
            date(),
            if (use_vst) "TMM -> VST" else "QN + within-cohort INT",
            ancestry, opt_modality))
if (opt_modality == "isoform_expression" && !use_vst) {
  cat("  NOTE: no --counts-input given — isoform_expression uses the INT path",
      "(fallback for the 29b-29f diagnostic gate).\n")
}

if (opt_nonparametric || opt_skip_combat) {
  cat("  NOTE: --nonparametric/--skip-combat are deprecated no-ops; ComBat was",
      "removed from this pipeline (within-cohort INT schema).\n")
}

# ---- Helper functions (match normalize_expression_hcp.R) ----

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
  # (Assigning result[ranked, j] <- mean_dist here would apply the
  # INVERSE permutation and scramble values across features.)
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
  # apply() simplifies to a vector when mat has a single row (single-sample
  # cohort): restore the matrix shape before setting dimnames
  if (!is.matrix(result)) {
    result <- matrix(result, nrow = nrow(mat), ncol = ncol(mat))
  }
  dimnames(result) <- dimnames(mat)
  result
}

# TMM -> VST bridge (PsychENCODE/isoTWAS convention; same implementation as
# normalize_expression_hcp.R — duplicated per repo convention).
# Input: samples x features count matrix (non-integer counts are fine — TMM is
# computed on the unrounded values; rounding happens only for
# DESeqDataSetFromMatrix, which requires integers).
# CRITICAL: never pre-divide counts by TMM factors and hand the result to
# DESeq2::vst() — that double-normalizes. Instead, TMM effective library
# sizes are injected as DESeq2 sizeFactors and vst() runs on raw counts:
#   eff_lib <- colSums(counts) * edgeR::calcNormFactors(counts, "TMM")
#   sizeFactors(dds) <- eff_lib / exp(mean(log(eff_lib)))   # geometric mean 1
# Returns the VST matrix in the same samples x features orientation.
tmm_vst <- function(counts_mat) {
  counts_gxs <- t(counts_mat)                       # features x samples
  lib_size <- colSums(counts_gxs)
  if (any(!is.finite(lib_size)) || any(lib_size <= 0)) {
    bad <- names(lib_size)[!is.finite(lib_size) | lib_size <= 0]
    stop("tmm_vst: ", length(bad), " sample(s) with zero/non-finite library size ",
         "(log(0) in the size-factor geometric mean): ",
         paste(utils::head(bad, 10), collapse = ", "))
  }
  tmm <- edgeR::calcNormFactors(counts_gxs, method = "TMM")
  eff_lib <- lib_size * tmm
  sf <- eff_lib / exp(mean(log(eff_lib)))           # DESeq2 size-factor scale
  dds <- DESeq2::DESeqDataSetFromMatrix(
    round(counts_gxs),
    colData = data.frame(row.names = colnames(counts_gxs)),
    design = ~ 1)
  DESeq2::sizeFactors(dds) <- sf
  # vst() stops when nrow < nsub (default 1000); cap nsub so small fixtures
  # (tests) run. Production feature counts are >> 1000, so nsub stays 1000.
  t(SummarizedExperiment::assay(
    DESeq2::vst(dds, blind = TRUE, nsub = min(1000L, nrow(counts_gxs)))))
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

# ---- Resolve cohort labels (for within-cohort INT) ----
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

# No batch merging: with ComBat dropped, every cohort is INT'd on its own
# ranks regardless of size (a single-sample cohort INTs to a constant 0 row).
# The per-cohort >= 2 non-NA guard below uses the raw cohort structure.
batch_table <- table(cohort_labels)
cat(sprintf("[%s] Cohort structure:\n", date()))
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
# Drop features where any cohort with >= 2 samples has <2 non-NA values
# (retained from the ComBat era as a conservative per-cohort detection guard;
# within-cohort INT on a single observed value would be uninformative)
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

# ---- Normalization (path-dependent) ----
if (use_vst) {
  # TMM -> VST on the pooled counts BED (PsychENCODE/isoTWAS convention;
  # terminal VST, no QN, no INT, no NA imputation). The TPM matrix above drove
  # the detection filter only — feature-set parity with the INT path.
  cat(sprintf("[%s] Loading counts BED: %s\n", date(), opt_counts_input))
  hdr_c <- readLines(opt_counts_input, n = 1)
  ncols_c <- length(strsplit(hdr_c, "\t")[[1]])
  counts_bed <- read.delim(opt_counts_input, sep = "\t", check.names = FALSE,
                           colClasses = c("character", "integer", "integer",
                                          "character", rep("numeric", ncols_c - 4)))
  missing_samples_c <- setdiff(sample_cols, colnames(counts_bed))
  if (length(missing_samples_c) > 0) {
    stop("Counts BED is missing ", length(missing_samples_c), " samples (e.g. ",
         paste(utils::head(missing_samples_c, 5), collapse = ", "), ")")
  }
  gidx <- match(colnames(data_mat), counts_bed$phenotype_id)
  if (any(is.na(gidx))) {
    stop("Counts BED is missing ", sum(is.na(gidx)), " of the ", ncol(data_mat),
         " TPM-filtered features (e.g. ",
         paste(utils::head(colnames(data_mat)[is.na(gidx)], 5), collapse = ", "),
         "). tximport_counts.R emits ALL features — regenerate the counts BEDs ",
         "rather than subsetting them upstream.")
  }
  counts_data <- t(as.matrix(counts_bed[gidx, sample_cols, drop = FALSE]))
  rownames(counts_data) <- sample_cols
  colnames(counts_data) <- colnames(data_mat)
  rm(counts_bed)
  # Guard: drop features with zero-variance (or NA) counts — VST cannot fit a
  # dispersion trend through degenerate rows.
  cv <- apply(counts_data, 2, var)
  zero_var_counts <- !is.finite(cv) | cv == 0
  if (any(zero_var_counts)) {
    cat(sprintf("  Dropping %d features with zero-variance/NA counts (VST guard)\n",
                sum(zero_var_counts)))
    counts_data <- counts_data[, !zero_var_counts, drop = FALSE]
    data_mat <- data_mat[, !zero_var_counts, drop = FALSE]
  }
  cat(sprintf("[%s] TMM -> VST (edgeR TMM size factors -> DESeq2 vst, blind)\n", date()))
  data_norm <- tmm_vst(counts_data)
  rm(counts_data)
  cat(sprintf("[%s] Normalized matrix: %d samples x %d features (TMM -> VST, terminal; no INT)\n",
              date(), nrow(data_norm), ncol(data_norm)))
} else {
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

  # ---- Quantile normalization (pooled across all stratum samples) ----
  cat(sprintf("[%s] Quantile normalization (pooled)\n", date()))
  data_qn <- quantile_normalize_rows(data_mat)

  # ---- Rank-based INT per feature, separately within each cohort ----
  # Every cohort is transformed on its own ranks regardless of size: each
  # feature is exactly N(0,1) within each cohort, so cohorts enter the pooled
  # scan on a common scale. Within-cohort sample ordering is preserved.
  cat(sprintf("[%s] Rank-based INT per feature within cohort\n", date()))
  data_int <- matrix(NA_real_, nrow = nrow(data_qn), ncol = ncol(data_qn),
                     dimnames = dimnames(data_qn))
  for (co in names(batch_table)) {
    idx <- names(cohort_labels)[cohort_labels == co]
    data_int[idx, ] <- inverse_normal_transform(data_qn[idx, , drop = FALSE])
  }
  cat(sprintf("  INT done: pooled mean=%.4f, sd=%.4f\n",
              mean(data_int, na.rm = TRUE), sd(data_int, na.rm = TRUE)))
  for (co in names(batch_table)) {
    idx <- names(cohort_labels)[cohort_labels == co]
    cat(sprintf("  INT within %s (n=%d): mean=%.4f, sd=%.4f (should be ~0, ~1)\n",
                co, length(idx),
                mean(data_int[idx, ], na.rm = TRUE), sd(data_int[idx, ], na.rm = TRUE)))
  }

  # ComBat stays dropped: within-cohort INT above already equilibrates
  # per-feature location and scale across cohorts, and a cohort fixed effect
  # in the tensorQTL covariate model (25_build_covariates.py) absorbs
  # remaining between-cohort mean structure.
  data_norm <- data_int
  cat(sprintf("[%s] Normalized matrix: %d samples x %d features (QN + within-cohort INT; no ComBat)\n",
              date(), nrow(data_norm), ncol(data_norm)))
}

# ---- Write normalized BED ----
feat_idx <- match(colnames(data_norm), bed$phenotype_id)
out_bed <- data.frame(
  chr = bed[["#chr"]][feat_idx],
  start = bed$start[feat_idx],
  end = bed$end[feat_idx],
  phenotype_id = colnames(data_norm),
  t(data_norm),
  check.names = FALSE
)
colnames(out_bed)[1] <- "#chr"
norm_suffix <- if (use_vst) "vst" else "int"
bed_path <- file.path(opt_output_dir,
                      sprintf("%s_%s_%s.bed", ancestry, opt_modality, norm_suffix))
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
                       sprintf("%s_%s_diagnostics.pdf", ancestry, opt_modality))
pdf(diag_path, width = 12, height = 5)
par(mfrow = c(1, 2))

cohort_colors <- as.numeric(as.factor(cohort_labels[rownames(data_norm)]))

# Panel 1: PCA of the final normalized matrix, colored by cohort. INT path:
# no cohort mean/scale separation is expected by construction. VST path
# (terminal VST, no INT): cohort structure MAY remain visible — expected;
# cohort effects are absorbed by covariates downstream and cohort variance is
# re-verified by the 29b-29f diagnostic gate.
tryCatch({
  pca_final <- prcomp(data_norm, scale. = FALSE, center = TRUE)
  pca_main <- if (use_vst) {
    sprintf("PCA after TMM->VST (%s; cohort structure may remain - expected)", opt_modality)
  } else {
    sprintf("PCA after QN + within-cohort INT (%s)", opt_modality)
  }
  plot(pca_final$x[, 1:2], col = cohort_colors, pch = 19, cex = 0.8,
       xlab = sprintf("PC1 (%.1f%%)", summary(pca_final)$importance[2, 1] * 100),
       ylab = sprintf("PC2 (%.1f%%)", summary(pca_final)$importance[2, 2] * 100),
       main = pca_main)
  legend("topright", legend = levels(as.factor(cohort_labels)),
         col = seq_along(levels(as.factor(cohort_labels))), pch = 19, cex = 0.6)
}, error = function(e) {
  cat(sprintf("  WARN: PCA failed: %s\n", e$message))
  plot.new(); text(0.5, 0.5, "PCA failed")
})

# Panel 2: per-cohort sd of the normalized values. INT path: ~1 in every
# cohort by construction (smaller for tiny cohorts due to INT quantile
# spacing). VST path: NOT equilibrated by construction — this panel is a
# diagnostic, and heterogeneity here is what the 29b-29f gate re-verifies.
tryCatch({
  cohort_sds <- sapply(names(batch_table), function(co) {
    idx <- names(cohort_labels)[cohort_labels == co]
    sd(data_norm[idx, ], na.rm = TRUE)
  })
  sd_main <- if (use_vst) {
    sprintf("Per-cohort scale after VST (%s; not equilibrated by construction)", opt_modality)
  } else {
    sprintf("Per-cohort scale after INT (%s)", opt_modality)
  }
  barplot(cohort_sds, ylab = "sd of normalized values",
          main = sd_main,
          col = seq_along(cohort_sds), las = 2,
          ylim = c(0, max(1.2, max(cohort_sds, na.rm = TRUE) * 1.1)))
  abline(h = 1, lty = 2, col = "grey35")
  legend("topright", legend = names(cohort_sds),
         col = seq_along(cohort_sds), pch = 15, cex = 0.6)
}, error = function(e) {
  cat(sprintf("  WARN: per-cohort sd panel failed: %s\n", e$message))
  plot.new(); text(0.5, 0.5, "Per-cohort sd panel failed")
})

dev.off()
cat(sprintf("  Written: %s\n", diag_path))

cat(sprintf("[%s] Done: %s / %s\n", date(), ancestry, opt_modality))
cat(sprintf("  BED:     %s\n", bed_path))
cat(sprintf("  Groups:  %s\n", groups_path))
cat(sprintf("  Diag:    %s\n", diag_path))
