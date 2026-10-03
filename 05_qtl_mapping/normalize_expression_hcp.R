#!/usr/bin/env Rscript
#
# normalize_expression_hcp.R — per-ancestry expression normalization + HCP
# estimation (formerly combat_normalize_hcp.R)
#
# TWO MODES (--mode):
#
#   vst (DEFAULT) — TMM -> VST on pooled counts, the PsychENCODE/isoTWAS
#     convention (Bhattacharya et al., Nat Genet 2023, PMC10703692):
#     Salmon counts -> tximport -> edgeR TMM -> DESeq2 VST, terminal (NO
#     per-gene INT). The TPM matrix is used ONLY for the devBrain detection
#     filter; the filtered feature set is then read off the pooled counts BED
#     (--counts; per-cohort counts from tximport_counts.R, pooled by
#     pool_expression_within_ancestry.py). Cohort is handled by k-1 dummy
#     covariates: as KNOWN covariates in the HCP Z matrix (appended after QC
#     correlation pruning) and in the tensorQTL covariate model
#     (25_build_covariates.py). NOTE: unlike within-cohort INT, terminal VST
#     does NOT equilibrate cohort scales by construction — cohort variance
#     heterogeneity is re-verified by the 29b-29f diagnostic gate. Writes
#     {ANC}_vst_expression.bed.
#     TMM bridge (do NOT pre-divide counts by TMM and hand to vst() — that is
#     double normalization): sf <- edgeR::calcNormFactors(counts,"TMM") *
#     colSums(counts); sizeFactors(dds) <- sf / exp(mean(log(sf))); vst(dds).
#
#   int — pooled quantile normalization, then per-gene rank-based INT
#     separately within each cohort (within-cohort INT schema, Oct 2026;
#     ComBat DROPPED — the last ComBat-schema commit is
#     ea2b344ab9cd7cdf116b742d97f072d3f7597186, see the 05 README note).
#     Every gene is exactly N(0,1) within each cohort, so cohorts enter the
#     pooled regression on a common scale. Writes {ANC}_int_expression.bed.
#
#   Shared steps:
#   1. Load pooled unnorm expression BED (TPM)
#   2. Gene filter: TPM > 0.1 in > 25% of stratum samples (devBrain §3.3)
#      + no-variance removal
#   3. Normalize (mode-dependent, above)
#   4. Expression-outlier removal: sample connectivity (signed biweight
#      midcorrelation network) z < -3 (devBrain §3.3; pure-R bicor, no WGCNA
#      dependency). Outlier list is written for the isoform modality
#      (script 20 --exclude-samples).
#   5. Write normalized expression BED
#   6. Standardize (center + unit SS) for HCP
#   7. HCP: hcp(Z = QC metrics [+ cohort dummies in vst mode], Y = normalized
#      expression, k, lambda1-3)
#   8. Extract hidden covariates W (n_samples x k) as tensorQTL covariate table
#
# Note (int mode): because INT is applied within cohort, between-cohort mean
# expression differences are removed by construction. Do NOT reuse the
# written BED for between-cohort differential-expression comparisons.
#
# Usage:
#   Rscript normalize_expression_hcp.R \
#       --expression EUR_pooled_expression.bed \
#       --counts EUR_pooled_expression_counts.bed \
#       --mode vst \
#       --qc-metrics all_qc_metrics.tsv \
#       --ancestry-map pooled_sample_ancestry_RNAseq.tsv \
#       --ancestry EUR \
#       --k 15 \
#       --output-dir hcp_output/
#
# Dependencies:
#   Rhcpp   (GitHub: mvaniterson/Rhcpp) — HCP
#   optparse (CRAN) — CLI parsing
#   edgeR, DESeq2 (Bioconductor) — vst mode only (TMM + VST)
#
# Install on seadragon:
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
              help = "Pooled unnorm expression BED file (TPM; drives the detection filter)"),
  make_option("--counts", type = "character", default = NULL,
              help = "Pooled counts BED (from tximport_counts.R via the pooler). REQUIRED when --mode vst."),
  make_option("--mode", type = "character", default = "vst",
              help = "Normalization mode: 'vst' (TMM -> VST on counts, terminal; default) or 'int' (pooled QN + within-cohort INT)"),
  make_option("--qc-metrics", type = "character",
              help = "PicardTools QC metrics TSV (from picard_qc.py)"),
  make_option("--ancestry-map", type = "character",
              help = "Ancestry map TSV: sample_id, assigned_ancestry, cohort"),
  make_option("--ancestry", type = "character",
              help = "Ancestry stratum to process (e.g. EUR, EAS, AFR, HIS)"),
  make_option("--k", type = "integer", default = 15,
              help = "Number of HCP hidden factors (default: 15)"),
  make_option("--lambda1", type = "double", default = 0.5,
              help = "HCP prior strength (default: 0.5, matching scripts 19/25a/25b; the callers pass --lambda1 explicitly, this default only affects standalone use)"),
  make_option("--lambda2", type = "double", default = 1,
              help = "HCP effect regularization (default: 1)"),
  make_option("--lambda3", type = "double", default = 1,
              help = "HCP coefficient regularization (default: 1)"),
  make_option("--output-dir", type = "character", default = ".",
              help = "Output directory"),
  make_option("--nonparametric", action = "store_true", default = FALSE,
              help = "DEPRECATED no-op: ComBat was removed from this script (within-cohort INT schema). Accepted for backward compatibility."),
  make_option("--skip-combat", action = "store_true", default = FALSE,
              help = "DEPRECATED no-op: ComBat was removed from this script (within-cohort INT schema). Accepted for backward compatibility."),
  make_option("--min-samples", type = "integer", default = 10,
              help = "Minimum samples required to run a stratum; smaller strata are skipped gracefully (default: 10)"),
  make_option("--skip-hcp", action = "store_true", default = FALSE,
              help = "Skip HCP factor estimation (normalized expression BED is still written). Use when HCP is run separately (e.g. the 25a optimization module)."),
  make_option("--qc-cor-threshold", type = "double", default = 0.9,
              help = "Drop QC metrics with |correlation| above this threshold to avoid singular Z'Z in HCP (default: 0.9)"),
  make_option("--tpm-min", type = "double", default = 0.1,
              help = "Gene filter: minimum TPM (default: 0.1, devBrain §3.3)"),
  make_option("--tpm-min-prop", type = "double", default = 0.25,
              help = "Gene filter: required fraction of samples above --tpm-min (default: 0.25, devBrain §3.3)"),
  make_option("--outlier-z", type = "double", default = -3,
              help = "Connectivity z-score threshold for expression-outlier removal (default: -3, devBrain §3.3)")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Use bracket notation for hyphenated option names
opt_expression    <- opt[["expression"]]
opt_counts        <- opt[["counts"]]
opt_mode          <- tolower(opt[["mode"]])
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
opt_tpm_min       <- opt[["tpm-min"]]
opt_tpm_min_prop  <- opt[["tpm-min-prop"]]
opt_outlier_z     <- opt[["outlier-z"]]

if (is.null(opt_expression) || is.null(opt_qc_metrics) ||
    is.null(opt_ancestry_map) || is.null(opt_ancestry)) {
  stop("--expression, --qc-metrics, --ancestry-map, and --ancestry are required")
}
if (!opt_mode %in% c("vst", "int")) {
  stop("--mode must be 'vst' or 'int', got: ", opt_mode)
}
if (opt_mode == "vst") {
  if (is.null(opt_counts)) {
    stop("--mode vst requires --counts (pooled counts BED from tximport_counts.R via the pooler)")
  }
  for (pkg in c("edgeR", "DESeq2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("--mode vst requires R package '%s' (not found in .libPaths(): %s)",
                   pkg, paste(.libPaths(), collapse = ", ")))
    }
  }
}

dir.create(opt_output_dir, showWarnings = FALSE, recursive = TRUE)
ancestry <- opt_ancestry
cat(sprintf("[%s] normalize_expression_hcp.R mode=%s for ancestry: %s\n",
            date(), opt_mode, ancestry))

if (opt_nonparametric || opt_skip_combat) {
  cat("  NOTE: --nonparametric/--skip-combat are deprecated no-ops; ComBat was",
      "removed from this pipeline (within-cohort INT schema).\n")
}

# ---- Helper functions ----

# Quantile normalization: normalize samples (rows) to the same distribution.
# Input mat is samples x genes. Transposes to genes x samples, normalizes
# columns (samples), transposes back. Matches PANTRY normalize_phenotypes.py.
quantile_normalize_rows <- function(mat) {
  # Impute NAs with column (gene) medians so that apply(,2,sort) below
  # returns a proper matrix. Without this, sort() drops NAs and columns with
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
  # apply() simplifies to a vector when mat has a single row (single-sample
  # cohort): restore the matrix shape before setting dimnames
  if (!is.matrix(result)) {
    result <- matrix(result, nrow = nrow(mat), ncol = ncol(mat))
  }
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

# TMM -> VST bridge (PsychENCODE/isoTWAS convention).
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

# Biweight midcorrelation between samples (WGCNA::bicor equivalent), pure R.
# mat: samples x genes. Returns the samples x samples bicor matrix.
#   u = (x - median(x)) / (9 * mad(x));  w = (1 - u^2)^2 for |u| < 1, else 0
#   bicor(x,y) = sum((x-mx)*wx * (y-my)*wy) /
#                sqrt( sum(((x-mx)*wx)^2) * sum(((y-my)*wy)^2) )
bicor_samples <- function(mat) {
  X <- as.matrix(mat)  # samples x genes
  med <- apply(X, 1, median, na.rm = TRUE)
  mad1 <- apply(X, 1, function(x) mad(x, constant = 1, na.rm = TRUE))
  mad1[!is.finite(mad1) | mad1 == 0] <- NA
  U <- (X - med) / (9 * mad1)          # length(med) == nrow(X): row-wise
  W <- (1 - U^2)^2
  W[!is.finite(W) | abs(U) >= 1] <- 0
  Y <- (X - med) * W
  Y[!is.finite(Y)] <- 0
  num <- tcrossprod(Y)                 # samples x samples
  den2 <- rowSums(Y^2)
  bc <- num / sqrt(outer(den2, den2))
  bc[!is.finite(bc)] <- 0
  dimnames(bc) <- list(rownames(X), rownames(X))
  bc
}

# Connectivity outlier detection (devBrain §3.3: WGCNA signed network,
# connectivity z < -3). Signed adjacency a = (1 + bicor)/2, power 1;
# connectivity = row sum excluding self.
connectivity_outliers <- function(mat, z_thresh) {
  bc <- bicor_samples(mat)
  adj <- (1 + bc) / 2
  diag(adj) <- 0
  conn <- rowSums(adj, na.rm = TRUE)
  conn_sd <- sd(conn)
  if (!is.finite(conn_sd) || conn_sd == 0) {
    return(list(outliers = character(0), z = rep(0, nrow(mat))))
  }
  z <- (conn - mean(conn)) / conn_sd
  names(z) <- rownames(mat)
  list(outliers = rownames(mat)[z < z_thresh], z = z)
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

# ---- Gene filter: TPM > --tpm-min in > --tpm-min-prop of samples (devBrain
# §3.3), plus no-variance removal ----
detect_frac <- colMeans(expr_data > opt_tpm_min, na.rm = TRUE)
no_var <- apply(expr_data, 2, function(x) length(unique(x[!is.na(x)])) <= 1)
keep_genes <- (detect_frac > opt_tpm_min_prop) & !no_var
cat(sprintf("  Gene filtering: %d -> %d (removed %d with TPM <= %g in >= %d%% of samples, %d no-variance)\n",
            ncol(expr_data), sum(keep_genes),
            sum(detect_frac <= opt_tpm_min_prop), opt_tpm_min,
            as.integer(opt_tpm_min_prop * 100),
            sum(no_var & detect_frac > opt_tpm_min_prop)))
expr_data <- expr_data[, keep_genes, drop = FALSE]

# ---- Resolve cohort labels (needed for within-cohort INT) ----
cohort_labels <- stratum_samples$cohort[match(expr_samples, stratum_samples$sample_id)]
names(cohort_labels) <- expr_samples
if (any(is.na(cohort_labels))) {
  stop("Missing cohort labels for samples: ",
       paste(expr_samples[is.na(cohort_labels)], collapse = ", "))
}

batch_table <- table(cohort_labels)
cat(sprintf("[%s] Cohort structure:\n", date()))
for (b in names(batch_table)) {
  cat(sprintf("    %s: %d samples\n", b, batch_table[b]))
}

# ---- Normalization (mode-dependent) ----
if (opt_mode == "vst") {
  # TMM -> VST on the pooled counts BED (PsychENCODE/isoTWAS convention;
  # terminal VST, no per-gene INT). The TPM matrix above drove the detection
  # filter only; the filtered feature set is read off the counts BED.
  cat(sprintf("[%s] Loading counts BED: %s\n", date(), opt_counts))
  bed_hdr_c <- strsplit(readLines(opt_counts, n = 1, warn = FALSE), "\t", fixed = TRUE)[[1]]
  bed_cc_c <- rep("numeric", length(bed_hdr_c))
  bed_cc_c[1:4] <- c("character", "integer", "integer", "character")
  counts_bed <- read.delim(opt_counts, sep = "\t", check.names = FALSE,
                           colClasses = bed_cc_c)
  missing_samples_c <- setdiff(expr_samples, colnames(counts_bed))
  if (length(missing_samples_c) > 0) {
    stop("Counts BED is missing ", length(missing_samples_c), " stratum samples (e.g. ",
         paste(utils::head(missing_samples_c, 5), collapse = ", "), ")")
  }
  gidx <- match(colnames(expr_data), counts_bed$phenotype_id)
  if (any(is.na(gidx))) {
    stop("Counts BED is missing ", sum(is.na(gidx)), " of the ", ncol(expr_data),
         " TPM-filtered genes (e.g. ",
         paste(utils::head(colnames(expr_data)[is.na(gidx)], 5), collapse = ", "),
         "). tximport_counts.R emits ALL features — regenerate the counts BEDs ",
         "rather than subsetting them upstream.")
  }
  counts_data <- t(as.matrix(counts_bed[gidx, expr_samples, drop = FALSE]))
  rownames(counts_data) <- expr_samples
  colnames(counts_data) <- colnames(expr_data)
  rm(counts_bed)
  # Guard: drop genes with zero-variance (or NA) counts — VST cannot fit a
  # dispersion trend through degenerate rows.
  cv <- apply(counts_data, 2, var)
  zero_var_counts <- !is.finite(cv) | cv == 0
  if (any(zero_var_counts)) {
    cat(sprintf("  Dropping %d genes with zero-variance/NA counts (VST guard)\n",
                sum(zero_var_counts)))
    counts_data <- counts_data[, !zero_var_counts, drop = FALSE]
    expr_data <- expr_data[, !zero_var_counts, drop = FALSE]
  }
  cat(sprintf("[%s] TMM -> VST (edgeR TMM size factors -> DESeq2 vst, blind)\n", date()))
  expr_norm <- tmm_vst(counts_data)
  rm(counts_data)
  cat(sprintf("[%s] Normalized matrix: %d samples x %d genes (TMM -> VST, terminal; no INT)\n",
              date(), nrow(expr_norm), ncol(expr_norm)))
} else {
  # ---- Quantile normalization (pooled across all stratum samples) ----
  cat(sprintf("[%s] Quantile normalization (pooled)\n", date()))
  expr_qn <- quantile_normalize_rows(expr_data)

  # ---- Rank-based INT per gene, separately within each cohort ----
  # Every cohort is transformed on its own ranks regardless of size: each gene
  # is exactly N(0,1) within each cohort, so cohorts enter the pooled scan on a
  # common scale. Within-cohort sample ordering is preserved.
  cat(sprintf("[%s] Rank-based INT per gene within cohort\n", date()))
  expr_int <- matrix(NA_real_, nrow = nrow(expr_qn), ncol = ncol(expr_qn),
                     dimnames = dimnames(expr_qn))
  for (co in names(batch_table)) {
    idx <- names(cohort_labels)[cohort_labels == co]
    expr_int[idx, ] <- inverse_normal_transform(expr_qn[idx, , drop = FALSE])
  }
  cat(sprintf("  INT done: pooled mean=%.4f, sd=%.4f\n",
              mean(expr_int, na.rm = TRUE), sd(expr_int, na.rm = TRUE)))
  for (co in names(batch_table)) {
    idx <- names(cohort_labels)[cohort_labels == co]
    cat(sprintf("  INT within %s (n=%d): mean=%.4f, sd=%.4f (should be ~0, ~1)\n",
                co, length(idx),
                mean(expr_int[idx, ], na.rm = TRUE), sd(expr_int[idx, ], na.rm = TRUE)))
  }
  # ComBat stays dropped: within-cohort INT above already equilibrates
  # per-gene location and scale across cohorts, and a cohort fixed effect in
  # the tensorQTL covariate model (25_build_covariates.py) absorbs remaining
  # between-cohort mean structure.
  expr_norm <- expr_int
  cat(sprintf("[%s] Normalized matrix: %d samples x %d genes (QN + within-cohort INT; no ComBat)\n",
              date(), nrow(expr_norm), ncol(expr_norm)))
}

# ---- Expression-outlier removal: connectivity z < --outlier-z (devBrain §3.3) ----
cat(sprintf("[%s] Connectivity outlier detection (signed bicor network, z < %g)\n",
            date(), opt_outlier_z))
conn <- connectivity_outliers(expr_norm, opt_outlier_z)
outlier_samples <- conn$outliers
outlier_path <- file.path(opt_output_dir, sprintf("%s_expression_outliers.tsv", ancestry))
write.table(data.frame(sample_id = outlier_samples), outlier_path,
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Outliers: %d of %d samples (z < %g) -> %s\n",
            length(outlier_samples), nrow(expr_norm), opt_outlier_z, outlier_path))
if (length(outlier_samples) > 0) {
  cat(sprintf("  Removed: %s\n", paste(outlier_samples, collapse = ", ")))
  keep_samples <- setdiff(rownames(expr_norm), outlier_samples)
  expr_norm <- expr_norm[keep_samples, , drop = FALSE]
  expr_samples <- keep_samples
  cohort_labels <- cohort_labels[expr_samples]
}

# ---- Write normalized expression BED ----
# Match gene metadata from the original BED by phenotype_id.
gene_idx <- match(colnames(expr_norm), expr_bed$phenotype_id)
expr_out <- data.frame(
  chr = expr_bed[["#chr"]][gene_idx],
  start = expr_bed$start[gene_idx],
  end = expr_bed$end[gene_idx],
  phenotype_id = colnames(expr_norm),
  t(expr_norm),
  check.names = FALSE
)
colnames(expr_out)[1] <- "#chr"
expr_out_path <- file.path(opt_output_dir,
                           sprintf("%s_%s_expression.bed", ancestry, opt_mode))
write.table(expr_out, expr_out_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Written: %s (%d genes x %d samples)\n",
            expr_out_path, nrow(expr_out), ncol(expr_out) - 4))

# ---- Skip HCP if requested (--skip-hcp) ----
if (opt_skip_hcp) {
  cat(sprintf("\n[%s] --skip-hcp: skipping HCP estimation. Normalized expression BED written.\n", date()))
  cat(sprintf("  Run the 25a optimization module / peer_factors.py separately to extract latent factors.\n"))
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

# In VST mode, append k-1 cohort dummies as KNOWN covariates for HCP
# (PsychENCODE/isoTWAS design: cohort handled by dummies in the mapping model
# AND as known covariates in HCP estimation). Appended AFTER the QC
# correlation pruning so the dummies are protected from it. Reference level =
# the largest cohort. In INT mode no dummies are needed: within-cohort INT
# already removes cohort mean/scale structure by construction.
if (opt_mode == "vst") {
  co <- factor(cohort_labels[expr_samples])
  ref_cohort <- names(sort(table(co), decreasing = TRUE))[1]
  co <- relevel(co, ref = ref_cohort)
  if (nlevels(co) > 1) {
    dummies <- model.matrix(~ co)[, -1, drop = FALSE]
    colnames(dummies) <- sub("^co", "cohort_", colnames(dummies))
    rownames(dummies) <- expr_samples
    cat(sprintf("  HCP known covariates: + %d cohort dummies (reference: %s)\n",
                ncol(dummies), ref_cohort))
    qc_subset <- cbind(qc_subset, dummies)
  } else {
    cat("  HCP known covariates: single-cohort stratum — no cohort dummies added\n")
  }
}

# Standardize for HCP (center + unit SS)
qc_std <- standardize_for_hcp(qc_subset)
expr_std <- standardize_for_hcp(expr_norm)

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

# 3. Connectivity z-scores (outlier diagnostic)
plot(sort(conn$z), pch = 19, cex = 0.6,
     xlab = "Sample (rank)", ylab = "Connectivity z-score",
     main = "Sample connectivity (outliers removed)")
abline(h = opt_outlier_z, col = "red", lty = 2)

# 4. PCA of the final normalized matrix, colored by cohort. In INT mode, no
# cohort mean or scale separation is expected by construction. In VST mode
# (terminal VST, no INT), cohort structure MAY remain visible — expected:
# VST does not equilibrate cohort scales; cohort effects are absorbed by the
# cohort dummies in HCP Z and in the tensorQTL covariate model, and cohort
# variance heterogeneity is re-verified by the 29b-29f diagnostic gate.
pca_final <- prcomp(expr_norm, scale. = FALSE, center = TRUE)
cohort_colors <- as.numeric(as.factor(cohort_labels[rownames(expr_norm)]))

pca_main <- if (opt_mode == "vst") {
  "PCA after TMM->VST (final; cohort structure may remain - expected)"
} else {
  "PCA after QN + within-cohort INT (final)"
}
plot(pca_final$x[, 1:2], col = cohort_colors, pch = 19, cex = 0.8,
     xlab = sprintf("PC1 (%.1f%%)", summary(pca_final)$importance[2, 1] * 100),
     ylab = sprintf("PC2 (%.1f%%)", summary(pca_final)$importance[2, 2] * 100),
     main = pca_main)
legend("topright", legend = levels(as.factor(cohort_labels)),
       col = seq_along(levels(as.factor(cohort_labels))), pch = 19, cex = 0.6)

dev.off()
cat(sprintf("  Written: %s\n", diag_path))

cat(sprintf("[%s] Done: ancestry %s\n", date(), ancestry))
cat(sprintf("  HCP factors: %s\n", hcp_out_path))
cat(sprintf("  Expression: %s\n", expr_out_path))
cat(sprintf("  Outliers: %s\n", outlier_path))
cat(sprintf("  Diagnostics: %s\n", diag_path))
