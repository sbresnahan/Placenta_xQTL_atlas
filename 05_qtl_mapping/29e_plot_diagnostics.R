#!/usr/bin/env Rscript
# =============================================================================
# 29e_plot_diagnostics.R — Diagnostic figures + validation summary for the
# within-cohort-INT schema gate (module 05; runs after 29b-29d, gates 31).
#
# Revision of the three ad-hoc diagnostic plotters (expression heterogeneity,
# Choi-vs-pooled gene-top/exact-pair, Choi significant retention), generalized
# over ancestry and cohort labels, plus the baseline-vs-rerun validation
# summary that gates fine-mapping.
#
# All inputs are optional individually; sections skip gracefully. Per ancestry:
#   --pcs            {ANC}_expression_sample_PCs.tsv                    (29b)
#   --variance       {ANC}_gene_residual_variance_by_cohort.tsv.gz      (29b)
#   --hcp            {ANC}_expression_hcp15_vs_hcp45_exact_pairs.tsv.gz (29b)
#   --gene-top       Choi_vs_pooled_{ANC}.gene_top_comparison.tsv.gz    (29c)
#   --exact-pairs    Choi_vs_pooled_{ANC}.exact_{ANC}_lead_comparison.tsv.gz
#   --retention      Choi_vs_pooled_{ANC}.Choi_significant_retention.tsv.gz
#   --catalog        Choi_vs_pooled_{ANC}.common_tested_catalog.tsv.gz  (29c)
#   --heterogeneity  {ANC}_expression_cohort_heterogeneity.tsv.gz       (29d)
#   --extra-pairs    {ANC}_expression_extra_pairs_withincohort.tsv.gz   (29d)
#   --top-table      {ANC}_expression_cisqtl_top.tsv                    (29)
#   --ancestry       label (default EAS)
#   --cohort-a / --cohort-b  comparison cohorts (default: two largest)
#   --baseline-summary  archived diagnostic_summary_stats.tsv (pre-fix run)
#   --baseline-metrics  archived baseline_metrics.tsv (metric/value; pre-fix
#                       retention, within-cohort, I2/Q/GxC, eGene baselines)
#   --fdr            FDR threshold (default 0.05)
#   --outdir         output directory
#
# Outputs: per-section PNG/PDF figures, diagnostic_summary_stats.tsv (same
# metric names as the archived baseline where cohorts match), and
# validation_summary.tsv (metric, baseline, rerun, target, pass).
# =============================================================================

.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value after ", flag)
  args[[i + 1]]
}

ancestry   <- get_arg("--ancestry", "EAS")
pcs_file   <- get_arg("--pcs")
var_file   <- get_arg("--variance")
hcp_file   <- get_arg("--hcp")
top_file   <- get_arg("--gene-top")
exact_file <- get_arg("--exact-pairs")
ret_file   <- get_arg("--retention")
cat_file   <- get_arg("--catalog")
het_file   <- get_arg("--heterogeneity")
extra_file <- get_arg("--extra-pairs")
top_table  <- get_arg("--top-table")
base_sum_f <- get_arg("--baseline-summary")
base_met_f <- get_arg("--baseline-metrics")
cohort_a   <- get_arg("--cohort-a")
cohort_b   <- get_arg("--cohort-b")
fdr        <- as.numeric(get_arg("--fdr", "0.05"))
outdir     <- get_arg("--outdir", ".")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------- style --------------------------------------

base_family <- tryCatch({
  if (requireNamespace("systemfonts", quietly = TRUE) &&
      "Liberation Sans" %in% systemfonts::system_fonts()$family)
    "Liberation Sans" else ""
}, error = function(e) "")

theme_diag <- function(base_size = 12) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size),
      axis.title = element_text(face = "bold"),
      legend.title = element_blank(),
      legend.position = "right",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

save_plot <- function(plot, stem, width = 7, height = 5.5) {
  ggsave(file.path(outdir, paste0(stem, ".png")), plot,
         width = width, height = height, units = "in", dpi = 300,
         bg = "white")
  ggsave(file.path(outdir, paste0(stem, ".pdf")), plot,
         width = width, height = height, units = "in",
         device = cairo_pdf, bg = "white")
}

known_colors <- c("GUSTO" = "#0072B2", "SNUH" = "#D55E00",
                  "NIEHS_RICHS" = "#009E73", "NIGMS" = "#CC79A7")
cohort_palette <- function(levels) {
  cols <- known_colors[levels]
  missing <- is.na(cols)
  if (any(missing)) {
    cols[missing] <- scales::hue_pal()(sum(missing))
  }
  names(cols) <- levels
  cols
}

# metric collectors
summary_metrics <- list()   # diagnostic_summary_stats.tsv
add_metric <- function(name, value) {
  summary_metrics[[name]] <<- value
}
validation <- list()        # validation_summary.tsv rows
add_validation <- function(metric, baseline, rerun, target, pass) {
  validation[[length(validation) + 1]] <<- data.frame(
    metric = metric, baseline = baseline, rerun = rerun,
    target = target, pass = pass, stringsAsFactors = FALSE)
}

read_gz_tsv <- function(f) read_tsv(f, show_col_types = FALSE,
                                    name_repair = "unique")

# ============================================================================
# Section 1: expression heterogeneity (29b inputs)
# ============================================================================

if (!is.null(pcs_file) && !is.null(var_file) &&
    file.exists(pcs_file) && file.exists(var_file)) {
  cat("\n== Section 1: expression heterogeneity ==\n")

  pcs <- read_gz_tsv(pcs_file)
  if (!"sample_id" %in% names(pcs)) names(pcs)[1] <- "sample_id"
  stopifnot(all(c("sample_id", "PC1", "PC2", "cohort") %in% names(pcs)))

  var_df <- read_gz_tsv(var_file)
  stopifnot(all(c("phenotype_id", "max_min_variance_ratio",
                  "log2_variance_ratio") %in% names(var_df)))

  # comparison cohorts: explicit flags, else the two largest in the PCs file
  cohort_counts <- pcs %>% count(cohort, name = "n") %>% arrange(desc(n))
  if (is.null(cohort_a)) cohort_a <- as.character(cohort_counts$cohort[1])
  if (is.null(cohort_b)) cohort_b <- as.character(cohort_counts$cohort[2])
  cat("  comparison cohorts:", cohort_a, "vs", cohort_b, "\n")
  col_a <- paste0("var_", cohort_a)
  col_b <- paste0("var_", cohort_b)
  stopifnot(all(c(col_a, col_b) %in% names(var_df)))

  cohort_levels <- c(cohort_a, cohort_b,
                     setdiff(cohort_counts$cohort, c(cohort_a, cohort_b)))
  pcs <- pcs %>%
    mutate(cohort = factor(cohort, levels = cohort_levels))
  ccolors <- cohort_palette(cohort_levels)

  var_df <- var_df %>%
    mutate(ratio_ba = .data[[col_b]] / .data[[col_a]],
           log2_ratio_ba = log2(ratio_ba))

  pct_b_gt_a   <- mean(var_df$ratio_ba > 1, na.rm = TRUE) * 100
  pct_b_gt_2x  <- mean(var_df$ratio_ba > 2, na.rm = TRUE) * 100
  median_ba    <- median(var_df$ratio_ba, na.rm = TRUE)
  p90_ba       <- quantile(var_df$ratio_ba, 0.90, na.rm = TRUE)
  p95_ba       <- quantile(var_df$ratio_ba, 0.95, na.rm = TRUE)
  allcoh_med   <- median(var_df$max_min_variance_ratio, na.rm = TRUE)
  allcoh_p90   <- quantile(var_df$max_min_variance_ratio, 0.90, na.rm = TRUE)
  allcoh_p95   <- quantile(var_df$max_min_variance_ratio, 0.95, na.rm = TRUE)

  add_metric("n_genes_residual_variance", nrow(var_df))
  add_metric(sprintf("fraction_%s_variance_gt_%s", cohort_b, cohort_a),
             pct_b_gt_a / 100)
  add_metric(sprintf("fraction_%s_variance_gt_2x_%s", cohort_b, cohort_a),
             pct_b_gt_2x / 100)
  add_metric(sprintf("median_%s_to_%s_variance_ratio", cohort_b, cohort_a),
             median_ba)
  add_metric(sprintf("p90_%s_to_%s_variance_ratio", cohort_b, cohort_a),
             p90_ba)
  add_metric(sprintf("p95_%s_to_%s_variance_ratio", cohort_b, cohort_a),
             p95_ba)
  add_metric("median_allcohort_max_min_variance_ratio", allcoh_med)
  add_metric("p90_allcohort_max_min_variance_ratio", allcoh_p90)
  add_metric("p95_allcohort_max_min_variance_ratio", allcoh_p95)

  # Fig 1: PC1 vs PC2 by cohort
  major_df <- pcs %>% filter(cohort %in% c(cohort_a, cohort_b))
  p1 <- ggplot(pcs, aes(PC1, PC2, color = cohort)) +
    stat_ellipse(data = major_df, aes(group = cohort), type = "norm",
                 level = 0.80, linewidth = 0.8, alpha = 0.9,
                 show.legend = FALSE) +
    geom_point(size = 2.1, alpha = 0.65) +
    scale_color_manual(values = ccolors, drop = FALSE) +
    labs(title = sprintf("Expression PC1 vs PC2 by cohort (%s)", ancestry),
         subtitle = paste(paste0(as.character(cohort_counts$cohort),
                                 " n=", cohort_counts$n),
                          collapse = "   |   "),
         x = "Expression PC1", y = "Expression PC2") +
    theme_diag()
  save_plot(p1, "01_expression_PC1_PC2_by_cohort", 7.2, 5.6)

  # Fig 2: per-cohort residual variance scatter
  scatter_df <- var_df %>%
    filter(is.finite(.data[[col_a]]), is.finite(.data[[col_b]]),
           .data[[col_a]] > 0, .data[[col_b]] > 0)
  variance_lim <- range(c(scatter_df[[col_a]], scatter_df[[col_b]]),
                        finite = TRUE)
  p2 <- ggplot(scatter_df, aes(x = .data[[col_a]], y = .data[[col_b]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.8, color = "grey35") +
    geom_point(alpha = 0.20, size = 1.2) +
    scale_x_log10(labels = label_number(accuracy = 0.01)) +
    scale_y_log10(labels = label_number(accuracy = 0.01)) +
    coord_equal(xlim = variance_lim, ylim = variance_lim) +
    labs(title = sprintf("Per-cohort residual expression variance (%s)",
                         ancestry),
         subtitle = sprintf("%s genes; dashed line = equal variance",
                            comma(nrow(scatter_df))),
         x = sprintf("Residual variance in %s (log10)", cohort_a),
         y = sprintf("Residual variance in %s (log10)", cohort_b)) +
    theme_diag()
  save_plot(p2, "02_residual_variance_by_cohort", 9.0, 6.1)

  # Fig 3: variance-ratio distribution
  ratio_df <- var_df %>% filter(is.finite(log2_ratio_ba))
  median_log2 <- median(ratio_df$log2_ratio_ba, na.rm = TRUE)
  p3 <- ggplot(ratio_df, aes(x = log2_ratio_ba)) +
    geom_histogram(bins = 60, boundary = 0, linewidth = 0.2,
                   color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8,
               color = "grey35") +
    geom_vline(xintercept = median_log2, linewidth = 0.9) +
    scale_x_continuous(breaks = log2(c(0.25, 0.5, 1, 2, 4, 8)),
                       labels = c("0.25x", "0.5x", "1x", "2x", "4x", "8x")) +
    labs(title = sprintf("%s:%s residual-variance ratio (%s)",
                         cohort_b, cohort_a, ancestry),
         subtitle = "Per-gene ratio of covariate-residual variance",
         x = sprintf("%s:%s residual-variance ratio (log2)",
                     cohort_b, cohort_a),
         y = "Genes") +
    theme_diag()
  save_plot(p3, "03_residual_variance_ratio_distribution", 7.2, 5.4)

  # Fig 4: HCP k-low vs k-high exact pairs
  if (!is.null(hcp_file) && file.exists(hcp_file)) {
    hcp <- read_gz_tsv(hcp_file)
    k_cols <- grep("^abs_t_k", names(hcp), value = TRUE)
    lead_cols <- grep("^lead_at_k", names(hcp), value = TRUE)
    stopifnot(length(k_cols) == 2, length(lead_cols) == 2)
    k_low <- sub("abs_t_k", "", k_cols[1])
    k_high <- sub("abs_t_k", "", k_cols[2])
    hcp <- hcp %>%
      filter(is.finite(.data[[k_cols[1]]]), is.finite(.data[[k_cols[2]]])) %>%
      mutate(lead_status = case_when(
        .data[[lead_cols[1]]] & .data[[lead_cols[2]]] ~ "Lead at both",
        .data[[lead_cols[1]]] ~ paste0("Lead at k=", k_low),
        .data[[lead_cols[2]]] ~ paste0("Lead at k=", k_high),
        TRUE ~ "Other"),
        lead_status = factor(lead_status,
                             levels = c("Lead at both",
                                        paste0("Lead at k=", k_low),
                                        paste0("Lead at k=", k_high),
                                        "Other")))
    median_tlow  <- median(hcp[[k_cols[1]]], na.rm = TRUE)
    median_thigh <- median(hcp[[k_cols[2]]], na.rm = TRUE)
    pct_weaker   <- mean(hcp[[k_cols[2]]] < hcp[[k_cols[1]]],
                         na.rm = TRUE) * 100
    wilcox_p <- suppressWarnings(
      wilcox.test(hcp[[k_cols[1]]], hcp[[k_cols[2]]], paired = TRUE,
                  exact = FALSE)$p.value)

    add_metric("n_HCP_exact_pairs", nrow(hcp))
    add_metric(sprintf("median_abs_t_k%s", k_low), median_tlow)
    add_metric(sprintf("median_abs_t_k%s", k_high), median_thigh)
    add_metric(sprintf("fraction_weaker_at_k%s", k_high), pct_weaker / 100)
    add_metric("paired_wilcoxon_p_HCP", wilcox_p)

    lead_colors <- c("Lead at both" = "#000000",
                     setNames("#0072B2", paste0("Lead at k=", k_low)),
                     setNames("#D55E00", paste0("Lead at k=", k_high)),
                     "Other" = "grey60")
    hcp_lim <- range(c(hcp[[k_cols[1]]], hcp[[k_cols[2]]]), finite = TRUE)
    p4 <- ggplot(hcp, aes(x = .data[[k_cols[1]]], y = .data[[k_cols[2]]],
                          color = lead_status)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  linewidth = 0.8, color = "grey35") +
      geom_point(alpha = 0.35, size = 1.5) +
      scale_color_manual(values = lead_colors, drop = FALSE) +
      coord_equal(xlim = hcp_lim, ylim = hcp_lim) +
      labs(title = sprintf("Exact-pair signal at HCP k=%s vs k=%s (%s)",
                           k_low, k_high, ancestry),
           subtitle = sprintf("%s chr1 SNP-gene pairs refit under both k",
                              comma(nrow(hcp))),
           x = sprintf("|t| at HCP k=%s", k_low),
           y = sprintf("|t| at HCP k=%s", k_high)) +
      theme_diag()
    save_plot(p4, sprintf("04_HCP_k%s_vs_k%s_exact_pair_signal",
                          k_low, k_high), 7.2, 6.0)
  } else {
    cat("  (HCP k-comparison input not provided; panel 4 skipped)\n")
  }
} else {
  cat("\n== Section 1 skipped: --pcs/--variance not provided ==\n")
}

# Fallback cohort resolution for later sections when section 1 did not run:
# infer from the heterogeneity file's beta_<cohort> columns.
if (is.null(cohort_a) || is.null(cohort_b)) {
  if (!is.null(het_file) && file.exists(het_file)) {
    bcols <- grep("^beta_", names(read_gz_tsv(het_file)), value = TRUE)
    bcols <- setdiff(sub("^beta_", "", bcols), c("meta", "gxc"))
    if (length(bcols) >= 2) {
      cohort_a <- bcols[1]
      cohort_b <- bcols[2]
      cat("  cohorts inferred from heterogeneity file:", cohort_a, "vs",
          cohort_b, "\n")
    }
  }
}
if (is.null(cohort_a)) cohort_a <- "GUSTO"
if (is.null(cohort_b)) cohort_b <- "SNUH"

# ============================================================================
# Section 2: external-reference (Choi) comparison (29c inputs)
# ============================================================================

if (!is.null(top_file) && !is.null(exact_file) &&
    file.exists(top_file) && file.exists(exact_file)) {
  cat("\n== Section 2: Choi vs pooled", ancestry, "==\n")

  top <- read_gz_tsv(top_file)
  exact <- read_gz_tsv(exact_file)
  t_ours   <- paste0("t_", ancestry)
  abs_ours <- paste0("abs_t_", ancestry)
  ml_ours  <- paste0("mlogp_", ancestry)
  stopifnot(all(c(t_ours, abs_ours, ml_ours) %in% names(top)))
  stopifnot(all(c(t_ours, abs_ours) %in% names(exact)))

  # Fig C1: per-gene best-variant signal
  top2 <- top[is.finite(top[[ml_ours]]) & is.finite(top$mlogp), ]
  rho1 <- cor(top2$mlogp, top2[[ml_ours]], method = "spearman",
              use = "complete.obs")
  top2$same_lead_variant <- factor(as.character(top2$same_lead_variant),
                                   levels = c("TRUE", "FALSE"))
  pc1 <- ggplot(top2, aes(x = mlogp, y = .data[[ml_ours]],
                          color = same_lead_variant)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey45") +
    geom_point(alpha = 0.55, size = 1.7) +
    scale_color_manual(values = c("TRUE" = "#0072B2", "FALSE" = "grey65"),
                       labels = c("TRUE" = "Same lead SNP",
                                  "FALSE" = "Different lead SNP")) +
    labs(title = sprintf("Gene-level cis-eQTL signal: Choi SNUH vs pooled %s",
                         ancestry),
         subtitle = sprintf("Each axis uses that study's strongest cis variant; Spearman rho = %.2f",
                            rho1),
         x = expression("Choi best cis signal  (" * -log[10](P) * ")"),
         y = sprintf("Pooled %s best cis signal  (-log10 P)", ancestry)) +
    theme_diag()
  save_plot(pc1, "05_gene_level_signal_Choi_vs_pooled", 7, 5.8)

  # Fig C2: exact-pair signed t
  ex <- exact[is.finite(exact[[t_ours]]) & is.finite(exact$t_Choi), ]
  rho2 <- cor(ex$t_Choi, ex[[t_ours]], method = "spearman",
              use = "complete.obs")
  dir_agree <- mean(sign(ex$t_Choi) == sign(ex[[t_ours]]), na.rm = TRUE)
  lim <- max(abs(c(ex$t_Choi, ex[[t_ours]])), na.rm = TRUE)
  pc2 <- ggplot(ex, aes(x = t_Choi, y = .data[[t_ours]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey45") +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey75") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey75") +
    geom_point(alpha = 0.45, size = 1.6) +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(title = "Exact same SNP-gene effects",
         subtitle = sprintf("Pooled-%s lead SNP evaluated in Choi; Spearman rho = %.2f; direction agreement = %.1f%%",
                            ancestry, rho2, 100 * dir_agree),
         x = expression("Choi " * beta / SE),
         y = sprintf("Pooled %s beta/SE", ancestry)) +
    theme_diag()
  save_plot(pc2, "06_exact_pair_signed_t_Choi_vs_pooled", 6.5, 6.0)

  # Fig C3: exact-pair absolute strength
  rho3 <- cor(ex$abs_t_Choi, ex[[abs_ours]], method = "spearman",
              use = "complete.obs")
  lim_abs <- max(c(ex$abs_t_Choi, ex[[abs_ours]]), na.rm = TRUE)
  pc3 <- ggplot(ex, aes(x = abs_t_Choi, y = .data[[abs_ours]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey45") +
    geom_point(alpha = 0.45, size = 1.6) +
    coord_equal(xlim = c(0, lim_abs), ylim = c(0, lim_abs)) +
    labs(title = "Exact-pair association strength",
         subtitle = sprintf("Spearman rho = %.2f; points below the line are attenuated in pooled %s",
                            rho3, ancestry),
         x = expression("Choi |" * beta / SE * "|"),
         y = sprintf("Pooled %s |beta/SE|", ancestry)) +
    theme_diag()
  save_plot(pc3, "07_exact_pair_abs_t_Choi_vs_pooled", 6.5, 6.0)

  # Fig C4: attenuation distribution
  ex$log2_t_ratio <- log2(ex[[abs_ours]] / ex$abs_t_Choi)
  med_ratio <- median(ex[[abs_ours]] / ex$abs_t_Choi, na.rm = TRUE)
  frac_weaker <- mean(ex[[abs_ours]] < ex$abs_t_Choi, na.rm = TRUE)
  pc4 <- ggplot(ex[is.finite(ex$log2_t_ratio), ], aes(x = log2_t_ratio)) +
    geom_histogram(bins = 60, color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
    labs(title = sprintf("Exact-pair signal ratio, pooled %s vs Choi",
                         ancestry),
         subtitle = sprintf("Median |t| ratio = %.2f; %.1f%% weaker in pooled %s",
                            med_ratio, 100 * frac_weaker, ancestry),
         x = sprintf("log2(|t| %s / |t| Choi)", ancestry),
         y = "SNP-gene pairs") +
    theme_diag()
  save_plot(pc4, "08_exact_pair_attenuation_distribution", 7, 5.2)

  add_metric("choi_exact_pairs", nrow(ex))
  add_metric("choi_exact_pair_direction_agreement", dir_agree)
  add_metric("choi_exact_pair_spearman_t", rho2)
  add_metric("choi_exact_pair_median_abs_t_ratio", med_ratio)
} else {
  cat("\n== Section 2 skipped: --gene-top/--exact-pairs not provided ==\n")
}

# ---- retention panels ----
retention_fraction <- NA_real_
lost_genes <- character(0)
if (!is.null(ret_file) && !is.null(cat_file) &&
    file.exists(ret_file) && file.exists(cat_file)) {
  cat("\n== Section 2b: Choi significant retention ==\n")
  d <- read_gz_tsv(ret_file)
  catg <- read_gz_tsv(cat_file)
  ml_choi <- "mlog10_pbeta_Choi"
  ml_ours <- paste0("mlog10_pbeta_", ancestry)
  stopifnot(all(c(ml_choi, ml_ours, "retention") %in% names(d)))

  retained_lbl <- paste0("Retained in pooled ", ancestry)
  lost_lbl <- paste0("Lost in pooled ", ancestry)
  retention_fraction <- mean(d$retention == retained_lbl)
  lost_genes <- d$gene_base[d$retention == lost_lbl]
  add_metric("choi_significant_tested_both", nrow(d))
  add_metric("choi_retention_fraction", retention_fraction)

  x <- d[is.finite(d[[ml_choi]]) & is.finite(d[[ml_ours]]), ]
  rho <- cor(x[[ml_choi]], x[[ml_ours]], method = "spearman",
             use = "complete.obs")
  ret_cols <- setNames(c("#0072B2", "#D55E00"), c(retained_lbl, lost_lbl))
  pr1 <- ggplot(x, aes(x = .data[[ml_choi]], y = .data[[ml_ours]],
                       color = retention)) +
    geom_point(alpha = 0.55, size = 1.7) +
    scale_color_manual(values = ret_cols) +
    labs(title = sprintf("Choi-significant eGenes in pooled %s", ancestry),
         subtitle = sprintf("Choi-significant genes tested in both analyses; Spearman rho = %.2f",
                            rho),
         x = expression("Choi SNUH gene-level  " * -log[10](P[beta])),
         y = parse(text = sprintf('"Pooled %s gene-level  " * -log[10](P[beta])',
                                  ancestry))) +
    theme_diag()
  save_plot(pr1, "09_Choi_eGenes_signal_retention", 7.2, 5.8)

  counts <- catg %>% filter(category != "Neither significant") %>%
    count(category)
  cat_levels <- c("Shared significant", "Choi-only significant",
                  sprintf("Pooled-%s-only significant", ancestry))
  counts$category <- factor(counts$category, levels = cat_levels)
  pr2 <- ggplot(counts, aes(category, n)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = comma(n)), vjust = -0.4, size = 4) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = "Expression eGene overlap",
         subtitle = sprintf("Genes tested in both Choi and pooled %s",
                            ancestry),
         x = NULL, y = "eGenes") +
    theme_diag() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
  save_plot(pr2, "10_eGene_overlap_counts", 6.8, 5.2)
} else {
  cat("\n== Section 2b skipped: --retention/--catalog not provided ==\n")
}

# ============================================================================
# Section 3: cohort heterogeneity (29d input)
# ============================================================================

median_i2 <- NA_real_
n_q_fdr <- NA_integer_
n_gxc_fdr <- NA_integer_
if (!is.null(het_file) && file.exists(het_file)) {
  cat("\n== Section 3: cohort heterogeneity ==\n")
  het <- read_gz_tsv(het_file)
  stopifnot(all(c("i2", "q_cochran_q") %in% names(het)))
  median_i2 <- median(het$i2, na.rm = TRUE)
  n_q_fdr <- sum(het$q_cochran_q <= fdr, na.rm = TRUE)
  n_gxc_fdr <- if ("q_gxc" %in% names(het))
    sum(het$q_gxc <= fdr, na.rm = TRUE) else NA_integer_
  add_metric("heterogeneity_n_pairs", nrow(het))
  add_metric("median_i2", median_i2)
  add_metric("n_cochran_q_fdr", n_q_fdr)
  add_metric("n_gxc_fdr", n_gxc_fdr)

  ph1 <- ggplot(het %>% filter(is.finite(i2)), aes(x = i2)) +
    geom_histogram(bins = 50, boundary = 0, linewidth = 0.2,
                   color = "white") +
    geom_vline(xintercept = median_i2, linewidth = 0.9) +
    labs(title = sprintf("Cross-cohort effect heterogeneity (%s)", ancestry),
         subtitle = sprintf("%s pooled-significant pairs; median I2 = %.2f; %d Cochran-Q and %s GxC FDR<%.2f hits",
                            comma(nrow(het)), median_i2, n_q_fdr,
                            ifelse(is.na(n_gxc_fdr), "NA",
                                   as.character(n_gxc_fdr)), fdr),
         x = "I2 (inconsistency across cohorts)", y = "Pairs") +
    theme_diag()
  save_plot(ph1, "11_cohort_heterogeneity_i2", 7.2, 5.0)
} else {
  cat("\n== Section 3 skipped: --heterogeneity not provided ==\n")
}

# ============================================================================
# Validation summary (baseline vs rerun)
# ============================================================================

cat("\n== Validation summary ==\n")

baseline <- list()
if (!is.null(base_sum_f) && file.exists(base_sum_f)) {
  bs <- read_gz_tsv(base_sum_f)
  baseline <- as.list(setNames(bs$value, bs$metric))
  cat("  loaded", length(baseline), "baseline metrics from", base_sum_f,
      "\n")
}
if (!is.null(base_met_f) && file.exists(base_met_f)) {
  bm <- read_gz_tsv(base_met_f)
  baseline <- modifyList(baseline, as.list(setNames(bm$value, bm$metric)))
  cat("  loaded", length(bm$metric), "baseline metrics from", base_met_f,
      "\n")
}
get_baseline <- function(name) {
  v <- baseline[[name]]
  if (is.null(v)) NA_real_ else as.numeric(v)
}

# variance-ratio targets
m_ratio <- sprintf("median_%s_to_%s_variance_ratio", cohort_b, cohort_a)
m_frac <- sprintf("fraction_%s_variance_gt_%s", cohort_b, cohort_a)
rerun_ratio <- as.numeric(summary_metrics[[m_ratio]] %||% NA)
rerun_frac <- as.numeric(summary_metrics[[m_frac]] %||% NA)
rerun_p90 <- as.numeric(
  summary_metrics[["p90_allcohort_max_min_variance_ratio"]] %||% NA)
add_validation(m_ratio, get_baseline(m_ratio), rerun_ratio,
               "[0.9, 1.1]",
               ifelse(is.na(rerun_ratio), NA,
                      rerun_ratio >= 0.9 && rerun_ratio <= 1.1))
add_validation(m_frac, get_baseline(m_frac), rerun_frac, "[0.4, 0.6]",
               ifelse(is.na(rerun_frac), NA,
                      rerun_frac >= 0.4 && rerun_frac <= 0.6))
add_validation("p90_allcohort_max_min_variance_ratio",
               get_baseline("p90_allcohort_max_min_variance_ratio"),
               rerun_p90, "<= 1.5",
               ifelse(is.na(rerun_p90), NA, rerun_p90 <= 1.5))

# Choi retention
add_validation("choi_retention_fraction",
               get_baseline("choi_retention_fraction"), retention_fraction,
               "> baseline",
               ifelse(is.na(retention_fraction) |
                      is.na(get_baseline("choi_retention_fraction")), NA,
                      retention_fraction >
                      get_baseline("choi_retention_fraction")))

# within-cohort-B nominal-P<0.05 fraction on lost eGenes (29d pair set B)
within_b_lost <- NA_real_
if (!is.null(extra_file) && file.exists(extra_file) &&
    length(lost_genes) > 0) {
  ex <- read_gz_tsv(extra_file)
  pcol_b <- paste0("p_", cohort_b)
  if (pcol_b %in% names(ex)) {
    # lost_genes are version-stripped base ids (29c); match on stripped ids
    ex_lost <- ex[sub("\\..*$", "", ex$phenotype_id) %in% lost_genes, ]
    within_b_lost <- mean(ex_lost[[pcol_b]] < 0.05, na.rm = TRUE)
    add_metric(sprintf("within_%s_lost_nominal_p05", cohort_b),
               within_b_lost)
  }
}
add_validation(sprintf("within_%s_lost_nominal_p05", cohort_b),
               get_baseline(sprintf("within_%s_lost_nominal_p05",
                                    cohort_b)),
               within_b_lost, "> baseline",
               ifelse(is.na(within_b_lost) |
                      is.na(get_baseline(sprintf(
                        "within_%s_lost_nominal_p05", cohort_b))), NA,
                      within_b_lost > get_baseline(sprintf(
                        "within_%s_lost_nominal_p05", cohort_b))))

# heterogeneity
add_validation("median_i2", get_baseline("median_i2"), median_i2,
               "< baseline",
               ifelse(is.na(median_i2) | is.na(get_baseline("median_i2")),
                      NA, median_i2 < get_baseline("median_i2")))
add_validation("n_cochran_q_fdr", get_baseline("n_cochran_q_fdr"), n_q_fdr,
               "< baseline",
               ifelse(is.na(n_q_fdr) | is.na(get_baseline("n_cochran_q_fdr")),
                      NA, n_q_fdr < get_baseline("n_cochran_q_fdr")))
add_validation("n_gxc_fdr", get_baseline("n_gxc_fdr"), n_gxc_fdr,
               "< baseline",
               ifelse(is.na(n_gxc_fdr) | is.na(get_baseline("n_gxc_fdr")),
                      NA, n_gxc_fdr < get_baseline("n_gxc_fdr")))

# discovery + calibration from the top table
if (!is.null(top_table) && file.exists(top_table)) {
  tt <- read_gz_tsv(top_table)
  if ("qval" %in% names(tt)) {
    n_eg <- sum(tt$qval <= fdr, na.rm = TRUE)
    add_metric("n_egenes", n_eg)
    add_validation("n_egenes", get_baseline("n_egenes"), n_eg,
                   "> baseline",
                   ifelse(is.na(get_baseline("n_egenes")), NA,
                          n_eg > get_baseline("n_egenes")))
  }
  if ("pval_nominal" %in% names(tt)) {
    pv <- tt$pval_nominal[is.finite(tt$pval_nominal)]
    lambda <- median(qchisq(pv, df = 1, lower.tail = FALSE),
                     na.rm = TRUE) / qchisq(0.5, df = 1)
    add_metric("lambda_gc", lambda)
    add_validation("lambda_gc", NA_real_, lambda, "[0.95, 1.05]",
                   lambda >= 0.95 && lambda <= 1.05)
  }
}

# ---- write tables ----
if (length(summary_metrics) > 0) {
  sm <- tibble(metric = names(summary_metrics),
               value = vapply(summary_metrics, function(x)
                 ifelse(is.null(x) || length(x) == 0, NA_real_,
                        as.numeric(x)), numeric(1)))
  write_tsv(sm, file.path(outdir, "diagnostic_summary_stats.tsv"))
  cat("  Written:", file.path(outdir, "diagnostic_summary_stats.tsv"), "\n")
}
if (length(validation) > 0) {
  vv <- bind_rows(validation)
  write_tsv(vv, file.path(outdir, "validation_summary.tsv"))
  cat("  Written:", file.path(outdir, "validation_summary.tsv"), "\n")
  print(as.data.frame(vv), row.names = FALSE)
}

cat("\n29e done. Figures in:", normalizePath(outdir), "\n")
