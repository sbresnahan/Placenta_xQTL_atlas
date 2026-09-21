# =============================================================================
# xqtl_report_functions.R — analysis functions for the placenta multi-ancestry
# xQTL report (PANTRY-style MAF 0.05 results: grouped, ungrouped, combined
# layers; honest-null stance; NO ACAT — dropped per user decision).
#
# Inputs (per ancestry ANC, modality MOD), all complete sorted TSVs:
#   grouped:   {ANC}_{MOD}_cisqtl_top.tsv        (one row per gene; carries
#              group_id, group_size; phenotype_id = driver phenotype)
#   ungrouped: {ANC}_{MOD}_ungrouped_cisqtl_top.tsv (one row per phenotype)
#   combined:  {ANC}_combined_cisqtl_top.tsv     (cross-modality gene groups;
#              phenotype_id is modality-namespaced: {mod}__{original})
#   expression:{ANC}_expression_cisqtl_top.tsv   (gene-level; no group_id)
#   independent: {ANC}_{MOD}_cisqtl_independent_top.tsv (may be EMPTY — 1 byte)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

MODALITY_LEVELS <- c("expression", "isoforms", "splicing", "intron_retention",
                     "alt_TSS", "alt_polyA", "RNA_editing", "stability")

# Colorblind-friendly modality palette
MODALITY_COLORS <- c(
  expression       = "#000000",
  isoforms         = "#0279EE",
  splicing         = "#FF9400",
  intron_retention = "#75A025",
  alt_TSS          = "#FD9BED",
  alt_polyA        = "#E9ED4C",
  RNA_editing      = "#7B3294",
  stability        = "#8C8C8C"
)

theme_report <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(text = element_text(family = "Liberation Sans"),
          strip.background = element_rect(fill = "grey92", color = NA),
          strip.text = element_text(face = "bold"),
          legend.key.size = unit(0.45, "cm"))
}

strip_version <- function(x) sub("\\..*$", "", x)

# ---- ID parsing -------------------------------------------------------------
# Handles bare IDs (per-modality layers) and namespaced IDs (combined layer:
# {modality}__{original_id}). Gene = first '__' token of the ORIGINAL id,
# version stripped. Returns df with modality_from_id, gene_id, original_id.
parse_phenotype_ids <- function(phenotype_id, layer_modality) {
  tok1 <- str_split_fixed(phenotype_id, "__", 2)[, 1]
  is_namespaced <- tok1 %in% MODALITY_LEVELS
  tibble(
    phenotype_id = phenotype_id,
    modality_from_id = if_else(is_namespaced, tok1, layer_modality),
    original_id = if_else(is_namespaced,
                          str_split_fixed(phenotype_id, "__", 2)[, 2],
                          phenotype_id),
    gene_id = strip_version(str_split_fixed(
      if_else(is_namespaced, str_split_fixed(phenotype_id, "__", 2)[, 2],
              phenotype_id), "__", 2)[, 1])
  )
}

# ---- Low-level TSV reader (empty-aware) ------------------------------------
read_result_tsv <- function(path) {
  if (!file.exists(path)) stop("file not found: ", path)
  if (file.size(path) < 10) {
    # empty independent output (no significant groups) -> 0-row tibble
    return(tibble())
  }
  read_tsv(path, col_types = cols(.default = col_double(),
                                  phenotype_id = col_character(),
                                  variant_id = col_character(),
                                  group_id = col_character()),
           progress = FALSE)
}

# ---- Load one layer ---------------------------------------------------------
# layer: "grouped" | "ungrouped" | "combined" | "expression"
load_layer <- function(data_dir, layer,
                       ancestries = c("EAS", "EUR"),
                       modalities = setdiff(MODALITY_LEVELS, "expression")) {
  out <- list()
  for (anc in ancestries) {
    mods <- switch(layer,
                   grouped = c("expression", modalities),
                   ungrouped = modalities,
                   combined = "combined",
                   expression = "expression")
    for (mod in mods) {
      path <- switch(layer,
                     grouped = file.path(data_dir, sprintf("%s_%s_cisqtl_top.tsv", anc, mod)),
                     ungrouped = file.path(data_dir, sprintf("%s_%s_ungrouped_cisqtl_top.tsv", anc, mod)),
                     combined = file.path(data_dir, sprintf("%s_combined_cisqtl_top.tsv", anc)),
                     expression = file.path(data_dir, sprintf("%s_expression_cisqtl_top.tsv", anc)))
      if (!file.exists(path)) {
        message(sprintf("  NOTE: missing %s/%s [%s] — skipped", anc, mod, layer))
        next
      }
      df <- read_result_tsv(path)
      if (nrow(df) == 0) next
      ids <- parse_phenotype_ids(df$phenotype_id, mod)
      df$modality_from_id <- ids$modality_from_id
      df$original_id <- ids$original_id
      # gene identity: prefer group_id column (grouped layers), else ID parse
      if ("group_id" %in% names(df)) {
        df$gene_id <- strip_version(df$group_id)
      } else {
        df$gene_id <- ids$gene_id
      }
      # variant position + phenotype coordinates from distance columns
      df$variant_pos <- suppressWarnings(
        as.integer(str_split_fixed(df$variant_id, ":", 4)[, 2]))
      df$pheno_start <- df$variant_pos - df$start_distance
      df$pheno_end   <- df$variant_pos - df$end_distance
      df$ancestry <- anc
      df$modality <- mod
      df$layer <- layer
      out[[paste(anc, mod, layer)]] <- df
    }
  }
  res <- bind_rows(out)
  if (nrow(res) == 0) return(res)
  res %>% mutate(modality = factor(modality, levels = c(MODALITY_LEVELS, "combined")),
                 ancestry = factor(ancestry, levels = ancestries))
}

# Load independent (stepwise) outputs; empty-aware. Returns 0-row tibble if
# all are empty (expected for this result set).
load_independent <- function(data_dir,
                             ancestries = c("EAS", "EUR"),
                             modalities = c("expression", setdiff(MODALITY_LEVELS, "expression"), "combined")) {
  out <- list()
  for (anc in ancestries) for (mod in modalities) {
    path <- file.path(data_dir, sprintf("%s_%s_cisqtl_independent_top.tsv", anc, mod))
    if (!file.exists(path)) next
    df <- read_result_tsv(path)
    if (nrow(df) == 0) next
    df$ancestry <- anc; df$modality <- mod
    out[[paste(anc, mod)]] <- df
  }
  bind_rows(out)
}

# ---- Diagnostics ------------------------------------------------------------
# Storey pi1 (fraction of true alternatives), lambda = 0.5
storey_pi1 <- function(p, lambda = 0.5) {
  p <- p[is.finite(p)]
  if (length(p) == 0) return(NA_real_)
  pi0 <- min(1, mean(p >= lambda) / (1 - lambda))
  max(0, 1 - pi0)
}

# QQ-plot data: expected vs observed -log10 p
qq_data <- function(p) {
  p <- p[is.finite(p) & p > 0]
  n <- length(p)
  tibble(expected = -log10(ppoints(n)),
         observed = -log10(sort(p)))
}

# ---- Top hits ---------------------------------------------------------------
top_hits <- function(df, n = 20) {
  df %>% arrange(qval, pval_beta) %>% head(n)
}

# ---- Driver modality (combined layer) ---------------------------------------
# The grouped-output row label is the phenotype with the best nominal
# association in the group (tensorQTL grouped map_cis) — i.e. the driver.
driver_modality <- function(combined_df) {
  combined_df %>%
    mutate(driver_modality = factor(modality_from_id, levels = MODALITY_LEVELS))
}

# ---- Gene-body-relative variant position (PANTRY Fig 2d style) --------------
# Requires gene bodies (Ensembl: chrom/start/end/strand). Strand-aware:
# 0 = gene start (5' end), 1 = gene end (3' end); <0 upstream, >1 downstream.
# NOTE: phenotype BED coordinates are NOT usable for this (expression uses TSS
# points, splicing uses junction spans) — gene bodies must come from annotation.
add_gene_relative_position <- function(df, gene_bodies) {
  gb <- gene_bodies %>% select(gene_id, gene_start = start, gene_end = end,
                               strand)
  df %>%
    left_join(gb, by = "gene_id") %>%
    mutate(
      gene_len = pmax(gene_end - gene_start, 1L),
      rel_pos = dplyr::case_when(
        is.na(gene_start) ~ NA_real_,
        strand == -1L ~ (gene_end - variant_pos) / gene_len,
        TRUE ~ (variant_pos - gene_start) / gene_len
      )
    )
}

# ---- Cross-ancestry concordance (comparative only; no pooled inference) -----
# Slope concordance: join the two ancestries' lead associations per gene
# (same layer/modality), report Spearman correlation of slopes where the SAME
# variant is lead in both, plus rank concordance of -log10(pval_beta) over
# genes tested in both.
cross_ancestry_concordance <- function(df_a, df_b, min_overlap = 20) {
  ja <- df_a %>% select(gene_id, variant_id, slope, pval_beta, qval)
  jb <- df_b %>% select(gene_id, variant_id, slope, pval_beta, qval)
  j <- ja %>% inner_join(jb, by = "gene_id", suffix = c("_a", "_b"))
  if (nrow(j) < min_overlap) {
    return(tibble(n_genes = nrow(j), same_variant = NA_integer_,
                  slope_spearman = NA_real_, rank_spearman = NA_real_,
                  top100_overlap = NA_integer_))
  }
  same <- j %>% filter(variant_id_a == variant_id_b)
  top_a <- j %>% slice_min(pval_beta_a, n = 100, with_ties = FALSE)
  top_b <- j %>% slice_min(pval_beta_b, n = 100, with_ties = FALSE)
  tibble(
    n_genes = nrow(j),
    same_variant = nrow(same),
    slope_spearman = if (nrow(same) >= min_overlap)
      cor(same$slope_a, same$slope_b, method = "spearman") else NA_real_,
    rank_spearman = cor(-log10(j$pval_beta_a), -log10(j$pval_beta_b),
                        method = "spearman"),
    top100_overlap = length(intersect(top_a$gene_id, top_b$gene_id))
  )
}

# ---- Power analysis ---------------------------------------------------------
# Minimum detectable r2 (variance explained) for a single-variant linear
# association test at sample size n, significance threshold alpha, given
# power. Noncentrality lambda = n * r2 / (1 - r2).
min_detectable_r2 <- function(n, alpha, power = 0.8) {
  fcrit <- qf(1 - alpha, df1 = 1, df2 = n - 2)
  f <- function(r2) {
    lambda <- n * r2 / (1 - r2)
    pf(fcrit, df1 = 1, df2 = n - 2, ncp = lambda) - (1 - power)
  }
  uniroot(f, interval = c(1e-6, 0.95))$root
}

# Power to detect a variant explaining r2 at sample size n, threshold alpha
power_at_n <- function(r2, n, alpha) {
  lambda <- n * r2 / (1 - r2)
  fcrit <- qf(1 - alpha, df1 = 1, df2 = n - 2)
  1 - pf(fcrit, df1 = 1, df2 = n - 2, ncp = lambda)
}

# Approximate r2 of an observed association on an INT-normalized phenotype:
# r2 ~ 2 * af * (1 - af) * slope^2  (phenotype variance ~ 1 after INT)
approx_r2 <- function(slope, af) 2 * af * (1 - af) * slope^2
