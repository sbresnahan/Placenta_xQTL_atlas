# =============================================================================
# plot_genotype_pcs_scree.R — Scree plots for PC-AiR genotype PCs
# =============================================================================
# Run interactively in RStudio on seadragon. Reads the PC TSV files produced
# by 22_genotype_pcair.sh, computes per-PC variance (proportional to
# eigenvalues), and draws scree plots with cumulative variance.
#
# Usage:
#   .libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))
#   source("/rsrch5/.../scripts/plot_genotype_pcs_scree.R")
#
# Or set PCAIR_DIR and ANCESTRIES at the top and source.
# =============================================================================

.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))

library(ggplot2)

# ---- Parameters (edit as needed) ----
PCAIR_DIR <- "/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/genotype_pcs"
ANCESTRIES <- c("EAS", "EUR")
N_PCS <- 20  # number of PCs in the TSV (must match what was computed)

# ---- Helper: compute variance-based scree from PC scores ----
# PC-AiR returns PC scores (vectors). The variance of each PC column across
# samples is proportional to the eigenvalue — larger variance = more
# population structure captured. This is the standard proxy when the
# singular values themselves aren't accessible.
compute_scree <- function(pcs, n_pcs) {
  pc_cols <- paste0("PC", 1:n_pcs)
  pc_cols <- intersect(pc_cols, colnames(pcs))
  variances <- sapply(pc_cols, function(c) var(pcs[[c]], na.rm = TRUE))
  pct_var <- variances / sum(variances) * 100
  cum_pct <- cumsum(pct_var)
  data.frame(
    PC = seq_along(variances),
    variance = variances,
    pct_var = pct_var,
    cum_pct = cum_pct
  )
}

# ---- Helper: draw scree plot ----
make_scree <- function(df, ancestry) {
  n <- nrow(df)
  ggplot(df, aes(x = PC, y = pct_var)) +
    geom_col(fill = "#0279EE", alpha = 0.7) +
    geom_line(aes(y = cum_pct), color = "#FF9400", linewidth = 0.8) +
    geom_point(aes(y = cum_pct), color = "#FF9400", size = 1.5) +
    scale_y_continuous(
      name = "Variance explained (%)",
      sec.axis = sec_axis(trans = ~., name = "Cumulative (%)")
    ) +
    scale_x_continuous(breaks = 1:n) +
    labs(
      title = sprintf("Genotype PC scree plot — %s", ancestry),
      subtitle = sprintf("%d PCs | variance computed from PC-AiR scores", n)
    ) +
    theme_bw(base_size = 11) +
    theme(
      text = element_text(family = "Liberation Sans"),
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# ---- Main loop ----
for (anc in ANCESTRIES) {
  pcs_path <- file.path(PCAIR_DIR, sprintf("%s_genotype_pcs.tsv", anc))
  if (!file.exists(pcs_path)) {
    cat(sprintf("  SKIP: %s — file not found: %s\n", anc, pcs_path))
    next
  }
  
  pcs <- read.csv(pcs_path, sep = "\t", check.names = FALSE)
  cat(sprintf("\n[%s] %d samples, %d PCs\n", anc, nrow(pcs),
              ncol(pcs) - 1))  # minus sample_id column
  
  scree_df <- compute_scree(pcs, N_PCS)
  
  # Print summary table
  cat("\n  PC   Variance    %% Var   Cum %%\n")
  cat("  ---  ----------  -------  ------\n")
  for (i in 1:min(nrow(scree_df), 10)) {
    cat(sprintf("  PC%-2d %.6f   %5.2f%%   %5.2f%%\n",
                scree_df$PC[i], scree_df$variance[i],
                scree_df$pct_var[i], scree_df$cum_pct[i]))
  }
  if (nrow(scree_df) > 10) {
    cat(sprintf("  ...  (%d more PCs)\n", nrow(scree_df) - 10))
  }
  
  # Suggest elbow
  # Simple heuristic: first PC where pct_var drops below 1/3 of PC1's value
  elbow <- which(scree_df$pct_var < scree_df$pct_var[1] / 3)[1]
  if (!is.na(elbow)) {
    cat(sprintf("\n  Suggested elbow (heuristic): PC%d (%.2f%% variance)\n",
                elbow, scree_df$pct_var[elbow]))
    cat(sprintf("  Cumulative variance at PC%d: %.2f%%\n",
                elbow, scree_df$cum_pct[elbow]))
  }
  
  # Draw and save
  p <- make_scree(scree_df, anc)
  print(p)  # display in RStudio Plots pane
  
  out_path <- file.path(PCAIR_DIR, sprintf("%s_genotype_pcs_scree.png", anc))
  ggsave(out_path, p, width = 8, height = 5, dpi = 150)
  cat(sprintf("  Saved: %s\n", out_path))
}

cat("\nDone. Adjust N_PCS at the top if fewer PCs were computed.\n")
