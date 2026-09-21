## =============================================================================
## explore_deconvolution_outputs.R
##
## Interactive exploration of the pooled MuSiC deconvolution outputs:
##   1. all_cohorts_cell_proportions_full.tsv      (samples x 27 cell types)
##   2. all_cohorts_cell_proportions_collapsed.tsv (samples x 8 collapsed types)
##   3. all_cohorts_maternal_flag.tsv              (sample, cohort, maternal_fraction, flagged)
##
## Produces (saved as PNGs in <deconv_dir>/plots/ and printed interactively):
##   - Stacked bar charts per cohort (x = sample, y = proportion, fill = cell type),
##     one set for the full 27-type table and one for the collapsed 8-type table
##   - Boxplot + jittered points of maternal_fraction by cohort, with the 10%
##     flag threshold marked
##
## Run interactively in RStudio. Requires: ggplot2, tidyr, dplyr.
##   install.packages(c("ggplot2", "tidyr", "dplyr"))
## =============================================================================

.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))

library(ggplot2)
library(tidyr)
library(dplyr)

## ---- Config -----------------------------------------------------------------
deconv_dir <- "/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY/deconvolution"
out_dir    <- file.path(deconv_dir, "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

MATERNAL_THRESHOLD <- 0.10  # flag threshold used by run_music_deconvolution.R

## Cohort display names — edit this list to rename cohorts in plots, tables,
## and output filenames. Cohorts not listed here keep their original IDs.
cohort_labels <- c(cohort1 = "NIEHS_RICHS",
                   cohort2 = "GUSTO",
                   cohort3 = "SNUH",
                   cohort4 = "NIGMS")

## ---- Read inputs -------------------------------------------------------------
prop_full <- read.delim(file.path(deconv_dir, "all_cohorts_cell_proportions_full.tsv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
prop_coll <- read.delim(file.path(deconv_dir, "all_cohorts_cell_proportions_collapsed.tsv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
mat_flag  <- read.delim(file.path(deconv_dir, "all_cohorts_maternal_flag.tsv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
mat_flag$flagged <- mat_flag$flagged == "True"


## ---- Apply cohort display names ----------------------------------------------
rename_cohorts <- function(x) {
  unknown <- setdiff(unique(x), names(cohort_labels))
  if (length(unknown) > 0)
    message("Note: no display name for cohort(s): ", paste(unknown, collapse = ", "),
            " — keeping original IDs")
  ifelse(x %in% names(cohort_labels), cohort_labels[x], x)
}
prop_full$cohort <- rename_cohorts(prop_full$cohort)
prop_coll$cohort <- rename_cohorts(prop_coll$cohort)
mat_flag$cohort  <- rename_cohorts(mat_flag$cohort)

## Quick sanity summary
cat("Full proportions:     ", nrow(prop_full), "samples x",
    ncol(prop_full) - 2, "cell types\n")
cat("Collapsed proportions:", nrow(prop_coll), "samples x",
    ncol(prop_coll) - 2, "cell types\n")
cat("Maternal flags:       ", nrow(mat_flag), "samples,",
    sum(mat_flag$flagged), "flagged (>", MATERNAL_THRESHOLD * 100, "% maternal)\n\n")
print(table(prop_full$cohort))
cat("\nFlagged by cohort:\n")
print(mat_flag %>% group_by(cohort) %>% summarise(n = n(), n_flagged = sum(flagged),
                                                  median_maternal = round(median(maternal_fraction), 4)))

## ---- Color palettes -----------------------------------------------------------
## 8-type palette: Okabe-Ito (colorblind-friendly)
pal8 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
          "#0072B2", "#D55E00", "#CC79A7", "#999999")

## 27-type palette: extended qualitative set (no palette is truly
## colorblind-safe at 27 categories; hues chosen for maximal separation)
pal27 <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
           "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
           "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
           "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
           "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
           "#3182bd", "#31a354")

## Save every plot as both PNG (150 dpi) and PDF (vector)
save_both <- function(p, path_stub, width, height) {
  ggsave(paste0(path_stub, ".png"), p, width = width, height = height, dpi = 150)
  ggsave(paste0(path_stub, ".pdf"), p, width = width, height = height)
}

## ---- Stacked bar chart function ------------------------------------------------
## df: pooled proportion table (sample, cohort, <cell type columns>)
## One plot per cohort; samples ordered by the dominant (first) cell type so
## trophoblast-rich vs. contaminated samples separate visually.
plot_stacked_bars <- function(df, palette, set_label, file_stub) {
  
  cell_types <- setdiff(colnames(df), c("sample", "cohort"))
  
  for (co in sort(unique(df$cohort))) {
    
    sub <- df %>% filter(cohort == co)
    
    ## Order samples by proportion of the first (dominant) cell type
    ord <- order(-sub[[cell_types[1]]])
    sub <- sub[ord, ]
    sub$sample <- factor(sub$sample, levels = sub$sample)  # lock order
    
    long <- sub %>%
      pivot_longer(all_of(cell_types), names_to = "cell_type", values_to = "proportion") %>%
      mutate(cell_type = factor(cell_type, levels = cell_types))  # keep column order
    
    p <- ggplot(long, aes(x = sample, y = proportion, fill = cell_type)) +
      geom_col(width = 1) +
      scale_fill_manual(values = palette, drop = FALSE) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(title = paste(co, "-", set_label),
           x = "sample", y = "proportion", fill = "cell type") +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_blank(),   # too many samples to label
            axis.ticks.x = element_blank(),
            panel.grid = element_blank(),
            legend.key.size = unit(0.35, "cm"))
    
    print(p)  # show in RStudio Plots pane
    safe_co <- gsub("[^A-Za-z0-9_.-]", "_", co)  # filesystem-safe label
    save_both(p, file.path(out_dir, paste0("stacked_", file_stub, "_", safe_co)),
              width = 10, height = 6)
  }
}

## ---- 1. Full 27-type stacked bars ---------------------------------------------
plot_stacked_bars(prop_full, pal27, "full 27 cell types", "full")

## ---- 2. Collapsed 8-type stacked bars ------------------------------------------
plot_stacked_bars(prop_coll, pal8, "collapsed 8 cell types", "collapsed")

## ---- 3. Maternal fraction boxplot ----------------------------------------------
p_mat <- ggplot(mat_flag, aes(x = cohort, y = maternal_fraction)) +
  geom_boxplot(outlier.shape = NA, fill = "grey85", color = "grey30") +
  geom_jitter(aes(color = flagged), width = 0.2, size = 1.2, alpha = 0.6) +
  geom_hline(yintercept = MATERNAL_THRESHOLD, linetype = "dashed",
             color = "red", linewidth = 0.6) +
  scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "red"),
                     name = paste0(">", MATERNAL_THRESHOLD * 100, "% maternal")) +
  labs(title = "Maternal contamination by cohort",
       x = "cohort", y = "maternal fraction") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(p_mat)
save_both(p_mat, file.path(out_dir, "maternal_fraction_by_cohort"),
          width = 7, height = 5)

cat("\nPlots saved to:", out_dir, "\n")
