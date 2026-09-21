.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))
library(GENESIS)
library(SNPRelate)
library(GWASTools)
library(data.table)
library(ggplot2)
library(dplyr)

# ============================================================
# Configuration
# ============================================================

COHORT   <- "NIGMS"

OUTDIR     <- paste0("/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/",COHORT,"/genotypes/imputed/qc")
REPORT_DIR <- file.path(OUTDIR, "report"); dir.create(REPORT_DIR, showWarnings = FALSE)
KG_DIR   <- "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/1kGP"
KG_LABELS <- file.path(KG_DIR, "1kg_superpop.tsv")

N_PCS    <- 20  # number of PCs to compute
PCA_PCS  <- 10  # PCs to use for ancestry projection

# Merged cohort + 1KG pgen from shell prep script
MERGED_PGEN <- file.path(OUTDIR, paste0(COHORT, "_plus_kg.for_pcair"))

# ============================================================
# SECTION 1: Convert BED to GDS for GENESIS
# ============================================================

cat("\n=== Step 1: Converting BED to GDS ===\n")

gds_file <- file.path(OUTDIR, paste0(COHORT, "_plus_kg.gds"))

SNPRelate::snpgdsBED2GDS(
  bed.fn    = paste0(MERGED_PGEN, ".bed"),
  fam.fn    = paste0(MERGED_PGEN, ".fam"),
  bim.fn    = paste0(MERGED_PGEN, ".bim"),
  out.gdsfn = gds_file
)
cat("GDS written to:", gds_file, "\n")

gds <- snpgdsOpen(gds_file)

# ============================================================
# SECTION 2: KING kinship estimation
# ============================================================

cat("\n=== Step 2: KING kinship estimation ===\n")

king_file <- file.path(OUTDIR, paste0(COHORT, "_king"))

king <- snpgdsIBDKING(gds, num.thread = 12)
king_mat <- king$kinship
rownames(king_mat) <- king$sample.id
colnames(king_mat) <- king$sample.id

saveRDS(king_mat, paste0(king_file, ".rds"))
cat("KING kinship matrix saved.\n")

# ============================================================
# SECTION 3: PC-AiR
# ============================================================

cat("\n=== Step 3: Running PC-AiR ===\n")

# Load 1KG sample IDs to define reference panel
kg_labels <- fread(KG_LABELS)
kg_sample_ids <- kg_labels[[1]]  # assumes first column is sample ID

all_sample_ids <- read.gdsn(index.gdsn(gds, "sample.id"))
cohort_ids     <- setdiff(all_sample_ids, kg_sample_ids)

cat("Cohort samples:", length(cohort_ids), "\n")
cat("1KG reference samples:", sum(all_sample_ids %in% kg_sample_ids), "\n")

pcair_result <- pcair(
  gdsobj     = gds,
  kinobj     = king_mat,
  divobj     = king_mat,
  sample.include = all_sample_ids,
  unrel.set  = kg_sample_ids,  # use 1KG as the unrelated reference
  num.cores = 12
)

cat("PC-AiR complete.\n")
saveRDS(pcair_result, file.path(OUTDIR, paste0(COHORT, "_pcair.rds")))

# Extract PCs for cohort samples only
pcs <- as.data.frame(pcair_result$vectors[cohort_ids, 1:N_PCS])
colnames(pcs) <- paste0("PC", 1:N_PCS)
pcs$sample_id <- rownames(pcs)

# ============================================================
# SECTION 4: Ancestry assignment via 1KG reference projection
# ============================================================

cat("\n=== Step 4: Ancestry assignment ===\n")

# Get 1KG PCs and superpop labels
kg_pcs <- as.data.frame(pcair_result$vectors[kg_sample_ids, 1:PCA_PCS])
colnames(kg_pcs) <- paste0("PC", 1:PCA_PCS)
kg_pcs$sample_id <- rownames(kg_pcs)

kg_labels_merged <- merge(kg_pcs, kg_labels, by.x = "sample_id", by.y = colnames(kg_labels)[1])
superpop_col <- colnames(kg_labels)[2]  # assumes second column is superpop

# Compute per-superpop centroid in PC space
pc_cols <- paste0("PC", 1:PCA_PCS)
centroids <- kg_labels_merged %>%
  dplyr::group_by(.data[[superpop_col]]) %>%
  dplyr::summarise(across(all_of(pc_cols), mean), .groups = "drop")

# Assign each cohort sample to nearest centroid (Euclidean distance)
cohort_pcs <- as.data.frame(pcair_result$vectors[cohort_ids, 1:PCA_PCS])
colnames(cohort_pcs) <- pc_cols

assign_ancestry <- function(sample_pcs, centroids, pc_cols, pop_col) {
  dists <- apply(centroids[, pc_cols], 1, function(centroid) {
    sqrt(sum((sample_pcs - centroid)^2))
  })
  centroids[[pop_col]][which.min(dists)]
}

ancestry_assignments <- data.frame(
  sample_id         = cohort_ids,
  assigned_ancestry = sapply(cohort_ids, function(sid) {
    assign_ancestry(cohort_pcs[sid, ], centroids, pc_cols, superpop_col)
  }),
  stringsAsFactors = FALSE
)

cat("\nAncestry assignments:\n")
print(table(ancestry_assignments$assigned_ancestry))

# Write ancestry assignments (input for scripts 3-4 and Rmd report)
ancestry_out <- file.path(OUTDIR, paste0(COHORT, "_ancestry_assignments.tsv"))
fwrite(ancestry_assignments, ancestry_out, sep = "\t")

ancestry_out <- file.path(REPORT_DIR, paste0(COHORT, "_ancestry_assignments.tsv"))
fwrite(ancestry_assignments, ancestry_out, sep = "\t")

cat("Ancestry assignments written to:", ancestry_out, "\n")

# ============================================================
# SECTION 5: Write cohort PCs for use as eQTL covariates
# ============================================================

cat("\n=== Step 5: Writing PCs for covariate use ===\n")

pcs_out <- merge(pcs, ancestry_assignments, by = "sample_id")
fwrite(pcs_out, file.path(OUTDIR, paste0(COHORT, "_pcair_pcs.tsv")), sep = "\t")
cat("PCs written to:", file.path(OUTDIR, paste0(COHORT, "_pcair_pcs.tsv")), "\n")

# ============================================================
# SECTION 6: PCA plots
# ============================================================

cat("\n=== Step 6: Plotting ===\n")

# Combine cohort + 1KG PCs for plotting
kg_plot <- kg_labels_merged[, c("sample_id", pc_cols, superpop_col)]
colnames(kg_plot)[colnames(kg_plot) == superpop_col] <- "ancestry"
kg_plot$group <- "1KG reference"

cohort_plot <- merge(pcs[, c("sample_id", pc_cols)], ancestry_assignments,
                     by = "sample_id")
colnames(cohort_plot)[colnames(cohort_plot) == "assigned_ancestry"] <- "ancestry"
cohort_plot$group <- COHORT

plot_df <- rbind(kg_plot[, c("sample_id", "PC1", "PC2", "PC3", "PC4", "ancestry", "group")],
                 cohort_plot[, c("sample_id", "PC1", "PC2", "PC3", "PC4", "ancestry", "group")])

library(patchwork)

make_pcair_plot <- function(df, pc_x, pc_y) {
  anc_levels <- sort(unique(df$ancestry))
  anc_colors <- setNames(scales::hue_pal()(length(anc_levels)), anc_levels)
  df$ancestry <- factor(df$ancestry, levels = anc_levels)
  
  kg_panel <- ggplot(df[df$group == "1KG reference", ],
                     aes(x = .data[[pc_x]], y = .data[[pc_y]], color = ancestry)) +
    geom_point(size = 1.2) +
    scale_color_manual(values = anc_colors, drop = FALSE) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    labs(title = "1KG reference") +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")
  
  cohort_panel <- ggplot(df[df$group == COHORT, ],
                         aes(x = .data[[pc_x]], y = .data[[pc_y]], color = ancestry)) +
    geom_point(size = 1.2) +
    scale_color_manual(values = anc_colors, drop = FALSE) +
    labs(title = COHORT) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
  
  kg_panel + cohort_panel + plot_layout(guides = "collect")
}

p1 <- make_pcair_plot(plot_df, "PC1", "PC2") +
  plot_annotation(title = paste(COHORT, "— PC-AiR: PC1 vs PC2"),
                  subtitle = paste("Cohort n =", length(cohort_ids),
                                   "| 1KG reference n =", sum(all_sample_ids %in% kg_sample_ids)))
ggsave(file.path(REPORT_DIR, paste0(COHORT, "_pcair_PC1_PC2.png")),
       p1, width = 12, height = 5, dpi = 150)

p2 <- make_pcair_plot(plot_df, "PC3", "PC4") +
  plot_annotation(title = paste(COHORT, "— PC-AiR: PC3 vs PC4"))
ggsave(file.path(REPORT_DIR, paste0(COHORT, "_pcair_PC3_PC4.png")),
       p2, width = 12, height = 5, dpi = 150)

snpgdsClose(gds)
cat("\n=== Script 2 complete ===\n")