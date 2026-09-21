.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))
library(data.table)
library(ggplot2)
library(rtracklayer)

# ============================================================
# Configuration
# ============================================================

POOLED_DIR    <- "/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/pooled/genotypes"
GTF_FILE      <- "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/combined_GTF/lr_assembly_GENCODE_v45_combined.gtf"
REPORT_DIR    <- file.path(POOLED_DIR, "report")

ANC_KEEP      <- c("AFR", "AMR", "EAS", "EUR", "SAS")
COHORTS       <- c("GUSTO", "NIEHS_RICHS", "SNUH", "MALI_G3A", "NIGMS")

RSQ_THRESH    <- 0.4
MAF_THRESH    <- 0.01
HWE_THRESH    <- 1e-6
GENO_THRESH   <- 0.05
CIS_WINDOW_BP <- 1e6
N_CIS_MIN     <- 50L

CHR_LENGTHS <- c(
  chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
  chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
  chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
  chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
  chr17=83257441,  chr18=80373285, chr19=58617616,  chr20=64444167,
  chr21=46709983,  chr22=50818468
)

dir.create(REPORT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# SECTION 1: Sample counts by cohort and ancestry
# ============================================================

cat("\n=== Section 1: Sample counts ===\n")

anc_file <- file.path(POOLED_DIR, "pooled_sample_ancestry.tsv")
anc_all  <- NULL

if (file.exists(anc_file)) {
  anc_all <- fread(anc_file)
  anc_all <- anc_all[assigned_ancestry %in% ANC_KEEP]
  anc_all[, assigned_ancestry := factor(assigned_ancestry, levels = ANC_KEEP)]
  anc_all[, cohort := factor(cohort, levels = COHORTS)]
  
  # Wide summary table
  anc_tab <- dcast(anc_all, cohort ~ assigned_ancestry,
                   value.var = "assigned_ancestry",
                   fun.aggregate = length, drop = FALSE)
  anc_tab[, Total := rowSums(.SD), .SDcols = ANC_KEEP]
  total_row        <- as.list(colSums(anc_tab[, -1, with = FALSE]))
  total_row$cohort <- "Total"
  anc_tab <- rbind(anc_tab, total_row, fill = TRUE)
  
  fwrite(anc_tab,
         file.path(REPORT_DIR, "pooled_sample_counts.tsv"),
         sep = "\t")
  cat("  Sample count table written.\n")
  
  # Plot
  p_samples <- ggplot(anc_all, aes(x = cohort, fill = assigned_ancestry)) +
    geom_bar(position = "stack") +
    labs(x = NULL, y = "N samples", fill = "Ancestry",
         title = "Sample composition by cohort and ancestry stratum") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "right")
  ggsave(file.path(REPORT_DIR, "pooled_sample_counts.png"),
         p_samples, width = 9, height = 4, dpi = 150)
  cat("  Sample count plot saved.\n")
} else {
  warning("pooled_sample_ancestry.tsv not found: ", anc_file)
}

# ============================================================
# SECTION 2: Variant counts passing pooled QC
# ============================================================

cat("\n=== Section 2: Variant counts ===\n")

summary_file <- file.path(POOLED_DIR, "pooled_variant_summary.tsv")
var_sum      <- NULL

if (file.exists(summary_file)) {
  var_sum <- fread(summary_file)
  var_sum <- var_sum[ancestry %in% ANC_KEEP]
  var_sum[, ancestry := factor(ancestry, levels = ANC_KEEP)]
  setorder(var_sum, ancestry)
  
  fwrite(var_sum,
         file.path(REPORT_DIR, "pooled_variant_summary.tsv"),
         sep = "\t")
  cat("  Variant summary table written.\n")
} else {
  warning("pooled_variant_summary.tsv not found: ", summary_file)
}

# Per-chromosome density
chr_list <- lapply(ANC_KEEP, function(anc) {
  f <- file.path(POOLED_DIR, paste0(anc, "_pooled_QC_summary_per_chr.tsv"))
  if (!file.exists(f)) return(NULL)
  d <- fread(f)
  d[, ancestry := anc]
  d
})
chr_all <- rbindlist(chr_list, use.names = TRUE, fill = TRUE)

if (nrow(chr_all) > 0) {
  setnames(chr_all, "CHROM", "chr", skip_absent = TRUE)
  chr_all[, chr := paste0("chr", sub("^chr", "", chr))]
  chr_all[, chr := factor(chr, levels = paste0("chr", 1:22))]
  chr_all[, ancestry := factor(ancestry, levels = ANC_KEEP)]
  chr_all[, chr_len_mb  := CHR_LENGTHS[as.character(chr)] / 1e6]
  chr_all[, snps_per_mb := round(n_pass / chr_len_mb, 1)]
  
  fwrite(chr_all,
         file.path(REPORT_DIR, "pooled_variants_per_chr.tsv"),
         sep = "\t")
  
  p_chr <- ggplot(chr_all, aes(x = chr, y = snps_per_mb, fill = ancestry)) +
    geom_col(position = position_dodge(width = 0.85)) +
    facet_wrap(~ ancestry, ncol = 1, scales = "free_y") +
    labs(x = "Chromosome", y = "Variants / Mb", fill = "Ancestry",
         title = "Pooled passing variants per Mb per chromosome, by ancestry stratum",
         caption = paste0("MAF >= ", MAF_THRESH,
                          " | HWE P > ", HWE_THRESH,
                          " | geno <= ", GENO_THRESH * 100, "%")) +
    theme_bw(base_size = 11) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  ggsave(file.path(REPORT_DIR, "pooled_variants_per_chr.png"),
         p_chr, width = 10, height = 3 * length(ANC_KEEP), dpi = 150)
  cat("  Per-chromosome density plot saved.\n")
} else {
  warning("No per-chromosome QC summary files found in ", POOLED_DIR)
}

# ============================================================
# SECTION 3: Cis-window variant density
# ============================================================

cat("\n=== Section 3: Cis-window variant density ===\n")
cat("Loading GTF:", GTF_FILE, "\n")

gtf <- import(GTF_FILE)
tx  <- as.data.frame(gtf[gtf$type == "transcript"])
tx  <- tx[as.character(tx$seqnames) %in% paste0("chr", 1:22), ]

tss <- data.table(
  transcript_id = tx$transcript_id,
  gene_id       = tx$gene_id,
  chr           = as.character(tx$seqnames),
  TSS           = ifelse(tx$strand == "+", tx$start, tx$end)
)
tss[, win_start := pmax(0L, TSS - CIS_WINDOW_BP)]
tss[, win_end   := TSS + CIS_WINDOW_BP]
tss[, chr       := factor(chr, levels = paste0("chr", 1:22))]
setkey(tss, chr, win_start, win_end)
cat("  Transcripts loaded:", nrow(tss), "\n")

cis_list <- lapply(ANC_KEEP, function(anc) {
  pvar_file    <- file.path(POOLED_DIR, paste0(anc, "_pooled.pvar"))
  snplist_file <- file.path(POOLED_DIR, paste0(anc, "_pooled_passing.snplist"))
  
  if (file.exists(pvar_file)) {
    pvar <- fread(cmd = paste0("grep -v '^##' ", pvar_file), sep = "\t")
    setnames(pvar, 1, "CHROM")
    pvar[, chr := factor(paste0("chr", sub("^chr", "", CHROM)),
                         levels = paste0("chr", 1:22))]
    vars <- pvar[, .(chr, POS)]
  } else if (file.exists(snplist_file)) {
    ids   <- fread(snplist_file, header = FALSE)$V1
    parts <- tstrsplit(ids, ":", fixed = TRUE)
    vars  <- data.table(
      chr = factor(paste0("chr", sub("^chr", "", parts[[1]])),
                   levels = paste0("chr", 1:22)),
      POS = as.integer(parts[[2]])
    )
  } else {
    warning("No pvar or snplist found for ancestry: ", anc)
    return(NULL)
  }
  
  vars <- vars[!is.na(chr) & !is.na(POS)]
  vars[, POS_end := POS]
  setkey(vars, chr, POS, POS_end)
  
  hits <- foverlaps(vars, tss,
                    by.x = c("chr", "POS", "POS_end"),
                    by.y = c("chr", "win_start", "win_end"),
                    type = "within", nomatch = 0L)
  
  counts <- hits[, .(n_cis_variants = .N), by = .(transcript_id, gene_id, chr)]
  all_tx <- unique(tss[, .(transcript_id, gene_id, chr)])
  counts <- merge(all_tx, counts,
                  by  = c("transcript_id", "gene_id", "chr"), all.x = TRUE)
  counts[is.na(n_cis_variants), n_cis_variants := 0L]
  counts[, ancestry := anc]
  cat("  ", anc, ": cis-window counts computed\n")
  counts
})

cis_all <- rbindlist(cis_list, use.names = TRUE, fill = TRUE)
cis_all[, ancestry     := factor(ancestry, levels = ANC_KEEP)]
cis_all[, low_coverage := n_cis_variants < N_CIS_MIN]

fwrite(cis_all,
       file.path(REPORT_DIR, "pooled_cis_window_variant_counts.tsv"),
       sep = "\t")
cat("  Full cis-window counts written.\n")

# Summary table
cis_sum <- cis_all[, .(
  n_transcripts = .N,
  median_cis    = as.numeric(median(n_cis_variants, na.rm = TRUE)),
  mean_cis      = round(mean(n_cis_variants, na.rm = TRUE), 1),
  n_low         = sum(low_coverage, na.rm = TRUE),
  pct_low       = round(mean(low_coverage, na.rm = TRUE) * 100, 1)
), by = ancestry]
setorder(cis_sum, ancestry)

fwrite(cis_sum,
       file.path(REPORT_DIR, "pooled_cis_window_summary.tsv"),
       sep = "\t")
cat("  Cis-window summary table written.\n")

# Per-chromosome median
cis_chr <- cis_all[, .(
  median_cis = median(as.numeric(n_cis_variants), na.rm = TRUE)
), by = .(ancestry, chr)]

fwrite(cis_chr,
       file.path(REPORT_DIR, "pooled_cis_window_per_chr.tsv"),
       sep = "\t")

# Plot: histogram
p_cis_hist <- ggplot(cis_all, aes(x = n_cis_variants, fill = ancestry)) +
  geom_histogram(bins = 60, colour = NA) +
  geom_vline(xintercept = N_CIS_MIN, linetype = "dashed",
             colour = "black", linewidth = 0.5) +
  facet_wrap(~ ancestry, ncol = 1, scales = "free_y") +
  scale_x_log10() +
  labs(
    title   = "Passing cis-window variants per transcript (pooled stratum)",
    x       = paste0("Variants within +/- ", CIS_WINDOW_BP / 1e6, " Mb of TSS (log10)"),
    y       = "Number of transcripts",
    caption = paste0("Dashed line: N = ", N_CIS_MIN, " variant minimum")
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")
ggsave(file.path(REPORT_DIR, "pooled_cis_window_hist.png"),
       p_cis_hist, width = 8, height = 3 * length(ANC_KEEP), dpi = 150)

# Plot: per-chromosome median
p_cis_chr <- ggplot(cis_chr, aes(x = chr, y = median_cis, fill = ancestry)) +
  geom_col(position = position_dodge(width = 0.85)) +
  facet_wrap(~ ancestry, ncol = 1, scales = "free_y") +
  labs(
    title   = "Median cis-window variant count per chromosome (pooled stratum)",
    x       = "Chromosome",
    y       = "Median variants per cis-window",
    caption = paste0("MAF >= ", MAF_THRESH,
                     " | HWE P > ", HWE_THRESH,
                     " | geno <= ", GENO_THRESH * 100, "%")
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave(file.path(REPORT_DIR, "pooled_cis_window_per_chr.png"),
       p_cis_chr, width = 10, height = 3 * length(ANC_KEEP), dpi = 150)

# Plot: % power-limited transcripts
p_cis_pct <- ggplot(cis_sum, aes(x = ancestry, y = pct_low, fill = ancestry)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(pct_low, "%")), vjust = -0.4, size = 3.5) +
  labs(
    x     = NULL,
    y     = paste0("% transcripts < ", N_CIS_MIN, " cis-variants"),
    title = "Proportion of transcripts with limited xQTL mapping power, by ancestry stratum"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
ggsave(file.path(REPORT_DIR, "pooled_cis_pct_low.png"),
       p_cis_pct, width = 6, height = 4, dpi = 150)

cat("\n=== pooled_imputation_QC_generate.R complete ===\n")
cat("All outputs written to:", REPORT_DIR, "\n")