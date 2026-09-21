.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))
library(data.table)
library(ggplot2)
library(rtracklayer)

# ============================================================
# Configuration
# ============================================================

GTF_FILE <- "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/Placenta_LRRNAseq/Placenta_LRRNAseq/STB/combined_GTF/lr_assembly_GENCODE_v45_combined.gtf"

COHORT_PVARS <- c(
  MALI    = "/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/MALI/genotypes/imputed/qc/MALI.rsq_pass.pvar",
  NIEHS_RICHS = "/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/NIEHS_RICHS/genotypes/imputed/qc/NIEHS_RICHS.rsq_pass.pvar",
  NIGMS       = "/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/NIGMS/genotypes/imputed/qc/NIGMS.rsq_pass.pvar"
)

COHORT_PVARS <- c(
  MALI    = "/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/MALI/genotypes/imputed/qc/MALI.rsq_pass.pvar"
) # for checking just MALI

COHORTS <- names(COHORT_PVARS)

CIS_WINDOW_BP <- 1e6
N_CIS_MIN     <- 50L

# ============================================================
# Build transcript TSS windows from GTF
# ============================================================

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
tss[, win_start := as.integer(pmax(0, TSS - CIS_WINDOW_BP))]
tss[, win_end   := as.integer(TSS + CIS_WINDOW_BP)]
tss[, chr       := factor(chr, levels = paste0("chr", 1:22))]
setkey(tss, chr, win_start, win_end)
cat("  Transcripts loaded:", nrow(tss), "\n")

# ============================================================
# Cis-window variant counts per cohort
# ============================================================

cis_list <- lapply(COHORTS, function(coh) {
  cat("Processing cohort:", coh, "\n"); flush.console()
  pvar_file <- COHORT_PVARS[[coh]]
  
  if (!file.exists(pvar_file)) {
    warning("Pvar file not found for cohort: ", coh, " (", pvar_file, ")")
    return(NULL)
  }
  
  pvar <- fread(cmd = paste0("grep -v '^##' ", pvar_file), sep = "\t")
  cat("  pvar rows read:", nrow(pvar), "| columns:", paste(names(pvar), collapse = ", "), "\n")
  setnames(pvar, 1, "CHROM")
  pvar[, chr := factor(paste0("chr", sub("^chr", "", CHROM)),
                       levels = paste0("chr", 1:22))]
  vars <- pvar[, .(chr, POS)]
  vars <- vars[!is.na(chr) & !is.na(POS)]
  cat("  vars rows after chr/POS NA filter:", nrow(vars), "\n")
  if (nrow(vars) == 0L) {
    warning("No valid variants after filtering for cohort: ", coh,
            " -- check CHROM value formatting (unexpected chr labels or non-numeric POS)")
    return(NULL)
  }
  vars[, POS_end := POS]
  setkey(vars, chr, POS, POS_end)
  
  chr_levels <- levels(vars$chr)
  counts_by_chr <- lapply(chr_levels, function(cc) {
    v_sub <- vars[chr == cc]
    t_sub <- tss[chr == cc]
    if (nrow(v_sub) == 0L || nrow(t_sub) == 0L) return(NULL)
    setkey(v_sub, chr, POS, POS_end)
    setkey(t_sub, chr, win_start, win_end)
    h <- foverlaps(v_sub, t_sub,
                   by.x = c("chr", "POS", "POS_end"),
                   by.y = c("chr", "win_start", "win_end"),
                   type = "within", nomatch = 0L)
    h[, .(n_cis_variants = .N), by = .(transcript_id, gene_id, chr)]
  })
  counts <- rbindlist(counts_by_chr, use.names = TRUE, fill = TRUE)
  
  all_tx <- unique(tss[, .(transcript_id, gene_id, chr)])
  counts <- merge(all_tx, counts,
                  by  = c("transcript_id", "gene_id", "chr"), all.x = TRUE)
  counts[is.na(n_cis_variants), n_cis_variants := 0L]
  counts[, cohort := coh]
  cat("  ", coh, ": cis-window counts computed\n")
  counts
})

cis_all <- rbindlist(cis_list, use.names = TRUE, fill = TRUE)
cis_all[, cohort := factor(cohort, levels = COHORTS)]

# ============================================================
# Plot: histogram
# ============================================================

p_cis_hist <- ggplot(cis_all, aes(x = n_cis_variants, fill = cohort)) +
  geom_histogram(bins = 60, colour = NA) +
  geom_vline(xintercept = N_CIS_MIN, linetype = "dashed",
             colour = "black", linewidth = 0.5) +
  facet_wrap(~ cohort, ncol = 1, scales = "free_y") +
  scale_x_log10() +
  labs(
    title   = "Passing cis-window variants per transcript (AFR strata)",
    x       = paste0("Variants within +/- ", CIS_WINDOW_BP / 1e6, " Mb of TSS (log10)"),
    y       = "Number of transcripts",
    caption = paste0("Dashed line: N = ", N_CIS_MIN, " variant minimum")
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")
print(p_cis_hist)

# ============================================================
# Plot: per-chromosome median
# ============================================================

cis_chr <- cis_all[, .(
  median_cis = median(as.numeric(n_cis_variants), na.rm = TRUE)
), by = .(cohort, chr)]

p_cis_chr <- ggplot(cis_chr, aes(x = chr, y = median_cis, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.85)) +
  facet_wrap(~ cohort, ncol = 1, scales = "free_y") +
  labs(
    title = "Median cis-window variant count per chromosome (AFR strata)",
    x     = "Chromosome",
    y     = "Median variants per cis-window"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        legend.position = "none")
print(p_cis_chr)