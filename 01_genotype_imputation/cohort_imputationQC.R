.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))
library(data.table)
library(ggplot2)
library(dplyr)

# ============================================================
# Configuration
# ============================================================

COHORT  <- "NIGMS"
OUTDIR <- paste("/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL",COHORT,
                "genotypes/imputed/qc",sep="/")
REPORTDIR <- paste0(OUTDIR,"/report")

RSQ_THRESH <- 0.4   # documentation only — filter was applied upstream

# ============================================================
# SECTION 1: Load ancestry assignments
# ============================================================

cat("\n=== Step 1: Loading PC-AiR ancestry assignments ===\n")
ancestry_df <- fread(file.path(OUTDIR, paste0(COHORT, "_ancestry_assignments.tsv")))
cat("Samples by assigned ancestry:\n")
print(table(ancestry_df$assigned_ancestry))

# ============================================================
# SECTION 2: Load Rsq-passing pvar (geometry only, no dosages)
# ============================================================

cat("\n=== Step 2: Loading Rsq-passing pvar ===\n")

pvar_file <- file.path(OUTDIR, paste0(COHORT, ".rsq_pass.pvar"))
pvar <- fread(cmd = paste0("grep -v '^##' ", pvar_file), sep = "\t")
setnames(pvar, 1, "CHROM")
pvar$CHROM <- paste0("chr", pvar$CHROM)

chr_order <- paste0("chr", 1:22)
pvar$CHROM <- factor(pvar$CHROM, levels = chr_order[chr_order %in% unique(pvar$CHROM)])

cat("Rsq-passing variants loaded:", nrow(pvar), "\n")

# ============================================================
# SECTION 3: Per-chromosome Rsq-passing variant counts
# ============================================================

CHR_LENGTHS <- c(
  chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
  chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
  chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
  chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
  chr17=83257441, chr18=80373285, chr19=58617616, chr20=64444167,
  chr21=46709983, chr22=50818468
)

chr_summary <- pvar[, .(n_rsq_pass = .N), by = CHROM]
chr_summary[, chr_len_mb  := CHR_LENGTHS[as.character(CHROM)] / 1e6]
chr_summary[, snps_per_mb := n_rsq_pass / chr_len_mb]
chr_summary[, cohort := COHORT]

fwrite(chr_summary,
       file.path(REPORTDIR, paste0(COHORT, "_rsq_pass_per_chr.tsv")),
       sep = "\t")
cat("Per-chromosome Rsq-passing summary written.\n")

# ============================================================
# SECTION 4: Plot — Rsq-passing SNPs per Mb per chromosome
# ============================================================

p1 <- ggplot(chr_summary, aes(x = CHROM, y = snps_per_mb)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title   = paste(COHORT, "— Rsq-passing SNPs per Mb per chromosome"),
    x       = "Chromosome",
    y       = "SNPs / Mb",
    caption = paste("Filter: Rsq >=", RSQ_THRESH)
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(REPORTDIR, paste0(COHORT, "_rsq_pass_snps_per_mb.png")),
       p1, width = 10, height = 4, dpi = 150)

cat("=== Per-cohort QC script complete ===\n")
cat("Outputs written to:", OUTDIR, "\n")
cat("  - ", COHORT, "_ancestry_assignments.tsv\n", sep = "")
cat("  - ", COHORT, "_rsq_pass_per_chr.tsv\n",     sep = "")
cat("  - ", COHORT, "_rsq_pass_snps_per_mb.png\n",  sep = "")
cat("\nNext step: run merge_and_filter_pooled.sh to merge within ancestry\n")
cat("stratum across cohorts and apply MAF/HWE/missingness filters.\n")