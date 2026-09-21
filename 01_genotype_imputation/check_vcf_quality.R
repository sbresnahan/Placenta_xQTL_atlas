.libPaths(c( "/home/stbresnahan/R/ubuntu/4.3.1" , .libPaths()))
library(ggplot2)
library(dplyr)
library(vcfR)



# ============================================================
# Before phasing/imputation
# ============================================================

# Load VCF file
vcf <- read.vcfR("/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/NIGMS/genotypes/raw/imputation_ready/NIGMS/NIGMS_filtered_norm.vcf.gz")

# Extract chromosome information
chr_data <- getCHROM(vcf)
pos_data <- getPOS(vcf)

# Create data frame
df <- data.frame(chr = chr_data, pos = pos_data)
df <- df[!grepl("random|chrUn|chrM|chrX|chrY", df$chr), ]

# Count SNPs per chromosome
snps_per_chr <- as.data.frame(table(df$chr))
colnames(snps_per_chr) <- c("chr", "count")

# Plot 1: SNPs per chromosome
p1 <- ggplot(snps_per_chr, aes(x = chr, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "SNPs per Chromosome",
       x = "Chromosome",
       y = "Number of SNPs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate chromosome lengths (max position per chr)
chr_lengths <- aggregate(pos ~ chr, data = df, FUN = max)
colnames(chr_lengths) <- c("chr", "length")

# Merge counts with lengths
snps_per_Mb <- merge(snps_per_chr, chr_lengths, by = "chr")
snps_per_Mb$snps_per_Mb <- (snps_per_Mb$count / snps_per_Mb$length) * 1000000

# Plot 2: SNPs per kb per chromosome
p2 <- ggplot(snps_per_Mb, aes(x = chr, y = snps_per_Mb)) +
  geom_bar(stat = "identity", fill = "coral") +
  theme_minimal() +
  labs(title = "SNPs per Mb per Chromosome",
       x = "Chromosome",
       y = "SNPs per Mb") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display plots
print(p1)
print(p2)


# ============================================================
# After phasing/imputation on SEAD
# ============================================================

# Load info file
info <- read.table("/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/mapqtl/SNUH/genotypes/imputed/chr1.SEAD.GRCh38.info", header = TRUE, stringsAsFactors = FALSE)
info$ALT_Frq <- as.numeric(info$ALT_Frq)
info$MAF <- as.numeric(info$MAF)
info$AvgCall <- as.numeric(info$AvgCall)
info$Rsq <- as.numeric(info$Rsq)

# Print QC summary
cat("=== Post-Imputation QC Summary ===\n")
cat("Total variants:", nrow(info), "\n")

# Imputation quality (Rsq)
cat("\n--- Imputation Quality (Rsq) ---\n")
cat("Mean Rsq:", mean(info$Rsq, na.rm = TRUE), "\n")
cat("Median Rsq:", median(info$Rsq, na.rm = TRUE), "\n")
cat("Min Rsq:", min(info$Rsq, na.rm = TRUE), "\n")
cat("Variants with Rsq < 0.3:", sum(info$Rsq < 0.3, na.rm = TRUE), "\n")
cat("Variants with Rsq < 0.5:", sum(info$Rsq < 0.5, na.rm = TRUE), "\n")
cat("Variants with Rsq < 0.8:", sum(info$Rsq < 0.8, na.rm = TRUE), "\n")

# MAF distribution
cat("\n--- Minor Allele Frequency ---\n")
cat("Mean MAF:", mean(info$MAF, na.rm = TRUE), "\n")
cat("Variants with MAF < 0.01:", sum(info$MAF < 0.01, na.rm = TRUE), "\n")
cat("Variants with MAF < 0.05:", sum(info$MAF < 0.05, na.rm = TRUE), "\n")

# Genotyped vs Imputed
cat("\n--- Genotyped vs Imputed ---\n")
cat("Genotyped variants:", sum(info$Genotyped == "Genotyped"), "\n")
cat("Imputed variants:", sum(info$Genotyped == "Imputed"), "\n")

# Generate QC plots
par(mfrow = c(2, 2))

# Plot 1: Rsq distribution
hist(info$Rsq, breaks = 50, main = "Imputation Quality (Rsq)", 
     xlab = "Rsq", col = "steelblue", xlim = c(0, 1))
abline(v = 0.8, col = "red", lty = 2, lwd = 2)
abline(v = 0.3, col = "orange", lty = 2, lwd = 2)
legend("topleft", legend = c("Rsq = 0.8", "Rsq = 0.3"), 
       col = c("red", "orange"), lty = 2, lwd = 2)

# Plot 2: Rsq by MAF
plot(info$MAF, info$Rsq, pch = 16, cex = 0.3, col = rgb(0, 0, 1, 0.1),
     main = "Rsq vs MAF", xlab = "MAF", ylab = "Rsq")
abline(h = 0.8, col = "red", lty = 2)

# Plot 3: MAF distribution
hist(info$MAF, breaks = 50, main = "Minor Allele Frequency",
     xlab = "MAF", col = "darkgreen")
abline(v = 0.01, col = "red", lty = 2)
abline(v = 0.05, col = "orange", lty = 2)

# Plot 4: Rsq by MAF bins
info$maf_bin <- cut(info$MAF, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 0.5))
boxplot(Rsq ~ maf_bin, data = info, 
        main = "Rsq by MAF Bins",
        xlab = "MAF Bin", ylab = "Rsq",
        col = "coral", las = 2)
abline(h = 0.8, col = "red", lty = 2)




# ============================================================
# After phasing/imputation on TOPMED
# ============================================================

library(data.table)
library(parallel)

read_topmed_info <- function(file, n_cores = parallel::detectCores() - 1, chunk_size = 500000) {
  
  # Get total line count (excluding header/comments)
  n <- as.integer(system(paste0("zcat ", file, " | grep -cv '^#'"), intern = TRUE))
  starts <- seq(1, n, by = chunk_size)
  
  process_chunk <- function(skip_n, nrows) {
    library(data.table)
    m <- fread(file, sep = "\t", header = FALSE, skip = skip_n,
               nrows = nrows, showProgress = FALSE)
    info_col <- m[[8]]
    info_split <- tstrsplit(info_col, ";", fixed = TRUE)
    
    extract_field <- function(key) {
      prefix <- paste0(key, "=")
      result <- rep(NA_real_, length(info_col))
      for (col in info_split) {
        hits <- !is.na(col) & startsWith(col, prefix)
        result[hits] <- as.numeric(substr(col[hits], nchar(prefix) + 1, nchar(col[hits])))
      }
      result
    }
    
    data.frame(
      CHROM     = m[[1]],
      POS       = m[[2]],
      ID        = m[[3]],
      REF       = m[[4]],
      ALT       = m[[5]],
      ALT_Frq   = extract_field("AF"),
      MAF       = extract_field("MAF"),
      AvgCall   = extract_field("AVG_CS"),
      Rsq       = extract_field("R2"),
      Genotyped = ifelse(grepl("TYPED", info_col, fixed = TRUE), "Genotyped", "Imputed"),
      stringsAsFactors = FALSE
    )
  }
  
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, { .libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths())); library(data.table) })
  
  result <- clusterMap(cl, fun = process_chunk,
                       skip_n = starts,
                       nrows  = pmin(chunk_size, n - starts + 1),
                       SIMPLIFY = FALSE)
  stopCluster(cl)
  
  rbindlist(result)
}

info <- read_topmed_info("/rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/imputed/MALI/chr1.info.gz",12)

# QC Summary

cat("\n=== Post-Imputation QC Summary ===\n")
cat("Total variants:", nrow(info), "\n")

cat("\n--- Imputation Quality (Rsq) ---\n")
cat("Mean Rsq:               ", mean(info$Rsq, na.rm = TRUE), "\n")
cat("Median Rsq:             ", median(info$Rsq, na.rm = TRUE), "\n")
cat("Min Rsq:                ", min(info$Rsq, na.rm = TRUE), "\n")
cat("Variants with Rsq > 0.5:", sum(info$Rsq > 0.5, na.rm = TRUE), "\n")
cat("Variants with Rsq > 0.8:", sum(info$Rsq > 0.8, na.rm = TRUE), "\n")

cat("\n--- Minor Allele Frequency ---\n")
cat("Mean MAF:                  ", mean(info$MAF, na.rm = TRUE), "\n")
cat("Variants with MAF > 0.05:  ", sum(info$MAF > 0.05, na.rm = TRUE), "\n")
cat("Variants with MAF > 0.01:  ", sum(info$MAF > 0.01, na.rm = TRUE), "\n")

cat("\n--- Genotyped vs Imputed ---\n")
cat("Genotyped variants:", sum(info$Genotyped == "Genotyped"), "\n")
cat("Imputed variants:  ", sum(info$Genotyped == "Imputed"), "\n")

cat("\n--- Variants passing QC thresholds (Rsq >= 0.8, MAF >= 0.05) ---\n")
pass <- sum(info$Rsq >= 0.8 & info$MAF >= 0.05, na.rm = TRUE)
cat("Passing variants:", pass, "(", round(100 * pass / nrow(info), 1), "% )\n")

# QC Plots

# Plot 1: Rsq distribution
hist(info$Rsq, breaks = 50, main = "Imputation Quality (Rsq)",
     xlab = "Rsq", col = "steelblue", xlim = c(0, 1))
abline(v = 0.8, col = "red",    lty = 2, lwd = 2)
abline(v = 0.3, col = "orange", lty = 2, lwd = 2)
legend("topleft", legend = c("Rsq = 0.8", "Rsq = 0.3"),
       col = c("red", "orange"), lty = 2, lwd = 2)

# Plot 2: MAF distribution
hist(info$MAF, breaks = 50, main = "Minor Allele Frequency",
     xlab = "MAF", col = "darkgreen")
abline(v = 0.01, col = "red",    lty = 2)
abline(v = 0.05, col = "orange", lty = 2)
legend("topright", legend = c("MAF = 0.01", "MAF = 0.05"),
       col = c("red", "orange"), lty = 2, lwd = 2)

# Check per LD block

blocks <- fread("/rsrch5/home/epi/stbresnahan/scratch/Placenta_QTL/imputation_pipeline/AFR_LD.bed")
info_pass <- info[info$Rsq >= 0.8 & info$MAF >= 0.01, ]
info_dt <- as.data.table(info_pass)
chr1_blocks <- blocks[blocks$chr == "chr1", ]
chr1_blocks <- chr1_blocks[order(chr1_blocks$start), ]
info_dt[, block := findInterval(POS, chr1_blocks$start)]
info_dt[, block := findInterval(POS, chr1_blocks$start)]
snps_per_block <- info_dt[, .N, by = block]
summary(snps_per_block$N)
