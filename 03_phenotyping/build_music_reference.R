#!/usr/bin/env Rscript
# Script 1: Build MuSiC SingleCellExperiment reference from Campbell et al. 2023
# placental scRNA-seq data (GEO GSE182381).
#
# Input:  GSE182381_reference_sample.txt.gz (16,003 genes x 40,494 cells, 27 cell types)
# Output: placenta_music_reference.rds (SingleCellExperiment with pseudo-subjects)
#
# The reference matrix has cell-type labels as column headers but no biological-
# subject metadata (the integrated Seurat object with biorep IDs was not deposited
# publicly). MuSiC's cross-subject variance weighting requires >=2 subjects, so we
# create K=5 pseudo-subjects by balanced random split within each cell type.
#
# Usage:
#   Rscript build_music_reference.R \
#     --input  GSE182381_reference_sample.txt.gz \
#     --output placenta_music_reference.rds \
#     --n-pseudo-subjects 5 \
#     --seed 42

.libPaths(c("/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1", .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(Matrix)
  library(data.table)
})

option_list <- list(
  make_option("--input", type="character", default=NULL,
              help="Path to GSE182381_reference_sample.txt.gz [required]"),
  make_option("--output", type="character", default="placenta_music_reference.rds",
              help="Output .rds path for SingleCellExperiment"),
  make_option("--n-pseudo-subjects", type="integer", default=5,
              help="Number of pseudo-subjects to create per cell type (default 5)"),
  make_option("--seed", type="integer", default=42,
              help="Random seed for pseudo-subject assignment (default 42)")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$`input`)) {
  stop("--input is required. Provide path to GSE182381_reference_sample.txt.gz")
}

set.seed(opt$`seed`)
n_ps <- opt$`n-pseudo-subjects`

cat("=== Building MuSiC reference ===\n")
cat("Input: ", opt$`input`, "\n", sep="")
cat("Output:", opt$output, "\n")
cat("Pseudo-subjects:", n_ps, "| Seed:", opt$`seed`, "\n\n")

# --- Read reference matrix, preserving duplicate cell-type column names ---
# data.table::fread deduplicates column names, so we scan the header first
# to capture the true cell-type labels, then restore them after reading.
cat("Reading header (cell-type labels)...\n")
hdr <- scan(opt$`input`, what="", sep="\t", nlines=1, quiet=TRUE)
cell_types <- hdr[-1]  # drop first element ("GeneSymbol")
cat("  Cells:", length(cell_types), "| Unique cell types:", length(unique(cell_types)), "\n")

cat("Reading expression matrix (this takes ~30s)...\n")
ref <- fread(cmd=paste("zcat", opt$`input`), header=TRUE, sep="\t", data.table=FALSE)
genes <- ref[[1]]
ref <- as.matrix(ref[, -1])
rownames(ref) <- genes
colnames(ref) <- cell_types  # restore true labels (with duplicates)
cat("  Matrix:", paste(dim(ref), collapse=" x "), "\n")

# --- Create pseudo-subject assignments ---
# For each cell type, assign cells to n_ps pseudo-subjects by balanced random split.
cat("Assigning pseudo-subjects (K=", n_ps, ")...\n", sep="")
subject <- character(ncol(ref))
ps_labels <- paste0("ps", seq_len(n_ps))
for (ct in unique(cell_types)) {
  idx <- which(cell_types == ct)
  subject[idx] <- sample(rep(ps_labels, length.out=length(idx)))
}
subject <- factor(subject)
cat("  Subjects:", paste(levels(subject), collapse=", "), "\n")
cat("  Subject x cell-type table (first 5 types):\n")
print(table(subject, factor(cell_types))[, 1:5])

# --- Build SingleCellExperiment ---
cat("\nBuilding SingleCellExperiment...\n")
sce <- SingleCellExperiment(
  assays = list(counts = as(ref, "dgCMatrix")),
  colData = DataFrame(
    cell.type = factor(cell_types),
    subject   = subject
  )
)
cat("  SCE dim:", paste(dim(sce), collapse=" x "), "\n")
cat("  Cell types:", nlevels(sce$cell.type), "\n")
cat("  Subjects:", nlevels(sce$subject), "\n")
cat("  Cell-type counts:\n")
print(table(sce$cell.type))

# --- Save ---
cat("\nSaving to", opt$output, "...\n")
saveRDS(sce, opt$output)
cat("Done. Size:", round(file.size(opt$output) / 1e6, 1), "MB\n")
cat("\n=== Reference build complete ===\n")
