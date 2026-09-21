# =============================================================================
# prep_PANTRY_gtf.R
# -----------------------------------------------------------------------------
# Normalize a SQANTI3-annotated GTF to be compatible with PANTRY.
#
# PANTRY was designed for GENCODE-style GTFs. The long-read assembly
# annotated by SQANTI3 only has `transcript` and `exon` features with attributes
# gene_id, transcript_id, EAS_term, EUR_preterm, EUR_term, exon_number. It is
# missing several things PANTRY scripts require:
#
#   1. `gene` feature lines          (map_edit_sites_to_genes.py, preprocess_gtf.py)
#   2. `gene_biotype` attribute      (preprocess_gtf.py, prepareAnnotations.R)
#   3. `transcript_biotype` attribute(prepareAnnotations.R)
#   4. `gene_name` attribute         (gtf_to_majiq_gff3.py, nice to have)
#   5. `CDS` feature lines           (prepareAnnotations.R cdsBy)
#   6. `tag "basic"` attribute       (extractTranscriptTags.py)
#
# This script adds all six using the SQANTI3 classification table as the source
# of biotype (coding/non_coding), gene name (associated_gene), and CDS coordinates
# (CDS_genomic_start/end).
#
# Usage: Rscript prep_PANTRY_gtf.R --gtf <raw.gtf> \
#          --classification <classification.txt> \
#          --reference <reference.gtf> --out <normalized.gtf> \
#          [--source-prefix HPLR]
# =============================================================================

# ---- Library path (seadragon R package library) ----
.libPaths(c("/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1", .libPaths()))

# ---- Command-line arguments ----
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--gtf", type = "character", default = NULL,
              help = "Raw SQANTI3-annotated GTF (input)"),
  make_option("--classification", type = "character", default = NULL,
              help = "SQANTI3 classification table (input)"),
  make_option("--reference", type = "character", default = NULL,
              help = "Reference GTF for Ensembl gene records (input)"),
  make_option("--out", type = "character", default = NULL,
              help = "Normalized PANTRY-compatible GTF to write (output)"),
  make_option("--source-prefix", type = "character", default = "HPLR",
              help = "Source label written into GTF column 2 for gene feature lines whose input line has no source, and into the ##source header comment. Default: HPLR")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("gtf", "classification", "reference", "out")
missing <- required[vapply(required, function(k) is.null(opt[[k]]), logical(1))]
if (length(missing) > 0) {
  stop(sprintf("Missing required argument(s): %s. Required: --gtf, --classification, --reference, --out",
               paste0("--", missing, collapse = ", ")))
}

gtf_file <- opt$gtf
classification_file <- opt$classification
output_file <- opt$out
REFERENCE_GTF <- opt$reference
SOURCE_PREFIX <- opt[["source-prefix"]]

# ---- Parameters ----
CODING_BIOTYPE <- "protein_coding"
NONCODING_BIOTYPE <- "lncRNA"
DEFAULT_BIOTYPE <- "protein_coding"

# ---- Helper functions ----

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

parse_attributes <- function(attr_str) {
  attr_str <- sub(";\\s*$", "", trimws(attr_str))
  if (attr_str == "") return(character(0))
  parts <- strsplit(attr_str, ";\\s*")[[1]]
  attrs <- character(0)
  for (p in parts) {
    p <- trimws(p)
    if (p == "") next
    sp <- regmatches(p, regexec("^(\\S+)\\s+(.*)$", p))[[1]]
    if (length(sp) == 3) {
      key <- sp[2]
      val <- gsub('^"|"$', "", sp[3])
      attrs[key] <- val
    } else {
      attrs[p] <- ""
    }
  }
  attrs
}

format_attributes <- function(attrs) {
  attrs <- attrs[!is.na(attrs)]
  if (length(attrs) == 0) return("")
  
  # Enforce Ensembl attribute order for the keys we emit
  ensembl_order <- c("gene_id", "gene_name", "gene_biotype",
                     "transcript_id", "transcript_biotype", "tag")
  present <- names(attrs)
  ordered_keys <- c(intersect(ensembl_order, present),
                    setdiff(present, ensembl_order))
  attrs <- attrs[ordered_keys]
  
  parts <- character(length(attrs))
  for (i in seq_along(attrs)) {
    key <- names(attrs)[i]
    val <- attrs[i]
    if (val == "") {
      parts[i] <- key
    } else {
      parts[i] <- paste0(key, ' "', val, '"')
    }
  }
  paste0(paste(parts, collapse = "; "), ";")
}

find_column <- function(df, candidates, description) {
  for (c in candidates) {
    if (c %in% colnames(df)) return(c)
  }
  stop(sprintf("Could not find %s column. Looked for: %s. Available: %s",
               description, paste(candidates, collapse=", "),
               paste(colnames(df), collapse=", ")))
}

cat("================================================================\n")
cat("prep_PANTRY_gtf.R\n")
cat("================================================================\n")
cat("GTF:               ", gtf_file, "\n")
cat("Classification:    ", classification_file, "\n")
cat("Reference GTF:     ", REFERENCE_GTF, "\n")
cat("Output:            ", output_file, "\n")
cat("Source prefix:     ", SOURCE_PREFIX, "\n")
cat("================================================================\n\n")

# =============================================================================
# 1. Read SQANTI3 classification table
# =============================================================================
cat("[1/7] Reading SQANTI3 classification table...\n")
classification <- read.delim(classification_file, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, quote = "")
cat("  Rows:", nrow(classification), "\n")

isoform_col <- find_column(classification, c("isoform", "transcript_id"), "isoform ID")
coding_col <- find_column(classification, c("coding", "coding_potential"), "coding status")
gene_name_col <- find_column(classification, c("associated_gene", "associatedGene"), "associated gene")
cds_start_col <- find_column(classification, c("CDS_genomic_start", "cds_genomic_start", "CDS_start"), "CDS genomic start")
cds_end_col <- find_column(classification, c("CDS_genomic_end", "cds_genomic_end", "CDS_end"), "CDS genomic end")

cat("  Columns identified:\n")
cat("    isoform:        ", isoform_col, "\n")
cat("    coding:         ", coding_col, "\n")
cat("    associated_gene:", gene_name_col, "\n")
cat("    CDS_start:      ", cds_start_col, "\n")
cat("    CDS_end:        ", cds_end_col, "\n")

classification$isoform_id <- classification[[isoform_col]]

coding_status <- classification[[coding_col]]
transcript_biotype <- ifelse(tolower(coding_status) == "coding",
                             CODING_BIOTYPE, NONCODING_BIOTYPE)
transcript_biotype[is.na(transcript_biotype) | transcript_biotype == ""] <- DEFAULT_BIOTYPE
names(transcript_biotype) <- classification$isoform_id

gene_name_map <- classification[[gene_name_col]]
names(gene_name_map) <- classification$isoform_id

cds_start_map <- classification[[cds_start_col]]
cds_end_map <- classification[[cds_end_col]]
names(cds_start_map) <- classification$isoform_id
names(cds_end_map) <- classification$isoform_id

cat("  Coding transcripts:    ", sum(transcript_biotype == CODING_BIOTYPE), "\n")
cat("  Non-coding transcripts:", sum(transcript_biotype == NONCODING_BIOTYPE), "\n\n")

# =============================================================================
# 2. Read GTF
# =============================================================================
cat("[2/7] Reading GTF...\n")
gtf_lines <- readLines(gtf_file)
cat("  Total lines:", length(gtf_lines), "\n")

header_mask <- grepl("^#", gtf_lines)
header_lines <- gtf_lines[header_mask]
feature_lines <- gtf_lines[!header_mask]
cat("  Header lines:", length(header_lines), "\n")
cat("  Feature lines:", length(feature_lines), "\n\n")

# =============================================================================
# 3. Parse all feature lines
# =============================================================================
cat("[3/7] Parsing GTF feature lines...\n")
fields_list <- strsplit(feature_lines, "\t")
n_fields <- sapply(fields_list, length)
if (any(n_fields != 9)) {
  bad <- which(n_fields != 9)
  stop(sprintf("%d feature lines do not have 9 tab-separated fields. First bad line: %s",
               length(bad), feature_lines[bad[1]]))
}

gtf_df <- data.frame(
  seqname = sapply(fields_list, `[`, 1),
  source = sapply(fields_list, `[`, 2),
  feature = sapply(fields_list, `[`, 3),
  start = as.integer(sapply(fields_list, `[`, 4)),
  end = as.integer(sapply(fields_list, `[`, 5)),
  score = sapply(fields_list, `[`, 6),
  strand = sapply(fields_list, `[`, 7),
  frame = sapply(fields_list, `[`, 8),
  stringsAsFactors = FALSE
)
attr_raw <- sapply(fields_list, `[`, 9)
gtf_df$attrs <- lapply(attr_raw, parse_attributes)

cat("  Feature types:\n")
ft <- table(gtf_df$feature)
for (f in names(ft)) cat("    ", f, ":", ft[f], "\n")
cat("\n")

gtf_df$gene_id <- sapply(gtf_df$attrs, function(a) a["gene_id"])
gtf_df$gene_id[is.na(gtf_df$gene_id)] <- ""
gtf_df$transcript_id <- sapply(gtf_df$attrs, function(a) a["transcript_id"])
gtf_df$transcript_id[is.na(gtf_df$transcript_id)] <- ""

# =============================================================================
# 4. Derive gene boundaries and gene-level biotype
# =============================================================================
cat("[4/7] Deriving gene boundaries and biotypes...\n")

# --- Load reference GTF gene records -----------------------------------------
cat("  Loading reference gene records from:", REFERENCE_GTF, "\n")
ref_gene_lines <- grep("\tgene\t", readLines(REFERENCE_GTF), value = TRUE)
ref_split <- strsplit(ref_gene_lines, "\t", fixed = TRUE)
ref_v <- do.call(rbind, ref_split)
get_attr <- function(x, key) {
  pat <- paste0(key, ' "([^"]*)"')
  out <- rep(NA_character_, length(x))
  hit <- grepl(pat, x)
  out[hit] <- sub(paste0(".*", pat, ".*"), "\\1", x[hit])
  out
}
ref_gene <- data.frame(
  gene_id   = get_attr(ref_v[, 9], "gene_id"),
  seqname   = ref_v[, 1],
  start     = as.integer(ref_v[, 4]),
  end       = as.integer(ref_v[, 5]),
  strand    = ref_v[, 7],
  gene_type = get_attr(ref_v[, 9], "gene_type"),
  gene_name = get_attr(ref_v[, 9], "gene_name"),
  stringsAsFactors = FALSE
)
ref_gene <- ref_gene[!duplicated(ref_gene$gene_id), ]
# Collapse GENCODE gene_type onto the two backbone biotypes accepted by
# txrevise (complete_transcripts). protein_coding -> CODING, all else ->
# NONCODING. NONCODING_BIOTYPE must be one of the noncoding labels your
# txrevise accepts (lincRNA / lncRNA / lnc_RNA).
ref_gene$gene_biotype <- ifelse(ref_gene$gene_type == "protein_coding",
                                CODING_BIOTYPE, NONCODING_BIOTYPE)
cat("  Reference genes loaded:", nrow(ref_gene), "\n")

# --- Derive boundaries/biotype/name for all genes (overridden below for ENSG) -
tx_mask <- gtf_df$feature == "transcript"
tx_df <- gtf_df[tx_mask, ]
gene_min_start <- tapply(tx_df$start, tx_df$gene_id, min)
gene_max_end <- tapply(tx_df$end, tx_df$gene_id, max)
gene_chrom <- tapply(tx_df$seqname, tx_df$gene_id, function(x) x[1])
gene_strand <- tapply(tx_df$strand, tx_df$gene_id, function(x) x[1])
gene_ids <- names(gene_min_start)
gene_df <- data.frame(
  gene_id = gene_ids,
  seqname = gene_chrom[gene_ids],
  start = as.integer(gene_min_start[gene_ids]),
  end = as.integer(gene_max_end[gene_ids]),
  strand = gene_strand[gene_ids],
  stringsAsFactors = FALSE
)
rownames(gene_df) <- NULL
gene_to_txs <- split(tx_df$transcript_id, tx_df$gene_id)
gene_biotype <- sapply(gene_ids, function(gid) {
  txs <- gene_to_txs[[gid]]
  bt <- transcript_biotype[txs]
  bt <- bt[!is.na(bt)]
  if (length(bt) == 0) return(DEFAULT_BIOTYPE)
  if (any(bt == CODING_BIOTYPE)) return(CODING_BIOTYPE)
  return(NONCODING_BIOTYPE)
})
gene_df$gene_biotype <- gene_biotype[gene_ids]
gene_name <- sapply(gene_ids, function(gid) {
  txs <- gene_to_txs[[gid]]
  gnames <- gene_name_map[txs]
  gnames <- gnames[!is.na(gnames) & gnames != "" & gnames != "novel"]
  if (length(gnames) > 0) return(gnames[1])
  return(gid)
})
gene_df$gene_name <- gene_name[gene_ids]

# --- Override standard Ensembl genes with reference GTF values ----------------
# A standard Ensembl gene_id is a single versioned ENSG accession with no extra
# suffix. Fusion products (e.g. ENSGxxx_ENSGyyy) contain an underscore and are
# intentionally excluded, so their boundaries/biotype/name are derived above.
is_ensembl <- grepl("^ENSG[0-9]+\\.[0-9]+$", gene_df$gene_id)

ens_missing <- setdiff(gene_df$gene_id[is_ensembl], ref_gene$gene_id)
if (length(ens_missing) > 0) {
  stop(sprintf(
    "[4/7] %d Ensembl gene_id(s) not found in reference GTF (%s). Examples: %s",
    length(ens_missing), REFERENCE_GTF,
    paste(head(ens_missing, 10), collapse = ", ")
  ))
}

ref_idx <- match(gene_df$gene_id[is_ensembl], ref_gene$gene_id)
gene_df$seqname[is_ensembl]      <- ref_gene$seqname[ref_idx]
gene_df$start[is_ensembl]        <- ref_gene$start[ref_idx]
gene_df$end[is_ensembl]          <- ref_gene$end[ref_idx]
gene_df$strand[is_ensembl]       <- ref_gene$strand[ref_idx]
gene_df$gene_biotype[is_ensembl] <- ref_gene$gene_biotype[ref_idx]
gene_df$gene_name[is_ensembl]    <- ref_gene$gene_name[ref_idx]

cat("  Genes derived:", nrow(gene_df), "\n")
cat("  Ensembl genes (from reference):", sum(is_ensembl), "\n")
cat("  Novel/fusion genes (derived):", sum(!is_ensembl), "\n")
cat("  Protein-coding genes:", sum(gene_df$gene_biotype == CODING_BIOTYPE), "\n")
cat("  lncRNA genes:", sum(gene_df$gene_biotype == NONCODING_BIOTYPE), "\n\n")

# =============================================================================
# 5. Build CDS feature lines from SQANTI3 coordinates
# =============================================================================
cat("[5/7] Building CDS features...\n")

cds_records <- list()
exon_mask <- gtf_df$feature == "exon"
exon_df <- gtf_df[exon_mask, ]

for (i in seq_len(nrow(tx_df))) {
  tx_id <- tx_df$transcript_id[i]
  if (tx_id == "") next
  
  cds_s <- cds_start_map[tx_id]
  cds_e <- cds_end_map[tx_id]
  if (is.na(cds_s) || is.na(cds_e)) next
  cds_s <- as.integer(cds_s)
  cds_e <- as.integer(cds_e)
  if (is.na(cds_s) || is.na(cds_e)) next
  if (cds_s >= cds_e) next
  
  tx_exons <- exon_df[exon_df$transcript_id == tx_id, ]
  if (nrow(tx_exons) == 0) next
  
  for (j in seq_len(nrow(tx_exons))) {
    exon_start <- tx_exons$start[j]
    exon_end <- tx_exons$end[j]
    cds_seg_start <- max(cds_s, exon_start)
    cds_seg_end <- min(cds_e, exon_end)
    if (cds_seg_start <= cds_seg_end) {
      cds_records[[length(cds_records) + 1]] <- list(
        seqname = tx_exons$seqname[j],
        source = tx_exons$source[j],
        start = cds_seg_start,
        end = cds_seg_end,
        strand = tx_exons$strand[j],
        gene_id = tx_exons$gene_id[j],
        transcript_id = tx_id
      )
    }
  }
}

cat("  CDS segments created:", length(cds_records), "\n\n")

# =============================================================================
# 6. Assemble the normalized GTF
# =============================================================================
cat("[6/7] Assembling normalized GTF...\n")

# Precompute gene lookup (named indices) to replace repeated which()
gene_idx_lookup <- setNames(seq_len(nrow(gene_df)), gene_df$gene_id)

# 6b. Gene feature lines (vectorized)
gene_attr_strings <- vapply(seq_len(nrow(gene_df)), function(i) {
  format_attributes(c(
    gene_id = gene_df$gene_id[i],
    gene_biotype = gene_df$gene_biotype[i],
    gene_name = gene_df$gene_name[i]
  ))
}, character(1))

gene_lines <- paste(
  gene_df$seqname,
  gene_df$source %||% SOURCE_PREFIX,
  "gene",
  gene_df$start,
  gene_df$end,
  ".",
  gene_df$strand,
  ".",
  gene_attr_strings,
  sep = "\t"
)
cat("  Gene lines:", nrow(gene_df), "\n")

# 6c. Transcript and exon lines
keep <- gtf_df$feature != "gene"
sub <- gtf_df[keep, ]
n_rows <- nrow(sub)
feat_lines <- character(n_rows)

for (i in seq_len(n_rows)) {
  feat <- sub$feature[i]
  attrs <- sub$attrs[[i]]
  gid <- sub$gene_id[i]
  tid <- sub$transcript_id[i]
  
  idx <- gene_idx_lookup[[gid]]  # NULL if not present
  if (!is.null(idx)) {
    attrs["gene_biotype"] <- gene_df$gene_biotype[idx]
    attrs["gene_name"] <- gene_df$gene_name[idx]
  }
  
  if (feat == "transcript") {
    if (tid != "") {
      bt <- transcript_biotype[tid]
      if (is.na(bt) || bt == "") bt <- DEFAULT_BIOTYPE
      attrs["transcript_biotype"] <- bt
    }
    attrs["tag"] <- "basic"
  } else if (feat == "exon") {
    if (tid != "") {
      bt <- transcript_biotype[tid]
      if (is.na(bt) || bt == "") bt <- DEFAULT_BIOTYPE
      attrs["transcript_biotype"] <- bt
    }
  }
  
  feat_lines[i] <- paste(c(
    sub$seqname[i], sub$source[i], sub$feature[i],
    sub$start[i], sub$end[i], sub$score[i],
    sub$strand[i], sub$frame[i],
    format_attributes(attrs)
  ), collapse = "\t")
}
n_tx <- sum(sub$feature == "transcript")
n_exon <- sum(sub$feature == "exon")
cat("  Transcript lines:", n_tx, "\n")
cat("  Exon lines:", n_exon, "\n")

# 6d. CDS feature lines
n_cds <- length(cds_records)
cds_lines <- character(n_cds)
for (j in seq_len(n_cds)) {
  rec <- cds_records[[j]]
  attrs <- c(gene_id = rec$gene_id, transcript_id = rec$transcript_id)
  idx <- gene_idx_lookup[[rec$gene_id]]
  if (!is.null(idx)) {
    attrs["gene_biotype"] <- gene_df$gene_biotype[idx]
    attrs["gene_name"] <- gene_df$gene_name[idx]
  }
  if (rec$transcript_id != "") {
    bt <- transcript_biotype[rec$transcript_id]
    if (is.na(bt) || bt == "") bt <- DEFAULT_BIOTYPE
    attrs["transcript_biotype"] <- bt
  }
  cds_lines[j] <- paste(c(
    rec$seqname, rec$source, "CDS",
    rec$start, rec$end, ".",
    rec$strand, ".",
    format_attributes(attrs)
  ), collapse = "\t")
}
cat("  CDS lines:", n_cds, "\n")

# Build a sortable table of all feature lines, then order Ensembl-style.
# Sort keys: seqname, gene start, gene_id, transcript start, transcript_id,
# feature-type rank, feature start.

feature_rank <- c(gene = 1L, transcript = 2L, exon = 3L, CDS = 4L)

# Per-gene sort anchors
gene_start_lookup  <- setNames(gene_df$start,   gene_df$gene_id)
gene_seq_lookup    <- setNames(gene_df$seqname, gene_df$gene_id)

# Per-transcript sort anchors (min start of each transcript's rows in gtf_df)
tx_start_lookup <- tapply(sub$start, sub$transcript_id, min)

get_gene_start <- function(gid) {
  v <- gene_start_lookup[[gid]]; if (is.null(v)) NA_integer_ else v
}
get_gene_seq <- function(gid) {
  v <- gene_seq_lookup[[gid]]; if (is.null(v)) NA_character_ else v
}
get_tx_start <- function(tid) {
  if (tid == "" || is.na(tid)) return(NA_integer_)
  v <- tx_start_lookup[[tid]]; if (is.null(v)) NA_integer_ else v
}

# Gene records
gene_keys <- data.frame(
  line       = gene_lines,
  seqname    = gene_df$seqname,
  gene_start = gene_df$start,
  gene_id    = gene_df$gene_id,
  tx_start   = -1L,               # gene sorts before its transcripts
  tx_id      = "",
  frank      = feature_rank[["gene"]],
  fstart     = gene_df$start,
  stringsAsFactors = FALSE
)

# Transcript/exon records
feat_keys <- data.frame(
  line       = feat_lines,
  seqname    = sub$seqname,
  gene_start = vapply(sub$gene_id, get_gene_start, integer(1)),
  gene_id    = sub$gene_id,
  tx_start   = vapply(sub$transcript_id, get_tx_start, integer(1)),
  tx_id      = sub$transcript_id,
  frank      = feature_rank[sub$feature],
  fstart     = sub$start,
  stringsAsFactors = FALSE
)

# CDS records
cds_gene_id <- vapply(cds_records, function(r) r$gene_id,       character(1))
cds_tx_id   <- vapply(cds_records, function(r) r$transcript_id, character(1))
cds_start   <- vapply(cds_records, function(r) r$start,         integer(1))
cds_seqname <- vapply(cds_records, function(r) r$seqname,       character(1))

cds_keys <- data.frame(
  line       = cds_lines,
  seqname    = cds_seqname,
  gene_start = vapply(cds_gene_id, get_gene_start, integer(1)),
  gene_id    = cds_gene_id,
  tx_start   = vapply(cds_tx_id, get_tx_start, integer(1)),
  tx_id      = cds_tx_id,
  frank      = feature_rank[["CDS"]],
  fstart     = cds_start,
  stringsAsFactors = FALSE
)

all_keys <- rbind(gene_keys, feat_keys, cds_keys)

ord <- order(
  all_keys$seqname,
  all_keys$gene_start,
  all_keys$gene_id,
  all_keys$tx_start,
  all_keys$tx_id,
  all_keys$frank,
  all_keys$fstart
)

output_lines <- c(
  "##gff-version 2",
  paste0("##source: ", SOURCE_PREFIX, " normalized for PANTRY (", date(), ")"),
  all_keys$line[ord]
)

# =============================================================================
# 7. Write output
# =============================================================================
cat("[7/7] Writing normalized GTF...\n")
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
writeLines(output_lines, output_file)
cat("  Written:", output_file, "\n")
cat("  Total lines:", length(output_lines), "\n\n")

cat("================================================================\n")
cat("DONE. Normalized GTF written to:\n")
cat("  ", output_file, "\n")
cat("================================================================\n")
