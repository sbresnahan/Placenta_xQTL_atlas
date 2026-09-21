.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))

## =============================================================================
## build_rnaseq_to_array_id_map.R
##
## Builds the RNA-seq -> genotype array sample ID map for the Placenta QTL
## project across 4 cohorts (NIEHS_RICHS, GUSTO, SNUH, NIGMS).
##
## Output columns: rnaseq_id, array_id, ancestry, cohort
##   - rnaseq_id : sample ID as it appears in the PANTRY transcriptomics BED files
##   - array_id  : sample ID as it appears in the pooled variant call files
##                 (psam/pgen), i.e. pooled_sample_ancestry.tsv$sample_id
##   - ancestry  : assigned_ancestry from pooled_sample_ancestry.tsv
##   - cohort    : NIEHS_RICHS / GUSTO / SNUH / NIGMS
##
## Verified linkage logic (against the actual files, 2026-09):
##   NIEHS_RICHS: key file Run (SRR) -> SUBJECT_ID_Array (S1, S2, ...) = pooled ID
##   GUSTO      : BED IDs are "J" + covars$ID (J1001-J1200);
##                pooled ID = paste0(SubjectID, "_B", SubjectID)
##   SNUH       : SRA metadata Run (SRR) -> isolate (OGF###) = pooled ID
##   NIGMS      : RNA-Seq runs with body_site == "placenta" (excludes DEC decidua;
##                keeps both fetal sampling sites -1 and -2) -> submitted_subject_id
##                -> WXS row -> WXS Run (SRR) = pooled ID
##
## Runs interactively in RStudio; base R only, no package dependencies.
## =============================================================================

## ---- Config: edit paths here -------------------------------------------------
pooled_file <- "/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/pooled/genotypes/pooled_sample_ancestry.tsv"

richs_key_file   <- "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/mapqtl/NIEHS_RICHS/genotypes/RICHS_GT_to_Seq_sampleKey.txt"
gusto_covars_file <- "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/mapqtl/GUSTO/covariates/20220823-Full_200_RNAseq_covars_v3.csv"
snuh_meta_file   <- "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/mapqtl/SNUH/RNA/SNUH_SRA_metadata.csv"
nigms_meta_file  <- "/rsrch5/home/epi/stbresnahan/bhattacharya_lab/data/mapqtl/NIGMS_SexDifferences_Placentas/NIGMS_SRA_metadata.csv"

## PANTRY cohort directories (for optional BED-header validation)
pantry_dir <- "/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/PANTRY"
cohort_dirs <- c(NIEHS_RICHS = file.path(pantry_dir, "cohort1"),
                 GUSTO       = file.path(pantry_dir, "cohort2"),
                 SNUH        = file.path(pantry_dir, "cohort3"),
                 NIGMS       = file.path(pantry_dir, "cohort4"))

out_file <- "/rsrch9/home/epi/bhattacharya_lab/data/Placenta_QTL/pooled/rnaseq_to_array_id_map.csv"

## ---- Helpers -----------------------------------------------------------------

## Join a per-cohort data.frame(rnaseq_id, array_key) to the pooled table.
## Prints match diagnostics; returns data.frame(rnaseq_id, array_id, ancestry, cohort).
join_to_pooled <- function(df, cohort_name, pooled) {
  pool <- pooled[pooled$cohort == cohort_name, ]
  idx  <- match(df$array_key, pool$sample_id)
  n_hit <- sum(!is.na(idx))
  cat(sprintf("[%s] RNA-seq samples: %d | matched to pooled genotypes: %d | unmatched: %d\n",
              cohort_name, nrow(df), n_hit, nrow(df) - n_hit))
  if (any(is.na(idx))) {
    cat(sprintf("[%s]   unmatched rnaseq_id(s): %s\n",
                cohort_name, paste(df$rnaseq_id[is.na(idx)], collapse = ", ")))
  }
  out <- data.frame(rnaseq_id = df$rnaseq_id[!is.na(idx)],
                    array_id  = pool$sample_id[idx[!is.na(idx)]],
                    ancestry  = pool$assigned_ancestry[idx[!is.na(idx)]],
                    cohort    = cohort_name,
                    stringsAsFactors = FALSE)
  if (anyDuplicated(out$rnaseq_id))
    stop(sprintf("[%s] duplicate rnaseq_id after join - investigate before proceeding", cohort_name))
  if (any(is.na(out$array_id)) || any(is.na(out$ancestry)))
    stop(sprintf("[%s] NA array_id/ancestry after join - investigate before proceeding", cohort_name))
  out
}

## Read the sample IDs from the header of a PANTRY expression.bed.gz and
## report what fraction are covered by the map.
check_bed_coverage <- function(cohort_name, bed_path, map_ids) {
  if (!file.exists(bed_path)) {
    cat(sprintf("[%s] BED not found (%s) - skipping coverage check\n", cohort_name, bed_path))
    return(invisible(NULL))
  }
  con <- gzfile(bed_path, "rt")
  hdr <- readLines(con, n = 1)
  close(con)
  samples <- strsplit(hdr, "\\s+")[[1]]
  samples <- setdiff(samples, c("#chr", "start", "end", "phenotype_id"))
  in_map <- samples %in% map_ids
  cat(sprintf("[%s] BED header: %d samples | covered by map: %d (%.1f%%)\n",
              cohort_name, length(samples), sum(in_map), 100 * mean(in_map)))
  if (any(!in_map))
    cat(sprintf("[%s]   BED samples missing from map: %s\n",
                cohort_name, paste(samples[!in_map], collapse = ", ")))
  invisible(samples)
}

## ---- Load pooled genotype sample table ---------------------------------------
pooled <- read.delim(pooled_file, stringsAsFactors = FALSE)
stopifnot(all(c("sample_id", "assigned_ancestry", "cohort") %in% names(pooled)))
cat(sprintf("Pooled genotype table: %d samples across cohorts: %s\n\n",
            nrow(pooled), paste(names(table(pooled$cohort)), as.integer(table(pooled$cohort)),
                                sep = "=", collapse = ", ")))

pooled <- pooled[!pooled$cohort=="MALI_G3A",]

## ---- Cohort 1: NIEHS_RICHS ---------------------------------------------------
## Key file: dbGaP_Subject_ID  SUBJECT_ID_Array  SUBJECT_ID_Seq  Run
richs_key <- read.table(richs_key_file, header = TRUE, stringsAsFactors = FALSE)
stopifnot(all(c("SUBJECT_ID_Array", "Run") %in% names(richs_key)))
richs_df <- data.frame(rnaseq_id = richs_key$Run,
                       array_key = richs_key$SUBJECT_ID_Array,
                       stringsAsFactors = FALSE)
map_richs <- join_to_pooled(richs_df, "NIEHS_RICHS", pooled)

## ---- Cohort 2: GUSTO ----------------------------------------------------------
## BED sample IDs are "J" + covars$ID (J1001-J1200).
## Pooled genotype ID = SubjectID_B<SubjectID> (e.g. 010-04002_B010-04002).
gusto <- read.csv(gusto_covars_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
stopifnot(all(c("SubjectID", "ID") %in% names(gusto)))
gusto_df <- data.frame(rnaseq_id = paste0("J", gusto$ID),
                       array_key = paste0(gusto$SubjectID, "_B", gusto$SubjectID),
                       stringsAsFactors = FALSE)
map_gusto <- join_to_pooled(gusto_df, "GUSTO", pooled)

## ---- Cohort 3: SNUH -----------------------------------------------------------
## SRA metadata: Run (SRR) -> isolate (OGF###) = pooled ID.
snuh <- read.csv(snuh_meta_file, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("Run", "isolate") %in% names(snuh)))
snuh <- snuh[snuh[["Assay Type"]] == "RNA-Seq", ]
snuh_df <- data.frame(rnaseq_id = snuh$Run,
                      array_key = snuh$isolate,
                      stringsAsFactors = FALSE)
map_snuh <- join_to_pooled(snuh_df, "SNUH", pooled)

## ---- Cohort 4: NIGMS ----------------------------------------------------------
## Fetal placenta RNA-Seq runs only (body_site == "placenta"; excludes DEC decidua).
## Link: RNA submitted_subject_id (OBG...) -> WXS row -> WXS Run = pooled ID.
## Note: each subject has two placenta sampling sites (PLAC-RNA1/-2), so each
## array_id legitimately appears twice in this cohort's map.
nigms <- read.csv(nigms_meta_file, stringsAsFactors = FALSE, check.names = FALSE)
nigms_rna <- nigms[nigms[["Assay Type"]] == "RNA-Seq" & nigms$body_site == "placenta", ]
nigms_wxs <- nigms[nigms[["Assay Type"]] == "WXS", c("submitted_subject_id", "Run")]
stopifnot(!anyDuplicated(nigms_wxs$submitted_subject_id))
nigms_df <- data.frame(rnaseq_id = nigms_rna$Run,
                       array_key = nigms_wxs$Run[match(nigms_rna$submitted_subject_id,
                                                       nigms_wxs$submitted_subject_id)],
                       stringsAsFactors = FALSE)
no_wxs <- is.na(nigms_df$array_key)
if (any(no_wxs))
  cat(sprintf("[NIGMS] %d placenta RNA run(s) dropped - subject has no WXS exome: %s\n",
              sum(no_wxs), paste(nigms_df$rnaseq_id[no_wxs], collapse = ", ")))
nigms_df <- nigms_df[!no_wxs, ]
map_nigms <- join_to_pooled(nigms_df, "NIGMS", pooled)

## ---- Combine, validate, write -------------------------------------------------
id_map <- rbind(map_richs, map_gusto, map_snuh, map_nigms)
id_map <- id_map[, c("rnaseq_id", "array_id", "ancestry", "cohort")]

stopifnot(!anyDuplicated(id_map$rnaseq_id))

cat("\nFinal map:\n")
print(table(id_map$cohort))
cat("\nAncestry x cohort:\n")
print(table(id_map$cohort, id_map$ancestry))

write.csv(id_map, out_file, row.names = FALSE, quote = FALSE)
cat(sprintf("\nWrote %d rows to:\n  %s\n", nrow(id_map), out_file))

## ---- Optional: validate against PANTRY BED headers ----------------------------
cat("\n--- BED header coverage ---\n")
for (co in names(cohort_dirs)) {
  bed <- file.path(cohort_dirs[co], "output", "expression.bed.gz")
  check_bed_coverage(co, bed, id_map$rnaseq_id[id_map$cohort == co])
}

cat("\nNote: NIGMS array_ids each appear twice (two fetal placenta sampling sites\n",
    "per subject, PLAC-RNA1 and PLAC-RNA2). Choose one run per subject downstream\n",
    "before QTL analysis (e.g. the run present in the cohort4 BED files).\n", sep = "")
