#!/usr/bin/env Rscript
# =============================================================================
# 44_colocboost.R — multi-trait colocBoost worker for one shard of regions
# =============================================================================
# Objective 1.6: joint multi-trait colocalization per region x ancestry.
# Regions and their outcome manifests come from 43_prepare_colocboost.py.
# Each region stacks every overlapping fine-mapped xQTL phenotype (all 9
# modalities) with every ancestry-matched GWAS and runs colocBoost in
# summary-statistic mode.
#
# LD is supplied as reference genotype matrices (X_ref mode): in-sample
# genotypes for xQTL outcomes, 1KG superpopulation genotypes for GWAS
# outcomes (plink2 --export A dosages, column-mean imputed). X_ref avoids
# materializing P x P LD matrices for large merged regions (N << P).
#
# Outputs (per region, under {outdir}/results/{ancestry}/):
#   {region_id}.clusters.tsv   colocBoost CoS summary (colocalized outcome
#                              sets, purity, top variant, VCP) + metadata
#   {region_id}.vcp.tsv        per-variant VCP
#   {region_id}.done           sentinel
# Per shard:
#   {outdir}/diagnostics/{ancestry}.shard-{idx}.diagnostics.tsv
#
# Usage:
#   Rscript 44_colocboost.R --regions EAS.regions.tsv --outcomes EAS.outcomes.tsv \
#     --shard-index 1 --n-shards 20 --outdir $COLOC_DIR/colocboost \
#     --ld-xqtl-pgen $QTL_DIR/EAS_qtl --ld-gwas-pgen $KG_PGEN \
#     --ld-gwas-keep $COLOC_DIR/loci/EAS.1kg.keep
# =============================================================================

R_LIB <- "/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"
if (!dir.exists(R_LIB)) stop("R package library does not exist: ", R_LIB)
.libPaths(c(R_LIB, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(colocboost)
})

option_list <- list(
  make_option("--regions", type = "character", help = "{ANC}.regions.tsv"),
  make_option("--outcomes", type = "character", help = "{ANC}.outcomes.tsv"),
  make_option("--shard-index", type = "integer", default = 1L),
  make_option("--n-shards", type = "integer", default = 1L),
  make_option("--outdir", type = "character"),
  make_option("--plink2", type = "character", default = "plink2"),
  make_option("--ld-xqtl-pgen", type = "character"),
  make_option("--ld-gwas-pgen", type = "character"),
  make_option("--ld-gwas-keep", type = "character", default = ""),
  make_option("--min-variants", type = "integer", default = 50L),
  make_option("--min-outcomes", type = "integer", default = 2L),
  make_option("--M", type = "integer", default = 500L,
              help = "max boosting rounds per outcome [default %default]"),
  make_option("--tmp-dir", type = "character", default = tempdir())
)
opt <- parse_args(OptionParser(option_list = option_list))

.args_all <- commandArgs(trailingOnly = FALSE)
.file_arg <- sub("^--file=", "", grep("^--file=", .args_all, value = TRUE))
.scripts_dir <- Sys.getenv("SCRIPTS_DIR",
                           unset = if (length(.file_arg)) dirname(.file_arg[1]) else getwd())
source(file.path(.scripts_dir, "coloc_common.R"))

dir.create(file.path(opt$outdir, "diagnostics"), showWarnings = FALSE, recursive = TRUE)

regions <- fread(opt$regions)
outcomes <- fread(opt$outcomes)
ancestry <- unique(regions$ancestry)
if (length(ancestry) != 1L) stop("regions file mixes ancestries")
diag_path <- file.path(opt$outdir, "diagnostics",
                       sprintf("%s.shard-%04d.diagnostics.tsv", ancestry,
                               opt$`shard-index`))

shard_rows <- which(((seq_len(nrow(regions)) - 1L) %% opt$`n-shards`) + 1L ==
                    opt$`shard-index`)
regions <- regions[shard_rows, ]
message(sprintf("[shard %d/%d] %d regions (%s)", opt$`shard-index`,
                opt$`n-shards`, nrow(regions), ancestry))

# --- helpers -----------------------------------------------------------------

# slice one outcome's summary stats to the region; standard colocboost
# sumstat columns: variant, beta, sebeta, n
slice_outcome <- function(oc, chrom, start, end, n_xqtl) {
  dt <- tabix_slice(oc$file, chrom, start, end)
  if (is.null(dt)) return(NULL)
  if (oc$type == "xqtl") {
    dt <- dt[phenotype_id == oc$phenotype_id]
    if (nrow(dt) == 0L) return(NULL)
    out <- data.table(variant = dt$variant_id, beta = dt$slope,
                      sebeta = dt$slope_se, n = n_xqtl)
  } else {
    out <- data.table(variant = dt$var_id, beta = dt$beta,
                      sebeta = dt$se, n = as.numeric(dt$n))
  }
  out <- out[is.finite(beta) & is.finite(sebeta) & sebeta > 0]
  out <- out[!duplicated(variant)]
  if (nrow(out) == 0L) return(NULL)
  out
}

# --- per-region worker -------------------------------------------------------

run_region <- function(reg, reg_outcomes, tmp_dir) {
  t0 <- Sys.time()
  res_dir <- file.path(opt$outdir, "results", ancestry)
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
  clu_path <- file.path(res_dir, paste0(reg$region_id, ".clusters.tsv"))
  vcp_path <- file.path(res_dir, paste0(reg$region_id, ".vcp.tsv"))
  done_path <- file.path(res_dir, paste0(reg$region_id, ".done"))

  diag <- data.table(region_id = reg$region_id, ancestry = ancestry,
                     status = "error", message = "", n_outcomes = NA_integer_,
                     n_variants_xqtl = NA_integer_, n_variants_gwas = NA_integer_,
                     n_clusters = NA_integer_, n_nonconverged = NA_integer_,
                     walltime_sec = NA_real_)
  finish <- function(status, message = "") {
    diag$status <- status
    diag$message <- gsub("[\t\n\r]", " ", substr(message, 1, 300))
    diag$walltime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    terminal <- c("ok", "too_few_variants", "too_few_outcomes",
                  "dosage_failed_xqtl", "dosage_failed_gwas",
                  "colocboost_fail")
    if (status %in% terminal) writeLines("ok", done_path)
    diag
  }
  if (file.exists(done_path)) return(finish("ok", "already done"))

  n_xqtl <- psam_n(opt$`ld-xqtl-pgen`)
  if (is.na(n_xqtl) || n_xqtl < 10)
    return(finish("error", "could not read xQTL N from psam"))

  # 1. slice all outcomes -----------------------------------------------------
  sumstats <- vector("list", nrow(reg_outcomes))
  names(sumstats) <- reg_outcomes$outcome_name
  for (i in seq_len(nrow(reg_outcomes))) {
    oc <- reg_outcomes[i, ]
    sumstats[[i]] <- tryCatch(
      slice_outcome(oc, reg$chrom, reg$start, reg$end, n_xqtl),
      error = function(e) NULL)
  }
  ok <- !vapply(sumstats, is.null, logical(1))
  sumstats <- sumstats[ok]
  reg_outcomes <- reg_outcomes[ok, ]
  diag$n_outcomes <- length(sumstats)
  if (length(sumstats) < opt$`min-outcomes`)
    return(finish("too_few_outcomes"))
  if (!any(reg_outcomes$type == "gwas"))
    return(finish("too_few_outcomes", "no GWAS outcome with variants"))

  # 2. reference dosages per side --------------------------------------------
  vars_xqtl <- unique(unlist(lapply(sumstats[reg_outcomes$type == "xqtl"],
                                    `[[`, "variant")))
  vars_gwas <- unique(unlist(lapply(sumstats[reg_outcomes$type == "gwas"],
                                    `[[`, "variant")))
  Xx <- export_dosages(opt$plink2, opt$`ld-xqtl-pgen`, NULL, vars_xqtl,
                       file.path(tmp_dir, "xref_xqtl"))
  if (is.null(Xx)) return(finish("dosage_failed_xqtl"))
  Xg <- export_dosages(opt$plink2, opt$`ld-gwas-pgen`, opt$`ld-gwas-keep`, vars_gwas,
                       file.path(tmp_dir, "xref_gwas"))
  if (is.null(Xg)) return(finish("dosage_failed_gwas"))
  diag$n_variants_xqtl <- ncol(Xx)
  diag$n_variants_gwas <- ncol(Xg)
  if (ncol(Xx) < opt$`min-variants` || ncol(Xg) < opt$`min-variants`)
    return(finish("too_few_variants"))

  # restrict each outcome to variants present in its side's panel
  for (i in seq_along(sumstats)) {
    ref <- if (reg_outcomes$type[i] == "xqtl") colnames(Xx) else colnames(Xg)
    sumstats[[i]] <- sumstats[[i]][variant %in% ref]
  }
  ok <- vapply(sumstats, nrow, integer(1)) >= 10
  sumstats <- sumstats[ok]
  reg_outcomes <- reg_outcomes[ok, ]
  if (length(sumstats) < opt$`min-outcomes` ||
      !any(reg_outcomes$type == "gwas"))
    return(finish("too_few_outcomes", "after panel restriction"))

  # 3. colocBoost -------------------------------------------------------------
  dict <- cbind(seq_along(sumstats),
                ifelse(reg_outcomes$type == "xqtl", 1L, 2L))
  n_nonconv <- 0L
  res <- tryCatch(
    withCallingHandlers(
      colocboost(sumstat = unname(sumstats),
                 X_ref = list(Xx, Xg),
                 dict_sumstatLD = dict,
                 outcome_names = names(sumstats),
                 M = opt$M),
      warning = function(w) {
        if (grepl("did not coverage|did not converge", conditionMessage(w))) {
          n_nonconv <<- n_nonconv + 1L
          invokeRestart("muffleWarning")
        }
      }),
    error = function(e) e)
  diag$n_nonconverged <- n_nonconv
  if (inherits(res, "error"))
    return(finish("colocboost_fail", conditionMessage(res)))

  # 4. outputs ----------------------------------------------------------------
  summ <- tryCatch(as.data.table(get_colocboost_summary(res)$cos_summary),
                   error = function(e) NULL)
  n_clusters <- if (is.null(summ)) 0L else nrow(summ)
  diag$n_clusters <- n_clusters
  meta <- data.table(region_id = reg$region_id, ancestry = ancestry,
                     chrom = reg$chrom, start = reg$start, end = reg$end,
                     n_outcomes = length(sumstats),
                     n_xqtl_outcomes = sum(reg_outcomes$type == "xqtl"),
                     n_gwas_traits = sum(reg_outcomes$type == "gwas"),
                     n_variants_xqtl = ncol(Xx), n_variants_gwas = ncol(Xg),
                     n_nonconverged = n_nonconv)
  if (n_clusters > 0L)
    fwrite(cbind(meta[rep(1L, nrow(summ))], summ), clu_path, sep = "\t")
  else
    fwrite(cbind(meta, data.table(cos_id = NA_character_)), clu_path, sep = "\t")

  vcp <- res$vcp
  if (!is.null(vcp)) {
    vdt <- data.table(variant = names(vcp), vcp = as.numeric(vcp))
    vdt[, `:=`(region_id = reg$region_id, ancestry = ancestry)]
    fwrite(vdt, vcp_path, sep = "\t")
  }
  finish("ok")
}

# --- run shard ---------------------------------------------------------------

tmp_dir <- tempfile(pattern = "colocboost_", tmpdir = opt$`tmp-dir`)
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

diags <- vector("list", nrow(regions))
for (i in seq_len(nrow(regions))) {
  reg <- regions[i, ]
  reg_outcomes <- outcomes[region_id == reg$region_id]
  message(sprintf("[%d/%d] %s  %s:%d-%d  (%d outcomes)", i, nrow(regions),
                  reg$region_id, reg$chrom, reg$start, reg$end,
                  nrow(reg_outcomes)))
  diags[[i]] <- tryCatch(run_region(reg, reg_outcomes, tmp_dir),
                         error = function(e) data.table(
                           region_id = reg$region_id, ancestry = ancestry,
                           status = "error",
                           message = gsub("[\t\n\r]", " ",
                                          substr(conditionMessage(e), 1, 300)),
                           n_outcomes = NA_integer_, n_variants_xqtl = NA_integer_,
                           n_variants_gwas = NA_integer_, n_clusters = NA_integer_,
                           n_nonconverged = NA_integer_, walltime_sec = NA_real_))
}

diag_all <- rbindlist(diags, fill = TRUE)
fwrite(diag_all, diag_path, sep = "\t")
message("status counts:")
print(diag_all[, .N, by = status][order(-N)])
message("diagnostics -> ", diag_path)
