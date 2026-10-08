#!/usr/bin/env Rscript
# =============================================================================
# 41_susie_coloc.R — SuSiE-coloc worker for one shard of coloc tasks
# =============================================================================
# Objective 1.6: pairwise colocalization between placental xQTL signals and
# GWAS of birth weight / gestational duration / glycemic traits / childhood
# adiposity, using ancestry-matched LD on both sides.
#
# Each row of the task list (40_prepare_coloc_loci.py ->
# {COLOC_DIR}/loci/{MOD}.tasks.tsv) is one
#   modality x phenotype_id (fine-mapped locus) x ancestry x GWAS trait
# combination. For every task this worker:
#   1. slices the merged genome-wide nominal xQTL stats (39_run_nominal.py)
#      and the harmonized GWAS stats (38_harmonize_gwas.py) to the locus
#      window with tabix and merges them on var_id (chr:pos:ref:alt);
#   2. computes in-sample LD for the xQTL side (plink2 --r-unphased square on
#      the intersected {ANC}_qtl pgen) and reference LD for the GWAS side
#      (1KG pgen restricted to the matching superpopulation via --keep);
#   3. drops variants missing from either LD source and requires
#      >= --min-variants shared variants;
#   4. fits coloc::runsusie() on each dataset (xQTL L from the fine-mapping
#      independent-signal count, GWAS L = --gwas-L) and runs
#      coloc::coloc.susie();
#   5. writes per-task outputs and a .done sentinel.
#
# Convergence failures are FLAGGED (status column), not rescued — the aims
# eCAVIAR fallback is intentionally not implemented.
#
# Outputs (per task, under {outdir}/results/{modality}/):
#   {anc}_{phenotype}_{trait}.coloc.tsv     one row per credible-set pair
#                                           (PP.H0..PP.H4) + task metadata
#   {anc}_{phenotype}_{trait}.variants.tsv  merged per-variant stats +
#                                           SNP.PP.H4 from the best pair
#   {anc}_{phenotype}_{trait}.done          sentinel (touch)
# Per shard:
#   {outdir}/diagnostics/{modality}.shard-{idx}.diagnostics.tsv
#
# Usage:
#   Rscript 41_susie_coloc.R --tasks {MOD}.tasks.tsv --shard-index 1 \
#       --n-shards 40 --outdir $RESULTS_DIR/coloc
# =============================================================================

R_LIB <- "/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"
if (!dir.exists(R_LIB)) stop("R package library does not exist: ", R_LIB)
.libPaths(c(R_LIB, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(coloc)
})

option_list <- list(
  make_option("--tasks", type = "character", help = "task list TSV from 40_prepare_coloc_loci.py"),
  make_option("--shard-index", type = "integer", default = 1L,
              help = "1-based shard index (e.g. $LSB_JOBINDEX)"),
  make_option("--n-shards", type = "integer", default = 1L,
              help = "total number of shards for this task list"),
  make_option("--outdir", type = "character", help = "coloc output dir (COLOC_DIR)"),
  make_option("--plink2", type = "character", default = "plink2",
              help = "plink2 binary [default %default]"),
  make_option("--min-variants", type = "integer", default = 50L,
              help = "minimum shared variants with LD on both sides [default %default]"),
  make_option("--gwas-L", type = "integer", default = 10L,
              help = "max credible sets for the GWAS-side SuSiE fit [default %default]"),
  make_option("--max-L", type = "integer", default = 20L,
              help = "cap on the xQTL-side L taken from the locus list [default %default]"),
  make_option("--pp-h4", type = "double", default = 0.7,
              help = "PP.H4 threshold for the colocalization call [default %default]"),
  make_option("--tmp-dir", type = "character", default = tempdir(),
              help = "scratch dir for plink2 LD output [default tempdir()]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# shared helpers (tabix_slice, run_ld, psam_n, sanitize_id)
.args_all <- commandArgs(trailingOnly = FALSE)
.file_arg <- sub("^--file=", "", grep("^--file=", .args_all, value = TRUE))
.scripts_dir <- Sys.getenv("SCRIPTS_DIR",
                           unset = if (length(.file_arg)) dirname(.file_arg[1]) else getwd())
source(file.path(.scripts_dir, "coloc_common.R"))

dir.create(file.path(opt$outdir, "diagnostics"), showWarnings = FALSE, recursive = TRUE)

tasks <- fread(opt$tasks)
modality <- unique(tasks$modality)
if (length(modality) != 1L) stop("task file mixes modalities: ", paste(modality, collapse = ","))
diag_path <- file.path(opt$outdir, "diagnostics",
                       sprintf("%s.shard-%04d.diagnostics.tsv", modality, opt$`shard-index`))

# round-robin shard assignment: row i -> shard ((i-1) %% n_shards) + 1
shard_rows <- which(((seq_len(nrow(tasks)) - 1L) %% opt$`n-shards`) + 1L == opt$`shard-index`)
tasks <- tasks[shard_rows, ]
message(sprintf("[shard %d/%d] %d tasks from %s", opt$`shard-index`, opt$`n-shards`,
                nrow(tasks), opt$tasks))

sanitize <- sanitize_id
get_n_xqtl <- psam_n

# --- per-task worker ---------------------------------------------------------

run_task <- function(task, tmp_dir) {
  t0 <- Sys.time()
  base <- sprintf("%s_%s_%s", task$ancestry, sanitize(task$phenotype_id), task$trait_id)
  res_dir <- file.path(opt$outdir, "results", task$modality)
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
  coloc_path <- file.path(res_dir, paste0(base, ".coloc.tsv"))
  var_path   <- file.path(res_dir, paste0(base, ".variants.tsv"))
  done_path  <- file.path(res_dir, paste0(base, ".done"))

  diag <- data.table(modality = task$modality, phenotype_id = task$phenotype_id,
                     trait_id = task$trait_id, ancestry = task$ancestry,
                     status = "error", message = "", nsnps_merged = NA_integer_,
                     nsnps_final = NA_integer_, pp_h4_max = NA_real_,
                     converged_xqtl = NA, converged_gwas = NA,
                     walltime_sec = NA_real_)

  finish <- function(status, message = "") {
    diag$status <- status
    diag$message <- gsub("[\t\n\r]", " ", substr(message, 1, 300))
    diag$walltime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    # .done for terminal states (incl. flagged convergence failures); transient
    # I/O failures (ld_failed_*, unexpected errors) are retried on rerun
    terminal <- c("ok", "no_cs_xqtl", "no_cs_gwas", "no_cs_pair",
                  "susie_fail_xqtl", "susie_fail_gwas", "coloc_fail",
                  "too_few_variants", "no_xqtl_variants", "no_gwas_variants",
                  "no_shared_variants")
    if (status %in% terminal) writeLines("ok", done_path)
    diag
  }

  if (file.exists(done_path)) return(finish("ok", "already done"))

  # 1. slice summary statistics ---------------------------------------------
  xq <- tabix_slice(task$xqtl_file, task$chrom, task$start, task$end)
  if (is.null(xq)) return(finish("no_xqtl_variants"))
  xq <- xq[phenotype_id == task$phenotype_id]
  if (nrow(xq) == 0L) return(finish("no_xqtl_variants", "phenotype absent from slice"))

  gw <- tabix_slice(task$gwas_file, task$chrom, task$start, task$end)
  if (is.null(gw)) return(finish("no_gwas_variants"))

  # 2. merge on var_id --------------------------------------------------------
  xq <- xq[, .(var_id = variant_id, pos_x = pos, slope, slope_se, af,
               pval_nominal)]
  gw <- gw[, .(var_id, pos_g = pos, beta, se, eaf, pval, n)]
  m <- merge(xq, gw, by = "var_id")
  if (nrow(m) == 0L) return(finish("no_shared_variants"))
  # dedupe (paranoid; ids should be unique) keeping the min xQTL p
  setorder(m, pval_nominal)
  m <- m[!duplicated(var_id)]
  # usable rows only
  m <- m[is.finite(slope) & is.finite(slope_se) & slope_se > 0 &
         is.finite(beta) & is.finite(se) & se > 0 &
         is.finite(af) & af > 0 & af < 1 &
         is.finite(eaf) & eaf > 0 & eaf < 1]
  diag$nsnps_merged <- nrow(m)
  if (nrow(m) < opt$`min-variants`) return(finish("too_few_variants", "after merge"))
  setorder(m, pos_x)
  var_ids <- m$var_id

  # 3. LD on both sides -------------------------------------------------------
  ld_x <- run_ld(opt$plink2, task$ld_xqtl_pgen, NULL, var_ids, file.path(tmp_dir, "ld_xqtl"))
  if (is.null(ld_x)) return(finish("ld_failed_xqtl"))
  ld_g <- run_ld(opt$plink2, task$ld_gwas_pgen, task$ld_gwas_keep, var_ids,
                 file.path(tmp_dir, "ld_gwas"))
  if (is.null(ld_g)) return(finish("ld_failed_gwas"))

  keep <- intersect(var_ids, intersect(ld_x$present, ld_g$present))
  m <- m[var_id %in% keep]
  # drop variants monomorphic (NA diagonal) in either panel
  dx <- diag(ld_x$LD)[match(m$var_id, dimnames(ld_x$LD)[[1]])]
  dg <- diag(ld_g$LD)[match(m$var_id, dimnames(ld_g$LD)[[1]])]
  m <- m[is.finite(dx) & is.finite(dg)]
  diag$nsnps_final <- nrow(m)
  if (nrow(m) < opt$`min-variants`) return(finish("too_few_variants", "after LD filter"))

  LDx <- ld_x$LD[m$var_id, m$var_id, drop = FALSE]
  LDg <- ld_g$LD[m$var_id, m$var_id, drop = FALSE]

  # 4. coloc datasets ---------------------------------------------------------
  n_xqtl <- get_n_xqtl(task$ld_xqtl_pgen)
  if (is.na(n_xqtl) || n_xqtl < 10) return(finish("error", "could not read xQTL N from psam"))
  n_gwas <- as.numeric(stats::median(m$n, na.rm = TRUE))
  if (!is.finite(n_gwas)) n_gwas <- NA_real_

  d1 <- list(beta = m$slope, varbeta = m$slope_se^2, snp = m$var_id,
             position = m$pos_x, type = "quant", N = n_xqtl, MAF = m$af,
             LD = LDx)
  trait_type <- if (is.na(task$trait_type) || task$trait_type == "") "quantitative" else task$trait_type
  d2 <- list(beta = m$beta, varbeta = m$se^2, snp = m$var_id,
             position = m$pos_x,
             type = if (trait_type == "cc") "cc" else "quant",
             N = n_gwas, MAF = m$eaf, LD = LDg)
  if (trait_type == "cc") {
    s <- suppressWarnings(as.numeric(task$prop_cases))
    if (is.na(s) || s <= 0 || s >= 1)
      return(finish("error", "cc trait without valid prop_cases"))
    d2$s <- s
  }

  # 5. SuSiE fits + coloc -----------------------------------------------------
  L_x <- min(max(as.integer(task$L), 5L), opt$`max-L`)
  fit1 <- tryCatch(
    runsusie(d1, suffix = 1, repeat_until_convergence = FALSE, L = L_x),
    error = function(e) e)
  if (inherits(fit1, "error"))
    return(finish("susie_fail_xqtl", conditionMessage(fit1)))
  diag$converged_xqtl <- isTRUE(fit1$converged)
  if (length(fit1$sets$cs) == 0L) return(finish("no_cs_xqtl"))

  fit2 <- tryCatch(
    runsusie(d2, suffix = 2, repeat_until_convergence = FALSE, L = opt$`gwas-L`),
    error = function(e) e)
  if (inherits(fit2, "error"))
    return(finish("susie_fail_gwas", conditionMessage(fit2)))
  diag$converged_gwas <- isTRUE(fit2$converged)
  if (length(fit2$sets$cs) == 0L) return(finish("no_cs_gwas"))

  res <- tryCatch(coloc.susie(fit1, fit2), error = function(e) e)
  if (inherits(res, "error"))
    return(finish("coloc_fail", conditionMessage(res)))
  # coloc.susie returns a bare data.table(nsnps = NA) when no credible-set
  # pair can be tested
  if (is.null(res$summary) || !("PP.H4.abf" %in% names(res$summary)))
    return(finish("no_cs_pair"))

  # 6. outputs ----------------------------------------------------------------
  summ <- as.data.table(res$summary)
  meta <- data.table(modality = task$modality, phenotype_id = task$phenotype_id,
                     trait_id = task$trait_id, ancestry = task$ancestry,
                     chrom = task$chrom, start = task$start, end = task$end,
                     nsnps = nrow(m), n_xqtl = n_xqtl, n_gwas = n_gwas,
                     converged_xqtl = diag$converged_xqtl,
                     converged_gwas = diag$converged_gwas,
                     n_cs_xqtl = length(fit1$sets$cs),
                     n_cs_gwas = length(fit2$sets$cs))
  summ <- cbind(meta[rep(1L, nrow(summ))], summ)
  pp_col <- grep("PP.H4", names(summ), value = TRUE)[1]
  summ[, coloc_call := get(pp_col) >= opt$`pp-h4`]
  diag$pp_h4_max <- max(summ[[pp_col]], na.rm = TRUE)
  fwrite(summ, coloc_path, sep = "\t")

  best <- summ[which.max(get(pp_col))]
  # res$results: data.table(snp, SNP.PP.H4.abf) for a single credible-set
  # pair, or data.table(snp, SNP.PP.H4.row1..rowK) for K pairs (row k matches
  # summary row k)
  per_var <- as.data.table(res$results)
  setnames(per_var, "snp", "var_id")
  pp_var_cols <- grep("^SNP\\.PP\\.H4", names(per_var), value = TRUE)
  best_row <- which.max(summ[[pp_col]])
  best_col <- if ("SNP.PP.H4.abf" %in% pp_var_cols) "SNP.PP.H4.abf" else
    sprintf("SNP.PP.H4.row%d", best_row)
  if (best_col %in% names(per_var))
    setnames(per_var, best_col, "SNP.PP.H4_best")
  per_var <- merge(m[, .(var_id, chrom = task$chrom, pos = pos_x, slope,
                         slope_se, pval_nominal, beta, se, pval, af, eaf)],
                   per_var, by = "var_id", all.x = TRUE)
  per_var[, `:=`(modality = task$modality, phenotype_id = task$phenotype_id,
                 trait_id = task$trait_id, ancestry = task$ancestry,
                 best_hit1 = best$hit1, best_hit2 = best$hit2,
                 pp_h4_best = best[[pp_col]])]
  fwrite(per_var, var_path, sep = "\t")

  finish("ok")
}

# --- run shard ---------------------------------------------------------------

tmp_dir <- tempfile(pattern = "coloc_", tmpdir = opt$`tmp-dir`)
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

diags <- vector("list", nrow(tasks))
for (i in seq_len(nrow(tasks))) {
  task <- tasks[i, ]
  message(sprintf("[%d/%d] %s | %s | %s | %s:%d-%d", i, nrow(tasks),
                  task$ancestry, task$phenotype_id, task$trait_id,
                  task$chrom, task$start, task$end))
  diags[[i]] <- tryCatch(run_task(task, tmp_dir),
                         error = function(e) data.table(
                           modality = task$modality, phenotype_id = task$phenotype_id,
                           trait_id = task$trait_id, ancestry = task$ancestry,
                           status = "error",
                           message = gsub("[\t\n\r]", " ", substr(conditionMessage(e), 1, 300)),
                           nsnps_merged = NA_integer_, nsnps_final = NA_integer_,
                           pp_h4_max = NA_real_, converged_xqtl = NA, converged_gwas = NA,
                           walltime_sec = NA_real_))
}

diag_all <- rbindlist(diags, fill = TRUE)
fwrite(diag_all, diag_path, sep = "\t")
message("status counts:")
print(diag_all[, .N, by = status][order(-N)])
message("diagnostics -> ", diag_path)
