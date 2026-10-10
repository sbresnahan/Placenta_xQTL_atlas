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
#   1. reads summary-statistics and LD artifacts prepared by
#      41_prepare_susie_coloc_inputs.py in the host LSF environment;
#   2. drops variants missing from either LD source and requires
#      >= --min-variants shared variants;
#   3. fits coloc::runsusie() on each dataset (xQTL L from the fine-mapping
#      independent-signal count, GWAS L = --gwas-L) and runs
#      coloc::coloc.susie();
#   4. writes per-task outputs and a .done sentinel.
#
# This R worker intentionally never invokes command-line tools.
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
#   Rscript 41_susie_coloc.R --prepared-tasks prepared.tasks.tsv \
#       --shard-index 1 --n-shards 40 --outdir $RESULTS_DIR/coloc
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
  make_option("--prepared-tasks", type = "character",
              help = "prepared shard TSV from 41_prepare_susie_coloc_inputs.py"),
  make_option("--shard-index", type = "integer", default = 1L,
              help = "1-based shard index (e.g. $LSB_JOBINDEX)"),
  make_option("--n-shards", type = "integer", default = 1L,
              help = "total number of shards for this task list"),
  make_option("--outdir", type = "character", help = "coloc output dir (COLOC_DIR)"),
  make_option("--min-variants", type = "integer", default = 50L,
              help = "minimum shared variants with LD on both sides [default %default]"),
  make_option("--gwas-L", type = "integer", default = 10L,
              help = "max credible sets for the GWAS-side SuSiE fit [default %default]"),
  make_option("--max-L", type = "integer", default = 20L,
              help = "cap on the xQTL-side L taken from the locus list [default %default]"),
  make_option("--pp-h4", type = "double", default = 0.7,
              help = "PP.H4 threshold for the colocalization call [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

sanitize_id <- function(x) gsub("[^A-Za-z0-9._+-]", "_", x)
read_ld <- function(matrix_path, vars_path) {
  if (is.na(matrix_path) || is.na(vars_path) ||
      !file.exists(matrix_path) || !file.exists(vars_path)) return(NULL)
  ids <- readLines(vars_path, warn = FALSE)
  ids <- ids[!grepl("^#", ids)]
  if (length(ids) == 0L) return(NULL)
  mat <- as.matrix(data.table::fread(matrix_path, header = FALSE,
                                     showProgress = FALSE))
  storage.mode(mat) <- "double"
  if (nrow(mat) != length(ids) || ncol(mat) != length(ids)) return(NULL)
  dimnames(mat) <- list(ids, ids)
  list(LD = mat, present = ids)
}

dir.create(file.path(opt$outdir, "diagnostics"), showWarnings = FALSE, recursive = TRUE)

tasks <- fread(opt$`prepared-tasks`)
if (nrow(tasks) == 0L) stop("prepared task manifest has no rows: ", opt$`prepared-tasks`)
modality <- unique(tasks$modality)
if (length(modality) != 1L) stop("task file mixes modalities: ", paste(modality, collapse = ","))
diag_path <- file.path(opt$outdir, "diagnostics",
                       sprintf("%s.shard-%04d.diagnostics.tsv", modality, opt$`shard-index`))
message(sprintf("[shard %d/%d] %d prepared tasks from %s", opt$`shard-index`,
                opt$`n-shards`, nrow(tasks), opt$`prepared-tasks`))

sanitize <- sanitize_id

# --- per-task worker ---------------------------------------------------------

run_task <- function(task) {
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
    if (length(message) == 0L || is.na(message)) message <- ""
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

  # External I/O (tabix/plink2) has already run in the host-side preparer.
  prep_status <- if (is.na(task$prep_status)) "error" else task$prep_status
  prep_message <- if (is.na(task$prep_message)) "" else task$prep_message
  if (prep_status != "ok") return(finish(prep_status, prep_message))

  m <- tryCatch(fread(task$merged_file), error = function(e) NULL)
  if (is.null(m)) return(finish("error", "could not read prepared merged stats"))
  diag$nsnps_merged <- nrow(m)
  if (nrow(m) < opt$`min-variants`) return(finish("too_few_variants", "after merge"))
  setorder(m, pos_x)
  var_ids <- m$var_id

  # 1. LD files prepared by plink2 outside R ---------------------------------
  ld_x <- read_ld(task$ld_x_matrix, task$ld_x_vars)
  if (is.null(ld_x)) return(finish("ld_failed_xqtl"))
  ld_g <- read_ld(task$ld_g_matrix, task$ld_g_vars)
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

  # 2. coloc datasets ---------------------------------------------------------
  n_xqtl <- suppressWarnings(as.numeric(task$n_xqtl))
  if (is.na(n_xqtl) || n_xqtl < 10) return(finish("error", "invalid prepared xQTL N"))
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

  # 3. SuSiE fits + coloc -----------------------------------------------------
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

  # 4. outputs ----------------------------------------------------------------
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

diags <- vector("list", nrow(tasks))
for (i in seq_len(nrow(tasks))) {
  task <- tasks[i, ]
  message(sprintf("[%d/%d] %s | %s | %s | %s:%d-%d", i, nrow(tasks),
                  task$ancestry, task$phenotype_id, task$trait_id,
                  task$chrom, task$start, task$end))
  diags[[i]] <- tryCatch(run_task(task),
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
