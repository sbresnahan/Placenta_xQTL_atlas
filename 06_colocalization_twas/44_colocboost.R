#!/usr/bin/env Rscript
# =============================================================================
# 44_colocboost.R — pure-R colocBoost statistical worker
# =============================================================================
# External I/O (tabix slicing and plink2 dosage export) is performed before R
# by 44_prepare_colocboost_inputs.py on the host compute node. This script only
# reads prepared files, fits colocBoost, and writes result/diagnostic files.
# It intentionally contains no system(), system2(), pipe(), or fread(cmd=...).
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
  make_option("--prepared-manifest", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--min-variants", type = "integer", default = 50L),
  make_option("--min-outcomes", type = "integer", default = 2L),
  make_option("--M", type = "integer", default = 500L,
              help = "max boosting rounds per outcome [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!file.exists(opt$`prepared-manifest`))
  stop("prepared manifest not found: ", opt$`prepared-manifest`)
prepared <- fread(opt$`prepared-manifest`)
if (nrow(prepared) == 0L) stop("prepared manifest has no rows")

ancestries <- unique(prepared$ancestry)
if (length(ancestries) != 1L) stop("prepared manifest mixes ancestries")
ancestry <- ancestries[1]
dir.create(file.path(opt$outdir, "diagnostics"), showWarnings = FALSE, recursive = TRUE)

# The host preparer writes the PLINK .raw file. R only parses it.
read_dosage_raw <- function(path, var_ids) {
  if (!file.exists(path)) return(NULL)
  raw <- fread(path)
  if (!("IID" %in% names(raw))) return(NULL)
  sample_ids <- as.character(raw[["IID"]])
  meta_cols <- unique(c("#FID", "FID", "IID", "PAT", "MAT", "SEX",
                        "PHENOTYPE", grep("^PHENO", names(raw), value = TRUE)))
  drop <- intersect(names(raw), meta_cols)
  dat <- as.data.frame(raw)[, setdiff(names(raw), drop), drop = FALSE]
  ids <- names(dat)
  # PLINK can append the counted allele to a dosage column. Only remove that
  # suffix when doing so recovers an explicitly requested canonical ID.
  ids_fix <- ids
  no_suffix <- sub("_[^_]+$", "", ids)
  need_fix <- !(ids %in% var_ids) & (no_suffix %in% var_ids)
  ids_fix[need_fix] <- no_suffix[need_fix]
  names(dat) <- ids_fix
  keep <- intersect(var_ids, names(dat))
  if (length(keep) == 0L) return(NULL)
  dat <- dat[, keep, drop = FALSE]
  mat <- as.matrix(dat)
  storage.mode(mat) <- "double"
  for (j in seq_len(ncol(mat))) {
    miss <- is.na(mat[, j])
    if (any(miss)) {
      mu <- mean(mat[, j], na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0
      mat[miss, j] <- mu
    }
  }
  sds <- apply(mat, 2, stats::sd)
  mat <- mat[, is.finite(sds) & sds > 0, drop = FALSE]
  if (ncol(mat) == 0L) return(NULL)
  rownames(mat) <- sample_ids
  mat
}

run_region <- function(reg) {
  t0 <- Sys.time()
  res_dir <- file.path(opt$outdir, "results", reg$ancestry)
  dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
  clu_path <- file.path(res_dir, paste0(reg$region_id, ".clusters.tsv"))
  vcp_path <- file.path(res_dir, paste0(reg$region_id, ".vcp.tsv"))
  done_path <- file.path(res_dir, paste0(reg$region_id, ".done"))
  diag <- data.table(region_id = reg$region_id, ancestry = reg$ancestry,
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

  prep_status <- if (is.na(reg$prep_status)) "error" else reg$prep_status
  prep_message <- if (is.na(reg$prep_message)) "" else reg$prep_message
  if (file.exists(done_path) || prep_status == "already_done")
    return(finish("ok", "already done"))
  if (prep_status != "ok")
    return(finish(prep_status, prep_message))

  sumdt <- tryCatch(fread(reg$sumstats_file), error = function(e) NULL)
  if (is.null(sumdt) || nrow(sumdt) == 0L)
    return(finish("error", "prepared sumstats unreadable or empty"))
  required <- c("outcome_name", "type", "variant", "beta", "sebeta", "n")
  if (!all(required %in% names(sumdt)))
    return(finish("error", "prepared sumstats missing required columns"))

  outcome_names <- unique(sumdt$outcome_name)
  sumstats <- lapply(outcome_names, function(nm) {
    sumdt[outcome_name == nm, .(variant, beta, sebeta, n)]
  })
  names(sumstats) <- outcome_names
  outcome_type <- vapply(outcome_names, function(nm) {
    unique(sumdt[outcome_name == nm, type])[1]
  }, character(1))

  diag$n_outcomes <- length(sumstats)
  if (length(sumstats) < opt$`min-outcomes`)
    return(finish("too_few_outcomes"))
  if (!any(outcome_type == "gwas"))
    return(finish("too_few_outcomes", "no GWAS outcome with variants"))

  vars_x <- unique(sumdt[type == "xqtl", variant])
  vars_g <- unique(sumdt[type == "gwas", variant])
  Xx <- read_dosage_raw(reg$xqtl_raw, vars_x)
  if (is.null(Xx)) return(finish("dosage_failed_xqtl"))
  Xg <- read_dosage_raw(reg$gwas_raw, vars_g)
  if (is.null(Xg)) return(finish("dosage_failed_gwas"))
  diag$n_variants_xqtl <- ncol(Xx)
  diag$n_variants_gwas <- ncol(Xg)
  if (ncol(Xx) < opt$`min-variants` || ncol(Xg) < opt$`min-variants`)
    return(finish("too_few_variants"))

  for (i in seq_along(sumstats)) {
    ref <- if (outcome_type[i] == "xqtl") colnames(Xx) else colnames(Xg)
    sumstats[[i]] <- sumstats[[i]][variant %in% ref]
  }
  ok <- vapply(sumstats, nrow, integer(1)) >= 10L
  sumstats <- sumstats[ok]
  outcome_type <- outcome_type[ok]
  if (length(sumstats) < opt$`min-outcomes` || !any(outcome_type == "gwas"))
    return(finish("too_few_outcomes", "after panel restriction"))

  dict <- cbind(seq_along(sumstats), ifelse(outcome_type == "xqtl", 1L, 2L))
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

  summ <- tryCatch(as.data.table(get_colocboost_summary(res)$cos_summary),
                   error = function(e) NULL)
  n_clusters <- if (is.null(summ)) 0L else nrow(summ)
  diag$n_clusters <- n_clusters
  meta <- data.table(region_id = reg$region_id, ancestry = reg$ancestry,
                     chrom = reg$chrom, start = reg$start, end = reg$end,
                     n_outcomes = length(sumstats),
                     n_xqtl_outcomes = sum(outcome_type == "xqtl"),
                     n_gwas_traits = sum(outcome_type == "gwas"),
                     n_variants_xqtl = ncol(Xx), n_variants_gwas = ncol(Xg),
                     n_nonconverged = n_nonconv)
  if (n_clusters > 0L)
    fwrite(cbind(meta[rep(1L, nrow(summ))], summ), clu_path, sep = "\t")
  else
    fwrite(cbind(meta, data.table(cos_id = NA_character_)), clu_path, sep = "\t")
  vcp <- res$vcp
  if (!is.null(vcp)) {
    vdt <- data.table(variant = names(vcp), vcp = as.numeric(vcp))
    vdt[, `:=`(region_id = reg$region_id, ancestry = reg$ancestry)]
    fwrite(vdt, vcp_path, sep = "\t")
  }
  finish("ok")
}

diags <- vector("list", nrow(prepared))
for (i in seq_len(nrow(prepared))) {
  reg <- prepared[i]
  message(sprintf("[%d/%d] %s  %s:%s-%s", i, nrow(prepared), reg$region_id,
                  reg$chrom, reg$start, reg$end))
  diags[[i]] <- tryCatch(run_region(reg), error = function(e) data.table(
    region_id = reg$region_id, ancestry = reg$ancestry, status = "error",
    message = gsub("[\t\n\r]", " ", substr(conditionMessage(e), 1, 300)),
    n_outcomes = NA_integer_, n_variants_xqtl = NA_integer_,
    n_variants_gwas = NA_integer_, n_clusters = NA_integer_,
    n_nonconverged = NA_integer_, walltime_sec = NA_real_))
}

diag_all <- rbindlist(diags, fill = TRUE)
manifest_stem <- tools::file_path_sans_ext(basename(opt$`prepared-manifest`))
diag_path <- file.path(opt$outdir, "diagnostics",
                       paste0(ancestry, ".", manifest_stem, ".diagnostics.tsv"))
fwrite(diag_all, diag_path, sep = "\t")
message("status counts:")
print(diag_all[, .N, by = status][order(-N)])
message("diagnostics -> ", diag_path)
