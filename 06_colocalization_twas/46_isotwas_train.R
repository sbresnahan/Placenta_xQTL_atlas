#!/usr/bin/env Rscript
# =============================================================================
# 46_isotwas_train.R — pure-R isoTWAS/TWAS model training worker
# =============================================================================
# 46_prepare_isotwas_inputs.py performs all host-side external I/O first:
# bgzip BED extraction is done with Python gzip and cis dosages are exported
# with plink2. R only reads prepared TSV/.raw files, residualizes phenotypes and
# genotypes, fits glmnet models, and writes FUSION weight files.
#
# This script intentionally contains no system(), system2(), pipe(), or
# fread(cmd=...) calls.
# =============================================================================
R_LIB <- "/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"
if (!dir.exists(R_LIB)) stop("R package library does not exist: ", R_LIB)
.libPaths(c(R_LIB, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(glmnet)
})

option_list <- list(
  make_option("--weight-set", type = "character", help = "EAS, EUR, or pooled"),
  make_option("--shard-index", type = "integer", default = 1L),
  make_option("--prepared-manifest", type = "character"),
  make_option("--qtl-dir", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--alpha", type = "double", default = 0.5),
  make_option("--r2-min", type = "double", default = 0.01),
  make_option("--nfolds", type = "integer", default = 5L),
  make_option("--min-variants", type = "integer", default = 10L),
  make_option("--max-isoforms", type = "integer", default = 50L,
              help = "cap isoforms per gene (most variable kept)")
)
opt <- parse_args(OptionParser(option_list = option_list))

ws <- opt$`weight-set`
ancestries <- if (ws == "pooled") c("EAS", "EUR") else ws
if (!all(ancestries %in% c("EAS", "EUR", "AFR", "AMR", "SAS")))
  stop("unknown weight set: ", ws)
if (!file.exists(opt$`prepared-manifest`))
  stop("prepared manifest not found: ", opt$`prepared-manifest`)
prepared <- fread(opt$`prepared-manifest`)

dir.create(file.path(opt$outdir, "diagnostics"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$outdir, "weights", ws, "genes"),
           showWarnings = FALSE, recursive = TRUE)
diag_path <- file.path(opt$outdir, "diagnostics",
                        sprintf("isotwas_%s.shard-%04d.diagnostics.tsv", ws,
                                opt$`shard-index`))
pos_path <- file.path(opt$outdir, "weights", ws,
                      sprintf("shard-%04d.pos", opt$`shard-index`))

read_covariates <- function(path) {
  dt <- fread(path)
  cn <- names(dt)
  id_col <- cn[1]
  samples <- cn[-1]
  mat <- t(as.matrix(dt[, -1, with = FALSE]))
  storage.mode(mat) <- "double"
  rownames(mat) <- samples
  colnames(mat) <- dt[[id_col]]
  mat
}

read_bed_subset <- function(path) {
  dt <- fread(path)
  if (ncol(dt) < 4L) stop("prepared BED subset has <4 columns: ", path)
  setnames(dt, 1:4, c("chr", "start", "end", "phenotype_id"))
  dt
}

# Read a host-prepared PLINK .raw file. Variant IDs in this project are
# canonical chrom:pos:ref:alt strings; PLINK may append _ALLELE to dosage
# column labels, which is removed only for that canonical four-field pattern.
read_dosage_raw <- function(path) {
  if (!file.exists(path)) return(NULL)
  raw <- fread(path)
  if (!("IID" %in% names(raw))) return(NULL)
  sample_ids <- as.character(raw[["IID"]])
  meta_cols <- unique(c("#FID", "FID", "IID", "PAT", "MAT", "SEX",
                        "PHENOTYPE", grep("^PHENO", names(raw), value = TRUE)))
  drop <- intersect(names(raw), meta_cols)
  dat <- as.data.frame(raw)[, setdiff(names(raw), drop), drop = FALSE]
  ids <- names(dat)
  canon <- grepl("^.+:[0-9]+:[^:]+:[^:]+_[^_]+$", ids)
  ids[canon] <- sub("^(.+:[0-9]+:[^:]+:[^:]+)_[^_]+$", "\\1",
                    ids[canon], perl = TRUE)
  names(dat) <- ids
  dat <- dat[, !duplicated(names(dat)), drop = FALSE]
  if (ncol(dat) == 0L) return(NULL)
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

message("loading covariates for: ", paste(ancestries, collapse = ", "))
covs <- list()
for (anc in ancestries) {
  covs[[anc]] <- list(
    expression = read_covariates(file.path(
      opt$`qtl-dir`, sprintf("%s_covariates_expression.tsv", anc))),
    isoform_expression = read_covariates(file.path(
      opt$`qtl-dir`, sprintf("%s_covariates_isoform_expression.tsv", anc))))
}

train_gene <- function(rec) {
  gene <- rec$gene
  t0 <- Sys.time()
  wgt_dir <- file.path(opt$outdir, "weights", ws, "genes")
  rdat_path <- file.path(wgt_dir, paste0(gene, ".wgt.RDat"))
  diag <- data.table(gene = gene, weight_set = ws, status = "error",
                     message = "", n_variants = NA_integer_,
                     n_isoforms_tested = NA_integer_, n_models_kept = NA_integer_,
                     best_rsq = NA_real_, walltime_sec = NA_real_)
  finish <- function(status, message = "") {
    diag$status <- status
    diag$message <- gsub("[\t\n\r]", " ", substr(message, 1, 300))
    diag$walltime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    diag
  }

  prep_status <- if (is.na(rec$prep_status)) "error" else rec$prep_status
  prep_message <- if (is.na(rec$prep_message)) "" else rec$prep_message
  if (file.exists(rdat_path) || prep_status == "already_done")
    return(finish("ok", "already done"))
  if (prep_status != "ok") return(finish(prep_status, prep_message))

  data_dir <- rec$data_dir
  X_list <- list()
  bed_e <- bed_i <- list()
  for (anc in ancestries) {
    X <- read_dosage_raw(file.path(data_dir, paste0(anc, ".raw")))
    if (is.null(X)) return(finish("dosage_failed", anc))
    X_list[[anc]] <- X
    bed_e[[anc]] <- read_bed_subset(file.path(data_dir,
                                               paste0(anc, ".expression.tsv")))
    bed_i[[anc]] <- read_bed_subset(file.path(data_dir,
                                               paste0(anc, ".isoform_expression.tsv")))
  }

  # The original worker intersected pvar IDs before exporting pooled dosages.
  # Here export happens first; intersecting prepared dosage columns is the same
  # operation after additionally removing panel-specific monomorphic variants.
  if (length(X_list) > 1L) {
    common_vars <- Reduce(intersect, lapply(X_list, colnames))
    if (length(common_vars) < opt$`min-variants`)
      return(finish("too_few_variants"))
    X_list <- lapply(X_list, function(x) x[, common_vars, drop = FALSE])
  }

  Xtr_list <- Yiso_list <- list()
  yexpr_list <- list()
  for (anc in ancestries) {
    X <- X_list[[anc]]
    erow <- bed_e[[anc]][phenotype_id == gene]
    if (nrow(erow) == 0L) next
    e_samples <- setdiff(names(bed_e[[anc]]), c("chr", "start", "end", "phenotype_id"))
    e_val <- as.numeric(erow[1, ..e_samples])
    names(e_val) <- e_samples

    irows <- bed_i[[anc]][startsWith(phenotype_id, paste0(gene, "__"))]
    iso_mat <- NULL
    if (nrow(irows) > 0L) {
      i_samples <- setdiff(names(bed_i[[anc]]), c("chr", "start", "end", "phenotype_id"))
      iso_mat <- t(as.matrix(irows[, ..i_samples]))
      storage.mode(iso_mat) <- "double"
      colnames(iso_mat) <- irows$phenotype_id
      rownames(iso_mat) <- i_samples
    }

    Ce <- covs[[anc]]$expression
    Ci <- covs[[anc]]$isoform_expression
    common_sets <- list(rownames(X), names(e_val), rownames(Ce), rownames(Ci))
    if (!is.null(iso_mat)) common_sets <- c(common_sets, list(rownames(iso_mat)))
    common <- Reduce(intersect, common_sets)
    if (length(common) < 20L) next

    X <- X[common, , drop = FALSE]
    e_val <- e_val[common]
    if (!is.null(iso_mat)) iso_mat <- iso_mat[common, , drop = FALSE]
    Ce <- Ce[common, , drop = FALSE]
    Ci <- Ci[common, , drop = FALSE]
    Cu <- cbind(Ce, Ci[, setdiff(colnames(Ci), colnames(Ce)), drop = FALSE])
    Xr <- qr.resid(qr(cbind(1, Cu)), X)
    y_e <- as.numeric(qr.resid(qr(cbind(1, Ce)), matrix(e_val, ncol = 1)))
    names(y_e) <- common
    s <- stats::sd(y_e)
    if (!is.finite(s) || s == 0) next
    y_e <- (y_e - mean(y_e)) / s

    Yi <- NULL
    if (!is.null(iso_mat)) {
      Yi <- qr.resid(qr(cbind(1, Ci)), iso_mat)
      keep <- apply(Yi, 2, function(v) is.finite(stats::sd(v)) && stats::sd(v) > 0)
      Yi <- Yi[, keep, drop = FALSE]
      if (ncol(Yi) > 0L) Yi <- scale(Yi) else Yi <- NULL
    }
    Xtr_list[[anc]] <- Xr
    yexpr_list[[anc]] <- y_e
    Yiso_list[[anc]] <- Yi
  }

  if (length(Xtr_list) == 0L) return(finish("error", "no usable ancestry block"))
  # Ensure any post-residualization ancestry blocks still share columns/order.
  common_vars <- Reduce(intersect, lapply(Xtr_list, colnames))
  if (length(common_vars) < opt$`min-variants`)
    return(finish("too_few_variants"))
  Xtr_list <- lapply(Xtr_list, function(x) x[, common_vars, drop = FALSE])
  Xtr <- do.call(rbind, Xtr_list)
  yexpr <- unlist(yexpr_list)

  Yiso <- NULL
  if (!any(vapply(Yiso_list, is.null, logical(1)))) {
    iso_ids <- Reduce(intersect, lapply(Yiso_list, colnames))
    if (length(iso_ids) > 0L) {
      Yiso <- do.call(rbind, lapply(Yiso_list,
                                    function(m) m[, iso_ids, drop = FALSE]))
      if (nrow(Yiso) != nrow(Xtr)) Yiso <- NULL
    }
  }

  if (ws == "pooled" && length(Xtr_list) > 1L) {
    anc_dummy <- rep(seq_along(Xtr_list), vapply(Xtr_list, nrow, integer(1)))
    D <- cbind(1, model.matrix(~ factor(anc_dummy) - 1))
    Xtr <- qr.resid(qr(D), Xtr)
    yexpr <- as.numeric(qr.resid(qr(D), as.matrix(yexpr)))
    if (!is.null(Yiso)) Yiso <- qr.resid(qr(D), Yiso)
  }

  diag$n_variants <- ncol(Xtr)
  if (ncol(Xtr) < opt$`min-variants`) return(finish("too_few_variants"))
  Xs <- scale(Xtr)
  Xs <- Xs[, apply(Xs, 2, function(v) all(is.finite(v))), drop = FALSE]
  if (ncol(Xs) < opt$`min-variants`) return(finish("too_few_variants"))

  fit_column <- function(y) {
    fit <- tryCatch(cv.glmnet(Xs, y, family = "gaussian", alpha = opt$alpha,
                              nfolds = opt$nfolds, keep = TRUE,
                              standardize = FALSE),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    idx <- which(fit$lambda == fit$lambda.min)[1]
    pred <- fit$fit.preval[, idx]
    if (stats::sd(pred) == 0) return(NULL)
    ct <- suppressWarnings(cor.test(y, pred))
    rsq <- as.numeric(ct$estimate^2)
    cf <- coef(fit, s = "lambda.min")
    w <- as.numeric(cf)[-1]
    names(w) <- rownames(cf)[-1]
    list(w = w, rsq = rsq, pval = ct$p.value)
  }

  models <- list()
  me <- fit_column(as.numeric(yexpr))
  if (!is.null(me)) models[[gene]] <- me

  if (!is.null(Yiso) && ncol(Yiso) >= 1L) {
    if (ncol(Yiso) > opt$`max-isoforms`) {
      vv <- apply(Yiso, 2, stats::var)
      Yiso <- Yiso[, order(-vv)[1:opt$`max-isoforms`], drop = FALSE]
    }
    diag$n_isoforms_tested <- ncol(Yiso)
    if (ncol(Yiso) == 1L) {
      mi <- fit_column(Yiso[, 1])
      if (!is.null(mi)) models[[colnames(Yiso)[1]]] <- mi
    } else {
      fit <- tryCatch(cv.glmnet(Xs, Yiso, family = "mgaussian",
                                alpha = opt$alpha, nfolds = opt$nfolds,
                                keep = TRUE, standardize = FALSE),
                      error = function(e) NULL)
      if (!is.null(fit)) {
        idx <- which(fit$lambda == fit$lambda.min)[1]
        prev <- fit$fit.preval[, , idx, drop = FALSE]
        cfs <- coef(fit, s = "lambda.min")
        for (j in seq_len(ncol(Yiso))) {
          pred <- prev[, j, 1]
          if (stats::sd(pred) == 0) next
          ct <- suppressWarnings(cor.test(Yiso[, j], pred))
          w <- as.numeric(cfs[[j]])[-1]
          names(w) <- rownames(cfs[[j]])[-1]
          models[[colnames(Yiso)[j]]] <-
            list(w = w, rsq = as.numeric(ct$estimate^2), pval = ct$p.value)
        }
      }
    }
  }

  keep <- vapply(models, function(m) is.finite(m$rsq) && m$rsq > opt$`r2-min`,
                 logical(1))
  models <- models[keep]
  diag$n_models_kept <- length(models)
  if (length(models) == 0L) return(finish("no_heritable_model"))
  diag$best_rsq <- max(vapply(models, `[[`, numeric(1), "rsq"))

  var_keep <- colnames(Xs)
  parts <- do.call(rbind, strsplit(var_keep, ":", fixed = TRUE))
  if (ncol(parts) != 4L)
    return(finish("error", "non-canonical variant ID encountered"))
  snps_full <- data.frame(V1 = parts[, 1], V2 = var_keep, V3 = 0,
                          V4 = as.integer(parts[, 2]), V5 = parts[, 4],
                          V6 = parts[, 3], stringsAsFactors = FALSE)
  for (mn in names(models)) {
    wgt.matrix <- matrix(0, nrow = length(var_keep), ncol = 1,
                         dimnames = list(var_keep, mn))
    w <- models[[mn]]$w
    common <- intersect(var_keep, names(w))
    wgt.matrix[common, mn] <- w[common]
    cv.performance <- matrix(c(models[[mn]]$rsq, models[[mn]]$pval),
                             nrow = 2, dimnames = list(c("rsq", "pval"), mn))
    snps <- snps_full
    hsq <- NA_real_
    hsq.pv <- NA_real_
    N.tot <- nrow(Xs)
    save(wgt.matrix, snps, cv.performance, hsq, hsq.pv, N.tot,
         file = file.path(wgt_dir, paste0(mn, ".wgt.RDat")))
    pos_row <- data.frame(WGT = paste0(mn, ".wgt.RDat"), ID = mn, GENE = gene,
                          CHR = rec$chrom, P0 = rec$win_start, P1 = rec$win_end)
    fwrite(pos_row, pos_path, sep = "\t", append = file.exists(pos_path),
           col.names = !file.exists(pos_path))
  }
  finish("ok")
}

if (file.exists(pos_path)) file.remove(pos_path)
diags <- vector("list", nrow(prepared))
for (i in seq_len(nrow(prepared))) {
  rec <- prepared[i]
  if (i %% 25L == 1L || i == nrow(prepared))
    message(sprintf("[%d/%d] %s", i, nrow(prepared), rec$gene))
  diags[[i]] <- tryCatch(train_gene(rec), error = function(e) data.table(
    gene = rec$gene, weight_set = ws, status = "error",
    message = gsub("[\t\n\r]", " ", substr(conditionMessage(e), 1, 300)),
    n_variants = NA_integer_, n_isoforms_tested = NA_integer_,
    n_models_kept = NA_integer_, best_rsq = NA_real_, walltime_sec = NA_real_))
}

diag_all <- if (length(diags)) rbindlist(diags, fill = TRUE) else data.table(
  gene = character(), weight_set = character(), status = character(),
  message = character(), n_variants = integer(), n_isoforms_tested = integer(),
  n_models_kept = integer(), best_rsq = numeric(), walltime_sec = numeric())
fwrite(diag_all, diag_path, sep = "\t")
message("status counts:")
if (nrow(diag_all)) print(diag_all[, .N, by = status][order(-N)])
message("diagnostics -> ", diag_path)
