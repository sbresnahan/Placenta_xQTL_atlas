#!/usr/bin/env Rscript
# =============================================================================
# 46_isotwas_train.R — isoTWAS / TWAS weight training worker (one shard)
# =============================================================================
# Objective 1.6: train SNP-based prediction weights for placental gene
# expression and, jointly, for all isoforms of each gene (isoTWAS:
# multivariate elastic net, glmnet family = "mgaussian", alpha = 0.5).
#
# Per gene x weight set:
#   * cis window = gene body +/- --cis-window (default 1 Mb)
#   * genotypes: plink2 --export A dosages from the ancestry pgen(s); the
#     pooled weight set stacks the per-ancestry dosage blocks on their shared
#     variants (no merged pgen required)
#   * phenotypes: harmonized BEDs (expression = gene level;
#     isoform_expression = {gene_id}__{transcript_id} rows), residualized on
#     the modality covariates ({ANC}_covariates_{MOD}.tsv); for the pooled
#     weight set each ancestry is residualized/z-scored separately, stacked,
#     and an ancestry indicator is residualized out of X and Y
#   * model: cv.glmnet (5-fold CV), weights at lambda.min; genotypes are
#     standardized BEFORE fitting and glmnet standardize = FALSE so weights
#     live on the standardized-genotype scale expected by FUSION.assoc_test.R
#   * retention: per-column cross-validated R^2 (prevalidated predictions)
#     > --r2-min (default 0.01)
#
# Outputs (under {outdir}):
#   weights/{WS}/genes/{gene}.wgt.RDat   FUSION format: wgt.matrix, snps,
#                                        cv.performance, hsq, hsq.pv, N.tot
#   weights/{WS}/shard-{idx}.pos         FUSION .pos fragment for the shard
#   diagnostics/isotwas_{WS}.shard-{idx}.diagnostics.tsv
#
# Usage:
#   Rscript 46_isotwas_train.R --weight-set EAS --shard-index 1 --n-shards 100 \
#     --qtl-dir $QTL_DIR --outdir $RESULTS_DIR/isotwas
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(glmnet)
})

option_list <- list(
  make_option("--weight-set", type = "character",
              help = "EAS, EUR, or pooled"),
  make_option("--shard-index", type = "integer", default = 1L),
  make_option("--n-shards", type = "integer", default = 1L),
  make_option("--qtl-dir", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--cis-window", type = "double", default = 1e6),
  make_option("--alpha", type = "double", default = 0.5),
  make_option("--r2-min", type = "double", default = 0.01),
  make_option("--nfolds", type = "integer", default = 5L),
  make_option("--min-variants", type = "integer", default = 10L),
  make_option("--max-isoforms", type = "integer", default = 50L,
              help = "cap isoforms per gene (most variable kept)"),
  make_option("--gene-list", type = "character", default = "",
              help = "optional file of gene IDs to restrict training"),
  make_option("--plink2", type = "character", default = "plink2"),
  make_option("--tmp-dir", type = "character", default = tempdir())
)
opt <- parse_args(OptionParser(option_list = option_list))

.args_all <- commandArgs(trailingOnly = FALSE)
.file_arg <- sub("^--file=", "", grep("^--file=", .args_all, value = TRUE))
.scripts_dir <- Sys.getenv("SCRIPTS_DIR",
                           unset = if (length(.file_arg)) dirname(.file_arg[1]) else getwd())
source(file.path(.scripts_dir, "coloc_common.R"))

ws <- opt$`weight-set`
ancestries <- if (ws == "pooled") c("EAS", "EUR") else ws
if (!all(ancestries %in% c("EAS", "EUR", "AFR", "AMR", "SAS")))
  stop("unknown weight set: ", ws)

dir.create(file.path(opt$outdir, "diagnostics"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$outdir, "weights", ws, "genes"),
           showWarnings = FALSE, recursive = TRUE)
diag_path <- file.path(opt$outdir, "diagnostics",
                       sprintf("isotwas_%s.shard-%04d.diagnostics.tsv", ws,
                               opt$`shard-index`))
pos_path <- file.path(opt$outdir, "weights", ws,
                      sprintf("shard-%04d.pos", opt$`shard-index`))

# --- data loading (once per shard) -------------------------------------------

read_bed <- function(path) {
  dt <- fread(path)
  setnames(dt, 1:4, c("chr", "start", "end", "phenotype_id"))
  dt
}

# covariates: rows = covariates, cols = samples -> samples x covariates matrix
read_covariates <- function(path) {
  dt <- fread(path)
  cn <- names(dt)
  id_col <- cn[1]
  samples <- cn[-1]
  mat <- t(as.matrix(dt[, -1, with = FALSE]))
  rownames(mat) <- samples
  colnames(mat) <- dt[[id_col]]
  mat
}

message("loading BEDs and covariates for: ", paste(ancestries, collapse = ", "))
beds <- covs <- list()
bed_samples <- list()
pgens <- character(0)
for (anc in ancestries) {
  beds[[anc]] <- list(
    expression = read_bed(file.path(opt$`qtl-dir`, sprintf("%s_expression.bed.gz", anc))),
    isoform_expression = read_bed(file.path(opt$`qtl-dir`, sprintf("%s_isoform_expression.bed.gz", anc))))
  bed_samples[[anc]] <- setdiff(names(beds[[anc]]$expression),
                                c("chr", "start", "end", "phenotype_id"))
  covs[[anc]] <- list(
    expression = read_covariates(file.path(
      opt$`qtl-dir`, sprintf("%s_covariates_expression.tsv", anc))),
    isoform_expression = read_covariates(file.path(
      opt$`qtl-dir`, sprintf("%s_covariates_isoform_expression.tsv", anc))))
  pgens[anc] <- file.path(opt$`qtl-dir`, paste0(anc, "_qtl"))
}


# gene universe: expression phenotype_ids present in every ancestry
gene_sets <- lapply(ancestries, function(a) beds[[a]]$expression$phenotype_id)
genes <- Reduce(intersect, gene_sets)
if (nzchar(opt$`gene-list`)) {
  gl <- readLines(opt$`gene-list`)
  genes <- intersect(genes, gl)
}
genes <- sort(genes)
shard_rows <- which(((seq_along(genes) - 1L) %% opt$`n-shards`) + 1L ==
                    opt$`shard-index`)
genes <- genes[shard_rows]
message(sprintf("[shard %d/%d] %d genes, weight set %s", opt$`shard-index`,
                opt$`n-shards`, length(genes), ws))

# variant positions per pgen (cached)
pvar_cache <- new.env()
get_pvar <- function(pgen) {
  if (exists(pgen, envir = pvar_cache)) return(get(pgen, envir = pvar_cache))
  pv <- fread(paste0(pgen, ".pvar"), select = c(1, 2, 3),
              col.names = c("chrom", "pos", "id"))
  pv[, chrom := sub("^chr", "", as.character(chrom))]
  assign(pgen, pv, envir = pvar_cache)
  pv
}

# residualize M (rows = samples) on covariates C (rows = samples)
residualize <- function(M, C) {
  common <- intersect(rownames(M), rownames(C))
  M <- M[common, , drop = FALSE]
  C <- C[common, , drop = FALSE]
  qr_dec <- qr(cbind(1, C))
  qr.resid(qr_dec, M)
}

# --- per-gene worker ---------------------------------------------------------

train_gene <- function(gene, tmp_dir) {
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
  if (file.exists(rdat_path)) return(finish("ok", "already done"))

  # gene coordinates from the first ancestry's expression BED
  grow <- beds[[ancestries[1]]]$expression[phenotype_id == gene]
  if (nrow(grow) == 0L) return(finish("error", "gene missing from BED"))
  chrom <- sub("^chr", "", as.character(grow$chr[1]))
  win_start <- max(0, grow$start[1] - opt$`cis-window`)
  win_end <- grow$end[1] + opt$`cis-window`

  # cis variants per panel
  var_lists <- list()
  for (nm in names(pgens)) {
    pv <- get_pvar(pgens[nm])
    var_lists[[nm]] <- pv[pv$chrom == chrom & pv$pos >= win_start &
                          pv$pos <= win_end]$id
  }
  var_ids <- Reduce(intersect, var_lists)
  if (length(var_ids) < opt$`min-variants`)
    return(finish("too_few_variants"))

  # dosages per ancestry (pooled stacks both ancestry blocks below)
  X_list <- list()
  for (nm in names(pgens)) {
    X <- export_dosages(opt$plink2, pgens[nm], NULL, var_ids,
                        file.path(tmp_dir, paste0("twas_", nm)))
    if (is.null(X)) return(finish("dosage_failed", nm))
    X_list[[nm]] <- X
  }

  # assemble training matrices per ancestry, then combine ------------------
  # Per ancestry: samples = intersection of dosage, expression, isoform and
  # covariate samples. X is residualized on the UNION of the expression and
  # isoform covariate sets (projection is harmless for prediction); each
  # phenotype is residualized on its own modality covariates, then z-scored.
  Xtr_list <- Yiso_list <- list()
  yexpr_list <- list()
  for (anc in ancestries) {
    X <- X_list[[anc]]
    rownames(X) <- attr(X, "samples")

    bed_e <- beds[[anc]]$expression
    bed_i <- beds[[anc]]$isoform_expression
    erow <- bed_e[phenotype_id == gene]
    if (nrow(erow) == 0L) next
    e_samples <- setdiff(names(bed_e), c("chr", "start", "end", "phenotype_id"))
    e_val <- as.numeric(erow[1, ..e_samples])
    names(e_val) <- e_samples
    irows <- bed_i[startsWith(phenotype_id, paste0(gene, "__"))]
    iso_mat <- NULL
    if (nrow(irows) > 0L) {
      i_samples <- setdiff(names(bed_i), c("chr", "start", "end", "phenotype_id"))
      iso_mat <- t(as.matrix(irows[, ..i_samples]))
      colnames(iso_mat) <- irows$phenotype_id
      rownames(iso_mat) <- i_samples
    }

    Ce <- covs[[anc]]$expression
    Ci <- covs[[anc]]$isoform_expression
    common <- Reduce(intersect, c(list(rownames(X), names(e_val),
                                       rownames(Ce), rownames(Ci)),
                                  if (!is.null(iso_mat)) list(rownames(iso_mat))
                                  else list(NULL)))
    if (length(common) < 20) next
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

  # stack ancestries
  Xtr <- do.call(rbind, Xtr_list)
  yexpr <- unlist(yexpr_list)
  # isoTWAS needs every isoform measured in every ancestry block
  Yiso <- NULL
  if (!any(vapply(Yiso_list, is.null, logical(1)))) {
    iso_ids <- Reduce(intersect, lapply(Yiso_list, colnames))
    if (length(iso_ids) > 0L) {
      Yiso <- do.call(rbind, lapply(Yiso_list, function(m) m[, iso_ids, drop = FALSE]))
      if (nrow(Yiso) != nrow(Xtr)) Yiso <- NULL  # paranoid
    }
  }

  # pooled: residualize ancestry indicator out of X and Y
  if (ws == "pooled" && length(Xtr_list) > 1L) {
    anc_dummy <- rep(seq_along(Xtr_list), vapply(Xtr_list, nrow, integer(1)))
    D <- cbind(1, model.matrix(~ factor(anc_dummy) - 1))
    Xtr <- qr.resid(qr(D), Xtr)
    yexpr <- as.numeric(qr.resid(qr(D), as.matrix(yexpr)))
    if (!is.null(Yiso)) Yiso <- qr.resid(qr(D), Yiso)
  }

  diag$n_variants <- ncol(Xtr)
  if (ncol(Xtr) < opt$`min-variants`) return(finish("too_few_variants"))

  # standardize genotypes; weights must live on the standardized scale for
  # FUSION.assoc_test.R (which scale()s the reference panel)
  Xs <- scale(Xtr)
  Xs <- Xs[, apply(Xs, 2, function(v) all(is.finite(v))), drop = FALSE]
  if (ncol(Xs) < opt$`min-variants`) return(finish("too_few_variants"))

  fit_column <- function(y) {
    # returns list(w = named weight vector, rsq, pval) or NULL
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

  models <- list()  # name -> list(w, rsq, pval)

  # per-gene expression model
  me <- fit_column(as.numeric(yexpr))
  if (!is.null(me)) models[[gene]] <- me

  # isoTWAS: multivariate elastic net across isoforms
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
        # mgaussian fit.preval is [samples, responses, lambdas]
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

  # retention gate
  keep <- vapply(models, function(m) is.finite(m$rsq) && m$rsq > opt$`r2-min`,
                 logical(1))
  models <- models[keep]
  diag$n_models_kept <- length(models)
  if (length(models) == 0L) return(finish("no_heritable_model"))
  diag$best_rsq <- max(vapply(models, `[[`, numeric(1), "rsq"))

  # FUSION .wgt.RDat files: ONE PER MODEL (isoTWAS convention: each isoform
  # is tested separately by FUSION.assoc_test.R and aggregated at gene level
  # downstream). File name = model ID; the .pos GENE column links isoforms
  # back to their gene.
  var_keep <- colnames(Xs)
  parts <- do.call(rbind, strsplit(var_keep, ":"))
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
                          CHR = chrom, P0 = win_start, P1 = win_end)
    fwrite(pos_row, pos_path, sep = "\t", append = file.exists(pos_path),
           col.names = !file.exists(pos_path))
  }
  finish("ok")
}

# --- run shard ---------------------------------------------------------------

tmp_dir <- tempfile(pattern = "isotwas_", tmpdir = opt$`tmp-dir`)
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
if (file.exists(pos_path)) file.remove(pos_path)

diags <- vector("list", length(genes))
for (i in seq_along(genes)) {
  if (i %% 25 == 1L || i == length(genes))
    message(sprintf("[%d/%d] %s", i, length(genes), genes[i]))
  diags[[i]] <- tryCatch(train_gene(genes[i], tmp_dir),
                         error = function(e) data.table(
                           gene = genes[i], weight_set = ws, status = "error",
                           message = gsub("[\t\n\r]", " ",
                                          substr(conditionMessage(e), 1, 300)),
                           n_variants = NA_integer_, n_isoforms_tested = NA_integer_,
                           n_models_kept = NA_integer_, best_rsq = NA_real_,
                           walltime_sec = NA_real_))
}

diag_all <- rbindlist(diags, fill = TRUE)
fwrite(diag_all, diag_path, sep = "\t")
message("status counts:")
print(diag_all[, .N, by = status][order(-N)])
message("diagnostics -> ", diag_path)
