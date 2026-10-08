# =============================================================================
# coloc_common.R — shared helpers for 41_susie_coloc.R and 43_colocboost.R
# =============================================================================
# Sourced by both workers; expects data.table to be attached. All functions
# take explicit arguments (no reliance on caller globals).
# =============================================================================

sanitize_id <- function(x) gsub("[^A-Za-z0-9._+-]", "_", x)

# cached header reader for bgzipped TSVs
.make_header_cache <- function() {
  cache <- new.env()
  function(path) {
    if (exists(path, envir = cache)) return(get(path, envir = cache))
    h <- names(data.table::fread(cmd = sprintf("zcat %s | head -1",
                                               shQuote(path))))
    assign(path, h, envir = cache)
    h
  }
}
get_header <- .make_header_cache()

# tabix region slice; returns NULL when empty/unreadable
tabix_slice <- function(path, chrom, start, end) {
  hdr <- get_header(path)
  cmd <- sprintf("tabix %s %s:%s-%s", shQuote(path), chrom, start, end)
  dt <- tryCatch(data.table::fread(cmd = cmd, header = FALSE, col.names = hdr),
                 error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L) return(NULL)
  dt
}

# plink2 --r-unphased square -> correlation matrix with variant-ID dimnames.
# Returns list(LD = matrix, present = ids) or NULL on failure.
run_ld <- function(plink2, pgen, keep, var_ids, out_prefix) {
  list_file <- paste0(out_prefix, ".extract.txt")
  writeLines(var_ids, list_file)
  args <- c("--pfile", pgen,
            "--extract", list_file,
            "--r-unphased", "square",
            "--threads", "1", "--silent",
            "--out", out_prefix)
  if (!is.null(keep) && length(keep) == 1L && !is.na(keep) && nzchar(keep) &&
      file.exists(keep))
    args <- c(args, "--keep", keep)
  status <- system2(plink2, args, stdout = FALSE, stderr = FALSE)
  vcor <- paste0(out_prefix, ".unphased.vcor1")
  vvars <- paste0(out_prefix, ".unphased.vcor1.vars")
  if (status != 0L || !file.exists(vcor) || !file.exists(vvars)) return(NULL)
  ids <- readLines(vvars, warn = FALSE)
  ids <- ids[!grepl("^#", ids)]
  if (length(ids) == 0L) return(NULL)
  mat <- as.matrix(data.table::fread(vcor, header = FALSE))
  if (nrow(mat) != length(ids) || ncol(mat) != length(ids)) return(NULL)
  dimnames(mat) <- list(ids, ids)
  list(LD = mat, present = ids)
}

# plink2 --export A -> samples x variants dosage matrix with variant-ID
# colnames, restricted to `var_ids`; column-mean imputed; monomorphic columns
# dropped; NULL on failure
export_dosages <- function(plink2, pgen, keep, var_ids, out_prefix) {
  list_file <- paste0(out_prefix, ".extract.txt")
  writeLines(var_ids, list_file)
  args <- c("--pfile", pgen, "--extract", list_file,
            "--export", "A", "--threads", "1", "--silent",
            "--out", out_prefix)
  if (!is.null(keep) && length(keep) == 1L && !is.na(keep) && nzchar(keep) &&
      file.exists(keep))
    args <- c(args, "--keep", keep)
  status <- system2(plink2, args, stdout = FALSE, stderr = FALSE)
  raw_path <- paste0(out_prefix, ".raw")
  if (status != 0L || !file.exists(raw_path)) return(NULL)
  raw <- data.table::fread(raw_path)
  sample_ids <- raw[["IID"]]
  drop <- intersect(names(raw), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
  raw <- as.data.frame(raw)[, setdiff(names(raw), drop), drop = FALSE]
  ids <- names(raw)
  # plink2 may suffix the allele to duplicated IDs; map back to requested IDs
  ids_fix <- ids
  no_suffix <- sub("_[^_]+$", "", ids)
  need_fix <- !(ids %in% var_ids) & (no_suffix %in% var_ids)
  ids_fix[need_fix] <- no_suffix[need_fix]
  names(raw) <- ids_fix
  raw <- raw[, intersect(var_ids, names(raw)), drop = FALSE]
  if (ncol(raw) == 0L) return(NULL)
  mat <- as.matrix(raw)
  # column-mean imputation for missing dosages
  for (j in seq_len(ncol(mat))) {
    na_idx <- is.na(mat[, j])
    if (any(na_idx)) {
      mu <- mean(mat[, j], na.rm = TRUE)
      if (!is.finite(mu)) mu <- 0
      mat[na_idx, j] <- mu
    }
  }
  # drop monomorphic columns (zero variance breaks standardization)
  sds <- apply(mat, 2, stats::sd)
  mat <- mat[, sds > 0, drop = FALSE]
  if (ncol(mat) == 0L) return(NULL)
  attr(mat, "samples") <- sample_ids
  mat
}

# sample size of a pgen (non-comment rows of .psam), cached
.make_psam_n_cache <- function() {
  cache <- new.env()
  function(pgen) {
    if (exists(pgen, envir = cache)) return(get(pgen, envir = cache))
    n <- tryCatch({
      con <- file(paste0(pgen, ".psam"), "r")
      on.exit(close(con))
      cnt <- 0L
      while (length(l <- readLines(con, n = 1L)) > 0L)
        if (!startsWith(l, "#")) cnt <- cnt + 1L
      cnt
    }, error = function(e) NA_integer_)
    assign(pgen, n, envir = cache)
    n
  }
}
psam_n <- .make_psam_n_cache()
