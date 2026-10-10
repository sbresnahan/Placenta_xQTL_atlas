# =============================================================================
# coloc_common.R — pure-R shared helpers for colocalization/TWAS workers
# =============================================================================
# IMPORTANT: R code in module 06 must not invoke command-line tools. External
# region extraction, PLINK LD/dosage generation, compression/indexing, etc. are
# performed by host-side shell/Python preparation steps before R is launched.
# This file therefore contains only R-native helpers.
# =============================================================================

sanitize_id <- function(x) gsub("[^A-Za-z0-9._+-]", "_", x)

# Sample size of a pgen, obtained by reading its text .psam companion directly.
# No external command is invoked. Cached because many tasks share a panel.
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
