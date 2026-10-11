#!/usr/bin/env Rscript
# =============================================================================
# 48b_fusion_assoc.R — resilient wrapper around FUSION.assoc_test.R
# =============================================================================
# FUSION.assoc_test.R (gusevlab/fusion_twas clone, not part of this repo) loads
# plink2R/optparse at startup but does NOT set .libPaths(), and bash-level
# R_LIBS_* / R_LIBS_USER are not reliably propagated to R on this cluster.
# This wrapper points .libPaths() at the lab R 4.3.1 package library INSIDE R
# and then evaluates a minimally patched copy of FUSION.assoc_test.R.
#
# The patch is applied in memory only; the external FUSION clone is not edited.
# It wraps FUSION's per-weight/model loop in tryCatch(), so a numerical/runtime
# error in one phenotype is recorded in FUSION.STATUS/FUSION.ERROR and the
# remaining phenotypes continue. FUSION's own per-phenotype cur.FAIL/skips are
# also marked failed. Errors outside that loop (bad inputs, unreadable LD, etc.)
# remain fatal because they invalidate the whole shard.
#
# FUSION.assoc_test.R parses its flags with optparse::parse_args(), which reads
# commandArgs(trailingOnly = TRUE); arguments passed after this wrapper on the
# Rscript command line are therefore seen by FUSION unchanged.
# =============================================================================

R_LIB <- "/rsrch5/home/epi/bhattacharya_lab/software/R_package_library/ubuntu/4.3.1"
if (!dir.exists(R_LIB)) stop("R package library does not exist: ", R_LIB)
.libPaths(c(R_LIB, .libPaths()))

fusion_dir <- Sys.getenv("FUSION_DIR", unset = "")
if (!nzchar(fusion_dir)) stop("FUSION_DIR is not set in the environment")
fusion_assoc <- file.path(fusion_dir, "FUSION.assoc_test.R")
if (!file.exists(fusion_assoc)) stop("FUSION.assoc_test.R not found: ", fusion_assoc)

src <- readLines(fusion_assoc, warn = FALSE)

find_unique <- function(pattern, label) {
    idx <- grep(pattern, src, perl = TRUE)
    if (length(idx) != 1L) {
        stop(
            "Cannot safely patch FUSION.assoc_test.R: expected exactly one ",
            label, ", found ", length(idx), ". External FUSION source may have changed."
        )
    }
    idx
}

# Anchor to stable statements in the upstream script instead of line numbers.
idx_failctr <- find_unique("^FAIL\\.ctr\\s*=\\s*0\\s*$", "FAIL.ctr initializer")
idx_loop <- find_unique(
    "^for \\( w in 1:nrow\\(wgtlist\\) \\) \\{\\s*$",
    "per-weight for-loop"
)
idx_failinc <- find_unique(
    "if \\( cur\\.FAIL \\) FAIL\\.ctr = FAIL\\.ctr \\+ 1",
    "per-weight FAIL.ctr increment"
)
idx_next <- find_unique("^\\s*next\\s*$", "per-weight next statement")

# The loop's closing brace is the first nonblank line after the FAIL.ctr line.
idx_loop_close <- idx_failinc + 1L
while (idx_loop_close <= length(src) && !nzchar(trimws(src[idx_loop_close]))) {
    idx_loop_close <- idx_loop_close + 1L
}
if (idx_loop_close > length(src) || trimws(src[idx_loop_close]) != "}") {
    stop("Cannot safely patch FUSION.assoc_test.R: per-weight loop closing brace not found")
}
if (!(idx_failctr < idx_loop && idx_loop < idx_next && idx_next < idx_failinc && idx_failinc < idx_loop_close)) {
    stop("Cannot safely patch FUSION.assoc_test.R: unexpected anchor ordering")
}

patched <- character(0)
for (i in seq_along(src)) {
    # Insert the tryCatch closure immediately before the original loop close.
    if (i == idx_loop_close) {
        patched <- c(
            patched,
            "    }, error = function(e) {",
            "        raw.msg <- gsub('[\\t\\r\\n]+', ' ', conditionMessage(e))",
            "        is.expected.skip <- identical(raw.msg, '__FUSION_SKIP_NO_PREDICTIVE_MODEL__')",
            "        msg <- if (is.expected.skip) 'FUSION did not have a predictive model; skipped' else raw.msg",
            "        wgt.file.err <- paste(opt$weights_dir, '/', wgtlist$WGT[w], sep='')",
            "        if (sum(names(wgtlist) == 'PANEL') == 1) out.tbl$PANEL[w] <<- wgtlist$PANEL[w]",
            "        out.tbl$FILE[w] <<- wgt.file.err",
            "        out.tbl$ID[w] <<- wgtlist$ID[w]",
            "        out.tbl$CHR[w] <<- wgtlist$CHR[w]",
            "        out.tbl$P0[w] <<- wgtlist$P0[w]",
            "        out.tbl$P1[w] <<- wgtlist$P1[w]",
            "        out.tbl$TWAS.Z[w] <<- NA_real_",
            "        out.tbl$TWAS.P[w] <<- NA_real_",
            "        out.tbl$FUSION.STATUS[w] <<- 'failed'",
            "        out.tbl$FUSION.ERROR[w] <<- msg",
            "        if (!is.expected.skip) FAIL.ctr <<- FAIL.ctr + 1L",
            "        cat(if (is.expected.skip) 'WARNING :' else 'ERROR :', 'phenotype', as.character(wgtlist$ID[w]), 'failed:', msg, '-- continuing\\n')",
            "    })"
        )
    }

    # The upstream no-predictive-model branch uses `next`. Once the loop body
    # is wrapped by tryCatch(), convert that one control transfer into a tagged
    # condition handled below; this avoids changing R's loop-control context.
    line_i <- src[i]
    if (i == idx_next) line_i <- "        stop('__FUSION_SKIP_NO_PREDICTIVE_MODEL__')"
    patched <- c(patched, line_i)

    # Add durable status/error columns to the normal FUSION output table.
    if (i == idx_failctr) {
        patched <- c(
            patched,
            "out.tbl$FUSION.STATUS = rep('pending', nrow(wgtlist))",
            "out.tbl$FUSION.ERROR = rep(NA_character_, nrow(wgtlist))"
        )
    }

    # Start the per-phenotype error boundary inside the existing for-loop.
    if (i == idx_loop) {
        patched <- c(patched, "    tryCatch({")
    }

    # If FUSION itself handled the phenotype via cur.FAIL, record that too.
    if (i == idx_failinc) {
        patched <- c(
            patched,
            "        out.tbl$FUSION.STATUS[w] = if (isTRUE(cur.FAIL)) 'failed' else 'ok'",
            "        if (isTRUE(cur.FAIL) && is.na(out.tbl$FUSION.ERROR[w])) out.tbl$FUSION.ERROR[w] = 'FUSION marked phenotype as failed/skipped; see shard log'"
        )
    }

    # Safety net: no completed output row should remain pending. Convert any
    # such row to an explicit failure before the output is written.
    if (i == idx_loop_close) {
        patched <- c(
            patched,
            "fusion.pending <- which(out.tbl$FUSION.STATUS == 'pending')",
            "if (length(fusion.pending) > 0) {",
            "    out.tbl$FUSION.STATUS[fusion.pending] <- 'failed'",
            "    out.tbl$FUSION.ERROR[fusion.pending] <- 'FUSION skipped phenotype before completing its output row; see shard log'",
            "    for (ww in fusion.pending) {",
            "        if (sum(names(wgtlist) == 'PANEL') == 1) out.tbl$PANEL[ww] <- wgtlist$PANEL[ww]",
            "        out.tbl$FILE[ww] <- paste(opt$weights_dir, '/', wgtlist$WGT[ww], sep='')",
            "        out.tbl$ID[ww] <- wgtlist$ID[ww]",
            "        out.tbl$CHR[ww] <- wgtlist$CHR[ww]",
            "        out.tbl$P0[ww] <- wgtlist$P0[ww]",
            "        out.tbl$P1[ww] <- wgtlist$P1[ww]",
            "        out.tbl$TWAS.Z[ww] <- NA_real_",
            "        out.tbl$TWAS.P[ww] <- NA_real_",
            "    }",
            "}"
        )
    }
}

# Evaluate in the global environment so FUSION's top-level option parsing and
# object behavior remain the same as direct execution.
eval(parse(text = patched, keep.source = FALSE), envir = .GlobalEnv)
