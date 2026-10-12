#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit)) return(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = TRUE))
  normalizePath("prepare_book.R", mustWork = TRUE)
}

book_dir <- dirname(script_path())
repo_root <- normalizePath(file.path(book_dir, ".."), mustWork = TRUE)
reports_dir <- file.path(repo_root, "reports")
stage_dir <- file.path(book_dir, ".staging")
assets_dir <- file.path(book_dir, "book_assets")
manifest_path <- file.path(book_dir, "report_manifest.tsv")

if (!dir.exists(reports_dir)) {
  stop(
    "Local reports directory not found: ", reports_dir,
    "\nThe reports directory is intentionally git-ignored. Build the book from the checkout ",
    "where the reports and report inputs were generated."
  )
}
if (!file.exists(manifest_path)) stop("Missing report manifest: ", manifest_path)
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required. Install it with install.packages('rmarkdown').")
}

manifest <- utils::read.delim(
  manifest_path,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  encoding = "UTF-8"
)
required_cols <- c("part_order", "chapter_order", "part_title", "source", "kind", "title", "id")
missing_cols <- setdiff(required_cols, names(manifest))
if (length(missing_cols)) stop("Manifest missing columns: ", paste(missing_cols, collapse = ", "))
manifest <- manifest[order(manifest$part_order, manifest$chapter_order), , drop = FALSE]

# Resolve the current local report archive once, rather than baking a dated
# archive name into the website builder. This is also used to repair a few
# historical report defaults that pointed at paths used only during report
# development (for example data/extras or coloc_aggregated.zip).
archive_candidates <- list.dirs(reports_dir, recursive = FALSE, full.names = TRUE)
archive_candidates <- archive_candidates[grepl("^placenta_xqtl_report_inputs_[0-9]+$", basename(archive_candidates))]
archive_dir <- if (length(archive_candidates)) sort(archive_candidates)[length(archive_candidates)] else NA_character_
archive_rel <- if (!is.na(archive_dir)) basename(archive_dir) else NA_character_

param_overrides <- list()
if (!is.na(archive_rel)) {
  param_overrides[["xqtl-mapping"]] <- list(
    archive_dir = archive_rel
  )
  param_overrides[["fine-mapping"]] <- list(
    agg_dir = file.path(archive_rel, "data", "extras"),
    data_dir = file.path(archive_rel, "data", "results"),
    gene_map_path = file.path(archive_rel, "gene_map.tsv")
  )
  param_overrides[["colocalization"]] <- list(
    coloc_input = file.path("coloc", "aggregated"),
    colocboost_path = file.path("coloc", "aggregated", "colocboost_clusters.tsv.gz"),
    gene_map_path = file.path(archive_rel, "gene_map.tsv")
  )
  param_overrides[["twas"]] <- list(
    twas_input = file.path("coloc", "aggregated"),
    gene_map_path = file.path(archive_rel, "gene_map.tsv"),
    coloc_best_path = file.path("coloc", "aggregated", "coloc_best.tsv.gz"),
    colocboost_path = file.path("coloc", "aggregated", "colocboost_clusters.tsv.gz")
  )
}
param_overrides[["cell-type"]] <- list(results_dir = "celltype_diagnostics_bundle")
param_overrides[["gxe-ppbmi"]] <- list(gxe_dir = "gxe")
param_overrides[["integration"]] <- list(results_root = ".")

discovered_sources <- c(
  list.files(reports_dir, pattern = "^report_.*\\.Rmd$", full.names = FALSE),
  list.files(reports_dir, pattern = "^literature_.*\\.md$", full.names = FALSE)
)
discovered_sources <- sort(unique(discovered_sources))
unlisted_sources <- setdiff(discovered_sources, manifest$source)
if (length(unlisted_sources)) {
  stop(
    "Scientific report source(s) exist in reports/ but are not in report_manifest.tsv:\n  - ",
    paste(unlisted_sources, collapse = "\n  - "),
    "\nAdd them to the manifest so the published book remains complete."
  )
}

source_paths <- file.path(reports_dir, manifest$source)
missing_sources <- manifest$source[!file.exists(source_paths)]
if (length(missing_sources)) {
  stop(
    "The following book sources are missing from reports/:\n  - ",
    paste(missing_sources, collapse = "\n  - "),
    "\n\nThe integration report is expected at reports/report_integration.Rmd."
  )
}

# Remove only generated book inputs/resources; never touch reports/.
if (dir.exists(stage_dir)) unlink(stage_dir, recursive = TRUE, force = TRUE)
if (dir.exists(assets_dir)) unlink(assets_dir, recursive = TRUE, force = TRUE)
unlink(list.files(book_dir, pattern = "^zbook_.*\\.Rmd$", full.names = TRUE), force = TRUE)
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage_dir, ".intermediates"), recursive = TRUE, showWarnings = FALSE)
dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)

inline_r_matches <- function(lines) {
  # knitr inline expressions must stay on one physical line. A few historical
  # reports wrap immediately after `r, e.g. `r\nn_call_traits`. knitr leaves
  # those expressions untouched in the first Markdown render, after which the
  # bookdown pass tries to evaluate them without the report setup environment.
  text <- paste(lines, collapse = "\n")
  m <- gregexpr("(?s)`r\\s+(?:(?!`).)*`", text, perl = TRUE)
  if (length(m) == 0L || m[[1]][1] == -1L) return(character())
  regmatches(text, m)[[1]]
}

normalize_multiline_inline_r <- function(lines) {
  text <- paste(lines, collapse = "\n")
  m <- gregexpr("(?s)`r\\s+(?:(?!`).)*`", text, perl = TRUE)
  if (length(m) == 0L || m[[1]][1] == -1L) return(lines)

  hits <- regmatches(text, m)[[1]]
  fixed <- vapply(hits, function(hit) {
    if (!grepl("[\\r\\n]", hit)) return(hit)
    gsub("[\\r\\n]+[ \\t]*", " ", hit, perl = TRUE)
  }, character(1))
  regmatches(text, m) <- list(fixed)
  strsplit(text, "\n", fixed = TRUE)[[1]]
}

strip_front_matter <- function(x) {
  if (!length(x)) return(x)
  first <- which(nzchar(trimws(x)))[1]
  if (is.na(first) || trimws(x[first]) != "---") return(x)
  if (first == length(x)) return(x)
  rest <- which(trimws(x[(first + 1):length(x)]) %in% c("---", "..."))
  if (!length(rest)) return(x)
  last <- first + rest[1]
  x[-seq.int(first, last)]
}

front_matter_title <- function(path) {
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (!length(x) || trimws(x[1]) != "---") return(NA_character_)
  if (length(x) < 2) return(NA_character_)
  end <- which(trimws(x[-1]) %in% c("---", "..."))[1]
  if (is.na(end)) return(NA_character_)
  yaml <- x[2:end]
  title <- grep("^title\\s*:", yaml, value = TRUE)
  if (!length(title)) return(NA_character_)
  title <- sub("^title\\s*:\\s*", "", title[1])
  title <- sub("^['\"]", "", title)
  title <- sub("['\"]$", "", title)
  title
}

plain_heading <- function(x) {
  x <- sub("^#+\\s*", "", x)
  x <- sub("\\s*\\{#[^}]+\\}\\s*$", "", x)
  x <- gsub("[`*_]", "", x)
  x <- gsub("[^[:alnum:]]+", " ", tolower(x))
  trimws(x)
}

maybe_drop_title <- function(x, candidates = character(), drop_first_h1 = FALSE) {
  if (!length(x)) return(x)
  nz <- which(nzchar(trimws(x)))
  if (!length(nz)) return(x)
  i <- nz[1]

  if (grepl("^#\\s+", x[i])) {
    h <- plain_heading(x[i])
    cand <- plain_heading(candidates[!is.na(candidates) & nzchar(candidates)])
    matched <- any(vapply(cand, function(z) {
      nzchar(z) && (identical(h, z) || startsWith(h, z) || startsWith(z, h))
    }, logical(1)))
    if (drop_first_h1 || matched) return(x[-i])
  }

  if (i < length(x) && grepl("^=+\\s*$", x[i + 1])) {
    h <- plain_heading(x[i])
    cand <- plain_heading(candidates[!is.na(candidates) & nzchar(candidates)])
    matched <- any(vapply(cand, function(z) {
      nzchar(z) && (identical(h, z) || startsWith(h, z) || startsWith(z, h))
    }, logical(1)))
    if (drop_first_h1 || matched) return(x[-c(i, i + 1)])
  }
  x
}

strip_style_blocks <- function(x) {
  out <- character()
  in_style <- FALSE
  in_css_fence <- FALSE
  for (line in x) {
    t <- trimws(line)
    if (!in_style && grepl("^<style(?:\\s|>)", t, ignore.case = TRUE)) {
      in_style <- !grepl("</style>", t, ignore.case = TRUE)
      next
    }
    if (in_style) {
      if (grepl("</style>", t, ignore.case = TRUE)) in_style <- FALSE
      next
    }
    if (!in_css_fence && grepl("^```(?:\\{css|css)(?:[, }]|$)", t, ignore.case = TRUE)) {
      in_css_fence <- TRUE
      next
    }
    if (in_css_fence) {
      if (grepl("^```\\s*$", t)) in_css_fence <- FALSE
      next
    }
    out <- c(out, line)
  }
  out
}

demote_headings <- function(x) {
  out <- x
  in_fence <- FALSE
  fence_token <- NULL
  for (i in seq_along(x)) {
    t <- trimws(x[i])
    if (grepl("^(```|~~~)", t)) {
      token <- substr(t, 1, 3)
      if (!in_fence) {
        in_fence <- TRUE
        fence_token <- token
      } else if (identical(token, fence_token)) {
        in_fence <- FALSE
        fence_token <- NULL
      }
      next
    }
    if (!in_fence && grepl("^#{1,5}\\s+", x[i])) out[i] <- paste0("#", x[i])
  }
  out
}

is_external_resource <- function(x) {
  grepl("^(?:https?:|mailto:|data:|#)", x, ignore.case = TRUE) || grepl("^//", x)
}

safe_resource_rel <- function(x) {
  x <- sub("[?#].*$", "", x)
  x <- gsub("\\\\", "/", x)
  x <- gsub("^\\./", "", x)
  while (grepl("^\\.\\./", x)) x <- sub("^\\.\\./", "", x)
  if (!nzchar(x) || grepl("^/", x) || grepl("^[A-Za-z]:/", x)) x <- basename(x)
  x <- gsub("[^A-Za-z0-9._/-]", "_", x)
  if (!nzchar(x)) x <- "resource"
  x
}

relative_to <- function(path, root) {
  p <- normalizePath(path, mustWork = TRUE)
  r <- normalizePath(root, mustWork = TRUE)
  prefix <- paste0(r, .Platform$file.sep)
  if (!startsWith(p, prefix)) return(NA_character_)
  substring(p, nchar(prefix) + 1)
}

extract_resource_paths <- function(text) {
  paths <- character()

  md_match <- gregexpr("!\\[[^]]*\\]\\((?:<)?[^ )>]+", text, perl = TRUE)
  md_matches <- regmatches(text, md_match)[[1]]
  if (length(md_matches)) {
    paths <- c(paths, sub("^.*\\]\\((?:<)?", "", md_matches, perl = TRUE))
  }

  html_match <- gregexpr(
    "<(?:img|source)[^>]+(?:src|srcset)=['\"][^'\"]+",
    text,
    perl = TRUE,
    ignore.case = TRUE
  )
  html_matches <- regmatches(text, html_match)[[1]]
  if (length(html_matches)) {
    hp <- sub("^.*(?:src|srcset)=['\"]", "", html_matches, perl = TRUE, ignore.case = TRUE)
    hp <- sub("\\s+.*$", "", hp)
    paths <- c(paths, hp)
  }

  unique(paths[nzchar(paths)])
}

stage_resources <- function(x, raw_md_dir, source_dir, chapter_id) {
  if (!length(x)) return(x)
  text <- paste(x, collapse = "\n")
  paths <- extract_resource_paths(text)
  if (!length(paths)) return(x)

  for (old in paths) {
    clean <- gsub("^<|>$", "", old)
    if (is_external_resource(clean)) next
    clean_file <- sub("[?#].*$", "", clean)

    candidates <- c(file.path(raw_md_dir, clean_file), file.path(source_dir, clean_file))
    candidates <- candidates[file.exists(candidates) & !dir.exists(candidates)]
    if (!length(candidates)) next
    src <- normalizePath(candidates[1], mustWork = TRUE)

    rel_stage <- relative_to(src, stage_dir)
    if (!is.na(rel_stage)) {
      safe_rel <- file.path("generated", gsub("\\\\", "/", rel_stage))
    } else {
      safe_rel <- safe_resource_rel(clean_file)
    }

    dest <- file.path(assets_dir, chapter_id, safe_rel)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(src, dest, overwrite = TRUE)) stop("Could not stage resource: ", src)

    new <- paste0(
      "book_assets/", chapter_id, "/",
      gsub("\\\\", "/", safe_rel)
    )
    text <- gsub(old, new, text, fixed = TRUE)
  }

  strsplit(text, "\n", fixed = TRUE)[[1]]
}

chapter_filename <- function(part_order, chapter_order, id) {
  idx <- (as.integer(part_order) - 1L) * 100L + as.integer(chapter_order)
  sprintf("zbook_%03d_%s.Rmd", idx, id)
}

# Part dividers are generated at the book root so all local resource paths are
# resolved consistently by bookdown/rmarkdown.
parts <- unique(manifest[c("part_order", "part_title")])
parts <- parts[order(parts$part_order), , drop = FALSE]
for (i in seq_len(nrow(parts))) {
  part_file <- file.path(book_dir, sprintf("zbook_part_%02d.Rmd", as.integer(parts$part_order[i])))
  writeLines(c(sprintf("# (PART) %s {-}", parts$part_title[i]), ""), part_file, useBytes = TRUE)
}

status <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  src <- file.path(reports_dir, row$source)
  out_name <- chapter_filename(row$part_order, row$chapter_order, row$id)
  out_path <- file.path(book_dir, out_name)

  message(sprintf("[%d/%d] %s", i, nrow(manifest), row$title))

  if (tolower(row$kind) == "rmd") {
    raw_name <- paste0("raw_", row$id, ".md")
    intermediate_dir <- file.path(stage_dir, ".intermediates", row$id)
    dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

    fmt <- rmarkdown::md_document(variant = "markdown_github", preserve_yaml = FALSE)

    render_params <- param_overrides[[row$id]]
    if (is.null(render_params)) render_params <- list()

    # Normalize malformed historical inline-R wrapping in a temporary source
    # beside the original report. Keeping the temporary input in reports/
    # preserves source-relative behavior without modifying the ignored report.
    render_src <- src
    src_lines <- readLines(src, warn = FALSE, encoding = "UTF-8")
    normalized_src_lines <- normalize_multiline_inline_r(src_lines)
    if (!identical(src_lines, normalized_src_lines)) {
      render_src <- file.path(reports_dir, paste0(".bookdown_render_", row$id, ".Rmd"))
      writeLines(normalized_src_lines, render_src, useBytes = TRUE)
      message("  normalized multiline inline R in ", row$source)
    }

    rendered <- tryCatch(
      rmarkdown::render(
        input = render_src,
        output_format = fmt,
        output_file = raw_name,
        params = render_params,
        output_dir = stage_dir,
        intermediates_dir = intermediate_dir,
        knit_root_dir = reports_dir,
        envir = new.env(parent = globalenv()),
        clean = TRUE,
        quiet = FALSE
      ),
      finally = {
        if (!identical(render_src, src) && file.exists(render_src)) unlink(render_src, force = TRUE)
      }
    )
    raw_path <- normalizePath(rendered, mustWork = TRUE)
    body <- readLines(raw_path, warn = FALSE, encoding = "UTF-8")

    # The second bookdown pass must never receive executable inline R from a
    # report: its setup chunks have already been evaluated and removed by the
    # first render. Fail here with a targeted message instead of an undefined
    # object error deep in bookdown.
    leftover_inline_r <- inline_r_matches(body)
    if (length(leftover_inline_r)) {
      preview <- paste(
        utils::head(gsub("[\\r\\n]+", " ", leftover_inline_r), 5L),
        collapse = "\n  - "
      )
      stop(
        "Intermediate Markdown for ", row$source,
        " still contains unevaluated inline R after the first render:\n  - ", preview,
        "\nKeep each inline expression on one physical line in the source report."
      )
    }
    body <- strip_front_matter(body)
    body <- maybe_drop_title(body, c(front_matter_title(src), row$title), drop_first_h1 = FALSE)
    body <- strip_style_blocks(body)
    body <- demote_headings(body)
    body <- stage_resources(body, dirname(raw_path), reports_dir, row$id)
  } else if (tolower(row$kind) == "md") {
    body <- readLines(src, warn = FALSE, encoding = "UTF-8")
    body <- strip_front_matter(body)
    body <- maybe_drop_title(body, row$title, drop_first_h1 = TRUE)
    body <- strip_style_blocks(body)
    body <- demote_headings(body)
    body <- stage_resources(body, reports_dir, reports_dir, row$id)
  } else {
    stop("Unsupported manifest kind '", row$kind, "' for ", row$source)
  }

  source_note <- sprintf(
    '<p class="source-note">Book chapter rebuilt locally from <code>reports/%s</code>.</p>',
    row$source
  )
  chapter <- c(
    sprintf("# %s {#%s}", row$title, row$id),
    "",
    source_note,
    "",
    body,
    ""
  )
  writeLines(chapter, out_path, useBytes = TRUE)

  status[[i]] <- data.frame(
    source = row$source,
    chapter = out_name,
    kind = row$kind,
    stringsAsFactors = FALSE
  )
}

status <- do.call(rbind, status)
utils::write.table(
  status,
  file.path(stage_dir, "prepared_sources.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Raw Markdown is no longer needed after chapter normalization. Generated
# resource files have already been copied into book_assets/.
unlink(list.files(stage_dir, pattern = "^raw_.*\\.md$", full.names = TRUE), force = TRUE)

message("Prepared ", nrow(manifest), " book chapters in ", book_dir)
