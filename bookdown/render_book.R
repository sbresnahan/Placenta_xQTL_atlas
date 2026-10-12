#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit)) return(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = TRUE))
  normalizePath("render_book.R", mustWork = TRUE)
}

book_dir <- dirname(script_path())
assets_dir <- file.path(book_dir, "book_assets")
output_dir <- file.path(book_dir, "_book")

required <- c("bookdown", "rmarkdown", "knitr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Missing required package(s): ", paste(missing, collapse = ", "),
    "\nInstall with: install.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "), "))"
  )
}
if (!rmarkdown::pandoc_available()) {
  stop("Pandoc is required. Install RStudio or Pandoc and retry.")
}

# Always refresh intermediate book chapters from the ignored report sources.
source(file.path(book_dir, "prepare_book.R"), local = new.env(parent = globalenv()), chdir = TRUE)

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(book_dir)

if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE, force = TRUE)

message("Rendering bookdown site ...")
bookdown::render_book(
  input = "index.Rmd",
  output_format = "bookdown::gitbook",
  envir = new.env(parent = globalenv()),
  clean = TRUE
)

copy_tree <- function(src, dst) {
  src <- normalizePath(src, mustWork = TRUE)
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(src, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  if (!length(entries)) return(invisible(TRUE))
  for (p in entries) {
    rel <- substring(p, nchar(src) + 2L)
    target <- file.path(dst, rel)
    if (dir.exists(p)) {
      dir.create(target, recursive = TRUE, showWarnings = FALSE)
    } else {
      dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(p, target, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)) {
        stop("Failed to copy book resource: ", p)
      }
    }
  }
  invisible(TRUE)
}

# Ensure chapter figures and other local assets are present at exactly the URLs
# written into the generated chapter Markdown.
if (dir.exists(assets_dir)) {
  copy_tree(assets_dir, file.path(output_dir, "book_assets"))
}

# GitHub Pages should serve the generated files verbatim and not run Jekyll.
file.create(file.path(output_dir, ".nojekyll"))

if (!file.exists(file.path(output_dir, "index.html"))) {
  stop("Book render completed without producing _book/index.html")
}

html_pages <- list.files(output_dir, pattern = "\\.html$", recursive = FALSE)
message("Book complete: ", output_dir)
message("HTML pages: ", length(html_pages))
