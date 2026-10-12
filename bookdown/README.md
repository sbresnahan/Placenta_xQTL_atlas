# Placenta xQTL Atlas bookdown site

This directory is the **tracked website builder** for the Placenta xQTL Atlas.
The scientific report directory (`/reports/`) remains intentionally ignored by
Git. The website is rebuilt locally from the report **source files**, then the
static output is committed to a dedicated GitHub Pages branch.

The book does **not** iframe or otherwise embed the existing standalone report
HTML files. For each analysis report, `prepare_book.R` re-runs the corresponding
parameterized `.Rmd` against the local report inputs, renders an intermediate
Markdown representation, normalizes it into a book chapter, and then bookdown
assembles the chapters into one navigable site. Literature-review `.md` files
are incorporated directly as native chapters.

## Included chapters

`report_manifest.tsv` is the authoritative chapter manifest. It currently
includes:

1. Multi-ancestry xQTL mapping
2. xQTL mapping literature context
3. Cross-ancestry fine-mapping
4. Fine-mapping literature context
5. GWAS colocalization
6. Colocalization literature context
7. TWAS
8. TWAS literature context
9. Cell-type attribution
10. Genotype × pre-pregnancy BMI (G×ppBMI)
11. G×ppBMI literature context
12. Integrated locus prioritization

The expected local integration source is
`reports/report_integration.Rmd`.

## Requirements

You need a local checkout that contains the ignored `reports/` directory and all
inputs needed by the individual report Rmd files. You also need R, Pandoc, the
analysis packages required by the reports, and these website packages:

```r
install.packages(c("bookdown", "rmarkdown", "knitr"))
```

Using the same R environment that already renders the standalone reports is the
safest choice because the report dependencies are deliberately not duplicated
here.

## Build locally

From the repository root:

```bash
./bookdown/build_bookdown.sh
```

or:

```bash
make -C bookdown build
```

The static site is written to:

```text
bookdown/_book/
```

The generated `_book/`, `.staging/`, `book_assets/`, and `zbook_*.Rmd` files are ignored. The builder
itself is what should be committed on the analysis branch.

## Publish to GitHub Pages

The publication helper builds the book, creates or updates a `gh-pages`
worktree, replaces that branch with the static site, commits the result, and
pushes it to `origin`:

```bash
./bookdown/publish_gh_pages.sh
```

To use another branch:

```bash
./bookdown/publish_gh_pages.sh website
```

To create the local publication commit without pushing it:

```bash
PUBLISH_NO_PUSH=1 ./bookdown/publish_gh_pages.sh
```

If the Pages branch already contains a `CNAME`, the publishing script preserves
it. A `.nojekyll` file is always added so GitHub Pages serves bookdown assets
verbatim.

In the GitHub repository settings, configure **Pages → Deploy from a branch**
and select:

```text
Branch: gh-pages
Folder: /(root)
```

GitHub Pages serves the files already rendered locally; it does not need to run
R or access the ignored reports/results.

## Updating the book

After regenerating any scientific report locally:

```bash
./bookdown/publish_gh_pages.sh
```

The builder re-executes every `.Rmd` in the manifest. This makes the website a
fresh snapshot of the complete report suite rather than a mixture of stale
standalone HTML artifacts.

To add another report, add one row to `report_manifest.tsv`, choose `kind=rmd`
or `kind=md`, then add the generated chapter filename to `_bookdown.yml`. The
chapter filename is deterministic:

```text
(part_order - 1) * 100 + chapter_order
```

formatted with three digits, followed by the chapter `id`. For example,
part 3 / chapter 20 / `gxe-ppbmi` becomes:

```text
.staging/220_gxe-ppbmi.Rmd
```

## Important distinction

There are two kinds of Git history involved:

- **analysis branch (for example `main`)**: tracks `bookdown/`, but not
  `reports/` or generated website files;
- **publication branch (`gh-pages`)**: tracks only the rendered static website.

This keeps controlled/local report inputs out of Git while still making the
complete report narrative publicly browsable.
