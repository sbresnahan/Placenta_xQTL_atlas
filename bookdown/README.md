# bookdown site — Placenta xQTL Atlas

This directory is the [bookdown](https://bookdown.org) source for the project
website. The two analysis reports and the two literature companions are
rendered as **native chapters** of one gitbook site (no embedded iframes, no
separate HTML viewers).

## Structure

| File | Role |
|---|---|
| `index.Rmd` | Homepage (`docs/index.html`) with document cards |
| `01-report1-intro.Rmd` | Part divider + title page for Report 1 |
| `02-report1-body.Rmd` | Report 1 source (xQTL mapping; identical to `reports/report_xqtl_mapping.Rmd`) |
| `03-report2-intro.Rmd` | Part divider + title page for Report 2 |
| `04-report2-body.Rmd` | Report 2 source (fine-mapping; identical to `reports/report_finemap.Rmd`) |
| `05-literature-xqtl.md` | Literature companion to Report 1 (unnumbered chapters, manual section numbers preserved) |
| `06-literature-finemap.md` | Literature companion to Report 2 |
| `_bookdown.yml` | Book config: `new_session: true` (each report knits in its own R session so its YAML `params` are honored), explicit `rmd_files` order |
| `_output.yml` | gitbook format config |
| `style.css` | Site theme (Phylo palette, scrollable tables) |
| `gene_map.tsv`, `xqtl_report_functions.R` | Helper inputs sourced by both reports |

Each report H1 becomes its own chapter page; Report 1 sections keep their
original numbering (chapters 1–11), Report 2 chapters continue the sequence,
and the literature companions are unnumbered (their manual section numbers
are preserved as written).

## Building

Input data are **not** in the repository. Regenerate the report inputs on the
cluster with `reports/make_report_archive.sh`, then make them available here
as `bookdown/data/` (symlink or copy), alongside `gene_map.tsv`:

```bash
ln -s /path/to/report_inputs/data bookdown/data
```

Then render:

```r
# from the bookdown/ directory, with bookdown installed
bookdown::render_book("index.Rmd", "bookdown::gitbook")
```

Output is written to `../docs/` (per `_bookdown.yml`), ready for GitHub
Pages. Knit side outputs (`figures/`, `tables/`, `*_files/`) and `data/` are
gitignored.

## Publishing

Commit `bookdown/` and `docs/`, push, then in GitHub:
**Settings → Pages → Build and deployment → Deploy from a branch → `main` /
`/docs` folder**. The site is served at
`https://sbresnahan.github.io/Placenta_xQTL_atlas/` with `docs/index.html` as
the homepage. `docs/.nojekyll` is included so Pages serves the bookdown
`libs/` and figure directories without Jekyll processing.

Note: the rendered figure PNGs in `docs/*_files/figure-html/` are required by
the chapter pages — commit `docs/` in full.
