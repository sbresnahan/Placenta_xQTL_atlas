# reports — analysis summary reports

## Placenta xQTL report: GTEx-conventions mapping, EAS + EUR

- 📄 `report_placenta_xqtl.html` — rendered, self-contained report
  (open in a browser; all figures embedded)
- 📜 `report_placenta_xqtl.Rmd` — report source
- 📜 `xqtl_report_functions.R` — helper module sourced by the Rmd (result
  loaders, QQ/power/concordance helpers, report theme)
- 📄 `gene_map.tsv` — gene_id → symbol/biotype/description map used for gene
  annotation in the report

Headline result: two FDR ≤ 5% associations across 32 ancestry × modality ×
layer cells — an EAS isoform xQTL at *TSPAN3* (q = 0.021 ungrouped, 0.022
grouped) and a EUR RNA-editing xQTL at *NCOA4* (ungrouped q = 0.029) — with
clean calibration and a power-limited null elsewhere. See
[`../docs/round_history.md`](../docs/round_history.md) for design context.

## Re-rendering

The report reads Stage-5 mapping outputs (top tables + QC files). Point the
params at a results directory and render:

```r
rmarkdown::render(
  "report_placenta_xqtl.Rmd",
  params = list(
    data_dir         = "<dir with *_cisqtl_top.tsv>",
    qc_dir           = "<dir with QC files>",
    module_path      = "xqtl_report_functions.R",
    gene_map_path    = "gene_map.tsv",
    gene_bodies_path = "<gene_bodies.tsv>",
    fig_dir          = "<fig out>", table_dir = "<table out>"
  )
)
```

Requires R ≥ 4.3 with rmarkdown, knitr, and tidyverse packages.

## Validated reproduction (2026-09-21)

The report was re-rendered end-to-end (55/55 chunks) from the current EAS +
EUR results using the shipped Rmd + module + gene map. The re-render
reproduces both FDR hits exactly (EAS isoforms at *TSPAN3*, lead
15:77800023:G:A, q = 0.021; EUR RNA editing at *NCOA4*, lead
10:45322442:G:A, q = 0.029) and all embedded figures. One practical note:

1. **`gene_bodies.tsv` is an explicit parameter** (`gene_bodies_path`). The
   file needs columns `gene_id, start, end, strand` (read as
   `col_types = "cddiic"`; strand as integer ±1) and is used to compute
   driver-variant position relative to the gene body. It can be regenerated
   from GENCODE v45 gene rows (strip version suffixes from `gene_id`); note a
   GENCODE-derived file will not cover HPLRv2-novel genes (they drop out of
   that one density figure only).
