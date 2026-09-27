# reports — analysis summary reports

## Placenta xQTL report: GTEx-conventions mapping, EAS + EUR

- 📄 `report_placenta_xqtl.html` — rendered, self-contained report
  (open in a browser; all figures embedded)
- 📜 `report_placenta_xqtl.Rmd` — report source
- 📜 `xqtl_report_functions.R` — helper module sourced by the Rmd (result
  loaders, QQ/power/concordance helpers, report theme)
- 📄 `gene_map.tsv` — gene_id → symbol/biotype/description map used for gene
  annotation in the report

Headline result: 12,482 FDR ≤ 5% associations across the ancestry ×
modality × layer cells; the combined cross-modality eGene catalog carries
1,381 EAS and 1,602 EUR genes, with the *ERAP2*/*ERAP1* locus among the
strongest in both ancestries. Splicing is the one empirically null modality
(Storey π1 = 0, calibrated QQ).

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

## Validated reproduction (2026-09-27)

The report was re-rendered end-to-end (83/83 chunks) from the current EAS +
EUR results (9 modalities, all layers) using the shipped Rmd + module +
gene map. The re-render reproduces the headline catalog (12,482 FDR ≤ 5%
associations; 1,381 EAS / 1,602 EUR combined-layer eGenes; *ERAP2*/*ERAP1*
among the top loci in both ancestries) and all embedded figures. One
practical note:

1. **`gene_bodies.tsv` is an explicit parameter** (`gene_bodies_path`). The
   file needs columns `gene_id, start, end, strand` (read as
   `col_types = "cddiic"`; strand as integer ±1) and is used to compute
   driver-variant position relative to the gene body. It can be regenerated
   from GENCODE v45 gene rows (strip version suffixes from `gene_id`); note a
   GENCODE-derived file will not cover HPLRv2-novel genes (they drop out of
   that one density figure only).
