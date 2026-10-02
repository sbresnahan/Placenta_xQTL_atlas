# reports

Analysis summary reports.

## Contents

- `report_placenta_xqtl.html` — rendered, self-contained main report (open
  in a browser; all figures embedded)
- `report_placenta_xqtl.Rmd` — main report source
- `xqtl_report_functions.R` — helper module sourced by the Rmd (result
  loaders, QQ/power/concordance helpers, report theme)
- `gene_map.tsv` — gene_id -> symbol/biotype/description map used for gene
  annotation in the report
- `report_finemap.Rmd` — SuSHiE fine-mapping diagnostic report (render
  command in ../05_qtl_mapping/README.md, section 6.5)
- `make_report_archive.sh` — builds the report-input tarball
  (`placenta_xqtl_report_inputs_<date>.tar.gz`) from mapping outputs

## Current results

GTEx-conventions mapping is complete for the EAS (n = 280) and EUR
(n = 140) ancestry strata across 9 RNA modalities: 12,482 FDR <= 5%
associations across the ancestry x modality x layer cells, and a combined
cross-modality catalog of 1,381 EAS and 1,602 EUR eGenes, with the
ERAP2/ERAP1 preeclampsia locus among the strongest in both ancestries.
Splicing is the one empirically null modality (Storey pi1 = 0, calibrated
QQ). AFR/AMR/SAS pooled genotypes exist and are the next mapping targets.

The EAS/EUR mapping was run without collapsing technical replicates: both
NIGMS runs per individual entered normalization and the QTL inputs, and duplicate
individual columns were averaged downstream (scripts 23/25/26). Run-level
and individual-level sample counts therefore differ in the report (e.g.
EUR 145 individuals at normalization vs 161 runs at QTL input). The optional
collapse pre-step is documented in ../05_qtl_mapping/README.md.

## Re-rendering

The report reads stage-5 mapping outputs (top tables + QC files). Point
the params at a results directory and render:

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

Requires R >= 4.3 with rmarkdown, knitr, and tidyverse packages.

`gene_bodies.tsv` is an explicit parameter (`gene_bodies_path`). The file
needs columns `gene_id, start, end, strand` (read as
`col_types = "cddiic"`; strand as integer +/-1) and is used to compute
driver-variant position relative to the gene body. It can be regenerated
from GENCODE v45 gene rows (strip version suffixes from `gene_id`); a
GENCODE-derived file will not cover HPLRv2-novel genes (they drop out of
that one density figure only).
