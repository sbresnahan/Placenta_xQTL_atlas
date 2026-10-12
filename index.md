---
title: "Placenta xQTL Atlas"
subtitle: "Multi-ancestry regulatory genomics, fine-mapping, GWAS integration, cell-type attribution, and G×E"
author: "Bhattacharya Lab"
date: "2026-10-12"
site: bookdown::bookdown_site
documentclass: book
github-repo: sbresnahan/Placenta_xQTL_atlas
description: "Integrated documentation and results for the Placenta xQTL Atlas."
---

# About this atlas {-}

This website is the integrated results book for the **Placenta xQTL Atlas**. It is assembled locally from the analysis-report sources in the repository's ignored `reports/` directory, then rendered as a static bookdown site for publication with GitHub Pages.

The book contains the complete scientific report sequence rather than embedding pre-rendered report HTML files. Parameterized `.Rmd` analyses are re-executed locally against the report inputs, converted to native book chapters, and combined with the companion literature-review `.md` files. The final synthesis chapter integrates xQTL mapping, fine-mapping, GWAS colocalization, TWAS, cell-type attribution, and genotype × pre-pregnancy BMI evidence.

## Reading order {-}

The book is organized in four parts: regulatory architecture, GWAS integration, biological context, and integrated synthesis. The analytical report for each layer is followed by its literature context when a companion review exists.

## Reproducibility {-}

The published website is a static artifact. The source analyses, controlled-access inputs, and intermediate report products are **not** required by GitHub Pages. To reproduce the website, build it locally from a working checkout containing the ignored `reports/` directory and its inputs, then publish the generated `_book` directory to the Pages branch.
