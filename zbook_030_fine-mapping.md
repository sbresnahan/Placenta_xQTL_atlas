# Cross-ancestry fine-mapping {#fine-mapping}

## Overview

SuSHiE jointly fine-maps each significant xQTL locus across ancestry strata, leveraging ancestry-specific LD to refine 95% credible sets and estimate cross-ancestry effect correlations. This report summarizes run diagnostics, credible-set properties, cross-ancestry concordance, and functional annotation enrichment of high-PIP variants.

| Metric                                           |   Value |
|:-------------------------------------------------|--------:|
| Loci attempted                                   |    5529 |
| Loci fine-mapped                                 |    5269 |
| Converged (min-tolerance)                        |    5275 |
| Loci with \>=1 credible set                      |    4290 |
| Credible sets total                              |    6347 |
| Variants fine-mapped (unique locus-variant rows) | 3382937 |



## PIP diagnostics


<img src="book_assets/fine-mapping/generated/html_fallback/fine-mapping/media/a9179aed5ff6423b0f13b5f51b3753b568b3a641.png" width="1200" alt="Maximum joint PIP per locus." />

Maximum joint PIP per locus.



<img src="book_assets/fine-mapping/generated/html_fallback/fine-mapping/media/08a6db7c67cbd19dc2f7729d6eff41b8d4825ca2.png" width="1200" alt="Distribution of joint PIPs across all fine-mapped variants (log scale)." />

Distribution of joint PIPs across all fine-mapped variants (log scale).




## Credible sets


<img src="book_assets/fine-mapping/generated/html_fallback/fine-mapping/media/70b89648f7d81ea5980845492e4fdedbfe91117f.png" width="1200" alt="Credible-set size distribution (number of variants per 95% CS)." />

Credible-set size distribution (number of variants per 95% CS).



<img src="book_assets/fine-mapping/generated/html_fallback/fine-mapping/media/f9d7968dc06605186381f156449b8e7b5851a8b0.png" width="1200" alt="Credible sets per locus." />

Credible sets per locus.




## Cross-ancestry concordance


<img src="book_assets/fine-mapping/generated/html_fallback/fine-mapping/media/02b2361f5ebc1764bb8f3a11ff0e99b2586a1c66.png" width="1200" alt="Cross-ancestry effect correlation (rho) per credible set. Values near 1 indicate shared causal effects; low/negative values indicate ancestry-specific effects." />

Cross-ancestry effect correlation (rho) per credible set. Values near 1 indicate shared causal effects; low/negative values indicate ancestry-specific effects.



<img src="book_assets/fine-mapping/generated/html_fallback/fine-mapping/media/8f1c0f4c321e6a8a5cc4b2c1ce6cadbb447679f1.png" width="1200" alt="Per-ancestry posterior effect weights for credible-set variants. Concordant effects lie along the diagonal." />

Per-ancestry posterior effect weights for credible-set variants. Concordant effects lie along the diagonal.




## Top loci


<img src="book_assets/fine-mapping/generated/html_fallback/fine-mapping/media/1266b810c0089a10d9507cb9a9e52232c5e50bda.png" width="1200" alt="Joint PIP by position for the six loci with the highest lead-variant PIP. Points in a credible set are colored by CS index." />

Joint PIP by position for the six loci with the highest lead-variant PIP. Points in a credible set are colored by CS index.




## Functional-annotation enrichment

High-PIP variants (PIP \>= 0.9) tested for enrichment against fastVEP consequence classes, ENCODE SCREEN cCRE classes, and placenta open-chromatin peaks; background = all other fine-mapped variants. Odds ratios and 95% CIs are recomputed from the 2x2 cell counts; FDR is Benjamini-Hochberg across classes.


<img src="book_assets/fine-mapping/generated/html_fallback/fine-mapping/media/2b7aca491ee8aecf53ef9c6ec6461d1f468bf3ae.png" width="1200" alt="Annotation enrichment of high-PIP variants (Fisher exact odds ratio, 95% CI)." />

Annotation enrichment of high-PIP variants (Fisher exact odds ratio, 95% CI).




## Run diagnostics

| Status            | Loci |
|:------------------|-----:|
| failed            |  254 |
| ok                | 5269 |
| vcf_export_failed |    6 |

<img src="book_assets/fine-mapping/generated/html_fallback/fine-mapping/media/761ea82d1dcf0ddfffe63d7c989e1f079fc0352b.png" width="1200" />



