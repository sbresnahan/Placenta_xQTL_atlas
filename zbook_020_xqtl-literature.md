# Literature context: placental xQTL mapping {#xqtl-literature}


**Project:** Placenta xQTL atlas (multi-ancestry, multi-modality cis-QTL mapping)  
**Revision date:** 2026-10-10  
**Companion analysis:** `report_placenta_xqtl.html`  
**Current result tables:** `tables/top_ungrouped_hits.tsv`, `tables/xqtl_catalog_phenotype_level.tsv`, `tables/phenotype_feature_counts.tsv`, `tables/egene_catalog_combined.tsv`, and enrichment tables.

**Citation convention.** Legacy bracketed citations `[N]` retain the indices from the original `xqtl_top_hits_evidence_table.csv`; their wording is preserved only where the corresponding gene remains in the current top-hit set. New literature or reference annotations added during this review are labeled `[322]`–`[339]` and are listed in the additional evidence (§9). This avoids renumbering the original evidence table.

---

### 1. Background and scope

The atlas maps cis genetic regulation in placenta across gene expression plus eight multi-phenotype RNA-processing layers: isoform expression, isoform usage, splicing, intron retention, alternative transcription start sites (alt_TSS), alternative polyadenylation (alt_polyA), RNA editing, and transcript stability. The final mapped sample sets contain **281 EAS and 136 EUR individuals** after genotype intersection, outlier removal, and technical-replicate collapse.

The current mapping report contains **8,116 FDR-significant phenotype-level discovery associations** (3,606 EAS and 4,510 EUR), implicating 1,757 EAS and 2,243 EUR unique genes across the discovery layer. At the combined gene-level layer there are 1,404 EAS and 1,654 EUR significant rows (1,403 and 1,650 unique gene IDs, respectively); 572 unique eGenes are shared across ancestries.

This literature review focuses on the **eight ungrouped multi-phenotype scans**. For each modality and ancestry, the five lowest-q entries from the current `top_ungrouped_hits.tsv` are reviewed, yielding **80 top-table entries representing 56 unique gene symbols**. A feature is labeled **shared** here when that same phenotype ID is present at q ≤ 0.05 in the opposite-ancestry phenotype-level catalog. By that definition, **39 of 80 top-table entries are shared**.

The purpose is descriptive rather than causal: to ask whether the strongest placental regulatory loci intersect known placental biology, prior QTL findings, or pregnancy-disease genetics, and to flag novel signals for mechanistic follow-up.

### 3. Methods

- **Hit selection.** Within each ancestry × ungrouped modality table, entries were ordered by Storey q-value and nominal p-value, and the top five were retained. Expression and the cross-modality combined scan are not part of this top-80 review because they do not have the same ungrouped multi-feature structure.
- **Cross-ancestry sharing.** “Shared” means the exact `phenotype_id` is FDR-significant (q ≤ 0.05) in both ancestries. Lead variants may differ. This is replication of a molecular feature, not proof of a shared causal variant and not a colocalization result.
- **Atlas-wide counts.** Significant-feature counts come from `xqtl_catalog_phenotype_level.tsv`; tested feature counts use the normalized (VST/INT) stage in `phenotype_feature_counts.tsv`.
- **Literature synthesis.** existing literature sources and numeric references were retained for genes still present. targeted literature searches were performed for newly entering genes, prioritizing PubMed papers with direct placental/pregnancy relevance, followed by mechanistic primary literature and curated reference resources when direct pregnancy literature was absent.
- **Opposite-ancestry q-values.** The supplied current tables include the significant catalog and top-10 tables, but not every non-significant ungrouped row. Therefore the tables report a binary shared flag rather than inventing or carrying forward stale `q(other)` values.

### 4. Overview of the current top hits


| Modality | EAS sig./tested | EUR sig./tested | Shared among top 5×2 |
|---|---:|---:|---:|
| Isoform expression | 1,333 / 37,217 | 1,603 / 41,150 | 8 / 10 |
| Isoform usage | 364 / 21,063 | 251 / 22,756 | 5 / 10 |
| Splicing | 5 / 56,258 | 28 / 79,006 | 0 / 10 |
| Intron retention | 56 / 3,818 | 102 / 10,765 | 0 / 10 |
| Alternative TSS | 184 / 29,868 | 90 / 37,279 | 5 / 10 |
| Alternative polyA | 371 / 26,322 | 459 / 33,027 | 10 / 10 |
| RNA editing | 30 / 3,120 | 246 / 17,951 | 7 / 10 |
| Transcript stability | 140 / 4,501 | 203 / 7,687 | 4 / 10 |

Alternative polyadenylation remains strikingly cross-ancestry reproducible at the top of the distribution (10/10). Isoform expression (8/10), RNA editing (7/10), isoform usage and alt-TSS (5/10 each) also contain substantial shared signal. By contrast, none of the selected splicing or intron-retention entries is shared at the exact-feature level.

#### 4.1 Recurring multi-modality loci

Several loci recur across independent molecular layers:

- **ERAP1/ERAP2 (chr5q15).** ERAP2 variant 5:96943786:A:G is the top EUR isoform-expression hit and drives two EUR promoter groups in opposite directions, including an HPLRv2 promoter. ERAP1 variant 5:96773906:G:A dominates both EUR isoform usage and alternative polyadenylation. This remains the atlas's strongest convergence with established preeclampsia genetics.
- **PSCA (chr8q24).** Variant 8:142680513:C:T drives opposite effects on two PSCA isoforms and the corresponding polyA features; both polyA features and both isoform-usage features are cross-ancestry significant.
- **ERGIC3 (chr20q13).** Variant 20:35549038:C:T is a top EAS isoform-expression, intron-retention, and alt-TSS hit; the isoform-expression and alt-TSS features are shared across ancestries.
- **MRPL43 (chr10q24).** Variant 10:100963159:G:C is the top EAS intron-retention hit and a top EAS alt-polyA hit.
- **PECAM1 (chr17q23).** Variant 17:64322640:G:A appears in EUR isoform expression, isoform usage, and alternative polyadenylation, including an HPLRv2 polyA feature.
- **GBP3 (chr1p22).** Variant 1:89009135:A:G is a top EUR isoform-expression and transcript-stability signal, linking transcript abundance to post-transcriptional stability at an interferon-responsive locus.
- **PSG9 (chr19q13).** Variant 19:43303030:G:T is a top EAS intron-retention and stability signal in a placenta-specific pregnancy glycoprotein.
- **MYL12B (chr18p11).** Variant 18:3278477:A:G drives opposite effects on two EAS TSS groups, both shared across ancestries.
- **HSP90B1 (chr12q23).** Alternative-polyA features are strongly significant in both ancestries, with different lead variants in EAS and EUR.

#### 4.2 Switch-like transcript regulation

The clearest same-variant opposite-direction patterns are **PSCA** (isoform usage and polyA), **ERAP1** (isoform usage and polyA), **MYL12B** (alt-TSS), **ERAP2** (two EUR TSS groups), and **SCARB2** (two EUR TSS groups). These loci are attractive for mechanism follow-up because the genetic effect appears to redistribute transcript identity or end choice rather than simply scale total RNA abundance.

#### 4.3 Splicing: corrected interpretation after modality-matched HCP covariates

In the current top five, EAS splicing lead AFs are approximately 0.053–0.098 and EUR lead minor-allele frequencies are approximately 0.052–0.054 (the TSR1 lead has AF 0.948, hence MAF ≈ 0.052). EUR top hits are concentrated around a chromosome-19 signal in which 19:44433568:C:T is the lead for ZNF234, APOC1P1/APOC1, and ZNF180. EAS remains a weak layer overall: only five ungrouped phenotypes reach FDR ≤ 5%.

This is now best described as **ancestry-asymmetric, modest-frequency splicing regulation**, not as a set of extreme rare-variant associations. The concentration of several EUR rows under one lead and the shared CMYA5/TENT2 EAS splice cluster also reinforce the need for locus-level rather than naïve gene-by-gene interpretation.

#### 4.4 HPLRv2-specific features among the top 80

Five current top entries explicitly carry `HPLRT_*` identifiers: **TRIM5** and **DRAM2** isoforms, an **ERAP2** alt-TSS feature, an **HSP90B1** alt-polyA feature, and a **PECAM1** alt-polyA feature. These are concrete examples in which the long-read placental annotation exposes regulatory phenotypes absent from a standard transcript-only view.

The chr22 `ENSG00000290199` hit remains poorly characterized, but the current feature ID is `ENST00000703580.1`; therefore the statement that this transcript was necessarily “HPLRv2-only” has been removed.

#### 4.5 RNA-editing coverage asymmetry

The normalized RNA-editing feature universe is **17,951 in EUR versus 3,120 in EAS (~5.75×)**. Consequently, the 246 versus 30 significant editing phenotypes cannot be interpreted as a simple ancestry biology contrast. Shared sites such as CTNND1 and HSPD1 are particularly valuable because they replicate despite the strong coverage imbalance.

#### 4.6 Atlas-wide ancestry sharing and functional context

At the combined gene-level layer, **572 of 1,403 EAS and 1,650 EUR unique eGenes are shared**. The current report estimates very strong enrichment for sharing relative to the tested-in-both universe (Fisher OR 6.44). Functional enrichment of combined and shared eGene sets is dominated by metabolic pathways, with “Metabolic pathways” the leading KEGG term for the EAS combined and shared sets. This broader metabolic signal complements—but should not be conflated with—the disease-focused literature enrichment apparent among the hand-selected top loci.

---

### 5. Modality-by-modality results and literature context


#### 5.1 Isoform expression

The strongest isoform-expression signals retain ERAP2/ZSCAN9 as positive controls but now bring **RBM23** into the EAS top five and **GBP3** into the EUR top five. Eight of the ten top-table entries correspond to features that are FDR-significant in both ancestries.


| Anc | Rank | Gene | Feature | Lead variant | Slope | AF | q | Shared |
|---|---|---|---|---|---|---|---|---|
| EAS | 1 | PGM1 | ENST00000371084.8 | 1:63616853:A:G | -0.55 | 0.129 | 5.19e-28 | yes |
| EAS | 2 | ERGIC3 | ENST00000489071.5 | 20:35549038:C:T | +1.25 | 0.133 | 1.15e-20 | yes |
| EAS | 3 | RBM23 | ENST00000346528.9 | 14:22918419:T:C | +1.49 | 0.785 | 1.99e-20 | no |
| EAS | 4 | ENSG00000290199 | ENST00000703580.1 | 22:23915150:G:A | +2.54 | 0.448 | 2.15e-19 | yes |
| EAS | 5 | FBLN1 | ENST00000262722.11 | 22:45483551:C:T | +0.50 | 0.554 | 3.20e-19 | no |
| EUR | 1 | ERAP2 | ENST00000437043.8 | 5:96943786:A:G | +1.72 | 0.415 | 8.19e-29 | yes |
| EUR | 2 | ZSCAN9 | ENST00000527436.5 | 6:28233360:A:G | +0.92 | 0.423 | 5.38e-24 | yes |
| EUR | 3 | GBP3 | ENST00000370481.9 | 1:89009135:A:G | -1.02 | 0.261 | 3.20e-23 | yes |
| EUR | 4 | PECAM1 | ENST00000563924.6 | 17:64322640:G:A | -0.63 | 0.552 | 4.69e-23 | yes |
| EUR | 5 | SOHLH2 | ENST00000379881.8 | 13:36218513:A:G | +0.89 | 0.452 | 1.03e-21 | yes |



- **PGM1** (EAS1). *Placental biology:* phosphoglucomutase isozymes were purified from human placenta decades ago [199], but the gene's modern literature is dominated by PGM1-CDG, a congenital glycosylation disorder treatable with galactose [200, 201]. No prior placental QTL or pregnancy-disease genetics — novel in this context.


- **ERGIC3** (EAS2). *Biology:* an ER–Golgi intermediate compartment protein; the ERGIC is a stress-responsive hub integrating membrane trafficking and adaptation [155], and ERGIC3-dependent secretory cargo trafficking is regulated by the E3 ligase MARCH2 [156]. No placental literature — a novel, strongly shared, multi-modality locus.


- **RBM23** (EAS3). *Biology:* RBM23 (CAPERβ) is a U2AF65-related RNA-binding protein that couples steroid-receptor transcriptional coactivation with alternative-splicing regulation [322]. I found no direct placental QTL or pregnancy-disease literature in the targeted literature search; the EAS isoform-expression QTL is therefore novel in placental context.


- **ENSG00000290199** (EAS4). This chr22 locus is currently represented by Ensembl transcript ENST00000703580.1 and remains poorly characterized. The results support a large-effect EAS isoform-expression QTL that is also FDR-significant in EUR. The literature review called the feature “HPLRv2-novel”; because the current feature ID is an Ensembl transcript rather than an HPLRT identifier, that novelty claim has been removed.


- **FBLN1** (EAS5). *Placental biology:* fibulin-1 is expressed during mouse placental development, with elevated transcripts in abnormal placentation [121]; in human pregnancy, plasma fibulin-1 rises across gestation and is consistently lower in PPROM [124]. No prior QTL evidence.


- **ERAP2** (EUR1). *Placental biology/disease:* ERAP2 is one of the most replicated preeclampsia candidate genes: associated with preeclampsia in Australian and Norwegian cohorts [110], fetal ERAP2 rs2549782 associated with preeclampsia in African Americans [112], a Chilean haplotype study links rs2549782 to ERAP2 protein expression and PE risk [113], and ERAP1/ERAP2 show altered placental expression and lymphocyte activation in PE [107]; ERAP2 was also associated with PE (and ERAP1 with eclampsia) in Brazilian women [108]. *Prior QTL evidence:* ERAP2 was one of the notable eGenes in the Kikas et al. placental eQTL study [56] — the present hit is a direct cross-ancestry replication and extension to isoform- and promoter-level regulation.


- **ZSCAN9** (EUR2). *Prior QTL evidence:* ZSCAN9 was likewise a notable placental eGene in Kikas et al. [56], making this a second direct replication. *Biology:* it is a KRAB zinc-finger protein; KRAB-ZFPs repress transposable elements and shape mammalian regulatory evolution [316, 317], a mechanism of particular relevance in placenta. The atlas continues to place multiple chromosome-19 zinc-finger loci among strong regulatory signals, although the specific top splicing genes have changed.


- **GBP3** (EUR3). *Biology:* guanylate-binding proteins are interferon-induced pattern-recognition/antimicrobial effectors [114, 115]. Novel in placenta; joins OAS1 and TRIM5 in the innate-immunity theme.
 The same variant, 1:89009135:A:G, is also the EUR3 transcript-stability hit, making GBP3 a multi-modality innate-immune locus.


- **PECAM1** (EUR4). *Placental biology:* CD31/PECAM-1 marks placental endothelium; microvessel density stained for CD31 is markedly reduced in preeclamptic placentas [222] and PECAM-1 expression is decreased in PE placentae [218], though one study found no difference in CD31 expression or vascular growth [223], and PECAM-1 itself does not appear to drive trophoblast invasion or PE/FGR pathophysiology [219]. The QTL is shared across ancestries and acts on both isoform abundance and polyA-site usage.


- **SOHLH2** (EUR5). *Biology:* a germ-cell-specific bHLH transcription factor essential for oocyte and spermatogonial differentiation [202, 203, 204]. No prior placental QTL or pregnancy-disease association — a novel placental regulatory signal, plausibly reflecting its role in germline/early-lineage transcription programs.



#### 5.2 Isoform usage

The PSCA and ERAP1 switch-like loci remain dominant. The EUR fifth-ranked entry is now **PECAM1**, replacing ZNF880 from the literature review.


| Anc | Rank | Gene | Feature | Lead variant | Slope | AF | q | Shared |
|---|---|---|---|---|---|---|---|---|
| EAS | 1 | PSCA | ENST00000513264.1 | 8:142680513:C:T | -0.79 | 0.407 | 4.10e-18 | yes |
| EAS | 2 | PSCA | ENST00000301258.5 | 8:142680513:C:T | +0.78 | 0.407 | 5.48e-18 | yes |
| EAS | 3 | OAS1 | ENST00000202917.10 | 12:112943944:T:C | -0.93 | 0.748 | 5.48e-18 | no |
| EAS | 4 | METTL2B | ENST00000480046.5 | 7:128477173:G:A | +0.89 | 0.296 | 1.51e-16 | no |
| EAS | 5 | PARP2 | ENST00000250416.9 | 14:20344308:C:A | +1.04 | 0.176 | 3.80e-16 | yes |
| EUR | 1 | ERAP1 | ENST00000296754.7 | 5:96773906:G:A | +1.28 | 0.287 | 1.97e-24 | no |
| EUR | 2 | ERAP1 | ENST00000443439.7 | 5:96773906:G:A | -1.27 | 0.287 | 1.01e-20 | no |
| EUR | 3 | TRIM5 | HPLRT_f3660beda5b | 11:5938806:G:A | -1.09 | 0.430 | 3.14e-15 | yes |
| EUR | 4 | DRAM2 | HPLRT_a8773b56079 | 1:111100631:TC:T | +1.23 | 0.280 | 4.34e-15 | yes |
| EUR | 5 | PECAM1 | ENST00000563924.6 | 17:64322640:G:A | -0.74 | 0.552 | 4.70e-15 | no |



- **PSCA** (EAS1–2). A textbook isoform-switch QTL: the same allele lowers one PSCA isoform and raises the other. *Literature:* PSCA's published record is entirely oncologic — expression rises with prostate cancer grade, stage, and metastasis [213], with context-dependent tumor-promoting or -suppressing roles [214] and prognostic/therapeutic interest across cancers [215]. There is no placental or pregnancy literature: this is a novel, cross-ancestry-shared placental isoform QTL for a GPI-anchored cell-surface antigen.


- **OAS1** (EAS3). *Placental biology:* OAS1 is an innate immune sensor of dsRNA [95]; the placenta maintains baseline type I interferon/ISG signaling that preserves pregnancy homeostasis and antiviral readiness [97], and hemochorial placentas constitutively express type III interferon in trophoblasts [98]. A biologically coherent EAS-specific isoform QTL in the placental antiviral pathway.


- **METTL2B** (EAS4). *Biology:* an m3C tRNA methyltransferase with mutually exclusive substrate selection against METTL6 [283]; m3C32 modification controls serine codon-biased translation, cell cycle, and DNA-damage response [286]. Novel in placenta.


- **PARP2** (EAS5). *Disease relevance:* PARP-family activity is mechanistically tied to placental dysfunction — NAD+ depletion from NAD+-consuming enzymes (including PARPs) marks an inflammatory subclass of preeclampsia and is reversible with nicotinamide riboside in models [271, 272], and PARP activation is elevated in gestational diabetic pregnancies [268]. *Biology:* PARP2 has distinct and shared functions with PARP1 in inflammation, metabolism, and oxidative stress [270].


- **ERAP1** (EUR1–2). An isoform-switch QTL mirroring PSCA. *Disease:* ERAP1 expression is elevated in preeclamptic placentas and inducible by hypoxia in trophoblasts [111]; ERAP1/ERAP2 placental expression and lymphocyte activation are altered in PE [107]; ERAP1 is associated with eclampsia in Brazilian women [108]. *Prior QTL evidence:* the adjacent ERAP2 was a notable placental eGene [56]; the present data resolve the locus into isoform- and polyA-level ERAP1 regulation.


- **TRIM5** (EUR3). *Biology:* TRIM5α is a canonical retroviral restriction factor that recognizes and destabilizes incoming capsids [173], within the broader TRIM antiviral family [170]. No placental literature — novel; given the placenta's constitutive antiviral state (see OAS1 below), a plausible innate-immunity locus.


- **DRAM2** (EUR4). *Biology:* a lysosomal autophagy regulator; biallelic mutations cause retinal dystrophy [276], and despite divergence from DRAM it contributes to autophagy induction [277, 278]. Novel in placenta, where autophagy is central to trophoblast homeostasis.


- **PECAM1** (EUR5). *Placental biology:* CD31/PECAM-1 marks placental endothelium; microvessel density stained for CD31 is markedly reduced in preeclamptic placentas [222] and PECAM-1 expression is decreased in PE placentae [218], though one study found no difference in CD31 expression or vascular growth [223], and PECAM-1 itself does not appear to drive trophoblast invasion or PE/FGR pathophysiology [219]. The QTL is shared across ancestries and acts on both isoform abundance and polyA-site usage.
 This same PECAM1 transcript is also among the top EUR isoform-expression hits, while a novel HPLRv2 polyA feature at PECAM1 is EUR4 in alt_polyA.



#### 5.3 Splicing

Under modality-matched HCP covariates, EUR splicing is no longer null and EAS is weak but nonzero. The old extreme 1% minor-allele-frequency LILRB/CNN2/DOCK8/SREK1/ZNF439/440 top set is not present in the results.


| Anc | Rank | Gene | Feature | Lead variant | Slope | AF | q | Shared |
|---|---|---|---|---|---|---|---|---|
| EAS | 1 | ZFP62 | chr5_180850329_180850582_clu_18215_- | 5:180193795:T:G | +0.41 | 0.086 | 0.0054 | no |
| EAS | 2 | CMYA5 | chr5_79669028_79679579_clu_17479_+ | 5:80551240:G:C | +0.18 | 0.063 | 0.0054 | no |
| EAS | 3 | TENT2 | chr5_79669028_79679579_clu_17479_+ | 5:80551240:G:C | +0.18 | 0.063 | 0.0054 | no |
| EAS | 4 | PLAC4 | chr21_41182302_41186123_clu_14257_- | 21:41210825:A:G | +0.30 | 0.053 | 0.0192 | no |
| EAS | 5 | TRMT13 | chr1_100144143_100146527_clu_1689_+ | 1:99717801:C:A | +0.42 | 0.098 | 0.0192 | no |
| EUR | 1 | ZNF234 | chr19_44141733_44142171_clu_11815_+ | 19:44433568:C:T | +0.44 | 0.054 | 3.81e-06 | no |
| EUR | 2 | APOC1P1 | chr19_44916325_44931000_clu_11828_+ | 19:44433568:C:T | +0.50 | 0.054 | 4.41e-05 | no |
| EUR | 3 | APOC1 | chr19_44916325_44931000_clu_11828_+ | 19:44433568:C:T | +0.50 | 0.054 | 4.41e-05 | no |
| EUR | 4 | ZNF180 | chr19_44485157_44497284_clu_11058_- | 19:44433568:C:T | +0.45 | 0.054 | 0.001 | no |
| EUR | 5 | TSR1 | chr17_2334896_2335256_clu_9525_- | 17:3024477:G:C | -1.06 | 0.948 | 0.0018 | no |



- **ZFP62** (EAS1). *Biology:* ZFP62 is a poorly characterized human zinc-finger protein [333]. No direct placental or pregnancy literature was identified in the targeted literature search, so the EAS splicing association is currently best treated as an uncharacterized regulatory signal.


- **CMYA5** (EAS2). *Biology:* CMYA5/myospryn is a very large striated-muscle protein induced by MEF2A, localized to the costamere, and interacting with α-actinin-2 [323]. Direct placental evidence is lacking; because the same splice cluster is also annotated to TENT2, gene attribution at this locus should remain cautious.
 CMYA5 and TENT2 are assigned to the same splice cluster and lead variant; gene attribution is therefore ambiguous.


- **TENT2** (EAS3). *Biology:* TENT2/GLD2 is a non-canonical terminal nucleotidyl transferase involved in RNA 3′-end metabolism, a class of enzymes that regulates RNA stability and translation [324]. No direct placental disease evidence was found; the shared CMYA5/TENT2 splice-cluster annotation makes locus-level interpretation preferable to single-gene attribution.


- **PLAC4** (EAS4). *Placental biology/disease:* PLAC4 is highly placenta-enriched, localizes to syncytiotrophoblast, increases with syncytialization, and was reported to be upregulated in severe early-onset preeclamptic placenta [325]. This makes the EAS PLAC4 splicing QTL one of the most directly placenta-relevant new hits created by the corrected splicing analysis.


- **TRMT13** (EAS5). *Biology:* TRMT13 is a predicted tRNA methyltransferase involved in tRNA methylation [331]. Human Protein Atlas placental immunohistochemistry reports medium signal in decidual cells and low signal in trophoblastic cells [332]. No pregnancy-disease genetic literature was identified, so this is a new RNA-processing candidate rather than a known disease locus.


- **ZNF234** (EUR1). *Biology:* ZNF234 is a chromosome-19 C2H2 zinc-finger protein predicted to regulate DNA-templated transcription [334]. Direct placental or pregnancy evidence was not identified. It participates in the same EUR chr19 splicing hotspot as APOC1/APOC1P1 and ZNF180.


- **APOC1P1** (EUR2). *Annotation:* APOC1P1 is a transcribed/pseudogene locus adjacent to APOC1 [335]. Its top splice cluster is identical to the APOC1-annotated cluster, so the signal should be discussed as an APOC1/APOC1P1 locus rather than as independent evidence for the pseudogene.


- **APOC1** (EUR3). *Placental biology/disease:* Recent placental single-cell/spatial multi-omics implicated APOC1 among central lipid-transport genes dysregulated in preeclampsia, with subtype- and cell-type-specific changes in EVT and Hofbauer cells [326]. A 2026 study further linked macrophage APOC1 deficiency to ferroptosis, inflammatory polarization, impaired trophoblast invasion, and recurrent pregnancy loss [327]. The EUR splicing signal therefore intersects a newly emerging maternal–fetal-interface literature.
 The APOC1 and APOC1P1 rows share the same splice cluster and lead variant, so these should not be read as independent mechanistic signals.


- **ZNF180** (EUR4). *Biology:* ZNF180 is a Kruppel-related C2H2 zinc-finger protein with broadly expressed transcriptional-regulatory potential [336]. No direct placental disease literature was identified; its association shares the same EUR chr19 lead variant as ZNF234 and APOC1/APOC1P1.


- **TSR1** (EUR5). *Biology:* TSR1 is a nucleolar RNA-binding ribosome-maturation factor implicated in pre-rRNA processing and small-subunit maturation [337]. No direct placental literature was identified, but the hit adds ribosome biogenesis to the RNA-processing signals represented among top xQTLs.



#### 5.4 Intron retention

The ancestry split remains strong: none of the ten top entries is FDR-significant in both ancestries. EAS is led by MRPL43/ERGIC3/PSG9; EUR by GM2A/IQGAP1/LGALS14/CYP19A1. The fifth EUR row is a second annotation of the CYP19A1 retention event.


| Anc | Rank | Gene | Feature | Lead variant | Slope | AF | q | Shared |
|---|---|---|---|---|---|---|---|---|
| EAS | 1 | MRPL43 | IR_ENSG00000055950_chr10_100986573_100986587_-_14340 | 10:100963159:G:C | -0.75 | 0.530 | 5.62e-16 | no |
| EAS | 2 | ERGIC3 | IR_ENSG00000125991_chr20_35543140_35543657_+_69075 | 20:35549038:C:T | +0.92 | 0.133 | 2.21e-13 | no |
| EAS | 3 | PSG9 | IR_ENSG00000183668_chr19_43269136_43269367_-_56775 | 19:43303030:G:T | +0.75 | 0.368 | 3.48e-13 | no |
| EAS | 4 | LGALS8 | IR_ENSG00000116977_chr1_236543041_236543559_+_10681 | 1:236542978:G:T | -0.84 | 0.151 | 2.14e-10 | no |
| EAS | 5 | VAMP8 | IR_ENSG00000118640_chr2_85581717_85581753_+_62360 | 2:85585623:T:C | -0.70 | 0.398 | 4.54e-10 | no |
| EUR | 1 | GM2A | IR_ENSG00000196743_chr5_151267146_151267295_+_71349 | 5:151262704:A:G | +1.07 | 0.261 | 2.40e-14 | no |
| EUR | 2 | IQGAP1 | IR_ENSG00000140575_chr15_90439400_90440501_+_30450 | 15:90359318:C:T | +1.07 | 0.211 | 1.52e-12 | no |
| EUR | 3 | LGALS14 | IR_ENSG00000006659_chr19_39707582_39709196_+_43912 | 19:39709274:C:G | +0.79 | 0.607 | 1.52e-12 | no |
| EUR | 4 | CYP19A1 | IR_ENSG00000137869_chr15_51255762_51256467_-_28695 | 15:51253257:A:G | -0.72 | 0.330 | 1.93e-12 | no |
| EUR | 5 | CYP19A1 | IR_ENSG00000137869_chr15_51255762_51256467_-_28695 | 15:51253257:A:G | -0.72 | 0.330 | 1.93e-12 | no |



- **MRPL43** (EAS1). *Prior QTL evidence:* a single 97 kb cis-eQTL spanning MRPL43 has antagonistic effects across regulatory stages — raising mRNA and ribosome occupancy while lowering protein — via linked variants affecting TF binding, splicing, and miRNA targeting [265]; the present intron-retention and polyA QTLs fit that multi-layer architecture. *Placental relevance:* nuclear-encoded mitochondrial ribosomal proteins are required to initiate gastrulation [262], and placental mitochondrial dysfunction is a recurring feature of preeclampsia, IUGR/FGR, and preterm birth [263, 264].


- **ERGIC3** (EAS2). *Biology:* an ER–Golgi intermediate compartment protein; the ERGIC is a stress-responsive hub integrating membrane trafficking and adaptation [155], and ERGIC3-dependent secretory cargo trafficking is regulated by the E3 ligase MARCH2 [156]. No placental literature — a novel, strongly shared, multi-modality locus.


- **PSG9** (EAS3). *Placental biology/disease:* pregnancy-specific glycoproteins are placenta-secreted immunomodulators [3]; PSG9 (with PSG7) rises before preeclampsia onset [4]; endogenous retroviral elements LTR8B/MER65 rewire PSG9 regulation to control trophoblast syncytialization and PE risk [5]; PSG9 induces FoxP3+ regulatory T cells via TGF-β1 [6] and enhances endothelial nitric-oxide release, with lower serum PSG9 in PE [7]. One of the most placenta-specific genes in the entire hit set.


- **LGALS8** (EAS4). *Placental biology:* galectins collectively shape early pregnancy and pregnancy pathologies [18]. Note that the well-studied chr19 placental galectin cluster (LGALS13/14/16) does not include LGALS8 (chr1); LGALS8-specific placental evidence is limited — largely novel.


- **VAMP8** (EAS5). *Biology:* a SNARE for regulated exocytosis, selective for granule-to-granule fusion [181] and general across exocrine secretion [182, 183]. No placental literature — novel; secretory machinery is nonetheless central to syncytiotrophoblast hormone release.


- **GM2A** (EUR1). *Biology:* the GM2-activator protein is a lysosomal lipid-transport cofactor for GM2 ganglioside hydrolysis [224, 226, 228]. Novel in placenta.


- **IQGAP1** (EUR2). *Placental biology:* IQGAP1 is a scaffold promoting cell motility and invasion via Cdc42/Rac1 [165, 168] and specifically promotes trophoblast motility and invasion [167] — directly relevant to extravillous trophoblast function.


- **LGALS14** (EUR3). *Placental biology:* galectin-14 is a placenta-specific, primate galectin that promotes trophoblast migration and invasion via MMP-9/N-cadherin [15], belongs to a chr19 cluster evolved at the maternal-fetal interface to induce maternal immune-cell death [16], modulates immune-cell survival and cytokines [19], and protects trophoblasts from oxidative stress while promoting extravillous differentiation [20]. *Disease:* the chr19 galectin cluster is complexly dysregulated in preeclampsia [17]. A high-priority placenta-specific hit.


- **CYP19A1** (EUR4–5). *Placental biology/disease:* aromatase drives placental estrogen synthesis; its downregulation is proposed to play a dual (causal and predictive) role in preeclampsia [43], placental estradiol is markedly decreased in PE [47], and oxygen-dependent CYP19A1 induction during trophoblast differentiation is mediated by ERRγ/HIF-1 [50], with tissue-specific promoter usage long established [48]. Intron retention is a plausible post-transcriptional layer on this heavily regulated gene.



#### 5.5 Alternative TSS

The promoter results replace the prior rare ACBD3/MMP14/RAF1/PRKCH-AS1/TEX10 top hits. Strong switch-like patterns now include MYL12B in EAS and ERAP2/SCARB2 in EUR.


| Anc | Rank | Gene | Feature | Lead variant | Slope | AF | q | Shared |
|---|---|---|---|---|---|---|---|---|
| EAS | 1 | ERGIC3 | grp_2_upstream_ENST00000489071.5 | 20:35549038:C:T | +0.92 | 0.133 | 5.94e-15 | yes |
| EAS | 2 | MYL12B | grp_1_upstream_ENST00000400175.9 | 18:3278477:A:G | +0.60 | 0.566 | 3.33e-12 | yes |
| EAS | 3 | MYL12B | grp_1_upstream_ENST00000237500.10 | 18:3278477:A:G | -0.51 | 0.566 | 3.33e-12 | yes |
| EAS | 4 | LGALS8 | grp_2_upstream_ENST00000489586.6 | 1:236537086:G:C | -0.74 | 0.152 | 2.06e-09 | no |
| EAS | 5 | RPL14 | grp_1_upstream_ENST00000396203.7 | 3:40453373:C:A | -0.87 | 0.119 | 3.72e-08 | no |
| EUR | 1 | ERAP2 | grp_2_upstream_HPLRT_039d234c78d | 5:96943786:A:G | -0.86 | 0.415 | 2.68e-16 | yes |
| EUR | 2 | SCARB2 | grp_1_upstream_ENST00000638295.1 | 4:76257226:A:G | +0.91 | 0.604 | 9.28e-08 | no |
| EUR | 3 | ENOSF1 | grp_2_upstream_ENST00000584259.6 | 18:691173:G:A | -0.49 | 0.548 | 1.53e-07 | yes |
| EUR | 4 | ERAP2 | grp_2_upstream_ENST00000437043.8 | 5:96943786:A:G | +0.82 | 0.415 | 2.40e-07 | no |
| EUR | 5 | SCARB2 | grp_1_upstream_ENST00000264896.8 | 4:76257226:A:G | -0.88 | 0.604 | 2.40e-07 | no |



- **ERGIC3** (EAS1). *Biology:* an ER–Golgi intermediate compartment protein; the ERGIC is a stress-responsive hub integrating membrane trafficking and adaptation [155], and ERGIC3-dependent secretory cargo trafficking is regulated by the E3 ligase MARCH2 [156]. No placental literature — a novel, strongly shared, multi-modality locus.


- **MYL12B** (EAS2–3). A promoter-switch QTL (opposite slopes on two TSS groups). *Placental biology:* myosin regulatory light chains maintain myosin II stability and cellular integrity [230]; human placental stem villi contain ultraslow myosin motors [233] whose contraction uses actin–myosin cross-bridges under NO–cGMP control [234]. A contractility-relevant gene with shared cross-ancestry promoter regulation.


- **LGALS8** (EAS4). *Placental biology:* galectins collectively shape early pregnancy and pregnancy pathologies [18]. Note that the well-studied chr19 placental galectin cluster (LGALS13/14/16) does not include LGALS8 (chr1); LGALS8-specific placental evidence is limited — largely novel.


- **RPL14** (EAS5). *Biology:* RPL14 encodes a large-subunit ribosomal protein [339]. The literature review already noted that mTORC1 controls ribosome biogenesis and protein synthesis in primary human trophoblasts [259]; the EAS alt-TSS QTL places a specific ribosomal gene under cis promoter regulation in placenta.


- **ERAP2** (EUR1–4). *Placental biology/disease:* ERAP2 is one of the most replicated preeclampsia candidate genes: associated with preeclampsia in Australian and Norwegian cohorts [110], fetal ERAP2 rs2549782 associated with preeclampsia in African Americans [112], a Chilean haplotype study links rs2549782 to ERAP2 protein expression and PE risk [113], and ERAP1/ERAP2 show altered placental expression and lymphocyte activation in PE [107]; ERAP2 was also associated with PE (and ERAP1 with eclampsia) in Brazilian women [108]. *Prior QTL evidence:* ERAP2 was one of the notable eGenes in the Kikas et al. placental eQTL study [56] — the present hit is a direct cross-ancestry replication and extension to isoform- and promoter-level regulation.
 The same 5:96943786:A:G allele drives opposite directions across two ERAP2 promoter groups in the EUR top five (novel HPLRv2 promoter negative; ENST00000437043.8 promoter positive), while the ERAP2 isoform-expression association is also top-ranked.


- **SCARB2** (EUR2–5). *Placental context:* SCARB2 is a lysosomal membrane protein. Human Protein Atlas placental staining is high in syncytiotrophoblast cell body/microvilli and Hofbauer cells and medium in endothelium [329]. I found no direct pregnancy-disease QTL literature, making the EUR promoter-switch signal a plausible but largely novel placental regulatory locus.
 The same 4:76257226:A:G allele has opposite effects on two SCARB2 TSS groups (+0.91 and −0.88), a promoter-switch pattern.


- **ENOSF1** (EUR3). *Biology:* ENOSF1 can encode a mitochondrial L-fuconate-pathway enzyme and also produces transcripts historically characterized as antisense regulators of adjacent TYMS [330]. No direct placental or pregnancy literature was identified. Its EUR alt-TSS feature is cross-ancestry significant despite limited prior biological context.



#### 5.6 Alternative polyA

Alternative polyadenylation remains the most cross-ancestry-robust top-hit modality: all ten top-table entries correspond to features significant in both ancestries.


| Anc | Rank | Gene | Feature | Lead variant | Slope | AF | q | Shared |
|---|---|---|---|---|---|---|---|---|
| EAS | 1 | PSCA | grp_2_downstream_ENST00000513264.1 | 8:142680513:C:T | -0.87 | 0.407 | 6.63e-23 | yes |
| EAS | 2 | PSCA | grp_2_downstream_ENST00000301258.5 | 8:142680513:C:T | +0.86 | 0.407 | 6.63e-23 | yes |
| EAS | 3 | ELP5 | grp_2_downstream_ENST00000574255.6 | 17:7254884:T:G | -0.89 | 0.655 | 2.89e-21 | yes |
| EAS | 4 | MRPL43 | grp_2_downstream_ENST00000299179.9 | 10:100963159:G:C | +0.81 | 0.530 | 3.28e-21 | yes |
| EAS | 5 | HSP90B1 | grp_1_downstream_HPLRT_270ba7eb207 | 12:103943444:GTAAT:G | +0.72 | 0.434 | 6.05e-21 | yes |
| EUR | 1 | ERAP1 | grp_1_downstream_ENST00000443439.7 | 5:96773906:G:A | -1.20 | 0.287 | 1.76e-35 | yes |
| EUR | 2 | ERAP1 | grp_1_downstream_ENST00000296754.7 | 5:96773906:G:A | +1.20 | 0.287 | 2.10e-34 | yes |
| EUR | 3 | HSP90B1 | grp_2_downstream_ENST00000550595.2 | 12:103946200:G:A | +0.79 | 0.342 | 1.37e-28 | yes |
| EUR | 4 | PECAM1 | grp_1_downstream_HPLRT_683224049f0 | 17:64322640:G:A | +0.87 | 0.552 | 3.96e-28 | yes |
| EUR | 5 | ERAP1 | grp_2_downstream_ENST00000443439.7 | 5:96773906:G:A | -1.29 | 0.287 | 5.18e-28 | yes |



- **PSCA** (EAS1–2). A textbook isoform-switch QTL: the same allele lowers one PSCA isoform and raises the other. *Literature:* PSCA's published record is entirely oncologic — expression rises with prostate cancer grade, stage, and metastasis [213], with context-dependent tumor-promoting or -suppressing roles [214] and prognostic/therapeutic interest across cancers [215]. There is no placental or pregnancy literature: this is a novel, cross-ancestry-shared placental isoform QTL for a GPI-anchored cell-surface antigen.


- **ELP5** (EAS3). *Biology:* Elongator-complex subunit; the complex's two subcomplexes have divergent neurodevelopmental roles [290, 291] and are linked to human neurological disease [292]. No placental literature — a novel, strongly shared (q ≈ 1e-20 in both ancestries) polyA QTL.


- **MRPL43** (EAS4). *Prior QTL evidence:* a single 97 kb cis-eQTL spanning MRPL43 has antagonistic effects across regulatory stages — raising mRNA and ribosome occupancy while lowering protein — via linked variants affecting TF binding, splicing, and miRNA targeting [265]; the present intron-retention and polyA QTLs fit that multi-layer architecture. *Placental relevance:* nuclear-encoded mitochondrial ribosomal proteins are required to initiate gastrulation [262], and placental mitochondrial dysfunction is a recurring feature of preeclampsia, IUGR/FGR, and preterm birth [263, 264].


- **HSP90B1** (EAS5, EUR3). *Placental biology:* gp96/GRP94 is expressed at the materno-fetal interface, most intensely in trophoblast and glandular epithelium [134], and modulates decidual antigen-presenting and NK cell activation [135]. *Disease:* unfolded-protein-response activation, including HSP90B1-related signaling, distinguishes early-onset from late-onset preeclampsia placentae [136].


- **ERAP1** (EUR1–2–5). An isoform-switch QTL mirroring PSCA. *Disease:* ERAP1 expression is elevated in preeclamptic placentas and inducible by hypoxia in trophoblasts [111]; ERAP1/ERAP2 placental expression and lymphocyte activation are altered in PE [107]; ERAP1 is associated with eclampsia in Brazilian women [108]. *Prior QTL evidence:* the adjacent ERAP2 was a notable placental eGene [56]; the present data resolve the locus into isoform- and polyA-level ERAP1 regulation.


- **PECAM1** (EUR4). *Placental biology:* CD31/PECAM-1 marks placental endothelium; microvessel density stained for CD31 is markedly reduced in preeclamptic placentas [222] and PECAM-1 expression is decreased in PE placentae [218], though one study found no difference in CD31 expression or vascular growth [223], and PECAM-1 itself does not appear to drive trophoblast invasion or PE/FGR pathophysiology [219]. The QTL is shared across ancestries and acts on both isoform abundance and polyA-site usage.



#### 5.7 RNA editing

The principal shared editing sites remain CTNND1 and HSPD1. **PGAP6** enters the EAS top five, replacing CNOT12 in the literature review.


| Anc | Rank | Gene | Feature | Lead variant | Slope | AF | q | Shared |
|---|---|---|---|---|---|---|---|---|
| EAS | 1 | ANXA1 | chr9_73160319 | 9:73163667:AG:A | +0.91 | 0.185 | 1.28e-12 | no |
| EAS | 2 | LCMT1 | chr16_25147778 | 16:25114437:T:C | +0.75 | 0.349 | 2.39e-12 | yes |
| EAS | 3 | HSPD1 | chr2_197497294 | 2:197446138:A:G | +0.71 | 0.541 | 7.10e-12 | yes |
| EAS | 4 | CTNND1 | chr11_57761957 | 11:57790196:C:T | +1.04 | 0.120 | 1.47e-11 | yes |
| EAS | 5 | PGAP6 | chr16_377516 | 16:377657:G:A | +0.75 | 0.308 | 2.67e-11 | yes |
| EUR | 1 | CTNND1 | chr11_57761957 | 11:57826368:G:A | -1.16 | 0.566 | 2.96e-26 | yes |
| EUR | 2 | CTSB | chr8_11844424 | 8:11844866:G:A | +1.43 | 0.257 | 7.28e-24 | no |
| EUR | 3 | RBBP4 | chr1_32680141 | 1:32669173:A:C | +1.42 | 0.752 | 1.85e-22 | no |
| EUR | 4 | HSPD1 | chr2_197497294 | 2:197400449:T:A | +1.25 | 0.701 | 2.11e-22 | yes |
| EUR | 5 | CTSB | chr8_11845033 | 8:11847833:G:C | -1.19 | 0.634 | 3.50e-19 | yes |



- **ANXA1** (EAS1). *Disease:* annexin A1 inhibits trophoblast ferroptosis in preeclampsia via KISS1 downregulation [71]; silencing ANXA1 suppresses apoptosis and inflammation in PE rat trophoblasts [72]; plasma ANXA1 is elevated in PE, especially early-onset [78]; placental ANXA1 influences DNA-damage response in gestational diabetes models [77]; the annexin family is broadly implicated at the maternal-fetal interface [75]. Rich disease literature; the editing QTL is EAS-specific.


- **LCMT1** (EAS2). *Biology:* leucine carboxyl methyltransferase-1 methylates PP4/PP6 and regulates phosphatase holoenzyme assembly [160]; global loss causes fetal-liver hematopoietic failure and embryonic lethality [158]; protein is expressed in placental trophoblasts [161]. Novel as an editing QTL.


- **HSPD1** (EAS3, EUR4). *Disease:* Hsp60 is overexpressed in preeclamptic placentas [152], circulating Hsp60 is elevated in PE [148], Hsp60 lactylation promotes mitochondrial dysfunction and trophoblast apoptosis in PE [150], and heat-shock proteins broadly are candidate biomarkers of placental ischemic disease [149].


- **CTNND1** (EAS4, EUR1). *Placental biology:* p120-catenin regulates trophoblast proliferation and lineage specification during human embryo implantation [79], and catenin expression (including p120) is tightly regulated during cytotrophoblast fusion into syncytium [81, 83]. The strongest RNA-editing QTL in the atlas, replicating across ancestries at the same molecular site.


- **PGAP6** (EAS5). *Biology/development:* PGAP6 is a GPI-specific phospholipase A2. It sheds the GPI-anchored Nodal co-receptor CRIPTO, a process important in early embryogenesis [328], and NCBI expression data show relatively high placental expression [338]. Direct pregnancy-disease evidence remains limited, but the shared EAS RNA-editing QTL connects a developmental signaling enzyme to post-transcriptional regulation.


- **CTSB** (EUR2–5). *Placental biology/disease:* cysteine cathepsins, including cathepsin B, are important for placental development and trophoblast invasion-related proteolysis [67], and serum cathepsin B is a candidate severity marker in preeclampsia [63]. Two genetically regulated editing sites in one gene.
 Two distinct CTSB editing sites remain in the EUR top five; one is FDR-significant in both ancestries and one is EUR-specific.


- **RBBP4** (EUR3). *Biology:* a chromatin-remodeling factor essential for preimplantation development — loss causes inner-cell-mass defects and lethality [185, 186] — and an epigenetic barrier to totipotency transitions [187]. Novel in placenta, but squarely in early-lineage chromatin biology.



#### 5.8 Transcript stability

The stability top five is essentially unchanged in gene composition. Four of ten top entries are shared, and the EUR list continues to include the canonical maternal–fetal immune locus HLA-C.


| Anc | Rank | Gene | Feature | Lead variant | Slope | AF | q | Shared |
|---|---|---|---|---|---|---|---|---|
| EAS | 1 | COX6C | ENSG00000164919.11 | 8:99848918:G:A | -0.95 | 0.119 | 2.91e-15 | no |
| EAS | 2 | PSG9 | ENSG00000183668.18 | 19:43303030:G:T | -0.72 | 0.368 | 1.18e-13 | no |
| EAS | 3 | PPIE | ENSG00000084072.17 | 1:39747206:G:A | -0.75 | 0.235 | 2.62e-12 | no |
| EAS | 4 | SDR39U1 | ENSG00000100445.18 | 14:24440269:C:G | +0.75 | 0.819 | 3.32e-12 | no |
| EAS | 5 | SH3YL1 | ENSG00000035115.22 | 2:254215:G:A | -0.67 | 0.224 | 3.93e-12 | yes |
| EUR | 1 | APIP | ENSG00000149089.13 | 11:34925165:A:G | -1.21 | 0.317 | 1.41e-25 | no |
| EUR | 2 | SNX19 | ENSG00000120451.11 | 11:130914721:T:C | -1.20 | 0.193 | 2.56e-21 | yes |
| EUR | 3 | GBP3 | ENSG00000117226.12 | 1:89009135:A:G | -1.07 | 0.261 | 3.47e-20 | yes |
| EUR | 4 | HLA-C | ENSG00000204525.18 | 6:31268539:G:T | -1.31 | 0.723 | 9.20e-17 | yes |
| EUR | 5 | RPS17 | ENSG00000182774.13 | 15:82550590:G:A | -1.02 | 0.367 | 5.61e-16 | no |



- **COX6C** (EAS1). *Placental disease:* cytochrome-c-oxidase activity is decreased in preeclamptic trophoblast [254], COX I mRNA is reduced in PE placentas [257], trophoblast mitochondrial dysfunction correlates with sFlt-1 expression [253], and mitochondrial dysfunction/oxidative stress are central to both early- and late-onset PE [252] and to adverse perinatal outcomes generally [263]. A stability QTL in a mitochondrial respiratory subunit fits this disease biology directly.


- **PSG9** (EAS2). *Placental biology/disease:* pregnancy-specific glycoproteins are placenta-secreted immunomodulators [3]; PSG9 (with PSG7) rises before preeclampsia onset [4]; endogenous retroviral elements LTR8B/MER65 rewire PSG9 regulation to control trophoblast syncytialization and PE risk [5]; PSG9 induces FoxP3+ regulatory T cells via TGF-β1 [6] and enhances endothelial nitric-oxide release, with lower serum PSG9 in PE [7]. One of the most placenta-specific genes in the entire hit set.


- **PPIE** (EAS3). *Biology:* a nuclear cyclophilin; cyclophilins regulate chromatin, transcription, and pre-mRNA splicing [235], and spliceosomal immunophilins act as isomerases promoting splicing efficiency [240]. Novel in placenta.


- **SDR39U1** (EAS4). *Biology:* a short-chain dehydrogenase/reductase-family epimerase [311] with a solved NADPH-bound structure [313]; placental steroidogenic pathways are gestational-age-dependent and trophoblast-centered [310], providing indirect context. Mostly novel.


- **SH3YL1** (EAS5). *Placental biology:* SH3YL1 binds phosphoinositides to drive dorsal ruffle formation [246] and cooperates with ESCRT-I in EGFR sorting/degradation [250]; placental RNA-landscape analyses report altered SH3YL1 expression in preeclampsia and fetal growth restriction [251]. A shared stability QTL with existing placental differential-expression evidence.


- **APIP** (EUR1). *Biology:* APIP suppresses hypoxic cell death via sustained AKT/ERK activation [143] and Apaf-1/caspase-9 inhibition [144], independently of its methionine-salvage enzymatic activity [146]; trophoblast hypoxia adaptation is a balance whose dysfunction drives gestational complications [145]. A coherent hypoxia-survival hit for the placenta.


- **SNX19** (EUR2). *Prior QTL evidence:* genetically predicted placental expression of SNX19 (among other genes) is associated with adult blood-pressure traits in a placental TWAS [194] — the clearest example in this set of a placental QTL already linked to a postnatal phenotype, supporting the developmental-origins hypothesis [52].


- **GBP3** (EUR3). *Biology:* guanylate-binding proteins are interferon-induced pattern-recognition/antimicrobial effectors [114, 115]. Novel in placenta; joins OAS1 and TRIM5 in the innate-immunity theme.


- **HLA-C** (EUR4). *Disease genetics:* maternal KIR–fetal HLA-C combinations are among the most established genetic interactions in preeclampsia — activating KIR (KIR2DS1) protects against reproductive failure when the fetus carries HLA-C2 [27], with supporting systematic reviews [22, 26], broader HLA immunogenetics at the feto-maternal interface [25], and central NK-cell involvement in PE [21]. The stability QTL adds a post-transcriptional regulatory layer to this canonical locus.


- **RPS17** (EUR5). *Placental biology:* mTORC1 transcriptionally regulates ribosome biogenesis and protein synthesis in primary human trophoblasts, linking ribosomal protein gene regulation to placental function and fetal growth [259]; RPS17 itself has mapped nuclear/nucleolar localization signals [261]. Novel as a stability QTL.



---

### 6. Synthesis

#### 6.1 Agreement with prior placental QTL and pregnancy literature

The strongest current loci continue to include direct placental-QTL positive controls. **ERAP2** and **ZSCAN9** were notable placental eGenes in the earlier Kikas et al. work [56], and the present atlas extends ERAP2 into transcript- and promoter-level regulation. **SNX19** remains notable because genetically predicted placental expression has previously been associated with adult blood-pressure traits in a placental TWAS [194].

The disease-context signal is also robust to the hit list. Genes with substantial independent preeclampsia or maternal–fetal-interface evidence include **ERAP1/ERAP2**, **HLA-C**, **LGALS14**, **CYP19A1**, **ANXA1**, **HSPD1**, **COX6C**, **PSG9**, **PECAM1**, **HSP90B1**, and **CTSB**. The splicing layer adds **PLAC4**, a particularly compelling placenta-specific gene with direct severe early-onset preeclampsia and syncytialization evidence [325], and **APOC1**, which now has both preeclampsia multi-omics and recurrent-pregnancy-loss mechanistic support [A5, A6].

These overlaps are biologically coherent but should not be presented as disease colocalization. The atlas maps molecular traits in placenta; formal inference about mediation of pregnancy-disease GWAS requires the dedicated colocalization/TWAS analyses.

#### 6.2 Recurring biological themes

1. **Immune regulation at the maternal–fetal interface:** HLA-C, ERAP1/2, OAS1, TRIM5, GBP3, LGALS14, APOC1, and CTNND1.
2. **Placenta-specific/endocrine programs:** PSG9, CYP19A1, PLAC4, and the galectin locus.
3. **Mitochondrial and cellular-stress biology:** COX6C, MRPL43, HSPD1, HSP90B1, PARP2, APIP.
4. **Trophoblast structure, invasion, and membrane organization:** IQGAP1, MYL12B, CTNND1, PECAM1, SCARB2, ERGIC3.
5. **RNA processing itself as a regulated phenotype:** RBM23, TENT2, TRMT13, TSR1, RPL14, PPIE, the strong isoform/polyA switch loci.
6. **Lipid/immune-metabolic regulation:** APOC1, PGAP6, SCARB2, GM2A, and the atlas-wide metabolic-pathway enrichment.

#### 6.3 Novelty

A substantial fraction of the 56 current top-hit genes still lack direct placental QTL or pregnancy-disease literature. The strongest “new biology” candidates include **RBM23, ERGIC3, PSCA, TRIM5, DRAM2, ELP5, PGAP6, SCARB2, ENOSF1, ZFP62, ZNF234, ZNF180, TSR1, METTL2B, GM2A, VAMP8, SDR39U1, PPIE, SOHLH2**, and the poorly characterized chr22 locus **ENSG00000290199**. Novelty here means that the targeted literature review did not identify direct prior placental/pregnancy evidence; it is not proof that no such paper exists.

The five explicit HPLRv2 top features are additionally novel at the phenotype-definition level and demonstrate the value of the long-read reference.

#### 6.4 A changed view of ancestry specificity

The results separate three phenomena that were conflated previously:

- **true feature-level replication**, abundant for alt-polyA, isoform expression, and RNA editing;
- **detection/phenotyping asymmetry**, especially intron retention and RNA editing;
- **covariate-specification sensitivity**, dramatically illustrated by EUR splicing.

Accordingly, absence of a shared top hit should not automatically be labeled an ancestry-specific biological effect. It can reflect phenotype observability, power, feature definition, or different LD/lead variants. The most defensible ancestry claims should come from matched-feature effect comparisons, cross-ancestry fine-mapping, and formal heterogeneity analyses.

### 7. Limitations

- **Rank-selected review.** Top five per ancestry/modality is a descriptive lens, not a genome-wide enrichment test and not a cross-ancestry meta-analysis.
- **Shared does not mean same causal variant.** The shared flag requires only FDR significance of the same molecular feature in both ancestries. Different lead variants can tag different causal configurations.
- **Gene attribution can be ambiguous.** This is especially important for the CMYA5/TENT2 splice cluster, APOC1P1/APOC1 cluster, CYP19A1 overlapping intron-retention annotations, and dense chromosome-19 zinc-finger regions.
- **Splicing interpretation is covariate-sensitive.** The corrected modality-matched HCP analysis resolves the prior EUR null, but EAS remains a low-yield layer and deserves cautious interpretation rather than overclaiming biological absence.
- **Coverage differs by ancestry and modality.** RNA editing is the clearest example (17,951 versus 3,120 normalized features).
- **Effect estimates are discovery estimates.** Large slopes, especially at lower-frequency variants, remain susceptible to winner's curse.
- **Abstract/reference-level literature review.** The literature review was primarily abstract-level. New additions were also targeted rather than systematic. “No direct literature found” should be read literally as a search result, not as proof of novelty.
- **No stale opposite-ancestry q-values.** Exact `q(other)` values for non-significant features cannot be reconstructed from the attached top-10 and FDR-only tables; this review intentionally omits them.

### 8. Recommended follow-up from the literature review

The strongest candidates for immediate integration with fine-mapping/colocalization are: **ERAP1/ERAP2** (positive control and pregnancy genetics), **PSCA** (multi-layer switch), **PECAM1** (endothelial and multi-layer), **PSG9/LGALS14/PLAC4** (placenta-specialized biology), **HLA-C/APOC1** (maternal–fetal immune interface), **CTNND1/HSPD1** (cross-ancestry RNA-editing replication), and **SCARB2/PGAP6/RBM23** as mechanistically interpretable but comparatively novel loci.

For splicing specifically, follow-up should treat the chr19 19:44433568:C:T set as a **single locus with multiple annotated consequences** until conditional/fine-mapping evidence separates ZNF234, APOC1P1/APOC1, and ZNF180.

### 9. Evidence-table additions in this review

The entries below are new to the 2026-10-10 revision. They have now been appended to the companion evidence table as references **322–339**, while all original evidence-table indices remain unchanged.


- **[322] RBM23.** Dowhan DH, Hong EP, Auboeuf D, et al. Steroid hormone receptor coactivation and alternative RNA splicing by U2AF65-related proteins CAPERalpha and CAPERbeta. Mol Cell. 2005;17(3):429-439. PMID 15694343; DOI 10.1016/j.molcel.2004.12.025. RBM23/CAPERβ links transcriptional coactivation to alternative splicing. https://pubmed.ncbi.nlm.nih.gov/15694343/

- **[323] CMYA5.** Durham JT, Brand OM, Arnold M, et al. Myospryn is a direct transcriptional target for MEF2A that encodes a striated muscle, alpha-actinin-interacting, costamere-localized protein. J Biol Chem. 2006;281(10):6841-6849. PMID 16407236; DOI 10.1074/jbc.M510499200. CMYA5/myospryn is a MEF2A target and costamere-associated structural/signaling protein. https://pubmed.ncbi.nlm.nih.gov/16407236/

- **[324] TENT2.** Warkocki Z, Liudkovska V, Gewartowska O, Mroczek S, Dziembowski A. Terminal nucleotidyl transferases (TENTs) in mammalian RNA metabolism. Philos Trans R Soc Lond B Biol Sci. 2018;373(1762):20180162. PMID 30397099; DOI 10.1098/rstb.2018.0162. Functional context for TENT2/GLD2 and mammalian RNA 3′-end metabolism. https://pubmed.ncbi.nlm.nih.gov/30397099/

- **[325] PLAC4.** Tuohey L, Macintire K, Ye L, et al. PLAC4 is upregulated in severe early onset preeclampsia and upregulated with syncytialisation but not hypoxia. Placenta. 2013;34(3):256-260. PMID 23312075; DOI 10.1016/j.placenta.2012.12.009. PLAC4 is placenta-enriched, syncytiotrophoblast-localized, syncytialization-responsive, and increased in severe early-onset PE. https://pubmed.ncbi.nlm.nih.gov/23312075/

- **[326] APOC1.** Chen Y, Wu L, Huang D, et al. Subtype specific immune-metabolic reprogramming in preeclampsia revealed by multiomics and serum biomarkers. Hypertens Res. 2026;49(3):641-657. PMID 41419624; DOI 10.1038/s41440-025-02504-5. Placental single-cell/spatial multi-omics identify APOC1 among dysregulated lipid-transport regulators in PE, with cell- and subtype-specific patterns. https://pubmed.ncbi.nlm.nih.gov/41419624/

- **[327] APOC1.** Du YM, Zhang JX, Lu S, et al. APOC1 deficiency disrupts M1/M2 macrophage balance via ferroptosis and suppresses trophoblast invasion in recurrent pregnancy loss. Cell Biosci. 2026. DOI 10.1186/s13578-026-01646-4. Macrophage APOC1 deficiency drives ferroptosis/inflammation and impairs trophoblast invasion in RPL models. https://doi.org/10.1186/s13578-026-01646-4

- **[328] PGAP6.** Lee GH, Fujita M, Nakanishi H, et al. PGAP6, a GPI-specific phospholipase A2, has narrow substrate specificity against GPI-anchored proteins. J Biol Chem. 2020;295(42):14501-14509. PMID 32816994; DOI 10.1074/jbc.RA120.014643. PGAP6 sheds CRIPTO and modulates a developmental Nodal-signaling co-receptor. https://pubmed.ncbi.nlm.nih.gov/32816994/

- **[329] SCARB2.** Human Protein Atlas. SCARB2 tissue expression: placenta. . High protein staining in syncytiotrophoblast and Hofbauer cells; medium in endothelium. https://www.proteinatlas.org/ENSG00000138760-SCARB2/tissue/placenta

- **[330] ENOSF1.** NCBI Gene 55556: ENOSF1 enolase superfamily member 1. . Functional annotation includes mitochondrial L-fuconate metabolism and antisense/TYMS-regulatory transcripts. https://www.ncbi.nlm.nih.gov/gene/55556

- **[331] TRMT13.** NCBI Gene 54482: TRMT13 tRNA methyltransferase 13. . Predicted tRNA methyltransferase involved in tRNA methylation. https://www.ncbi.nlm.nih.gov/gene/54482

- **[332] TRMT13.** Human Protein Atlas. TRMT13 tissue expression: placenta. . Placental immunohistochemistry reports medium decidual and low trophoblastic protein signal. https://www.proteinatlas.org/ENSG00000122435-TRMT13/tissue/placenta

- **[333] ZFP62.** NCBI Gene 643836: ZFP62 zinc finger protein. . Reference annotation for poorly characterized ZFP62. https://www.ncbi.nlm.nih.gov/gene/643836

- **[334] ZNF234.** NCBI Gene 10780: ZNF234 zinc finger protein 234. . Predicted DNA-binding transcription regulator on 19q13.31. https://www.ncbi.nlm.nih.gov/gene/10780

- **[335] APOC1P1.** NCBI Gene 342: APOC1P1 apolipoprotein C1 pseudogene 1. . Reference annotation for APOC1P1 and its pseudogene status. https://www.ncbi.nlm.nih.gov/gene/342

- **[336] ZNF180.** NCBI Gene 7733: ZNF180 zinc finger protein 180. . Kruppel-related C2H2 zinc-finger protein. https://www.ncbi.nlm.nih.gov/gene/7733

- **[337] TSR1.** NCBI Gene 55720: TSR1 ribosome maturation factor. . RNA-binding/nucleolar factor implicated in pre-rRNA processing and SSU maturation. https://www.ncbi.nlm.nih.gov/gene/55720

- **[338] PGAP6.** NCBI Gene 58986: PGAP6 post-GPI attachment to proteins 6. . Reference annotation reports phospholipase A2/lipid-metabolism function and relatively high placental expression. https://www.ncbi.nlm.nih.gov/gene/58986

- **[339] RPL14.** NCBI Gene 9045: RPL14 ribosomal protein L14. . Reference annotation for RPL14 as a large-subunit ribosomal protein. https://www.ncbi.nlm.nih.gov/gene/9045


### 10. Companion files

- `literature_xqtl_mapping_revised_20261010.md` — this review.
- `xqtl_top_hits_current_20261010.tsv` — exact current top-five-per-ancestry/modality entries used here, with the recomputed shared flag.
- `xqtl_literature_revision_sources_20261010.tsv` — structured bibliography for `[322]`–`[339]`.
- Original `xqtl_top_hits_evidence_table.csv` — still required to resolve legacy numeric citations `[N]`.

