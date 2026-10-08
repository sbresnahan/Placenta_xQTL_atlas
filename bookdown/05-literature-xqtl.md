# (PART) Literature companions {-}

# Top placental xQTL hits in the context of the published literature {-}

**Project:** Placenta xQTL atlas (multi-ancestry, multi-modality cis-QTL mapping)
**Date:** 2026-10-05
**Companion report:** `report_xqtl_mapping.html` (Report 1: xQTL mapping) — this document is its literature-review companion; the in-report "Literature context" section points here.
**Companion files:** `xqtl_top_hits_ungrouped.tsv` (full 80-hit table), `xqtl_top_hits_evidence_table.csv` (bibliographic details and key-finding notes for every cited reference)

Bracketed numbers `[N]` are citation indices that resolve to full bibliographic records in the evidence table (`xqtl_top_hits_evidence_table.csv`, column `index`).

***

## 1. Background and scope {-}

The placenta xQTL atlas maps cis genetic regulation across eight RNA-processing modalities — gene expression, isoform expression, isoform ratios, splicing, intron retention, alternative transcription start sites (alt_TSS), alternative polyadenylation (alt_polyA), RNA editing, and transcript stability — phenotyped against the HPLRv2 long-read placental transcriptome annotation and tested with tensorQTL under GTEx conventions (Storey q-values) in two ancestry strata (EUR, n ≈ 136 mapped; EAS, n ≈ 272 mapped) drawn from the NIEHS_RICHS, GUSTO, SNUH (PRJNA820329), and NIGMS (PRJNA671171) cohorts.

This report examines the **top hits of the ungrouped top-variant scan** — the scan in which every molecular feature was tested individually and the single best-associated variant per feature is reported — in the context of the PubMed literature. For each of the eight modalities that have an ungrouped scan (isoform expression, isoforms, splicing, intron retention, alt_TSS, alt_polyA, RNA editing, stability), the **top 5 features by q-value in each ancestry** were selected (80 hits, 58 unique genes). Each gene is annotated with a short structured snapshot covering three angles: **placental biology**, **prior QTL evidence**, and **pregnancy-disease genetics**. Loci whose association is shared across ancestries are presented first within each modality.

## 2. Methods {-}

- **Hit selection.** From each ungrouped top-variant scan (one tensorQTL file per modality × ancestry), features were ranked by Storey q-value (ties broken by nominal p-value) and the top 5 per ancestry were retained. Expression and combined modalities have no ungrouped scan and are excluded. Two grouped stepwise-conditional runs (EAS splicing, EAS combined) were still in progress and are not part of this scan.
- **Cross-ancestry sharing.** A hit is flagged **shared** when the *same feature* reaches q < 0.05 in both the EUR and EAS ungrouped scans (the lead variants may differ). `q(other)` is the feature's q-value in the opposite ancestry; `NA` means the variant was not testable there (typically because it is too rare).
- **Literature synthesis.** Targeted PubMed searches were run for each of the 58 unique genes (gene name + placenta/trophoblast/pregnancy/preeclampsia/QTL terms) plus four umbrella queries (placental eQTL studies, preeclampsia GWAS, gestational diabetes genetics, fetal growth/placental weight GWAS). Synthesis is at the abstract level; the 175 references cited here are indexed in the evidence table. Claims about the atlas itself (sample sizes, cohorts, annotation) come from the project record, not the literature.

## 3. Overview of the 80 top hits {-}

| Modality | EUR sig. features (q<0.05) | EAS sig. features (q<0.05) | Shared among top 5×2 |
|---|---|---|---|
| Isoform expression | 1,509 / 39,554 tested | 1,519 / 36,441 | 7 / 10 |
| Isoforms | 255 | 408 | 6 / 10 |
| Splicing | 16,117 / 75,470 (21.4%) | 2,211 / 56,149 | 0 / 10 |
| Intron retention | 113 | 57 | 0 / 10 |
| alt_TSS | 988 | 471 | 4 / 10 |
| alt_polyA | 1,069 | 565 | 10 / 10 |
| RNA editing | 294 / 17,484 | 44 / 3,120 | 7 / 10 |
| Stability | 206 | 148 | 4 / 10 |

**38 of 80 hits are shared across ancestries** at the feature level. Sharing is concentrated in the dosage-like modalities (alt_polyA, isoform expression, RNA editing) and absent among splicing and intron-retention top hits, which are dominated by rare, ancestry-specific variants (see §4.3).

## 4. Cross-cutting patterns {-}

### 4.1 Recurring multi-modality loci {-}

The same lead variant (or same locus) recurs across independent modalities, marking regulatory hotspots:

- **ERAP1/ERAP2 locus (chr5q15).** Variant 5:96943786 is the top EUR hit for ERAP2 isoform expression *and* ERAP2 alt_TSS (novel HPLRv2 promoter); variant 5:96773906 dominates ERAP1 isoforms (EUR ranks 1–2) and ERAP1 alt_polyA (EUR ranks 1, 2, 4). This is the atlas's strongest convergence with prior pregnancy-disease genetics (§5.1, §5.2).
- **PSCA (chr8q24).** Variant 8:142680513 is an isoform-switch QTL in EAS (isoforms ranks 1–2, slopes −0.79/+0.78) and the top EAS alt_polyA locus (ranks 1–2), with feature-level sharing to EUR.
- **HSP90B1 (chr12q23).** Variant 12:103946200 is a top-5 EUR hit in both isoform expression and alt_polyA; the EAS alt_polyA hit is a nearby indel (12:103943444) — cross-ancestry replication of the same feature with different leads.
- **MRPL43 (chr10q24).** Variant 10:100963159 is the top EAS intron-retention hit and a top-5 EAS alt_polyA hit; the locus is a known antagonistic cis-eQTL [265].
- **ERGIC3 (chr20q13).** Variant 20:35549038 is shared EUR–EAS for isoform expression and alt_TSS, and is also the EAS intron-retention rank-2 hit.
- **PSG9 (chr19q13).** Variant 19:43303030 is a top EAS hit for both intron retention and stability — a pregnancy-specific glycoprotein under retroelement control [5].
- **MYL12B (chr18q21).** Variant 18:3278477 produces opposite-direction alt_TSS effects on two promoters (EAS ranks 4–5), shared with EUR.
- Others: **ZNF439/ZNF440** (one variant, splicing EUR ranks 4–5), **LILRB1/LILRB4** (one rare variant, splicing EAS ranks 1–4), **PRKCH-AS1** (alt_TSS EUR ranks 3–4, two nearby variants), **CTSB** (RNA editing EUR ranks 2 and 5, two sites), **CTNND1** and **HSPD1** (RNA editing shared across ancestries with different lead variants), **PECAM1** (isoform expression EUR3 + alt_polyA EUR5, same variant).

### 4.2 Isoform-switch QTLs {-}

Three loci show the same variant pushing two isoforms/promoters of one gene in opposite directions — the signature of isoform usage regulation rather than expression dosage: **PSCA** (isoforms EAS, −0.79/+0.78), **ERAP1** (isoforms EUR, +1.28/−1.27), and **MYL12B** (alt_TSS EAS, −0.51/+0.60). These are prime candidates for mechanism-follow-up because the molecular consequence is a change in transcript *identity*, not just abundance.

### 4.3 Splicing and intron-retention top hits are rare-variant and ancestry-specific {-}

All ten splicing top hits have lead-variant allele frequencies of 0.011–0.016 in the discovery stratum and were not testable in the other ancestry (`q(other) = NA`): CNN2, DOCK8, SREK1, ZNF439/ZNF440 (EUR); LILRB1/LILRB4, RDX (EAS). The same holds for all ten intron-retention hits despite common-variant leads (AF 0.15–0.61) — cross-ancestry replication there is limited by feature detection, not just variant frequency. These hits should be treated as ancestry-stratified signals requiring within-ancestry confirmation.

### 4.4 Novel HPLRv2 transcripts among the top hits {-}

Seven top hits are carried by transcripts that exist only in the HPLRv2 long-read annotation (feature IDs `HPLRT_*`): TRIM5 and DRAM2 novel isoforms (isoforms EUR ranks 3–4), an ERAP2 novel promoter (alt_TSS EUR1), PRKCH-AS1 and TEX10 promoters (alt_TSS EUR ranks 3 and 5), a PECAM1 polyA site (alt_polyA EUR5), an HSP90B1 polyA site (alt_polyA EAS4), and a novel chr22 lncRNA (ENSG00000290199; isoform expression EAS3, shared with EUR). These associations are invisible to short-read annotations and are a concrete yield of the long-read resource.

### 4.5 RNA-editing coverage asymmetry {-}

Five times more editing sites were testable in EUR than EAS (17,484 vs 3,120), so the RNA-editing modality's EUR dominance (294 vs 44 significant sites) partly reflects coverage rather than biology. Cross-ancestry shared editing QTLs (CTNND1, HSPD1, LCMT1, CNOT12, CTSB) are therefore the most robust signals in this modality.

***

## 5. Modality-by-modality results and gene snapshots {-}

### 5.1 Isoform expression {-}

| Anc | Rank | Gene | Feature | Lead variant | Slope | AF | q | q(other) | Shared |
|---|---|---|---|---|---|---|---|---|---|
| EAS | 1 | PGM1 | ENST00000371084.8 | 1:63616853:A:G | −0.55 | 0.129 | 9.9e-27 | 2.7e-19 | yes |
| EAS | 2 | ERGIC3 | ENST00000489071.5 | 20:35549038:C:T | 1.25 | 0.133 | 5.1e-20 | 1.6e-08 | yes |
| EAS | 3 | HPLRv2-novel-lncRNA | ENST00000703580.1 | 22:23915150:G:A | 2.54 | 0.448 | 2.9e-19 | 0.017 | yes |
| EAS | 4 | FBLN1 | ENST00000262722.11 | 22:45483551:C:T | 0.50 | 0.554 | 5.1e-19 | 0.23 | no |
| EAS | 5 | DRC9 | ENST00000478903.5 | 3:197947511:C:T | 4.44 | 0.043 | 5.9e-19 | 0.18 | no |
| EUR | 1 | ERAP2 | ENST00000437043.8 | 5:96943786:A:G | 1.72 | 0.415 | 2.6e-30 | 2.6e-15 | yes |
| EUR | 2 | ZSCAN9 | ENST00000527436.5 | 6:28233360:A:G | 0.92 | 0.423 | 5.3e-25 | 3.7e-07 | yes |
| EUR | 3 | PECAM1 | ENST00000563924.6 | 17:64322640:G:A | −0.63 | 0.552 | 2.7e-24 | 1.6e-09 | yes |
| EUR | 4 | SOHLH2 | ENST00000379881.8 | 13:36218513:A:G | 0.89 | 0.452 | 4.6e-23 | 2.2e-09 | yes |
| EUR | 5 | HSP90B1 | ENST00000550595.2 | 12:103946200:G:A | 0.88 | 0.342 | 9.5e-22 | NA | no |

**Shared loci**

- **ERAP2** (EUR1; also alt_TSS EUR1, same variant). *Placental biology/disease:* ERAP2 is one of the most replicated preeclampsia candidate genes: associated with preeclampsia in Australian and Norwegian cohorts [110], fetal ERAP2 rs2549782 associated with preeclampsia in African Americans [112], a Chilean haplotype study links rs2549782 to ERAP2 protein expression and PE risk [113], and ERAP1/ERAP2 show altered placental expression and lymphocyte activation in PE [107]; ERAP2 was also associated with PE (and ERAP1 with eclampsia) in Brazilian women [108]. *Prior QTL evidence:* ERAP2 was one of the notable eGenes in the Kikas et al. placental eQTL study [56] — the present hit is a direct cross-ancestry replication and extension to isoform- and promoter-level regulation.
- **ZSCAN9** (EUR2). *Prior QTL evidence:* ZSCAN9 was likewise a notable placental eGene in Kikas et al. [56], making this a second direct replication. *Biology:* it is a KRAB zinc-finger protein; KRAB-ZFPs repress transposable elements and shape mammalian regulatory evolution [316, 317], a mechanism of particular relevance in placenta (see §5.2, ZNF880; §5.3, ZNF439/440).
- **PECAM1** (EUR3; also alt_polyA EUR5, same variant). *Placental biology:* CD31/PECAM-1 marks placental endothelium; microvessel density stained for CD31 is markedly reduced in preeclamptic placentas [222] and PECAM-1 expression is decreased in PE placentae [218], though one study found no difference in CD31 expression or vascular growth [223], and PECAM-1 itself does not appear to drive trophoblast invasion or PE/FGR pathophysiology [219]. The QTL is shared across ancestries and acts on both isoform abundance and polyA-site usage.
- **SOHLH2** (EUR4). *Biology:* a germ-cell-specific bHLH transcription factor essential for oocyte and spermatogonial differentiation [202, 203, 204]. No prior placental QTL or pregnancy-disease association — a novel placental regulatory signal, plausibly reflecting its role in germline/early-lineage transcription programs.
- **PGM1** (EAS1). *Placental biology:* phosphoglucomutase isozymes were purified from human placenta decades ago [199], but the gene's modern literature is dominated by PGM1-CDG, a congenital glycosylation disorder treatable with galactose [200, 201]. No prior placental QTL or pregnancy-disease genetics — novel in this context.
- **ERGIC3** (EAS2; also intron retention EAS2 and alt_TSS EAS2, same variant). *Biology:* an ER–Golgi intermediate compartment protein; the ERGIC is a stress-responsive hub integrating membrane trafficking and adaptation [155], and ERGIC3-dependent secretory cargo trafficking is regulated by the E3 ligase MARCH2 [156]. No placental literature — a novel, strongly shared, multi-modality locus.
- **HPLRv2-novel-lncRNA (ENSG00000290199)** (EAS3). A novel lncRNA on chr22 present only in the HPLRv2 annotation; by definition there is no prior literature. Its large effect (slope 2.54) at intermediate allele frequency in both ancestries makes it a priority for functional characterization.

**Ancestry-specific loci**

- **FBLN1** (EAS4). *Placental biology:* fibulin-1 is expressed during mouse placental development, with elevated transcripts in abnormal placentation [121]; in human pregnancy, plasma fibulin-1 rises across gestation and is consistently lower in PPROM [124]. No prior QTL evidence.
- **DRC9** (EAS5). *Biology:* a dynein-regulatory-complex (N-DRC) subunit (aliases CFAP122/IQCG) required for ciliary/flagellar motility [279, 281]; the mouse ortholog Iqcg is essential for sperm flagellum formation [282]. Novel in placenta; notable for the largest effect size in this modality (slope 4.44) at low allele frequency (AF 0.043).
- **HSP90B1** (EUR5; shared with EAS in alt_polyA). *Placental biology:* gp96/GRP94 is expressed at the materno-fetal interface, most intensely in trophoblast and glandular epithelium [134], and modulates decidual antigen-presenting and NK cell activation [135]. *Disease:* unfolded-protein-response activation, including HSP90B1-related signaling, distinguishes early-onset from late-onset preeclampsia placentae [136].

### 5.2 Isoforms {-}

| Anc | Rank | Gene | Feature | Lead variant | Slope | AF | q | q(other) | Shared |
|---|---|---|---|---|---|---|---|---|---|
| EAS | 1 | PSCA | ENST00000513264.1 | 8:142680513:C:T | −0.79 | 0.407 | 1.5e-18 | 4.6e-08 | yes |
| EAS | 2 | PSCA | ENST00000301258.5 | 8:142680513:C:T | +0.78 | 0.407 | 1.5e-18 | 7.0e-06 | yes |
| EAS | 3 | OAS1 | ENST00000202917.10 | 12:112943944:T:C | −0.93 | 0.748 | 1.5e-18 | NA | no |
| EAS | 4 | METTL2B | ENST00000480046.5 | 7:128477173:G:A | 0.89 | 0.296 | 5.9e-17 | NA | no |
| EAS | 5 | PARP2 | ENST00000250416.9 | 14:20344308:C:A | 1.04 | 0.176 | 8.1e-17 | 1.0e-04 | yes |
| EUR | 1 | ERAP1 | ENST00000296754.7 | 5:96773906:G:A | +1.28 | 0.287 | 1.6e-25 | NA | no |
| EUR | 2 | ERAP1 | ENST00000443439.7 | 5:96773906:G:A | −1.27 | 0.287 | 2.5e-21 | 0.14 | no |
| EUR | 3 | TRIM5 | HPLRT_f3660beda5b | 11:5938806:G:A | −1.09 | 0.430 | 1.2e-15 | 4.5e-07 | yes |
| EUR | 4 | DRAM2 | HPLRT_a8773b56079 | 1:111100631:TC:T | 1.23 | 0.280 | 1.2e-15 | 3.5e-07 | yes |
| EUR | 5 | ZNF880 | ENST00000600321.5 | 19:52384174:A:C | −0.94 | 0.401 | 1.8e-15 | 4.2e-07 | yes |

**Shared loci**

- **PSCA** (EAS1–2; also alt_polyA EAS1–2, same variant). A textbook isoform-switch QTL: the same allele lowers one PSCA isoform and raises the other. *Literature:* PSCA's published record is entirely oncologic — expression rises with prostate cancer grade, stage, and metastasis [213], with context-dependent tumor-promoting or -suppressing roles [214] and prognostic/therapeutic interest across cancers [215]. There is no placental or pregnancy literature: this is a novel, cross-ancestry-shared placental isoform QTL for a GPI-anchored cell-surface antigen.
- **PARP2** (EAS5). *Disease relevance:* PARP-family activity is mechanistically tied to placental dysfunction — NAD+ depletion from NAD+-consuming enzymes (including PARPs) marks an inflammatory subclass of preeclampsia and is reversible with nicotinamide riboside in models [271, 272], and PARP activation is elevated in gestational diabetic pregnancies [268]. *Biology:* PARP2 has distinct and shared functions with PARP1 in inflammation, metabolism, and oxidative stress [270].
- **TRIM5** (EUR3; novel HPLRv2 isoform). *Biology:* TRIM5α is a canonical retroviral restriction factor that recognizes and destabilizes incoming capsids [173], within the broader TRIM antiviral family [170]. No placental literature — novel; given the placenta's constitutive antiviral state (see OAS1 below), a plausible innate-immunity locus.
- **DRAM2** (EUR4; novel HPLRv2 isoform). *Biology:* a lysosomal autophagy regulator; biallelic mutations cause retinal dystrophy [276], and despite divergence from DRAM it contributes to autophagy induction [277, 278]. Novel in placenta, where autophagy is central to trophoblast homeostasis.
- **ZNF880** (EUR5). *Biology:* a KRAB zinc-finger protein; KRAB-ZFPs repress transposable elements and drove regulatory innovations in mammalian evolution [316], show diverse homeostatic functions in overexpression screens [317], engage in arms races with endogenous retroelements [321], and include placenta-essential members (murine ZFP568 controls Igf2 for embryo-placental development [318]). Recurrent KRAB-ZFP hits (ZSCAN9, ZNF439/440) suggest transposable-element control is a heritable regulatory axis in placenta.

**Ancestry-specific loci**

- **OAS1** (EAS3). *Placental biology:* OAS1 is an innate immune sensor of dsRNA [95]; the placenta maintains baseline type I interferon/ISG signaling that preserves pregnancy homeostasis and antiviral readiness [97], and hemochorial placentas constitutively express type III interferon in trophoblasts [98]. A biologically coherent EAS-specific isoform QTL in the placental antiviral pathway.
- **METTL2B** (EAS4). *Biology:* an m3C tRNA methyltransferase with mutually exclusive substrate selection against METTL6 [283]; m3C32 modification controls serine codon-biased translation, cell cycle, and DNA-damage response [286]. Novel in placenta.
- **ERAP1** (EUR1–2; also alt_polyA EUR ranks 1, 2, 4, same variant). An isoform-switch QTL mirroring PSCA. *Disease:* ERAP1 expression is elevated in preeclamptic placentas and inducible by hypoxia in trophoblasts [111]; ERAP1/ERAP2 placental expression and lymphocyte activation are altered in PE [107]; ERAP1 is associated with eclampsia in Brazilian women [108]. *Prior QTL evidence:* the adjacent ERAP2 was a notable placental eGene [56]; the present data resolve the locus into isoform- and polyA-level ERAP1 regulation.

### 5.3 Splicing {-}

All ten hits are rare-variant (AF 0.011–0.016), ancestry-specific signals; none were testable in the other ancestry.

| Anc | Rank | Gene | Feature | Lead variant | Slope | AF | q | q(other) | Shared |
|---|---|---|---|---|---|---|---|---|---|
| EAS | 1 | LILRB1 | chr19_54635609_54636732_clu_11802_+ | 19:54557553:T:C | 0.54 | 0.013 | 3.6e-24 | NA | no |
| EAS | 2 | LILRB1 | chr19_54635609_54667873_clu_11802_+ | 19:54557553:T:C | 0.54 | 0.013 | 3.6e-24 | NA | no |
| EAS | 3 | LILRB4 | chr19_54635609_54636732_clu_11802_+ | 19:54557553:T:C | 0.54 | 0.013 | 3.6e-24 | NA | no |
| EAS | 4 | LILRB4 | chr19_54635609_54667873_clu_11802_+ | 19:54557553:T:C | 0.54 | 0.013 | 3.6e-24 | NA | no |
| EAS | 5 | RDX | chr11_110264234_110272536_clu_4344_- | 11:110158016:T:A | 0.56 | 0.013 | 5.6e-22 | NA | no |
| EUR | 1 | CNN2 | chr19_1032696_1036067_clu_11295_+ | 19:1128719:C:T | 1.56 | 0.012 | 2.2e-37 | NA | no |
| EUR | 2 | DOCK8 | chr9_215029_260162_clu_22024_+ | 9:722423:G:T | 1.79 | 0.011 | 2.0e-30 | NA | no |
| EUR | 3 | SREK1 | chr5_66175041_66191156_clu_17707_+ | 5:65792660:GT:G | 1.36 | 0.011 | 3.3e-30 | NA | no |
| EUR | 4 | ZNF439 | chr19_11947574_11978116_clu_11482_+ | 19:12090144:C:T | 1.29 | 0.016 | 1.8e-29 | NA | no |
| EUR | 5 | ZNF440 | chr19_11947574_11978116_clu_11482_+ | 19:12090144:C:T | 1.29 | 0.016 | 1.8e-29 | NA | no |

- **LILRB1/LILRB4** (EAS1–4, one variant). *Placental biology:* LILRB1 is the receptor through which trophoblast HLA-G modulates maternal antigen-presenting cells [101], with structural basis for LILRB recognition of HLA-G isoforms [102]; LILRB1/2 proteins are expressed in placental stromal and perivascular cells (not trophoblast) [106]. *Disease:* KIR/LILRB/ligand polymorphisms have been proposed as biomarkers in recurrent implantation failure [103]. A rare EAS variant rewiring splicing across this immune-checkpoint cluster is a strong maternal-fetal-tolerance candidate. (Note: LILRB1 and LILRB4 splice clusters overlap genomically; gene attribution is ambiguous at this locus.)
- **RDX** (EAS5). *Placental biology:* the ERM family (ezrin/radixin/moesin) builds placental microvilli — ezrin alone is ~5% of microvillus protein [212] — and ERMs integrate membrane-cortex signaling [211]; placental ezrin is increased in gestational diabetes [207]. RDX-specific placental literature is thin, so this is largely novel.
- **CNN2** (EUR1). *Biology:* calponin-2 regulates the actin cytoskeleton in smooth muscle and non-muscle cells [126] and is essential for vascular development and endothelial migration [129]. *Disease:* serum calponin-2 discriminates ectopic pregnancy from viable pregnancy (AUC ≈ 0.93) [125, 130]. A rare EUR splicing variant in a cytoskeletal/vascular gene with an existing pregnancy-biomarker link.
- **DOCK8** (EUR2). *Biology:* DOCK8 is a guanine-exchange factor central to innate and adaptive immunity [176, 177]; deficiency causes severe combined immunodeficiency with allergy, infection, and autoimmunity [178]. Novel in placenta.
- **SREK1** (EUR3). *Prior QTL/splicing evidence:* placental alternative splicing is extensive (≈149k events in decidual cells) and involves SREK1/SFRS12-related modulation [245]; placental splicing is partly genetically driven in cis and trans [242]. SREK1 mis-splicing (exon 10 inclusion) is oncogenic in hepatocellular carcinoma [243]. The hit extends genetically driven SREK1 splicing regulation to placenta.
- **ZNF439/ZNF440** (EUR4–5, one variant). KRAB-ZFP transposable-element repressors [316, 317, 321]; placenta-essential KRAB-ZFP precedent (ZFP568–Igf2) [318]. Same theme as ZNF880/ZSCAN9 — here acting through splicing of a chr19 KRAB-ZFP cluster.

### 5.4 Intron retention {-}

No top hit is shared at the feature level (all `q(other) = NA`), reflecting detection asymmetry between ancestries rather than variant rarity (most leads are common).

| Anc | Rank | Gene | Feature | Lead variant | Slope | AF | q | q(other) | Shared |
|---|---|---|---|---|---|---|---|---|---|
| EAS | 1 | MRPL43 | IR_ENSG00000055950 (chr10) | 10:100963159:G:C | −0.75 | 0.530 | 2.4e-16 | NA | no |
| EAS | 2 | ERGIC3 | IR_ENSG00000125991 (chr20) | 20:35549038:C:T | 0.92 | 0.133 | 1.7e-13 | NA | no |
| EAS | 3 | PSG9 | IR_ENSG00000183668 (chr19) | 19:43303030:G:T | 0.75 | 0.368 | 7.2e-13 | NA | no |
| EAS | 4 | LGALS8 | IR_ENSG00000116977 (chr1) | 1:236542978:G:T | −0.84 | 0.151 | 1.6e-10 | NA | no |
| EAS | 5 | VAMP8 | IR_ENSG00000118640 (chr2) | 2:85585623:T:C | −0.70 | 0.398 | 2.8e-10 | NA | no |
| EUR | 1 | GM2A | IR_ENSG00000196743 (chr5) | 5:151262704:A:G | 1.07 | 0.261 | 6.6e-15 | NA | no |
| EUR | 2 | IQGAP1 | IR_ENSG00000140575 (chr15) | 15:90359318:C:T | 1.07 | 0.211 | 8.0e-13 | NA | no |
| EUR | 3 | GLDN | IR_ENSG00000186417 (chr15) | 15:51377428:A:C | 0.87 | 0.222 | 8.7e-13 | NA | no |
| EUR | 4 | LGALS14 | IR_ENSG00000006659 (chr19) | 19:39709274:C:G | 0.79 | 0.607 | 9.1e-13 | NA | no |
| EUR | 5 | CYP19A1 | IR_ENSG00000137869 (chr15) | 15:51253257:A:G | −0.72 | 0.330 | 9.2e-13 | NA | no |

- **MRPL43** (EAS1; also alt_polyA EAS5, same variant). *Prior QTL evidence:* a single 97 kb cis-eQTL spanning MRPL43 has antagonistic effects across regulatory stages — raising mRNA and ribosome occupancy while lowering protein — via linked variants affecting TF binding, splicing, and miRNA targeting [265]; the present intron-retention and polyA QTLs fit that multi-layer architecture. *Placental relevance:* nuclear-encoded mitochondrial ribosomal proteins are required to initiate gastrulation [262], and placental mitochondrial dysfunction is a recurring feature of preeclampsia, IUGR/FGR, and preterm birth [263, 264].
- **ERGIC3** (EAS2) — see §5.1; the same variant also drives the EAS intron-retention signal.
- **PSG9** (EAS3; also stability EAS2, same variant). *Placental biology/disease:* pregnancy-specific glycoproteins are placenta-secreted immunomodulators [3]; PSG9 (with PSG7) rises before preeclampsia onset [4]; endogenous retroviral elements LTR8B/MER65 rewire PSG9 regulation to control trophoblast syncytialization and PE risk [5]; PSG9 induces FoxP3+ regulatory T cells via TGF-β1 [6] and enhances endothelial nitric-oxide release, with lower serum PSG9 in PE [7]. One of the most placenta-specific genes in the entire hit set.
- **LGALS8** (EAS4). *Placental biology:* galectins collectively shape early pregnancy and pregnancy pathologies [18]. Note that the well-studied chr19 placental galectin cluster (LGALS13/14/16) does not include LGALS8 (chr1); LGALS8-specific placental evidence is limited — largely novel.
- **VAMP8** (EAS5). *Biology:* a SNARE for regulated exocytosis, selective for granule-to-granule fusion [181] and general across exocrine secretion [182, 183]. No placental literature — novel; secretory machinery is nonetheless central to syncytiotrophoblast hormone release.
- **GM2A** (EUR1). *Biology:* the GM2-activator protein is a lysosomal lipid-transport cofactor for GM2 ganglioside hydrolysis [224, 226, 228]. Novel in placenta.
- **IQGAP1** (EUR2). *Placental biology:* IQGAP1 is a scaffold promoting cell motility and invasion via Cdc42/Rac1 [165, 168] and specifically promotes trophoblast motility and invasion [167] — directly relevant to extravillous trophoblast function.
- **GLDN** (EUR3). *Disease genetics:* gliomedin variants cause lethal congenital contracture syndrome 11 / fetal akinesia deformation sequence [298, 301] — a fetal-developmental gene with no prior placental regulatory evidence.
- **LGALS14** (EUR4). *Placental biology:* galectin-14 is a placenta-specific, primate galectin that promotes trophoblast migration and invasion via MMP-9/N-cadherin [15], belongs to a chr19 cluster evolved at the maternal-fetal interface to induce maternal immune-cell death [16], modulates immune-cell survival and cytokines [19], and protects trophoblasts from oxidative stress while promoting extravillous differentiation [20]. *Disease:* the chr19 galectin cluster is complexly dysregulated in preeclampsia [17]. A high-priority placenta-specific hit.
- **CYP19A1** (EUR5). *Placental biology/disease:* aromatase drives placental estrogen synthesis; its downregulation is proposed to play a dual (causal and predictive) role in preeclampsia [43], placental estradiol is markedly decreased in PE [47], and oxygen-dependent CYP19A1 induction during trophoblast differentiation is mediated by ERRγ/HIF-1 [50], with tissue-specific promoter usage long established [48]. Intron retention is a plausible post-transcriptional layer on this heavily regulated gene.

### 5.5 Alternative TSS {-}

| Anc | Rank | Gene | Feature | Lead variant | Slope | AF | q | q(other) | Shared |
|---|---|---|---|---|---|---|---|---|---|
| EAS | 1 | ACBD3 | grp_2_upstream_ENST00000366812.6 | 1:225333399:G:C | 0.77 | 0.011 | 8.7e-23 | 0.40 | no |
| EAS | 2 | ERGIC3 | grp_2_upstream_ENST00000489071.5 | 20:35549038:C:T | 0.92 | 0.133 | 1.6e-15 | 0.005 | yes |
| EAS | 3 | MMP14 | grp_1_upstream_ENST00000311852.11 | 14:22821187:C:G | −0.55 | 0.013 | 2.2e-13 | 0.71 | no |
| EAS | 4 | MYL12B | grp_1_upstream_ENST00000237500.10 | 18:3278477:A:G | −0.51 | 0.566 | 6.7e-13 | 3.3e-04 | yes |
| EAS | 5 | MYL12B | grp_1_upstream_ENST00000400175.9 | 18:3278477:A:G | +0.60 | 0.566 | 6.7e-13 | 1.2e-06 | yes |
| EUR | 1 | ERAP2 | grp_2_upstream_HPLRT_039d234c78d | 5:96943786:A:G | −0.86 | 0.415 | 1.7e-17 | 0.002 | yes |
| EUR | 2 | RAF1 | grp_1_upstream_ENST00000691643.1 | 3:11779825:C:T | 1.46 | 0.019 | 7.1e-11 | NA | no |
| EUR | 3 | PRKCH-AS1 | grp_1_upstream_HPLRT_85bb29b580c | 14:60360050:A:G | 1.24 | 0.019 | 2.0e-10 | 0.71 | no |
| EUR | 4 | PRKCH-AS1 | grp_1_upstream_ENST00000661303.1 | 14:60416948:TCA:T | −1.23 | 0.019 | 8.2e-10 | NA | no |
| EUR | 5 | TEX10 | grp_1_upstream_HPLRT_7a1fc1a5d6d | 9:99584699:T:C | 1.48 | 0.012 | 2.9e-09 | 0.65 | no |

**Shared loci**

- **ERAP2** (EUR1) — see §5.1; here the QTL acts on a novel HPLRv2 promoter, adding promoter-level regulation to the isoform-level signal.
- **ERGIC3** (EAS2) — see §5.1; third modality for this variant.
- **MYL12B** (EAS4–5). A promoter-switch QTL (opposite slopes on two TSS groups). *Placental biology:* myosin regulatory light chains maintain myosin II stability and cellular integrity [230]; human placental stem villi contain ultraslow myosin motors [233] whose contraction uses actin–myosin cross-bridges under NO–cGMP control [234]. A contractility-relevant gene with shared cross-ancestry promoter regulation.

**Ancestry-specific loci**

- **ACBD3** (EAS1). *Biology:* a Golgi scaffold organizing stacking proteins and Rab33b-GAP [137], a multifunctional membrane-domain organizer [138] required for Golgi integrity and glucosylceramide trafficking [142]. Novel in placenta; rare EAS lead (AF 0.011).
- **MMP14** (EAS3). *Placental biology:* MT1-MMP is produced in early human placenta and implicated in implantation [58], is developmentally regulated with other MMPs in first-trimester trophoblast [60], co-localizes with MMP-2 in extravillous cytotrophoblasts [61], is required for leptin-promoted extravillous trophoblast invasion via Notch1–PI3K/Akt crosstalk [62], and is down-regulated by endothelin-1, reducing trophoblast invasion [59]. One of the best-characterized invasion genes in the hit set.
- **RAF1** (EUR2). *Placental biology:* Raf-1 loss causes midgestation lethality with placental anomalies and fetal liver apoptosis [91]; the Raf/ERK pathway is essential for extraembryonic vascular development [87] and its components are expressed in preimplantation embryos and trophoblast stem cells [90]. Rare EUR lead (AF 0.019).
- **PRKCH-AS1** (EUR3–4). An antisense lncRNA to PRKCH with no direct published literature; the two hits include a novel HPLRv2 promoter. Analogous placental lncRNAs (e.g., PROX1-AS1) regulate trophoblast migration/invasion and preeclampsia-relevant phenotypes [287], making this a plausible but uncharacterized regulatory RNA.
- **TEX10** (EUR5). *Biology:* a pluripotency factor that fine-tunes Wnt signaling in primordial germ cell and male germline development [294, 297]; the related mammalian-specific Tex19.1 is essential for spermatogenesis and placenta-supported development [296]. Novel in placenta; rare EUR lead (AF 0.012).

### 5.6 Alternative polyA {-}

All ten hits are shared across ancestries — the most cross-ancestry-robust modality.

| Anc | Rank | Gene | Feature | Lead variant | Slope | AF | q | q(other) | Shared |
|---|---|---|---|---|---|---|---|---|---|
| EAS | 1 | PSCA | grp_2_downstream_ENST00000513264.1 | 8:142680513:C:T | −0.87 | 0.407 | 3.7e-23 | 2.0e-13 | yes |
| EAS | 2 | PSCA | grp_2_downstream_ENST00000301258.5 | 8:142680513:C:T | +0.86 | 0.407 | 3.7e-23 | 1.3e-13 | yes |
| EAS | 3 | ELP5 | grp_2_downstream_ENST00000574255.6 | 17:7254884:T:G | −0.89 | 0.655 | 1.7e-22 | 5.3e-20 | yes |
| EAS | 4 | HSP90B1 | grp_1_downstream_HPLRT_270ba7eb207 | 12:103943444:GTAAT:G | 0.72 | 0.434 | 2.8e-21 | 2.0e-22 | yes |
| EAS | 5 | MRPL43 | grp_2_downstream_ENST00000299179.9 | 10:100963159:G:C | 0.81 | 0.530 | 3.2e-21 | 3.7e-21 | yes |
| EUR | 1 | ERAP1 | grp_1_downstream_ENST00000443439.7 | 5:96773906:G:A | −1.20 | 0.287 | 1.3e-37 | 6.4e-06 | yes |
| EUR | 2 | ERAP1 | grp_1_downstream_ENST00000296754.7 | 5:96773906:G:A | +1.20 | 0.287 | 2.6e-36 | 5.2e-06 | yes |
| EUR | 3 | HSP90B1 | grp_2_downstream_ENST00000550595.2 | 12:103946200:G:A | 0.79 | 0.342 | 9.9e-30 | 3.2e-21 | yes |
| EUR | 4 | ERAP1 | grp_2_downstream_ENST00000443439.7 | 5:96773906:G:A | −1.29 | 0.287 | 6.5e-29 | 1.3e-04 | yes |
| EUR | 5 | PECAM1 | grp_1_downstream_HPLRT_683224049f0 | 17:64322640:G:A | 0.87 | 0.552 | 1.8e-28 | 3.2e-06 | yes |

- **PSCA** (EAS1–2) — see §5.2; the same variant switches both isoform ratios and polyA-site usage, suggesting the isoform switch is driven by alternative 3′-end choice.
- **ELP5** (EAS3). *Biology:* Elongator-complex subunit; the complex's two subcomplexes have divergent neurodevelopmental roles [290, 291] and are linked to human neurological disease [292]. No placental literature — a novel, strongly shared (q ≈ 1e-20 in both ancestries) polyA QTL.
- **HSP90B1** (EAS4, EUR3) — see §5.1. Cross-ancestry replication of the same polyA feature with different lead variants (EAS indel 12:103943444; EUR SNP 12:103946200), plus the EUR isoform-expression hit — a three-signal locus.
- **MRPL43** (EAS5) — see §5.4; same variant as the intron-retention hit, consistent with the known multi-layer antagonistic eQTL at this gene [265].
- **ERAP1** (EUR1, 2, 4) — see §5.2; three of the top four EUR alt_polyA hits are ERAP1 features under one variant, mirroring the isoform-switch pattern at the 3′ end.
- **PECAM1** (EUR5) — see §5.1; same variant as the isoform-expression hit, acting on a novel HPLRv2 polyA site.

### 5.7 RNA editing {-}

| Anc | Rank | Gene | Editing site | Lead variant | Slope | AF | q | q(other) | Shared |
|---|---|---|---|---|---|---|---|---|---|
| EAS | 1 | ANXA1 | chr9_73160319 | 9:73163667:AG:A | 0.91 | 0.185 | 1.1e-12 | 0.59 | no |
| EAS | 2 | LCMT1 | chr16_25147778 | 16:25114437:T:C | 0.75 | 0.349 | 1.3e-12 | 0.033 | yes |
| EAS | 3 | CTNND1 | chr11_57761957 | 11:57790196:C:T | 1.04 | 0.120 | 4.6e-12 | 3.3e-27 | yes |
| EAS | 4 | HSPD1 | chr2_197497294 | 2:197446138:A:G | 0.71 | 0.541 | 5.3e-12 | 2.1e-22 | yes |
| EAS | 5 | CNOT12 | chr11_57309990 | 11:57309990:T:C | 0.83 | 0.307 | 1.2e-11 | 2.3e-04 | yes |
| EUR | 1 | CTNND1 | chr11_57761957 | 11:57826368:G:A | −1.16 | 0.566 | 3.3e-27 | 4.6e-12 | yes |
| EUR | 2 | CTSB | chr8_11844424 | 8:11844866:G:A | 1.43 | 0.257 | 2.2e-24 | NA | no |
| EUR | 3 | RBBP4 | chr1_32680141 | 1:32669173:A:C | 1.42 | 0.752 | 5.1e-23 | NA | no |
| EUR | 4 | HSPD1 | chr2_197497294 | 2:197400449:T:A | 1.25 | 0.701 | 2.1e-22 | 5.3e-12 | yes |
| EUR | 5 | CTSB | chr8_11845033 | 8:11847833:G:C | −1.19 | 0.634 | 6.5e-20 | 1.2e-08 | yes |

**Shared loci**

- **CTNND1** (EAS3, EUR1 — same editing site, different lead variants). *Placental biology:* p120-catenin regulates trophoblast proliferation and lineage specification during human embryo implantation [79], and catenin expression (including p120) is tightly regulated during cytotrophoblast fusion into syncytium [81, 83]. The strongest RNA-editing QTL in the atlas, replicating across ancestries at the same molecular site.
- **HSPD1** (EAS4, EUR4 — same site, different leads). *Disease:* Hsp60 is overexpressed in preeclamptic placentas [152], circulating Hsp60 is elevated in PE [148], Hsp60 lactylation promotes mitochondrial dysfunction and trophoblast apoptosis in PE [150], and heat-shock proteins broadly are candidate biomarkers of placental ischemic disease [149].
- **LCMT1** (EAS2). *Biology:* leucine carboxyl methyltransferase-1 methylates PP4/PP6 and regulates phosphatase holoenzyme assembly [160]; global loss causes fetal-liver hematopoietic failure and embryonic lethality [158]; protein is expressed in placental trophoblasts [161]. Novel as an editing QTL.
- **CNOT12** (EAS5). *Biology:* a CCR4–NOT complex subunit; the complex couples deadenylation to translation and mRNA decay [304, 307]. Novel in placenta — a shared editing QTL in the machinery that itself controls RNA stability (cf. the stability modality).
- **CTSB** (EUR5; EUR2 is a second, EUR-specific site). *Placental biology/disease:* cysteine cathepsins, including cathepsin B, are important for placental development and trophoblast invasion-related proteolysis [67], and serum cathepsin B is a candidate severity marker in preeclampsia [63]. Two genetically regulated editing sites in one gene.

**Ancestry-specific loci**

- **ANXA1** (EAS1). *Disease:* annexin A1 inhibits trophoblast ferroptosis in preeclampsia via KISS1 downregulation [71]; silencing ANXA1 suppresses apoptosis and inflammation in PE rat trophoblasts [72]; plasma ANXA1 is elevated in PE, especially early-onset [78]; placental ANXA1 influences DNA-damage response in gestational diabetes models [77]; the annexin family is broadly implicated at the maternal-fetal interface [75]. Rich disease literature; the editing QTL is EAS-specific.
- **RBBP4** (EUR3). *Biology:* a chromatin-remodeling factor essential for preimplantation development — loss causes inner-cell-mass defects and lethality [185, 186] — and an epigenetic barrier to totipotency transitions [187]. Novel in placenta, but squarely in early-lineage chromatin biology.

### 5.8 Transcript stability {-}

| Anc | Rank | Gene | Lead variant | Slope | AF | q | q(other) | Shared |
|---|---|---|---|---|---|---|---|---|
| EAS | 1 | COX6C | 8:99848918:G:A | −0.95 | 0.119 | 2.3e-15 | 0.81 | no |
| EAS | 2 | PSG9 | 19:43303030:G:T | −0.72 | 0.368 | 8.9e-14 | 0.57 | no |
| EAS | 3 | PPIE | 1:39747206:G:A | −0.75 | 0.235 | 1.9e-12 | 0.22 | no |
| EAS | 4 | SDR39U1 | 14:24440269:C:G | 0.75 | 0.819 | 2.1e-12 | 0.65 | no |
| EAS | 5 | SH3YL1 | 2:254215:G:A | −0.67 | 0.224 | 3.3e-12 | 1.3e-09 | yes |
| EUR | 1 | APIP | 11:34925165:A:G | −1.21 | 0.317 | 8.4e-27 | 0.76 | no |
| EUR | 2 | SNX19 | 11:130914721:T:C | −1.20 | 0.193 | 2.7e-21 | 2.5e-09 | yes |
| EUR | 3 | GBP3 | 1:89009135:A:G | −1.07 | 0.261 | 2.1e-20 | 1.8e-09 | yes |
| EUR | 4 | HLA-C | 6:31268539:G:T | −1.31 | 0.723 | 1.5e-17 | 1.5e-11 | yes |
| EUR | 5 | RPS17 | 15:82550590:G:A | −1.02 | 0.367 | 2.4e-17 | NA | no |

**Shared loci**

- **SH3YL1** (EAS5). *Placental biology:* SH3YL1 binds phosphoinositides to drive dorsal ruffle formation [246] and cooperates with ESCRT-I in EGFR sorting/degradation [250]; placental RNA-landscape analyses report altered SH3YL1 expression in preeclampsia and fetal growth restriction [251]. A shared stability QTL with existing placental differential-expression evidence.
- **SNX19** (EUR2). *Prior QTL evidence:* genetically predicted placental expression of SNX19 (among other genes) is associated with adult blood-pressure traits in a placental TWAS [194] — the clearest example in this set of a placental QTL already linked to a postnatal phenotype, supporting the developmental-origins hypothesis [52].
- **GBP3** (EUR3). *Biology:* guanylate-binding proteins are interferon-induced pattern-recognition/antimicrobial effectors [114, 115]. Novel in placenta; joins OAS1 and TRIM5 in the innate-immunity theme.
- **HLA-C** (EUR4). *Disease genetics:* maternal KIR–fetal HLA-C combinations are among the most established genetic interactions in preeclampsia — activating KIR (KIR2DS1) protects against reproductive failure when the fetus carries HLA-C2 [27], with supporting systematic reviews [22, 26], broader HLA immunogenetics at the feto-maternal interface [25], and central NK-cell involvement in PE [21]. The stability QTL adds a post-transcriptional regulatory layer to this canonical locus.

**Ancestry-specific loci**

- **COX6C** (EAS1). *Placental disease:* cytochrome-c-oxidase activity is decreased in preeclamptic trophoblast [254], COX I mRNA is reduced in PE placentas [257], trophoblast mitochondrial dysfunction correlates with sFlt-1 expression [253], and mitochondrial dysfunction/oxidative stress are central to both early- and late-onset PE [252] and to adverse perinatal outcomes generally [263]. A stability QTL in a mitochondrial respiratory subunit fits this disease biology directly.
- **PSG9** (EAS2) — see §5.4; same variant as the intron-retention hit.
- **PPIE** (EAS3). *Biology:* a nuclear cyclophilin; cyclophilins regulate chromatin, transcription, and pre-mRNA splicing [235], and spliceosomal immunophilins act as isomerases promoting splicing efficiency [240]. Novel in placenta.
- **SDR39U1** (EAS4). *Biology:* a short-chain dehydrogenase/reductase-family epimerase [311] with a solved NADPH-bound structure [313]; placental steroidogenic pathways are gestational-age-dependent and trophoblast-centered [310], providing indirect context. Mostly novel.
- **APIP** (EUR1). *Biology:* APIP suppresses hypoxic cell death via sustained AKT/ERK activation [143] and Apaf-1/caspase-9 inhibition [144], independently of its methionine-salvage enzymatic activity [146]; trophoblast hypoxia adaptation is a balance whose dysfunction drives gestational complications [145]. A coherent hypoxia-survival hit for the placenta.
- **RPS17** (EUR5). *Placental biology:* mTORC1 transcriptionally regulates ribosome biogenesis and protein synthesis in primary human trophoblasts, linking ribosomal protein gene regulation to placental function and fetal growth [259]; RPS17 itself has mapped nuclear/nucleolar localization signals [261]. Novel as a stability QTL.

***

## 6. Synthesis {-}

**Agreement with prior placental QTL studies.** The atlas's top hits directly replicate two of the most notable placental eGenes from Kikas et al. — ERAP2 and ZSCAN9 [56] — and are consistent with the broader picture that placental eQTLs are numerous, often placenta-specific, and mediate GWAS signals for postnatal traits [52, 53], including in East Asian placentas [51], across gestational stages [55], and in pathological placentas [54]. Genetically driven splicing variation in placenta [242] is strongly confirmed — splicing is the largest modality by significant-feature count in EUR (21.4% of tested features). The SNX19 stability QTL connects to the one existing placental-TWAS result (adult blood pressure) [194].

**Convergence on pregnancy-disease biology.** Despite no disease phenotype being used in the mapping, the top hits are enriched for genes with independent preeclampsia/pregnancy evidence: ERAP1/ERAP2 [107, 108, 110, 111, 112, 113], HLA-C [22, 26, 27], LGALS14 [15, 17, 20], CYP19A1 [43, 47], MMP14 [58, 61, 62], ANXA1 [71, 78], HSPD1 [148, 150, 152], COX6C/mitochondrial function [252, 253, 254, 257], PSG9 [4, 5, 7], PECAM1 [218, 222], HSP90B1 [136], CTSB [63], and CNN2 [125, 130]. This mirrors the locus-level themes of preeclampsia GWAS (FLT1 and blood-pressure/endothelial pathways [30, 31]) and suggests cis-regulatory variation in placenta is a plausible mediating layer, as hypothesized by developmental-origins eQTL work [52] and placental multi-omics integration for birthweight [14]. Gestational-diabetes genetics (HKDC1, MTNR1B [35, 38]) and placental-weight GWAS [8] do not overlap the top hits directly, but PARP2/NAD+ biology links to both PE and GDM [268, 271, 272].

**Novelty.** Roughly half of the 58 genes have no prior placental or pregnancy literature (e.g., PSCA, TRIM5, DRAM2, ERGIC3, ACBD3, ELP5, CNOT12, RBBP4, GM2A, VAMP8, METTL2B, DOCK8, SDR39U1, PPIE, GLDN, SOHLH2, PGM1, DRC9, TEX10, and the novel chr22 lncRNA). These are the atlas's most original contributions: common-variant, often cross-ancestry-shared regulatory effects on genes whose placental role is uncharacterized. The seven hits carried by novel HPLRv2 transcripts (§4.4) are novel by construction.

**Recurring biological themes.** (i) *Immune regulation at the maternal-fetal interface* — HLA-C, LILRB1/4, OAS1, TRIM5, GBP3, ERAP1/2, LGALS14. (ii) *KRAB-ZFP/transposable-element control* — ZSCAN9, ZNF880, ZNF439/440, plus retroelement-wired PSG9 [5]. (iii) *Mitochondrial/oxidative-stress axis* — COX6C, MRPL43, HSPD1, HSP90B1, PARP2, APIP. (iv) *Trophoblast invasion and cytoskeleton* — MMP14, IQGAP1, CNN2, MYL12B, RDX, CTNND1. (v) *Secretory/membrane trafficking* — ERGIC3, ACBD3, VAMP8, SH3YL1, SNX19.

**Gaps.** Splicing and intron-retention top hits lack any cross-ancestry replication (rare variants and detection asymmetry, respectively); RNA-editing coverage is 5.6× deeper in EUR; and the two ongoing grouped conditional runs (EAS splicing, EAS combined) may reorder which loci are independent within genes.

## 7. Limitations {-}

- **Selection by rank, not genome-wide novelty.** Top-5-per-stratum is a descriptive cut of the strongest associations, not an FDR-controlled cross-ancestry meta-analysis; shared flags reflect feature-level q < 0.05 in both scans and do not require the same lead variant or colocalization.
- **Winner's curse and rare variants.** Effect sizes (slopes) for the rare-variant splicing hits are upwardly biased; several leads (AF ≈ 0.01) need within-ancestry confirmation.
- **Gene attribution ambiguity.** At gene-cluster loci (LILRB1/LILRB4, ZNF439/ZNF440, chr19 galectins), splice/retention features span multiple genes; the assigned symbol is the best positional match.
- **Abstract-level literature synthesis.** Snapshots rely on titles/abstracts of targeted searches, not full-text review; absence of prior evidence means "none found by these searches," not proof of absence. Literature for paralogs (e.g., ezrin for RDX, Tex19.1 for TEX10) is flagged as family-level context only.
- **Ancestry asymmetries.** EAS has ~2× the samples of EUR, and editing-site coverage differs 5.6×; cross-ancestry comparisons of hit counts confound power with biology.

## 8. Files {-}

- `report_placenta_xqtl_top_hits_literature.md` — this report
- `xqtl_top_hits_ungrouped.tsv` — master table of all 80 hits (modality, ancestry, rank, gene, feature, variant, slope, AF, q, other-ancestry q, shared flag)
- `xqtl_top_hits_evidence_table.csv` — all 179 cited references with bibliographic metadata and key-finding notes; the `index` column matches the `[N]` citations above
