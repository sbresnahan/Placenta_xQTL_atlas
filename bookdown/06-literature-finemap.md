# Fine-mapped placental xQTL genes in the context of the published literature {-}

*Companion literature review to `report_finemap.html` (Report 2). It is the
fine-mapping counterpart to `literature_xqtl_mapping.md`, which reviews the
top variant-level association hits from Report 1.*

## 1. Background and scope {-}

Report 2 (`report_finemap.html`) summarizes cross-ancestry (EUR + EAS) SuSHiE
fine-mapping of every placental xQTL locus across nine molecular modalities
(gene expression, isoform expression, isoform ratios, splicing, intron
retention, alternative TSS, alternative polyA, RNA editing, and transcript
stability). Fine-mapping collapses each associated locus to one or more 95%
credible sets (CSs) of candidate causal variants and assigns every variant a
posterior inclusion probability (PIP). The high-confidence gene set used here
is the best credible set per locus with lead-variant PIP >= 0.9, which
nominates **5,236 genes** with at least one confidently resolved regulatory
signal.

Where Report 1's companion review asked "which genes sit under the strongest
association peaks," this review asks a sharper question enabled by
fine-mapping: **which genes have a confidently resolved, often single-variant,
causal regulatory signal, and what is already known about them in placental
biology, pregnancy outcomes, and intragenomic conflict?** Placental eQTL and
multi-omic studies have established that common regulatory variants shape the
placental transcriptome and connect to developmental and pregnancy phenotypes
[6, 7, 8, 9]. Because the placenta is the arena of maternal-fetal conflict
over nutrient transfer [1, 2, 3, 4], and because placental regulatory
variation is increasingly tied to preeclampsia, fetal growth, and gestational
duration [10, 11, 12], we read the fine-mapped catalog through three lenses:

- **Placental biology** — genes that define trophoblast identity,
  syncytialization, and placental function, including genes with little or no
  prior placental characterization (candidate novel biology).
- **Pregnancy and birth outcomes** — genes with established or candidate links
  to preeclampsia, fetal growth, parturition, and related complications.
- **Intragenomic conflict** — imprinted genes and other loci where
  parent-of-origin selection (Haig's kinship theory) is expected to leave
  strong, sometimes antagonistic, regulatory signatures [1, 2, 5].

## 2. Methods {-}

- **High-confidence set.** From `finemap_credible_sets.tsv.gz`, the credible
  set with the highest lead-variant PIP was taken per locus; a locus enters the
  high-confidence set when that lead PIP >= 0.9. This yields 5,236 genes. The
  lead-PIP distribution is strongly bimodal (median ~0.93), so PIP >= 0.9 alone
  is permissive; cross-modality recurrence and credible-set size are used as
  secondary evidence (Report 2, "Credible sets").
- **Candidate ranking.** Genes were ranked by (i) the number of modalities in
  which they reach the high-confidence threshold (cross-modality recurrence),
  (ii) maximum lead PIP, and (iii) minimum credible-set size, and were flagged
  against a priori lists: imprinted genes (geneimprint), preeclampsia GWAS
  candidates, syncytiotrophoblast (STB) syncytialization markers, HPA
  placenta-enriched genes, and KRAB-zinc-finger proteins.
- **Literature screen.** All 264 ranked/flagged candidate genes were screened
  against PubMed (E-utilities `esearch`, 0.34 s rate limit) with three
  lens-specific queries per gene (placenta/trophoblast; pregnancy/birth
  outcomes; imprinting/intragenomic conflict). Hit counts and the top PMIDs per
  lens were recorded (`finemap_lit_screen_candidates.tsv`).
- **Deep-dive selection.** Twenty genes were chosen to span the three lenses,
  balancing (a) cross-modality recurrence, (b) strength and specificity of the
  existing literature, and (c) novelty. For each, lens-specific PubMed queries
  were re-run with `sort=relevance` and the top records retrieved (`efetch`),
  with citation counts from the NIH iCite API. Synthesis is at the abstract
  level; the 59 references cited here are indexed in
  `finemap_hits_evidence_table.csv`, whose `index` column matches the `[N]`
  citations below. Claims about the atlas itself (PIP values, credible-set
  sizes, modality recurrence) come from the project data, not the literature.

## 3. Overview of the 20 deep-dive genes {-}

| Gene | Lens | Modalities | Max lead PIP | Min CS size | Prior flag |
|------|------|-----------|--------------|-------------|------------|
| WSB1 | Placental (novel) | 7 | 1.00 | 1 | none |
| GBP3 | Placental (novel) | 6 | 1.00 | 1 | none |
| CYP19A1 | Placental | 5 | 1.00 | 1 | STB marker |
| CAST | Placental (novel) | 4 | 1.00 | 1 | none |
| SIGLEC6 | Placental | 4 | 1.00 | 1 | — |
| GCM1 | Placental | 1 | 1.00 | 1 | STB marker |
| GSTM1 | Placental | 3 | 1.00 | 1 | — |
| ERAP1 | Pregnancy | 6 | 1.00 | 1 | PE GWAS |
| FLT1 | Pregnancy | 1 | 1.00 | 1 | PE GWAS |
| PAPPA | Pregnancy | 1 | 1.00 | 1 | — |
| HLA-C | Pregnancy | 1 | 1.00 | 1 | PE GWAS |
| ACVR2A | Pregnancy | 1 | 1.00 | 1 | PE GWAS |
| HSD11B2 | Pregnancy | 1 | 1.00 | 1 | STB marker |
| IL1R1 | Pregnancy | 4 | 1.00 | 1 | — |
| IGF2 | Conflict | 1 | 1.00 | 1 | imprinted |
| IGF2R | Conflict | 1 | 1.00 | 1 | imprinted |
| MEG3 | Conflict | 1 | 1.00 | 1 | imprinted |
| PEG3 | Conflict | 1 | 1.00 | 1 | imprinted |
| MEST | Conflict | 1 | 1.00 | 1 | imprinted |
| GNAS | Conflict | 1 | 1.00 | 1 | imprinted |

Every selected gene reaches the high-confidence threshold with a single-variant
credible set (min CS size = 1) and lead PIP ~1.0 in at least one modality, so
each represents a confidently resolved regulatory signal rather than a diffuse
association peak.

## 4. Cross-cutting patterns {-}

### 4.1 The most recurrent signals are the least characterized {-}

The genes with the broadest cross-modality regulatory control are precisely
those with the thinnest placental literature. **WSB1** (7 modalities) and
**GBP3** (6 modalities) each carry a single-variant credible set at lead PIP
= 1.0 across six to seven independent molecular layers, yet neither has any
prior placental characterization (Section 5.1). This inversion — strongest
genetic evidence, least biological annotation — is the clearest discovery
signal in the fine-mapped catalog and mirrors a pattern seen in Report 1, where
several top hits were novel HPLRv2 transcripts.

### 4.2 Splicing QTLs dominate the pregnancy and conflict lenses {-}

Eleven of the thirteen pregnancy-outcome and conflict genes are fine-mapped
**only** through splicing (FLT1, PAPPA, HLA-C, ACVR2A, HSD11B2, IGF2, MEG3,
PEG3, MEST, GNAS) or expression (IGF2R). This is consistent with Report 2's
headline enrichment — high-confidence splicing-QTL genes are enriched for
syncytiotrophoblast syncytialization markers (FDR ~ 8e-5) — and suggests that
post-transcriptional regulation, more than transcript abundance, is the
readout through which many pregnancy-relevant loci act in term placenta.

### 4.3 Fine-mapping recovers known pregnancy genes as positive controls {-}

The pregnancy lens is anchored by genes with strong prior evidence — FLT1 [10,
11, 31], PAPPA [33, 34], HLA-C [35, 36], ACVR2A [37, 38], ERAP1 [28, 29] —
each resolved to a single-variant credible set. That fine-mapping
re-nominates these established loci with high confidence cross-validates the
pipeline and calibrates interpretation of the novel signals.

### 4.4 Imprinted genes are enriched and cleanly resolved {-}

Six canonical imprinted genes (IGF2, IGF2R, MEG3, PEG3, MEST, GNAS) carry
single-variant credible sets, part of a broader pool of 27 imprinted genes in
the high-confidence set. Their clean resolution is notable because imprinting
places effectively all regulatory weight on one parental allele, so regulatory
variants at these loci can have outsized, parent-of-origin-specific effects on
placental and fetal growth [1, 2, 3].

## 5. Gene snapshots by lens {-}

### 5.1 Placental biology {-}

**WSB1 (WD repeat and SOCS box 1) — 7 modalities, lead PIP 1.0, single-variant
CS.** The single most recurrent high-confidence gene in the atlas, with
resolved signals in alt_TSS, expression, intron retention, isoform expression,
RNA editing, splicing, and stability. WSB1 encodes the substrate-recognition
subunit of an E3 ubiquitin ligase and is a HIF-1 target that promotes
degradation of VHL and HIPK2, positioning it at the center of the hypoxia
response that dominates early placentation [13]. Critically, WSB1 produces at
least three splice isoforms with distinct, partly opposing effects on cell
proliferation and apoptosis [14] — so the multi-layer splicing/isoform signal
we resolve is likely functionally consequential, not neutral. No placental
role has been reported; WSB1 is a priority candidate for functional follow-up.

**GBP3 (guanylate-binding protein 3) — 6 modalities, lead PIP 1.0,
single-variant CS.** Resolved across alt_polyA, expression, isoform
expression, isoforms, splicing, and stability. GBPs are interferon-inducible
GTPases that act as intracellular pattern-recognition proteins [15]; GBP3
specifically governs caspase-4 activation on cytosolic bacteria to trigger
pyroptotic cell death [16]. Cell-autonomous innate immunity at the
maternal-fetal interface is biologically plausible (defense against ascending
infection) but unstudied for GBP3 in placenta — a second high-confidence novel
nomination.

**CYP19A1 (aromatase) — 5 modalities, lead PIP 1.0, single-variant CS; STB
marker.** Aromatase catalyzes estrogen synthesis and is one of the most
placenta-specific genes known: its promoter drives trophoblast-exclusive
expression (the basis of the Cyp19a1-Cre placenta-specific mouse) [17], and
placental CYP19A1 also metabolizes xenobiotics such as nifedipine [18].
Aromatase/CYP activity is dysregulated in preeclampsia organoid models [19].
The strong multi-modality signal here is a reassuring positive control for a
bona fide placenta-enriched regulator.

**CAST (calpastatin) — 4 modalities, lead PIP 1.0, single-variant CS.**
Calpastatin is the endogenous inhibitor of the calcium-dependent calpain
proteases. The calpain/calpastatin system acts at mitochondria-associated ER
membranes to regulate calcium homeostasis, autophagy, and apoptosis, and its
components are altered in preeclamptic placentas [20]. CAST's cross-modality
signal (isoform expression, isoforms, RNA editing, splicing) is essentially
unstudied in human placenta and represents a third novel, mechanistically
grounded candidate.

**SIGLEC6 (sialic acid-binding Ig-like lectin 6) — 4 modalities, lead PIP 1.0,
single-variant CS.** SIGLEC6 is among the most placenta-specific genes and an
emerging preeclampsia biomarker: circulating SIGLEC6 is elevated ~9.5-fold in
preterm preeclampsia and tracks disease severity across seven cohorts [21].
Mechanistically, SIGLEC6 overexpression impairs trophoblast syncytialization
and promotes macrophage M1 polarization in early-onset disease [22], and on
syncytiotrophoblasts it mediates extracellular-vesicle uptake through a
noncanonical glycolipid-binding pocket [23]. The fine-mapped regulatory
variants may therefore influence both placental development and a clinically
measured circulating biomarker.

**GCM1 (glial cells missing 1) — splicing, lead PIP 1.0, single-variant CS;
STB marker.** GCM1 is the master transcription factor of syncytiotrophoblast
differentiation, activating syncytin and hCG-beta [24] and sitting at the
center of the DeltaNp63alpha/GCM1 switch that governs the
cytotrophoblast-to-syncytiotrophoblast transition [25]. A single-variant
splicing credible set at GCM1 directly implicates post-transcriptional control
of the syncytialization program — dovetailing with the report-level enrichment
of splicing QTLs for STB markers.

**GSTM1 (glutathione S-transferase mu 1) — 3 modalities, lead PIP 1.0,
single-variant CS.** GSTM1 is a central glutathione-conjugating enzyme in
cellular oxidative-stress defense [27]. Although the common GSTM1-null
deletion shows no robust association with preeclampsia [26], our signal is a
finely resolved expression/splicing QTL, suggesting regulatory (rather than
null-allele) modulation of placental glutathione metabolism — consistent with
the KEGG glutathione-metabolism enrichment among fine-mapped expression genes
(Report 2).

### 5.2 Pregnancy and birth outcomes {-}

**ERAP1 (endoplasmic reticulum aminopeptidase 1) — 6 modalities, lead PIP 1.0,
single-variant CS; PE GWAS candidate.** ERAP1 trims peptides for MHC-I loading
and degrades vasoactive peptides; it is expressed in placenta, changes in
preeclampsia from the first trimester, and associates with immune-cell
activation [29], and its variants are associated with eclampsia/preeclampsia
[28]. ERAP1 sits at the trophoblast-immune interface that governs spiral-artery
remodeling [30]. Its six-modality recurrence (expression, intron retention,
isoform expression, RNA editing, splicing, stability) makes it one of the most
heavily regulated pregnancy genes in the atlas; the adjacent ERAP2 (the
balancing-selection positive control of Report 1) shares the same haplotype
background.

**FLT1 (fms-like tyrosine kinase 1) — splicing, lead PIP 1.0, single-variant
CS; PE GWAS.** FLT1 encodes the VEGF receptor whose soluble splice isoform,
sFlt-1, is the anti-angiogenic driver of the maternal preeclampsia phenotype;
the sFlt-1/PlGF ratio is a frontline clinical biomarker [31, 32]. Fetal-genome
GWAS identified FLT1 as a robust preeclampsia risk locus [10, 11]. That we
resolve a **splicing** credible set at FLT1 is mechanistically apt: the
pathogenic species is itself a splice isoform, so genetically regulated
splicing at FLT1 is a direct candidate mechanism for modulating the
sFlt-1/PlGF balance.

**PAPPA (pregnancy-associated plasma protein A) — splicing, lead PIP 1.0,
single-variant CS.** PAPP-A is a metalloprotease that frees insulin-like
growth factors from their binding proteins, thereby promoting fetal growth;
low first-trimester PAPP-A is an established predictor of preeclampsia and
adverse outcomes [33] and a core component of first-trimester competing-risks
screening [34]. A single-variant splicing signal at PAPPA nominates a
regulatory handle on a clinically deployed biomarker.

**HLA-C (major histocompatibility complex class I, C) — splicing, lead PIP
1.0, single-variant CS; PE GWAS.** HLA-C on extravillous trophoblasts is the
dominant ligand for killer-cell immunoglobulin-like receptors (KIR) on uterine
NK cells, and this interaction calibrates trophoblast invasion and
spiral-artery remodeling [35, 36]; specific HLA-C/KIR combinations alter
preeclampsia risk. A resolved splicing QTL at HLA-C adds a regulatory layer to
a locus usually considered through protein-coding and copy-number variation.

**ACVR2A (activin A receptor type 2A) — splicing, lead PIP 1.0, single-variant
CS; PE GWAS candidate.** ACVR2A was nominated by linkage and association as a
preeclampsia positional candidate at 2q22 [37], though replication has been
population-dependent [38]. Functionally, ACVR2A promotes trophoblast migration
and invasion via TCF7/c-JUN and is downregulated in preeclamptic placentas
[39]. The fine-mapped splicing signal offers a candidate causal mechanism
reconciling the genetic association with the observed expression changes.

**HSD11B2 (hydroxysteroid 11-beta dehydrogenase 2) — splicing, lead PIP 1.0,
single-variant CS; STB marker.** HSD11B2 forms the placental glucocorticoid
barrier, inactivating maternal cortisol to protect the fetus [40]; its
placental expression is epigenetically tuned by maternal distress [41] and its
variants are associated with hypertensive pregnancy [42]. A single-variant
splicing credible set identifies a regulatory lever on a barrier central to
fetal programming.

**IL1R1 (interleukin 1 receptor type 1) — 4 modalities, lead PIP 1.0,
single-variant CS.** IL1R1 transduces IL-1 signaling, a core inflammatory
pathway in parturition: IL-1 receptor signaling is necessary for
infection-induced preterm labor in mice [44], IL1R1 is dynamically regulated
across rat parturition [43], and IL1R1 polymorphisms associate with
chorioamnionitis risk [45]. Its cross-modality signal (expression, isoform
expression, RNA editing, splicing) implicates regulated IL-1 responsiveness in
the timing of labor.

### 5.3 Intragenomic conflict {-}

The six imprinted genes below are the clearest molecular footprint of
maternal-fetal conflict in the catalog [1, 2, 3]. Because imprinting enforces
monoallelic expression, a regulatory variant at these loci acts on effectively
the only expressed copy, so parent-of-origin-specific effects on placental and
fetal growth are expected [2, 48].

**IGF2 (insulin-like growth factor 2) — splicing, lead PIP 1.0, single-variant
CS; imprinted (paternally expressed).** IGF2 is the canonical
growth-promoting imprinted gene: placental-endocrine Igf2 manipulates maternal
metabolism to partition nutrients to the fetus [46], and the IGF2-H19 complex
matches placental nutrient supply to fetal demand [47]. Chorionic-villus IGF2
expression correlates with fetal growth [48]. A resolved splicing QTL at IGF2
is a direct regulatory handle on the fetal side of the conflict.

**IGF2R (IGF2 receptor) — expression, lead PIP 1.0, single-variant CS;
imprinted (maternally expressed).** IGF2R is the functional antagonist of
IGF2 — it sequesters and degrades IGF2, restraining growth — and the
IGF2/IGF2R pair was the first described embodiment of parental conflict [49].
A cleanly resolved expression QTL at the maternally expressed receptor
complements the paternally expressed ligand signal above.

**MEG3 (maternally expressed 3) — splicing, lead PIP 1.0, single-variant CS;
imprinted.** MEG3 is a maternally expressed lncRNA in the Dlk1-Dio3 imprinted
domain; the Meg3 transcript is required in cis to repress paternal Dlk1 [50],
and its differentially methylated region is the imprinting-control element of
the locus, with loss causing embryonic death and placental loss of imprinting
[51]. A splicing credible set at MEG3 points to post-transcriptional control
of an imprinting-control lncRNA.

**PEG3 (paternally expressed 3) — splicing, lead PIP 1.0, single-variant CS;
imprinted.** PEG3 anchors a ~500-kb imprinted domain controlling fetal growth
and maternal-caring behavior [52], and imprinted DMRs (including this
pathway's H19/MEST) are sensitive to the conception environment [53]. The
resolved splicing signal adds PEG3 to the set of conflict loci under genetic
regulatory control in term placenta.

**MEST (mesoderm-specific transcript) — splicing, lead PIP 1.0, single-variant
CS; imprinted (paternally expressed).** Loss of paternally expressed Mest
causes embryonic growth retardation and abnormal maternal behavior in mice
[54], placental MEST is downregulated in early-onset preeclampsia [55], and
its placenta-conserved imprinting across mammals supports a conflict-driven
origin [56]. MEST thus links the conflict lens directly back to pregnancy
disease.

**GNAS (GNAS complex locus) — splicing, lead PIP 1.0, single-variant CS;
imprinted.** GNAS is a complex imprinted locus with oppositely imprinted
transcripts (Gs-alpha vs XL-alpha-s) whose parent-of-origin effects on growth
and metabolism follow the conflict model [57]; its disruption causes human
imprinting disorders [58], and altered GNAS NESP55 methylation is seen in
placental overgrowth (placental mesenchymal dysplasia) [59]. A single-variant
splicing credible set resolves regulatory variation at this dosage-sensitive
locus.

## 6. Synthesis {-}

Reading the high-confidence fine-mapped genes through the three lenses yields a
coherent picture. The **pregnancy** and **conflict** lenses are dominated by
single-variant **splicing** credible sets at genes with strong prior evidence
(FLT1, PAPPA, HLA-C, ACVR2A, HSD11B2) or established imprinting (IGF2, IGF2R,
MEG3, PEG3, MEST, GNAS) — these serve as positive controls that validate the
fine-mapping and point to post-transcriptional regulation as the operative
mechanism. The **placental-biology** lens, by contrast, surfaces the atlas's
most novel contribution: broadly cross-modality genes (WSB1, GBP3, CAST) with
single-variant resolution but no prior placental annotation, alongside
mechanistically anchored regulators of syncytialization (GCM1, SIGLEC6,
CYP19A1). The convergence of the splicing-QTL enrichment for
syncytiotrophoblast markers (Report 2) with the splicing-dominated pregnancy
and conflict signals here suggests that **regulation of RNA processing in the
syncytiotrophoblast is a principal route by which common genetic variation
shapes placental function and pregnancy outcomes.**

## 7. Limitations {-}

- **Abstract-level synthesis.** Key findings are drawn from abstracts (and, for
  framing papers, review-level knowledge); individual effect sizes and
  subgroup results were not re-extracted from full texts.
- **Selection is illustrative, not exhaustive.** Twenty genes were chosen to
  span the three lenses; the full 264-gene screen
  (`finemap_lit_screen_candidates.tsv`) and the 5,236-gene high-confidence set
  contain additional candidates (e.g., SGCE, ZDBF2, KISS1, the PSG family,
  GLS, PSCA) not given deep-dives here.
- **Lens assignment is not exclusive.** Several genes (e.g., SIGLEC6, MEST,
  CYP19A1) span placental biology and pregnancy disease; each is discussed
  under its dominant lens.
- **Credible-set causality is statistical.** A single-variant credible set at
  lead PIP ~1.0 is strong evidence of a causal regulatory variant but does not
  by itself prove the molecular mechanism or the direction of effect on the
  phenotype.
- **Literature recency.** PubMed and iCite were queried at the time of writing;
  citation counts and the literature landscape will drift.

## 8. Files {-}

- `literature_finemap.md` — this report
- `finemap_hits_evidence_table.csv` — all 59 cited references with
  bibliographic metadata and key-finding notes; the `index` column matches the
  `[N]` citations above
- `finemap_lit_screen_candidates.tsv` — the 264-gene PubMed screen (hit counts
  and top PMIDs per lens) used to select the deep-dive genes
- `report_finemap.html` — Report 2 (fine-mapping results), the companion this
  review supports
