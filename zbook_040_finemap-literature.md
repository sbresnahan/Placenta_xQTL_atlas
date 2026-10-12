# Literature context: fine-mapped loci {#finemap-literature}


**Project:** Placenta xQTL atlas  
**Revision:** 2026-10-10, rebuilt against the newly supplied `report_finemap.html`,
`finemap_locus_summary.tsv`, and `finemap_credible_sets.tsv.gz`.

*Companion literature review to `report_finemap.html`. It is the fine-mapping
counterpart to `literature_xqtl_mapping.md`, which reviews the strongest
variant-level association signals from the xQTL mapping report.*

### 1. Background and scope

The current fine-mapping run uses SuSHiE to resolve significant placental xQTL
loci across nine molecular modalities: gene expression, isoform expression,
isoform usage, splicing, intron retention, alternative transcription start
sites, alternative polyadenylation, RNA editing, and transcript stability.
The run report records **5,529 loci attempted, 5,269 successfully fine-mapped,
and 4,290 loci with at least one credible set**.

this review does **not** retain the literature review's claim of 5,236
high-confidence genes. Recalculation from the supplied results gives
**2,055 current loci whose best credible set has lead PIP ≥ 0.9**, spanning
**1,142 unique Ensembl genes**. Of those high-confidence loci, **1,899
(92.4%) have a single-variant credible set**. Thus the strongest conclusion of
the new run is not the sheer size of the gene catalog but the high degree of
variant-level resolution among loci that pass the PIP threshold.

The biological reading is also different from the literature review. The old
fine-mapping review emphasized a splicing-heavy set of established pregnancy
and imprinted genes. In the new run, many of those loci are absent or no longer
reach PIP ≥ 0.9. Instead, the leading pattern is **cross-modality recurrence**:
WSB1 is high-confidence in six modalities, while SCARB2 and OLAH are
high-confidence in five. Established placental and pregnancy genes remain in
the catalog — including CYP19A1, SIGLEC6, ERAP1, KISS1, LGALS14, PSG3,
PSG9, PECAM1, and HSPD1 — but the conflict/imprinting lens is much smaller,
with IGF2R and SGCE providing the clearest retained examples.

Existing numeric citations `[N]` below preserve the numbering of the original
fine-mapping literature review where the same literature-supported claim is
retained. Newly added references are labeled `[A1]`–`[A18]` and are listed in
Section 9. The companion `finemap_hits_evidence_table.csv` was not included
with the new fine-mapping upload, so the new references are intentionally not
assigned permanent numeric indices yet.

### 2. Methods and result reconciliation

#### 2.1 Current-run membership

The supplied `finemap_credible_sets.tsv.gz` contains **6,347 credible-set
rows**, matching the total displayed in the HTML report. However, **763 of
those rows (534 phenotype IDs)** do not occur in the newly supplied
`finemap_locus_summary.tsv`. The 4,290 current loci with credible sets in the
locus summary account for **5,584 credible sets**.

To prevent carry-over results from a previous mapping/fine-mapping run from
entering the literature interpretation, all gene rankings in this review use
only credible sets whose `(modality, phenotype_id)` pair is present in the new
locus summary. This filtering step should also be incorporated into the report
generation code, or the credible-set file should be regenerated cleanly.

#### 2.2 High-confidence definition

For each current `(modality, phenotype_id)` locus, the credible set with the
highest `lead_pip` was selected. A locus is called high-confidence when that
best credible set has `lead_pip >= 0.9`. Gene IDs were taken from the Ensembl
prefix of `phenotype_id`. This gives:

- **2,055 high-confidence loci**
- **1,142 unique genes**
- **1,899 / 2,055 (92.4%) single-variant credible sets**
- **405 genes** high-confidence in at least two modalities
- **79 genes** in at least three modalities
- **25 genes** in at least four modalities
- only **3 genes** in at least five modalities: WSB1, SCARB2, and OLAH.

| Modality | HC loci (lead PIP ≥0.9) | Unique genes | Single-variant HC loci |
|---|---:|---:|---:|
| RNA_editing | 136 | 46 | 123 |
| alt_TSS | 63 | 42 | 60 |
| alt_polyA | 213 | 110 | 196 |
| expression | 640 | 640 | 593 |
| intron_retention | 48 | 37 | 45 |
| isoform_expression | 714 | 570 | 669 |
| isoforms | 132 | 104 | 115 |
| splicing | 24 | 21 | 23 |
| stability | 85 | 85 | 75 |

#### 2.3 Deep-dive selection

Twenty genes were selected to represent the strongest new themes: recurrent
multi-modality regulation, established placental biology/pregnancy disease,
and the smaller retained imprinting/conflict signal. Every gene below has at
least one **current** lead-PIP ≥0.9 credible set and a minimum credible-set size
of one variant.

| Gene | Lens | Modalities | Current HC loci | Max lead PIP | Min CS size |
|---|---|---:|---:|---:|---:|
| WSB1 | Recurrent/novel | 6 | 7 | 1.00 | 1 |
| SCARB2 | Recurrent/placental | 5 | 8 | 1.00 | 1 |
| OLAH | Recurrent/pregnancy | 5 | 6 | 1.00 | 1 |
| CTSB | Recurrent/placental | 4 | 27 | 1.00 | 1 |
| GPNMB | Recurrent/placental | 4 | 12 | 1.00 | 1 |
| FOLR1 | Recurrent/placental | 4 | 5 | 1.00 | 1 |
| GBP3 | Recurrent/novel | 4 | 7 | 1.00 | 1 |
| HSP90B1 | Placental stress | 3 | 4 | 1.00 | 1 |
| CYP19A1 | Placental/pregnancy | 3 | 6 | 1.00 | 1 |
| SIGLEC6 | Placental/pregnancy | 2 | 6 | 1.00 | 1 |
| ERAP1 | Pregnancy | 4 | 10 | 1.00 | 1 |
| IL1R1 | Pregnancy | 3 | 7 | 0.999 | 1 |
| KISS1 | Pregnancy | 2 | 2 | 1.00 | 1 |
| LGALS14 | Placental/pregnancy | 3 | 4 | 1.00 | 1 |
| PSG3 | Placental | 3 | 7 | 1.00 | 1 |
| PSG9 | Placental/pregnancy | 2 | 6 | 1.00 | 1 |
| PECAM1 | Placental vascular | 2 | 17 | 1.00 | 1 |
| HSPD1 | Placental stress | 1 | 1 | 0.994 | 1 |
| IGF2R | Conflict/imprinting | 1 | 1 | 0.999 | 1 |
| SGCE | Conflict/imprinting | 1 | 1 | 1.00 | 1 |

### 3. Cross-cutting patterns

#### 3.1 Multi-modality recurrence replaces splicing dominance as the main signal

The earlier review identified WSB1 and GBP3 as the most recurrent genes and
then interpreted much of the pregnancy/conflict catalog through splicing. The
new run retains WSB1 but changes the ranking substantially. **WSB1 now spans
six modalities**, followed by **SCARB2 and OLAH at five**. Twenty-two additional
genes span four modalities. The recurrent genes are therefore distributed
across expression, isoform-level traits, RNA editing, RNA stability, and
alternative ends rather than being dominated by one RNA-processing mechanism.

WSB1 remains the clearest "genetics first, biology second" discovery. It has
high-confidence signals in RNA editing, expression, intron retention, isoform
expression, isoform usage, and stability. Its known biology — HIF-responsive
ubiquitin regulation and functionally distinct splice isoforms [13,14] — makes
the recurrence plausible in a hypoxia-responsive organ, but a direct placental
role remains poorly characterized.

SCARB2 and OLAH now join WSB1 as the two other five-modality genes. SCARB2 is a
lysosomal membrane protein with strong protein expression in
syncytiotrophoblast and Hofbauer cells `[A1]`, while OLAH has unusually direct
pregnancy relevance: placental OLAH changes across gestation and is altered in
preeclampsia and fetal growth restriction; OLAH knockdown changes IGF2,
oxidative-stress, and apoptosis-related transcription `[A2]`.

#### 3.2 The high-confidence catalog is extremely well resolved

A single-variant credible set is the minimum possible 95% CS. **More than 92%
of current PIP ≥0.9 loci are single-variant**, and all 20 deep-dive genes have
at least one such locus. This makes fine-mapping useful for variant
prioritization, but the statistical precision should not be confused with
mechanistic proof: a one-variant CS identifies a strong causal-regulatory
candidate, not the molecular route by which it acts.

#### 3.3 Established placental genes remain, but the positive controls changed

The current high-confidence catalog still contains biologically reassuring
placental genes, but not the exact set highlighted previously. CYP19A1,
SIGLEC6, ERAP1, KISS1, LGALS14, PSG3, PSG9, PECAM1, HSPD1 and HSP90B1 all
retain at least one single-variant PIP ≥0.9 signal. By contrast, **FLT1, PAPPA,
ACVR2A, HSD11B2, IGF2, PEG3, MEST and GNAS are not current high-confidence
genes**, and GCM1, HLA-C and MEG3 are present only below the 0.9 threshold.

The implication is important: the literature review's statement that fine-mapping
"recovered known pregnancy genes as positive controls" was too dependent on an
older result set. The new run should be interpreted from the loci it actually
resolves, not from biologically attractive loci that were present previously.

#### 3.4 Imprinting remains visible, but it is no longer a dominant catalog-wide theme

The literature review highlighted six canonical imprinted genes and stated that
27 imprinted genes were present in the high-confidence set. That claim is not
supported by the current files supplied here. Among the original six, **only
IGF2R remains PIP ≥0.9**. A second clear current imprinted signal is **SGCE**,
which has a single-variant isoform-expression credible set. SGCE is maternally
imprinted `[A17]`, is expressed across human gestation `[A16]`, and together
with IGF2R has been implicated in environmentally responsive placental
imprinted-gene expression `[A18]`.

Thus intragenomic conflict remains a useful biological lens, but it is now a
focused secondary theme rather than one of the dominant statistical patterns.

### 4. Gene snapshots

#### 4.1 Recurrent and placental-biology signals

**WSB1 — 6 modalities; 7 high-confidence loci; max lead PIP 1.00; min CS = 1.**  
WSB1 is the most recurrent gene in the current fine-mapped catalog. The
resolved modalities are RNA editing, expression, intron retention, isoform
expression, isoform usage, and stability. WSB1 encodes the
substrate-recognition component of an E3 ubiquitin ligase and participates in
hypoxia signaling [13]. Its splice isoforms can have distinct effects on
proliferation and apoptosis [14], making the breadth of RNA-level regulation
biologically credible. Direct placental characterization remains limited,
which keeps WSB1 among the atlas's strongest novel candidates.

**SCARB2 — 5 modalities; 8 high-confidence loci; max lead PIP 1.00; min CS = 1.**  
SCARB2/LIMP-2 is a lysosomal membrane receptor. Unlike WSB1, its placental
presence is directly supported by tissue proteomics: SCARB2 is high in
syncytiotrophoblast cell body and microvilli and in Hofbauer cells `[A1]`.
Fine-mapping across RNA editing, alternative TSS, expression, isoform
expression and isoform usage therefore nominates a lysosome/membrane-trafficking
axis operating in both trophoblast and placental macrophage compartments.
Mechanistic pregnancy literature is still sparse.

**OLAH — 5 modalities; 6 high-confidence loci; max lead PIP 1.00; min CS = 1.**  
OLAH is the strongest newly validated pregnancy-relevant recurrent signal. A
human placental study found that OLAH expression rises with gestation, is
altered in preterm preeclampsia and fetal growth restriction, and responds to
placental dysfunction; OLAH knockdown perturbed BAX, BCL2, IGF2 and NOX4
expression `[A2]`. The current fine-map spans RNA editing, alt-TSS, expression,
isoform expression and stability, converting that disease-expression
observation into a candidate genetically regulated placental pathway.

**CTSB — 4 modalities; 27 high-confidence loci; max lead PIP 1.00; min CS = 1.**  
CTSB has the largest number of high-confidence loci among the selected genes,
distributed across RNA editing, alt-polyA, expression and isoform expression.
Cathepsin B protein is high in placental trophoblast `[A3]`; experimental
trophoblast work places CTSB downstream of beta-catenin-dependent migration
and invasion in a hypoxia/preeclampsia context `[A4]`. The striking number of
resolved CTSB molecular phenotypes suggests that this locus is not a single
expression QTL but a dense regulatory hotspot.

**GPNMB — 4 modalities; 12 high-confidence loci; max lead PIP 1.00; min CS = 1.**  
GPNMB is resolved across alt-polyA, expression, isoform expression and isoform
usage. Placental single-cell reference data show particularly strong
expression in Hofbauer cells `[A5]`, pointing to placental macrophages rather
than trophoblast as the most likely cellular context. This is a useful example
of why cell-type attribution should accompany bulk-placenta QTL interpretation.

**FOLR1 — 4 modalities; 5 high-confidence loci; max lead PIP 0.99995; min CS = 1.**  
FOLR1 encodes folate receptor-alpha, a core component of maternal-fetal folate
transfer. FR-alpha localizes to the apical syncytiotrophoblast `[A7]`, and
placental FOLR1 expression/methylation has been studied in preterm birth
`[A6]`. The four-modality fine-map — alt-TSS, expression, isoform expression
and isoform usage — identifies several genetically controlled points at which
placental folate transport could vary among individuals.

**GBP3 — 4 modalities; 7 high-confidence loci; max lead PIP 1.00; min CS = 1.**  
GBP3 remains a recurrent innate-immunity signal, although its current
recurrence is four rather than six modalities: expression, isoform expression,
isoform usage and stability. GBPs are interferon-inducible antimicrobial
GTPases [15], and GBP3 can participate in caspase-4-dependent pyroptosis [16].
There is still little direct placental literature, so GBP3 remains a
high-confidence but biologically under-characterized maternal-fetal-interface
candidate.

**HSP90B1 — 3 modalities; 4 high-confidence loci; max lead PIP 0.99976; min CS = 1.**  
HSP90B1/GRP94 is resolved at alternative polyadenylation, expression and
isoform expression. Human Protein Atlas staining is high in both trophoblastic
and decidual cells `[A8]`. The locus therefore adds genetically regulated
endoplasmic-reticulum stress/chaperone biology to the same broad stress axis
represented by HSPD1 below.

#### 4.2 Pregnancy and trophoblast-function signals

**CYP19A1 — 3 modalities; 6 high-confidence loci; max lead PIP 1.00; min CS = 1.**  
CYP19A1 remains a strong placental positive control, although the new run
supports three modalities rather than five: alt-polyA, intron retention and
isoform expression. Aromatase drives placental estrogen synthesis and has
trophoblast-specific promoter biology [17]; placental CYP activity is altered
in preeclampsia-related models [19]. Fine-mapping indicates that genetically
regulated transcript architecture, not only total expression, is relevant at
this canonical placental gene.

**SIGLEC6 — 2 modalities; 6 high-confidence loci; max lead PIP 1.00; min CS = 1.**  
SIGLEC6 remains one of the most biologically specific placental hits. It is a
placenta-enriched molecule and circulating biomarker candidate in
preeclampsia [21], with mechanistic evidence connecting excess SIGLEC6 to
impaired syncytialization and inflammatory polarization [22,23]. The current
fine-map resolves RNA-editing and expression loci rather than the four
modalities stated previously.

**ERAP1 — 4 modalities; 10 high-confidence loci; max lead PIP 1.00; min CS = 1.**  
ERAP1 remains one of the most compelling pregnancy-disease loci in the
catalog. It participates in MHC-I peptide processing and vasoactive peptide
metabolism, has placental and preeclampsia evidence [28,29], and sits directly
at the trophoblast-immune interface [30]. The current resolved modalities are
RNA editing, intron retention, isoform expression and stability. This still
supports unusually broad regulatory control, but the old expression/splicing
claims should not be carried forward.

**IL1R1 — 3 modalities; 7 high-confidence loci; max lead PIP 0.9986; min CS = 1.**  
IL1R1 remains high-confidence across RNA editing, expression and isoform
expression. IL-1 receptor signaling is a well-established inflammatory route
into parturition and infection-induced preterm labor [43,44], and IL1R1
variation has been associated with chorioamnionitis risk [45]. The current
fine-map therefore supports genetically regulated IL-1 responsiveness without
requiring a splicing mechanism.

**KISS1 — 2 modalities; 2 high-confidence loci; max lead PIP ~1.00; min CS = 1.**  
KISS1 is resolved for expression and isoform expression. Kisspeptin is tightly
linked to trophoblast invasion: placental KISS1 is altered in preeclampsia,
and experimental overexpression suppresses trophoblast invasion `[A9,A10]`.
The QTL therefore provides a direct genetic route into a pathway already
implicated in defective placentation.

**LGALS14 — 3 modalities; 4 high-confidence loci; max lead PIP 1.00; min CS = 1.**  
LGALS14 is a primate placenta-associated galectin resolved for expression,
isoform expression and stability. Galectin-14 promotes trophoblast migration
and invasion and is reduced in early pregnancy loss and preeclampsia `[A11]`.
The combination of placental specificity, disease dysregulation and
single-variant QTL resolution makes LGALS14 a particularly interpretable
fine-mapped locus.

**PSG3 — 3 modalities; 7 high-confidence loci; max lead PIP 1.00; min CS = 1.**  
PSG3 is resolved for intron retention, isoform expression and isoform usage.
Pregnancy-specific glycoproteins are produced in large amounts by placental
trophoblast; PSG3 is among the genes strongly induced during trophoblast
differentiation `[A12]`. The multi-layer signal suggests genetic control over
both PSG3 abundance and transcript architecture.

**PSG9 — 2 modalities; 6 high-confidence loci; max lead PIP 1.00; min CS = 1.**  
PSG9 is resolved for intron retention and isoform expression. Recent work
identifies an LTR8B/PSG9 regulatory program controlling syncytialization and
links elevated PSG9 to preeclampsia `[A13]`. The fine-mapped locus therefore
intersects a newly defined, placenta-evolved regulatory mechanism.

**PECAM1 — 2 modalities; 17 high-confidence loci; max lead PIP 1.00; min CS = 1.**  
PECAM1 has a large number of resolved RNA-editing and isoform-usage loci.
Its placental role should be framed carefully: PECAM1 is a strong endothelial
marker, but direct placental-bed studies did not support trophoblast PECAM1 as
a driver of invasion or preeclampsia `[A14]`. The atlas signal may therefore
reflect placental vascular/endothelial regulation rather than a trophoblast
mechanism.

**HSPD1 — RNA editing; 1 high-confidence locus; lead PIP 0.994; CS = 1.**  
HSPD1/Hsp60 is a lower-recurrence but biologically strong stress candidate.
Recent work shows aberrant Hsp60 lactylation in preeclamptic trophoblasts,
promoting mitochondrial fission and apoptosis `[A15]`. A single-variant
RNA-editing QTL provides a distinct post-transcriptional layer on this
mitochondrial stress pathway.

#### 4.3 Intragenomic conflict and imprinting

**IGF2R — expression; 1 high-confidence locus; lead PIP 0.9992; CS = 1.**  
IGF2R is the only one of the literature review's six canonical conflict genes
that remains high-confidence in the new run. It antagonizes IGF2 signaling
and is a classic parent-of-origin growth-regulation locus [49]. Its clean
single-variant expression credible set remains a compelling molecular
connection between cis regulation and placental growth conflict.

**SGCE — isoform expression; 1 high-confidence locus; lead PIP ~1.00; CS = 1.**  
SGCE is a new conflict-lens replacement for the prior IGF2/MEG3/PEG3/MEST/GNAS
set. SGCE is maternally imprinted `[A17]`, is expressed throughout human
gestation `[A16]`, and appears among placental imprinted genes whose expression
tracks maternal one-carbon nutritional status `[A18]`. The resolved
isoform-expression QTL therefore provides a current-run example of genetic
regulation operating within an imprinted placental domain.

### 6. Synthesis

The fine-mapping catalog supports three conclusions.

First, **fine-mapping resolution is exceptionally sharp**: more than 92% of
PIP≥0.9 current loci have one-variant credible sets. This creates a strong
variant-prioritization resource even where biological annotation is sparse.

Second, the strongest recurring biology is **multi-layer RNA regulation rather
than splicing alone**. WSB1, SCARB2, OLAH, CTSB, GPNMB, FOLR1, GBP3 and ERAP1
recur across multiple molecular readouts, often mixing transcript abundance,
isoform regulation, RNA editing, intron retention, stability, and alternative
ends. This is precisely the kind of regulatory architecture that a
multi-modality xQTL atlas can reveal but a standard eQTL study cannot.

Third, the new run retains a coherent set of **placental and pregnancy
positive controls** — CYP19A1, SIGLEC6, KISS1, LGALS14, PSG3/PSG9, PECAM1,
HSPD1 and ERAP1 — while simultaneously weakening the prior claim that the
fine-mapped catalog is dominated by canonical pregnancy GWAS and imprinted
genes. The scientifically stronger interpretation is therefore that the atlas
recovers known placental pathways **and** highlights recurrent loci whose
placental role is incompletely characterized.

### 7. Limitations

- **Current-file mismatch.** The credible-set file contains 534 loci absent
  from the new locus summary. This review filters them out, but the report
  generation pipeline should be corrected so all fine-mapping products share
  the same run membership.
- **PIP threshold is a statistical filter.** Lead PIP ≥0.9 and a one-variant
  credible set are strong causal-regulatory evidence, not proof of molecular
  mechanism.
- **Gene-level recurrence is not locus independence.** A gene with many
  high-confidence molecular phenotypes can contain correlated features driven
  by the same underlying variant/haplotype.
- **Bulk placenta mixes cell types.** GPNMB, SCARB2 and PECAM1 illustrate why
  trophoblast, Hofbauer-cell and endothelial signals must be distinguished
  using cell-type-resolved annotations.
- **Literature selection is illustrative.** The 20 deep dives prioritize
  recurrence and placental interpretability; they are not the only biologically
  important genes among the 1,142 current high-confidence genes.
- **Citation integration is incomplete.** The original
  `finemap_hits_evidence_table.csv` was not supplied with this review. Legacy
  numeric citation indices are retained where possible; new sources remain
  `[A#]` until that evidence table is provided.

### 8. Companion files

- `literature_finemap_revised_20261010.md` — this literature review.
- `finemap_current_high_confidence_genes_20261010.tsv` — recalculated
  PIP≥0.9 gene summary using only loci present in the current locus summary.
- `finemap_old_deep_dive_audit_20261010.tsv` — status of each of the prior
  20 deep-dive genes in the current run.
- `finemap_literature_revision_sources_20261010.tsv` — new `[A1]`–`[A18]`
  literature/reference sources.
- `finemap_locus_summary.tsv` and `finemap_credible_sets.tsv.gz` — supplied
  current fine-mapping outputs.
- `report_finemap.html` — supplied current fine-mapping report.

### 9. New literature/reference addendum

The entries below were added during this review. They should be appended to
the original `finemap_hits_evidence_table.csv` once that file is available,
preserving its existing numeric indices and assigning these references new
indices after the current maximum.

- **[A1]** Human Protein Atlas: SCARB2 staining in placenta. SCARB2 protein is high in syncytiotrophoblast cell body/microvilli and Hofbauer cells in placenta. https://www.proteinatlas.org/ENSG00000138760-SCARB2/tissue/placenta
- **[A2]** Placental OLAH Levels Are Altered in Fetal Growth Restriction, Preeclampsia and Models of Placental Dysfunction (2022); PMID: 36139751. Placental OLAH rises with gestation, is altered in preterm preeclampsia/FGR, and OLAH knockdown changes apoptosis, IGF2 and oxidative-stress genes. https://pubmed.ncbi.nlm.nih.gov/36139751/
- **[A3]** Human Protein Atlas: CTSB staining in placenta. CTSB protein is high in trophoblastic cells and medium in decidual cells. https://v19.proteinatlas.org/ENSG00000164733-CTSB/tissue/placenta
- **[A4]** Effect of hypoxia on expression of placental trophoblast cells SATB1 and beta-catenin and its correlation with the pathogenesis of preeclampsia (2016). Beta-catenin perturbation reduces trophoblast migration/invasion and CTSB expression in a preeclampsia-related hypoxia model. https://www.sciencedirect.com/science/article/pii/S1995764516300542
- **[A5]** Human Protein Atlas: GPNMB single-cell expression. GPNMB is strongly expressed in placental Hofbauer cells, supporting a placental macrophage interpretation. https://www.proteinatlas.org/ENSG00000136235-GPNMB/single+cell
- **[A6]** Epigenetic regulation of folate receptor-alpha (FOLR1) in human placenta of preterm newborns (2020); PMID: 32421530. FOLR1 is a placental folate transporter; placental expression and methylation were studied in preterm birth. https://pubmed.ncbi.nlm.nih.gov/32421530/
- **[A7]** Recent Advances in Folates and Autoantibodies against Folate Receptors in Early Pregnancy and Miscarriage (2023). FR-alpha/FOLR1 localizes to the apical syncytiotrophoblast and contributes to maternal-fetal folate transfer. https://pmc.ncbi.nlm.nih.gov/articles/PMC10708193/
- **[A8]** Human Protein Atlas: HSP90B1 staining in placenta. HSP90B1 is highly expressed in trophoblastic and decidual cells. https://www.proteinatlas.org/ENSG00000166598-HSP90B1/tissue/placenta
- **[A9]** Elevated expression of KiSS-1 in placenta of preeclampsia and its effect on trophoblast (2011); PMID: 21804632. KISS1 was increased in preeclamptic placentas and overexpression reduced trophoblast invasion. https://pubmed.ncbi.nlm.nih.gov/21804632/
- **[A10]** Altered placental expression of kisspeptin and its receptor in pre-eclampsia (2012); PMID: 22526494. KISS1/KISS1R are expressed in placenta and altered in preeclampsia, implicating invasion/angiogenesis. https://pubmed.ncbi.nlm.nih.gov/22526494/
- **[A11]** Galectin-14 Promotes Trophoblast Migration and Invasion by Upregulating the Expression of MMP-9 and N-Cadherin (2021). LGALS14 is placenta/trophoblast-specific, promotes trophoblast migration/invasion, and is reduced in preeclampsia/early pregnancy loss. https://pmc.ncbi.nlm.nih.gov/articles/PMC8007908/
- **[A12]** Expression and transcriptional regulation of individual pregnancy-specific glycoprotein genes in differentiating trophoblast cells (2010); PMID: 20116096. PSG1/3/5 rise strongly during trophoblast differentiation; PSG transcripts/proteins are early differentiation markers. https://pubmed.ncbi.nlm.nih.gov/20116096/
- **[A13]** Endogenous retroviral elements LTR8B and MER65 rewire PSG9 regulation to control trophoblast syncytialization and pre-eclampsia risk (2026). LTR8B/PSG9 regulates syncytialization and PSG9 is strongly elevated in preeclampsia. https://pmc.ncbi.nlm.nih.gov/articles/PMC12969887/
- **[A14]** Human trophoblast invasion and spiral artery transformation: the role of PECAM-1 in normal pregnancy, preeclampsia, and fetal growth restriction (2001); PMID: 11337369. PECAM1 marks placental endothelium but was not detected in trophoblast; the study did not support a direct trophoblast-invasion role. https://pubmed.ncbi.nlm.nih.gov/11337369/
- **[A15]** Hsp60 lactylation promotes mitochondrial dysfunction and trophoblast apoptosis in preeclampsia (2025); DOI: 10.1016/j.bbrc.2025.152379. HSPD1/Hsp60 is abundant in trophoblast and aberrant Hsp60 lactylation promotes mitochondrial dysfunction and apoptosis in preeclampsia. https://www.sciencedirect.com/science/article/abs/pii/S0006291X25010940
- **[A16]** Temporal regulation of the expression of syncytin (HERV-W), maternally imprinted PEG10, and SGCE in human placenta (2003); PMID: 12620933. SGCE is expressed in human placenta and shows gestational regulation distinct from neighboring PEG10. https://pubmed.ncbi.nlm.nih.gov/12620933/
- **[A17]** The epsilon-sarcoglycan gene (SGCE), mutated in myoclonus-dystonia syndrome, is maternally imprinted (2003); PMID: 12634861. SGCE shows parent-of-origin-specific methylation/expression consistent with maternal imprinting. https://pubmed.ncbi.nlm.nih.gov/12634861/
- **[A18]** The roles of maternal one-carbon metabolism and placental imprinted gene expression in placental development and somatic growth in a longitudinal birth cohort (2025). Placental expression of imprinted SGCE and IGF2R was associated with maternal choline status. https://www.sciencedirect.com/science/article/abs/pii/S0143400425001286

