# Literature context: colocalization findings {#coloc-literature}

<p class="source-note">Book chapter rebuilt locally from <code>reports/literature_colocalization.md</code>.</p>


**Project:** Placenta xQTL atlas  
**Companion report:** `report_colocalization.Rmd` / rendered `report_colocalization.html`  
**Results basis:** `coloc_best.tsv.gz`, `coloc_results.tsv.gz`, `coloc_diagnostics_summary.tsv`, and `colocboost_clusters.tsv.gz` from the supplied colocalization aggregate.

Bracketed references **[1]–[14]** resolve to `coloc_hits_evidence_table.csv`, retained from the archived report convention. References **[15]–[17]** below cover the colocBoost additions and new interpretation.

### 1. Scope and interpretation

The colocalization analysis asks whether a placental molecular QTL and an
external GWAS association are statistically consistent with a **shared causal
variant** rather than merely occupying the same genomic neighborhood. The
pipeline uses the SuSiE extension of coloc, which permits multiple causal
signals within a locus and evaluates pairs of resolved signals rather than
assuming one causal variant per region [1].

The supplied nonredundant `coloc_best.tsv.gz` contains **3,192 xQTL-feature ×
GWAS-trait × ancestry comparison cells**. **194 are called colocalized** by the
supplied flag (which corresponds to PP4 ≥0.70 in these files), with **167 at
PP4 ≥0.80, 115 at PP4 ≥0.90, and 87 at PP4 ≥0.95**. The called set spans **139
molecular phenotypes, 66 primary Ensembl gene IDs, and 15 GWAS traits**. There
are **92 EAS and 102 EUR calls**.

The trait panel covers four related but genetically distinct early-growth
settings. The 2019 EGG analyses distinguish direct fetal, indirect maternal,
and own-birth-weight genetic effects [2]. The age-10 body-size GWAS is a UK
Biobank measure of comparative childhood body size [3]. The JECS data provide
dense early-childhood anthropometric GWAS in a Japanese birth cohort [4], and
the maternal gestational-duration GWAS comes from the 2023 EGG meta-analysis
[5]. These distinctions matter: a placental xQTL that colocalizes with fetal
birth weight and one that colocalizes with a maternal birth-weight effect do
not have the same causal interpretation.

This review focuses on **loci**, not raw row counts. Several molecular
phenotypes from one gene—and several neighboring genes—can be driven by the
same local haplotype. Repeated colocalization across transcripts, poly(A)
sites, editing sites, or modalities is therefore evidence of a coherent
regulatory hotspot, not proof of multiple independent causal mechanisms.

### 2. Overview of the result landscape

The 194 calls are concentrated in **alternative polyadenylation (53),
expression (44), and isoform expression (43)**, followed by RNA editing (16),
stability (12), isoform usage (11), intron retention (9), and alternative TSS
(6). No splicing comparison is called in the current aggregate. Because the
diagnostic funnel contains very few successfully resolved splicing
comparisons, this absence should not be interpreted as a biological exclusion
of splicing.

The GWAS side is dominated by developmental anthropometric phenotypes:
comparative body size at age 10 contributes 70 calls; fetal/direct birth
weight contributes 24; several JECS BMI timepoints contribute substantial
additional counts; own birth weight contributes 11; and the maternal
birth-weight and maternal gestational-duration resources contribute smaller
sets.

Only **31/194 called comparisons use the identical xQTL and GWAS lead variant**.
This is expected under LD-aware colocalization: H4 concerns a shared causal
signal, not equality of the marginal lead SNP.

The most important negative result is cross-ancestry replication. **No exact
molecular-feature–trait call, and no primary gene–trait call, reaches the
colocalization criterion in both EAS and EUR.** This makes the current atlas a
set of ancestry-stratified mechanistic hypotheses rather than internally
replicated trans-ancestry colocalizations.

### 2.1 What colocBoost adds to the pairwise catalog

The archived report was based on pairwise SuSiE-coloc. The new `colocboost_clusters.tsv.gz` adds a different view: **joint multi-trait colocalization**. ColocBoost uses a multi-task gradient-boosting framework that can accommodate multiple causal variants and discover configurations of several outcomes sharing one effect. Each event is represented by a colocalization confidence set (CoS) of candidate variants and a set of colocalized outcomes [15,16].

The supplied colocBoost aggregate contains **564 CoS rows**. Of these, **12 events jointly contain at least one placental xQTL and at least one GWAS outcome**. Those events yield **22 ancestry-matched gene–trait pairs involving 13 genes and 12 GWAS traits**. **Seven of the 22 pairs are also called by pairwise SuSiE-coloc in the same ancestry; 15 are colocBoost-only.** This is not an “agreement rate” in the conventional sense because the methods have different targets and signal decompositions. The overlap is best used as a high-confidence convergence tier; colocBoost-only events are complementary multi-trait hypotheses.

The 12 mixed events have high LD coherence overall (median CoS purity ~0.974), but **purity is not a posterior probability**. In the colocBoost output, `purity` is the minimum absolute correlation among variants in the CoS, `top_variable_vcp` is the top variant's colocalization probability, `cos_npc` summarizes normalized colocalization support for the event, and `min_npc_outcome` is the weakest normalized support among its outcomes [16].

Two events are especially informative for developmental trajectories. A **SH3YL1** EAS event jointly links seven placental xQTL outcomes across expression/isoform/poly(A)/stability layers with **four JECS BMI GWAS outcomes** (Dr1m, C1Y_1, C3Y and C4Y). A **PLCE1-AS1** EAS event links isoform/poly(A) regulation with **five JECS BMI outcomes** (Dr1m, C6m, C1Y_1, C1Y_2 and C3HY). The appropriate interpretation is a shared regulatory configuration across several growth timepoints, not four or five independent causal variants.

Other mixed events implicate **FAP, C3orf38, FES, FUBP1, GCLM, MRPL18/TCP1, BCKDHB, RBIS, ARF4, and NPEPL1**. These provide a compact shortlist for follow-up because the xQTL and GWAS outcomes were selected jointly rather than by repeated pairwise scanning.

### 3. The chr17p13 birth-weight hotspot: ELP5, CLDN7 and CTDNEP1

The strongest and best externally validated locus in the current results lies
at **17p13**, where ELP5, CLDN7 and CTDNEP1 sit within a compact regulatory
region.

#### 3.1 ELP5: the dominant molecular-feature signal

ELP5 accounts for **29 called comparison rows**, the largest count for any
primary gene ID. In EUR, the fetal/direct birth-weight GWAS colocalizes with
ELP5 expression, isoform expression, isoform usage, and multiple
alternative-polyadenylation phenotypes, many with PP4 ~0.995–1.00. Several of
these comparisons use the same xQTL and GWAS lead variant
`17:7259132:G:A`. In EAS, ELP5 poly(A), isoform-expression, and isoform-usage
signals colocalize with trans-ancestry own birth weight, JECS low birth weight,
and neonatal/early-childhood BMI codes.

ELP5 itself is a component of the Elongator complex and is annotated for tRNA
modification and cell-migration biology [8]. It does not yet have a strong
placenta-specific disease literature. Consequently, the strength of the
ELP5 signal comes primarily from **genetic convergence across molecular
modalities and growth traits**, not from prior mechanistic expectation.

#### 3.2 CLDN7: an independently recognizable birth-weight region

CLDN7 has four EUR fetal-birth-weight calls across expression and
alternative-polyadenylation phenotypes. The strongest has PP4 ~0.9998, and the
same `17:7259132:G:A` variant is the xQTL and GWAS lead in several comparisons.
CLDN7 encodes a tight-junction claudin at 17p13.1 [14].

Importantly, the region is not novel to birth-weight genetics. Fine-mapping in
the Black Women's Health Study previously evaluated a **SLC2A4/CLDN7 17p13
birth-weight locus** [7]. Thus the present placental xQTL colocalization places
a molecular regulatory readout directly under a previously recognized
birth-weight association region.

#### 3.3 CTDNEP1: direct replication of prior placental multi-omics

CTDNEP1 provides the strongest external validation of the atlas's biological
interpretation. In the current results, an EAS CTDNEP1 intron-retention
phenotype colocalizes with trans-ancestry own birth weight (PP4 ~0.967) and
with JECS Dr0m/Dr1m BMI signals.

A prior placental multi-omics study integrated birth-weight GWAS, placental DNA
methylation, and placental gene expression and independently highlighted
**CTDNEP1 at chr17p13**. That study found evidence for shared genetic
architecture linking birth weight, placental methylation, and CTDNEP1-related
expression at the locus, and identified CTDNEP1 among candidate functional
birth-weight genes [6]. The present atlas therefore does not merely rediscover
a nearby GWAS peak: it extends a previously implicated placental locus into
**RNA-processing regulation**, particularly intron retention.

#### 3.4 Interpretation of the hotspot

ELP5, CLDN7 and CTDNEP1 should not be presented as three independent causal
genes solely because all appear in the call table. Their physical proximity,
shared lead variants, and overlapping trait associations make the more
defensible interpretation a **chr17p13 placental regulatory hotspot**. The
current data nominate several candidate molecular consequences within that
hotspot; distinguishing the true effector gene(s) will require conditional
analyses, variant-level functional annotation, and experimental perturbation.

### 4. WSB1: five molecular routes to fetal birth weight

WSB1 is one of the cleanest multi-modality colocalizations outside chr17p13.
Five EUR molecular phenotypes—RNA editing, intron retention, isoform usage,
transcript stability, and isoform expression—colocalize with the EGG
fetal/direct birth-weight signal at PP4 ~0.997. In every one of these
comparisons, the same lead variant (`17:27310317:G:A`) is selected for the
xQTL and GWAS signal.

WSB1 has a plausible placental stress biology even though it has not been
established as a canonical birth-weight gene. It encodes the
substrate-recognition component of an E3 ubiquitin ligase, is a HIF-1 target,
regulates proteins involved in hypoxic responses, and produces multiple
functional transcripts; earlier tissue profiling reported WSB1 transcripts in
placenta [9]. This makes the locus particularly attractive: the colocalization
links a fetal-growth GWAS signal to **five distinct placental RNA phenotypes**
at a gene already positioned in hypoxia-responsive biology.

The data do not establish which RNA phenotype mediates fetal growth. Rather,
they show that a single fetal-birth-weight locus has unusually coherent effects
across WSB1 RNA processing.

### 5. PLCE1-AS1: a cross-trait developmental-growth locus

PLCE1-AS1 has **11 calls spanning seven traits, three molecular modalities, and
both ancestry strata**, making it the widest cross-trait locus in the current
catalog. EUR signals include fetal/direct birth weight, maternal birth weight,
and comparative body size at age 10; EAS signals include several JECS BMI
timepoints. The molecular phenotypes include alternative polyadenylation,
isoform expression, and stability.

PLCE1-AS1 is a validated antisense noncoding RNA with detectable expression in
placenta [10], but there is little direct literature tying it to placental
growth regulation. That combination—broad developmental-trait recurrence,
cross-ancestry appearance at different trait comparisons, and sparse prior
placental mechanism—makes it one of the atlas's stronger **novel
developmental-regulation candidates**.

The recurrence should not be read as proof of one causal pathway from birth
weight to later BMI. These are separate GWAS traits. The result instead
suggests that genetic variation at the PLCE1-AS1 locus intersects several
stages of growth and body-size biology.

### 6. Dense single-trait regulatory hotspots

#### 6.1 MRPL43 and JECS C1Y_2

MRPL43 contributes **16 EAS calls**, all to one JECS childhood-BMI code
(C1Y_2), but distributed across four molecular modalities: alternative
polyadenylation, intron retention, isoform expression, and isoform usage. The
maximum PP4 is ~0.954.

The concentration of many RNA features on one GWAS signal suggests a broad
post-transcriptional regulatory haplotype rather than 16 independent
mechanisms. MRPL43 encodes a mitochondrial ribosomal protein; the current
literature search did not identify direct evidence that MRPL43 itself is a
placental childhood-growth effector. This is therefore a strong genetic signal
with comparatively weak prior biological annotation.

#### 6.2 SH3YL1 and JECS C6m

SH3YL1 contributes **seven EAS calls to JECS C6m BMI** across five modalities:
alternative polyadenylation, expression, stability, isoform expression, and
isoform usage. PP4 reaches ~0.912.

As with MRPL43, the scientific value lies in cross-modality consistency at one
locus. The result nominates SH3YL1-associated RNA regulation during a period of
rapid infant growth, but direct placental mechanistic evidence remains limited.
It is a candidate for follow-up rather than a known positive control.

The colocBoost result materially strengthens the *locus-level* case for SH3YL1. Rather than connecting the locus to only one JECS timepoint at a time, the same EAS CoS includes four childhood-BMI GWAS outcomes and multiple placental molecular readouts. This does not prove one biological effect persists unchanged across childhood, but it supports a parsimonious shared developmental-growth signal.

### 7. RNA-editing colocalizations with childhood body size

RNA editing is a comparatively small component of the total call set, but it
contains several highly interpretable recurrent loci for the age-10 body-size
GWAS.

#### 7.1 TFPI

Four EUR TFPI editing sites colocalize with comparative body size at age 10;
the strongest PP4 is ~0.963. TFPI is biologically notable because it is highly
expressed in placenta and sits in the coagulation pathway. Circulating TFPI has
also been reported to rise before clinical preeclampsia and in established
disease [11].

The colocalization trait here is **childhood body size, not preeclampsia**.
Thus the pregnancy literature supports placental plausibility but does not
demonstrate the mechanism connecting TFPI editing to later body size. The
novel observation is that a childhood anthropometric GWAS signal shares
genetic architecture with specific placental TFPI RNA-editing phenotypes.

#### 7.2 CTSB

Four EUR CTSB editing sites also colocalize with the age-10 body-size GWAS,
with PP4 up to ~0.994. Cathepsin B has stronger trophoblast context than many
editing hits: experimental trophoblast work links CTSB expression to
beta-catenin-dependent migration/invasion under a preeclampsia-relevant
hypoxia framework [12].

Again, the result should not be overextended to disease causality. It provides
a candidate bridge between a placental invasion/proteolysis gene and later
body-size genetics through a **post-transcriptional editing phenotype**.

### 8. The IGF2-overlapping composite stability phenotype

One of the most biologically provocative EAS signals is a composite stability
phenotype labeled
`ENSG00000240801.1_ENSG00000284779.3`. It colocalizes with JECS C6m, Dr0m,
and Dr1m BMI; the strongest comparison has PP4 ~0.972 and uses the same xQTL
and GWAS lead variant (`11:2146751:G:A`).

This phenotype must be described carefully. The first Ensembl ID corresponds
to an IGF2-overlapping noncoding locus, but the phenotype ID is composite; it
should **not be renamed as an IGF2 QTL** without checking its construction.
Nevertheless, its genomic context is compelling because placental IGF2 is a
canonical regulator of placental and fetal growth. Experimental deletion of a
placenta-specific Igf2 transcript reduces placental growth and subsequently
causes fetal growth restriction [13].

The colocalization therefore merits targeted annotation of the underlying
stability phenotype and its transcript boundaries. It is a good example of
where the long-read/derived phenotype definition matters as much as the
gene-level label.

### 9. SLC6A2 and neonatal JECS BMI

SLC6A2 has four EAS alternative-polyadenylation calls with the JECS Dr0m BMI
trait, all with PP4 >0.96 and a common GWAS lead. This is a highly resolved
3'-end-regulation signal, but the targeted literature search did not identify
strong evidence for SLC6A2 as a known placental growth gene.

Accordingly, this locus belongs in the **novel-mechanism** category. The
specificity for alternative polyadenylation is useful: rather than simply
nominating SLC6A2 expression, the atlas points to 3'-end choice as the
molecular phenotype most consistent with the GWAS signal.

### 10. Gestational duration: one borderline signal

Only one called comparison involves the 2023 maternal gestational-duration
GWAS: EUR **NCKIPSD expression**, with PP4 ~0.7003. The EGG resource is a large
maternal meta-analysis of gestational duration [5], but the atlas posterior is
essentially at the inclusive call threshold.

This signal should therefore be reported as **suggestive/provisional**, not as
a core colocalization result. A higher PP4 threshold would remove it.

### 10.1 New colocBoost-prioritized loci

#### FES and maternal birth weight

The **FES** mixed event is one of the strongest external-validation additions. In EUR, FES expression and isoform expression colocalize jointly with the EGG maternal-birth-weight GWAS in a **one-variant CoS**. This is important because an independent placental multi-omics study previously used multi-trait colocalization to prioritize **FES** as a likely functional effector for birth weight [6]. The current result therefore recovers FES through a different placental molecular dataset and a different multi-trait colocalization framework.

Notably, FES is not significant in the current gene-level TWAS for maternal birth weight. That discrepancy is scientifically useful: it demonstrates that prediction-based association and shared-signal inference are complementary rather than interchangeable.

#### FAP and decidual remodeling

An EAS **FAP** expression/isoform event colocalizes with a JECS early-childhood BMI outcome. FAP has direct pregnancy context: fibroblast activation protein is expressed on decidual stromal fibroblasts during early pregnancy and has been implicated in physiological tissue remodeling [17]. The growth connection is not established in prior literature, so this remains a novel mechanistic hypothesis, but it is more interpretable than an unannotated locus because the gene maps to a placenta-relevant stromal process.

#### C3orf38 and age-10 body size

The EUR **C3orf38** event colocalizes with comparative body size at age 10 and is also supported by pairwise SuSiE-coloc. The same gene–trait pair is strongly significant in the current TWAS, creating three-way statistical convergence. Direct placental/growth biology for C3orf38 remains sparse; accordingly, it should be highlighted as a high-priority *novel* candidate rather than assigned a mechanism prematurely.

#### Remaining colocBoost-only configurations

**GCLM, MRPL18/TCP1, BCKDHB, RBIS, ARF4, and NPEPL1** occur in mixed events that are absent from the pairwise call set for at least some represented gene–trait pairs. These are best treated as a secondary discovery tier pending locus-level inspection. Their value lies in joint configuration support across outcomes, while the absence of pairwise PP4 support argues for additional sensitivity analysis before mechanistic claims.

### 11. Cross-cutting interpretation

#### 11.1 Placental RNA processing, not only total expression, links to growth GWAS

The largest call class is alternative polyadenylation, and many flagship loci
are supported simultaneously by isoform expression, isoform usage, editing,
intron retention, or stability. The key contribution of this atlas is
therefore not simply another placental eQTL–GWAS colocalization scan. It
resolves *which layer of RNA regulation* shares genetic architecture with a
developmental trait.

The chr17p13 locus is the clearest example: prior placental multi-omics already
implicated CTDNEP1 and nearby regulatory signals in birth weight [6], whereas
the present analysis reveals a much wider spectrum of placental
RNA-processing phenotypes across CTDNEP1, ELP5 and CLDN7.

#### 11.2 Birth weight and childhood size share recurrent loci, but not necessarily one pathway

ELP5/chr17p13 and PLCE1-AS1 recur across birth-weight and childhood-BMI/body-size
resources. This supports developmental continuity at the level of genetic
loci, but it does **not** show that placental regulation mediates a causal path
from birth weight to later adiposity. Formal mediation or multivariable causal
analyses would be needed for that conclusion.

#### 11.3 Exact cross-ancestry replication is absent

Despite 194 calls, no exact gene–trait pair is called in both EAS and EUR.
That absence may reflect real ancestry differences, but several technical and
statistical explanations are equally plausible: different LD structure,
allele frequency, xQTL power, molecular-feature detection, and GWAS ancestry.

The strongest way to use the multi-ancestry design at this stage is therefore
to ask whether *loci and pathways* recur across ancestries, while keeping
exact PP4 calls ancestry-specific.

### 12. Limitations

- **PP4 is evidence for shared genetic causality, not molecular mediation.**
  A colocalized xQTL can be a correlated readout of the causal regulatory
  process rather than the mediator itself.
- **Gene attribution is phenotype-based.** Neighboring genes at a compact
  locus can all inherit evidence from the same haplotype.
- **PP4 ≥0.70 is permissive.** The 0.90–0.95 subset is a more conservative
  shortlist.
- **Repeated molecular phenotypes are correlated.** Row counts should not be
  interpreted as counts of independent causal variants.
- **No exact cross-ancestry replication is present.**
- **The LD reference can influence SuSiE decomposition.** The project uses
  ancestry-matched 1000 Genomes LD rather than the original GWAS individual
  genotypes.
- **The aggregate does not expose all harmonization/QC fields.** Variant
  alignment and per-locus LD diagnostics must be audited upstream.
- **JECS labels are source codes.** A final manuscript should replace codes
  such as `Dr0m`, `C1Y_2`, and `C3HY` with definitions from the original
  phenotype dictionary.
- **Literature absence is search-dependent.** “Novel” here means no compelling
  direct placental/growth mechanism was identified in the targeted searches,
  not proof that none exists.

### 13. Synthesis

The addition of colocBoost changes the interpretation from a purely pairwise catalog to a **pairwise-plus-multi-trait atlas**. SuSiE-coloc remains the clearest per-feature/per-GWAS posterior framework, while colocBoost reveals configurations that span several molecular modalities and several GWAS timepoints. The strongest candidates are those that converge across methods; the most interesting novel candidates are colocBoost-only events with coherent placental biology or TWAS support.


Three findings stand out.

First, **chr17p13 is a compelling positive-validation locus**. Birth-weight
genetics previously implicated the SLC2A4/CLDN7 region [7], and independent
placental multi-omics prioritized CTDNEP1 at the same locus [6]. The current
atlas adds near-certain ELP5/CLDN7 colocalizations and CTDNEP1 RNA-processing
signals, strongly supporting a placenta-relevant regulatory hotspot.

Second, **multi-modality convergence identifies candidates that conventional
eQTL analysis would understate**. WSB1 is the cleanest example, with one fetal
birth-weight signal colocalizing across five RNA-processing layers. MRPL43 and
SH3YL1 show similarly broad molecular effects for individual JECS BMI
timepoints.

Third, the catalog suggests **genetic continuity between placental regulation
and postnatal anthropometry**, particularly at PLCE1-AS1, TFPI, CTSB, and
several JECS-linked loci. These observations are hypotheses about shared
genetic architecture—not proof of developmental mediation—but they create a
focused set of variants and molecular phenotypes for the next stage of causal
and functional work.

### 14. Companion files

- `report_colocalization.Rmd` — descriptive local report.
- `coloc_hits_evidence_table.csv` — literature/source table for `[N]` citations.
- `coloc_best.tsv.gz` — one maximum-PP4 signal pair per comparison cell.
- `coloc_results.tsv.gz` — all tested SuSiE credible-set pairs.
- `coloc_diagnostics_summary.tsv` — run/failure-status aggregation.


### Revised references for colocBoost additions

15. Cao X, Sun H, Feng R, et al. **Integrative multi-omics QTL colocalization maps regulatory architecture in aging human brain.** medRxiv 2025. doi:10.1101/2025.04.17.25326042. https://www.medrxiv.org/content/10.1101/2025.04.17.25326042v3
16. ColocBoost documentation. **Interpretation of ColocBoost output / colocalization confidence sets.** https://statfungen.github.io/colocboost/
17. Scholz M, et al. **Expression of the Peptidase “Fibroblast Activation Protein” on Decidual Stromal Cells Facilitating Tissue Remodeling.** 2020. PMID:33171480. https://pubmed.ncbi.nlm.nih.gov/33171480/

