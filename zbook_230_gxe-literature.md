# Literature context: G×ppBMI findings {#gxe-literature}

<p class="source-note">Book chapter rebuilt locally from <code>reports/literature_gxe_ppBMI.md</code>.</p>


**Placenta xQTL Atlas — Module 07 G×E (ppBMI)**  
**Review date:** 2026-10-12

### Scope

This review places the current placenta genotype × maternal pre-pregnancy BMI
(G×ppBMI) results in the context of prior work on maternal adiposity, placental
molecular biology, and BMI-dependent genetic regulation. The project results
are treated as the primary evidence; external literature is used to assess
biological plausibility and novelty.

A critical distinction is maintained throughout: most prior placental studies
test **maternal BMI or obesity main effects** on placental molecular traits,
whereas the present analysis asks whether **genetic effects on placental
molecular phenotypes vary as a function of ppBMI**. A BMI-associated placental
gene is therefore not automatically a G×BMI gene, and vice versa.

### Summary of the current ppBMI G×E results

The supplied Module-07 aggregates contain 444,582 ancestry-specific Tier-1
feature tests for the exact `ppBMI` exposure label. At q < 0.05, there are
2,067 significant ancestry-specific rows representing 2,064 unique molecular
features. Splicing accounts for 1,782 of these rows (86.2%), with a particularly
large EUR contribution.

Cross-ancestry synthesis is much more restrictive. Only two molecular features
are significant in both EAS and EUR, and neither uses the same ancestry-specific
lead variant. The features map to **TAF6** and **MRAS**. These should be
described as recurrent cross-ancestry molecular phenotypes rather than
variant-level replications. Across the broader Tier-1 meta table, 1,285
features pass the feature-level ACAT q-value threshold.

Tier 2 contains 35,090 genetically prioritized phenotypes. Although 33,949
(96.7%) have `q_meta < 0.05`, most are supported by only one ancestry and the
set was pre-enriched by external genetic prioritization. Requiring two
ancestries, concordant interaction direction, nominal P < 0.05 in both
ancestries, and P_heterogeneity >= 0.05 reduces the set to 664 molecular
features. That stricter subset is more useful for biological prioritization
than the raw Tier-2 q-value count.

### 1. Maternal adiposity measurably alters placental molecular state

There is strong prior evidence that maternal BMI and obesity affect the human
placental transcriptome, epigenome, and metabolic environment.

A 2024 genome-wide RNA-seq study of 450 placentas from the Gen3G cohort found
that higher maternal BMI in early pregnancy was associated with lower placental
**EPYC** expression and higher **IGFBP6**, **CHRDL1**, and **CXCL13**
expression after multiple-testing correction [1]. This establishes that
maternal BMI can be reflected in the term placental transcriptome even after
adjustment for major pregnancy covariates.

Earlier placental transcriptome studies reached broadly consistent pathway-level
conclusions. In a small case-control study, maternal pre-pregnancy obesity was
associated with 72 differentially regulated placental genes and with
inflammation/immune, lipid-metabolism, cancer, and angiogenesis pathways [2].
A separate RNA-seq study of placentas from lean versus obese pregnancies found
288 differentially expressed genes and evidence for altered lipid metabolism,
angiogenesis, hormone activity, cytokine activity, oxidative stress, and
lipotoxicity [3]. These studies are not G×E analyses, but they support ppBMI as
a biologically meaningful placental environment.

Placental epigenetic studies provide a parallel line of evidence. In the NICHD
Fetal Growth Studies, pre-pregnancy BMI was associated with placental DNA
methylation at loci including **EGFL7** and **VEZT**; both ppBMI and methylation
at the EGFL7-associated CpG were negatively correlated with placental EGFL7
expression [4]. Thus maternal adiposity is capable of engaging multiple
molecular layers of the placenta.

#### Implication for the present atlas

The present G×ppBMI results extend this literature from *mean molecular
differences by BMI* to *genotype-dependent molecular responses across BMI*.
They are therefore conceptually complementary rather than redundant.

### 2. BMI-dependent regulatory effects are biologically established

Human G×BMI regulatory mapping has precedent outside placenta.

Glastonbury and colleagues analyzed RNA-seq from 856 female twins across
multiple tissues and identified cis G×BMI effects on gene expression, with
significant interactions concentrated in adipose tissue and enriched in
adipocyte enhancers and metabolic/inflammatory programs [5]. This study is
important for two reasons: it shows that BMI can act as a continuous regulatory
environment, and it demonstrates that G×BMI effects can be strongly
tissue-specific.

More recently, Signer and colleagues used an interaction framework to map
BMI-dynamic eQTLs across tissues and showed that BMI-by-genotype models can
identify context-dependent regulatory relationships influenced by cell type,
putative enhancer activity, and inflammatory cytokines [6]. Their results
further support the general premise that a SNP's regulatory effect is not
necessarily fixed across an individual's physiological environment.

#### Relevance to placenta

Neither study establishes the specific placental interactions observed here.
Instead, they validate the **regulatory model class**: BMI is a plausible
continuous context that can modify genetic control of gene expression. The
placenta is especially appropriate for this question because it is directly
exposed to the maternal metabolic milieu while retaining a fetal genotype.

### 3. The strong splicing signal has emerging biological precedent

The most striking feature of the ppBMI Tier-1 results is the predominance of
splicing interactions. This pattern should be treated as a major hypothesis,
not as already-established biology, because direct human placental
genotype×ppBMI splicing studies remain sparse.

Two recent studies nevertheless make the observation biologically plausible.

First, a 2026 rat time-series placental RNA-seq study reported that maternal
high-fat diet disrupted placental spliceosome pathways during mid-to-late
gestation, including altered expression of multiple spliceosome genes [7].
Importantly, the authors observed spliceosome disruption even in
obesity-resistant animals, suggesting that metabolic/dietary context and
adiposity phenotype need not be interchangeable.

Second, a 2026 human study found widespread differential pre-mRNA splicing in
placentas from pregnancies affected by gestational diabetes mellitus and
implicated SRSF10 as a regulator of a subset of the missplicing [8].
Gestational diabetes is not ppBMI, and this study did not test genotype
interactions, but it independently demonstrates that placental RNA processing
is responsive to a metabolic pregnancy environment.

#### Interpretation

The atlas result therefore fits an emerging model in which placental splicing
is metabolically responsive. The novel contribution is the suggestion that
**the magnitude of cis-genetic control over splicing itself depends on
maternal ppBMI**. Replication in an independent placenta cohort would be
particularly valuable because the number and correlation structure of splice
events can magnify modality-level counts.

### 4. Tier-1 candidate loci

#### CFLAR

**CFLAR** is among the strongest/recurrent Tier-1 signals, with significant
ppBMI interactions across molecular modalities and ancestries. The strongest
current signal is an EAS alternative-polyadenylation phenotype.

CFLAR encodes FLICE-inhibitory protein (FLIP), an inhibitor of caspase-8.
Human placental work has demonstrated long and short FLIP isoforms in early
and late gestation placenta and in trophoblast cells, with expression patterns
varying by cell lineage and differentiation state [9]. The study interpreted
FLIP as an important component of trophoblast survival.

This provides direct placental biological support for CFLAR, although not for
a maternal-BMI interaction specifically. A ppBMI-dependent regulatory effect
at CFLAR is mechanistically plausible in a placenta exposed to obesity-related
inflammatory and oxidative-stress pathways, but functional validation is
needed.

**Evidence assessment:** strong placental relevance; G×ppBMI mechanism novel.

#### PRLR

**PRLR** has multiple significant Tier-1 splicing features, particularly in
EUR. The prolactin receptor is directly expressed in human fetal membranes,
decidua, placenta, and trophoblast compartments [10]. The broader placental
lactogenic axis is also sensitive to maternal obesity: obesity has been
associated with disrupted placental lactogen chromosomal architecture and
more than a 40% decrease in term human placental lactogen RNA in one study
[11].

These observations make PRLR one of the more compelling biologically
contextualized G×ppBMI candidates. They do not establish that PRLR splicing is
caused by obesity or that the current interacting variant is functional.

**Evidence assessment:** direct placental receptor evidence plus obesity-
sensitive lactogenic context; specific G×ppBMI splicing effect novel.

#### SKP2

**SKP2** shows significant G×ppBMI evidence across more than one molecular
phenotype in the current results. A 2026 placental study showed that SKP2 is
highly expressed in placenta and co-expressed in syncytiotrophoblast precursor,
syncytiotrophoblast, and cytotrophoblast populations; it also identified a
RhoBTB1/CUL3-dependent mechanism regulating SKP2 stability [12].

This is unusually direct recent placental evidence for a project hit. The
connection to maternal BMI remains untested experimentally.

**Evidence assessment:** strong contemporary placental evidence; ppBMI-
dependent genetic regulation novel.

#### SGMS1

**SGMS1**, encoding sphingomyelin synthase 1, has recurrent EUR splicing
signals in the current Tier-1 results. Maternal obesity is well established to
alter placental lipid metabolism and create a lipotoxic/oxidative-stress
environment [3]. This gives SGMS1 a coherent pathway-level rationale, but the
literature reviewed here does not provide direct evidence that SGMS1 itself is
a placental maternal-BMI target.

**Evidence assessment:** strong pathway plausibility; direct locus-level
placental G×ppBMI evidence currently limited.

#### SLC25A36

**SLC25A36**, a mitochondrial carrier, is another recurrent Tier-1 signal.
Maternal obesity has been linked to altered placental lipid handling,
oxidative stress, and mitochondrial biology [3]. As with SGMS1, this supports
the biological neighborhood rather than the specific gene–variant interaction.

**Evidence assessment:** mitochondrial/metabolic plausibility; direct
placental G×ppBMI evidence limited.

#### SYT8

**SYT8** is one of the most recurrent Tier-1 genes by number of significant
features and appears across ancestry-specific results. Outside placenta, SYT8
has been linked to glucose-regulated insulin-secretory biology, making it a
plausible metabolic-context candidate. However, the literature reviewed for
this report did not identify strong direct evidence connecting SYT8 to
maternal obesity in human placenta.

**Evidence assessment:** interesting metabolic candidate; placental G×ppBMI
role appears novel.

#### PBRM1

**PBRM1** appears across significant molecular features and ancestries in the
project. PBRM1 is a chromatin-remodeling component with reproductive biology
evidence, but direct human placental maternal-BMI evidence is currently thin.
Its value here is therefore primarily statistical/internal until validated.

**Evidence assessment:** mechanistically plausible chromatin regulator;
placental ppBMI specificity novel.

#### TAF6 and MRAS: cross-ancestry molecular recurrence

The only two Tier-1 molecular features significant in both EAS and EUR map to
**TAF6** and **MRAS**. This is important evidence of cross-ancestry recurrence
at the phenotype level, but it is **not same-variant replication**: both
features have different ancestry-specific lead variants. MRAS also has
different interaction directions at the two distinct leads.

These loci should therefore be prioritized for locus-level conditional or
joint fine-mapping rather than described as replicated allele-specific G×ppBMI
effects.

### 5. Tier-2 candidates and how to use them

Tier 2 is enriched by design because it begins with loci supported by
colocalization and/or TWAS-based criteria. Its very high rate of q_meta < 0.05
therefore cannot be compared directly with Tier-1 discovery yield.

A stricter two-ancestry sensitivity subset retains 664 features. Recurrent
genes in this set include **CCDC107, ACTG1, CUTA, CNGA1, DAXX, BTN3A2,
CLEC1A**, and **RPS18**, among others. These are useful targets for integrated
prioritization, especially when the same gene also has independent
fine-mapping, colocalization/colocBoost, TWAS, or cell-type evidence.

The literature review should not assign equal biological weight to every
Tier-2 q-value. In particular:

- one-ancestry Tier-2 rows are discovery/follow-up evidence, not replication;
- the selected variant is the best meta-analysis variant within a
  pre-prioritized phenotype;
- molecular phenotypes from the same gene are correlated; and
- for one-ancestry rows, the final normal-theory meta z-test can be more
  extreme than the finite-df ancestry nominal p-value.

The two-ancestry, directionally concordant subset is therefore the preferable
Tier-2 layer for headline biological interpretation.

### 6. Ancestry heterogeneity is a result, not merely a nuisance

The project shows major ancestry differences in the number and distribution of
G×ppBMI discoveries. This can arise from several nonexclusive sources:
sample-size/power differences, allele-frequency differences, genotype support
within the exposure distribution, LD structure, distinct ancestry-specific
lead variants, and true biological heterogeneity.

The current data do **not** justify concluding that ppBMI-dependent regulation
is intrinsically more prevalent in one ancestry. The scarcity of same-lead
cross-ancestry features makes that especially important. Ancestry-specific
results should remain primary, with cross-ancestry summaries serving as a
synthesis layer.

### 7. Integrated biological interpretation

Taken together, the current results and prior literature support four working
hypotheses.

1. **Maternal adiposity is a genuine placental molecular environment.** Human
   transcriptomic and epigenetic studies consistently demonstrate placental
   responses to maternal BMI and obesity [1-4].
2. **Genetic regulation can be BMI-dependent.** Human G×BMI eQTL studies show
   that adiposity modifies regulatory variant effects in a strongly
   tissue/context-specific manner [5,6].
3. **RNA processing may be an important placental metabolic-response layer.**
   The splicing predominance in this atlas is novel, but recent placental
   metabolic-stress studies independently implicate spliceosome and
   pre-mRNA-splicing regulation [7,8].
4. **The strongest candidates are those with orthogonal evidence.** Genes such
   as CFLAR, PRLR, and SKP2 combine strong internal interaction evidence with
   direct placental biology, whereas several other high-ranking loci remain
   attractive precisely because their placental G×ppBMI roles appear novel.

### 8. Recommended next analyses

For atlas-level prioritization, the most informative next step is to intersect
the G×ppBMI results with independent evidence rather than rank on G×E q-value
alone. In particular:

- annotate Tier-1 and robust Tier-2 features with Module-05 credible sets;
- intersect G×ppBMI genes with SuSiE-coloc and colocBoost loci;
- overlay TWAS gene-level evidence;
- use the cell-type interaction annotations to ask whether G×ppBMI loci are
  concentrated in specific placental cellular contexts;
- distinguish exact-variant overlap from gene-only overlap;
- perform locus plots for TAF6 and MRAS to understand the discordant
  ancestry-specific leads; and
- where individual-level data permit, pursue joint main+interaction
  fine-mapping rather than treating the top marginal interaction variant as
  the causal variant.

### Conclusions

The placenta ppBMI G×E atlas is consistent with a broader literature in which
maternal adiposity reshapes placental molecular state and BMI modifies genetic
regulatory effects in human tissues. Its most distinctive observation is the
large burden of ppBMI-dependent **splicing** regulation. The literature
provides emerging metabolic-splicing precedent but not direct replication,
making this a potentially important novel contribution.

The strongest locus-level biological candidates should be prioritized on
multiple lines of evidence rather than interaction significance alone.
CFLAR, PRLR, and SKP2 have particularly clear independent placental support;
TAF6 and MRAS merit focused cross-ancestry locus analysis; and recurrent
splicing candidates such as SGMS1, SLC25A36, SYT8, ASPRV1, and others provide
a broader discovery set for functional follow-up.

---

### References

1. Sordillo JE, et al. **Higher Maternal Body Mass Index Is Associated With
   Lower Placental Expression of EPYC: A Genome-Wide Transcriptomic Study.**
   *J Clin Endocrinol Metab.* 2024;109:e1159-e1166.
   doi:10.1210/clinem/dgad619.  
   https://pubmed.ncbi.nlm.nih.gov/37864851/

2. Altmäe S, et al. **Maternal Pre-Pregnancy Obesity Is Associated with
   Altered Placental Transcriptome.** 2017. PMID:28125591.  
   https://pubmed.ncbi.nlm.nih.gov/28125591/

3. Saben J, et al. **Maternal Obesity is Associated with a Lipotoxic Placental
   Environment.** 2014. PMID:24484739.  
   https://pubmed.ncbi.nlm.nih.gov/24484739/

4. Shrestha D, et al. **Placental DNA methylation changes associated with
   maternal pre-pregnancy BMI and gestational weight gain.** 2020.
   PMID:32071425.  
   https://pubmed.ncbi.nlm.nih.gov/32071425/

5. Glastonbury CA, et al. **Adiposity-Dependent Regulatory Effects on
   Multi-tissue Transcriptomes.** *Am J Hum Genet.* 2016;99:567-579.
   doi:10.1016/j.ajhg.2016.07.001.  
   https://pubmed.ncbi.nlm.nih.gov/27588447/

6. Signer R, et al. **BMI-genome interactions regulate global gene expression
   with emphasis in brain and gut.** *Cell Genomics.* 2026;6:101280.
   doi:10.1016/j.xgen.2026.101280. PMID:42314667.  
   https://pubmed.ncbi.nlm.nih.gov/42314667/

7. Du X, et al. **Maternal high-fat diet disrupts placental spliceosome
   pathways independent of obesity phenotype: time-series transcriptome
   analysis during mid-to-late gestation in rats.** *Biol Reprod.* 2026;
   114:1980-1990. doi:10.1093/biolre/ioag023. PMID:41593795.  
   https://pubmed.ncbi.nlm.nih.gov/41593795/

8. Engal E, et al. **Gestational Diabetes Mellitus Alters Placental Precursor
   mRNA Splicing.** *Diabetes.* 2026;75. doi:10.2337/db25-0333.
   PMID:41252199.  
   https://pubmed.ncbi.nlm.nih.gov/41252199/

9. Ka H, Hunt JS. **FLICE-inhibitory protein: expression in early and late
   gestation human placentas.** *Placenta.* 2006;27:626-634.
   doi:10.1016/j.placenta.2005.08.004. PMID:16174530.  
   https://pubmed.ncbi.nlm.nih.gov/16174530/

10. Maaskant RA, et al. **The human prolactin receptor in the fetal membranes,
    decidua, and placenta.** *J Clin Endocrinol Metab.* 1996;81:396-405.
    doi:10.1210/jcem.81.1.8550784. PMID:8550784.  
    https://pubmed.ncbi.nlm.nih.gov/8550784/

11. **Obesity and regulation of human placental lactogen production in
    pregnancy.** 2020. PMID:32500948.  
    https://pubmed.ncbi.nlm.nih.gov/32500948/

12. Tvina AL, et al. **Evidence that S-phase kinase associated protein 2
    (SKP2) is ubiquitinated and degraded via a Rho-related BTB domain
    containing 1 (RhoBTB1) and Cullin-3 mechanism in placenta.**
    *Physiol Rep.* 2026;14:e71039. doi:10.14814/phy2.71039.
    PMID:42552612.  
    https://pubmed.ncbi.nlm.nih.gov/42552612/

