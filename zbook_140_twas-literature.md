# Literature context: TWAS findings {#twas-literature}


**Project:** Placenta xQTL atlas  
**Companion report:** `report_twas.Rmd` / rendered `report_twas.html`  
**Results basis:** `twas_results.tsv.gz` and `twas_gene_results.tsv.gz`, with optional integration of `coloc_best.tsv.gz` and `colocboost_clusters.tsv.gz`.

### 1. Scope and interpretation

The TWAS layer asks whether the genetically predicted component of a placental molecular phenotype is associated with an external GWAS trait. The implementation uses FUSION-style predictive weights for both gene-level and long-read isoform phenotypes. FUSION explicitly frames TWAS as prediction of the genetic component of a functional phenotype followed by association testing against GWAS summary statistics [1].

The supplied aggregate contains **258,124 model-level tests** and **131,523 gene-level ACAT rows** across 16 growth-related GWAS traits. At the supplied multiple-testing thresholds, **3,041 model rows have `TWAS.Q < 0.05`** and **2,098 gene × trait × weight-set rows have `acat_q < 0.05`**. The gene-level layer uses the Cauchy combination test (ACAT), which provides a computationally efficient way to combine dependent p-values [2]. In this application, ACAT is useful because multiple gene and isoform prediction models can exist for the same primary gene.

These results require the standard TWAS caution: association of predicted placental expression or isoform abundance with a phenotype does **not** by itself establish mediation. A TWAS signal can reflect LD between causal variants, shared genetic regulation, or pleiotropy. Colocalization is therefore used as an orthogonal prioritization layer rather than as a redundant significance test [3].

### 2. Overall architecture of the current TWAS

The model-level result set contains 8,075 primary Ensembl genes. Prediction weights are unevenly distributed across ancestry panels: approximately two thirds of model tests use pooled weights, with a much larger EAS-specific than EUR-specific component among the remaining tests. The imbalance largely reflects the GWAS resources available for each ancestry: the Japanese JECS childhood traits have extensive EAS-weight testing but almost no EUR-weight testing, whereas the age-10 comparative body-size GWAS has broad EUR-weight testing but essentially no EAS-weight coverage.

This design means that raw EAS-versus-EUR discovery counts are **not a replication experiment**. Pooled weights provide the broadest cross-trait view, while ancestry-specific weights are most informative when the corresponding external GWAS is well matched.

The median cross-validated prediction R² across model-level tests is about **0.025**. Isoform models are numerous and contribute substantially to discovery. This is scientifically useful because long-read transcript architecture can localize a locus to a particular transcript model, but significant transcript models from the same gene are correlated and should not be counted as independent genes or mechanisms.

### 3. Recurrent body-weight genes provide biological calibration

#### 3.1 ADCY3

**ADCY3** is the broadest recurrent gene-level TWAS signal in the current aggregate: it is significant across ten GWAS traits and more than one weight set. This is a strong biological calibration point. Independent human genetics established that rare loss-of-function mutations in ADCY3 can cause severe monogenic obesity, while other loss-of-function variants increase risk of obesity and type 2 diabetes [4,5].

The relevant caution is tissue interpretation. ADCY3 has well-established roles in neuronal ciliary signaling and energy balance; a placental TWAS association therefore does not imply that the known obesity mechanism is placental. Instead, the result shows that cis-genetic predictors learned in placenta track anthropometric GWAS signal at a canonical body-weight locus. Tissue-comparison or colocalization analyses are needed before claiming a placenta-specific effector route.

#### 3.2 LEP and LEPR

**LEP** is significant for ten childhood-BMI traits in the supplied gene aggregate, and **LEPR** for nine. The leptin pathway is central to body-weight regulation, but it is also directly relevant to pregnancy: human placenta expresses leptin and leptin receptors, and placental leptin has long been proposed to participate in fetal growth, trophoblast biology, and maternal metabolic adaptation [6,7].

The recurrence of both ligand and receptor across childhood anthropometric phenotypes is therefore biologically coherent. It should nevertheless be described as convergence on a known metabolic pathway rather than proof that placental leptin signaling mediates later childhood BMI.

### 4. Placenta-specific precedent: CLDN7 and the 17p13 growth locus

A key strength of the atlas is that several TWAS signals fall in loci already implicated by placental functional genomics. Peng and colleagues used placental transcriptome data to show that genetically regulated placental expression underlies birth weight and childhood obesity risk; their results included the **SLC2A4/CLDN7 region** among birth-weight-linked placental expression signals [8].

In the present aggregate, **CLDN7** is significant for fetal/direct birth weight, JECS low birth weight, neonatal BMI, and a later infant BMI trait. The corresponding chr17p13 region is also one of the strongest xQTL–GWAS colocalization hotspots in the atlas. This convergence is more persuasive than TWAS alone because it links (i) prior placental transcriptomic evidence, (ii) the current genetically predicted placental phenotype, and (iii) shared-signal colocalization.

The same compact region contains **ELP5** and **CTDNEP1**, so gene-level assignment remains nontrivial. The locus is best presented as a regulatory hotspot whose candidate effector genes can be ranked by orthogonal evidence rather than as several independent causal genes.

### 5. CTDNEP1: convergence with prior placental multi-omics

**CTDNEP1** is significant in the current TWAS for fetal/direct birth weight and a neonatal BMI phenotype. A 2022 placental multi-omics study integrating birth-weight GWAS with placental methylation and expression independently identified **CTDNEP1** among candidate functional birth-weight genes [9]. That study also used multi-trait colocalization to prioritize placental effector genes.

This prior evidence substantially strengthens the current CTDNEP1 TWAS interpretation. The novel contribution of the present atlas is that CTDNEP1 sits within a much richer RNA-regulatory context, including xQTL signals beyond conventional total expression.

### 6. SLC38A1 connects genetic association to placental nutrient transport

**SLC38A1** encodes SNAT1, one of the three System A neutral amino-acid transporters expressed in human placenta. Placental System A/L activity is closely linked to fetal growth: reduced activity is associated with fetal growth restriction, while increased activity has been linked to fetal overgrowth [10,11].

The strongest current SLC38A1 TWAS association is with a neonatal/early-childhood BMI trait and is significant with both EAS and pooled weights. It also has SuSiE-coloc support in the atlas. This makes SLC38A1 a particularly interpretable candidate because the statistical evidence lands on a pathway with a direct physiological route—maternal-to-fetal nutrient transfer—rather than on a gene with only generic developmental annotation.

The inference should remain precise: the data support a genetically linked placental SLC38A1 molecular phenotype at a growth-associated locus; they do not establish that altered SNAT1 transport is the causal mediator in these individuals.

### 7. FAP: a novel growth association with decidual-remodeling context

**FAP** is one of the most recurrent gene-level TWAS signals, significant across nine traits and more than one weight set. It is additionally involved in a mixed EAS xQTL–GWAS colocBoost event for a JECS childhood-BMI phenotype.

FAP is not a canonical childhood-obesity gene, but there is pregnancy-specific biological context: fibroblast activation protein is strongly expressed on decidual stromal fibroblasts in early pregnancy and has been proposed to facilitate physiological tissue remodeling during decidualization [12]. The current TWAS/colocalization convergence therefore nominates a plausible stromal-remodeling route that would not have been obvious from anthropometric GWAS alone.

This is an example of an atlas-driven candidate where the placental biology is plausible but the specific link to childhood growth remains novel and requires functional follow-up.

### 8. PLCE1-AS1 and SH3YL1: recurrent developmental-growth loci strengthened by colocBoost

**PLCE1-AS1** is significant across ten GWAS traits in the current TWAS and is one of the broadest loci in the pairwise colocalization catalog. colocBoost adds an important multi-trait observation: one EAS colocalization confidence set jointly connects PLCE1-AS1 isoform/poly(A) regulation with five JECS childhood-BMI GWAS outcomes. Several of those gene–trait pairs are also TWAS significant.

**SH3YL1** is significant across nine traits. colocBoost identifies a single EAS multi-trait event connecting multiple SH3YL1 molecular modalities with four JECS BMI outcomes. This pattern is more parsimoniously interpreted as a shared developmental-growth regulatory locus than as four independent GWAS/TWAS mechanisms.

Direct placenta-specific literature for PLCE1-AS1 and SH3YL1 remains limited. Their priority therefore comes from **internal convergence across TWAS, multiple xQTL modalities, and multi-trait colocalization**, making them strong candidates for targeted functional annotation.

### 9. C3orf38: strong statistical triangulation with sparse prior biology

For comparative body size at age 10, **C3orf38** is a striking EUR result: the gene-level TWAS is highly significant and the same gene–trait pair is supported by both SuSiE-coloc and colocBoost. The colocBoost event is exceptionally compact, with a high-purity confidence set at the locus.

Unlike ADCY3 or SLC38A1, C3orf38 has little established placental or growth biology. That contrast is useful. It marks C3orf38 as a **novel statistical-triangulation candidate**—high priority for follow-up precisely because the genetics are coherent while prior mechanistic expectations are weak.

### 10. FES illustrates why TWAS and colocalization should not be forced to agree

The placental multi-omics literature previously prioritized **FES** as a likely functional effector at a birth-weight locus [9]. In the current colocBoost output, an EUR FES expression/isoform-expression event colocalizes with the maternal birth-weight GWAS in a one-variant confidence set. Yet FES is **not significant in the current gene-level TWAS** (`acat_q` is well above 0.05 for maternal birth weight).

This discordance is informative rather than contradictory. TWAS power depends on the predictive weight model and the linear combination of cis variants; colocalization asks whether association patterns share a causal signal. A locus may therefore colocalize convincingly without producing a significant genetically predicted-expression association. FES should be retained as a colocalization-supported placental birth-weight candidate, not discarded because it misses the TWAS threshold.

### 11. Isoform-aware TWAS is a major added value of this atlas

Most published placental TWAS or transcriptome-imputation studies have focused on total gene expression. The present pipeline carries long-read-derived isoform prediction models into FUSION and then combines them at the primary-gene level with ACAT. This has two advantages:

1. it can identify a transcript model with substantially stronger association than the corresponding total-expression model; and
2. it provides a bridge back to the xQTL atlas, where isoform expression, isoform usage, alternative TSS/poly(A), intron retention, editing, and stability are separately mapped.

The cost is increased correlation and multiplicity. Isoform hits should therefore be interpreted through gene-level aggregation and, where possible, through colocalization or fine-mapping rather than simply counting every significant transcript as an independent discovery.

### 12. Practical evidence tiers for manuscript prioritization

A useful hierarchy for the results is:

**Tier 1 — TWAS + colocalization + independent placental biology.** Examples include CLDN7/CTDNEP1 at chr17p13 and SLC38A1. These are the strongest candidates for mechanistic discussion.

**Tier 2 — TWAS + colocalization, limited prior placental literature.** PLCE1-AS1, SH3YL1, C3orf38, and FAP fit here. These are high-value novel candidates because multiple genetic analyses converge even though external literature is sparse.

**Tier 3 — recurrent TWAS with established general growth biology but incomplete placental causal evidence.** ADCY3, LEP, and LEPR provide strong biological coherence and can anchor the phenotype interpretation, but tissue-specific mediation remains unresolved.

**Tier 4 — TWAS alone.** These associations remain useful for discovery but need the most caution because shared LD and correlated prediction models can generate TWAS signals without a shared causal molecular/GWAS variant.

### 13. Limitations

First, TWAS is not a causal mediation test. Colocalization, conditional analysis, and experimental perturbation are needed to discriminate causal regulation from LD-mediated association [3].

Second, gene-level ACAT combines p-values but does not estimate a single signed effect. The sign of the strongest component TWAS model should not be reported as the direction of the ACAT association without qualification.

Third, ancestry-specific weight coverage is highly unequal across GWAS traits. The present EAS and EUR discovery counts cannot be used as a clean test of ancestry heterogeneity.

Fourth, prediction R² is modest for many models. Strong GWAS loci can generate statistically compelling TWAS results even when molecular prediction is imperfect, and weak prediction reduces power for true effects.

Fifth, several extreme p-values are stored as numerical zero because of floating-point underflow. Ranking among such values requires the underlying test statistic or higher-precision computation rather than the printed p-value.

Finally, placenta is a heterogeneous organ. Bulk molecular prediction can mix trophoblast, stromal, vascular, immune, and other cell-type programs. The atlas's downstream cell-type annotations should therefore be used when translating a gene-level TWAS into a cellular mechanism.

### 14. Synthesis

The TWAS adds a trait-facing layer to the placental xQTL atlas. It recovers canonical body-weight genetics, identifies placental-growth loci with strong prior support, and extends the atlas into isoform-aware predicted molecular phenotypes. The most defensible biological claims arise where **TWAS, xQTL–GWAS colocalization, and prior placental biology converge**.

Three themes stand out. First, **chr17p13 (CLDN7/ELP5/CTDNEP1)** remains the clearest birth-weight regulatory hotspot. Second, **SLC38A1** links the genetic results to a well-established placental nutrient-transport mechanism. Third, **PLCE1-AS1, SH3YL1, FAP, and C3orf38** exemplify the atlas's discovery potential: they gain priority from cross-method convergence rather than from being expected positive controls.

At the same time, FES shows why the report should resist a simple “TWAS-significant = causal” framework. Multi-trait colocalization can provide compelling locus evidence even when a prediction-based test is null. The combined Module 6 outputs are therefore most powerful when interpreted as complementary views of the same genetic architecture.

### References

1. Gusev Lab. **FUSION / Functional Summary-based Imputation.** GitHub repository and documentation. https://github.com/gusevlab/fusion_twas
2. Liu Y, Chen S, Li Z, et al. **ACAT: A Fast and Powerful p Value Combination Method for Rare-Variant Analysis in Sequencing Studies.** *Am J Hum Genet.* 2019;104:410–421. doi:10.1016/j.ajhg.2019.01.002. https://pubmed.ncbi.nlm.nih.gov/30849328/
3. Wainberg M, Sinnott-Armstrong N, Mancuso N, et al. **Opportunities and challenges for transcriptome-wide association studies.** *Nat Genet.* 2019;51:592–599. doi:10.1038/s41588-019-0385-z.
4. Saeed S, Bonnefond A, Tamanini F, et al. **Loss-of-function mutations in ADCY3 cause monogenic severe obesity.** *Nat Genet.* 2018;50:175–179. doi:10.1038/s41588-017-0023-6.
5. Grarup N, Moltke I, Andersen MK, et al. **Loss-of-function variants in ADCY3 increase risk of obesity and type 2 diabetes.** *Nat Genet.* 2018;50:172–174. doi:10.1038/s41588-017-0022-7.
6. Masuzaki H, Ogawa Y, Sagawa N, et al. **Nonadipose tissue production of leptin: leptin as a novel placenta-derived hormone in humans.** *Nat Med.* 1997;3:1029–1033.
7. Henson MC, Castracane VD. **Leptin in pregnancy: an update.** *Biol Reprod.* 2006;74:218–229. doi:10.1095/biolreprod.105.045120.
8. Peng S, Deyssenroth MA, Di Narzo AF, et al. **Genetic regulation of the placental transcriptome underlies birth weight and risk of childhood obesity.** *PLoS Genet.* 2018;14:e1007799. doi:10.1371/journal.pgen.1007799.
9. Tekola-Ayele F, Zeng X, Chatterjee S, et al. **Placental multi-omics integration identifies candidate functional genes for birthweight.** *Nat Commun.* 2022;13:2384. doi:10.1038/s41467-022-30007-1.
10. Lager S, Powell TL. **Regulation of nutrient transport across the placenta.** *J Pregnancy.* 2012;2012:179827. doi:10.1155/2012/179827.
11. Rosario FJ, Powell TL, Jansson T. **Regulation of Placental Amino Acid Transport in Health and Disease.** 2024 review, PMCID: PMC11162343. https://pmc.ncbi.nlm.nih.gov/articles/PMC11162343/
12. Scholz M, et al. **Expression of the Peptidase “Fibroblast Activation Protein” on Decidual Stromal Cells Facilitating Tissue Remodeling.** 2020. PMID:33171480. https://pubmed.ncbi.nlm.nih.gov/33171480/
13. Cao X, Sun H, Feng R, et al. **Integrative multi-omics QTL colocalization maps regulatory architecture in aging human brain.** medRxiv 2025. doi:10.1101/2025.04.17.25326042. (ColocBoost methodology.)

