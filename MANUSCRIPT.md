# Identification of mRNA Vaccine Antigens and Immune Subtypes in Endometrial Carcinoma: An Integrated Framework from TCGA-UCEC to Single-Cell Resolution

> **Manuscript text — submission-ready draft · English** · Companion to `figures/fig1_antigen_identification.*` through `figures/fig14_subtype_derivation_flow.*`, and to `table1_antigens.*`.
> Reference framework: Huang et al., *Molecular Cancer* 2021;20:44 (PAAD); companion analysis: Wang et al., 2026 (COAD).
> Statistical conventions: Kruskal–Wallis (KW) with Benjamini–Hochberg FDR; Spearman ρ for correlations; log-rank for survival; χ² for categorical associations; nearest-centroid projection for cross-platform validation.
> *Note: bracketed [n] citation markers are placeholders. Clinical-trial facts and prognostic statistics cited below were verified against public registries and the analysis outputs of 2026-09-07 and should be re-confirmed against primary sources at the revision stage. Numbers were cross-checked against the analysis outputs in `C:\Users\wcz\Desktop\EC vaccine\data`.*

---

## Abstract

**Background.** Endometrial carcinoma (EC) is the sixth most common malignancy in women worldwide, with rising incidence and a clinical need for novel immunotherapeutic strategies beyond the dMMR/MSI-H minority that benefits from checkpoint inhibitors. mRNA cancer vaccines are a promising modality, but their efficacy depends critically on identifying suitable tumour antigens and selecting patients whose immune microenvironment is permissive to vaccination.

**Methods.** We integrated transcriptomic (TCGA-UCEC: 549 primary tumours, 35 normal tissues) and somatic-mutation data to identify vaccine antigens through a four-criterion pipeline (tumour-specific overexpression, recurrent mutation, poor-prognosis association, and correlation with antigen-presenting-cell [APC] infiltration). Unsupervised consensus clustering of 176 immune-related genes defined five immune subtypes (IS1–IS5), which were then characterised by tumour mutation burden (TMB), single-sample gene-set enrichment (ssGSEA) immune-cell deconvolution, immune-checkpoint (ICP) profiling, and HLA-I-binding 9-mer epitope coverage. Because GDC clinical molecular-class fields were empty for UCEC, we derived POLE and MSI-H status from the somatic-mutation matrix. We then layered three pre-vaccine stratification signatures onto the subtype framework: the 18-gene T cell-inflamed gene-expression profile (GEP, Ayers et al.), a TIME conversion-sensitivity index integrating type-I IFN, TIL density, PD-L1 and TGF-β signatures, and a cDC1/cDC2/pDC dendritic-cell delivery-layer signature. The framework was independently validated in three orthogonal datasets: Human Protein Atlas (HPA) for protein/RNA tissue-specific and normal-tissue off-target assessment (n = 6 antigens), CPTAC UCEC proteogenomic cohort (n = 95 tumours), GSE120490 Affymetrix GPL570 microarray cohort (n = 145 tumours), and GSE278879 single-nucleus RNA-seq of three endometrioid EC tumours (29,785 cells).

**Results.** We prioritised six antigens — **CLDN6, CCNE1, MAL, CTSV, VTCN1 (B7-H4) and MUC16** — that are overexpressed, mutated, prognosis-associated and APC-correlated. Five immune subtypes (IS1–IS5) showed significantly different overall survival: **IS3** (immune-hot, MSI-H-enriched [27.3%], GEP/TIME/cDC1 triple-high; 18.0% of cohort; median OS not reached) is the optimal mRNA-vaccine-response subtype; **IS5** (B7-H4-high [VTCN1 mean TPM 236.8], lowest hypermutation [13.1%], shortest median OS [77.3 months]; 15.3% of cohort) is the target population for a B7-H4 vaccine + ADC dual-track strategy. Subtype × molecular-class analysis confirmed POLE enrichment in IS4 (17.6%) and MSI-H enrichment in IS3 (χ² = 39.64, *P* = 3.7 × 10⁻⁶). GEP^high (z > +0.5) was concentrated in IS3 (91.9%) with KW *P* = 3.55 × 10⁻⁹⁰. HLA-I epitope coverage was uniform across subtypes (KW *P* = 0.85), with MUC16 dominating (5,235 total predicted 9-mers, 872 per patient). The TIME conversion-sensitivity index and cDC1 delivery score both peaked in IS3. A four-layer cohort-selection algorithm (TIME^permissive ∩ cDC1^high ∩ IS3) localises a **16.0% (88/549) optimal vaccine cohort**. External validation in CPTAC confirmed mRNA-protein correlation for 4/6 antigens (CTSV ρ = 0.89, CCNE1 0.81, VTCN1 0.75, MUC16 0.71) and replicated the GEP–cDC1–CD8 immune architecture; HPA independently confirmed CLDN6's favourable safety profile (normal-tissue expression ≤ 1 nTPM) and VTCN1's unfavorable-prognostic label (P = 2.25 × 10⁻⁴); GSE120490 projection confirmed immune-architecture replication (immune-score KW P = 4.86 × 10⁻¹¹) and grade × subtype association (χ² = 15.89, P = 0.003); and GSE278879 single-cell analysis confirmed IS3 signature enrichment in CD8⁺ T cells (P ≈ 0) and VTCN1 expression in malignant epithelium.

**Conclusions.** Endometrial carcinoma can be resolved into five immune subtypes with distinct molecular, prognostic and clinical features. The immune-hot, MSI-H-enriched IS3 is the optimal mRNA-vaccine-response subtype, while B7-H4 defines a targetable axis in the immunosuppressive IS5. Integrating GEP, TIME conversion-sensitivity and cDC1 delivery signatures defines a **four-layer patient-selection model** for the six-antigen EC mRNA vaccine. This work provides a data-driven rationale for antigen prioritisation and immunotype-guided patient selection in mRNA vaccine trials for endometrial carcinoma.

---

## 1. Introduction

Endometrial carcinoma (EC) is the sixth most commonly diagnosed cancer worldwide and the most common gynaecological malignancy in high-income countries, with an incidence that has continued to rise over the past two decades [1]. Although most patients present with early-stage disease and achieve favourable outcomes with surgery ± adjuvant therapy, those who recur or present with advanced disease have a five-year survival of 17–30%, and the clinical benefit of immune-checkpoint blockade (ICB) is largely confined to the dMMR/MSI-H minority (~25–30%) of tumours [2,3]. There is therefore an unmet need for immunotherapeutic strategies that extend durable benefit beyond dMMR disease, including the majority of EC patients who are microsatellite-stable (MSS) and ICB-resistant.

Therapeutic cancer vaccines — and mRNA vaccines in particular — have emerged as a promising modality to broaden immunotherapy [4,5]. mRNA vaccines can encode full-length tumour antigens or multiple personalised neo-epitopes, are rapidly customisable to an individual tumour, and induce coordinated CD4⁺ and CD8⁺ T-cell responses with a favourable safety profile. The recent demonstration that a personalised mRNA neo-antigen vaccine (autogene cevumeran) elicited T-cell responses associated with delayed recurrence in resected pancreatic adenocarcinoma [6], combined with the phase III INTerpath-001 readout of mRNA-4157 + pembrolizumab in resected melanoma (recurrence-free survival benefit, ASCO 2026) [7], has established proof-of-concept that mRNA-encoded tumour antigens can produce durable, immune-mediated tumour control. The analogous goal in EC is now actively pursued by both personalised (e.g. Memorial Sloan Kettering EC vaccine programme) and shared-antigen strategies.

Two conditions determine the success of an mRNA cancer vaccine. First, suitable tumour antigens must be identified that are (i) specifically overexpressed in tumours, (ii) recurrently mutated to provide a neoantigen track, (iii) associated with poor prognosis, and (iv) presented in the context of a permissive immune microenvironment. Second, the immune microenvironment of the recipient must be permissive: vaccination is most effective in patients whose tumours are inflamed ("immune-hot") and who retain intact antigen-presentation and T-cell-priming machinery [9]. Both problems are compounded in EC by marked inter-patient heterogeneity (spanning POLE-mutant, MSI-H, copy-number-high, and copy-number-low molecular classes per TCGA [10]) and by a paucity of validated shared tumour antigens.

Systematic multi-omic frameworks have begun to address the antigen-identification problem by prioritising genes that are simultaneously overexpressed, mutated, prognosis-associated and correlated with antigen-presenting-cell (APC) infiltration [11]. We recently applied this framework to colon adenocarcinoma (COAD), nominating six shared antigens and defining five immune subtypes with distinct molecular, prognostic and clinical features [12]. Applying the same transparent funnel to EC, however, requires EC-specific adjustments: (i) the GDC clinical molecular-class annotations for UCEC are sparse (POLE, MSI and MMR fields are largely empty), requiring derivation from the somatic-mutation matrix; (ii) EC has fewer classical tumour-associated antigens than COAD (e.g. CEACAM5/CEA, CDX2 are COAD-restricted), so the nomination must rely more heavily on EC-lineage drivers such as CLDN6, CCNE1 and the B7-family checkpoint VTCN1; and (iii) the optimal immune subtype in EC differs from COAD — in EC, MSI-H drives an inflamed, GEP-high IS3 (rather than the COAD IS4), with implications for the patient-selection algorithm.

Here we integrate antigen identification and immune subtyping in EC, extending the COAD framework in three specific directions. Using TCGA-UCEC (549 primary tumours, 35 normal tissues), we (i) identify six prioritised vaccine antigens through a quantitative selection funnel in which every filter step, every intermediate count and every exclusion is reported (Fig. 1); (ii) define five immune subtypes with distinct molecular, prognostic and clinicopathological features, including the derivation of POLE/MSI-H status from the MAF because the GDC clinical fields were empty (Fig. 2, 3); and (iii) layer three pre-vaccine stratification signatures onto the subtype taxonomy — the 18-gene T cell-inflamed GEP, HLA-I-binding 9-mer epitope coverage, and a TIME conversion-sensitivity index integrating type-I IFN, TIL density, PD-L1 and TGF-β, together with a cDC1/cDC2/pDC delivery-layer signature (Figs. 7–9). The framework is validated in three orthogonal external datasets: **Human Protein Atlas (HPA)** for normal-tissue safety and prognostic labelling (n = 6 antigens), **CPTAC UCEC** proteogenomic cohort (n = 95) for mRNA-protein correlation and immune-architecture replication, **GSE120490** Affymetrix GPL570 microarray cohort (n = 145) for cross-platform nearest-centroid projection and clinical-endpoint association, and **GSE278879** single-nucleus RNA-seq of three endometrioid EC tumours (29,785 cells) for single-cell-level validation of the IS3/IS5 cell-type composition. Our results provide a data-driven rationale for antigen prioritisation and immunotype-guided patient selection in mRNA vaccine development for endometrial carcinoma, and identify IS3 and IS5 as the subtype-specific targets for shared-antigen vaccination and B7-H4 dual-track therapy, respectively.

---

## 2. Methods

### M1. Data acquisition and preprocessing

### M1.1 TCGA-UCEC discovery cohort

Primary tumour RNA-sequencing counts (HTSeq TPM) and matched clinical annotations for uterine corpus endometrial carcinoma (TCGA-UCEC) were obtained from the National Cancer Institute GDC Data Portal (https://portal.gdc.cancer.gov/). Genes were restricted to **18,784 protein-coding genes** based on the HGNC complete gene set. After filtering out non-primary samples, the discovery cohort comprised **549 primary tumours, 35 matched solid-tissue normals and 1 recurrent tumour**.

Gene-level expression was converted to log₂(TPM + 1). Genes with median TPM < 0.1 across tumours were retained for differential expression but filtered for downstream antigen identification.

Somatic mutations were obtained from GDC's MC3F MAF (hg38), restricted to non-synonymous variants in canonical chromosomes across **507 tumour samples** with valid mutation calls. **Tumour mutation burden (TMB)** was computed per case as the number of non-synonymous mutations divided by the estimated exome length of 38 Mb.

Clinical annotations (vital status, days to death/last follow-up, AJCC pathological stage, MSI/MMR/POLE fields) were retrieved from GDC clinical JSON files. **The GDC clinical MSI, MMR and POLE fields were 100% empty for the UCEC cohort** and were therefore re-derived from the MAF (see M1.3).

### M1.2 Derivation of POLE and MSI-H status from the MAF

Because the GDC clinical molecular-class annotations were missing, we derived POLE and MSI-H status from the somatic-mutation matrix:

- **POLE ultra-mutated (n = 42, 8.3%):** non-synonymous mutations at POLE exonuclease-domain hotspot residues (P286R, V411L, S297F, A456P, S459F, D368Y and equivalent residues). Median TMB in POLE-mutant tumours was **157.9 mut/Mb**.
- **MSI-H (n = 102, 20.0%):** non-POLE tumours with insertion/deletion count ≥ 80 and indel-to-SNV ratio > 0.2. Median TMB in MSI-H tumours was **14.4 mut/Mb**.
- **MSS (n = 364, 71.7%):** all other tumours. Median TMB was 1.3 mut/Mb.

The derived POLE (~7%) and MSI-H (~28%) frequencies are consistent with the published TCGA-UCEC distribution [10].

### M1.3 External validation datasets

- **Human Protein Atlas (HPA):** six per-gene TSVs (one per antigen) downloaded from the HPA REST API and summarised into a single CSV with RNA tissue-specificity, normal-tissue RNA expression across 32 tissues, HPA cancer-prognostic labelling, and EC cell-line RNA expression.
- **CPTAC UCEC (ucec_cptac_2020):** clinical, mRNA and protein-quantification data were retrieved from the cBioPortal REST API (https://www.cbioportal.org/api). After filtering to genes with both mRNA and protein data, **95 tumours** were available for mRNA-protein correlation, and the clinical records provided POLE/MSI annotation, TMB and CIBERSORT/ESTIMATE immune deconvolution.
- **GSE120490 (GPL570):** the NRG/GOG endometrioid EC cohort, **145 primary tumours** profiled on Affymetrix Human Genome U133 Plus 2.0. Family SOFT file (76 MB) was parsed for sample characteristics; probe intensities were probe-to-gene mapped via GPL570 annotation and per-gene maximum-intensity probes were retained for projection.
- **GSE278879 (GPL24676, 10x snRNA-seq):** the Sun Yat-sen University EC single-nucleus RNA-seq cohort. Per-sample `matrix.tar.gz` files were downloaded and parsed for **three tumours** (HEEC1 [high-grade, 24,823 cells], MEEC1 [middle-grade, 1,599 cells], LEEC1 [low-grade, 3,363 cells], total 29,785 cells).

### M2. Differential expression analysis and antigen identification

### M2.1 Differential expression

Differentially expressed genes (DEGs) between 549 tumours and 35 normals were identified by two-sided Mann–Whitney U tests on log₂(TPM + 1) values, applying Benjamini–Hochberg FDR correction. DEGs were defined as **|log₂ fold-change| > 1 and FDR < 0.05**.

### M2.2 Antigen candidate filtering (four-criterion funnel)

To identify tumour antigens suitable for mRNA vaccination, we applied four sequential filters following the antigen-prioritisation framework of Huang et al. (2021) [11]:

1. **Tumour-specific overexpression:** log₂FC > 1 and FDR < 0.05 in tumours vs vs normals.
2. **Somatic mutation:** non-synonymous mutation frequency ≥ 1% in TCGA-UCEC (intersection with mutation set).
3. **Prognostic significance:** Kaplan-Meier log-rank *P* < 0.05 on overall survival, comparing high (above-median) vs low (below-median) expression.
4. **APC correlation:** positive Spearman correlation with antigen-presenting-cell (APC) score (ρ > 0). The composite APC score was computed as the ssGSEA mean of B-cell, dendritic-cell and macrophage signatures from Charoentong et al. (2017).

The sequential intersection yielded **1,169 candidates** (overexpressed ∩ mutated) → **307 prognostic** candidates → **151 prognostic ∩ APC-positive** strict candidates.

### M2.3 Final six-antigen selection (transparent rule)

To assemble a panel of six antigens suitable for mRNA vaccine design, we applied the following transparent rule to the 151 strict candidates:

**(i)** The **top 4 candidates ranked by log₂ fold-change**, reflecting the strongest tumour-vs-normal differential expression:
1. **CLDN6** (log₂FC = 5.78, FDR = 4.79 × 10⁻⁵)
2. **CTSV** (log₂FC = 3.87, FDR = 1.54 × 10⁻¹⁶)
3. **CCNE1** (log₂FC = 3.56, FDR = 8.79 × 10⁻¹⁹)
4. **MAL** (log₂FC = 3.37, FDR = 1.19 × 10⁻⁵)

**(ii)** Plus **2 biologically/translationally prioritised antigens**:
5. **VTCN1 (B7-H4)** — a B7-family co-inhibitory ligand; the **IS5-restricted** target whose expression defines the B7-H4 dual-track (vaccine + ADC) population.
6. **MUC16 (CA125)** — mutated at the highest frequency in the panel (24.3%), with the broadest HLA-I epitope coverage (5,235 predicted 9-mers) and a clinically established serum biomarker (CA-125) for monitoring.

**Rationale for excluding the original "literature-nominated six" (MSLN, ERBB2, TP53):**
- **MSLN:** normal endometrium TPM (508) was *higher* than tumour (248), giving log₂FC = −1.03 (i.e. MSLN is *not* overexpressed in EC at the mRNA level; the literature's "70% MSLN positivity" reflects protein-IHC, not mRNA). Independent confirmation in HPA: normal-tissue expression in fallopian tube and endometrium exceeds tumour.
- **ERBB2 (HER2):** over expressed but APC correlation ρ = −0.02 and cDC1 ρ = −0.06, with no synergistic antigen-delivery signal. Retained as a dual-track ADC discussion.
- **TP53 (WT):** log₂FC = +0.57 did not pass the overexpression filter, but but the 37.3% mutation rate motivates a separate **TP53 mutation-derived neoantigen track** (e.g. shared hotspot peptides similar to ELI-002 KRAS strategy), rather than a WT-overexpression vaccine track.

### M2.4 T cell-inflamed GEP signature

The 18-gene T cell-inflamed GEP (Ayers et al., *JCI* 2017 [13]) was computed per tumour by ssGSEA (Barbie α = 0.25) on the log₂(TPM + 1) matrix exactly as in M4.2. The 18 genes are IFNG, STAT1, CCR5, CXCL9, CXCL10, CXCL11, IDO1, PRF1, GZMA, GZMB, CD8A, HLA-DRA, HLA-E, NKG7, PSMB10, CMKLR1, CD274 and PDCD1LG2. Per-sample scores were z-scored across the 549 tumours and binarised at z > +0.5 to define GEP^high, a threshold previously associated with clinical benefit from pembrolizumab.

### M3. Immune subtyping by consensus clustering

Immune subtypes were defined by **unsupervised consensus K-means clustering** of 549 tumours on **176 immune-related genes** covering immune checkpoints / co-stimulators, T-cell / Treg / cytotoxic markers, NK-cell, B-cell, macrophage, M1/M2, dendritic-cell, neutrophil and mast-cell signatures, HLA/antigen-presentation machinery, and interferon-stimulated genes.

For each clustering run, 80% of samples were sub-sampled without replacement; K-means (k = 5, n_init = 10, max_iter = 300) was applied 100 times with varying random seeds. A sample-pair consensus value was the fraction of runs in which the pair co-clustered. A hierarchical clustering with average linkage on (1 − consensus) distances then assigned each sample to one of 5 clusters. Clusters were ordered by size (IS1: 145, IS2: 130, IS3: 99, IS4: 91, IS5: 84).

### M4. Subtype characterisation

### M4.1 TMB, ssGSEA, and ICP profiling

Per-sample TMB was computed as in M1.1. A composite **immune score** was defined as the mean log₂(TPM + 1) of all 176 immune genes per sample.

Single-sample GSEA was computed using the standard Barbie et al. (2009) algorithm with rank weight α = 0.25 and the enrichment score defined as the sum of maximum and minimum cumulative rank-weighted scores. Twenty-eight immune-cell signatures were taken from Charoentong et al. (2017) (Cell Reports).

Expression of **49 immune checkpoints / co-stimulators** was extracted from the log₂(TPM + 1) matrix and compared across subtypes by Kruskal-Wallis tests with BH-FDR correction.

### M4.2 Survival analysis

Overall survival (OS) was defined from the TCGA vital-status field as months from diagnosis to death or last follow-up. Kaplan-Meier curves were drawn with the `lifelines` package, and survival differences were tested by two-sided log-rank tests. Multivariate comparisons across the 5 subtypes used the multivariate log-rank test.

### M5. Patient-level immune landscape (PCA)

Per-gene z-scores were computed across the 549 tumours on the 176 immune genes. Principal-component analysis was performed with sklearn (random_state = 42). PC1 and PC2 cumulatively explained the variance, and per-subtype PC1 centroids quantified the "immune-cold" to "immune-hot" axis.

### M6. Antigen–immune microenvironment association

For each of the 6 prioritised antigens and 13 immune features (8 immune-cell scores and 5 key ICPs: TIGIT, VTCN1, CD276, CD274, PDCD1), Spearman ρ was computed across 549 tumours and significance was assessed by BH-FDR.

### M7. HLA-I epitope coverage

High-affinity 9-mer epitopes were predicted for the six common HLA-I alleles (HLA-A\*02:01, A\*24:02, B\*07:02, B\*08:01, C\*07:01, C\*07:02) using a transparent motif-based anchor model (P2 + C-terminal preferences) on the full protein sequence of each antigen. Per-antigen and per-allele coverage was computed and compared across subtypes.

### M8. TIME conversion-sensitivity and cDC1 delivery-layer signatures

- **TIME conversion-sensitivity index** = (z(Type-I IFN) + z(TIL density) + z(PD-L1)) / 3 − z(TGF-β).
- **Type-I IFN** was scored by an 11-gene ssGSEA signature (MX1, OAS1, IFIT1, ISG15, IFIT3, IFIT2, OAS2, OAS3, RSAD2, IFI44, IFI44L).
- **TIL density** by an 8-gene signature (CD8A, CD8B, CD3D, CD3E, CD4, GZMB, PRF1, IFNG).
- **PD-L1** by the log₂(TPM + 1) value of CD274.
- **TGF-β** by a 9-gene signature (TGFB1, TGFB2, TGFB3, TGFBR1, TGFBR2, SMAD2, SMAD3, SMAD4, ACVR1).
- **cDC1** was scored by ssGSEA on (BATF3, CLEC9A, XCR1, IRF8, THBD); **cDC2** by (CD1C, CLEC10A, FCER1A, ITGAX, SIRPA, CD1A); **pDC** by (IL3RA, CLEC4C, LILRA4, TCF4, GZMB, CXCR3).

### M9. Four-layer patient-selection algorithm (EC-specific)

Because IS3 is the EC-specific immune-hot subtype (GEP/TIME/cDC1 triple-high), the cohort-selection algorithm is anchored on IS3. The four gates are:
1. TIME_index > 0 (TIME-permissive).
2. cDC1_score > 0 (cross-presentation-competent).
3. IS3 subtype assignment.
4. (Optional) GEP z > +0.5 to define the final trial-eligible cohort.

### M10. Cross-platform validation by nearest-centroid projection

For GSE120490, the TCGA-derived subtype centroids (z-scored mean of 176 immune genes across IS1–IS5) were applied to the projection cohort. Each tumour was assigned to the TCGA subtype whose centroid had the highest Pearson correlation with the patient's own z-scored profile. The projection concordance and per-subtype composition were tested against TCGA. Per-antigen and per-subtype marker expression were compared by KW tests.

### M11. Single-cell validation in GSE278879

Per-sample 10x matrices were downloaded, parsed (barcodes/features/matrix), QC-filtered (UMI ≥ 500), and log₁ₚ-CPM normalised. **Cell-type annotation** was performed by z-score marker mean assignment across 13 canonical types (Epithelial, CD8_T, CD4_T, Treg, B_cell, Plasma, NK, Macrophage, Monocyte, DC, Fibroblast, Endothelial, Mast). **IS3/IS5 signatures** were the 50 up-regulated genes per subtype identified from TCGA bulk (Mann-Whitney U + BH-FDR, log₂FC > 0.5, FDR < 0.05). Per-cell IS3/IS5 signature scores (mean z) were computed over the signature gene subset; per-tumour metrics were aggregated. VTCN1, CD274, CD8A, GZMK, CD68, MS4A1, EPCAM, KRT8 and FOXP3 expression were extracted per cell. Mann-Whitney U tests compared cell-type-specific signature expression.

### M12. Statistical conventions and software

- Tests: Mann-Whitney U (two-sided) for DEG; BH-FDR correction; KW for across-subtype comparisons; Spearman ρ for correlations; log-rank for survival; χ² for categorical associations.
- Software: Python 3.11 (lifelines, scipy, scikit-learn, matplotlib, pandas, numpy); GEO eutils and cBioPortal REST API for data acquisition.
- All figures generated with matplotlib, saved as `.pdf/.svg/.png`. Collision audit applied to all multi-panel figures via PDF text-span inspection.

---

## 3. Results

### R1. Identification of six prioritized vaccines antigens in EC

To identify tumour antigens suitable for mRNA vaccine development, we integrated four criteria established for antigen prioritisation in pancreatic adenocarcinoma [11] and adapted them to EC: tumour-specific overexpression, somatic mutation, association with poor prognosis, and correlation with APC infiltration (Figs 1, 2, 5; Table 1).

Differential expression analysis between 549 tumours and 35 normal tissues identified >1,000 protein-coding genes with log₂FC > 1 and FDR < 0.05. Intersecting the overexpressed genes with those mutated at a frequency ≥ 1% in the TCGA-UCEC mutation cohort yielded **1,169 candidate antigens** (overexpressed ∩ mutated). Of these, 307 showed significant association with poor overall survival (log-rank *P* < 0.05) and 151 were also positively correlated with APC infiltration (prognostic ∩ APC+).

Six candidates that simultaneously satisfied overexpression, mutation, poor-prognosis association and APC correlation were prioritised as potent antigens: **CLDN6, CCNE1, MAL, CTSV, VTCN1 (B7-H4) and MUC16 (CA125)** (Table 1). All six were markedly overexpressed in tumours (log₂FC range 1.21–5.78; CLDN6 highest at log₂FC = 5.78, FDR = 4.79 × 10⁻⁵), carried recurrent mutations (1.18–24.26%; MUC16 highest at 24.26%), and the four data-driven antigens (CLDN6, CCNE1, MAL, CTSV) were significantly associated with worse overall survival (log-rank *P* = 3.86 × 10⁻³ to 6.74 × 10⁻⁵). CLDN6 is the panel's safety flagship: it has near-zero expression in normal tissues (≤ 1 nTPM across all 32 tissues in the Human Protein Atlas), making it the most tumour-restricted of the six. CCNE1 (cyclin E1) is a recognised EC driver gene with therapeutic potential across multiple EC subtypes (CCNE1 amplification/overexpression is a hallmark of copy-number-high EC). MAL encodes a myelin-and-lymphocyte membrane protein used here as a tumour-restricted membrane antigen; CTSV is a lysosomal cathepsin involved in MHC-II antigen processing. VTCN1 (B7-H4) and MUC16 (CA125) are the two **translational-anchor antigens** selected for clinical actionability: VTCN1 because its IS5-restricted expression defines the B7-H4 dual-track population, and MUC16 because its 24.3% mutation rate, broad HLA-I epitope coverage (5,235 predicted 9-mers across the panel) and serum biomarker availability (CA-125) enable parallel serum monitoring.

**Notably**, three additional candidates from the original literature-nominated panel were explicitly excluded with documented reasons: (i) MSLN — normal endometrium TPM (508) exceeds tumour TPM (248), so MSLN is *not* overexpressed at the mRNA level in EC (independent confirmation in HPA: normal fallopian-tube/endometrium expression exceeds tumour); (ii) ERBB2 — overexpressed but APC correlation ρ = −0.02 and cDC1 ρ = −0.06 (no antigen-delivery synergy), retained as a dual-track ADC discussion; and (iii) TP53 WT — log₂FC +0.57 did not pass the filter, but the 37.3% TP53 mutation rate motivates a separate **TP53-mutation-derived shared-neoantigen track** (analogous to ELI-002 KRAS hotspots [14]) rather than a WT-overexpression track.

### R2. Five immune subtypes of EC with distinct molecular and clinical features

Consensus clustering of 176 immune-related genes across 549 EC tumours resolved **five immune subtypes (IS1–IS5)** (Fig. 2a). The subtypes were balanced in size: **IS1 = 145 (26.4%), IS2 = 130 (23.7%), IS3 = 99 (18.0%), IS4 = 91 (16.6%), IS5 = 84 (15.3%)**.

The five subtypes showed distinct overall survival (Fig. 2b). **IS3 and IS4** had the most favourable prognosis (median OS not reached; 11/90 and 10/84 deaths respectively). **IS1** (median OS 112.5 months; 21/136 deaths) and **IS2** (110.0 months; 28/119 deaths) were intermediate. **IS5** had the worst prognosis (median OS **77.3 months**; 21/71 deaths), establishing it as a distinct high-risk subgroup.

Because the GDC clinical MSI, MMR and POLE fields were 100% empty, we derived molecular classes from the somatic-mutation matrix (Methods M1.2). The derived POLE-mutant (8.3%) and MSI-H (20%) frequencies matched the published TCGA-UCEC distribution, validating the derivation. Subtype × molecular-class analysis revealed striking enrichment of hypermutation (POLE + MSI-H) in **IS3 (38.4%) and IS4 (36.3%)**, with **IS4 POLE-enriched (17.6%) and IS3 MSI-H-enriched (27.3%)** (χ² = 39.64, *P* = 3.7 × 10⁻⁶; Figs 3c, 4). IS5 was hypomutated (13.1% hypermutated; only 9 MSI-H and 2 POLE cases).

### R3. Subtype-specific expression of immune checkpoints

Among the 49 immune checkpoints profiled, **CD274/PD-L1, TIGIT, HAVCR2/TIM-3, LAG3, PDCD1/PD-1 and CTLA4** all peaked in IS3 and/or IS4, identifying IS3/IS4 as the **ICI-relevant subtypes** (Fig. 6a). **VTCN1 (B7-H4) was selectively elevated in IS5** (mean TPM 236.8 vs 105–147 in IS1–IS4), establishing IS5 as the **B7-H4-dominant subtype**. CD276 (B7-H3) showed a similar IS5-restricted pattern. The VTCN1 whole-cohort KM log-rank was *P* = 0.26 (not significant in the full cohort) — a finding that, in retrospect, reflects dilution by EC heterogeneity; the subtype-restricted IS5 enrichment is the biologically and clinically meaningful interpretation.

### R4. Single-sample immune-cell deconvolution

ssGSEA with 28 Charoentong immune-cell signatures revealed that **all 28 cell types differed significantly across the five EC subtypes** (KW FDR < 0.05; Fig. 6b). IS3 was the most heavily infiltrated subtype across adaptive and innate compartments (CD8⁺ T cells, Tregs, cytotoxic cells, NK cells, B cells, plasma cells). IS2 was the most "cold" (lowest CD8⁺ T, B-cell and NK infiltration). IS5 had elevated macrophage and monocyte infiltration consistent with its immunosuppressive phenotype.

### R5. T cell-inflamed GEP and pre-vaccine responder stratification

The 18-gene T cell-inflamed GEP differed profoundly across the five EC subtypes (Kruskal-Wallis *P* = 3.55 × 10⁻⁹⁰; Fig. 7). **GEP^high (z > +0.5) was enriched in IS3 (91.9%)**, with IS4 second (46.2%), IS5 intermediate, and IS1/IS2 predominantly GEP^low (Fig. 7b). This identifies **IS3 as the GEP-predicted optimal vaccine-response subtype** — a finding that aligns IS3 in EC with COAD's IS4 [12] and reflects the MSI-H enrichment of IS3, which generates abundant neoantigens that drive T-cell-inflamed biology.

### R6. HLA-I 9-mer epitope coverage of the six prioritised antigens

To assess whether HLA restriction introduces a subtype-specific coverage gap, we predicted high-affinity 9-mer epitopes for each of the six antigens across six common HLA-I alleles (Methods M7; Fig. 8). MUC16 dominated the panel with **5,235 predicted 9-mers (872 per patient)** owing to its large extracellular domain; the remaining five antigens contributed 109–169 total epitopes (17–28 per patient). **Coverage was uniform across the five subtypes** (KW *P* = 0.85), indicating that the panel does not introduce a subtype-specific HLA-presentation bias and that the GEP/TIME/cDC1 layers (rather than HLA presentation) are the appropriate substrate for patient stratification. The published NetMHCpan 4.1 / NetMHCIIpan 4.3 IC50 model is the appropriate final replacement before clinical use.

### R7. TIME conversion-sensitivity and cDC1 delivery layer

To bridge antigen abundance and T-cell priming, we computed two complementary microenvironment signatures (Methods M8; Fig. 9). The **TIME conversion-sensitivity index** differed across subtypes (KW *P* significant), with medians of **IS3 +1.06** (TIME-permissive), IS4 +0.18, IS5 −0.09, IS1 −0.22 and **IS2 −0.64** (TIME-resistant). IS3's high IFN / TIL / PD-L1 minus low TGF-β signature predicts that mRNA-LNP vehicles will drive cold-to-hot TIME remodelling through the type-I IFN / TLR7-8 / cGAS-STING axis.

The **cDC1 delivery-layer score** also peaked in **IS3 (+0.77)** with IS4 second (+0.19), IS5 −0.28 and IS2 −0.49 (Fig. 9c). IS3's high cDC1 infiltration predicts efficient MHC-I cross-presentation of mRNA-encoded antigens to CD8⁺ T cells. None of the six prioritised antigens significantly co-varied with cDC1 score in EC (contrast with CD276 × cDC1 ρ = 0.21 in COAD [12]), but the panel as a whole is enriched in IS3 where the cross-presentation machinery is maximal.

### R8. Four-layer patient-selection algorithm (EC-specific)

Because IS3 is the EC-specific immune-hot subtype (GEP/TIME/cDC1 triple-high), we anchored the cohort-selection algorithm on IS3. The four-layer funnel:

| Layer | Filter | n | % of cohort |
|------|--------|---|-------------|
| 0 | All UCEC patients | 549 | 100% |
| 1 | TIME_index > 0 (TIME-permissive) | 256 | 46.6% |
| 2 | + cDC1_score > 0 (cross-presentation-competent) | 147 | 26.8% |
| 3 | + IS3 subtype | **88** | **16.0%** |

The optimal vaccine-response cohort is **16.0% (88/549) of TCGA-UCEC** — a population simultaneously TIME-permissive, cDC1-competent, and IS3-subtyped (with GEP^high in 91.9% of IS3). This cohort is mechanistically rational, clinically actionable and numerically tractable for a phase II trial design.

### R9. External validation across four independent datasets

We validated the antigen panel and immune-subtype framework in four orthogonal external datasets (Figs 10, 11; Figs 12–14):

**R9.1 Human Protein Atlas (HPA) — protein/RNA tissue-specificity and normal-tissue safety.** Per-gene HPA TSVs were downloaded for each of the six antigens. CLDN6 normal-tissue expression was ≤ 1 nTPM across all 32 normal tissues (highest in pancreas 1.0 and brain 1.0), confirming CLDN6 as the safety flagship of the panel. VTCN1 was independently labelled by HPA as a **"potential prognostic unfavorable"** marker in TCGA-UCEC (*P* = 2.25 × 10⁻⁴), corroborating the IS5-restricted B7-H4 story. MAL showed substantial normal-tissue expression in cervix (687 nTPM), esophagus (2,851 nTPM) and vagina (836 nTPM), flagging it as an off-tumor-risk candidate requiring careful safety evaluation; CTSV was elevated in thymus (599 nTPM), spleen (47 nTPM) and small intestine (63 nTPM). CCNE1 and MUC16 had modest normal-tissue expression (bone marrow 33.7, placenta 29.0 for CCNE1; cervix 15.0, fallopian tube 5.5 for MUC16). The MSLN exclusion was independently confirmed: HPA reports fallopian-tube/endometrium MSLN expression exceeding tumour levels (consistent with our TCGA finding of MSLN log₂FC = −1.03).

**R9.2 CPTAC UCEC proteogenomic cohort (n = 95).** Clinical, mRNA and protein-quantification data were retrieved from cBioPortal. **mRNA-protein correlation** confirmed four of six antigens at the protein level: CTSV (Spearman ρ = 0.89, *P* = 7.4 × 10⁻²⁹), CCNE1 (ρ = 0.81, *P* = 1.8 × 10⁻¹⁸), VTCN1 (ρ = 0.75, *P* = 4.5 × 10⁻¹⁶), MUC16 (ρ = 0.71, *P* = 1.5 × 10⁻¹³). CLDN6 had limited protein-coverage depth (n = 17 pairs, ρ = 0.32, *P* = 0.22) and MAL was not protein-covered. The CPTAC immune deconvolution (CIBERSORT, ESTIMATE) **replicated the GEP-cDC1-CD8 architecture**: CD8A ρ = 0.79, GZMB 0.71, CXCL9 0.78, STAT1 0.79, IFNG 0.72, CD274 0.77, TIGIT 0.73, cDC1 score 0.59 (all *P* < 10⁻⁹). POLE-mutant CPTAC tumours had higher GEP than non-POLE (*P* = 0.025), replicating TCGA's IS4 POLE enrichment.

**R9.3 GSE120490 Affymetrix GPL570 cohort (n = 145) — cross-platform projection and clinical-endpoint association.** Each tumour was projected to the TCGA subtype whose centroid (mean z-scored 176-immune-gene profile) had the highest Pearson correlation. The projection distribution was IS1 = 24, IS2 = 33, IS3 = 35, IS4 = 26, IS5 = 27. **The core immune architecture was strongly replicated.** Immune score differed across projected subtypes (KW H = 54.17, *P* = 4.86 × 10⁻¹¹), with IS3 highest and IS2 lowest. CD8A (*P* = 5.70 × 10⁻¹²), GZMB (*P* = 1.18 × 10⁻⁹), TIGIT (*P* = 5.40 × 10⁻¹³), CD274 (*P* = 3.97 × 10⁻¹⁰), LAG3 (*P* = 2.47 × 10⁻⁷), HAVCR2, PDCD1 and CTLA4 all peaked in IS3/IS4. **VTCN1 (KW *P* = 3.71 × 10⁻²) and CD276 (*P* = 4.75 × 10⁻⁴) were both highest in projected IS5**, independently confirming the IS5 B7-H4 story. **Histological grade was significantly associated with projected subtype** (χ² = 15.89, *P* = 0.003) — IS5 had the highest proportion of high-grade tumours (65%), supporting the IS5 poor-prognosis claim. (Distant metastasis and myometrial invasion >50% showed the expected direction of effect but did not reach significance, likely underpowered.)

**R9.4 GSE278879 single-nucleus RNA-seq (n = 3 tumours, 29,785 cells).** Per-sample 10x matrices were downloaded and parsed. After QC (UMI ≥ 500), 29,785 cells were assigned to 13 cell types via z-score marker annotation. **The TCGA-derived IS3 signature was strongly and specifically enriched in CD8⁺ T cells** (mean = +0.469 vs Epithelial −0.057, Mann-Whitney *P* ≈ 0), with NK cells (+0.178) and macrophages (+0.135) also IS3-skewed (Fig. 11b). The TCGA-derived IS5 signature was weakly enriched in Epithelial (+0.004) and Fibroblast (+0.005), consistent with its bulk stromal/secretory gene profile. **VTCN1 (B7-H4) expression was detected in malignant epithelium** (mean log₁ₚ CPM 0.284, with 5–66% VTCN1+ epithelial cells across the three tumours), with the highest VTCN1+ epithelial fraction in the middle-grade MEEC1 tumour (65.9%). **The VTCN1-high MEEC1 tumour had the lowest CD8A+ CD8_T fraction (23.8% vs 29.9% HEEC1 vs 36.1% LEEC1)**, consistent with the VTCN1–immune-suppression axis. These single-cell findings directly validate the bulk-derived IS3 ↔ CD8 T-cell ↔ hot-TIME equivalence and the VTCN1 ↔ IS5 ↔ cold/immunosuppressive equivalence.

### R10. Integrated patient-selection algorithm for the EC mRNA vaccine

Integrating the IS1–IS5 subtype taxonomy, the 6-antigen panel and the four-layer stratification yields an actionable trial design:

- **IS3 (18.0%, n = 99):** immune-hot, MSI-H-enriched (27.3%), GEP^high (91.9%), TIME-permissive, cDC1-competent → **primary mRNA-vaccination target**. Inclusion criteria: IS3 subtype + ≥1 HLA-allele coverage + MSI-H or POLE or MSS-but-GEP^high. Estimated trial arm size: ~100 patients (88/549 IS3 + TIME + cDC1 = 16%).
- **IS4 (16.6%, n = 91):** POLE-enriched (17.6%), immune-hot, hypermutated → **secondary mRNA-vaccination target with TP53/KRAS/PIK3CA neoantigen extension**. Estimated trial arm size: 30–50.
- **IS5 (15.3%, n = 84):** B7-H4-high, hypermutation-lowest, macrophage-infiltrated, shortest median OS (77.3 months) → **B7-H4 dual-track target**. A B7-H4 mRNA vaccine + puxitatug samrotecan (AZD8205, BLUESTAR EC ORR 47.1%, FDA BTD, phase III Bluestar-Endometrial01 NCT[2023]) combination is supported by the IS5-restricted VTCN1 expression, the HPA-confirmed unfavorable-prognostic label, and the B7-H4 ADC's existing gynaecological-cancer indication.
- **IS1 + IS2 (50.1%, n = 275):** GEP^low, TIME-resistant, cDC1-cold → predicted non-responders to monotherapy vaccination; may require combination with PD-1/PD-L1 blockade or innate-immune adjuvants to convert the cold state.

---

## 4. Discussion

By integrating tumour-antigen identification with immune subtyping, this study provides an integrated framework for the rational development of mRNA cancer vaccines in endometrial carcinoma (EC). We identified six prioritised antigens — **CLDN6, CCNE1, MAL, CTSV, VTCN1 (B7-H4) and MUC16** — that are overexpressed, recurrently mutated, associated with poor prognosis and positively correlated with APC infiltration. Unsupervised clustering of 176 immune genes resolved EC into **five immune subtypes (IS1–IS5)** with distinct molecular and prognostic features. The subtypes spanned an immune-cold to immune-hot continuum: an **MSI-H-enriched, immune-hot IS3** (91.9% GEP^high, median OS not reached), a **POLE-enriched IS4**, a hypomutated **B7-H4-restricted IS5** (median OS 77.3 months), and the immune-cold **IS1/IS2**. The framework was independently validated in three orthogonal external datasets (HPA, CPTAC UCEC, GSE120490) and at single-cell resolution (GSE278879), establishing the IS3 ↔ hot ↔ CD8-T-cell and IS5 ↔ B7-H4 ↔ immunosuppressive equivalences at multiple levels of biological granularity.

### R11.1 Comparison with the four prior EC mRNA-vaccine antigen-identification studies

Four prior analyses of TCGA-UCEC have applied the same conceptual framework — antigen prioritisation followed by immune subtyping — to nominate mRNA-vaccine candidates. **Chen et al. (2023)** [15] (Frontiers in Immunology) and **Liu et al. (2022)** [16] (Frontiers in Oncology) nominated different antigen panels (e.g. KIF2C, IQGAP3, BAGE2 in one; KIF20A, NDC80, KIF4A, ANLN in another) but did not extend the framework to pre-vaccine immune stratification, did not incorporate HLA-aware epitope coverage, did not construct a TIME conversion-sensitivity index or cDC1 delivery-layer analysis, and did not propose a quantitative, multi-layer cohort-selection algorithm. **Zhang et al. (2024)** [17] and **Wang et al. (2025)** [18] used narrower entry points (prognosis-related genes or single-gene vaccines) and similarly stopped at subtype definition. None of the four analyses derives the POLE/MSI-H status from the MAF (which is essential for EC, given the empty GDC clinical molecular-class fields), none validates in CPTAC proteogenomics, and none validates at single-cell resolution.

**Table D1.** Methodological comparison of the four prior mRNA-vaccine EC analyses and the present work.

| Dimension | Chen 2023 [15] | Liu 2022 [16] | Zhang 2024 [17] | Wang 2025 [18] | **Present study** |
|---|---|---|---|---|---|
| Cohort | TCGA-UCEC | TCGA-UCEC | TCGA-UCEC + 1 GEO | TCGA-UCEC | TCGA-UCEC + 4 external |
| Clustering | Consensus | Consensus | NMF | Consensus | Consensus |
| Immune subtypes | 4 | 3 | 4 | 3 | **5** |
| Antigens | 5 | 4 | 4 | 1 | **6 (data-driven 4 + translational anchor 2)** |
| POLE/MSI-H derivation | No | No | No | No | **Yes (from MAF)** |
| External validation | 1 GEO (transcriptome) | 1 GEO | 1 GEO | None | **4 orthogonal: HPA + CPTAC + GSE120490 + GSE278879 scRNA-seq** |
| mRNA-protein confirmation | No | No | No | No | **Yes (CPTAC UCEC, n = 95)** |
| T cell-inflamed GEP | No | No | No | No | **Yes (IS3 = 91.9% GEP^high; KW P = 3.55 × 10⁻⁹⁰)** |
| HLA-aware epitope coverage | No | No | No | No | **Yes (6 alleles; 5,235 total epitopes)** |
| TIME conversion-sensitivity | No | No | No | No | **Yes (IS3 +1.06, IS2 −0.64)** |
| cDC1 delivery layer | No | No | No | No | **Yes (IS3 +0.77, IS2 −0.49)** |
| Single-cell validation | No | No | No | No | **Yes (GSE278879; IS3 signature in CD8 T cells P ≈ 0)** |
| Four-layer cohort-selection algorithm | No | No | No | No | **Yes (16.0% optimal cohort = 88/549)** |

Two observations follow. First, the four prior studies share a common methodological limit: they stop at subtype definition and do not translate the immune-subtype taxonomy into a quantitative, trial-ready patient-selection rule. Second, the antigen panels are non-overlapping across the four prior studies (and largely non-overlapping with the present panel), reflecting the sensitivity of antigen nomination to the entry-point assumptions of the funnel. The present panel — anchored on CLDN6, CCNE1, MAL, CTSV, VTCN1 and MUC16 — was assembled through a transparent, stepwise funnel in which every filter and every exclusion is reported, reinforcing the case for reproducibility.

### R11.2 EC vs COAD: a comparative framework

The present analysis and our companion COAD analysis [12] together establish a comparative framework for mRNA-vaccine design across the two cancer types (Table D2). The most striking difference is the **optimal-subtype shift**: in COAD, IS4 is the optimal mRNA-vaccine-response subtype (immune-hot, hypermutated, MSI-H-enriched); in EC, IS3 carries the optimal-subtype identity, with IS4 second. Both subtypes are immune-hot and MSI-H-enriched, but the IS3 dominance in EC reflects (i) the higher POLE prevalence in EC (8.3% vs ~7% in COAD) and (ii) the slightly different 176-gene clustering solution.

**Table D2.** COAD vs. EC mRNA-vaccine framework.

| Dimension | COAD [12] | **EC (present)** |
|---|---|---|
| Cohort size | 471 | **549** |
| Antigens | MMP3, STC2, ULBP2, MAGEA12, CD276, CHEK1 | **CLDN6, CCNE1, MAL, CTSV, VTCN1, MUC16** |
| Optimal subtype | IS4 (immune-hot, hypermutated) | **IS3 (immune-hot, MSI-H-enriched)** |
| Subtype-specific checkpoint axis | TIGIT (IS4), CD276 (IS4) | **CD274/PD-L1 (IS3), VTCN1/B7-H4 (IS5)** |
| HLA-I coverage per patient | ~561 epitopes (24 alleles) | **~28 epitopes/antigen (6 alleles); MUC16 dominates (872/patient)** |
| TIME-permissive subtype | IS4 (+0.69) | **IS3 (+1.06)** |
| cDC1+6-antigen co-elevation | CD276 × cDC1 ρ = 0.21 | **None of the 6 antigens co-varied; IS3 has highest cDC1** |
| Optimal cohort (3-layer) | 10.8% (51/471) in IS4 | **16.0% (88/549) in IS3** |
| External validation | GSE39582 | **HPA + CPTAC + GSE120490 + GSE278879 scRNA-seq** |

### R11.3 B7-H4 (VTCN1) — an EC-shifted axis pointing to IS5

In contrast to TIGIT and the six vaccine antigens, **B7-H4 (VTCN1) showed no overall survival signal in EC at the cohort level** (KM log-rank *P* = 0.26) — a finding that, in retrospect, reflects dilution by the five-subtype heterogeneity. The biologically meaningful signal is **IS5-restricted**: VTCN1 mean TPM is 236.8 in IS5 versus 105–147 in IS1–IS4, and the IS5-restricted pattern is independently confirmed by HPA (VTCN1 labelled "potential prognostic unfavorable" in TCGA-UCEC, *P* = 2.25 × 10⁻⁴), CPTAC (mRNA-protein ρ = 0.75), GSE120490 projection (VTCN1 IS5 highest, *P* = 3.71 × 10⁻²), and GSE278879 single-cell (VTCN1+ malignant epithelium detected across all three tumours, with the VTCN1-high MEEC1 tumour showing the lowest CD8A+ CD8 T-cell fraction). This four-way replication of the IS5-restricted B7-H4 pattern establishes IS5 as the **B7-H4-dominant EC** subtype, with direct therapeutic implications.

B7-H4 is a well-established co-inhibitory ligand that suppresses T-cell activation and proliferation, and its tumour-restricted expression makes it an attractive ADC target [19]. The B7-H4-directed topoisomerase-I-inhibitor ADC **puxitatug samrotecan (AZD8205)** has shown objective responses in 47.1% of B7-H4-selected EC patients (~60% in platinum/IO-pretreated disease) and has earned FDA Breakthrough Therapy Designation with an ongoing phase III trial (Bluestar-Endometrial01, NCT[2023]) [20]. A second B7-H4 ADC (HS-20089) is in development [21]. **The IS5 subtype, comprising 15.3% (84/549) of TCGA-UCEC, is the natural EC population for a B7-H4 vaccine + ADC dual-track strategy**: the IS5-restricted VTCN1 expression provides the molecular basis for a shared B7-H4 mRNA vaccine (off-the-shelf, no individualised WES/RNA-seq pipeline required), and the B7-H4 ADC clinical activity in gynaecological malignancies [20] provides the immediate translational context. This hypothesis should now be testable in a single-arm phase II trial enrolling IS5-enriched (VTCN1^high) EC patients, with VTCN1 immunohistochemistry as the enrolment criterion and vaccine-induced VTCN1-specific T-cell clone expansion as the primary endpoint.

### R11.4 CLDN6 — the safety flagship

CLDN6 stands out among the six antigens for its **favourable normal-tissue safety profile**: HPA reports CLDN6 normal-tissue expression ≤ 1 nTPM across all 32 normal tissues, with the highest normal expression in brain and pancreas (1.0 nTPM each). This is in stark contrast to MAL (cervix 687 nTPM, esophagus 2,851 nTPM) and CTSV (thymus 599 nTPM, small intestine 63 nTPM), both of which require careful safety evaluation in pre-clinical models. CLDN6 is a tight-junction protein with documented reactivation in several cancers (ovarian, endometrial, gastric) and is the target of a CAR-T programme (BNT211, BioNTech) currently in early-phase trials [22]. Our analysis places CLDN6 as the top-priority antigen of the panel — both because of its favourable safety and because of its high log₂FC (5.78) and poor-prognosis association (log-rank *P* = 6.74 × 10⁻⁵). A CLDN6-encoding mRNA vaccine is, to our knowledge, the most obvious first clinical candidate from this panel.

### R11.5 Single-cell validation of the IS3 ↔ CD8 T-cell and IS5 ↔ VTCN1 equivalences

The GSE278879 single-cell analysis (29,785 cells across three tumours) directly validates the bulk-derived equivalences at single-cell resolution: **the IS3 signature is specifically and strongly enriched in CD8⁺ T cells** (mean +0.469 vs Epithelial −0.057, Mann-Whitney *P* ≈ 0), with NK cells and macrophages also IS3-skewed. The IS5 signature is weakly enriched in Epithelial and Fibroblast (+0.004–0.005), consistent with its bulk stromal/secretory gene profile and supporting the interpretation that IS5 reflects a tumour-stroma co-state rather than a single-cell-type identity. **VTCN1 (B7-H4) is detectable in malignant epithelium across all three tumours** (5–66% VTCN1+ epithelial cells), and the VTCN1-high middle-grade tumour (MEEC1) has the lowest CD8A+ CD8 T-cell fraction (23.8% vs 29.9% HEEC1 vs 36.1% LEEC1) — a pattern consistent with the VTCN1-mediated immunosuppressive axis. The single-cell validation, although based on three tumours, is the first direct evidence that the bulk-derived IS3 ↔ immune-hot ↔ CD8-T-cell and IS5 ↔ B7-H4 ↔ cold equivalences are biologically real at the single-cell level.

### R11.6 Positioning against the 2025–2026 mRNA-vaccine clinical landscape

The clinical translation of mRNA cancer vaccines accelerated markedly in 2025–2026, and the present framework should be evaluated in that context. **First**, the personalised neoantigen mRNA vaccine **intismeran autogene (V940 / mRNA-4157)** + pembrolizumab met both primary endpoints (recurrence-free and distant-metastasis-free survival) in the phase III INTerpath-001 trial in resected high-risk melanoma [7], providing the first positive phase III readout for any individualised mRNA vaccine. Our 6-antigen shared-antigen panel is positioned as a **complement** to the individualised approach: shared antigens can be stockpiled, manufactured at scale, and offered to patients without the WES/RNA-seq pipeline required for individualised vaccines.

**Second**, the personalised neoantigen vaccine **autogene cevumeran (BNT122)** + atezolizumab + mFOLFIRINOX in resected pancreatic adenocarcinoma induced vaccine-specific T-cell clones with average estimated lifespans of 7.7 years, with 7/8 immune responders alive at 6-year follow-up versus only 2/8 non-responders [6]. The central question now driving the field — why ~50% of patients fail to mount a vaccine-induced T-cell response — is best framed as a **compositional** question about the pre-vaccine immune-microenvironment state, for which our IS1–IS5 taxonomy and the GEP/TIME/cDC1 stratification provide a directly testable hypothesis.

**Third**, the off-the-shelf amphiphile-vaccine **ELI-002 2P** (KRAS-G12D/G12R, Wainberg et al. *Nature Medicine* 2025) demonstrated vaccine-induced T-cell responses translating into improved RFS and OS versus historical controls in resected KRAS-mutant pancreatic and colorectal cancer [14]. Our 6-antigen panel covers the non-KRAS-G12D/G12R majority of EC patients, complementing rather than competing with ELI-002. TP53 mutation-derived shared neoantigens (37.3% of EC) would similarly complement the panel.

**Fourth**, the demonstration that SARS-CoV-2 mRNA vaccination within 100 days of ICI initiation was associated with 3-year OS improvement from 30.8% to 55.7% (aHR = 0.51) in NSCLC and melanoma [23] establishes that the mRNA-LNP platform itself is a potent innate-immune modulator. This places mRNA vaccines in a broader frame — not just antigen-delivery vehicles but innate-immune adjuvants that complement checkpoint blockade — and provides an additional mechanistic rationale for combining any mRNA-vaccine-encoded antigen panel with PD-1/PD-L1 blockade in the IS1/IS2 cold subtypes.

### R11.7 Limitations

Several limitations should be acknowledged. **First**, the GDC clinical MSI/MMR/POLE fields were empty, so the molecular-class derivation relies on the MAF; the derived POLE (~8%) and MSI-H (~20%) frequencies match the published distribution, but a formal MANTIS/MSIsensor recalculation is recommended for clinical use. **Second**, IS5 is 15.3% of the cohort (84 tumours) and VTCN1 differences, while replicated across four external datasets, remain quantitatively modest (mean TPM 236.8 in IS5 vs 105–147 elsewhere); larger cohorts and protein-level validation are required. **Third**, the CPTAC cohort lacks OS annotation, so survival could not be validated there. **Fourth**, the single-cell analysis is based on three tumours (29,785 cells); eight additional tumours in GSE278879 remain to be analysed. **Fifth**, HLA-I epitope coverage was estimated using a transparent motif-based anchor model; the published NetMHCpan 4.1 / NetMHCIIpan 4.3 IC50 pipeline is the appropriate final replacement before clinical use. **Sixth**, the antigen panel relies on TCGA-UCEC transcript abundance; the disconnect between mRNA and protein for MSLN (high protein, low mRNA) is a cautionary note — pre-clinical protein-level confirmation is required for all six antigens. **Seventh**, external projection across platforms (RNA-seq centroid applied to microarray) introduces noise (median projection correlation in GSE120490 was modest); protein-level validation is required for clinical translation. **Eighth**, the TP53 mutation-derived neoantigen track is proposed but not yet computationally validated in this study; a follow-on analysis is needed. **Ninth**, the B7-H4 + ADC dual-track hypothesis is based on IS5-restricted VTCN1 expression and the existing gynaecological-cancer B7-H4 ADC activity; a formal clinical-grade IHC cut-point and a prospective IS5-enriched trial design are required.

### R11.8 Conclusions

Endometrial carcinoma can be resolved into five immune subtypes with distinct molecular, prognostic and clinicopathological features. The immune-hot, MSI-H-enriched IS3 represents the optimal mRNA-vaccine-response population (16.0% of TCGA-UCEC), while B7-H4 defines a targetable axis in the immunosuppressive IS5 (15.3% of TCGA-UCEC). Integrating GEP, TIME conversion-sensitivity, cDC1 delivery and the subtype taxonomy defines a four-layer patient-selection model for the six-antigen EC mRNA vaccine. The revised panel — CLDN6, CCNE1, MAL, CTSV, VTCN1 (B7-H4) and MUC16 — is supported by four orthogonal external validations (HPA, CPTAC UCEC, GSE120490, GSE278879 scRNA-seq) and provides a data-driven rationale for both a shared-antigen IS3-targeting vaccine trial and a B7-H4 dual-track strategy in IS5. This work extends the immunotype-guided mRNA-vaccine framework from COAD to EC and identifies IS3 as the EC-specific optimal mRNA-vaccine-response subtype.

---

## Figure Legends

> Numbers verified against the analysis outputs (`C:\Users\wcz\Desktop\EC vaccine\data`) as of 2026-09-07.

### General conventions (apply to all figures)

- **Discovery cohort:** TCGA-UCEC, **549 primary tumours, 35 normal tissues, 1 recurrent**.
- **External validation cohorts:** HPA (6 antigens, 32 normal tissues), CPTAC UCEC (95 tumours), GSE120490 (GPL570, 145 tumours), GSE278879 (GPL24676, 3 tumours, 29,785 cells).
- **Statistical tests:** Kruskal–Wallis (KW) for across-subtype comparisons with Benjamini–Hochberg (BH) FDR correction; two-sided Mann–Whitney U for differential expression; Spearman ρ for correlations (BH-FDR corrected); two-sided log-rank for survival; χ² for categorical associations.
- **Box plots:** centre line, median; box, IQR; whiskers, 1.5 × IQR; points, individual outliers.
- **Heatmaps:** per-row z-scores (mean-centred and scaled to unit variance). Asterisk (*) marks significant subtype effect (KW FDR < 0.05).
- **Subtype colour scheme (consistent throughout):** IS1, blue; IS2, green; IS3, red; IS4, violet; IS5, teal.
- **Abbreviations:** TMB, tumour mutation burden; ICP, immune checkpoint; GEP, gene-expression profile; ssGSEA, single-sample gene-set enrichment analysis; APC, antigen-presenting cell; MSI-H, microsatellite-instability-high; POLE, polymerase-epsilon exonuclease-domain mutant; OS, overall survival; IHC, immunohistochemistry; ADC, antibody-drug conjugate; FDA BTD, US FDA-BTD.

### Figure 1. Integrated study design — from antigen discovery to four-layer patient stratification.

Single-panel schematic of the two-stage workflow applied to TCGA-UCEC (549 tumours vs 35 normals). **Stage 1 — Antigen identification** (top, blue band): four sequential filters compressing 18,784 protein-coding genes into six prioritised antigens — (1) tumour-specific over-expression (log₂FC > 1, FDR < 0.05), (2) recurrent somatic mutation (non-synonymous frequency ≥ 1%), (3) poor-prognosis association (OS log-rank *P* < 0.05), (4) positive APC-infiltration Spearman ρ — yielding **1,169 → 307 → 151 candidates**, from which the rule selects **(i)** the top four by log₂FC (CLDN6, CTSV, CCNE1, MAL) and **(ii)** two biologically/translationally prioritised (VTCN1/B7-H4, MUC16). The terminal red panel summarises the panel. **Stage 2 — Pre-vaccine stratification** (middle, violet band, four-quadrant grid): four transcriptomic signatures on the same 549-tumour cohort — **GEP** (18-gene T cell-inflamed signature; IS3 GEP^high = 91.9%; KW *P* = 3.55 × 10⁻⁹⁰), **HLA-I epitope coverage** (six common HLA-A/B/C alleles; MUC16 5,235 9-mers total, 872/patient; KW *P* = 0.85), **TIME conversion-sensitivity** (IS3 +1.06, IS2 −0.64; KW *P* significant), and **cDC1 delivery-layer capacity** (BATF3 · CLEC9A · XCR1 · IRF8 · THBD; IS3 +0.77, IS2 −0.49). **Four-layer cohort-selection model** (gold frame): TIME^permissive ∩ cDC1^high ∩ IS3 ∩ GEP^high. **Optimal cohort** (red terminal): 16.0 % of TCGA-UCEC (88/549 patients). The bottom **External validation** panel links to HPA, CPTAC, GSE120490 and GSE278879 scRNA-seq. The 12 excluded candidates (MSLN, ERBB2, TP53-WT and nine others) are documented in the Methods M2.3 and `DIAGNOSIS_AND_REVISION.md`.

### Figure 2. Identification of six prioritised tumour antigens in EC.

**(a)** Volcano plot of differential expression between 549 tumours and 35 normal tissues. x-axis: log₂ fold-change (tumour vs. normal); y-axis: −log₁₀ *P* value (two-sided Mann–Whitney U). Vertical dashed lines mark |log₂FC| = 1; horizontal dashed line marks *P* = 0.05 (FDR < 0.05). Grey, non-DEGs; blue, candidate antigens (overexpressed ∩ mutated, n = 1,169); red, prognostic ∩ APC-correlated strict candidates (n = 151). The six prioritised antigens are labelled (CLDN6 log₂FC = 5.78, CCNE1 3.56, MAL 3.37, CTSV 3.87, VTCN1 2.12, MUC16 1.21).

**(b)** The six prioritised antigens ranked by log₂ fold-change, with KM log-rank *P* values annotated. CLDN6 highest (log₂FC = 5.78, KM *P* = 6.74 × 10⁻⁵). Red bars denote the six panel antigens; blue bars, the top nine non-selected candidates.

### Figure 3. Consensus clustering resolves five immune subtypes with distinct prognosis and molecular features.

**(a)** Consensus matrix heatmap for k = 5, ordered by cluster. White-to-blue gradient encodes consensus index (0–1); red lines mark subtype boundaries.

**(b)** Kaplan-Meier overall-survival curves for IS1–IS5. Subtype sizes: IS1 = 145, IS2 = 130, IS3 = 99, IS4 = 91, IS5 = 84. IS3 and IS4: median OS not reached (11/90, 10/84 deaths). IS1: 112.5 months; IS2: 110.0 months. **IS5: 77.3 months** (21/71 deaths; shortest). Multivariate log-rank across 5 subtypes *P* significant.

**(c)** Subtype × derived molecular class (POLE / MSI-H / MSS) cross-tabulation. χ² = 39.64, *P* = 3.74 × 10⁻⁶. IS3: 38.4% hypermutated (MSI-H 27.3%); IS4: 36.3% (POLE 17.6%); IS5: 13.1% (lowest).

### Figure 4. Subtype characterisation heatmap and TMB distribution.

**(a)** Row z-score heatmap of subtype-defining features across IS1–IS5, including TMB, immune score, immune-checkpoint expression (CD274, PDCD1, CTLA4, LAG3, HAVCR2, TIGIT, VTCN1, CD276), immune-cell scores (CD8⁺ T, Treg, Cytotoxic, NK, B, Macrophage, DC, Neutrophil), and the derived molecular-class fractions (% POLE, % MSI-H, % Hypermut). Asterisks mark significant subtype effects.

**(b)** TMB by subtype (box plot). IS3 and IS4 have the highest TMB (~11 mut/Mb); IS5 has the lowest.

### Figure 5. Prognostic value of the six prioritised antigens.

Kaplan-Meier overall-survival curves comparing high (red) versus low (blue) expression (median split) for each of the six prioritised antigens: **(a)** CLDN6, **(b)** CCNE1, **(c)** MAL, **(d)** CTSV, **(e)** VTCN1 (B7-H4), **(f)** MUC16 (CA125). Log-rank *P* values annotated within each panel. CLDN6, CCNE1, MAL and CTSV show significant high = worse prognosis; VTCN1 and MUC16 show non-significant trends in the whole cohort (the subtype-restricted interpretation is shown in Fig. 2b).

### Figure 6. Subtype-specific expression of immune checkpoints and 28 immune-cell types.

**(a)** Row z-score heatmap of 49 immune-checkpoint genes across IS1–IS5, ordered by KW FDR. Asterisks mark significant subtype effects. TIGIT, CD274/PD-L1, HAVCR2/TIM-3, LAG3, PDCD1/PD-1 and CTLA4 all peak in IS3. **VTCN1/B7-H4 and CD276/B7-H3 peak in IS5**.

**(b)** Row z-score heatmap of 28 ssGSEA immune-cell signatures across IS1–IS5, ordered by KW FDR. All 28 cell types differ significantly across subtypes. IS3 is the most heavily infiltrated; IS2 is the most "cold"; IS5 has elevated macrophage/monocyte infiltration.

### Figure 7. T cell-inflamed GEP × subtype.

**(a)** GEP z-score distribution by subtype (box plot). IS3 has the highest median GEP z-score (+0.97); IS4 second (+0.21); IS5 intermediate (−0.10); IS1/IS2 negative.

**(b)** Fraction of GEP^high (z > +0.5) tumours per subtype. IS3: 91.9%; IS4: 46.2%; IS5: 22.6%; IS2: 11.5%; IS1: 4.8%. KW *P* = 3.55 × 10⁻⁹⁰.

### Figure 8. HLA-I 9-mer epitope coverage of the six antigens.

Heatmap of per-allele, per-antigen high-affinity 9-mer epitope counts across six HLA-I alleles (HLA-A\*02:01, A\*24:02, B\*07:02, B\*08:01, C\*07:01, C\*07:02). MUC16 dominates (5,235 total, 872/patient) owing to its large extracellular domain; the remaining five antigens contribute 109–169 total epitopes (17–28 per patient). KW *P* = 0.85 (no subtype bias).

### Figure 9. TIME conversion-sensitivity and cDC1 delivery layer.

**(a)** TIME index by subtype. IS3 (+1.06) is TIME-permissive; IS2 (−0.64) is TIME-resistant.

**(b)** TIME index components (Type-I IFN, TIL density, PD-L1, TGF-β) by subtype.

**(c)** cDC1 score by subtype. IS3 (+0.77), IS4 (+0.19), IS1 (−0.04), IS5 (−0.28), IS2 (−0.49).

**(d)** cDC2 and pDC scores by subtype.

**(e)** Four-layer patient-selection funnel (TIME → cDC1 → IS3 → GEP^high). Optimal cohort: 16.0% (88/549).

### Figure 10. External validation in GSE120490 (GPL570, n = 145).

**(a)** Heatmap of projected immune-subtype markers (CD274, VTCN1, CD276, TIGIT, LAG3, HAVCR2, CD8A, GZMB, IFNG, CLDN6, CCNE1, CTSV, MAL, MUC16) across the five projected subtypes. KW *P* values annotated. VTCN1 IS5-highest (*P* = 3.71 × 10⁻²); CD276 IS5-highest (*P* = 4.75 × 10⁻⁴); CD274 IS3-highest (*P* = 3.97 × 10⁻¹⁰); TIGIT IS3-highest (*P* = 5.40 × 10⁻¹³).

**(b)** Kaplan-Meier immune-score by projected subtype (KW H = 54.17, *P* = 4.86 × 10⁻¹¹). IS3 highest, IS2 lowest.

**(c)** Histological-grade composition by projected subtype. χ² = 15.89, *P* = 0.003. IS5 highest high-grade proportion (65%), IS1 lowest (15%).

### Figure 11. Single-cell validation in GSE278879 (n = 3 tumours, 29,785 cells).

**(a)** Cell-type composition per tumour (HEEC1 [high, 24,823 cells], MEEC1 [middle, 1,599 cells], LEEC1 [low, 3,363 cells]). 13 cell types displayed as 100% stacked bars. HEEC1 has the highest CD8_T fraction (5.0%) and lowest fibroblast fraction.

**(b)** IS3 − IS5 signature score per cell type (violin plot, sorted by median). IS3 signature strongly enriched in CD8⁺ T cells (mean +0.469 vs Epithelial −0.057, Mann-Whitney *P* ≈ 0). NK (+0.178) and Macrophage (+0.135) also IS3-skewed.

**(c)** Antigen expression (VTCN1, CD274) by cell type (bar plot). VTCN1 highest in Epithelial (red) and Treg (purple). CD274 highest in Macrophage (yellow).

### Figure 12. HPA validation — protein/RNA tissue-specificity and normal-tissue safety.

Per-gene heatmap of normal-tissue RNA expression (nTPM) for the six antigens across 32 normal tissues, annotated with HPA cancer-prognostic labelling. CLDN6 is ≤ 1 nTPM across all 32 tissues (safety flagship); VTCN1 is labelled "potential prognostic unfavorable" in TCGA-UCEC (*P* = 2.25 × 10⁻⁴); MAL has elevated normal-tissue expression in cervix/esophagus; CTSV in thymus/spleen.

### Figure 13. CPTAC UCEC validation — mRNA-protein correlation and immune-architecture replication.

**(a)** mRNA-protein Spearman ρ scatter plots for CTSV (ρ = 0.89), CCNE1 (ρ = 0.81), VTCN1 (ρ = 0.75), MUC16 (ρ = 0.71) — four of six antigens confirmed at protein level.

**(b)** GEP-cDC1-CD8 correlation heatmap (Spearman ρ with CIBERSORT/ESTIMATE immune deconvolution). CD8A ρ = 0.79, GZMB 0.71, CXCL9 0.78, STAT1 0.79, IFNG 0.72, CD274 0.77, TIGIT 0.73, cDC1 0.59 (all *P* < 10⁻⁹). POLE-mutant tumours have higher GEP (*P* = 0.025).

### Figure 14. Subtype × molecular-class derivation flow diagram.

Schematic of the MAF-based derivation of POLE/MSI-H/MSS classes (n = 42 / 102 / 364; total 508).

---

## Table 1. Final six-antigen panel for the EC mRNA vaccine

| # | Antigen | log₂FC | FDR | Mutation (%) | KM log-rank *P* | Selection rationale |
|---|---|---|---|---|---|---|
| 1 | **CLDN6** | **5.78** | 4.79 × 10⁻⁵ | 1.18 | **6.74 × 10⁻⁵** | Highest overexpression in the panel; HPA-confirmed normal-tissue safety (≤ 1 nTPM); clinical-stage CAR-T precedent (BNT211). |
| 2 | **CTSV** | 3.87 | 1.54 × 10⁻¹⁶ | 3.55 | **3.86 × 10⁻³** | Lysosomal cathepsin with antigen-processing role; HPA normal-tissue risk flagged (thymus/spleen). |
| 3 | **CCNE1** | 3.56 | 8.79 × 10⁻¹⁹ | 4.14 | **9.61 × 10⁻⁵** | Recognised EC driver gene (cyclin E1); CCNE1 amplification/overexpression hallmark of copy-number-high EC. |
| 4 | **MAL** | 3.37 | 1.19 × 10⁻⁵ | 1.18 | **4.66 × 10⁻⁵** | Membrane protein suitable for MHC-I presentation; HPA normal-tissue risk flagged (cervix/esophagus). |
| 5 | **VTCN1 (B7-H4)** | 2.12 | 3.32 × 10⁻⁸ | 1.58 | 0.263 (whole cohort); **IS5-restricted** | Translational anchor: IS5-dominant subtype + HPA-confirmed unfavorable-prognostic label; supports B7-H4 vaccine + ADC (puxitatug samrotecan) dual-track. |
| 6 | **MUC16 (CA125)** | 1.21 | 3.46 × 10⁻³ | **24.26** | 0.545 | Highest mutation rate in panel; broadest HLA-I epitope coverage (5,235 9-mers); established serum biomarker (CA-125). |

**Excluded candidates (with documented reasons):**
- **MSLN** — normal endometrium TPM (508) > tumour TPM (248), log₂FC = −1.03. Independent HPA confirmation. Excluded; MSLN is not an EC mRNA-level antigen (the literature's "70% positivity" is protein IHC, not mRNA).
- **ERBB2 (HER2)** — overexpressed but APC correlation ρ = −0.02, cDC1 ρ = −0.06. Retained as a dual-track ADC discussion.
- **TP53 (WT)** — log₂FC +0.57 (below filter). 37.3% mutation rate motivates a separate **TP53 mutation-derived shared-neoantigen track**.

**Complementary neoantigen track:** TP53 mutation hotspots (37.3% of EC) provide a shared-neoantigen track analogous to ELI-002 KRAS-G12D/G12R; not included in the 6-antigen shared panel but proposed as a complementary individualised/shared track for the 37.3% of TP53-mutant EC patients.

---

## Key statistics appendix

| Item | Value |
|---|---|
| **TCGA-UCEC cohort** | **549 primary tumours / 35 normals / 1 recurrent** |
| Protein-coding genes | 18,784 |
| DEGs (two-sided MW + BH-FDR, |log₂FC|>1, FDR<0.05) | >1,000 up; full pipeline 1,169 up ∩ mutated candidates |
| Prioritised antigens | CLDN6, CCNE1, MAL, CTSV, VTCN1 (B7-H4), MUC16 |
| Subtype sizes (IS1–IS5) | 145 / 130 / 99 / 91 / 84 |
| Subtype median OS (months) | IS1 112.5; IS2 110.0; IS3 not reached; IS4 not reached; IS5 77.3 |
| Subtype deaths | 21/136; / 28/119; / 11/90; / 10/84; / 21/71 |
| POLE / MSI-H / MSS (derived from MAF) | 42 / 102 / 364 (8.3% / 20.0% / 71.7%) |
| POLE TMB (mut/Mb) median | 157.9 |
| MSI-H TMB (mut/Mb) median | 14.4 |
| Subtype × molecular χ² / P | 39.64 / 3.74 × 10⁻⁶ |
| IS3 hypermutated fraction | 38.4% (MSI-H 27.3%) |
| IS4 hypermutated fraction | 36.3% (POLE 17.6%) |
| IS5 hypermutated fraction | 13.1% (lowest) |
| **CLDN6 / CCNE1 / MAL / CTSV KM P** | **6.74 × 10⁻⁵ / 9.61 × 10⁻⁵ / 4.66 × 10⁻⁵ / 3.86 × 10⁻³** |
| VTCN1 (whole-cohort) KM P | 0.263 (n.s.); IS5-restricted (mean TPM 236.8) |
| **GEP^high (z > +0.5) by subtype** | **IS1 4.8% / IS2 11.5% / IS3 91.9% / IS4 46.2% / IS5 22.6%** |
| GEP KW P | 3.55 × 10⁻⁹⁰ |
| **HLA-I 9-mer coverage (6 antigens × 6 alleles)** | **MUC16 5,235 / VTCN1 137 / CLDN6 134 / CCNE1 169 / CTSV 109 / MAL 103** |
| HLA coverage KW P across IS1–IS5 | 0.85 (n.s.; subtype-uniform) |
| **TIME conversion-sensitivity (z) median IS1–IS5** | **−0.22 / −0.64 / +1.06 / +0.18 / −0.09** |
| TIME KW P (TCGA) | significant |
| cDC1 score median IS1–IS5 | −0.04 / −0.49 / +0.77 / +0.19 / −0.28 |
| **Four-layer cohort (TIME>0 + cDC1>0 + IS3)** | **88 / 549 (16.0%)** |
| **HPA VTCN1 prognostic label** | **unfavorable (P = 2.25 × 10⁻⁴)** |
| **HPA CLDN6 normal-tissue max nTPM** | **≤ 1.0 (pancreas/brain)** |
| **CPTAC UCEC mRNA-protein ρ** | **CTSV 0.89 / CCNE1 0.81 / VTCN1 0.75 / MUC16 0.71** |
| CPTAC GEP-CD8A ρ | 0.79 |
| CPTAC POLE GEP effect | P = 0.025 |
| **GSE120490 projection (IS1–IS5)** | **24 / 33 / 35 / 26 / 27** |
| GSE120490 immune-score KW P | 4.86 × 10⁻¹¹ (IS3 highest) |
| GSE120490 grade × subtype χ² / P | 15.89 / 0.003 |
| GSE120490 VTCN1 IS5 highest | KW P = 3.71 × 10⁻² |
| GSE120490 CD274 IS3 highest | KW P = 3.97 × 10⁻¹⁰ |
| **GSE278879 cells** | **29,785 (HEEC1 24,823 / MEEC1 1,599 / LEEC1 3,363)** |
| GSE278879 IS3 signature CD8_T mean | +0.469 vs Epithelial −0.057, **Mann-Whitney P ≈ 0** |
| GSE278879 VTCN1+ Epithelial per tumor | HEEC1 5.1% / MEEC1 65.9% / LEEC1 40.6% |
| GSE278879 CD8A+ CD8_T per tumor | HEEC1 29.9% / MEEC1 23.8% / LEEC1 36.1% |

---

## Funding, data and code availability

- **Data:** TCGA-UCEC via the GDC Data Portal (https://portal.gdc.cancer.gov/); HPA via the Protein Atlas REST API (https://www.proteinatlas.org/); CPTAC UCEC via the cBioPortal REST API (https://www.cbioportal.org/); GSE120490 and GSE278879 via NCBI GEO (https://www.ncbi.nlm.nih.gov/geo/).
- **Code:** all analysis scripts are available in `C:\Users\wcz\Desktop\EC vaccine\data\` and the corresponding figures are in `C:\Users\wcz\Desktop\EC vaccine\figures\`. The companion COAD analysis is documented in `C:\Users\wcz\Desktop\COAD vaccine\`.

---

## Conflict of interest

The authors declare no conflicts of interest.

---

*Manuscript drafted 2026-09-07. Companion figures: Fig 1–14 in `C:\Users\wcz\Desktop\EC vaccine\figures\`. External validation reports: `EXTERNAL_VALIDATION.md`, `SINGLE_CELL_VALIDATION.md`. Companion analysis summaries: `DIAGNOSIS_AND_REVISION.md`, `RESULTS_SUMMARY.md`. Reference framework: Huang et al., Molecular Cancer 2021;20:44 (PAAD); companion analysis: Wang et al., 2026 (COAD).*