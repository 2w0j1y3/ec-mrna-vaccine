# Identification of mRNA Vaccine Antigens and Immune Subtypes in Endometrial Carcinoma

> Companion repository for the manuscript "Identification of mRNA Vaccine Antigens and Immune Subtypes in Endometrial Carcinoma: An Integrated Framework from TCGA-UCEC to Single-Cell Resolution".

This project extends the COAD mRNA-vaccine framework (companion analysis at [`COAD-vaccine`](https://github.com/2w0j1y3/COAD-vaccine)) to **endometrial carcinoma (EC)**, integrating TCGA-UCEC transcriptomics with somatic-mutation, proteomic (CPTAC), microarray (GSE120490) and single-cell RNA-seq (GSE278879) data to nominate a 6-antigen mRNA-vaccine panel, define 5 immune subtypes, and propose a 4-layer patient-selection algorithm.

---

## Key findings

| Item | Value |
|---|---|
| TCGA-UCEC discovery cohort | 549 primary tumours / 35 normal / 1 recurrent |
| **Revised 6-antigen panel** | **CLDN6 · CCNE1 · MAL · CTSV · VTCN1 (B7-H4) · MUC16 (CA125)** |
| 5 immune subtypes (IS1–IS5) | 145 / 130 / 99 / 91 / 84 tumours |
| **Optimal mRNA-vaccine subtype** | **IS3** — MSI-H enriched (27.3%), GEP^high 91.9%, median OS not reached |
| **B7-H4 dual-track subtype** | **IS5** — VTCN1 mean TPM 236.8, median OS 77.3 mo |
| 4-layer optimal cohort | TIME^permissive ∩ cDC1^high ∩ IS3 = **88 / 549 (16.0%)** |
| External validation | HPA · CPTAC UCEC (n = 95) · GSE120490 (n = 145) · GSE278879 scRNA-seq (n = 3 tumours, 29,785 cells) |

The panel excludes MSLN (normal endometrium TPM > tumour), ERBB2 (APC ρ = −0.02) and TP53-WT (log₂FC +0.57). TP53 mutation hotspots (37.3% of EC) are proposed as a separate shared-neoantigen track.

> ⚠️ Note: Figure numbering is **sequential 1–14** throughout the manuscript and tables. The physical figure filenames in `figures/` retain their historical naming (e.g. `fig13_*.png` is referenced as **Figure 7**). See [Figure renumbering map](#figure-renumbering-map) below.

---

## Repository layout

```
EC-mRNA-vaccine/
├── MANUSCRIPT.md                       # Submission-ready manuscript (Markdown)
├── EC_mRNA_Vaccine_Manuscript.docx    # Generated Word document (Stage-3 docx)
├── DIAGNOSIS_AND_REVISION.md          # Diagnosis of "unsatisfactory" first run + revision strategy
├── RESULTS_SUMMARY.md                 # Complete analysis result summary
├── EXTERNAL_VALIDATION.md             # HPA / CPTAC / GSE120490 reports
├── SINGLE_CELL_VALIDATION.md          # GSE278879 single-cell validation report
│
├── data/                              # Analysis scripts + per-step CSV outputs
│   ├── analysis_antigen.py            #   4-filter antigen identification funnel
│   ├── analysis_subtype.py            #   Consensus K-means (k = 5) on 176 immune genes
│   ├── analysis_gep.py                #   T cell-inflamed GEP (18 genes)
│   ├── analysis_hla_epitope.py        #   HLA-I 9-mer epitope coverage (6 alleles)
│   ├── analysis_time_cdc.py           #   TIME conversion-sensitivity + cDC1 delivery layer
│   ├── download_expr.py / merge_expr.py   #   GDC RNA-seq acquisition
│   ├── process_maf.py / derive_molecular.py / finalize_molecular.py
│   │                                 #   MAF-based POLE / MSI-H / MSS derivation
│   ├── parse_geo_ec.py                #   GSE120490 family SOFT parsing
│   ├── validate_geo_ec.py             #   Nearest-centroid projection + clinical validation
│   ├── fetch_cptac.py / validate_cptac.py
│   ├── check_hpa.py / summarize_hpa.py / hpa_antigen_check.csv
│   ├── derive_subtype_signatures.py   #   TCGA bulk → IS3/IS5 gene signatures
│   ├── parse_and_annotate_sc.py       #   10x snRNA-seq parser + cell-type annotation
│   ├── analyze_sc_validation.py       #   Single-cell validation + Fig 17
│   ├── make_figures.py                #   Master figure-generation script (fig1–6, 13–17)
│   ├── *_signature.csv / *_validation.csv    # Per-stage intermediate CSVs (committed)
│   └── *_family.soft.gz / *.npz / *.tar.gz   # Raw downloads (gitignored, re-fetchable)
│
├── figures/                           # Publication-grade figures (PDF / SVG / PNG)
│   ├── fig1_antigen_identification.*  #   Volcano + revised 6-antigen bar
│   ├── fig2_immune_subtypes.*         #   Consensus matrix + KM (median OS annotated)
│   ├── fig3_subtype_features.*       #   Heatmap (TMB/immune/ICP + % POLE/% MSI-H/% Hypermut) + TMB box
│   ├── fig4_antigen_survival.*        #   2×3 KM grid for revised 6 antigens
│   ├── fig5_*                         #   (skipped per sequential renumbering)
│   ├── fig6_immune_checkpoint.*       #   49-ICP + 28-cell ssGSEA heatmaps
│   ├── fig13_gep.*                    #   Fig 7 — GEP × subtype
│   ├── fig14_hla_epitopes.*           #   Fig 8 — HLA-I coverage
│   ├── fig15_time_cdc1.*              #   Fig 9 — TIME/cDC1
│   ├── fig16_external_validation.*    #   Fig 10 — GSE120490
│   ├── fig17_sc_validation.*          #   Fig 11 — GSE278879
│   └── _fig*_preview*.png             #   Preview thumbnails (gitignored or removed)
│
└── output/                            # tencent-docx pipeline intermediate state (gitignored)
```

---

## Methods at a glance

1. **Data acquisition** — TCGA-UCEC RNA-seq (549 tumours / 35 normals / 1 recurrent) from GDC, plus CPTAC UCEC (n = 95), GSE120490 GPL570 (n = 145) and GSE278879 10x snRNA-seq (3 tumours, 29,785 cells).
2. **POLE / MSI-H derivation** — GDC clinical fields are 100% empty; derive from MC3F MAF:
   - **POLE** (n = 42, 8.3%): exonuclease-domain hotspots → median TMB 157.9 mut/Mb
   - **MSI-H** (n = 102, 20.0%): non-POLE with n_indel ≥ 80 & indel/SNV > 0.2 → TMB 14.4 mut/Mb
   - **MSS** (n = 364, 71.7%): remainder → TMB 1.3 mut/Mb
3. **Antigen identification** — 4-filter funnel on 18,784 protein-coding genes:
   `overexpressed (log₂FC > 1, FDR < 0.05) → mutated (≥ 1%) → prognostic (log-rank P < 0.05) → APC-positive (Spearman ρ > 0)` → 1169 → 307 → 151 strict candidates; final 6 selected with transparent rule (top 4 by log₂FC + 2 translational anchors).
4. **Immune subtyping** — Consensus K-means (k = 5, 80% subsampling, 100 runs) on 176 immune genes.
5. **Pre-vaccine stratification** — 18-gene T cell-inflamed GEP (Ayers 2017), HLA-I 9-mer motif model, TIME conversion-sensitivity index (IFN + TIL + PD-L1 − TGF-β), cDC1 delivery-layer (BATF3 · CLEC9A · XCR1 · IRF8 · THBD).
6. **Single-cell validation** — Per-sample matrix.tar.gz from GEO, 10x parsing, QC (UMI ≥ 500), z-score marker annotation, IS3/IS5 signature scoring.
7. **External validation** — HPA tissue-safety + prognostic label, CPTAC mRNA-protein Spearman, GSE120490 nearest-centroid projection + clinical endpoints, GSE278879 per-tumour IS3/IS5 + VTCN1/CD8A composition.

---

## Reproducibility

```bash
# 1. Create venv (managed Python 3.11+)
python -m venv .venv
source .venv/bin/activate   # or .venv\Scripts\activate on Windows

# 2. Install dependencies (numpy, pandas, scipy, scikit-learn, lifelines, matplotlib, requests)
pip install -r requirements.txt

# 3. Re-fetch TCGA / GEO / CPTAC data (≥ 6 GB raw, not in repo)
python data/download_expr.py
python data/process_maf.py
python data/parse_geo_ec.py
# ... etc.

# 4. Run analysis pipeline (each is independent)
python data/analysis_antigen.py
python data/analysis_subtype.py
python data/analysis_gep.py
python data/analysis_hla_epitope.py
python data/analysis_time_cdc.py
python data/analyze_sc_validation.py

# 5. Render figures
python data/make_figures.py

# 6. Render Word document (10x pipeline via tencent-docx skill)
#    Stage 1: local_md passthrough
#    Stage 2: design-token + doc-typeset + html-review
#    Stage 3: html-to-docx
```

---

## Figure renumbering map

Per the user's request to "renumber figures starting from 1 in sequential order":

| New # | Old filename (in `figures/`) | Title |
|---|---|---|
| 1 | fig1_antigen_identification | Integrated study design |
| 2 | fig2_immune_subtypes | Six prioritised antigens identification (R1 part) |
| 3 | fig3_subtype_features | Consensus clustering (k = 5) |
| 4 | fig4_antigen_survival | Subtype characterisation heatmap + TMB |
| 5 | fig5_* | Prognostic value of six antigens |
| 6 | fig6_immune_checkpoint | Subtype-specific ICP + 28 immune cells |
| **7** | fig13_gep | T cell-inflamed GEP × subtype |
| **8** | fig14_hla_epitopes | HLA-I 9-mer epitope coverage |
| **9** | fig15_time_cdc1 | TIME conversion-sensitivity + cDC1 |
| **10** | fig16_external_validation | GSE120490 cross-platform validation |
| **11** | fig17_sc_validation | GSE278879 single-cell validation |
| **12** | (was Fig S1) | HPA tissue-safety + prognostic label |
| **13** | (was Fig S2) | CPTAC mRNA-protein correlation |
| **14** | (was Fig S3) | Subtype × molecular-class derivation flow |

All cross-references in the manuscript text and figure legends have been updated to the new sequential numbers.

---

## Companion analyses

- **[COAD vaccine](https://github.com/2w0j1y3/COAD-vaccine)** — Colon adenocarcinoma. Same framework but with COAD-specific 6 antigens (MMP3/STC2/ULBP2/MAGEA12/CD276/CHEK1), IS4 as optimal-vaccine subtype.

---

## Citation

If you use this framework, please cite the underlying methodology:

> Huang X. et al. (2021). *Identifying therapeutic targets in cancer using genomic data.* Molecular Cancer 20:44.

For EC-specific adaptations and the B7-H4 dual-track IS5 strategy, see the manuscript draft at `MANUSCRIPT.md` (submission date TBA).

---

## License

TBA — academic non-commercial research use unless otherwise noted.
