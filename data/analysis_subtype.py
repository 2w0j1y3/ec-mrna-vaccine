"""EC (UCEC) 免疫亚型鉴定。

EC 特殊考虑：
1. TCGA 4 分子亚型：POLE (超突变) / MSI-H / CN-low (MSS) / CN-high (serous-like)
2. 我们的免疫亚型（IS1-IS5）作为新分层维度，独立于分子亚型
3. 比较免疫亚型与分子亚型的对应关系

输出：
- ucec_consensus_matrix.csv
- ucec_subtypes.csv
- ucec_subtype_features.csv
- ucec_subtype_summary.csv
- ucec_subtype_vs_molecular.csv
"""
import os
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.cluster import KMeans
from collections import Counter

DATA = os.path.dirname(os.path.abspath(__file__))

# 同 COAD 的免疫基因列表
IMMUNE_GENES = [
    "CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1", "IDO2",
    "CD80", "CD86", "CD28", "ICOS", "ICOSLG", "TNFRSF9", "TNFRSF4", "TNFSF4", "TNFRSF18",
    "TNFSF9", "CD40", "CD40LG", "CD27", "CD70", "VTCN1", "BTLA", "CD244", "CD160", "VSIR",
    "ENTPD1", "NT5E", "LILRB1", "LILRB2", "SIGLEC7", "LGALS9", "CD276", "NCR3LG1",
    "CD3D", "CD3E", "CD3G", "CD247", "CD4", "CD8A", "CD8B", "CD2", "CD5", "CD7",
    "ZAP70", "LCK", "ITK", "TRAC", "TRBC1", "TRBC2",
    "FOXP3", "IL2RA", "IL7R", "CCR4", "CCR8", "IKZF2",
    "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "GNLY", "FASLG", "NKG7", "IFNG", "TNF",
    "KLRK1", "KLRD1", "KLRC1", "KLRC2", "KLRB1", "NCR1", "NCR2", "NCR3", "FCGR3A",
    "KIR2DL1", "KIR2DL3", "KIR3DL1",
    "CD19", "CD79A", "CD79B", "MS4A1", "CD22", "CR2", "PAX5", "BANK1", "BLK", "FCRL1",
    "SDC1", "TNFRSF17", "MZB1", "XBP1",
    "CD68", "CD163", "MRC1", "CSF1R", "CD14", "ITGAM", "ITGAX", "MSR1", "MARCO",
    "NOS2", "IL12A", "IL12B", "ARG1", "IL10", "TGFB1", "CCL18", "CHI3L1",
    "CD1C", "CLEC9A", "CLEC10A", "ITGAE", "BATF3", "FLT3", "THBD", "LAMP3", "NRP1",
    "FCGR3B", "CSF3R", "S100A8", "S100A9", "CEACAM8",
    "KIT", "TPSAB1", "CPA3", "MS4A2", "HDC",
    "IL2", "IL4", "IL6", "IL7", "IL13", "IL15", "IL17A", "IL21", "IL23A",
    "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CCL2", "CCL3", "CCL4", "CCL5", "CCL19", "CCL21",
    "CXCR3", "CCR5", "CCR7", "IL2RG",
    "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
    "HLA-DQA1", "HLA-DQB1", "B2M", "TAP1", "TAP2", "TAPBP", "CD74", "CIITA",
    "STAT1", "STAT2", "IRF1", "IRF7", "IRF8", "MX1", "OAS1", "OAS2", "ISG15",
    "IFI27", "IFI44", "IFIT1", "IFIT2", "IFIT3", "RSAD2", "GBP1", "GBP2",
]


def main():
    print("载入表达矩阵 ...")
    m = pd.read_csv(os.path.join(DATA, "ucec_tpm_matrix.csv"), index_col=0)
    smap = pd.read_csv(os.path.join(DATA, "ucec_sample_map.csv"))
    tumor = [s for s in smap[smap["sample_type"] == "Primary Tumor"]["sample_id"] if s in m.columns]
    print(f"肿瘤样本: {len(tumor)}")

    immune_found = [g for g in IMMUNE_GENES if g in m.index]
    print(f"免疫基因可用: {len(immune_found)}/{len(IMMUNE_GENES)}")
    pd.DataFrame({"gene": immune_found}).to_csv(os.path.join(DATA, "ucec_immune_genes.csv"), index=False)

    X = np.log2(m.loc[immune_found, tumor].astype(float) + 1).T.values
    X = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-8)
    print(f"聚类矩阵: {X.shape}")

    # ---------- 一致性聚类 (k=5) ----------
    print("一致性聚类 (k=5, 100 次, 80% 子抽样) ...")
    N = X.shape[0]
    K = 5
    n_iter = 100
    consensus = np.zeros((N, N))
    counts = np.zeros((N, N))
    subsample = int(N * 0.8)

    rng = np.random.RandomState(42)
    for it in range(n_iter):
        idx = rng.choice(N, subsample, replace=False)
        Xsub = X[idx]
        km = KMeans(n_clusters=K, n_init=10, random_state=it, max_iter=300)
        labels = km.fit_predict(Xsub)
        mask = (labels[:, None] == labels[None, :])
        consensus[np.ix_(idx, idx)] += mask
        counts[np.ix_(idx, idx)] += 1

    consensus = consensus / np.maximum(counts, 1)
    np.fill_diagonal(consensus, 1.0)
    consensus_df = pd.DataFrame(consensus, index=tumor, columns=tumor)
    consensus_df.to_csv(os.path.join(DATA, "ucec_consensus_matrix.csv"))
    print("共识矩阵已保存")

    # 层次聚类
    dist = 1 - consensus
    Z = linkage(squareform(dist, checks=False), method="average")
    labels = fcluster(Z, t=K, criterion="maxclust")
    order = [c for c, _ in Counter(labels).most_common()]
    remap = {old: new + 1 for new, old in enumerate(order)}
    subtype = np.array([remap[l] for l in labels])

    subtype_df = pd.DataFrame({
        "sample_id": tumor,
        "subtype": [f"IS{s}" for s in subtype],
    }).sort_values("subtype")
    subtype_df.to_csv(os.path.join(DATA, "ucec_subtypes.csv"), index=False)
    print("亚型分布:")
    print(subtype_df["subtype"].value_counts().sort_index())

    # ---------- TMB ----------
    print("\n计算 TMB ...")
    # Use the MAF freq file
    mf = pd.read_csv(os.path.join(DATA, "ucec_mut_freq.csv"))
    # Need per-patient TMB - load from raw MAF if available
    # For now, aggregate per-patient from MAF files (similar to process_maf.py)
    import gzip, glob
    case_mut_counts = {}
    for fpath in sorted(glob.glob(os.path.join(DATA, "maf", "*.maf.gz"))):
        try:
            with gzip.open(fpath, "rt", encoding="utf-8", errors="ignore") as fh:
                header = None
                i_class = i_patient = None
                for line in fh:
                    if line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if header is None:
                        header = parts
                        i_class = header.index("Variant_Classification") if "Variant_Classification" in header else None
                        i_patient = header.index("Tumor_Sample_Barcode")
                        continue
                    if i_class is not None and parts[i_class] in ("Silent", "RNA", "5'UTR", "3'UTR", "5'Flank", "3'Flank", "IGR", "Intron"):
                        continue
                    tb = parts[i_patient][:12]
                    case_mut_counts[tb] = case_mut_counts.get(tb, 0) + 1
        except Exception:
            continue
    # EC WES covered region ~38 Mb (TCGA standard)
    tmb = {c: n / 38.0 for c, n in case_mut_counts.items()}
    sample_to_case = dict(zip(smap["sample_id"], smap["case_id"]))

    # ---------- 亚型特征 ----------
    print("计算亚型特征 ...")
    checkpoints = ["CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1", "VTCN1", "CD276", "PDCD1LG2"]
    checkpoints = [g for g in checkpoints if g in m.index]

    cell_markers = {
        "CD8 T": ["CD8A", "CD8B"],
        "Treg": ["FOXP3", "IL2RA", "IKZF2"],
        "Cytotoxic": ["GZMA", "GZMB", "PRF1", "IFNG"],
        "NK": ["KLRK1", "KLRD1", "NCR1", "NKG7"],
        "B cell": ["CD19", "MS4A1", "CD79A"],
        "Macrophage": ["CD68", "CD163", "MRC1"],
        "DC": ["CD1C", "CLEC9A", "ITGAX", "BATF3"],
        "Neutrophil": ["FCGR3B", "CSF3R", "S100A8"],
    }

    def gene_score(sample_ids, genes):
        genes = [g for g in genes if g in m.index]
        return np.log2(m.loc[genes, sample_ids].astype(float) + 1).mean(axis=0)

    features = subtype_df.copy()
    features["TMB"] = [tmb.get(sample_to_case.get(s, ""), np.nan) for s in features["sample_id"]]
    features["immune_score"] = gene_score(features["sample_id"], immune_found).values
    for cp in checkpoints:
        features[cp] = m.loc[cp, features["sample_id"]].astype(float).values
    for cell, genes in cell_markers.items():
        features[f"{cell}_score"] = gene_score(features["sample_id"], genes).values

    features.to_csv(os.path.join(DATA, "ucec_subtype_features.csv"), index=False)

    # 汇总各亚型特征均值
    summ = features.groupby("subtype").agg({
        "TMB": "median", "immune_score": "mean",
        **{cp: "mean" for cp in checkpoints},
        **{f"{c}_score": "mean" for c in cell_markers},
    }).round(3)
    summ.to_csv(os.path.join(DATA, "ucec_subtype_summary.csv"))
    print("\n各亚型特征汇总:")
    print(summ.to_string())

    # ---------- 免疫亚型 vs TCGA 分子亚型 ----------
    print("\n免疫亚型 vs TCGA 分子亚型（msi_status / mmr_status / pole_mutation）...")
    clin = pd.read_csv(os.path.join(DATA, "ucec_clinical.csv"))
    # Sample to case
    features["case_id"] = features["sample_id"].map(sample_to_case)
    merged = features.merge(clin[["submitter_id", "msi_status", "mmr_status", "pole_mutation", "tumor_grade", "ajcc_pathologic_stage"]],
                             left_on="case_id", right_on="submitter_id", how="left")

    crosstab = pd.crosstab(merged["subtype"], merged["msi_status"].fillna("NA"), margins=True)
    print("亚型 × MSI:")
    print(crosstab)

    crosstab_mmr = pd.crosstab(merged["subtype"], merged["mmr_status"].fillna("NA"), margins=True)
    print("\n亚型 × MMR:")
    print(crosstab_mmr)

    crosstab_grade = pd.crosstab(merged["subtype"], merged["tumor_grade"].fillna("NA"), margins=True)
    print("\n亚型 × 分级:")
    print(crosstab_grade)

    crosstab_grade.to_csv(os.path.join(DATA, "ucec_subtype_vs_grade.csv"))
    crosstab_mmr.to_csv(os.path.join(DATA, "ucec_subtype_vs_mmr.csv"))
    crosstab.to_csv(os.path.join(DATA, "ucec_subtype_vs_msi.csv"))

    print("\n=== 免疫亚型分析完成 ===")


if __name__ == "__main__":
    main()
