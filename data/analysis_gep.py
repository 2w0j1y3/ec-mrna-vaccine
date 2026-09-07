"""EC GEP (T-cell-inflamed gene expression profile) analysis。

复现 Ayers JCI 2017: 18-gene IFN-γ signature predicting ICI response.

EC-specific: 与免疫亚型、POLE/MSI 关联分析。

输出：
- ucec_gep.csv (每个样本的 GEP z-score)
- ucec_gep_by_subtype.csv
"""
import os, numpy as np, pandas as pd
from scipy import stats

DATA = os.path.dirname(os.path.abspath(__file__))

# Ayers JCI 2017 - 18-gene IFN-γ signature
# Original: IFNG, STAT1, CCR5, CXCL9, CXCL10, CXCL11, IDO1, PRF1, GZMA, GZMB,
# CD8A, CD8B, PDCD1, PDCD1LG2 (PD-L2), CD274 (PD-L1), LAG3, HAVCR2 (TIM-3), TIGIT
GEP_GENES = ["IFNG", "STAT1", "CCR5", "CXCL9", "CXCL10", "CXCL11", "IDO1",
             "PRF1", "GZMA", "GZMB", "CD8A", "CD8B", "PDCD1", "PDCD1LG2",
             "CD274", "LAG3", "HAVCR2", "TIGIT"]


def main():
    m = pd.read_csv(os.path.join(DATA, "ucec_tpm_matrix.csv"), index_col=0)
    smap = pd.read_csv(os.path.join(DATA, "ucec_sample_map.csv"))
    tumor = [s for s in smap[smap["sample_type"] == "Primary Tumor"]["sample_id"] if s in m.columns]
    print(f"肿瘤样本: {len(tumor)}")

    found = [g for g in GEP_GENES if g in m.index]
    missing = [g for g in GEP_GENES if g not in m.index]
    print(f"GEP 基因可用: {len(found)}/{len(GEP_GENES)}")
    if missing:
        print(f"  缺失: {missing}")

    # 计算 GEP z-score
    log2_tpm = np.log2(m.loc[found, tumor].astype(float) + 1)
    gep_per_gene_mean = log2_tpm.mean(axis=1)
    gep_per_gene_std = log2_tpm.std(axis=1)
    gep_z = log2_tpm.subtract(gep_per_gene_mean, axis=0).divide(gep_per_gene_std + 1e-8, axis=0)
    gep = gep_z.mean(axis=0)

    out = pd.DataFrame({"sample_id": tumor, "GEP_z": gep.values})
    out.to_csv(os.path.join(DATA, "ucec_gep.csv"), index=False)
    print(f"\nGEP 统计:")
    print(f"  Mean: {gep.mean():.3f}, Std: {gep.std():.3f}")
    print(f"  Min: {gep.min():.3f}, Max: {gep.max():.3f}")
    print(f"  z > +0.5: {(gep > 0.5).sum()} ({(gep > 0.5).mean()*100:.1f}%)")
    print(f"  z > +1.0: {(gep > 1.0).sum()} ({(gep > 1.0).mean()*100:.1f}%)")

    # 按免疫亚型
    subtypes = pd.read_csv(os.path.join(DATA, "ucec_subtypes.csv"))
    out = out.merge(subtypes, on="sample_id")
    by_sub = out.groupby("subtype")["GEP_z"].agg(["mean", "median", "std", "count"]).round(3)
    print(f"\nGEP × 免疫亚型:")
    print(by_sub)
    by_sub.to_csv(os.path.join(DATA, "ucec_gep_by_subtype.csv"))

    # Kruskal-Wallis
    groups = [g["GEP_z"].values for _, g in out.groupby("subtype")]
    h, p = stats.kruskal(*groups)
    print(f"\nKruskal-Wallis H={h:.2f}, P={p:.2e}")

    # 按 TCGA 分子亚型
    clin = pd.read_csv(os.path.join(DATA, "ucec_clinical.csv"))
    sample_to_case = dict(zip(smap["sample_id"], smap["case_id"]))
    out["case_id"] = out["sample_id"].map(sample_to_case)
    merged = out.merge(clin[["submitter_id", "msi_status", "mmr_status", "pole_mutation"]],
                       left_on="case_id", right_on="submitter_id", how="left")

    print(f"\nGEP × MSI 状态:")
    for msi in ["MSS", "MSI-H", "MSI-L"]:
        sub = merged[merged["msi_status"] == msi]
        if len(sub) > 5:
            print(f"  {msi}: n={len(sub)}, mean GEP_z = {sub['GEP_z'].mean():.3f}")

    print(f"\nGEP × POLE 突变:")
    for pole in ["yes", "no"]:
        sub = merged[merged["pole_mutation"] == pole]
        if len(sub) > 5:
            print(f"  POLE {pole}: n={len(sub)}, mean GEP_z = {sub['GEP_z'].mean():.3f}")

    print("\n=== GEP 分析完成 ===")


if __name__ == "__main__":
    main()
