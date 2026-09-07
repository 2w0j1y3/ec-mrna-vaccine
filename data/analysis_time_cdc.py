"""EC TIME (tumor immune microenvironment) + cDC1 cross-presentation analysis.

定义：
1. TIME 4 组分：
   - IFN (11 基因): IFNG, STAT1, IRF1, IRF7, MX1, OAS1, ISG15, IFI27, IFIT1, GBP1, CXCL10
   - TIL (8 基因): CD8A, CD8B, GZMA, GZMB, PRF1, GNLY, NKG7, CD3E
   - PD-L1 (1 基因): CD274
   - TGF-β (9 基因): TGFB1, TGFB2, TGFB3, ACVR1, ACVR1B, BMPR1A, BMPR1B, SMAD3, SMAD4
2. TIME 转换敏感指数 = (z(IFN) + z(TIL) + z(PD-L1))/3 - z(TGF-β)
3. cDC1 markers (5 基因): BATF3, CLEC9A, XCR1, IRF8, THBD

输出：
- ucec_time.csv (4 组分 + TIME index)
- ucec_cdc.csv (cDC1 score)
- ucec_ag_cdc1_corr.csv (抗原 × cDC1 相关性)
"""
import os, numpy as np, pandas as pd
from scipy import stats

DATA = os.path.dirname(os.path.abspath(__file__))


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    sorted_p = pvals[order]
    q = sorted_p * n / (np.arange(n) + 1.0)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    out = np.empty(n)
    out[order] = q
    return out


# 4 组分
IFN_GENES = ["IFNG", "STAT1", "IRF1", "IRF7", "MX1", "OAS1", "ISG15", "IFI27", "IFIT1", "GBP1", "CXCL10"]
TIL_GENES = ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "GNLY", "NKG7", "CD3E"]
PDL1_GENES = ["CD274"]
TGFB_GENES = ["TGFB1", "TGFB2", "TGFB3", "ACVR1", "ACVR1B", "BMPR1A", "BMPR1B", "SMAD3", "SMAD4"]
CDC1_GENES = ["BATF3", "CLEC9A", "XCR1", "IRF8", "THBD"]


def main():
    m = pd.read_csv(os.path.join(DATA, "ucec_tpm_matrix.csv"), index_col=0)
    smap = pd.read_csv(os.path.join(DATA, "ucec_sample_map.csv"))
    tumor = [s for s in smap[smap["sample_type"] == "Primary Tumor"]["sample_id"] if s in m.columns]
    print(f"肿瘤样本: {len(tumor)}")

    log2 = np.log2(m[tumor].astype(float) + 1)

    # Z-score per gene (across samples)
    z = log2.subtract(log2.mean(axis=1), axis=0).divide(log2.std(axis=1) + 1e-8, axis=0)

    def score(genes):
        found = [g for g in genes if g in z.index]
        return z.loc[found].mean(axis=0), len(found)

    # 4 组分
    ifn_s, ifn_n = score(IFN_GENES)
    til_s, til_n = score(TIL_GENES)
    pdl1_s, pdl1_n = score(PDL1_GENES)
    tgfb_s, tgfb_n = score(TGFB_GENES)
    cdc1_s, cdc1_n = score(CDC1_GENES)

    print(f"可用基因: IFN={ifn_n}/{len(IFN_GENES)}, TIL={til_n}/{len(TIL_GENES)}, PD-L1={pdl1_n}/{len(PDL1_GENES)}, TGF-β={tgfb_n}/{len(TGFB_GENES)}, cDC1={cdc1_n}/{len(CDC1_GENES)}")

    time_index = (ifn_s + til_s + pdl1_s) / 3 - tgfb_s

    out = pd.DataFrame({
        "sample_id": tumor,
        "IFN": ifn_s.values,
        "TIL": til_s.values,
        "PDL1": pdl1_s.values,
        "TGFB": tgfb_s.values,
        "TIME_index": time_index.values,
    })
    out.to_csv(os.path.join(DATA, "ucec_time.csv"), index=False)
    print(f"\nTIME 统计:")
    print(out[["IFN", "TIL", "PDL1", "TGFB", "TIME_index"]].describe().round(3))

    # 按免疫亚型
    subtypes = pd.read_csv(os.path.join(DATA, "ucec_subtypes.csv"))
    out = out.merge(subtypes, on="sample_id")
    by_sub = out.groupby("subtype").agg({
        "IFN": "mean", "TIL": "mean", "PDL1": "mean",
        "TGFB": "mean", "TIME_index": ["mean", "median"],
    }).round(3)
    print(f"\nTIME × 免疫亚型:")
    print(by_sub)
    by_sub.to_csv(os.path.join(DATA, "ucec_time_by_subtype.csv"))

    # Kruskal-Wallis
    groups = [g["TIME_index"].values for _, g in out.groupby("subtype")]
    h, p = stats.kruskal(*groups)
    print(f"\nTIME_index Kruskal-Wallis: H={h:.2f}, P={p:.2e}")

    # cDC1
    cdc1_df = pd.DataFrame({"sample_id": tumor, "cDC1_score": cdc1_s.values})
    cdc1_df.to_csv(os.path.join(DATA, "ucec_cdc.csv"), index=False)
    print(f"\ncDC1 统计: mean={cdc1_s.mean():.3f}, std={cdc1_s.std():.3f}")

    cdc1_df = cdc1_df.merge(subtypes, on="sample_id")
    print(f"\ncDC1 × 免疫亚型:")
    print(cdc1_df.groupby("subtype")["cDC1_score"].agg(["mean", "median", "count"]).round(3))

    cdc1_df.to_csv(os.path.join(DATA, "ucec_cdc_by_subtype.csv"))

    # 抗原 × cDC1 相关性
    print(f"\n抗原 × cDC1 相关性:")
    # Need to wait for final antigens
    final_antigens_path = os.path.join(DATA, "ucec_final_antigens.csv")
    if os.path.exists(final_antigens_path):
        ag_df = pd.read_csv(final_antigens_path)
        antigens = ag_df["gene"].tolist()
    else:
        antigens = ["VTCN1", "MSLN", "ERBB2", "MUC16", "CLDN6", "TP53"]

    cdc1_vec = cdc1_s.values
    rows = []
    for ag in antigens:
        if ag not in m.index:
            print(f"  {ag}: not in expression matrix")
            continue
        ag_expr = log2.loc[ag].values
        r, p = stats.spearmanr(ag_expr, cdc1_vec)
        rows.append({"antigen": ag, "spearman_r": r, "pvalue": p, "n": len(ag_expr)})
    ag_cdc = pd.DataFrame(rows)
    if len(ag_cdc) > 0:
        ag_cdc["FDR"] = bh_fdr(ag_cdc["pvalue"].values)
    ag_cdc.to_csv(os.path.join(DATA, "ucec_ag_cdc1_corr.csv"), index=False)
    print(ag_cdc.round(4).to_string(index=False))

    print("\n=== TIME + cDC1 分析完成 ===")


if __name__ == "__main__":
    main()
