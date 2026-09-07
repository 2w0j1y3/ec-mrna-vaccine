# -*- coding: utf-8 -*-
"""最终分子分型 + EC 的 cDC1 共变抗原鉴定。

1. 修正 MSI-H 分类：非 POLE 且 n_indel >= 80（indel/SNV 比验证）
2. 免疫亚型 × 分子分型（POLE / MSI-H / MSS）
3. 寻找 EC 中与 cDC1 共变的严格候选抗原（COAD CD276 故事的 EC 对应）
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

DATA = os.path.dirname(os.path.abspath(__file__))

# ---------- 1. 最终分子分型 ----------
df = pd.read_csv(os.path.join(DATA, "ucec_molecular_class.csv"))
df["is_pole"] = df["POLE_hotspot"].notna() & (df["POLE_hotspot"] != "")
df["indel_snv_ratio"] = df["n_indel"] / (df["n_snv"] + 1)

def classify(r):
    if r["is_pole"]:
        return "POLE"
    if r["n_indel"] >= 80 and r["indel_snv_ratio"] > 0.2:
        return "MSI-H"
    return "MSS"

df["molecular_class"] = df.apply(classify, axis=1)
print("最终分子分型分布:")
print(df["molecular_class"].value_counts())
print(f"\n各类别 TMB (SNV+indel / 38Mb):")
df["TMB_calc"] = (df["n_snv"] + df["n_indel"]) / 38.0
print(df.groupby("molecular_class")["TMB_calc"].agg(["median", "count"]).round(1))
df.to_csv(os.path.join(DATA, "ucec_molecular_class.csv"), index=False)

# ---------- 2. 免疫亚型 × 分子分型 ----------
subtype = pd.read_csv(os.path.join(DATA, "ucec_subtypes.csv"))
smap = pd.read_csv(os.path.join(DATA, "ucec_sample_map.csv"))
sample_to_case = dict(zip(smap["sample_id"], smap["case_id"]))
sub_cases = subtype.copy()
sub_cases["case_id"] = sub_cases["sample_id"].map(sample_to_case)
merged = sub_cases.merge(df[["case_id", "molecular_class", "TMB_calc"]], on="case_id", how="left")
merged["molecular_class"] = merged["molecular_class"].fillna("NA")

print(f"\n免疫亚型 × 分子分型:")
ct = pd.crosstab(merged["subtype"], merged["molecular_class"])
# 排序列
for col in ["POLE", "MSI-H", "MSS", "NA"]:
    if col not in ct.columns:
        ct[col] = 0
ct = ct[["POLE", "MSI-H", "MSS", "NA"]]
print(ct)
ct.to_csv(os.path.join(DATA, "ucec_subtype_vs_molecular.csv"))

chi2, p, dof, _ = stats.chi2_contingency(ct[["POLE", "MSI-H", "MSS"]].values)
print(f"\n卡方检验: chi2={chi2:.2f}, P={p:.2e}")

print(f"\n各亚型内超突变 (POLE+MSI-H) 比例:")
for s in sorted(merged["subtype"].unique()):
    sub = merged[merged["subtype"] == s]
    hot = (sub["molecular_class"].isin(["POLE", "MSI-H"])).sum()
    pole_n = (sub["molecular_class"] == "POLE").sum()
    msi_n = (sub["molecular_class"] == "MSI-H").sum()
    print(f"  {s}: POLE={pole_n}, MSI-H={msi_n}, 合计 {hot}/{len(sub)} ({hot/len(sub)*100:.1f}%)")

# ---------- 3. cDC1 共变抗原（EC 的 CD276 对应故事） ----------
print(f"\n" + "=" * 60)
print("3. 严格候选抗原中与 cDC1 显著共变的基因")
print("=" * 60)
final = pd.read_csv(os.path.join(DATA, "ucec_final_antigens.csv"))
ag_cdc = pd.read_csv(os.path.join(DATA, "ucec_ag_cdc1_corr.csv"))[["antigen", "spearman_r", "FDR"]]
ag_cdc = ag_cdc.rename(columns={"spearman_r": "cdc1_rho", "FDR": "cdc1_FDR"})
strict = final[(final["prognostic"]) & (final["apc_positive"])]
strict_cdc = strict.merge(ag_cdc, left_on="gene", right_on="antigen", how="left")
sig = strict_cdc[strict_cdc["cdc1_FDR"] < 0.05].sort_values("cdc1_rho", ascending=False)
print(f"严格候选中 cDC1 共变显著 (FDR<0.05): {len(sig)} 个")
print(sig[["gene", "log2FC", "mut_freq_pct", "km_logrank_p", "cdc1_rho", "cdc1_FDR"]].head(20).to_string(index=False))
sig.to_csv(os.path.join(DATA, "ucec_strict_cdc1_covarying.csv"), index=False)
