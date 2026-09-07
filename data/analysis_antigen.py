import os
# -*- coding: utf-8 -*-
"""EC (UCEC) 肿瘤抗原鉴定（适配 TCGA-UCEC）。

EC 特殊考虑：
1. 高 dMMR/MSI-H 比例（~30%），存在超突变子群
2. CN-high 亚型大量 TP53 突变
3. POLE 超突变子群 TMB 极高（>100 mut/Mb）
4. B7-H4 (VTCN1) 在 EC 高表达 —— 用户专长直接命中

输出：ucec_degs.csv, ucec_mut_freq.csv, ucec_antigen_candidates.csv, ucec_final_antigens.csv
"""
import numpy as np
import pandas as pd
from scipy import stats
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")


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


# ---------- 1. 载入数据 ----------
print("载入表达矩阵（样本级）...")
m = pd.read_csv(f"{DATA}\\ucec_tpm_matrix.csv", index_col=0)
smap = pd.read_csv(f"{DATA}\\ucec_sample_map.csv")

tumor = smap[smap["sample_type"] == "Primary Tumor"]["sample_id"].tolist()
normal = smap[smap["sample_type"] == "Solid Tissue Normal"]["sample_id"].tolist()
tumor = [s for s in tumor if s in m.columns]
normal = [s for s in normal if s in m.columns]
print(f"肿瘤 {len(tumor)} / 正常 {len(normal)}")

# 蛋白编码基因
pc_genes = set(open(f"{DATA}\\protein_coding_genes.txt").read().split("\n"))
m_pc = m[m.index.isin(pc_genes)]
print(f"蛋白编码基因: {m_pc.shape[0]}")

tumor_m = m_pc[tumor]
normal_m = m_pc[normal]

# ---------- 2. 差异表达 ----------
print("差异表达分析（Wilcoxon，蛋白编码基因）...")
tpm_t = tumor_m.values.astype(float)
tpm_n = normal_m.values.astype(float)
genes = tumor_m.index.tolist()

# EC 调整：normal 样本较少，使用更宽松的过滤
keep = (tpm_t.mean(axis=1) > 0.5) | (tpm_n.mean(axis=1) > 0.1)
genes = [g for g, k in zip(genes, keep) if k]
tpm_t = tpm_t[keep]
tpm_n = tpm_n[keep]
print(f"过滤后基因数: {len(genes)}")

log2fc = np.empty(len(genes))
pvals = np.empty(len(genes))
for i in range(len(genes)):
    log2fc[i] = np.log2(tpm_t[i].mean() + 1) - np.log2(tpm_n[i].mean() + 1)
    _, p = stats.mannwhitneyu(np.log2(tpm_t[i] + 1), np.log2(tpm_n[i] + 1))
    pvals[i] = p

fdr = bh_fdr(pvals)
deg = pd.DataFrame({
    "gene": genes, "log2FC": log2fc,
    "mean_tumor_tpm": tpm_t.mean(axis=1),
    "mean_normal_tpm": tpm_n.mean(axis=1),
    "pvalue": pvals, "FDR": fdr,
})
deg = deg.sort_values("FDR")
over = deg[(deg["log2FC"] > 1) & (deg["FDR"] < 0.05)]
under = deg[(deg["log2FC"] < -1) & (deg["FDR"] < 0.05)]
print(f"过表达基因: {len(over)}，下调基因: {len(under)}")
deg.to_csv(f"{DATA}\\ucec_degs.csv", index=False)

# ---------- 3. 突变频次 ----------
print("载入突变数据 ...")
mut_freq = pd.read_csv(f"{DATA}\\ucec_mut_freq.csv")
mut_freq = mut_freq[mut_freq["gene"].isin(pc_genes)]
# EC 有大量超突变患者，调整阈值为 1% 仍合适
mutated_genes = set(mut_freq[mut_freq["mut_freq_pct"] >= 1.0]["gene"])
print(f"突变频率 >=1% 的蛋白编码基因: {len(mutated_genes)}")

# ---------- 4. 抗原候选 = 过表达 ∩ 突变 ----------
over_genes = set(over["gene"])
candidates = over_genes & mutated_genes
print(f"过表达 ∩ 突变 候选基因: {len(candidates)}")

cand_df = deg[deg["gene"].isin(candidates)].copy()
cand_df = cand_df.merge(mut_freq, on="gene", how="left")
cand_df = cand_df.sort_values("log2FC", ascending=False)
cand_df.to_csv(f"{DATA}\\ucec_antigen_candidates.csv", index=False)
print("候选基因（前 30）:")
print(cand_df[["gene", "log2FC", "FDR", "mut_freq_pct"]].head(30).to_string())

# ---------- 5. 预后筛选（KM 生存） ----------
print("\n预后筛选（KM log-rank）...")
clin = pd.read_csv(f"{DATA}\\ucec_clinical.csv")
clin = clin.dropna(subset=["vital_status"])
clin = clin[clin["vital_status"].isin(["Alive", "Dead"])]
clin["event"] = (clin["vital_status"] == "Dead").astype(int)
clin["time"] = np.where(clin["event"] == 1, clin["days_to_death"], clin["days_to_last_follow_up"])
clin = clin[(clin["time"] > 0)]
clin["time_mon"] = clin["time"] / 30.44  # 转为月
surv = clin.set_index("submitter_id")[["time_mon", "event"]]
print(f"可用于生存分析病例: {len(surv)}")

# 肿瘤样本 -> 病例
sample_to_case = dict(zip(smap["sample_id"], smap["case_id"]))

def km_logrank(gene):
    """按基因表达中位数分组，log-rank 检验。返回 p 值。"""
    exp = tumor_m.loc[gene].astype(float)
    med = exp.median()
    high_samples = exp[exp > med].index.tolist()
    low_samples = exp[exp <= med].index.tolist()
    high_cases = [sample_to_case[s] for s in high_samples if s in sample_to_case]
    low_cases = [sample_to_case[s] for s in low_samples if s in sample_to_case]
    high_cases = [c for c in high_cases if c in surv.index]
    low_cases = [c for c in low_cases if c in surv.index]
    if len(high_cases) < 5 or len(low_cases) < 5:
        return np.nan
    t_high = surv.loc[high_cases, "time_mon"].values
    e_high = surv.loc[high_cases, "event"].values
    t_low = surv.loc[low_cases, "time_mon"].values
    e_low = surv.loc[low_cases, "event"].values
    res = logrank_test(t_high, t_low, event_observed_A=e_high, event_observed_B=e_low)
    return res.p_value

km_p = {}
for g in candidates:
    km_p[g] = km_logrank(g)

# ---------- 6. APC 浸润关联（Spearman） ----------
print("APC 浸润关联分析 ...")
# 抗原呈递细胞标记基因（DC + 巨噬细胞 + B 细胞）
APC_MARKERS = [
    # 树突状细胞
    "CD1C", "CLEC9A", "CLEC10A", "ITGAX", "ITGAE", "BATF3", "LAMP3", "FLT3", "THBD", "NRP1",
    # 巨噬细胞
    "CD68", "CD163", "MRC1", "CD14", "ITGAM", "CSF1R",
    # B 细胞
    "CD19", "CD79A", "CD79B", "MS4A1",
]
apc_markers_found = [g for g in APC_MARKERS if g in tumor_m.index]
print(f"APC 标记基因可用: {len(apc_markers_found)}/{len(APC_MARKERS)}")
apc_score = np.log2(tumor_m.loc[apc_markers_found].astype(float) + 1).mean(axis=0)

def spearman_with_apc(gene):
    exp = tumor_m.loc[gene].astype(float)
    r, p = stats.spearmanr(exp.values, apc_score.values)
    return r, p

apc_corr = {}
for g in candidates:
    apc_corr[g] = spearman_with_apc(g)

# ---------- 汇总最终抗原 ----------
print("\n汇总最终抗原 ...")
final = []
for g in sorted(candidates, key=lambda x: -deg.loc[deg["gene"] == x, "log2FC"].values[0]):
    row = deg[deg["gene"] == g].iloc[0]
    mf = mut_freq[mut_freq["gene"] == g]["mut_freq_pct"].values[0]
    p_km = km_p.get(g, np.nan)
    r_apc, p_apc = apc_corr.get(g, (np.nan, np.nan))
    final.append({
        "gene": g, "log2FC": row["log2FC"], "FDR": row["FDR"],
        "mean_tumor_tpm": row["mean_tumor_tpm"], "mut_freq_pct": mf,
        "km_logrank_p": p_km, "apc_spearman_r": r_apc, "apc_spearman_p": p_apc,
    })

final_df = pd.DataFrame(final)
# EC 调整：考虑到 EC 异质性，亚型间 KM 差异大
final_df["prognostic"] = final_df["km_logrank_p"] < 0.05
final_df["apc_positive"] = final_df["apc_spearman_r"] > 0
final_df.to_csv(f"{DATA}\\ucec_final_antigens.csv", index=False)

print(f"\n最终抗原候选: {len(final_df)} 个")
print(f"其中预后显著: {final_df['prognostic'].sum()}")
print(f"其中 APC 正相关: {final_df['apc_positive'].sum()}")
print(final_df.to_string(index=False))

# ---------- EC 特异：标注 B7-H4 ----------
print("\n=== EC 重点关注：B7-H4 (VTCN1) ===")
vtcn1_row = deg[deg["gene"] == "VTCN1"]
if len(vtcn1_row) > 0:
    print(f"VTCN1: log2FC={vtcn1_row['log2FC'].values[0]:.2f}, FDR={vtcn1_row['FDR'].values[0]:.2e}")
    print(f"  Tumor TPM: {vtcn1_row['mean_tumor_tpm'].values[0]:.2f}")
    print(f"  Normal TPM: {vtcn1_row['mean_normal_tpm'].values[0]:.2f}")
    print(f"  KM logrank P: {km_p.get('VTCN1', np.nan):.2e}")
    r_apc, p_apc = apc_corr.get('VTCN1', (np.nan, np.nan))
    print(f"  APC Spearman rho={r_apc:.3f}, P={p_apc:.2e}")

print("\n=== 抗原鉴定完成 ===")
