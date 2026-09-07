# -*- coding: utf-8 -*-
"""诊断 EC 分析中"不理想"的部分：
1. MSLN/TP53 为何没进候选（看它们的真实表达数据）
2. CLDN6 的 KM 方向（高表达是危险还是保护？）
3. VTCN1 (B7-H4) 在亚型/分级内的 KM（整体不显著可能是被稀释）
4. 各亚型的中位生存（确认 IS5 是否最差）
"""
import os
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

DATA = os.path.dirname(os.path.abspath(__file__))

m = pd.read_csv(os.path.join(DATA, "ucec_tpm_matrix.csv"), index_col=0)
smap = pd.read_csv(os.path.join(DATA, "ucec_sample_map.csv"))
clin = pd.read_csv(os.path.join(DATA, "ucec_clinical.csv"))
subtype = pd.read_csv(os.path.join(DATA, "ucec_subtypes.csv"))
deg = pd.read_csv(os.path.join(DATA, "ucec_degs.csv"))

tumor = [s for s in smap[smap["sample_type"] == "Primary Tumor"]["sample_id"] if s in m.columns]
normal = [s for s in smap[smap["sample_type"] == "Solid Tissue Normal"]["sample_id"] if s in m.columns]
sample_to_case = dict(zip(smap["sample_id"], smap["case_id"]))

# 生存
clin = clin[clin["vital_status"].isin(["Alive", "Dead"])].dropna(subset=["vital_status"])
clin["event"] = (clin["vital_status"] == "Dead").astype(int)
clin["time"] = np.where(clin["event"] == 1, clin["days_to_death"], clin["days_to_last_follow_up"])
clin = clin[clin["time"] > 0]
clin["time_mon"] = clin["time"] / 30.44
surv = clin.set_index("submitter_id")[["time_mon", "event"]]

print("=" * 70)
print("【1】6 个提名抗原的真实表达数据（不管是否过了过滤）")
print("=" * 70)
targets = ["VTCN1", "MSLN", "ERBB2", "MUC16", "CLDN6", "TP53"]
for g in targets:
    if g in m.index:
        t_mean = m.loc[g, tumor].astype(float).mean()
        n_mean = m.loc[g, normal].astype(float).mean() if normal else np.nan
        l2fc = np.log2(t_mean + 1) - np.log2(n_mean + 1)
        # 检测率
        det_t = (m.loc[g, tumor].astype(float) > 1).mean() * 100
        det_n = (m.loc[g, normal].astype(float) > 1).mean() * 100 if normal else np.nan
        print(f"  {g:8s}: tumor TPM={t_mean:8.2f} (检出{det_t:5.1f}%), normal TPM={n_mean:8.2f} (检出{det_n:5.1f}%), 手算log2FC={l2fc:+.2f}")
    else:
        print(f"  {g:8s}: 不在矩阵中")

print()
print("=" * 70)
print("【2】各亚型中位 OS（确认预后排序）")
print("=" * 70)
sub_surv = subtype.merge(pd.DataFrame({"sample_id": tumor}), on="sample_id")
sub_surv["case_id"] = sub_surv["sample_id"].map(sample_to_case)
sub_surv = sub_surv[sub_surv["case_id"].isin(surv.index)]
for s in sorted(sub_surv["subtype"].unique()):
    cases = sub_surv[sub_surv["subtype"] == s]["case_id"]
    t = surv.loc[cases, "time_mon"]
    e = surv.loc[cases, "event"]
    km = KaplanMeierFitter().fit(t, e)
    med = km.median_survival_time_
    print(f"  {s}: n={len(cases)}, 死亡{int(e.sum())}例, 中位OS = {med:.1f} 月" if not np.isnan(med) else f"  {s}: n={len(cases)}, 中位OS 未达到")

print()
print("=" * 70)
print("【3】CLDN6 KM 方向（高危还是保护）")
print("=" * 70)
def km_direction(gene, cases_mask=None, label=""):
    exp = m.loc[gene, tumor].astype(float)
    med = exp.median()
    hi_cases = [sample_to_case[s] for s in exp[exp > med].index if sample_to_case[s] in surv.index]
    lo_cases = [sample_to_case[s] for s in exp[exp <= med].index if sample_to_case[s] in surv.index]
    if cases_mask is not None:
        hi_cases = [c for c in hi_cases if c in cases_mask]
        lo_cases = [c for c in lo_cases if c in cases_mask]
    if len(hi_cases) < 10 or len(lo_cases) < 10:
        print(f"  {gene}{label}: 样本不足 (hi={len(hi_cases)}, lo={len(lo_cases)})")
        return
    t_h = surv.loc[hi_cases, "time_mon"].values; e_h = surv.loc[hi_cases, "event"].values
    t_l = surv.loc[lo_cases, "time_mon"].values; e_l = surv.loc[lo_cases, "event"].values
    res = logrank_test(t_h, t_l, event_observed_A=e_h, event_observed_B=e_l)
    km_h = KaplanMeierFitter().fit(t_h, e_h).median_survival_time_
    km_l = KaplanMeierFitter().fit(t_l, e_l).median_survival_time_
    med_h = f"{km_h:.1f}" if not np.isnan(km_h) else "NR"
    med_l = f"{km_l:.1f}" if not np.isnan(km_l) else "NR"
    direction = "高危(高表达预后差)" if (np.isnan(km_h) and not np.isnan(km_l)) or (not np.isnan(km_h) and not np.isnan(km_l) and km_h < km_l) else "保护(高表达预后好)"
    print(f"  {gene}{label}: P={res.p_value:.2e}, 高表达中位OS={med_h}月 vs 低表达={med_l}月 → {direction}")

km_direction("CLDN6")

print()
print("=" * 70)
print("【4】VTCN1 (B7-H4) 分层 KM（整体不显著 → 看亚型内）")
print("=" * 70)
km_direction("VTCN1", label=" [全队列]")
sub_case = sub_surv.set_index("case_id")["subtype"]
for s in ["IS5", "IS3", "IS2"]:
    mask = set(sub_case[sub_case == s].index)
    km_direction("VTCN1", cases_mask=mask, label=f" [{s}内]")

# G3 内
clin_g = clin.set_index("submitter_id")
g3_cases = set(clin_g[clin_g["tumor_grade"] == "G3"].index) & set(surv.index)
km_direction("VTCN1", cases_mask=g3_cases, label=" [G3内]")

print()
print("=" * 70)
print("【5】MSLN / TP53 细看：是不是 threshold 问题")
print("=" * 70)
mf = pd.read_csv(os.path.join(DATA, "ucec_mut_freq.csv"))
for g in ["MSLN", "TP53", "MUC16", "VTCN1", "CLDN6", "ERBB2"]:
    row = mf[mf["gene"] == g]
    mut = row["mut_freq_pct"].values[0] if len(row) else 0
    drow = deg[deg["gene"] == g]
    if len(drow):
        print(f"  {g:8s}: mut_freq={mut:5.2f}%, log2FC={drow['log2FC'].values[0]:+.2f}, FDR={drow['FDR'].values[0]:.2e}")
    else:
        # 查原始 DEG
        print(f"  {g:8s}: mut_freq={mut:5.2f}%, 不在DEG表中（表达量低于过滤线）")
