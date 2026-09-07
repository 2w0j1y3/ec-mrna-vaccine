# -*- coding: utf-8 -*-
"""CPTAC UCEC 验证分析：
1. 6 抗原 mRNA-蛋白一致性（Spearman）
2. VTCN1/B7-H4 蛋白丰度 vs 免疫浸润（CIBERSORT）+ 亚型相关
3. GEP/TIME/cDC1 相关架构 vs TCGA（用 mRNA 算 GEP/TIME/cDC1，与 CIBERSORT/ESTIMATE 交叉验证）
4. 6 抗原 mRNA vs TMB/MSI/POLE 亚型
输出: cptac_validation.csv
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

DATA = os.path.dirname(os.path.abspath(__file__))
GENES = ["CLDN6", "CCNE1", "MAL", "CTSV", "VTCN1", "MUC16"]
IMMUNE = ["CD274", "PDCD1", "CTLA4", "TIGIT", "LAG3", "HAVCR2", "CD276",
          "CD8A", "GZMB", "IFNG", "STAT1", "CXCL9", "FOXP3",
          "BATF3", "CLEC9A", "XCR1", "IRF8", "TGFB1", "CD68"]

# entrez → symbol 映射（重新解析）
import requests
HDRS = {"Accept": "application/json"}
entrez_to_sym = {}
for g in GENES + IMMUNE:
    try:
        r = requests.get(f"https://www.cbioportal.org/api/genes/{g}", headers=HDRS, timeout=30)
        if r.status_code == 200:
            entrez_to_sym[r.json()["entrezGeneId"]] = g
    except Exception:
        pass
print(f"Resolved {len(entrez_to_sym)} entrez → symbol mappings")

# ---------- 载入 ----------
mrna = pd.read_csv(os.path.join(DATA, "cptac_mrna.csv"))
prot = pd.read_csv(os.path.join(DATA, "cptac_protein.csv"))
clin = pd.read_csv(os.path.join(DATA, "cptac_clinical.csv"))
mrna["gene"] = mrna["entrezGeneId"].map(entrez_to_sym)
prot["gene"] = prot["entrezGeneId"].map(entrez_to_sym)

m_piv = mrna.pivot_table(index="gene", columns="sampleId", values="value")
p_piv = prot.pivot_table(index="gene", columns="sampleId", values="value")
print(f"mRNA matrix: {m_piv.shape}, protein matrix: {p_piv.shape}")
print(f"Protein genes available: {sorted(p_piv.index.dropna().tolist())}")

# ---------- 1. mRNA-蛋白一致性 ----------
print("\n" + "=" * 60)
print("1. 6 抗原 mRNA-蛋白一致性")
print("=" * 60)
results = []
for g in GENES:
    if g in m_piv.index and g in p_piv.index:
        common = [c for c in m_piv.columns if c in p_piv.columns]
        m_v = m_piv.loc[g, common].astype(float)
        p_v = p_piv.loc[g, common].astype(float)
        ok = ~(m_v.isna() | p_v.isna())
        m_v, p_v = m_v[ok], p_v[ok]
        if ok.sum() > 10:
            r_, p_ = stats.spearmanr(m_v, p_v)
            results.append({"gene": g, "rho_mrna_protein": r_, "p": p_, "n": ok.sum()})
            print(f"  {g:8s}: rho={r_:+.3f}, P={p_:.2e} (n={ok.sum()})")
        else:
            print(f"  {g:8s}: n={ok.sum()} 不足")
    else:
        print(f"  {g:8s}: mRNA={'Y' if g in m_piv.index else 'N'}, protein={'Y' if g in p_piv.index else 'N'}")

# ---------- 2. VTCN1 蛋白 vs 免疫浸润（CIBERSORT/ESTIMATE） ----------
print("\n" + "=" * 60)
print("2. VTCN1 蛋白/抗原 vs CPTAC 免疫指标")
print("=" * 60)
clin_idx = clin.drop_duplicates(subset="sampleId").set_index("sampleId")
imm_cols = {
    "CD8 T (CIBERSORT)": "CIBERSORT_T _CELLS _CD8",
    "Treg (CIBERSORT)": "CIBERSORT_T _CELLS _REGULATORY _(TREGS)",
    "Macrophage (CIBERSORT)": "CIBERSORT_MACROPHAGES _M0",
    "M1 (CIBERSORT)": "CIBERSORT_MACROPHAGES _M1",
    "M2 (CIBERSORT)": "CIBERSORT_MACROPHAGES _M2",
    "ImmuneScore (ESTIMATE)": "ESTIMATE_IMMUNESCORE",
    "TMB": "TMB_NONSYNONYMOUS",
}
for label, col in imm_cols.items():
    if col not in clin_idx.columns:
        print(f"  [skip] {label}: col not found")
        continue
    vals = pd.to_numeric(clin_idx[col], errors="coerce")
    for g in ["VTCN1", "CLDN6", "CD276"]:
        src = p_piv if g in p_piv.index else m_piv
        if g not in src.index:
            continue
        common = [c for c in src.columns if c in vals.index]
        gv = src.loc[g, common].astype(float)
        iv = vals.reindex(common)
        ok = ~(gv.isna() | iv.isna())
        if ok.sum() > 10:
            r_, p_ = stats.spearmanr(gv[ok], iv[ok])
            layer = "protein" if src is p_piv else "mRNA"
            star = " *" if p_ < 0.05 else ""
            print(f"  {g:8s}({layer}) vs {label:24s}: rho={r_:+.3f}, P={p_:.2e}{star}")

# ---------- 3. GEP / TIME / cDC1 计算与验证 ----------
print("\n" + "=" * 60)
print("3. GEP/TIME/cDC1 架构在 CPTAC 的计算")
print("=" * 60)
GEP_GENES = ["IFNG", "STAT1", "CCR5", "CXCL9", "CXCL10", "CXCL11", "IDO1",
             "PRF1", "GZMA", "GZMB", "CD8A", "CD8B", "PDCD1", "PDCD1LG2",
             "CD274", "LAG3", "HAVCR2", "TIGIT"]
CDC1_GENES = ["BATF3", "CLEC9A", "XCR1", "IRF8", "THBD"]

def gene_score(piv, genes):
    found = [g for g in genes if g in piv.index]
    print(f"    available {len(found)}/{len(genes)}")
    if len(found) < 2:
        return None
    z = piv.loc[found].astype(float)
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1) + 1e-8, axis=0)
    return z.mean(axis=0)

gep_c = gene_score(m_piv, GEP_GENES)
cdc1_c = gene_score(m_piv, CDC1_GENES)

if gep_c is not None and cdc1_c is not None:
    common = [c for c in gep_c.index if c in cdc1_c.index]
    r_, p_ = stats.spearmanr(gep_c[common], cdc1_c[common])
    print(f"  GEP vs cDC1: rho={r_:+.3f}, P={p_:.2e} (n={len(common)})")

# GEP vs CIBERSORT CD8 / ESTIMATE immune score
if gep_c is not None:
    for label, col in [("CD8 (CIBERSORT)", "CIBERSORT_T _CELLS _CD8"),
                       ("ImmuneScore (ESTIMATE)", "ESTIMATE_IMMUNESCORE"),
                       ("TMB", "TMB_NONSYNONYMOUS")]:
        if col not in clin_idx.columns:
            continue
        vals = pd.to_numeric(clin_idx[col], errors="coerce")
        common = [c for c in gep_c.index if c in vals.index]
        g_v = gep_c[common]
        i_v = vals.reindex(common)
        ok = ~(g_v.isna() | i_v.isna())
        if ok.sum() > 10:
            r_, p_ = stats.spearmanr(g_v[ok], i_v[ok])
            print(f"  GEP vs {label:24s}: rho={r_:+.3f}, P={p_:.2e}")

# GEP by POLE / MSI 亚型
if gep_c is not None:
    for subcol, name in [("POLE_SUBTYPE", "POLE 亚型"), ("MSI_STATUS", "MSI 状态")]:
        if subcol not in clin_idx.columns:
            continue
        print(f"\n  GEP by {name}:")
        groups = []
        subvals = clin_idx[subcol].reindex(gep_c.index)
        for grp in subvals.dropna().unique():
            g = gep_c[subvals == grp].dropna()
            if len(g) >= 3:
                groups.append(g.values)
                print(f"    {grp}: n={len(g)}, mean GEP z={g.mean():+.3f}")
        if len(groups) >= 2:
            h, p_ = stats.kruskal(*groups)
            print(f"    Kruskal-Wallis: H={h:.2f}, P={p_:.2e}")

# ---------- 4. 6 抗原 mRNA vs TMB/MSI/POLE ----------
print("\n" + "=" * 60)
print("4. 6 抗原 mRNA vs 分子亚型 (POLE/MSI)")
print("=" * 60)
for subcol in ["POLE_SUBTYPE", "MSI_STATUS"]:
    if subcol not in clin_idx.columns:
        continue
    subvals = clin_idx[subcol].reindex(m_piv.columns)
    print(f"\n  by {subcol}:")
    for g in GENES:
        if g not in m_piv.index:
            continue
        gv = m_piv.loc[g].astype(float)
        groups = []
        for grp in subvals.dropna().unique():
            v = gv[subvals == grp].dropna()
            if len(v) >= 3:
                groups.append(v.values)
        if len(groups) >= 2:
            h, p_ = stats.kruskal(*groups)
            if p_ < 0.1:
                print(f"    {g:8s}: H={h:.2f}, P={p_:.2e}")

# 保存汇总
pd.DataFrame(results).to_csv(os.path.join(DATA, "cptac_validation.csv"), index=False)
print("\nSaved cptac_validation.csv")
