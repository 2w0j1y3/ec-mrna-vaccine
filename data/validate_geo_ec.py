# -*- coding: utf-8 -*-
"""GSE120490 最近中心投影验证 TCGA-UCEC 免疫亚型（复用 COAD GSE39582 验证模式）。

验证点:
1. 免疫亚型分布与投影置信度
2. 免疫架构复现（immune score / checkpoints 按投影亚型）
3. VTCN1 (B7-H4) 是否 IS5 最高
4. 6 抗原表达按投影亚型
5. 转移/分级/肌层浸润 × 投影亚型关联

输出: gse120490_projection.csv, gse120490_validation_stats.txt
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

DATA = os.path.dirname(os.path.abspath(__file__))

# ---------- 1. TCGA 免疫基因中心 ----------
print("Loading TCGA reference ...")
m = pd.read_csv(os.path.join(DATA, "ucec_tpm_matrix.csv"), index_col=0)
smap = pd.read_csv(os.path.join(DATA, "ucec_sample_map.csv"))
subtype = pd.read_csv(os.path.join(DATA, "ucec_subtypes.csv"))
tumor = [s for s in smap[smap["sample_type"] == "Primary Tumor"]["sample_id"] if s in m.columns]
sub_map = dict(zip(subtype["sample_id"], subtype["subtype"]))
tumor_sub = [s for s in tumor if s in sub_map]

immune_genes_df = pd.read_csv(os.path.join(DATA, "ucec_immune_genes.csv"))
immune_genes = immune_genes_df["gene"].tolist()
immune_found = [g for g in immune_genes if g in m.index]
print(f"TCGA immune genes: {len(immune_found)}")

# TCGA z-score（用于中心）
tcga_log = np.log2(m.loc[immune_found, tumor_sub].astype(float) + 1)
tcga_z = tcga_log.sub(tcga_log.mean(axis=1), axis=0).div(tcga_log.std(axis=1) + 1e-8, axis=0)
tcga_labels = pd.Series([sub_map[s] for s in tumor_sub], index=tcga_z.columns)
centroids = tcga_z.T.groupby(tcga_labels).mean().T  # gene x subtype
print(f"Centroids: {centroids.shape}")

# ---------- 2. GSE120490 投影 ----------
print("\nProjecting GSE120490 ...")
geo = pd.read_csv(os.path.join(DATA, "gse120490_expr.csv"), index_col=0)
gclin = pd.read_csv(os.path.join(DATA, "gse120490_clinical.csv"))
print(f"GEO matrix: {geo.shape}")

geo_genes = [g for g in immune_found if g in geo.index]
print(f"Common immune genes: {len(geo_genes)}/{len(immune_found)}")

# GEO 已是 log2 尺度（Affy MAS5/RMA log2）；再 z-score
geo_sub = geo.loc[geo_genes].astype(float)
geo_z = geo_sub.sub(geo_sub.mean(axis=1), axis=0).div(geo_sub.std(axis=1) + 1e-8, axis=0)

# Pearson 相关投影
cent_sub = centroids.loc[geo_genes]
corr = pd.DataFrame(index=geo_z.columns, columns=cent_sub.columns, dtype=float)
for s in geo_z.columns:
    v = geo_z[s].values
    for c in cent_sub.columns:
        u = cent_sub[c].values
        ok = ~(np.isnan(v) | np.isnan(u))
        if ok.sum() > 10:
            corr.loc[s, c] = np.corrcoef(v[ok], u[ok])[0, 1]

assigned = corr.idxmax(axis=1)
max_corr = corr.max(axis=1)
proj = pd.DataFrame({
    "sample_id": geo_z.columns,
    "projected_subtype": assigned,
    "max_corr": max_corr,
})
print(f"\n投影分布:")
print(proj["projected_subtype"].value_counts().sort_index())
print(f"\nmax_corr: median={max_corr.median():.3f}, "
      f"IQR=({max_corr.quantile(0.25):.3f}-{max_corr.quantile(0.75):.3f})")

# ---------- 3. 验证 ----------
report = []
def log(s=""):
    print(s)
    report.append(s)

log("\n" + "=" * 60)
log("3. 免疫架构复现（投影亚型）")
log("=" * 60)
# immune score
imm_score = geo_z.mean(axis=0)
proj["immune_score"] = proj["sample_id"].map(imm_score)
groups = [g["immune_score"].values for _, g in proj.groupby("projected_subtype")]
h, p = stats.kruskal(*groups)
log(f"immune score by projected subtype: KW H={h:.2f}, P={p:.2e}")
log(proj.groupby("projected_subtype")["immune_score"].mean().round(3).to_string())

# checkpoints + VTCN1 + 6 抗原（用原始 log2 表达）
log("\n关键基因按投影亚型（均值 log2）:")
key_genes = ["VTCN1", "CD274", "PDCD1", "CTLA4", "TIGIT", "LAG3",
             "CLDN6", "CCNE1", "MAL", "CTSV", "MUC16", "CD276", "CD8A", "GZMB"]
expr_map = {}
for g in key_genes:
    if g in geo.index:
        vals = geo.loc[g].astype(float)
        proj[g] = proj["sample_id"].map(vals)
        expr_map[g] = True
    else:
        expr_map[g] = False
        log(f"  [missing] {g} not in GEO matrix")

for g in key_genes:
    if not expr_map.get(g):
        continue
    groups = [x[g].dropna().values for _, x in proj.groupby("projected_subtype") if len(x[g].dropna()) > 3]
    if len(groups) >= 2:
        h, p = stats.kruskal(*groups)
        means = proj.groupby("projected_subtype")[g].mean().round(2).to_dict()
        star = " *" if p < 0.05 else ""
        log(f"  {g:8s}: KW P={p:.2e}{star} | " + " ".join(f"{k}={v}" for k, v in sorted(means.items())))

# ---------- 4. 转移/分级/肌层浸润关联 ----------
log("\n" + "=" * 60)
log("4. 临床终点 × 投影亚型")
log("=" * 60)
gclin_idx = gclin.set_index("sample_id")
for field in ["matastasis", "grade", "mi50pluspercent"]:
    vals = gclin_idx[field].reindex(proj["sample_id"]).values
    proj[field] = vals
    ct = pd.crosstab(proj["projected_subtype"], proj[field])
    log(f"\n{field} × projected subtype:")
    log(ct.to_string())
    try:
        chi2, p, _, _ = stats.chi2_contingency(ct.values)
        log(f"chi2={chi2:.2f}, P={p:.3f}")
    except Exception as e:
        log(f"chi2 failed: {e}")

# ---------- 保存 ----------
proj.to_csv(os.path.join(DATA, "gse120490_projection.csv"), index=False)
corr.to_csv(os.path.join(DATA, "gse120490_corr.csv"))
with open(os.path.join(DATA, "gse120490_validation_stats.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(report))
print("\nSaved gse120490_projection.csv + gse120490_validation_stats.txt")
