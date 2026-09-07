# -*- coding: utf-8 -*-
"""Fig 16: GSE120490 外部队列验证图（3 panel，COAD fig11 模式）。

(a) 免疫评分 + 关键检查点按投影亚型（热图）
(b) VTCN1/CD276 + CD8A/GZMB 按投影亚型（条形图）
(c) 分级 × 投影亚型（堆叠条形图）
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

DATA = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.abspath(os.path.join(DATA, "..", "figures"))
os.makedirs(FIGDIR, exist_ok=True)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 7
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['legend.frameon'] = False

PALETTE = {
    "blue_main": "#0F4D92", "red_strong": "#B64342", "teal": "#42949E",
    "neutral_mid": "#767676", "neutral_dark": "#4D4D4D", "neutral_black": "#272727",
}
SUBTYPE_COLORS = ["#0F4D92", "#8BCF8B", "#B64342", "#9A4D8E", "#42949E"]

proj = pd.read_csv(os.path.join(DATA, "gse120490_projection.csv"))
subs = ["IS1", "IS2", "IS3", "IS4", "IS5"]
proj["projected_subtype"] = pd.Categorical(proj["projected_subtype"], categories=subs, ordered=True)

fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.2))

# ---------- (a) 免疫评分 + 检查点热图 ----------
ax = axes[0]
genes_order = ["immune_score", "CD8A", "GZMB", "TIGIT", "CD274", "PDCD1", "LAG3", "CTLA4", "VTCN1", "CD276"]
genes_avail = [g for g in genes_order if g in proj.columns]
mat = proj.groupby("projected_subtype", observed=True)[genes_avail].mean().reindex(subs)
mat_z = (mat - mat.mean(axis=0)) / (mat.std(axis=0) + 1e-8)
im = ax.imshow(mat_z.T.values, cmap="RdBu_r", aspect="auto", vmin=-2, vmax=2)
ax.set_xticks(range(len(subs)))
ax.set_xticklabels(subs, fontsize=7)
ax.set_yticks(range(len(genes_avail)))
ax.set_yticklabels(genes_avail, fontsize=6.5)
# KW P 标注
for j, g in enumerate(genes_avail):
    groups = [x[g].dropna().values for _, x in proj.groupby("projected_subtype", observed=True) if len(x[g].dropna()) > 3]
    if len(groups) >= 2:
        _, p = stats.kruskal(*groups)
        if p < 0.05:
            ax.text(len(subs) - 0.3, j, "*", ha="left", va="center", fontsize=8, color="black")
ax.set_title("Immune features by projected subtype\n(GSE120490, n=145)", fontsize=8, fontweight="bold", pad=3)
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cbar.set_label("Z-score", fontsize=6.5)
cbar.ax.tick_params(labelsize=6)
ax.text(-0.15, 1.08, "a", transform=ax.transAxes, fontsize=9, fontweight='bold',
        ha='left', va='bottom')

# ---------- (b) 关键基因条形图（VTCN1/CD276 + CD8A/GZMB）----------
ax = axes[1]
sel = ["CD8A", "GZMB", "VTCN1", "CD276"]
means = proj.groupby("projected_subtype", observed=True)[sel].mean().reindex(subs)
x = np.arange(len(subs))
w = 0.2
colors4 = ["#0F4D92", "#42949E", "#B64342", "#9A4D8E"]
for k, g in enumerate(sel):
    ax.bar(x + (k - 1.5) * w, means[g], width=w, color=colors4[k], label=g)
ax.set_xticks(x)
ax.set_xticklabels(subs)
ax.set_ylabel("Mean log$_2$ expression")
ax.set_title("Key antigens & immune markers", fontsize=8, fontweight="bold", pad=3)
ax.legend(fontsize=6, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.13))
ax.text(-0.15, 1.08, "b", transform=ax.transAxes, fontsize=9, fontweight='bold',
        ha='left', va='bottom')

# ---------- (c) 分级 × 投影亚型 ----------
ax = axes[2]
ct = pd.crosstab(proj["projected_subtype"], proj["grade"]).reindex(subs)
ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100
bottom = np.zeros(len(subs))
for cat, color in [(c, col) for c, col in zip(ct.columns, ["#8BCF8B", "#B64342"])]:
    if cat in ct_pct.columns:
        vals = ct_pct[cat].values
        ax.bar(subs, vals, bottom=bottom, color=color, label=f"Grade {'high (Yes)' if cat == 'Yes' else 'low/no (No)'}", width=0.6)
        for i, (v, b) in enumerate(zip(vals, bottom)):
            ax.text(i, b + v / 2, f"{int(ct.values[i][list(ct.columns).index(cat)])}", ha="center", va="center", fontsize=6.5, color="white")
        bottom += vals
chi2, p, _, _ = stats.chi2_contingency(ct.values)
ax.set_ylabel("% of patients")
ax.set_title(f"Tumor grade × projected subtype\n$\\chi^2$={chi2:.1f}, P={p:.3f}", fontsize=8, fontweight="bold", pad=3, linespacing=1.3)
ax.legend(fontsize=6, loc="lower center", bbox_to_anchor=(0.5, -0.28), ncol=2)
ax.set_ylim(0, 100)
ax.text(-0.15, 1.08, "c", transform=ax.transAxes, fontsize=9, fontweight='bold',
        ha='left', va='bottom')

fig.tight_layout(pad=1.0)
fig.savefig(os.path.join(FIGDIR, "fig16_external_validation.svg"), bbox_inches="tight")
fig.savefig(os.path.join(FIGDIR, "fig16_external_validation.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(FIGDIR, "fig16_external_validation.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Saved fig16_external_validation to {FIGDIR}")
