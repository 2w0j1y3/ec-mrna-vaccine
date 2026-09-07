"""Deeper single-cell validation: compare cell-type-specific expression of
IS3 vs IS5 signature genes, VTCN1 (B7-H4), and key checkpoints across
HEEC1 (high-grade), MEEC1 (middle), LEEC1 (low) tumors.

Key validations:
1. IS3 signature should be elevated in CD8 T cells and macrophages of hot tumors
2. IS5 signature should be elevated in malignant epithelium of cold/IS5 tumors
3. VTCN1 (B7-H4) should be expressed primarily in malignant epithelium
4. CD8A / GZMK should be enriched in CD8_T cells, CD68 in macrophages
5. Per-sample: IS3-IS5 score should correlate with T cell fraction

Output:
  data/sc_ucec/cell_level_validation.csv (per-cell IS3/IS5/VTCN1 by cell type)
  figures/fig17_sc_validation.png (3 panels)
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd

SC_DIR = Path(__file__).parent / "sc_ucec"
FIG_DIR = Path(__file__).parent.parent / "figures"
os.makedirs(FIG_DIR, exist_ok=True)

print("Loading cell assignments ...")
cells = pd.read_csv(SC_DIR / "cell_assignments.csv")
print(f"  {len(cells)} cells total")

# Add grade info per sample
SAMPLE_GRADE = {"HEEC1": "High", "HEEC2": "High", "HEEC3": "High",
                "MEEC1": "Middle", "MEEC2": "Middle",
                "LEEC1": "Low", "LEEC2": "Low"}
cells["grade"] = cells["sample"].map(SAMPLE_GRADE)
print(cells["grade"].value_counts())

# ====================== Per-cell-type expression analysis ======================
print("\n=== Per-cell-type mean gene expression ===")
key_genes = ["VTCN1", "CD274", "CD8A", "GZMK", "CD68", "MS4A1", "EPCAM", "KRT8", "FOXP3"]
cell_types_main = ["Epithelial", "CD8_T", "CD4_T", "Treg", "B_cell", "Plasma",
                  "Macrophage", "Monocyte", "DC", "NK", "Fibroblast", "Endothelial"]

ct_expr = (cells[~cells["cell_type"].isin(["Unknown"])]
           .groupby("cell_type")[[f"expr_{g}" for g in key_genes]].mean())
ct_expr.columns = [c.replace("expr_", "") for c in ct_expr.columns]
print(ct_expr.round(2).to_string())

# Save
ct_expr.to_csv(SC_DIR / "cell_type_expression.csv")

# ====================== IS3/IS5 score per cell type ======================
print("\n=== IS3/IS5 score distribution per cell type ===")
ct_score = cells.groupby("cell_type")[["IS3_score", "IS5_score"]].agg(["mean", "std"])
print(ct_score.round(3).to_string())

# ====================== Per-tumor summary ======================
print("\n=== Per-tumor metrics ===")
per_tumor = (cells.groupby(["sample", "grade"])
             .agg(n_cells=("barcode", "count"),
                  IS3_mean=("IS3_score", "mean"),
                  IS5_mean=("IS5_score", "mean"),
                  VTCN1_mean=("expr_VTCN1", "mean"),
                  CD8A_mean=("expr_CD8A", "mean"),
                  CD68_mean=("expr_CD68", "mean"))
             .reset_index())
per_tumor["IS3_minus_IS5"] = per_tumor["IS3_mean"] - per_tumor["IS5_mean"]
per_tumor["frac_IS3_high"] = (cells.groupby(["sample"])["IS3_score"]
                               .apply(lambda x: (x > 0).mean()).values)
per_tumor["frac_IS5_high"] = (cells.groupby(["sample"])["IS5_score"]
                               .apply(lambda x: (x > 0).mean()).values)
print(per_tumor.to_string(index=False))
per_tumor.to_csv(SC_DIR / "per_tumor_metrics.csv", index=False)

print("\n=== Per-tumor fraction of IS3-high and IS5-high cells ===")
for _, row in per_tumor.iterrows():
    print(f"  {row['sample']:8s} ({row['grade']:6s}): IS3+ {row['frac_IS3_high']:.1%}, "
          f"IS5+ {row['frac_IS5_high']:.1%}, VTCN1 mean={row['VTCN1_mean']:.3f}")

# ====================== Validate: IS3 signature in CD8 cells ==================
print("\n=== Validation 1: IS3 signature enrichment in CD8 T cells ===")
ep_v_cd8 = cells[cells["cell_type"].isin(["Epithelial", "CD8_T"])].copy()
from scipy.stats import mannwhitneyu
for sig in ["IS3_score", "IS5_score"]:
    e = ep_v_cd8[ep_v_cd8["cell_type"] == "Epithelial"][sig].dropna()
    c = ep_v_cd8[ep_v_cd8["cell_type"] == "CD8_T"][sig].dropna()
    if len(e) > 5 and len(c) > 5:
        u, p = mannwhitneyu(e, c, alternative="two-sided")
        print(f"  {sig}: Epithelial mean={e.mean():.3f}, CD8_T mean={c.mean():.3f}, P={p:.2e}")

print("\n=== Validation 2: VTCN1 expression per cell type ===")
v_by_ct = cells.groupby("cell_type")["expr_VTCN1"].mean().sort_values(ascending=False)
print(v_by_ct.round(3).to_string())

# ====================== Generate Fig 17 =========================================
print("\nGenerating Figure 17 (single-cell validation) ...")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PALETTE = {
    "blue_main": "#3B6FB6",
    "red_strong": "#C73E3A",
    "green_strong": "#3A8C5F",
    "purple_strong": "#7E4FA0",
    "orange_strong": "#E08A3C",
    "yellow_strong": "#D4B53C",
    "neutral_mid": "#999999",
    "neutral_dark": "#444444",
    "neutral_black": "#000000",
}
CT_COLOR = {
    "Epithelial": PALETTE["red_strong"],
    "CD8_T": PALETTE["blue_main"],
    "CD4_T": PALETTE["green_strong"],
    "Treg": PALETTE["purple_strong"],
    "B_cell": PALETTE["orange_strong"],
    "Plasma": "#C97B63",
    "Macrophage": PALETTE["yellow_strong"],
    "Monocyte": "#A3823C",
    "DC": "#5C9C9C",
    "NK": "#7B6FA0",
    "Fibroblast": "#808080",
    "Endothelial": "#406060",
    "Mast": "#A04060",
    "Unknown": "#CCCCCC",
}

fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.0))

# Panel labels added after tight_layout
panel_labels = {}

# ----- Panel a: cell-type fractions per sample -----
ax = axes[0]
cts_to_plot = ["Epithelial", "CD8_T", "CD4_T", "Treg", "B_cell", "Plasma",
               "Macrophage", "Monocyte", "DC", "NK", "Fibroblast"]
per_tumor_long = (cells.groupby(["sample", "grade", "cell_type"])
                  .size().reset_index(name="n"))
per_tumor_long = per_tumor_long[per_tumor_long["cell_type"].isin(cts_to_plot)]
per_tumor_long["frac"] = per_tumor_long.groupby("sample")["n"].transform(lambda x: x / x.sum())
samples_order = ["LEEC1", "MEEC1", "HEEC1"]
x = np.arange(len(samples_order))
bottom = np.zeros(len(samples_order))
for ct in cts_to_plot:
    vals = []
    for s in samples_order:
        sub = per_tumor_long[(per_tumor_long["sample"] == s) & (per_tumor_long["cell_type"] == ct)]
        vals.append(sub["frac"].values[0] if len(sub) else 0)
    ax.bar(x, vals, bottom=bottom, color=CT_COLOR[ct], label=ct, width=0.55)
    bottom = bottom + np.array(vals)
ax.set_xticks(x)
ax.set_xticklabels([f"{s}\n({SAMPLE_GRADE[s]})" for s in samples_order], fontsize=7)
ax.set_ylabel("Cell-type fraction", fontsize=8)
ax.set_title("Composition across grade groups", fontsize=8.5, fontweight="bold", pad=3)
ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=5.3, frameon=False, ncol=1)
ax.tick_params(axis="y", labelsize=7)
panel_labels[ax] = "a"

# ----- Panel b: IS3 vs IS5 score by cell type -----
ax = axes[1]
ct_main = ["CD8_T", "NK", "Macrophage", "CD4_T", "B_cell", "Treg", "Epithelial", "Fibroblast"]
data = []
n_counts = []
for ct in ct_main:
    sub = cells[cells["cell_type"] == ct]
    data.append(sub["IS3_score"].values - sub["IS5_score"].values)
    n_counts.append(len(sub))
parts = ax.violinplot(data, showmedians=True, showextrema=False, widths=0.7)
for pc, ct in zip(parts["bodies"], ct_main):
    pc.set_facecolor(CT_COLOR[ct])
    pc.set_alpha(0.7)
parts["cmedians"].set_color(PALETTE["neutral_black"])
parts["cmedians"].set_linewidth(1.4)
ax.axhline(0, color=PALETTE["neutral_mid"], linestyle="--", linewidth=0.8)
ax.set_xticks(range(1, len(ct_main) + 1))
labels = [ct.replace("_", " ") for ct in ct_main]
ax.set_xticklabels(labels, fontsize=6.8, rotation=30, ha="right")
# Add n counts below x-axis
for i, n in enumerate(n_counts):
    ax.text(i + 1, -3.8, f"n={n:,}", ha="center", va="top", fontsize=5.5,
            color=PALETTE["neutral_dark"])
ax.set_ylim(bottom=-4.2)
ax.set_ylabel("IS3 score − IS5 score", fontsize=8)
ax.set_title("IS3-skew by cell type", fontsize=8.5, fontweight="bold", pad=3)
ax.tick_params(axis="y", labelsize=7)
panel_labels[ax] = "b"

# ----- Panel c: key antigen expression per cell type -----
ax = axes[2]
plot_genes = ["VTCN1", "CD274"]
ct_comp = ["Epithelial", "Macrophage", "CD8_T", "Treg", "Fibroblast"]
x = np.arange(len(plot_genes))
width = 0.13
for i, ct in enumerate(ct_comp):
    vals = []
    for g in plot_genes:
        sub = cells[cells["cell_type"] == ct]
        vals.append(sub[f"expr_{g}"].mean() if len(sub) > 10 else 0)
    ax.bar(x + (i - 2) * width, vals, width, color=CT_COLOR[ct], label=ct)
ax.set_xticks(x)
ax.set_xticklabels(plot_genes, fontsize=8)
ax.set_ylabel("Mean log1p CPM", fontsize=8)
ax.set_title("Antigen expression by cell type", fontsize=8.5, fontweight="bold", pad=3)
ax.legend(loc="upper right", fontsize=5.8, frameon=False, ncol=1)
ax.tick_params(axis="y", labelsize=7)
panel_labels[ax] = "c"

plt.tight_layout()

# Add panel labels after tight_layout to avoid overlap
for ax, label in panel_labels.items():
    pos = ax.get_position()
    fig.text(pos.x0 - 0.015, pos.y1 + 0.02, label, fontsize=12,
             fontweight="bold", ha="center", va="center")

plt.savefig(FIG_DIR / "fig17_sc_validation.png", dpi=300, bbox_inches="tight")
plt.savefig(FIG_DIR / "fig17_sc_validation.pdf", bbox_inches="tight")
plt.savefig(FIG_DIR / "fig17_sc_validation.svg", bbox_inches="tight")
plt.close()
print(f"  Saved fig17_sc_validation.{{png,pdf,svg}}")