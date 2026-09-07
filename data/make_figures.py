# -*- coding: utf-8 -*-
"""UCEC (EC) 疫苗抗原 + 免疫亚型 发表级配图（nature-figure Python 后端，适配 COAD pipeline）。

图1 抗原鉴定：火山图 + 候选抗原排名
图2 免疫亚型：一致性聚类共识矩阵 + 亚型生存曲线
图3 亚型特征：特征热图 + TMB/免疫评分箱线图
图4 抗原生存：top6 抗原 KM 曲线
图13 GEP: GEP z-score × 免疫亚型
图14 HLA: HLA 表位覆盖热图
图15 TIME + cDC1: 4-layer 算法 + 相关性
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.dirname(SCRIPT_DIR) if os.path.basename(SCRIPT_DIR) == "scripts" else SCRIPT_DIR
FIGDIR = os.path.join(os.path.dirname(DATA) if os.path.basename(DATA) == "data" else DATA, "figures")
if not os.path.exists(FIGDIR):
    FIGDIR = os.path.join(DATA, "..", "figures")
FIGDIR = os.path.abspath(FIGDIR)
os.makedirs(FIGDIR, exist_ok=True)

# ── 发表级样式 ─────────────────────────────────────────
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 7
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['legend.frameon'] = False

PALETTE = {
    "blue_main": "#0F4D92", "green_3": "#8BCF8B", "red_strong": "#B64342",
    "teal": "#42949E", "violet": "#9A4D8E", "gold": "#FFD700",
    "neutral_mid": "#767676", "neutral_dark": "#4D4D4D", "neutral_black": "#272727",
}
SUBTYPE_COLORS = ["#0F4D92", "#8BCF8B", "#B64342", "#9A4D8E", "#42949E"]


def add_panel_label(ax, label, x=-0.06, y=1.08):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=9,
            fontweight='bold', ha='left', va='bottom', color='black')


# ── 载入数据 ─────────────────────────────────────────
print("载入分析数据 ...")
deg = pd.read_csv(os.path.join(DATA, "ucec_degs.csv"))
final = pd.read_csv(os.path.join(DATA, "ucec_final_antigens.csv"))
m = pd.read_csv(os.path.join(DATA, "ucec_tpm_matrix.csv"), index_col=0)
smap = pd.read_csv(os.path.join(DATA, "ucec_sample_map.csv"))
clin = pd.read_csv(os.path.join(DATA, "ucec_clinical.csv"))
subtype = pd.read_csv(os.path.join(DATA, "ucec_subtypes.csv"))
features = pd.read_csv(os.path.join(DATA, "ucec_subtype_features.csv"))
consensus = pd.read_csv(os.path.join(DATA, "ucec_consensus_matrix.csv"), index_col=0)

# 生存数据
clin = clin[clin["vital_status"].isin(["Alive", "Dead"])].dropna(subset=["vital_status"])
clin["event"] = (clin["vital_status"] == "Dead").astype(int)
clin["time"] = np.where(clin["event"] == 1, clin["days_to_death"], clin["days_to_last_follow_up"])
clin = clin[clin["time"] > 0]
clin["time_mon"] = clin["time"] / 30.44
surv = clin.set_index("submitter_id")[["time_mon", "event"]]
sample_to_case = dict(zip(smap["sample_id"], smap["case_id"]))
tumor = [s for s in smap[smap["sample_type"] == "Primary Tumor"]["sample_id"] if s in m.columns]


def km_data(gene):
    exp = m.loc[gene, tumor].astype(float)
    med = exp.median()
    hi = exp[exp > med].index.map(sample_to_case).dropna()
    lo = exp[exp <= med].index.map(sample_to_case).dropna()
    hi = [c for c in hi if c in surv.index]
    lo = [c for c in lo if c in surv.index]
    th = surv.loc[hi, "time_mon"].values; eh = surv.loc[hi, "event"].values
    tl = surv.loc[lo, "time_mon"].values; el = surv.loc[lo, "event"].values
    p = logrank_test(th, tl, event_observed_A=eh, event_observed_B=el).p_value
    return (th, eh), (tl, el), p


def draw_km(ax, gene):
    (th, eh), (tl, el), p = km_data(gene)
    kmh = KaplanMeierFitter(); kmh.fit(th, eh, label="High")
    kml = KaplanMeierFitter(); kml.fit(tl, el, label="Low")
    kmh.plot_survival_function(ax=ax, color=PALETTE["red_strong"], linewidth=1.6, ci_show=False)
    kml.plot_survival_function(ax=ax, color=PALETTE["blue_main"], linewidth=1.6, ci_show=False)
    ax.set_title(gene, fontsize=8, fontweight='bold', pad=3)
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Overall survival")
    ax.set_ylim(0, 1.05)
    ax.text(0.97, 0.05, f"log-rank P = {p:.3g}", transform=ax.transAxes,
            ha='right', va='bottom', fontsize=6.5)
    ax.legend(fontsize=6.5, loc="lower left")


# ════════════════ 图1：抗原鉴定 ════════════════
print("绘制图1：抗原鉴定 ...")
fig1, axes1 = plt.subplots(1, 2, figsize=(7.0, 3.1))

# 修订后 6 抗原（数据驱动 4 + 转化锚定 2）
REVISED_PANEL = ["CLDN6", "CCNE1", "MAL", "CTSV", "VTCN1", "MUC16"]

# (a) 火山图
ax = axes1[0]
strict_genes = set(final[(final["prognostic"]) & (final["apc_positive"])]["gene"])
cand_genes = set(final["gene"])
x = deg["log2FC"].values
y = -np.log10(deg["pvalue"].values)
colors = np.where(
    [g in strict_genes for g in deg["gene"]], PALETTE["red_strong"],
    np.where([g in cand_genes for g in deg["gene"]], PALETTE["blue_main"], "#C8C8C8"))
ax.scatter(x, y, s=4, c=colors, linewidths=0, alpha=0.7, rasterized=True)
ax.axvline(1, color=PALETTE["neutral_mid"], linestyle='--', linewidth=0.8)
ax.axvline(-1, color=PALETTE["neutral_mid"], linestyle='--', linewidth=0.8)
ax.axhline(-np.log10(0.05), color=PALETTE["neutral_mid"], linestyle='--', linewidth=0.8)
# 标注修订 6 抗原（手工排布 + 细引线，参照 COAD fig1 的防碰撞经验）
_LABEL_POS = {
    "CLDN6":  (5.3, 21.0, "left"),   # log2FC=5.78（最右）
    "CCNE1":  (2.7, 27.5, "left"),
    "CTSV":   (3.5, 32.5, "left"),
    "MAL":    (2.5, 38.0, "left"),
    "VTCN1":  (0.2, 34.5, "left"),
    "MUC16":  (-0.5, 25.5, "right"),
}
for g in REVISED_PANEL:
    row = deg[deg["gene"] == g]
    if len(row):
        x_g = row["log2FC"].values[0]
        y_g = -np.log10(row["pvalue"].values[0])
        lx, ly, ha = _LABEL_POS[g]
        ax.annotate(g, xy=(x_g, y_g), xytext=(lx, ly), textcoords='data',
                    ha=ha, va='center', fontsize=6.5, color=PALETTE["neutral_black"],
                    arrowprops=dict(arrowstyle='-', color=PALETTE["neutral_mid"], lw=0.5))
ax.set_xlabel("$\\log_{2}$ fold change (tumor vs. normal)", fontsize=8)
ax.set_ylabel("$-\\log_{10}$(P value)", fontsize=8)
ax.set_title("Differential expression", fontsize=8, fontweight='bold', pad=3)
add_panel_label(ax, "a")

# (b) 修订 6 抗原 log2FC 条形图
ax = axes1[1]
panel_rows = final[final["gene"].isin(REVISED_PANEL)].set_index("gene").loc[REVISED_PANEL]
genes = panel_rows.index.values[::-1]
vals = panel_rows["log2FC"].values[::-1]
bar_colors = [PALETTE["red_strong"] if g in ["CLDN6", "CCNE1", "MAL", "CTSV"] else PALETTE["teal"] for g in genes]
ax.barh(range(len(genes)), vals, color=bar_colors, height=0.7)
ax.set_yticks(range(len(genes)))
ax.set_yticklabels(genes, fontsize=6.5)
ax.set_xlabel("$\\log_{2}$ fold change", fontsize=8)
ax.set_title("Final 6-antigen panel", fontsize=8, fontweight='bold', pad=3)
ax.axvline(0, color=PALETTE["neutral_black"], linewidth=0.7)
# 标注 KM P 值
for i, g in enumerate(genes):
    p_km = panel_rows.loc[g, "km_logrank_p"]
    if not np.isnan(p_km):
        ax.text(vals[i] + 0.08, i, f"P={p_km:.1g}",
                va='center', fontsize=5.5, color=PALETTE["neutral_dark"])
add_panel_label(ax, "b")

fig1.tight_layout(pad=1.0)
fig1.savefig(f"{FIGDIR}\\fig1_antigen_identification.svg", bbox_inches="tight")
fig1.savefig(f"{FIGDIR}\\fig1_antigen_identification.pdf", bbox_inches="tight")
fig1.savefig(f"{FIGDIR}\\fig1_antigen_identification.png", dpi=300, bbox_inches="tight")
plt.close(fig1)
print("  图1 完成")

# ════════════════ 图2：免疫亚型 ════════════════
print("绘制图2：免疫亚型 ...")
fig2, axes2 = plt.subplots(1, 2, figsize=(7.0, 3.1))

# (a) 共识矩阵热图
ax = axes2[0]
subtype_sorted = subtype.sort_values("subtype")
order = subtype_sorted["sample_id"].tolist()
cons = consensus.loc[order, order].values
cmap = LinearSegmentedColormap.from_list("cons", ["#FFFFFF", PALETTE["blue_main"]])
im = ax.imshow(cons, cmap=cmap, aspect='auto', vmin=0, vmax=1)
boundaries = np.cumsum(subtype_sorted["subtype"].value_counts().sort_index().values)
for b in boundaries[:-1]:
    ax.axhline(b - 0.5, color=PALETTE["red_strong"], linewidth=0.8)
    ax.axvline(b - 0.5, color=PALETTE["red_strong"], linewidth=0.8)
ax.set_xticks([]); ax.set_yticks([])
ax.set_title("Consensus matrix (k = 5)", fontsize=8, fontweight='bold', pad=3)
cbar = fig2.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cbar.set_label("Consensus", fontsize=6.5)
cbar.ax.tick_params(labelsize=6)
add_panel_label(ax, "a")

# (b) 亚型生存曲线
ax = axes2[1]
sub_surv = {}
for isub, g in subtype.groupby("subtype"):
    cases = g["sample_id"].map(sample_to_case).dropna()
    cases = [c for c in cases if c in surv.index]
    t = surv.loc[cases, "time_mon"].values
    e = surv.loc[cases, "event"].values
    sub_surv[isub] = (t, e)
for i, isub in enumerate(sorted(sub_surv, key=lambda s: int(s[2:]))):
    t, e = sub_surv[isub]
    km = KaplanMeierFitter()
    med_os = km.fit(t, e).median_survival_time_
    med_str = f", med {med_os:.0f} mo" if not np.isnan(med_os) else ""
    km.plot_survival_function(ax=ax, label=f"{isub} (n={len(t)}{med_str})",
                              color=SUBTYPE_COLORS[i], linewidth=1.5, ci_show=False)
ax.set_xlabel("Time (months)")
ax.set_ylabel("Overall survival")
ax.set_ylim(0, 1.05)
ax.set_title("Survival by immune subtype", fontsize=8, fontweight='bold', pad=3)
ax.legend(fontsize=6, loc="lower left")
t_all = np.concatenate([sub_surv[s][0] for s in sorted(sub_surv, key=lambda x: int(x[2:]))])
e_all = np.concatenate([sub_surv[s][1] for s in sorted(sub_surv, key=lambda x: int(x[2:]))])
g_all = np.concatenate([[s] * len(sub_surv[s][0]) for s in sorted(sub_surv, key=lambda x: int(x[2:]))])
p_mlr = multivariate_logrank_test(t_all, g_all, e_all).p_value
ax.text(0.97, 0.05, f"log-rank P = {p_mlr:.3g}", transform=ax.transAxes,
        ha='right', va='bottom', fontsize=6.5)
add_panel_label(ax, "b")

fig2.tight_layout(pad=1.0)
fig2.savefig(f"{FIGDIR}\\fig2_immune_subtypes.svg", bbox_inches="tight")
fig2.savefig(f"{FIGDIR}\\fig2_immune_subtypes.pdf", bbox_inches="tight")
fig2.savefig(f"{FIGDIR}\\fig2_immune_subtypes.png", dpi=300, bbox_inches="tight")
plt.close(fig2)
print("  图2 完成")

# ════════════════ 图3：亚型特征 ════════════════
print("绘制图3：亚型特征 ...")
fig3, axes3 = plt.subplots(1, 2, figsize=(7.0, 3.1))

# (a) 特征热图（加入 POLE / MSI-H 比例行）
ax = axes3[0]
feat_cols = ["TMB", "immune_score", "CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT",
             "VTCN1", "CD276",
             "CD8 T_score", "Treg_score", "Cytotoxic_score", "NK_score",
             "B cell_score", "Macrophage_score", "DC_score"]
feat_cols = [c for c in feat_cols if c in features.columns]
mat = features.groupby("subtype")[feat_cols].mean()
mat = mat.loc[[f"IS{i}" for i in range(1, 6)]]

# 追加分子分型比例行（% POLE, % MSI-H, % hypermutated）
mol_path = os.path.join(DATA, "ucec_molecular_class.csv")
if os.path.exists(mol_path):
    mol = pd.read_csv(mol_path)
    sub2case = dict(zip(smap["sample_id"], smap["case_id"]))
    feats_cases = features.copy()
    feats_cases["case_id"] = feats_cases["sample_id"].map(sub2case)
    feats_cases = feats_cases.merge(mol[["case_id", "molecular_class"]], on="case_id", how="left")
    pct_pole = feats_cases.groupby("subtype")["molecular_class"].apply(lambda s: (s == "POLE").mean() * 100)
    pct_msih = feats_cases.groupby("subtype")["molecular_class"].apply(lambda s: (s == "MSI-H").mean() * 100)
    pct_hyper = feats_cases.groupby("subtype")["molecular_class"].apply(
        lambda s: s.isin(["POLE", "MSI-H"]).mean() * 100)
    mat["% POLE"] = pct_pole
    mat["% MSI-H"] = pct_msih
    mat["% Hypermut"] = pct_hyper
    feat_cols = feat_cols + ["% POLE", "% MSI-H", "% Hypermut"]
    mat = mat[feat_cols].loc[[f"IS{i}" for i in range(1, 6)]]

mat_z = (mat - mat.mean()) / (mat.std() + 1e-8)
im = ax.imshow(mat_z.T.values, cmap="RdBu_r", aspect='auto', vmin=-2, vmax=2)
ax.set_xticks(range(5))
ax.set_xticklabels(mat.index, fontsize=7)
ax.set_yticks(range(len(feat_cols)))
ax.set_yticklabels(feat_cols, fontsize=6)
ax.set_title("Immune subtype features", fontsize=8, fontweight='bold', pad=3)
cbar = fig3.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cbar.set_label("Z-score", fontsize=6.5)
cbar.ax.tick_params(labelsize=6)
add_panel_label(ax, "a")

# (b) TMB 箱线图
ax = axes3[1]
subs = [f"IS{i}" for i in range(1, 6)]
tmb_data = [features[features["subtype"] == s]["TMB"].dropna().values for s in subs]
bp = ax.boxplot(tmb_data, tick_labels=subs, patch_artist=True, widths=0.6,
                medianprops=dict(color="black", linewidth=1.2),
                flierprops=dict(marker='o', markersize=2, alpha=0.5))
for patch, color in zip(bp['boxes'], SUBTYPE_COLORS):
    patch.set_facecolor(color); patch.set_alpha(0.7)
ax.set_ylabel("TMB (mutations/Mb)")
ax.set_title("Tumor mutation burden", fontsize=8, fontweight='bold', pad=3)
add_panel_label(ax, "b")

fig3.tight_layout(pad=1.0)
fig3.savefig(f"{FIGDIR}\\fig3_subtype_features.svg", bbox_inches="tight")
fig3.savefig(f"{FIGDIR}\\fig3_subtype_features.pdf", bbox_inches="tight")
fig3.savefig(f"{FIGDIR}\\fig3_subtype_features.png", dpi=300, bbox_inches="tight")
plt.close(fig3)
print("  图3 完成")

# ════════════════ 图4：抗原生存（修订 6 抗原） ════════════════
print("绘制图4：抗原生存 ...")
top6 = [g for g in ["CLDN6", "CCNE1", "MAL", "CTSV", "VTCN1", "MUC16"] if g in m.index]
fig4, axes4 = plt.subplots(2, 3, figsize=(7.0, 4.8))
for ax, g, lab in zip(axes4.flat, top6, "abcdef"):
    draw_km(ax, g)
    add_panel_label(ax, lab)
fig4.tight_layout(pad=1.0)
fig4.savefig(f"{FIGDIR}\\fig4_antigen_survival.svg", bbox_inches="tight")
fig4.savefig(f"{FIGDIR}\\fig4_antigen_survival.pdf", bbox_inches="tight")
fig4.savefig(f"{FIGDIR}\\fig4_antigen_survival.png", dpi=300, bbox_inches="tight")
plt.close(fig4)
print("  图4 完成")

# ════════════════ 图13：GEP × 免疫亚型 ════════════════
if os.path.exists(os.path.join(DATA, "ucec_gep.csv")):
    print("绘制图13：GEP ...")
    gep_df = pd.read_csv(os.path.join(DATA, "ucec_gep.csv"))
    gep_df = gep_df.merge(subtype, on="sample_id")

    fig13, axes13 = plt.subplots(1, 2, figsize=(7.0, 3.1))

    # (a) GEP 分布（小提琴图）
    ax = axes13[0]
    subs = [f"IS{i}" for i in range(1, 6)]
    gep_data = [gep_df[gep_df["subtype"] == s]["GEP_z"].values for s in subs]
    parts = ax.violinplot(gep_data, positions=range(len(subs)), widths=0.7,
                          showmeans=False, showmedians=True)
    for i, body in enumerate(parts['bodies']):
        body.set_facecolor(SUBTYPE_COLORS[i])
        body.set_alpha(0.7)
    ax.set_xticks(range(len(subs)))
    ax.set_xticklabels(subs)
    ax.set_ylabel("GEP z-score")
    ax.set_title("T-cell inflamed GEP by subtype", fontsize=8, fontweight='bold', pad=3)

    groups = [g["GEP_z"].values for _, g in gep_df.groupby("subtype")]
    h, p = stats.kruskal(*groups)
    ax.text(0.97, 0.05, f"KW P = {p:.2e}", transform=ax.transAxes,
            ha='right', va='bottom', fontsize=6.5)
    add_panel_label(ax, "a")

    # (b) GEP 阈值比例堆叠图
    ax = axes13[1]
    high_pct = []
    for s in subs:
        d = gep_df[gep_df["subtype"] == s]["GEP_z"]
        high_pct.append((d > 0.5).mean() * 100)
    ax.bar(range(len(subs)), high_pct, color=SUBTYPE_COLORS)
    ax.set_xticks(range(len(subs)))
    ax.set_xticklabels(subs)
    ax.set_ylabel("% GEP^high (z > +0.5)")
    ax.set_title("GEP-responsive fraction", fontsize=8, fontweight='bold', pad=3)
    for i, v in enumerate(high_pct):
        ax.text(i, v + 1, f"{v:.1f}%", ha='center', fontsize=7)
    add_panel_label(ax, "b")

    fig13.tight_layout(pad=1.0)
    fig13.savefig(f"{FIGDIR}\\fig13_gep.svg", bbox_inches="tight")
    fig13.savefig(f"{FIGDIR}\\fig13_gep.pdf", bbox_inches="tight")
    fig13.savefig(f"{FIGDIR}\\fig13_gep.png", dpi=300, bbox_inches="tight")
    plt.close(fig13)
    print("  图13 完成")

# ════════════════ 图14：HLA 表位覆盖 ════════════════
if os.path.exists(os.path.join(DATA, "ucec_epitope_pivot.csv")):
    print("绘制图14：HLA 表位覆盖 ...")
    pivot = pd.read_csv(os.path.join(DATA, "ucec_epitope_pivot.csv"), index_col=0)
    fig14, ax = plt.subplots(figsize=(5.5, 3.1))
    im = ax.imshow(pivot.values, cmap="YlOrRd", aspect='auto')
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha='right', fontsize=7)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_title("HLA-I 9-mer epitope coverage", fontsize=8, fontweight='bold', pad=3)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            if v > 0:
                ax.text(j, i, str(int(v)), ha='center', va='center', fontsize=6,
                        color="white" if v > pivot.values.max() * 0.5 else "black")
    cbar = fig14.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("# predicted epitopes", fontsize=6.5)
    cbar.ax.tick_params(labelsize=6)
    fig14.tight_layout(pad=1.0)
    fig14.savefig(f"{FIGDIR}\\fig14_hla_epitopes.svg", bbox_inches="tight")
    fig14.savefig(f"{FIGDIR}\\fig14_hla_epitopes.pdf", bbox_inches="tight")
    fig14.savefig(f"{FIGDIR}\\fig14_hla_epitopes.png", dpi=300, bbox_inches="tight")
    plt.close(fig14)
    print("  图14 完成")

# ════════════════ 图15：TIME + cDC1 ════════════════
if os.path.exists(os.path.join(DATA, "ucec_time.csv")):
    print("绘制图15：TIME + cDC1 ...")
    time_df = pd.read_csv(os.path.join(DATA, "ucec_time.csv"))
    cdc_df = pd.read_csv(os.path.join(DATA, "ucec_cdc.csv"))
    merged = time_df.merge(cdc_df, on="sample_id").merge(subtype, on="sample_id")

    fig15, axes15 = plt.subplots(1, 3, figsize=(10.5, 3.1))

    # (a) TIME_index × 亚型
    ax = axes15[0]
    subs = [f"IS{i}" for i in range(1, 6)]
    time_data = [merged[merged["subtype"] == s]["TIME_index"].values for s in subs]
    parts = ax.violinplot(time_data, positions=range(len(subs)), widths=0.7,
                          showmeans=False, showmedians=True)
    for i, body in enumerate(parts['bodies']):
        body.set_facecolor(SUBTYPE_COLORS[i])
        body.set_alpha(0.7)
    ax.axhline(0, color=PALETTE["neutral_mid"], linestyle='--', linewidth=0.8)
    ax.set_xticks(range(len(subs)))
    ax.set_xticklabels(subs)
    ax.set_ylabel("TIME conversion index")
    ax.set_title("TIME sensitivity", fontsize=8, fontweight='bold', pad=3)
    add_panel_label(ax, "a")

    # (b) cDC1 × 亚型
    ax = axes15[1]
    cdc_data = [merged[merged["subtype"] == s]["cDC1_score"].values for s in subs]
    parts = ax.violinplot(cdc_data, positions=range(len(subs)), widths=0.7,
                          showmeans=False, showmedians=True)
    for i, body in enumerate(parts['bodies']):
        body.set_facecolor(SUBTYPE_COLORS[i])
        body.set_alpha(0.7)
    ax.set_xticks(range(len(subs)))
    ax.set_xticklabels(subs)
    ax.set_ylabel("cDC1 cross-presentation score")
    ax.set_title("Antigen delivery layer", fontsize=8, fontweight='bold', pad=3)
    add_panel_label(ax, "b")

    # (c) 4-layer optimal cohort Sankey-like display
    ax = axes15[2]
    # EC-specific: IS3 is the hot/optimal subtype (highest GEP, TIME, cDC1)
    n0 = len(merged)
    n1 = int((merged["TIME_index"] > 0).sum())  # TIME-permissive
    n2 = int(((merged["TIME_index"] > 0) & (merged["cDC1_score"] > 0)).sum())
    n3 = int(((merged["TIME_index"] > 0) & (merged["cDC1_score"] > 0) & (merged["subtype"] == "IS3")).sum())
    layers = ["All\npatients", "TIME\npermissive", "TIME ×\ncDC1", "TIME ×\ncDC1 ×\nIS3"]
    counts = [n0, n1, n2, n3]
    pcts = [n0/n0*100, n1/n0*100, n2/n0*100, n3/n0*100]
    bar_colors = [PALETTE["neutral_mid"]] + SUBTYPE_COLORS[:3]
    bars = ax.bar(range(4), pcts, color=bar_colors)
    ax.set_xticks(range(4))
    ax.set_xticklabels(layers, fontsize=7)
    ax.set_ylabel("% of cohort")
    ax.set_title("4-layer cohort enrichment", fontsize=8, fontweight='bold', pad=3)
    for i, (b, c) in enumerate(zip(bars, counts)):
        ax.text(b.get_x() + b.get_width()/2, b.get_height() + 1,
                f"{c} ({pcts[i]:.1f}%)", ha='center', fontsize=7)
    add_panel_label(ax, "c")

    fig15.tight_layout(pad=1.0)
    fig15.savefig(f"{FIGDIR}\\fig15_time_cdc1.svg", bbox_inches="tight")
    fig15.savefig(f"{FIGDIR}\\fig15_time_cdc1.pdf", bbox_inches="tight")
    fig15.savefig(f"{FIGDIR}\\fig15_time_cdc1.png", dpi=300, bbox_inches="tight")
    plt.close(fig15)
    print("  图15 完成")

print("\n=== 所有图完成 ===")
print(f"输出目录: {FIGDIR}")
