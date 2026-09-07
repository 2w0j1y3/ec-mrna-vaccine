"""Derive IS3 and IS5 up-regulated gene signatures from TCGA-UCEC bulk expression,
to be used for single-cell projection validation.

For each of IS3 vs rest and IS5 vs rest:
  - Mann-Whitney U test (non-parametric) per gene
  - Rank by log2FC and BH-FDR
  - Take top 50 up-regulated in IS3 and top 50 up-regulated in IS5 as signatures

Output:
  data/is3_signature.csv (gene, log2FC, pval, fdr)
  data/is5_signature.csv (gene, log2FC, pval, fdr)
  data/is3_is5_signatures.json (combined dict for downstream use)
"""
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def bh_fdr(pvals):
    """Manual Benjamini-Hochberg FDR (no statsmodels dependency)."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    fdr = ranked * n / (np.arange(n) + 1)
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    fdr = np.clip(fdr, 0, 1)
    out = np.empty_like(fdr)
    out[order] = fdr
    return out

DATA = Path(__file__).parent

# Load bulk expression matrix
print("Loading TCGA bulk expression ...")
tpm = pd.read_csv(DATA / "ucec_tpm_matrix.csv", index_col=0)
print(f"  {tpm.shape[1]} samples × {tpm.shape[0]} genes")

# Load subtype assignments (from analysis_subtype.py output)
sub = pd.read_csv(DATA / "ucec_subtypes.csv")
print(f"  subtype assignments: {sub.shape}")
# Map sample -> subtype
smap = dict(zip(sub["sample_id"], sub["subtype"]))
tumor_cols = [c for c in tpm.columns if smap.get(c)]
print(f"  tumor samples with subtype: {len(tumor_cols)}")

tpm_tumor = tpm[tumor_cols]
# Subtypes
labels = pd.Series([smap[c] for c in tumor_cols], index=tumor_cols)
print("  subtype distribution:")
print(labels.value_counts().sort_index())

# Compute per-gene IS3 vs rest and IS5 vs rest
print("\nComputing per-gene differential expression ..")
results = {}
for sub_target in ["IS3", "IS5"]:
    pos = labels[labels == sub_target].index
    neg = labels[labels != sub_target].index
    print(f"  {sub_target}: {len(pos)} vs {len(neg)}")
    rows = []
    for gene in tpm_tumor.index:
        v_pos = tpm_tumor.loc[gene, pos].values
        v_neg = tpm_tumor.loc[gene, neg].values
        # log2FC using mean TPM
        m_pos = v_pos.mean()
        m_neg = v_neg.mean()
        # Add small pseudocount for log safety
        log2fc = np.log2((m_pos + 0.1) / (m_neg + 0.1))
        try:
            u, p = mannwhitneyu(v_pos, v_neg, alternative="two-sided")
        except ValueError:
            u, p = np.nan, 1.0
        rows.append((gene, m_pos, m_neg, log2fc, p, u))
    df = pd.DataFrame(rows, columns=["gene", "mean_pos", "mean_neg", "log2FC", "pvalue", "u_stat"])
    # BH-FDR
    mask = df["pvalue"].notna()
    df["FDR"] = 1.0
    df.loc[mask, "FDR"] = bh_fdr(df.loc[mask, "pvalue"].values)
    df = df.sort_values(["log2FC", "FDR"], ascending=[False, True]).reset_index(drop=True)
    results[sub_target] = df
    df.to_csv(DATA / f"{sub_target.lower()}_de_all.csv", index=False)
    up = df[(df["log2FC"] > 0.5) & (df["FDR"] < 0.05)].head(50)
    print(f"  {sub_target}: top 50 up-regulated (log2FC>0.5 & FDR<0.05):")
    print(up[["gene", "log2FC", "FDR"]].head(10).to_string(index=False))
    up.to_csv(DATA / f"{sub_target.lower()}_signature.csv", index=False)

# Combined signature dict
sigs = {
    "IS3_up": results["IS3"][(results["IS3"]["log2FC"] > 0.5) & (results["IS3"]["FDR"] < 0.05)].head(50)["gene"].tolist(),
    "IS5_up": results["IS5"][(results["IS5"]["log2FC"] > 0.5) & (results["IS5"]["FDR"] < 0.05)].head(50)["gene"].tolist(),
    "IS3_down": results["IS3"][(results["IS3"]["log2FC"] < -0.5) & (results["IS3"]["FDR"] < 0.05)].head(50)["gene"].tolist(),
    "IS5_down": results["IS5"][(results["IS5"]["log2FC"] < -0.5) & (results["IS5"]["FDR"] < 0.05)].head(50)["gene"].tolist(),
}
with open(DATA / "is3_is5_signatures.json", "w") as f:
    json.dump(sigs, f, indent=2)

print(f"\nSignatures saved:")
print(f"  IS3_up: {len(sigs['IS3_up'])} genes")
print(f"  IS5_up: {len(sigs['IS5_up'])} genes")
print(f"  IS3_down: {len(sigs['IS3_down'])} genes")
print(f"  IS5_down: {len(sigs['IS5_down'])} genes")