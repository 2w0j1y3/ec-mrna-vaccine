# -*- coding: utf-8 -*-
"""解析 GSE120490 family.soft.gz（145 例子宫内膜样 EC，NRG/GOG 队列）。

输出:
- gse120490_expr.csv (gene x sample, log2)
- gse120490_clinical.csv
"""
import os, gzip
import numpy as np
import pandas as pd

DATA = os.path.dirname(os.path.abspath(__file__))

samples = []          # [(gsm, {field: value})]
sample_expr = {}      # gsm -> {probe: value}
cur_gsm = None
cur_meta = {}
in_table = False

print("Parsing family.soft.gz ...")
with gzip.open(os.path.join(DATA, "GSE120490_family.soft.gz"), "rt",
               encoding="utf-8", errors="ignore") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("^SAMPLE"):
            # 保存上一个样本
            if cur_gsm is not None:
                samples.append((cur_gsm, cur_meta))
            cur_gsm = line.split("=")[1].strip()
            cur_meta = {}
            in_table = False
        elif line.startswith("!sample_table_begin"):
            in_table = True
            sample_expr[cur_gsm] = {}
        elif line.startswith("!sample_table_end"):
            in_table = False
        elif in_table:
            parts = line.split("\t")
            if len(parts) >= 2 and parts[0] not in ("ID_REF",):
                try:
                    sample_expr[cur_gsm][parts[0]] = float(parts[1])
                except (ValueError, IndexError):
                    pass
        elif line.startswith("!Sample_characteristics_ch1"):
            val = line.split("=", 1)[1].strip() if "=" in line else ""
            if ":" in val:
                field, v = val.split(":", 1)
                cur_meta[field.strip()] = v.strip()
        elif line.startswith("!Sample_geo_accession") and cur_gsm is None:
            pass

if cur_gsm is not None:
    samples.append((cur_gsm, cur_meta))

print(f"Parsed {len(samples)} samples")
clin = pd.DataFrame([{"sample_id": g, **m} for g, m in samples])
print(f"Clinical fields: {[c for c in clin.columns if c != 'sample_id']}")
for c in clin.columns:
    if c != "sample_id":
        vc = clin[c].value_counts()
        if len(vc) <= 6:
            print(f"  {c}: {dict(vc)}")

# ---------- 探针矩阵 ----------
print("\nBuilding probe matrix ...")
all_probes = set()
for gsm in sample_expr:
    all_probes.update(sample_expr[gsm].keys())
probes = sorted(all_probes)
gsm_list = sorted(sample_expr.keys())
print(f"Probes: {len(probes)}, Samples: {len(gsm_list)}")

probe_idx = {p: i for i, p in enumerate(probes)}
mat = np.full((len(probes), len(gsm_list)), np.nan, dtype=np.float32)
for j, gsm in enumerate(gsm_list):
    d = sample_expr[gsm]
    for p, v in d.items():
        mat[probe_idx[p], j] = v

# ---------- GPL570 注释 ----------
print("Reading GPL570 annotation ...")
annot = pd.read_csv(os.path.join(DATA, "GPL570.annot.gz"), sep="\t", skiprows=27)
annot = annot.rename(columns={"ID": "probe", "Gene symbol": "symbol"})
annot = annot[["probe", "symbol"]].dropna(subset=["symbol"])
annot = annot[annot["symbol"] != ""]
annot["probe"] = annot["probe"].astype(str)
annot["symbols"] = annot["symbol"].str.split("///").apply(
    lambda x: [s.strip() for s in x if s.strip()])
probe_to_symbols = dict(zip(annot["probe"], annot["symbols"]))

common = [p for p in probes if p in probe_to_symbols]
print(f"Annotated probes: {len(common)}/{len(probes)}")

# 每基因取均值最高探针
mat_common = mat[[probe_idx[p] for p in common]]
probe_mean = np.nanmean(mat_common, axis=1)
gene_probes = {}
for k, p in enumerate(common):
    for g in probe_to_symbols[p]:
        gene_probes.setdefault(g, []).append(k)

rows = []
for g, idxs in gene_probes.items():
    best = max(idxs, key=lambda i: probe_mean[i])
    rows.append((g, best))
gene_arr = np.array([r[0] for r in rows])
sel = np.array([r[1] for r in rows])
gene_mat = mat_common[sel]

df = pd.DataFrame(gene_mat, index=gene_arr, columns=gsm_list)
df = df[~df.index.duplicated(keep="first")]
print(f"Gene-level matrix: {df.shape}")

df.to_csv(os.path.join(DATA, "gse120490_expr.csv"))
clin.to_csv(os.path.join(DATA, "gse120490_clinical.csv"), index=False)
print("Saved gse120490_expr.csv + gse120490_clinical.csv")
