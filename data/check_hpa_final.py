# -*- coding: utf-8 -*-
"""HPA 最终核查：6 抗原（正确列代码）。

- t_RNA_endometrium_1: 正常子宫内膜 nTPM
- t_RNA_<tissue>: 各正常组织（疫苗安全性扫描）
- cell_RNA_uterine_cancer: EC 细胞系表达
- prognostic_UCEC: TCGA 预后分型
- rnats/rnatsm: 组织特异性
"""
import os, gzip
import requests
import pandas as pd
from io import StringIO

DATA = os.path.dirname(os.path.abspath(__file__))
GENES = ["CLDN6", "CCNE1", "MAL", "CTSV", "VTCN1", "MUC16"]

# 正常组织列表（HPA consensus，含 GTEx）
TISSUES = [
    "endometrium_1", "ovary", "fallopian_tube", "breast", "cervix", "placenta",
    "colon", "small_intestine", "stomach_1", "liver", "pancreas", "kidney",
    "lung", "heart_muscle", "skeletal_muscle", "brain_cortex",
    "cerebral_cortex", "brain", "bone_marrow", "lymph_node", "spleen", "thymus",
    "blood", "testis", "prostate", "urinary_bladder", "skin_1", "thyroid_gland",
    "adrenal_gland", "pituitary_gland",
]

cols = ["g", "gd", "rnats", "rnatsm",
        "prognostic_Uterine_Corpus_Endometrial_Carcinoma_(TCGA)",
        "cell_RNA_uterine_cancer"] + [f"t_RNA_{t}" for t in TISSUES]

def hpa_query(genes):
    params = {
        "search": " ".join(genes),
        "format": "tsv",
        "columns": ",".join(cols),
    }
    r = requests.get("https://www.proteinatlas.org/api/search_download.php",
                     params=params, timeout=60)
    try:
        text = gzip.decompress(r.content).decode("utf-8")
    except Exception:
        text = r.content.decode("utf-8", errors="ignore")
    return pd.read_csv(StringIO(text), sep="\t")

frames = []
for gene in GENES:
    params = {
        "search": gene,
        "format": "tsv",
        "columns": ",".join(cols),
    }
    r = requests.get("https://www.proteinatlas.org/api/search_download.php",
                     params=params, timeout=60)
    try:
        text = gzip.decompress(r.content).decode("utf-8")
    except Exception:
        text = r.content.decode("utf-8", errors="ignore")
    from io import StringIO
    sub = pd.read_csv(StringIO(text), sep="\t")
    sub = sub[sub["Gene"] == gene]
    frames.append(sub)
    print(f"  {gene}: {len(sub)} row(s)")

df = pd.concat(frames, ignore_index=True)
df.to_csv(os.path.join(DATA, "hpa_antigen_check.csv"), index=False)

pd.set_option("display.max_columns", 40)
pd.set_option("display.width", 300)
print("=== 正常组织 RNA (nTPM) ===")
tissue_cols = [c for c in df.columns if c.startswith("Tissue RNA - ")]
key_view = df[["Gene"] + [c for c in tissue_cols if any(k in c.lower() for k in
              ["endometri", "ovary", "breast", "colon", "liver", "kidney", "lung",
               "bone marrow", "lymph node", "spleen", "testis", "brain", "heart", "skeletal"])]]
print(key_view.to_string(index=False))

print("\n=== EC 细胞系 RNA ===")
cell_col = [c for c in df.columns if "uterine" in c.lower() and "cell" in c.lower()]
if cell_col:
    print(df[["Gene"] + cell_col].to_string(index=False))

print("\n=== TCGA 预后 ===")
prog_col = [c for c in df.columns if "prognostic" in c.lower() or "Endometrial Carcinoma" in c]
if prog_col:
    print(df[["Gene"] + prog_col].to_string(index=False))

print("\n=== 正常组织安全性（max nTPM 除子宫内膜外）===")
for _, row in df.iterrows():
    vals = {}
    for c in tissue_cols:
        if "endometr" in c.lower():
            continue
        v = pd.to_numeric(row[c], errors="coerce")
        if not pd.isna(v) and v == v:
            vals[c.replace("Tissue RNA - ", "")] = v
    if vals:
        top3 = sorted(vals.items(), key=lambda kv: -kv[1])[:3]
        print(f"  {row['Gene']:8s}: " + ", ".join(f"{t}={v:.0f}" for t, v in top3))

print("\nSaved hpa_antigen_check.csv")
