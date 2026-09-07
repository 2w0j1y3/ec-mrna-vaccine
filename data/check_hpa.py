# -*- coding: utf-8 -*-
"""HPA 单基因 TSV 核查 6 抗原：EC IHC + 正常组织 + RNA。

流程：Ensembl REST 解析 symbol→ENSG → HPA per-gene TSV → 提取 endometrial 相关行。
"""
import os
import requests
import pandas as pd
from io import StringIO

DATA = os.path.dirname(os.path.abspath(__file__))
GENES = ["CLDN6", "CCNE1", "MAL", "CTSV", "VTCN1", "MUC16"]

def resolve_ensembl(gene):
    r = requests.get(f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}",
                     headers={"Accept": "application/json"}, timeout=30)
    if r.status_code == 200:
        return r.json().get("id")
    return None

summary_rows = []

for gene in GENES:
    print(f"\n===== {gene} =====")
    ensg = resolve_ensembl(gene)
    print(f"Ensembl ID: {ensg}")
    if not ensg:
        continue
    r = requests.get(f"https://www.proteinatlas.org/{ensg}.tsv", timeout=30)
    try:
        text = r.content.decode("utf-8")
    except UnicodeDecodeError:
        import gzip
        text = gzip.decompress(r.content).decode("utf-8")
    df = pd.read_csv(StringIO(text), sep="\t", low_memory=False)
    print(f"Columns: {list(df.columns)[:12]} ... ({len(df.columns)} total, {len(df)} rows)")
    df.to_csv(os.path.join(DATA, f"hpa_{gene}.tsv"), sep="\t", index=False)
