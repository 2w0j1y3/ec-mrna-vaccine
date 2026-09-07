# -*- coding: utf-8 -*-
"""汇总 HPA 6 抗原核查结果：EC IHC 染色 + RNA + 正常组织表达。"""
import os
import pandas as pd

DATA = os.path.dirname(os.path.abspath(__file__))
GENES = ["CLDN6", "CCNE1", "MAL", "CTSV", "VTCN1", "MUC16"]

frames = []
for gene in GENES:
    df = pd.read_csv(os.path.join(DATA, f"hpa_{gene}.tsv"), sep="\t")
    rec = df.iloc[0].to_dict()
    # 找 endometrial / uterus 相关列
    ec_cols = [c for c in df.columns if "endometr" in c.lower() or "uterus" in c.lower() or "uterine" in c.lower()]
    row = {"Gene": gene}
    for c in ec_cols:
        row[c] = rec.get(c)
    # RNA consensus / tissue specificity
    for c in df.columns:
        lc = c.lower()
        if lc in ("tissue specificity", "tissue enriched (RNA)", "rna tissue specificity",
                  "tissue specificity (RNA)", "rna cancer specificity"):
            row[c] = rec.get(c)
    # 血液/免疫相关（VTCN1 风险）
    for c in df.columns:
        if "immune" in c.lower() and "cell" in c.lower() and len(str(rec.get(c))) < 120:
            row[c] = rec.get(c)
    frames.append(row)

out = pd.DataFrame(frames)
pd.set_option("display.max_columns", None)
pd.set_option("display.width", 250)
print("EC/uterus 相关列:")
print(out.to_string())

out.to_csv(os.path.join(DATA, "hpa_summary.csv"), index=False)
print(f"\nSaved hpa_summary.csv")

# 单独打印每个基因的 EC IHC 详情
print("\n" + "=" * 70)
print("详细：每基因的 endometrial cancer IHC + 正常子宫内膜")
print("=" * 70)
for gene in GENES:
    df = pd.read_csv(os.path.join(DATA, f"hpa_{gene}.tsv"), sep="\t")
    rec = df.iloc[0].to_dict()
    print(f"\n--- {gene} ---")
    for c in df.columns:
        lc = c.lower()
        if ("endometr" in lc or "uterus" in lc or "uterine" in lc
                or lc in ("tissue specificity", "rna tissue specificity")):
            v = str(rec.get(c, ""))
            if v and v != "nan":
                print(f"  {c}: {v[:150]}")
