# -*- coding: utf-8 -*-
"""从 MAF 直接推导分子分型标记：
1. POLE 超突变（外切酶结构域热点突变：P286R, V411L, S297F, A456P, S459F 等）
2. MSI 评分（微卫星区 indel 负荷 — 简化 MANTIS 风格）
3. 与免疫亚型交叉验证（修复 clinical 字段全空的问题）
"""
import os, gzip, glob
import numpy as np
import pandas as pd
from scipy import stats

DATA = os.path.dirname(os.path.abspath(__file__))

POLE_HOTSPOTS = {
    "P286R": ("P", 286, "R"), "V411L": ("V", 411, "L"), "S297F": ("S", 297, "F"),
    "A456P": ("A", 456, "P"), "S459F": ("S", 459, "F"), "D368Y": ("D", 368, "Y"),
    "M297R": ("M", 297, "R"), "L424V": ("L", 424, "V"), "C534Y": ("C", 534, "Y"),
    "P436R": ("P", 436, "R"), "G420S": ("G", 420, "S"), "R443Q": ("R", 443, "Q"),
}

files = sorted(glob.glob(os.path.join(DATA, "maf", "*.maf.gz")))
print(f"Processing {len(files)} MAF files ...")

pole_cases = {}       # case -> hotspot mutation
case_indels = {}      # case -> total indels
case_ms_indels = {}   # case -> indels in repeat/homopolymer context
case_snvs = {}

for i, fpath in enumerate(files):
    case = os.path.basename(fpath).replace(".maf.gz", "")
    try:
        with gzip.open(fpath, "rt", encoding="utf-8", errors="ignore") as fh:
            header = None
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if header is None:
                    header = parts
                    i_gene = header.index("Hugo_Symbol")
                    i_class = header.index("Variant_Classification")
                    i_type = header.index("Variant_Type") if "Variant_Type" in header else None
                    i_aa = None
                    for cand in ["Protein_Change", "HGVSp_Short", "HGVSp"]:
                        if cand in header:
                            i_aa = header.index(cand)
                            break
                    i_ref = header.index("Reference_Allele")
                    i_alt = header.index("Tumor_Seq_Allele2") if "Tumor_Seq_Allele2" in header else header.index("Tumor_Seq_Allele1")
                    i_start = header.index("Start_Position")
                    continue
                if len(parts) <= i_alt:
                    continue
                gene = parts[i_gene]
                vclass = parts[i_class]
                vtype = parts[i_type] if i_type is not None else ""
                ref, alt = parts[i_ref], parts[i_alt]

                # SNV 计数（非同义）
                if vtype == "SNP" and vclass not in ("Silent", "RNA", "Intron", "IGR", "3'UTR", "5'UTR", "3'Flank", "5'Flank"):
                    case_snvs[case] = case_snvs.get(case, 0) + 1

                # Indel 计数
                if vtype == "INS" or vtype == "DEL":
                    case_indels[case] = case_indels.get(case, 0) + 1
                    # 简化 MSI 判定：indel 长度是 1-5bp 且位于同聚物/重复上下文
                    # MAF 无上下文序列时，用 indel 长度 <= 2 作为代理（微卫星 indel 多为 1-2bp）
                    try:
                        indel_len = abs(len(ref) - len(alt))
                    except Exception:
                        indel_len = 1
                    if 1 <= indel_len <= 2:
                        case_ms_indels[case] = case_ms_indels.get(case, 0) + 1

                # POLE 热点
                if gene == "POLE" and i_aa is not None and len(parts) > i_aa:
                    aa = parts[i_aa]
                    for hs_name, (wt, pos, mut) in POLE_HOTSPOTS.items():
                        # 匹配 p.P286R / P286R 等格式
                        import re
                        mm = re.search(r"p\.([A-Z])(\d+)([A-Z])", aa)
                        if mm:
                            if int(mm.group(2)) == pos and mm.group(3) == mut:
                                pole_cases[case] = hs_name
    except Exception as e:
        continue
    if (i + 1) % 100 == 0:
        print(f"  {i+1}/{len(files)}", flush=True)

print(f"\nPOLE 热点突变病例: {len(pole_cases)}")
for c, hs in sorted(pole_cases.items()):
    print(f"  {c}: {hs}")

# MSI 评分：微卫星 indel 数 / 总 indel 数，加上绝对 indel 负荷
rows = []
all_cases = set(case_indels.keys()) | set(case_snvs.keys())
for c in all_cases:
    indels = case_indels.get(c, 0)
    ms_indels = case_ms_indels.get(c, 0)
    snvs = case_snvs.get(c, 0)
    ms_frac = ms_indels / indels if indels > 0 else 0
    rows.append({
        "case_id": c,
        "n_snv": snvs,
        "n_indel": indels,
        "n_ms_indel": ms_indels,
        "ms_indel_frac": round(ms_frac, 3),
        "POLE_hotspot": pole_cases.get(c, ""),
    })
msi_df = pd.DataFrame(rows)

# MSI 分类（EC 阈值：MSI-H 通常 >100 个 MS-indel 或 indel/SNV > 0.3）
def classify(r):
    if r["POLE_hotspot"]:
        return "POLE"
    if r["n_indel"] >= 100 and r["ms_indel_frac"] >= 0.5:
        return "MSI-H"
    if r["n_indel"] >= 30 and r["ms_indel_frac"] >= 0.4:
        return "MSI-H"
    return "MSS"
msi_df["molecular_class"] = msi_df.apply(classify, axis=1)

print(f"\n分子分型分布:")
print(msi_df["molecular_class"].value_counts())

msi_df.to_csv(os.path.join(DATA, "ucec_molecular_class.csv"), index=False)

# 与免疫亚型交叉
subtype = pd.read_csv(os.path.join(DATA, "ucec_subtypes.csv"))
smap = pd.read_csv(os.path.join(DATA, "ucec_sample_map.csv"))
sample_to_case = dict(zip(smap["sample_id"], smap["case_id"]))
sub_cases = subtype[subtype["sample_id"].isin(
    smap[smap["sample_type"] == "Primary Tumor"]["sample_id"])]
sub_cases["case_id"] = sub_cases["sample_id"].map(sample_to_case)
merged = sub_cases.merge(msi_df, on="case_id", how="left")
merged["molecular_class"] = merged["molecular_class"].fillna("NA")

print(f"\n免疫亚型 × 分子分型:")
ct = pd.crosstab(merged["subtype"], merged["molecular_class"])
print(ct)
ct.to_csv(os.path.join(DATA, "ucec_subtype_vs_molecular.csv"))

# 卡方检验
try:
    chi2, p, dof, _ = stats.chi2_contingency(ct.values)
    print(f"\n卡方检验: chi2={chi2:.2f}, P={p:.2e}")
except Exception as e:
    print(f"卡方失败: {e}")

# TMB by molecular class
print(f"\n各分子类别 TMB (SNV+indel / 38 Mb):")
merged["TMB_calc"] = (merged["n_snv"] + merged["n_indel"]) / 38.0
print(merged.groupby("molecular_class")["TMB_calc"].agg(["median", "mean", "count"]).round(1))

# 各亚型的 POLE/MSI 比例
print(f"\n各亚型内 POLE+MSI-H 比例:")
for s in sorted(merged["subtype"].unique()):
    sub = merged[merged["subtype"] == s]
    hot = (sub["molecular_class"].isin(["POLE", "MSI-H"])).sum()
    print(f"  {s}: {hot}/{len(sub)} ({hot/len(sub)*100:.1f}%)")
