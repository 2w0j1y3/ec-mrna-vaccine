# -*- coding: utf-8 -*-
"""CPTAC UCEC 验证（via cBioPortal API）：
1. 临床数据（OS/PFS + 分期/分级）
2. 6 抗原 mRNA（RSEM）+ 蛋白丰度（log2 ratio）
3. 验证：VTCN1 蛋白预后、6 抗原 mRNA-蛋白一致性
"""
import os, json, requests
import pandas as pd

DATA = os.path.dirname(os.path.abspath(__file__))
BASE = "https://www.cbioportal.org/api"
HDRS = {"Accept": "application/json"}
STUDY = "ucec_cptac_2020"
GENES = ["CLDN6", "CCNE1", "MAL", "CTSV", "VTCN1", "MUC16"]

def get(url, params=None):
    r = requests.get(url, headers=HDRS, params=params, timeout=60)
    r.raise_for_status()
    return r.json()

# ---------- 1. 样本列表 ----------
samples = get(f"{BASE}/studies/{STUDY}/samples")
sample_ids = [s["sampleId"] for s in samples]
print(f"CPTAC UCEC samples: {len(sample_ids)}")

# ---------- 2. 临床数据（批量端点） ----------
print("Fetching clinical data ...")
clin_rows = []
for dtype in ["PATIENT", "SAMPLE"]:
    data = get(f"{BASE}/studies/{STUDY}/clinical-data",
               params={"clinicalDataType": dtype, "projection": "DETAILED"})
    recs = {}
    for item in data:
        key = item.get("patientId") or item.get("sampleId")
        recs.setdefault(key, {"patientId": item.get("patientId"),
                              "sampleId": item.get("sampleId")})[item["clinicalAttributeId"]] = item["value"]
    clin_rows.extend(recs.values())
clin = pd.DataFrame(clin_rows)
# 合并 patient-level 和 sample-level（按 patientId）
if "sampleId" in clin.columns:
    pat = clin[clin["sampleId"].isna() | (clin["sampleId"] == clin["patientId"])]
    clin = clin.drop_duplicates(subset=["patientId"])
print(f"Clinical rows: {len(clin)}, columns: {list(clin.columns)[:30]}")
clin.to_csv(os.path.join(DATA, "cptac_clinical.csv"), index=False)

# 生存字段
surv_cols = [c for c in clin.columns if "SURVIVAL" in c.upper() or "OS_" in c.upper() or "STATUS" in c.upper() or "DFS" in c.upper() or "PFS" in c.upper()]
print(f"\nSurvival columns: {surv_cols}")

# ---------- 3. mRNA + 蛋白表达（6 抗原 + 免疫基因） ----------
print("\nFetching mRNA expression ...")
molecular_profile = f"{STUDY}_mrna"
gene_list = GENES + ["CD274", "PDCD1", "CTLA4", "TIGIT", "LAG3", "HAVCR2", "CD276",
                     "CD8A", "GZMB", "IFNG", "STAT1", "CXCL9", "FOXP3",
                     "BATF3", "CLEC9A", "XCR1", "IRF8", "TGFB1", "CD68"]
# symbol → entrez ID
entrez_ids = {}
sym_errors = []
for g in gene_list:
    try:
        gd = get(f"{BASE}/genes/{g}")
        entrez_ids[g] = gd["entrezGeneId"]
    except Exception:
        sym_errors.append(g)
print(f"Resolved {len(entrez_ids)}/{len(gene_list)} entrez IDs; missing: {sym_errors}")

payload = {"sampleListId": f"{STUDY}_all", "entrezGeneIds": list(entrez_ids.values())}
url = f"{BASE}/molecular-profiles/{molecular_profile}/molecular-data/fetch"
r = requests.post(url, headers={**HDRS, "Content-Type": "application/json"},
                  json=payload, timeout=120)
r.raise_for_status()
mrna = pd.DataFrame(r.json())
mrna.to_csv(os.path.join(DATA, "cptac_mrna.csv"), index=False)
print(f"mRNA rows: {len(mrna)}, genes: {mrna['entrezGeneId'].nunique()}")

print("\nFetching protein quantification ...")
prot_profile = f"{STUDY}_protein_quantification"
url = f"{BASE}/molecular-profiles/{prot_profile}/molecular-data/fetch"
r = requests.post(url, headers={**HDRS, "Content-Type": "application/json"},
                  json=payload, timeout=120)
r.raise_for_status()
prot = pd.DataFrame(r.json())
prot.to_csv(os.path.join(DATA, "cptac_protein.csv"), index=False)
print(f"Protein rows: {len(prot)}, genes: {prot['entrezGeneId'].nunique()}")

# ---------- 4. mRNA-蛋白一致性 ----------
print("\n=== 6 抗原 mRNA-蛋白一致性（Spearman）===")
from scipy import stats
m_piv = mrna.pivot_table(index="hugoGeneSymbol", columns="sampleId", values="value")
p_piv = prot.pivot_table(index="hugoGeneSymbol", columns="sampleId", values="value")
for g in GENES:
    if g in m_piv.index and g in p_piv.index:
        common = [c for c in m_piv.columns if c in p_piv.columns]
        m_v = m_piv.loc[g, common].astype(float)
        p_v = p_piv.loc[g, common].astype(float)
        if len(common) > 10:
            r_, p_ = stats.spearmanr(m_v, p_v)
            print(f"  {g:8s}: rho={r_:+.3f}, P={p_:.2e} (n={len(common)})")
        else:
            print(f"  {g:8s}: common samples insufficient ({len(common)})")
    else:
        in_m = g in m_piv.index
        in_p = g in p_piv.index
        print(f"  {g:8s}: mRNA={'Y' if in_m else 'N'}, protein={'Y' if in_p else 'N'}")
