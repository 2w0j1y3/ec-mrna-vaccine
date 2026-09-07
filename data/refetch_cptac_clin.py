# -*- coding: utf-8 -*-
"""重新抓取 CPTAC UCEC 临床数据（patient 层 + sample 层分开）。"""
import os
import requests
import pandas as pd

DATA = os.path.dirname(os.path.abspath(__file__))
BASE = "https://www.cbioportal.org/api"
HDRS = {"Accept": "application/json"}
STUDY = "ucec_cptac_2020"

def get(url, params=None):
    r = requests.get(url, headers=HDRS, params=params, timeout=60)
    r.raise_for_status()
    return r.json()

# PATIENT 层
data = get(f"{BASE}/studies/{STUDY}/clinical-data", params={"clinicalDataType": "PATIENT"})
rows = {}
for item in data:
    pid = item["patientId"]
    rows.setdefault(pid, {"patientId": pid})[item["clinicalAttributeId"]] = item["value"]
pat = pd.DataFrame(rows.values())
print(f"PATIENT rows: {len(pat)}, cols: {len(pat.columns)}")

# SAMPLE 层
data = get(f"{BASE}/studies/{STUDY}/clinical-data", params={"clinicalDataType": "SAMPLE"})
rows = {}
for item in data:
    sid = item["sampleId"]
    rows.setdefault(sid, {"sampleId": sid, "patientId": item["patientId"]})[item["clinicalAttributeId"]] = item["value"]
samp = pd.DataFrame(rows.values())
print(f"SAMPLE rows: {len(samp)}, cols: {len(samp.columns)}")

# 合并
merged = samp.merge(pat, on="patientId", how="left", suffixes=("_samp", ""))
print(f"Merged: {merged.shape}")
merged.to_csv(os.path.join(DATA, "cptac_clinical.csv"), index=False)

for c in ["POLE_SUBTYPE", "MSI_STATUS", "TMB_NONSYNONYMOUS", "ESTIMATE_IMMUNESCORE", "CIBERSORT_T _CELLS _CD8"]:
    if c in merged.columns:
        print(f"  {c}: non-null {merged[c].notna().sum()}, e.g. {merged[c].dropna().unique()[:4]}")
    else:
        # 可能在 patient 层
        found = [x for x in pat.columns if c.split("_")[0] in x]
        print(f"  {c}: NOT in merged (patient-layer candidates: {found[:3]})")
