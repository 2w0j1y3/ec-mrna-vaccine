# -*- coding: utf-8 -*-
"""下载 TCGA-UCEC 临床数据（GDC cases API），提取生存/分期/分级/年龄/性别等字段。
"""
import json
import csv
import os
import requests
from collections import Counter

OUT_JSON = r"C:\Users\wcz\Desktop\EC vaccine\data\ucec_clinical.json"
OUT_CSV = r"C:\Users\wcz\Desktop\EC vaccine\data\ucec_clinical.csv"

BASE = "https://api.gdc.cancer.gov/cases"

FIELDS = [
    "submitter_id",
    "demographic.sex_at_birth",
    "demographic.race",
    "demographic.ethnicity",
    "demographic.vital_status",
    "demographic.days_to_death",
    "demographic.age_at_index",
    "demographic.year_of_death",
    "diagnoses.age_at_diagnosis",
    "diagnoses.days_to_last_follow_up",
    "diagnoses.days_to_death",
    "diagnoses.ajcc_pathologic_stage",
    "diagnoses.tumor_stage",
    "diagnoses.tumor_grade",
    "diagnoses.ajcc_pathologic_t",
    "diagnoses.ajcc_pathologic_n",
    "diagnoses.ajcc_pathologic_m",
    "diagnoses.site_of_resection_or_biopsy",
    "diagnoses.primary_diagnosis",
    "diagnoses.morphology",
    "diagnoses.days_to_last_known_disease_status",
    "diagnoses.last_known_disease_status",
    "diagnoses.microsatellite_instability_test_result",  # EC has MSI status!
    "diagnoses.mismatch_repair_protein_result",           # EC has MMR!
    "diagnoses.pole_mutation_status",                     # EC POLE!
]

filters = {
    "op": "in",
    "content": {"field": "project.project_id", "value": ["TCGA-UCEC"]},
}

params = {
    "filters": json.dumps(filters),
    "fields": ",".join(FIELDS),
    "format": "JSON",
    "size": 1000,
}

all_cases = []
print("Fetching GDC UCEC cases ...")
r = requests.get(BASE, params=params, timeout=120)
r.raise_for_status()
data = r.json()
hits = data["data"]["hits"]
print(f"总病例数: {data['data']['pagination']['total']}")

for h in hits:
    rec = {"submitter_id": h.get("submitter_id")}
    demo = h.get("demographic") or {}
    rec["gender"] = demo.get("sex_at_birth") or demo.get("gender")
    rec["race"] = demo.get("race")
    rec["ethnicity"] = demo.get("ethnicity")
    rec["vital_status"] = demo.get("vital_status")
    rec["days_to_death"] = demo.get("days_to_death")
    rec["age_at_index"] = demo.get("age_at_index")
    rec["year_of_death"] = demo.get("year_of_death")
    diag_list = h.get("diagnoses") or []
    if diag_list:
        d = diag_list[0]
        rec["age_at_diagnosis"] = d.get("age_at_diagnosis")
        rec["days_to_last_follow_up"] = d.get("days_to_last_follow_up")
        rec["days_to_death_dx"] = d.get("days_to_death")
        rec["ajcc_pathologic_stage"] = d.get("ajcc_pathologic_stage")
        rec["tumor_stage"] = d.get("tumor_stage")
        rec["tumor_grade"] = d.get("tumor_grade")
        rec["ajcc_pathologic_t"] = d.get("ajcc_pathologic_t")
        rec["ajcc_pathologic_n"] = d.get("ajcc_pathologic_n")
        rec["ajcc_pathologic_m"] = d.get("ajcc_pathologic_m")
        rec["site_of_resection_or_biopsy"] = d.get("site_of_resection_or_biopsy")
        rec["primary_diagnosis"] = d.get("primary_diagnosis")
        rec["morphology"] = d.get("morphology")
        rec["msi_status"] = d.get("microsatellite_instability_test_result")
        rec["mmr_status"] = d.get("mismatch_repair_protein_result")
        rec["pole_mutation"] = d.get("pole_mutation_status")
    all_cases.append(rec)

with open(OUT_JSON, "w", encoding="utf-8") as f:
    json.dump(all_cases, f, ensure_ascii=False, indent=2)
print(f"已保存 {OUT_JSON} ({len(all_cases)} 条)")

if all_cases:
    keys = list(all_cases[0].keys())
    with open(OUT_CSV, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(all_cases)
    print(f"已保存 {OUT_CSV}")

print("vital_status:", Counter(c.get("vital_status") for c in all_cases))
print("gender:", Counter(c.get("gender") for c in all_cases))
print("stage:", Counter(c.get("ajcc_pathologic_stage") for c in all_cases))
print("grade:", Counter(c.get("tumor_grade") for c in all_cases))
print("MSI status:", Counter(c.get("msi_status") for c in all_cases))
print("MMR status:", Counter(c.get("mmr_status") for c in all_cases))
print("POLE mutation:", Counter(c.get("pole_mutation") for c in all_cases))

n_os = 0
for c in all_cases:
    dtd = c.get("days_to_death")
    dtf = c.get("days_to_last_follow_up")
    if dtd or dtf:
        n_os += 1
print(f"可计算 OS 时间病例数: {n_os}")