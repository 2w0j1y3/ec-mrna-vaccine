# -*- coding: utf-8 -*-
"""GDC UCEC 表达数据查询 - 获取 TCGA-UCEC (endometrial) RNA-seq 文件清单。

字段：file_id, file_name, file_size, case submitter_id, sample submitter_id,
      sample_type, aliquot submitter_id

保存到 EC vaccine/data/ucec_expr_files.json
"""
import os
import json
import requests
from collections import Counter

OUT = r"C:\Users\wcz\Desktop\EC vaccine\data\ucec_expr_files.json"

BASE = "https://api.gdc.cancer.gov/files"
filters = {
    "op": "and",
    "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-UCEC"]}},
        {"op": "in", "content": {"field": "files.experimental_strategy", "value": ["RNA-Seq"]}},
        {"op": "in", "content": {"field": "files.data_type", "value": ["Gene Expression Quantification"]}},
        {"op": "in", "content": {"field": "files.analysis.workflow_type", "value": ["STAR - Counts"]}},
        {"op": "in", "content": {"field": "files.access", "value": ["open"]}},
    ],
}
fields = "file_id,file_name,file_size,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type,cases.samples.portions.analytes.aliquots.submitter_id"

all_hits = []
params = {"filters": json.dumps(filters), "format": "JSON", "size": 100, "fields": fields}
while True:
    r = requests.get(BASE, params=params, timeout=120)
    r.raise_for_status()
    j = r.json()
    hits = j["data"]["hits"]
    all_hits.extend(hits)
    pagination = j["data"]["pagination"]
    total = pagination["total"]
    params["from"] = len(all_hits)
    if len(all_hits) >= total:
        break

print(f"查询到 {len(all_hits)} 个文件（总 {total}）")

# 解析
records = []
for h in all_hits:
    case_sid = h["cases"][0]["submitter_id"] if h.get("cases") else None
    samples = h["cases"][0].get("samples", []) if h.get("cases") else []
    sample_sid = samples[0].get("submitter_id") if samples else None
    sample_type = samples[0].get("sample_type") if samples else None
    aliquot_sid = None
    if samples:
        portions = samples[0].get("portions", [])
        if portions:
            analytes = portions[0].get("analytes", [])
            if analytes:
                aliquots = analytes[0].get("aliquots", [])
                if aliquots:
                    aliquot_sid = aliquots[0].get("submitter_id")
    records.append({
        "file_id": h["file_id"],
        "file_name": h["file_name"],
        "file_size": h.get("file_size"),
        "case_id": case_sid,
        "sample_id": sample_sid,
        "aliquot_id": aliquot_sid,
        "sample_type": sample_type,
    })

with open(OUT, "w", encoding="utf-8") as f:
    json.dump(records, f, ensure_ascii=False, indent=2)

# 统计
print("sample_type 分布:", Counter(r["sample_type"] for r in records))
print("唯一 case_id:", len(set(r["case_id"] for r in records)))
print("唯一 sample_id:", len(set(r["sample_id"] for r in records)))
print("唯一 aliquot_id:", len(set(r["aliquot_id"] for r in records)))
print(f"已保存 {OUT}")