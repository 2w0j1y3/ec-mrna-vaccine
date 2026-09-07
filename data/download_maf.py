# -*- coding: utf-8 -*-
"""下载 TCGA-UCEC MAF 文件并合并为单一 ucec_maf.tsv。"""
import json
import gzip
import os
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

DATA_DIR = r"C:\Users\wcz\Desktop\EC vaccine\data"
MAF_DIR = os.path.join(DATA_DIR, "maf")
OUT_MANIFEST = os.path.join(DATA_DIR, "ucec_maf_files.json")
OUT_MERGED = os.path.join(DATA_DIR, "ucec_maf.tsv")

os.makedirs(MAF_DIR, exist_ok=True)

KEEP_COLS = [
    "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
    "Variant_Classification", "Variant_Type", "Reference_Allele",
    "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "HGVSp_Short",
    "t_ref_count", "t_alt_count", "Tumor_Sample_Barcode",
]

BASE = "https://api.gdc.cancer.gov/files"
filters = {"op": "and", "content": [
    {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-UCEC"]}},
    {"op": "in", "content": {"field": "files.data_type", "value": ["Masked Somatic Mutation"]}},
    {"op": "in", "content": {"field": "files.data_format", "value": ["MAF"]}},
]}

manifest = []
print("查询 MAF 文件清单 ...")
params = {"filters": json.dumps(filters), "format": "JSON", "size": 1000,
          "fields": "file_id,file_name,file_size,cases.submitter_id"}
r = requests.get(BASE, params=params, timeout=120)
r.raise_for_status()
hits = r.json()["data"]["hits"]
for h in hits:
    sid = h["cases"][0]["submitter_id"] if h.get("cases") else "UNKNOWN"
    manifest.append({"file_id": h["file_id"], "file_name": h["file_name"],
                     "file_size": h.get("file_size"), "submitter_id": sid})
print(f"共 {len(manifest)} 个 MAF 文件")

with open(OUT_MANIFEST, "w", encoding="utf-8") as f:
    json.dump(manifest, f, ensure_ascii=False, indent=2)


def download_one(item):
    fid = item["file_id"]
    sid = item["submitter_id"]
    local = os.path.join(MAF_DIR, f"{sid}.maf.gz")
    if os.path.exists(local) and os.path.getsize(local) > 200:
        return local, sid, "cached"
    url = f"https://api.gdc.cancer.gov/data/{fid}"
    for attempt in range(3):
        try:
            resp = requests.get(url, timeout=120)
            resp.raise_for_status()
            with open(local, "wb") as f:
                f.write(resp.content)
            return local, sid, "ok"
        except Exception as e:
            if attempt == 2:
                return None, sid, f"FAIL:{e}"
    return None, sid, "FAIL"


print("并发下载 MAF ...")
ok = 0
fail = []
with ThreadPoolExecutor(max_workers=10) as ex:
    futs = {ex.submit(download_one, it): it for it in manifest}
    for i, fut in enumerate(as_completed(futs), 1):
        local, sid, status = fut.result()
        if status in ("ok", "cached"):
            ok += 1
        else:
            fail.append(sid)
        if i % 50 == 0:
            print(f"  进度 {i}/{len(manifest)}")
print(f"下载完成: {ok}/{len(manifest)} 成功, {len(fail)} 失败")
if fail:
    print("失败样本:", fail[:20])

# 合并
print("合并 MAF ...")
header = None
rows = []
colidx = None
for it in manifest:
    local = os.path.join(MAF_DIR, f"{it['submitter_id']}.maf.gz")
    if not os.path.exists(local):
        continue
    try:
        with gzip.open(local, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if header is None:
                    header = cols
                    colidx = {c: i for i, c in enumerate(header)}
                    missing = [c for c in KEEP_COLS if c not in colidx]
                    if missing:
                        print("警告: MAF 缺少列", missing)
                    continue
                rec = [cols[colidx[c]] if c in colidx and colidx[c] < len(cols) else "" for c in KEEP_COLS]
                rows.append("\t".join(rec) + "\n")
    except Exception as e:
        print(f"解析失败 {it['submitter_id']}: {e}")

with open(OUT_MERGED, "w", encoding="utf-8") as f:
    f.write("\t".join(KEEP_COLS) + "\n")
    f.writelines(rows)

print(f"合并完成: {len(rows)} 条突变记录 -> {OUT_MERGED}")
print(f"文件大小: {os.path.getsize(OUT_MERGED)/1e6:.1f} MB")