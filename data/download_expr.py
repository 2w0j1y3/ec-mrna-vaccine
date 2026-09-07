# -*- coding: utf-8 -*-
"""下载 TCGA-UCEC RNA-seq 表达数据 (STAR-Counts TPM)。

从 ucec_expr_files.json 读取 file_id，下载 GDC STAR-Counts 文件，
提取 TPM 列输出到 expr/<sample_id>.tpm.tsv。
"""
import os
import json
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import Counter

DATA = r"C:\Users\wcz\Desktop\EC vaccine\data"
files = json.load(open(os.path.join(DATA, "ucec_expr_files.json")))
os.makedirs(os.path.join(DATA, "expr"), exist_ok=True)

meta = []
for f in files:
    fid = f["file_id"]
    sid = f["sample_id"]
    st = f["sample_type"]
    meta.append((fid, sid, st))

print(f"总文件数: {len(meta)}")


def download_and_extract(item):
    fid, sid, st = item
    out = os.path.join(DATA, "expr", sid + ".tpm.tsv")
    if os.path.exists(out) and os.path.getsize(out) > 1000:
        return (sid, "cached", st)
    url = "https://api.gdc.cancer.gov/data/" + fid
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    try:
        r = urllib.request.urlopen(req, timeout=180)
        raw = r.read().decode("utf-8", "ignore")
    except Exception as e:
        return (sid, "ERROR:" + str(e)[:80], st)
    lines = raw.split("\n")
    header_idx = None
    for i, l in enumerate(lines):
        if l.startswith("gene_id\tgene_name"):
            header_idx = i
            break
    if header_idx is None:
        return (sid, "ERROR:no-header", st)
    hdr = lines[header_idx].split("\t")
    try:
        tpm_idx = hdr.index("tpm_unstranded")
    except ValueError:
        tpm_idx = hdr.index("fpkm_unstranded")
    out_lines = ["gene_id\tgene_name\ttpm"]
    for l in lines[header_idx + 1:]:
        if not l.strip():
            continue
        parts = l.split("\t")
        if len(parts) <= tpm_idx:
            continue
        gid = parts[0]
        gname = parts[1]
        tpm = parts[tpm_idx]
        if gid.startswith("ENSG"):
            out_lines.append(f"{gid}\t{gname}\t{tpm}")
    with open(out, "w") as fo:
        fo.write("\n".join(out_lines))
    return (sid, "ok", st)


results = []
with ThreadPoolExecutor(max_workers=20) as ex:
    futs = {ex.submit(download_and_extract, m): m for m in meta}
    done = 0
    for fut in as_completed(futs):
        sid, status, st = fut.result()
        results.append((sid, status, st))
        done += 1
        if done % 50 == 0:
            print(f"进度 {done}/{len(meta)}")

ok = [r for r in results if r[1] == "ok"]
cached = [r for r in results if r[1] == "cached"]
errs = [r for r in results if r[1].startswith("ERROR")]
print(f"成功: {len(ok)} | 缓存: {len(cached)} | 错误: {len(errs)}")
for e in errs[:10]:
    print("  ERR", e)

# 保存样本类型映射
meta_out = {}
for fid, sid, st in meta:
    meta_out[sid] = st
json.dump(meta_out, open(os.path.join(DATA, "ucec_sample_type.json"), "w"))
print("样本类型映射已保存")
print("Sample types:", Counter(meta_out.values()))