import requests, json, time, os
from concurrent.futures import ThreadPoolExecutor, as_completed

files = json.load(open('ucec_expr_files.json'))
existing = os.listdir('expr/') if os.path.exists('expr/') else []
print(f'Expression: {len(existing)}/{len(files)} downloaded')

maf_files = json.load(open('ucec_maf_files.json'))
maf_existing = os.listdir('maf/') if os.path.exists('maf/') else []
print(f'MAF: {len(maf_existing)}/{len(maf_files)} downloaded')

# Show last 5 downloaded
expr_ids = sorted([f for f in existing if f.endswith('.tsv')])
print(f'Last 5 expression files:')
for f in expr_ids[-5:]:
    sz = os.path.getsize(f'expr/{f}') / 1e6
    print(f'  {f}  {sz:.1f}MB')

print(f'Last 5 MAF files:')
for f in sorted(maf_existing)[-5:]:
    sz = os.path.getsize(f'maf/{f}') / 1e6
    print(f'  {f}  {sz:.1f}MB')
