"""合并 UCEC MAF 文件，计算基因突变频次。

输出：
- ucec_mut_freq.csv (gene, mut_freq_pct, n_mut_samples)
"""
import os, gzip, glob, re
import pandas as pd

DATA = os.path.dirname(os.path.abspath(__file__))


def main():
    files = sorted(glob.glob(os.path.join(DATA, "maf", "*.maf.gz")))
    print(f"Total MAF files: {len(files)}")

    # Build patient-level mutation set: gene -> set(patient_barcodes)
    gene_patients = {}
    n_processed = 0

    # MAF column 1 (after #version): Hugo_Symbol
    # Column 16 (Tumor_Sample_Barcode): TCGA-XX-XXXX (12 chars)
    for fpath in files:
        try:
            with gzip.open(fpath, "rt", encoding="utf-8", errors="ignore") as fh:
                header = None
                for line in fh:
                    if line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if header is None:
                        # First non-# line is header
                        header = parts
                        # Find indexes
                        try:
                            i_gene = header.index("Hugo_Symbol")
                            i_patient = header.index("Tumor_Sample_Barcode")
                            i_class = header.index("Variant_Classification") if "Variant_Classification" in header else None
                        except ValueError as e:
                            print(f"  Missing column in {fpath}: {e}")
                            break
                        continue
                    if len(parts) <= max(i_gene, i_patient):
                        continue
                    gene = parts[i_gene]
                    if not gene or gene == "Unknown":
                        continue
                    # Skip synonymous / silent mutations
                    if i_class is not None and parts[i_class] in ("Silent", "RNA", "5'UTR", "3'UTR", "5'Flank", "3'Flank", "IGR", "Intron"):
                        continue
                    # Patient barcode: first 12 chars (TCGA-XX-XXXX)
                    tb = parts[i_patient][:12]
                    gene_patients.setdefault(gene, set()).add(tb)
            n_processed += 1
            if n_processed % 50 == 0:
                print(f"  Processed {n_processed}/{len(files)} files")
        except Exception as e:
            print(f"  Error in {fpath}: {e}")
            continue

    n_patients = len(set().union(*gene_patients.values())) if gene_patients else 0
    # Better: count unique patients across all files
    all_patients = set()
    for ps in gene_patients.values():
        all_patients.update(ps)
    n_patients = len(all_patients)
    print(f"\nProcessed {n_processed}/{len(files)} files")
    print(f"Unique patients with mutations: {n_patients}")
    print(f"Unique mutated genes: {len(gene_patients)}")

    # Build frequency table
    rows = []
    for gene, patients in gene_patients.items():
        n = len(patients)
        freq = n / n_patients * 100 if n_patients else 0
        rows.append({"gene": gene, "mut_freq_pct": freq, "n_mut_samples": n})
    df = pd.DataFrame(rows).sort_values("mut_freq_pct", ascending=False)

    out_path = os.path.join(DATA, "ucec_mut_freq.csv")
    df.to_csv(out_path, index=False)
    print(f"Saved {out_path}")
    print(f"Top 20 mutated genes:")
    print(df.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
