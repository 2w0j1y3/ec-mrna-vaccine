"""合并 UCEC 表达 TSV 文件（numpy 版，内存最优）。

输出：
- ucec_tpm_matrix.csv (genes × samples) - 内存允许时输出
- ucec_tpm_matrix.npz (numpy 压缩) - 主输出
- ucec_sample_map.csv
"""
import os, glob
import pandas as pd
import numpy as np

DATA = os.path.dirname(os.path.abspath(__file__))


def parse_case_from_filename(fname):
    base = os.path.basename(fname).replace(".tpm.tsv", "")
    parts = base.split("-")
    if len(parts) < 4:
        return None, None, None
    case_id = "-".join(parts[:3])
    sample_part = parts[3]
    type_code = sample_part[:2]
    sample_type = {
        "01": "Primary Tumor",
        "02": "Recurrent",
        "07": "Recurrent",
        "11": "Solid Tissue Normal",
        "14": "Bone Marrow Normal",
    }.get(type_code, "Unknown")
    return case_id, sample_type, base


def main():
    expr_dir = os.path.join(DATA, "expr")
    tsv_files = sorted(glob.glob(os.path.join(expr_dir, "*.tpm.tsv")))
    print(f"Found {len(tsv_files)} files", flush=True)

    # Parse all files: collect per-sample {gene: value} dicts
    sample_data = {}
    sample_info = []

    for i, fpath in enumerate(tsv_files):
        case_id, sample_type, sample_id = parse_case_from_filename(fpath)
        if case_id is None:
            continue
        try:
            df = pd.read_csv(fpath, sep="\t", skiprows=1, header=None,
                             names=["gene_id", "gene_name", "value"], usecols=[1, 2])
        except Exception as e:
            print(f"  ERR {fpath}: {e}", flush=True)
            continue
        gene_expr = dict(zip(df["gene_name"].astype(str), df["value"].astype(float)))
        sample_data[sample_id] = gene_expr
        sample_info.append({"sample_id": sample_id, "case_id": case_id, "sample_type": sample_type})
        if (i + 1) % 100 == 0:
            print(f"  Read {i+1}/{len(tsv_files)}", flush=True)

    print(f"Loaded {len(sample_data)} samples", flush=True)

    all_samples = sorted(sample_data.keys())

    # Find unique genes
    all_genes_set = set()
    for sid, gd in sample_data.items():
        all_genes_set.update(gd.keys())
    all_genes = sorted(all_genes_set)
    n_genes = len(all_genes)
    n_samples = len(all_samples)
    print(f"Genes: {n_genes}, Samples: {n_samples}", flush=True)

    # Build numpy matrix directly
    gene_to_idx = {g: i for i, g in enumerate(all_genes)}
    matrix = np.zeros((n_genes, n_samples), dtype=np.float32)

    for j, sid in enumerate(all_samples):
        for g, v in sample_data[sid].items():
            i = gene_to_idx[g]
            matrix[i, j] = v
        if (j + 1) % 50 == 0:
            print(f"  Filled {j+1}/{n_samples}", flush=True)

    print(f"Matrix shape: {matrix.shape}", flush=True)
    print(f"Matrix dtype: {matrix.dtype}, size: {matrix.nbytes/1e6:.1f} MB", flush=True)

    # Save as npz (compressed, fast)
    np.savez_compressed(os.path.join(DATA, "ucec_tpm_matrix.npz"),
                        matrix=matrix, genes=np.array(all_genes), samples=np.array(all_samples))
    print(f"Saved ucec_tpm_matrix.npz", flush=True)

    # Also save as CSV for compatibility (smaller, gene subset to save time)
    # Build DataFrame only for protein-coding genes (use saved file)
    pc_genes_path = os.path.join(DATA, "protein_coding_genes.txt")
    if os.path.exists(pc_genes_path):
        pc_genes = set(open(pc_genes_path).read().split("\n"))
        keep_idx = [i for i, g in enumerate(all_genes) if g in pc_genes]
        keep_genes = [all_genes[i] for i in keep_idx]
        sub_matrix = matrix[keep_idx, :]
        df = pd.DataFrame(sub_matrix, index=keep_genes, columns=all_samples)
        df.index.name = "gene"
        df.to_csv(os.path.join(DATA, "ucec_tpm_matrix.csv"))
        print(f"Saved ucec_tpm_matrix.csv ({len(keep_genes)} genes × {n_samples} samples)", flush=True)
    else:
        print(f"protein_coding_genes.txt not found, skipping CSV output", flush=True)

    # Sample map
    smap = pd.DataFrame(sample_info).drop_duplicates(subset="sample_id")
    smap.to_csv(os.path.join(DATA, "ucec_sample_map.csv"), index=False)
    print(f"\nSample types:")
    print(smap["sample_type"].value_counts())
    print("\n=== 合并完成 ===")


if __name__ == "__main__":
    main()
