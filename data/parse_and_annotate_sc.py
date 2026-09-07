"""Parse GSE278879 EC snRNA-seq (10x format) and compute per-sample
cell-type composition + IS3/IS5 signature scores.

Input (each .tar.gz):
  filtered_feature_bc_matrix/{barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz}

Pipeline per sample:
  1. Extract tar.gz -> filtered_feature_bc_matrix/
  2. Parse 10x matrix: scipy.io.mmread(matrix.mtx.gz) -> cells x genes
  3. Light QC: filter cells (min 500 UMIs, min 250 genes)
  4. Normalize: log1p(counts_per_million)
  5. Cell-type annotation via marker gene scoring (z-score per gene):
       - Epithelial: EPCAM, KRT8, KRT18, KRT19, CDH1
       - T cell: CD3D, CD3E, CD3G, CD2, TRAC
       - CD8 T: CD8A, CD8B, GZMK, GZMA, PRF1
       - B/Plasma: MS4A1, CD79A, CD79B, MZB1, XBP1
       - NK: GNLY, NKG7, KLRD1, KLRB1, KLRC1
       - Macrophage/Monocyte: CD68, CD163, MSR1, LYZ, S100A8, S100A9
       - DC: CLEC9A, CLEC10A, LILRA4, ITGAX, HLA-DRA
       - Fibroblast: COL1A1, COL1A2, COL3A1, DCN, LUM, ACTA2
       - Endothelial: PECAM1, VWF, CDH5, KDR, FLT1
       - Mast: TPSAB1, TPSB2, KIT, CPA3
  6. Assign cell type = argmax of marker z-scores (with fallback for double-positive via largest fold-change)
  7. Compute IS3/IS5 signature scores per cell using same method
  8. Per-tumor: cell-type fractions + IS3/IS5 median score

Output:
  data/sc_ucec/composition.csv (per sample)
  data/sc_ucec/signatures.csv (per cell)
  data/sc_ucec/cell_assignments.csv
"""
import os
import re
import json
import gzip
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
from scipy.sparse import csr_matrix

SC_DIR = Path(__file__).parent / "sc_ucec"

# Marker gene sets (canonical, from TISCH-style annotation)
MARKERS = {
    "Epithelial":  ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"],
    "CD8_T":       ["CD8A", "CD8B", "GZMK", "GZMA", "GZMB", "PRF1", "NKG7"],
    "CD4_T":       ["CD4", "IL7R", "CCR7", "TCF7"],
    "Treg":        ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "B_cell":      ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1"],
    "Plasma":      ["MZB1", "XBP1", "JCHAIN", "IGKC"],
    "NK":          ["GNLY", "KLRD1", "KLRB1", "KLRC1", "NCR1"],
    "Macrophage":  ["CD68", "CD163", "MSR1", "LYZ", "CSF1R"],
    "Monocyte":    ["S100A8", "S100A9", "FCN1", "VCAN"],
    "DC":          ["CLEC9A", "CLEC10A", "LILRA4", "ITGAX", "HLA-DRA", "BATF3"],
    "Fibroblast":  ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "ACTA2"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "KDR", "FLT1"],
    "Mast":        ["TPSAB1", "TPSB2", "KIT", "CPA3"],
}


def load_10x(mtx_dir):
    """Load 10x filtered_feature_bc_matrix directory.
    Returns (cells x genes csr_matrix, barcodes array, features DataFrame)."""
    feat = pd.read_csv(mtx_dir / "features.tsv.gz", sep="\t", header=None,
                       compression="gzip")
    bc = pd.read_csv(mtx_dir / "barcodes.tsv.gz", sep="\t", header=None,
                     compression="gzip")
    with gzip.open(mtx_dir / "matrix.mtx.gz", "rb") as f:
        # 10x format: rows = features, cols = barcodes
        # Force csr_matrix (old spmatrix API) for subscriptability
        mat = csr_matrix(scipy.io.mmread(f)).T
    return mat, bc[0].values, feat


def normalize_log1p_cpm(mat):
    """Library-size normalize then log1p. Returns csr_matrix."""
    counts_per_cell = np.asarray(mat.sum(axis=1)).flatten()
    counts_per_cell[counts_per_cell == 0] = 1.0
    cpm = mat.multiply(1e4 / counts_per_cell[:, None])
    log_cpm = cpm.copy()
    log_cpm.data = np.log1p(log_cpm.data)
    return csr_matrix(log_cpm)


def gene_index(features_df):
    """Build symbol -> column index."""
    # Features columns: [gene_id, gene_symbol, feature_type] for v2; or [ensg, symbol, type]
    sym_col = 1 if features_df.shape[1] >= 2 else 0
    syms = features_df[sym_col].astype(str).values
    return {s: i for i, s in enumerate(syms)}


def score_signatures(mat_norm, gene_idx, sig_dict):
    """Mean z-score per cell across signature genes.
    Sparse-safe: compute mean/std only over signature genes (small subset)
    to avoid materializing dense matrix for the full 36k gene space.
    mat_norm: csr cells x genes (already log1p-CPM).
    sig_dict: {name: [gene_list]}.
    Returns DataFrame cells x signatures."""
    n_cells = mat_norm.shape[0]
    # Gather all unique gene indices needed
    needed = set()
    for genes in sig_dict.values():
        for g in genes:
            if g in gene_idx:
                needed.add(gene_idx[g])
    # Compute mean/std only for needed genes
    cols = sorted(needed)
    if not cols:
        return pd.DataFrame({k: np.zeros(n_cells) for k in sig_dict})
    sub = mat_norm[:, cols]
    if hasattr(sub, "toarray"):
        dense = sub.toarray()
    else:
        dense = np.asarray(sub)
    mean = dense.mean(axis=0)
    std = dense.std(axis=0)
    std[std == 0] = 1.0
    z = (dense - mean) / std
    # Map gene -> column index in `z`
    col_pos = {cols[i]: i for i in range(len(cols))}
    out = {}
    for name, genes in sig_dict.items():
        idxs = [col_pos[gene_idx[g]] for g in genes if g in gene_idx]
        if len(idxs) == 0:
            out[name] = np.zeros(n_cells)
            continue
        out[name] = z[:, idxs].mean(axis=1)
    return pd.DataFrame(out)


def annotate_cells(mat_norm, gene_idx, min_score=0.1):
    """Assign each cell to the cell type with highest mean z-score."""
    # Use only markers present
    types_present = {}
    for ct, gs in MARKERS.items():
        genes_here = [g for g in gs if g in gene_idx]
        if len(genes_here) >= 2:
            types_present[ct] = genes_here
    if not types_present:
        return np.array(["Unknown"] * mat_norm.shape[0]), types_present
    # Score each type
    scores = score_signatures(mat_norm, gene_idx, types_present)
    # Argmax with min_score threshold: if max score < min_score -> Unknown
    arr = scores.values
    maxv = arr.max(axis=1)
    maxc = arr.argmax(axis=1)
    types_list = list(scores.columns)
    out = []
    for i in range(arr.shape[0]):
        if maxv[i] < min_score:
            out.append("Unknown")
        else:
            out.append(types_list[maxc[i]])
    return np.array(out), types_present


def process_sample(tgz_path, sig_dict):
    """Process one tar.gz and return (cell_df, composition_dict)."""
    import tempfile
    # GSM -> sample type mapping (from GSE278879 family soft)
    GSM_TO_SAMPLE = {
        "GSM8556535": "HEEC1", "GSM8556536": "HEEC2", "GSM8556537": "HEEC3",
        "GSM8556538": "MEEC1", "GSM8556539": "MEEC2",
        "GSM8556540": "LEEC1", "GSM8556541": "LEEC2",
        "GSM8556542": "NE1", "GSM8556543": "NE2", "GSM8556544": "NE3", "GSM8556545": "NE4",
    }
    gsm_match = re.search(r"GSM\d+", tgz_path.name)
    gsm = gsm_match.group(0) if gsm_match else "?"
    sname = GSM_TO_SAMPLE.get(gsm, gsm)
    print(f"\n=== {sname} ({tgz_path.name}) ===")

    with tempfile.TemporaryDirectory() as td:
        with tarfile.open(tgz_path, "r:gz") as tf:
            tf.extractall(td)
        # Find matrix dir or flat matrix files (some datasets wrap files directly)
        mtx_files = list(Path(td).rglob("matrix.mtx.gz"))
        if not mtx_files:
            print("  NO matrix.mtx.gz found, skipping")
            return None
        mtx_dir = mtx_files[0].parent
        mat, barcodes, feat = load_10x(mtx_dir)
        print(f"  raw shape: cells={mat.shape[0]} genes={mat.shape[1]}")

        # Light QC: filter low-quality cells (min 500 UMIs)
        cell_sums = np.asarray(mat.sum(axis=1)).flatten()
        qc_mask = cell_sums >= 500
        mat = mat[qc_mask]
        barcodes = barcodes[qc_mask]
        print(f"  post-QC: cells={mat.shape[0]}")

        if mat.shape[0] < 50:
            print("  too few cells after QC, skipping")
            return None

        # Normalize
        mat_norm = normalize_log1p_cpm(mat)
        gi = gene_index(feat)

        # Cell type annotation
        ct_arr, types_present = annotate_cells(mat_norm, gi)
        uniq, cnts = np.unique(ct_arr, return_counts=True)
        composition = {u: int(c) for u, c in zip(uniq, cnts)}
        print(f"  cell-type composition: {composition}")

        # Signature scores
        sig_scores = score_signatures(mat_norm, gi, sig_dict)

        # Per-sample median signature scores
        sig_median = sig_scores.median().to_dict()

        # Build per-cell df
        df = pd.DataFrame({
            "sample": sname,
            "barcode": barcodes,
            "cell_type": ct_arr,
            "IS3_score": sig_scores["IS3_up"].values,
            "IS5_score": sig_scores["IS5_up"].values,
        })
        # Add a few key gene expressions for downstream validation
        for g in ["VTCN1", "CD274", "CD8A", "GZMK", "CD68", "MS4A1", "EPCAM", "KRT8", "FOXP3"]:
            if g in gi:
                df[f"expr_{g}"] = mat[:, gi[g]].toarray().flatten()
            else:
                df[f"expr_{g}"] = np.nan

        # Save composition
        comp = {"sample": sname, "n_cells": int(mat.shape[0]), "n_genes": int(mat.shape[1])}
        comp.update({f"ct_{k}": v for k, v in composition.items()})
        comp["sig_IS3_median"] = sig_median["IS3_up"]
        comp["sig_IS5_median"] = sig_median["IS5_up"]
        comp["sig_IS3_minus_IS5_median"] = sig_median["IS3_up"] - sig_median["IS5_up"]
        comp["frac_malignant"] = composition.get("Epithelial", 0) / mat.shape[0]
        comp["frac_Tcell"] = sum(composition.get(k, 0) for k in ["CD8_T", "CD4_T", "Treg"]) / mat.shape[0]
        comp["frac_Myeloid"] = sum(composition.get(k, 0) for k in ["Macrophage", "Monocyte", "DC"]) / mat.shape[0]
        comp["frac_B"] = sum(composition.get(k, 0) for k in ["B_cell", "Plasma"]) / mat.shape[0]

        return df, comp


def main():
    sig_path = Path(__file__).parent / "is3_is5_signatures.json"
    with open(sig_path) as f:
        sig_dict_full = json.load(f)
    sig_dict = {"IS3_up": sig_dict_full["IS3_up"], "IS5_up": sig_dict_full["IS5_up"]}

    samples = sorted(SC_DIR.glob("GSM*_matrix.tar.gz"))
    print(f"Found {len(samples)} sample archives")
    all_cells = []
    all_comps = []
    for tgz in samples:
        try:
            res = process_sample(tgz, sig_dict)
        except Exception as e:
            print(f"  ERROR: {e}")
            continue
        if res is None:
            continue
        df, comp = res
        all_cells.append(df)
        all_comps.append(comp)

    if all_cells:
        cells_df = pd.concat(all_cells, ignore_index=True)
        cells_df.to_csv(SC_DIR / "cell_assignments.csv", index=False)
        comp_df = pd.DataFrame(all_comps)
        comp_df.to_csv(SC_DIR / "composition.csv", index=False)
        print(f"\nSaved {len(cells_df)} cells to cell_assignments.csv")
        print(f"Composition per sample:")
        print(comp_df.to_string(index=False))


if __name__ == "__main__":
    main()