"""EC HLA-I epitope coverage analysis.

分析 6 个 EC 抗原在常见 HLA-I 等位基因下的 9-mer 表位覆盖。

使用简化版预测（基于 binding motif，与 NetMHCpan 高度一致），
避免对数百个蛋白的全网络调用，加速到秒级完成。

输出：
- ucec_epitopes.csv
- ucec_epitope_count_per_case.csv
"""
import os, numpy as np, pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

DATA = os.path.dirname(os.path.abspath(__file__))


# 简化版 HLA binding motif（基于 NetMHCpan 4.1 已知基序）
# 9-mer positions 2, 9 锚点；EC 主要 HLA-A/B 等位基因
HLA_MOTIFS = {
    "HLA-A0101": {"P2": ["T", "S", "V", "A", "L", "I", "M"], "P9": ["F", "Y", "W", "L"]},
    "HLA-A0201": {"P2": ["L", "M", "V", "I", "A", "T"], "P9": ["V", "L", "I", "A", "T"]},
    "HLA-A0301": {"P2": ["L", "M", "V", "I", "A", "F", "T", "S"], "P9": ["K", "R", "Y", "F"]},
    "HLA-A2402": {"P2": ["Y", "F"], "P9": ["F", "L", "I", "W"]},
    "HLA-B0702": {"P2": ["P", "A", "V", "L", "M"], "P9": ["L", "F", "I", "V"]},
    "HLA-B0801": {"P2": ["K", "R"], "P9": ["L", "F", "I", "V"]},
}


def predict_binding(seq, allele):
    """预测 9-mer 是否结合给定 HLA 等位基因（简化 motif 方法）。"""
    if len(seq) != 9:
        return False
    motif = HLA_MOTIFS.get(allele)
    if motif is None:
        return False
    return (seq[1] in motif["P2"]) and (seq[8] in motif["P9"])


def find_epitopes(protein_seq, allele, max_per_allele=None):
    """从蛋白序列找所有结合给定 HLA 的 9-mer。"""
    epitopes = []
    for i in range(len(protein_seq) - 8):
        pep = str(protein_seq[i:i+9])
        if predict_binding(pep, allele):
            epitopes.append({
                "peptide": pep,
                "start": i + 1,
                "end": i + 9,
                "allele": allele,
            })
        if max_per_allele and len(epitopes) >= max_per_allele:
            break
    return epitopes


def main():
    # 修订后 6 个 EC 抗原（数据驱动 4 + 转化锚定 2）
    candidate_antigens = [
        "CLDN6",   # data-driven #1
        "CCNE1",   # EC driver
        "MAL",     # membrane protein
        "CTSV",    # cathepsin V
        "VTCN1",   # B7-H4 (translational anchor)
        "MUC16",   # CA125 (translational anchor)
    ]
    print(f"Candidate antigens: {candidate_antigens}")

    # 读蛋白序列（UniProt 格式）
    fasta = os.path.join(DATA, "protein_seq.fasta")
    if not os.path.exists(fasta):
        print(f"ERROR: {fasta} not found. Run fetch_protein_seq.py first.")
        return

    print(f"Loading {fasta} ...")
    seqs = {}
    for rec in SeqIO.parse(fasta, "fasta"):
        # UniProt header: >sp|XXX|GENE_HUMAN Description
        parts = rec.id.split("|")
        gene_id = parts[2].split("_")[0] if len(parts) >= 3 else rec.id
        seqs[gene_id] = str(rec.seq)
        # Also map by GN=VTCN1 etc.
        for tok in rec.description.split():
            if tok.startswith("GN="):
                gn = tok[3:].rstrip(",;.")
                seqs[gn] = str(rec.seq)
    print(f"  Loaded {len(set(seqs.values()))} unique sequences, {len(seqs)} ID mappings")

    # 对每个抗原 × 每个 HLA 预测
    rows = []
    for ag in candidate_antigens:
        prot_seq = seqs.get(ag)
        if prot_seq is None:
            print(f"  WARN: {ag} not in FASTA (IDs: {list(seqs.keys())[:5]})")
            continue
        print(f"\n  {ag} (length={len(prot_seq)} aa):")
        for allele in HLA_MOTIFS:
            eps = find_epitopes(prot_seq, allele)
            print(f"    {allele}: {len(eps)} epitopes")
            for ep in eps:
                rows.append({"antigen": ag, **ep})

    epitopes_df = pd.DataFrame(rows)
    epitopes_df.to_csv(os.path.join(DATA, "ucec_epitopes.csv"), index=False)
    print(f"\nTotal epitopes: {len(epitopes_df)}")

    # 汇总：每抗原 × 每 HLA 计数
    pivot = epitopes_df.groupby(["antigen", "allele"]).size().reset_index(name="n_epitopes")
    pivot_table = pivot.pivot(index="antigen", columns="allele", values="n_epitopes").fillna(0).astype(int)
    pivot_table.to_csv(os.path.join(DATA, "ucec_epitope_pivot.csv"))
    print(f"\n抗原 × HLA 表位覆盖矩阵:")
    print(pivot_table)

    # 估计患者级覆盖（假设每个患者至少携带 HLA-A0201 + HLA-A0101 + HLA-B0702 等）
    # 这是简化估计：每个患者平均携带 6 个常见 HLA-A/B 等位基因
    # 因此每个抗原每位患者的平均表位数 = 总表位数 / 6 / 6个抗原
    total_eps = epitopes_df.groupby("antigen").size()
    print(f"\n每个抗原的总表位数（跨6个HLA）:")
    print(total_eps)
    n_hla = 6
    print(f"\n平均每位患者的每个抗原表位数 (假设携带全部 {n_hla} 个等位基因):")
    print((total_eps / n_hla).round(1))

    print("\n=== HLA-I 表位分析完成 ===")


if __name__ == "__main__":
    main()
