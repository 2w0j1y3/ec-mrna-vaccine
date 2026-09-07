"""Fetch protein sequences for the REVISED 6 EC antigens + TP53 from UniProt.

Output: protein_seq.fasta (overwrites)
"""
import os, urllib.request

DATA = os.path.dirname(os.path.abspath(__file__))

# 修订后 6 抗原（数据驱动 4 + 转化锚定 2）
ANTIGENS = [
    "CLDN6",   # Claudin 6 (data-driven #1)
    "CCNE1",   # Cyclin E1 (EC driver)
    "MAL",     # MAL (membrane protein)
    "CTSV",    # Cathepsin V
    "VTCN1",   # B7-H4 (translational anchor, IS5)
    "MUC16",   # CA125 (translational anchor)
]


def fetch_uniprot_fasta(gene):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene}+AND+organism_id:9606+AND+reviewed:true&format=fasta&size=1&fields=accession,sequence"
    req = urllib.request.Request(url, headers={"User-Agent": "EC-vaccine-research"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            return resp.read().decode("utf-8")
    except Exception as e:
        print(f"  UniProt search failed for {gene}: {e}")
        return None


def main():
    out_fasta = os.path.join(DATA, "protein_seq.fasta")
    if os.path.exists(out_fasta):
        os.remove(out_fasta)

    with open(out_fasta, "w") as out:
        for ag in ANTIGENS:
            print(f"Fetching {ag} from UniProt...")
            fasta = fetch_uniprot_fasta(ag)
            if fasta:
                lines = fasta.strip().split("\n")
                out.write(lines[0] + "\n")
                for line in lines[1:]:
                    if line.startswith(">"):
                        break
                    out.write(line + "\n")
                print(f"  OK")
            else:
                print(f"  FAILED for {ag}")

    # Verify
    with open(out_fasta) as f:
        content = f.read()
    n = content.count(">")
    genes = [l.split("GN=")[1].split()[0] for l in content.split("\n") if "GN=" in l]
    print(f"\n{out_fasta}: {n} sequences, genes: {genes}")


if __name__ == "__main__":
    main()
