"""Why did a candidate fail verification? Answer it with numbers.

For every chain in a candidate's structure this prints two different
measures of "is this the target protein":

  exact 8-mer identity   -- what verification currently uses
  aligned identity       -- Smith-Waterman, with coverage

The gap between those two columns is the diagnosis. An exact k-mer match
needs all k positions identical, so at sequence identity p the expected
match rate is p^k: a true ortholog at 90% identity scores 0.9^8 = 43%, and
a different viral strain of the same protein scores far lower still. A
threshold applied to the first column therefore rejects proteins that are
unambiguously the target.

Run: python scripts/diagnose_catalog.py [KEY ...]
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from plab.align import smith_waterman  # noqa: E402
from plab.catalog import CANDIDATES_BY_KEY  # noqa: E402
from plab.epitope import kmer_identity  # noqa: E402
from plab.pdbio import parse_structure  # noqa: E402

STRUCT_DIR = ROOT / "data" / "structures"
SEQ_DIR = ROOT / "data" / "sequences"

FAILED = ["hiv_gp120", "flu_ha", "il6r_ra", "rankl_osteoporosis", "il5_asthma"]


def read_fasta(acc: str) -> str | None:
    path = SEQ_DIR / f"{acc}.fasta"
    if not path.exists():
        return None
    return "".join(l.strip() for l in path.read_text(encoding="utf-8").splitlines()
                   if l.strip() and not l.startswith(">"))


def report(key: str) -> None:
    entry = CANDIDATES_BY_KEY.get(key)
    if entry is None:
        print(f"{key}: not a candidate")
        return
    target = read_fasta(entry.uniprot)
    path = STRUCT_DIR / f"{entry.pdb_complex}.pdb"
    print("=" * 78)
    print(f"{key}  gene {entry.gene}  acc {entry.uniprot} "
          f"({len(target) if target else '?'} aa)  pdb {entry.pdb_complex}")
    if target is None or not path.exists():
        print("  missing sequence or structure on disk; run build_catalog.py first")
        return
    struct = parse_structure(path)
    print(f"  {struct.title[:70]}")
    print(f"  {'chain':<6}{'len':>5}{'8mer':>7}{'aligned':>9}{'cover':>7}  sequence head")
    for cid in struct.chain_ids:
        seq = struct.sequence(cid)
        if len(seq) < 6:
            continue
        km = kmer_identity(seq, target)
        aln = smith_waterman(seq, target)
        cover = sum(1 for a, b in zip(aln.query_aligned, aln.target_aligned)
                    if a != "-" and b != "-") / len(seq)
        print(f"  {cid:<6}{len(seq):>5}{km:>6.0%}{aln.identity:>9.0%}"
              f"{cover:>7.0%}  {seq[:34]}")


def main() -> int:
    keys = sys.argv[1:] or FAILED
    for k in keys:
        report(k)
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
