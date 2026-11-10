"""Gene symbol in, reviewed human UniProt accession out.

The companion to rcsb_lookup. That script exists because recalling PDB ids
from memory put a bacterial encapsulin in this catalogue under the name of
interleukin-5. Accessions are the same hazard one step earlier: P01375 and
P01374 are TNF and lymphotoxin-alpha, one digit apart and both plausible, and
nothing downstream would notice the swap -- the structure would verify, the
panel would build, and the target would simply be the wrong protein.

So the research loop is: symbol -> accession here, accession -> real PDB
entries in rcsb_lookup, and nothing is typed from memory at either step. It
earned this on its first run: the symbol SELP returns selenoprotein P before
it returns P-selectin, and MET returns a transcription modulator alongside
the kinase.

  python scripts/resolve_gene.py TNF IL6R CD38
  python scripts/resolve_gene.py --pick TNF      (also list usable structures)
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from plab.lookup import resolve_symbol  # noqa: E402


def main() -> int:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    also_pick = "--pick" in sys.argv
    if not args:
        print(__doc__)
        return 1
    found: list[str] = []
    for symbol in args:
        hits = resolve_symbol(symbol)
        if not hits:
            print(f"{symbol:<10} NO REVIEWED HUMAN ENTRY")
            continue
        for h in hits[:3]:
            flag = "  <-- ambiguous" if len(hits) > 1 else ""
            print(f"{symbol:<10} {h['accession']:<8} {str(h['length']):>5} aa  "
                  f"{h['name'][:52]}{flag}")
        found.append(hits[0]["accession"])
    if also_pick and found:
        print()
        from rcsb_lookup import print_pick  # noqa: E402
        print_pick(found, top=3)
    return 0


if __name__ == "__main__":
    sys.exit(main())
