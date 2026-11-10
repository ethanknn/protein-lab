"""Ask RCSB what a structure actually contains, and which structures contain a protein.

Written after five catalogue candidates were rejected and the cause turned out
to be, in three cases, that the PDB id I had paired with an accession did not
hold that protein -- 4O9H is interleukin-6, not its receptor; 7KQ5 is a
bacterial encapsulin and has nothing to do with interleukin-5. Guessing PDB
ids from memory is how that happens. The fix is to stop guessing: RCSB will
say which UniProt accession every chain maps to, and which entries contain a
given accession alongside an antibody.

Three queries:

  entities PDBID...     what each chain in these entries is
  find ACCESSION        entries containing this protein, antibody ones first
  pick ACCESSION...     the entries most usable as a design target

The lookups themselves now live in src/plab/lookup.py, because the web
application's "design your own" tab runs the same two steps for whatever
protein someone types in, and a second copy of this logic would be a second
chance to get it wrong.

Run: python scripts/rcsb_lookup.py entities 4FQI 4GIQ
     python scripts/rcsb_lookup.py pick P05113
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from plab.lookup import (  # noqa: E402
    entity_map, looks_like_antibody, rank_structures, structures_for,
    uniprot_header,
)


def print_entities(pdb_ids: list[str]) -> None:
    for rec in entity_map(pdb_ids):
        res = rec["resolution"]
        print("=" * 78)
        print(f"{rec['pdb_id']}  {rec['method']}  "
              + (f"{res:.2f} A" if res else ""))
        print(f"  {rec['title'][:74]}")
        for e in rec["entities"]:
            acc = ",".join(a for a in e["accessions"] if a) or "-"
            print(f"  chains {','.join(e['chains']):<10} {str(e['length']):>5} aa  "
                  f"{acc:<12} {e['description'][:40]}")
            if e["organism"]:
                print(f"         {e['organism']}")


def print_find(accession: str, limit: int = 25) -> None:
    ids = structures_for(accession, limit)
    if not ids:
        print(f"{accession}: no structures found")
        return
    print(f"{accession}  {uniprot_header(accession)[:70]}")
    for rec in entity_map(ids):
        ab = any(looks_like_antibody(e["description"]) for e in rec["entities"])
        print(f"  {rec['pdb_id']}  {rec['resolution']} A  "
              f"{'antibody' if ab else 'partner '}  {rec['title'][:50]}")


def print_pick(accessions: list[str], limit: int = 60, top: int = 3) -> None:
    """For each accession, the structures most usable as a design target.

    This is the step that replaces recalling a PDB id.
    """
    for acc in accessions:
        print("=" * 78)
        print(f"{acc}  {uniprot_header(acc)[:72]}")
        ranked = rank_structures(acc, limit=limit, top=top)
        if not ranked:
            print("  no entries found with this protein and a bound partner")
            continue
        for r in ranked:
            kind = "antibody" if r["has_antibody"] else "partner"
            print(f"  {r['pdb_id']}  {r['resolution']} A  {kind}"
                  + ("  [large assembly]" if r["large_assembly"] else ""))
            print(f"    {r['title'][:70]}")
            print(f"    target  {','.join(r['target_chains'])}  "
                  f"{r['target_length']} aa  {r['target_description'][:38]}")
            for o in r["partners"][:3]:
                print(f"    bound   {','.join(o['chains']):<8} "
                      f"{str(o['length']):>4} aa  {o['description'][:38]}")
        print()


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    cmd = sys.argv[1]
    if cmd == "entities":
        print_entities(sys.argv[2:])
    elif cmd == "find":
        for acc in sys.argv[2:]:
            print_find(acc)
            print()
    elif cmd == "pick":
        print_pick(sys.argv[2:])
    else:
        print(__doc__)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
