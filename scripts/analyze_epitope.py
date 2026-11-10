"""Derive the real CLEC12A epitope from the 8W9J crystal structure.

8W9J is human CLEC12A bound to 50C1, an inhibitory antibody that blocks
ligand binding. Rather than guessing where a designed binder should aim, we
read the actual contacts the antibody makes and use those residues as the
target surface.

This also extracts the antibody's own contact residues (its paratope), which
become the benchmark: a designed candidate is measured against what a
molecule that demonstrably engages this target actually does.

Output: data/targets/epitope_8W9J.json
Run:    python scripts/analyze_epitope.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from plab.pdbio import (  # noqa: E402
    find_contacts, interface_residues, parse_structure, radius_of_gyration,
)

STRUCT_DIR = ROOT / "data" / "structures"
TARGET_DIR = ROOT / "data" / "targets"

CONTACT_CUTOFF = 4.5

# Immunoglobulin J-segment signatures, used to tell heavy from light chains.
HEAVY_MOTIFS = ("WGQG", "WGRG", "WGAG", "WGKG", "WGPG", "WGMG")
LIGHT_MOTIFS = ("FGQG", "FGGG", "FGSG", "FGPG", "FGTG")


def load_fasta(path: Path) -> str:
    lines = [l.strip() for l in path.read_text(encoding="utf-8").splitlines()]
    return "".join(l for l in lines if l and not l.startswith(">"))


def kmer_identity(query: str, reference: str, k: int = 8) -> float:
    """Fraction of query k-mers found verbatim in reference.

    A crude but decisive test of 'is this chain the same protein'. Crystal
    constructs have tags and gaps, so exact alignment is overkill here.
    """
    if len(query) < k:
        return 0.0
    ref_kmers = {reference[i:i + k] for i in range(len(reference) - k + 1)}
    hits = sum(1 for i in range(len(query) - k + 1) if query[i:i + k] in ref_kmers)
    return hits / (len(query) - k + 1)


def classify_chain(seq: str, clec_seq: str) -> str:
    ident = kmer_identity(seq, clec_seq)
    if ident > 0.5:
        return "CLEC12A"
    if any(m in seq for m in HEAVY_MOTIFS):
        return "Fab_heavy"
    if any(m in seq for m in LIGHT_MOTIFS):
        return "Fab_light"
    return "unknown"


def main() -> int:
    pdb_path = STRUCT_DIR / "8W9J.pdb"
    if not pdb_path.exists():
        print(f"Missing {pdb_path}. Run scripts/fetch_data.py first.")
        return 1

    clec_seq = load_fasta(TARGET_DIR / "Q5QGZ9.fasta")
    print(f"CLEC12A reference (UniProt Q5QGZ9): {len(clec_seq)} aa\n")

    struct = parse_structure(pdb_path)
    print(struct.summary())
    print()

    # --- Classify every chain in the asymmetric unit ----------------------
    roles: dict[str, str] = {}
    print("Chain assignment:")
    for cid in struct.chain_ids:
        seq = struct.sequence(cid)
        if len(seq) < 20:
            roles[cid] = "too_short"
            print(f"  {cid}: {len(seq):>4} aa  -> skipped (too short)")
            continue
        role = classify_chain(seq, clec_seq)
        roles[cid] = role
        ident = kmer_identity(seq, clec_seq)
        print(f"  {cid}: {len(seq):>4} aa  -> {role:<10} (CLEC12A 8-mer identity {ident:.0%})")

    target_chains = [c for c, r in roles.items() if r == "CLEC12A"]
    binder_chains = [c for c, r in roles.items() if r in ("Fab_heavy", "Fab_light")]

    if not target_chains or not binder_chains:
        print("\nCould not identify both a CLEC12A chain and Fab chains.")
        print("Inspect the chain assignment above before continuing.")
        return 1

    print(f"\nTarget chains: {target_chains}")
    print(f"Binder chains: {binder_chains}")

    # --- Find the interface for each target copy, keep the largest ---------
    best = None
    for tc in target_chains:
        t_res = [r for r in struct.chain(tc) if r.is_standard]
        b_res = [r for r in struct.residues
                 if r.chain in binder_chains and r.is_standard]
        contacts = find_contacts(t_res, b_res, CONTACT_CUTOFF)
        if not contacts:
            continue
        epi, para = interface_residues(t_res, b_res, CONTACT_CUTOFF)
        print(f"\n  copy {tc}: {len(contacts)} residue contacts, "
              f"{len(epi)} epitope / {len(para)} paratope residues")
        if best is None or len(contacts) > len(best["contacts"]):
            best = {"target_chain": tc, "contacts": contacts,
                    "epitope": epi, "paratope": para, "target_residues": t_res}

    if best is None:
        print("\nNo inter-chain contacts found at the chosen cutoff.")
        return 1

    epi, para, contacts = best["epitope"], best["paratope"], best["contacts"]

    print(f"\n{'=' * 68}")
    print(f"EPITOPE on CLEC12A (chain {best['target_chain']}) "
          f"-- {len(epi)} residues contacted by 50C1")
    print("=" * 68)
    for r in epi:
        partners = sorted({c.res_b.label for c in contacts if c.res_a.key == r.key})
        closest = min(c.min_distance for c in contacts if c.res_a.key == r.key)
        print(f"  {r.label:<12} {r.resname}  closest {closest:.2f} A  "
              f"<- {', '.join(partners[:4])}")

    print(f"\n{'=' * 68}")
    print(f"PARATOPE on 50C1 Fab -- {len(para)} residues")
    print("=" * 68)
    by_chain: dict[str, list] = {}
    for r in para:
        by_chain.setdefault(r.chain, []).append(r)
    for cid, rs in sorted(by_chain.items()):
        role = roles.get(cid, "?")
        seq_str = "".join(r.one_letter for r in rs)
        nums = ", ".join(f"{r.one_letter}{r.resseq}" for r in rs)
        print(f"  chain {cid} ({role}): {seq_str}")
        print(f"    {nums}")

    # --- Interface descriptors we will reuse as design constraints ---------
    epi_seq = "".join(r.one_letter for r in epi)
    para_seq = "".join(r.one_letter for r in para)
    n_salt_bridge = sum(
        1 for c in contacts
        if {c.res_a.one_letter, c.res_b.one_letter} & set("DE")
        and {c.res_a.one_letter, c.res_b.one_letter} & set("KR")
        and c.min_distance <= 4.0
    )
    n_hydrophobic = sum(
        1 for c in contacts
        if c.res_a.one_letter in "AVILMFWCY" and c.res_b.one_letter in "AVILMFWCY"
    )

    print(f"\n{'=' * 68}")
    print("INTERFACE CHARACTER (these become the design objective)")
    print("=" * 68)
    print(f"  epitope composition      : {epi_seq}")
    print(f"  paratope composition     : {para_seq}")
    print(f"  residue-residue contacts : {len(contacts)}")
    print(f"  candidate salt bridges   : {n_salt_bridge}")
    print(f"  hydrophobic contacts     : {n_hydrophobic}")
    print(f"  epitope radius of gyration: {radius_of_gyration(epi):.2f} A")
    print(f"  tightest contact         : {contacts[0].describe()}")

    # --- Persist -----------------------------------------------------------
    payload = {
        "source_pdb": "8W9J",
        "source_title": struct.title,
        "resolution_angstrom": struct.resolution,
        "reference_binder": "50C1 Fab (inhibitory anti-CLEC12A antibody)",
        "contact_cutoff_angstrom": CONTACT_CUTOFF,
        "chain_roles": roles,
        "target_chain": best["target_chain"],
        "binder_chains": binder_chains,
        "epitope": [
            {
                "chain": r.chain,
                "resseq": r.resseq,
                "icode": r.icode.strip(),
                "resname": r.resname,
                "aa": r.one_letter,
                "label": r.label,
                "centroid": [round(float(v), 3) for v in r.centroid],
                "min_distance_to_binder": round(
                    min(c.min_distance for c in contacts if c.res_a.key == r.key), 3),
            }
            for r in epi
        ],
        "paratope": [
            {
                "chain": r.chain,
                "role": roles.get(r.chain, "?"),
                "resseq": r.resseq,
                "resname": r.resname,
                "aa": r.one_letter,
                "label": r.label,
            }
            for r in para
        ],
        "epitope_sequence": epi_seq,
        "paratope_sequence": para_seq,
        "descriptors": {
            "n_residue_contacts": len(contacts),
            "n_salt_bridges": n_salt_bridge,
            "n_hydrophobic_contacts": n_hydrophobic,
            "epitope_radius_of_gyration": round(radius_of_gyration(epi), 3),
        },
        "contacts": [
            {
                "target": c.res_a.label,
                "binder": c.res_b.label,
                "min_distance": round(c.min_distance, 3),
                "atom_pairs": c.n_atom_pairs,
            }
            for c in contacts
        ],
    }

    out = TARGET_DIR / "epitope_8W9J.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
