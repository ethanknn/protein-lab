"""Download the public experimental data this project is built on.

Everything fetched here is open data from the RCSB Protein Data Bank and
UniProt. Nothing is generated or guessed -- if a file is missing, the
pipeline should fail loudly rather than fall back to invented values.

Run:  python scripts/fetch_data.py
"""

from __future__ import annotations

import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
STRUCT_DIR = ROOT / "data" / "structures"
TARGET_DIR = ROOT / "data" / "targets"

USER_AGENT = "protein-lab/0.1 (educational in-silico design sandbox)"

# PDB entries. Provenance for each one is recorded so the dashboard can cite it.
PDB_ENTRIES = {
    "8W8T": {
        "title": "Human CLEC12A carbohydrate-recognition domain (apo)",
        "role": "target_apo",
        "note": "Unliganded CLEC12A CRD. Used as the design target surface.",
    },
    "8W9J": {
        "title": "Human CLEC12A ectodomain complexed with 50C1 Fab",
        "role": "reference_complex",
        "note": (
            "CLEC12A bound to the validated inhibitory antibody 50C1 at 3.5 A. "
            "Used to derive the real epitope and as the benchmark binder."
        ),
    },
    "8W9H": {
        "title": "Anti-human CLEC12A antibody 50C1 Fab (unbound)",
        "role": "reference_binder_apo",
        "note": "50C1 Fab alone at 2.0 A. Higher-resolution source for CDR geometry.",
    },
    "1Q2N": {
        "title": "Z domain of staphylococcal protein A (NMR)",
        "role": "scaffold_template",
        "note": (
            "The three-helix bundle used as the binder scaffold. Real "
            "coordinates, so designed candidates are displayed on a genuine "
            "fold rather than an idealised cartoon."
        ),
    },
}

# UniProt accessions. CLEC12A plus the off-target controls we score against.
UNIPROT_ENTRIES = {
    "Q5QGZ9": {
        "name": "CLEC12A_HUMAN",
        "expect_gene": "CLEC12A",
        "label": "CLEC12A / CLL-1 (human) - the AML target",
        "role": "target",
    },
    # The off-target panel. A binder that engages any of these is a problem,
    # and the closer the relative the more likely the cross-reaction.
    "Q2HXU8": {
        "name": "CLC12B_HUMAN",
        "expect_gene": "CLEC12B",
        "label": "CLEC12B - immediate neighbour paralog, same gene cluster",
        "role": "off_target_paralog",
    },
    "Q6UXN8": {
        "name": "CLEC9A_HUMAN",
        "expect_gene": "CLEC9A",
        "label": "CLEC9A / DNGR-1 - C-type lectin on dendritic cells",
        "role": "off_target_paralog",
    },
    "Q9BXN2": {
        "name": "CLC7A_HUMAN",
        "expect_gene": "CLEC7A",
        "label": "CLEC7A / Dectin-1 - myeloid C-type lectin, same family",
        "role": "off_target_paralog",
    },
    "Q8NC01": {
        "name": "CLC1A_HUMAN",
        "expect_gene": "CLEC1A",
        "label": "CLEC1A - C-type lectin, NK gene complex",
        "role": "off_target_paralog",
    },
    "Q12918": {
        "name": "KLRB1_HUMAN",
        "expect_gene": "KLRB1",
        "label": "KLRB1 / CD161 - C-type lectin on NK and T cells",
        "role": "off_target_paralog",
    },
    "P20138": {
        "name": "CD33_HUMAN",
        "expect_gene": "CD33",
        "label": "CD33 - myeloid antigen with known HSC toxicity",
        "role": "off_target_reference",
    },
}


def _get(url: str, binary: bool = False):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=60) as resp:
        raw = resp.read()
    return raw if binary else raw.decode("utf-8")


def fetch_pdb(pdb_id: str) -> dict:
    """Fetch a structure, preferring .pdb and falling back to .cif."""
    result = {"id": pdb_id, "files": {}}
    for ext in ("pdb", "cif"):
        dest = STRUCT_DIR / f"{pdb_id}.{ext}"
        if dest.exists() and dest.stat().st_size > 0:
            print(f"  [cached] {dest.name} ({dest.stat().st_size:,} bytes)")
            result["files"][ext] = dest.name
            continue
        url = f"https://files.rcsb.org/download/{pdb_id}.{ext}"
        try:
            text = _get(url)
        except urllib.error.HTTPError as exc:
            print(f"  [skip]   {pdb_id}.{ext} -> HTTP {exc.code}")
            continue
        except Exception as exc:  # network trouble: report, do not invent
            print(f"  [fail]   {pdb_id}.{ext} -> {exc}")
            continue
        dest.write_text(text, encoding="utf-8")
        print(f"  [ok]     {dest.name} ({len(text):,} bytes)")
        result["files"][ext] = dest.name
    return result


def fetch_uniprot(acc: str) -> dict:
    """Fetch a FASTA sequence from UniProt."""
    dest = STRUCT_DIR.parent / "targets" / f"{acc}.fasta"
    if dest.exists() and dest.stat().st_size > 0:
        text = dest.read_text(encoding="utf-8")
        print(f"  [cached] {dest.name}")
    else:
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        try:
            text = _get(url)
        except Exception as exc:
            print(f"  [fail]   {acc}.fasta -> {exc}")
            return {"accession": acc, "sequence": None}
        dest.write_text(text, encoding="utf-8")
        print(f"  [ok]     {dest.name}")
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    header = lines[0] if lines and lines[0].startswith(">") else ""
    seq = "".join(ln for ln in lines if not ln.startswith(">"))
    return {"accession": acc, "header": header, "sequence": seq, "length": len(seq)}


def verify_gene(entry: dict, expected: str | None) -> bool:
    """Confirm the downloaded sequence is the protein we meant to fetch.

    An accession typo would silently substitute a different protein and the
    selectivity comparison would then be against the wrong thing -- worse
    than not comparing at all, because the number would still look fine.
    """
    if not expected:
        return True
    header = entry.get("header", "")
    ok = f"GN={expected}" in header or f"{expected}_HUMAN" in header
    if not ok:
        print(f"    !! MISMATCH: expected gene {expected}, header says:")
        print(f"       {header[:110]}")
    return ok


def main() -> int:
    STRUCT_DIR.mkdir(parents=True, exist_ok=True)
    TARGET_DIR.mkdir(parents=True, exist_ok=True)

    manifest = {"pdb": {}, "uniprot": {}}

    print("Structures (RCSB PDB):")
    for pdb_id, meta in PDB_ENTRIES.items():
        print(f"  {pdb_id}: {meta['title']}")
        entry = fetch_pdb(pdb_id)
        entry.update(meta)
        manifest["pdb"][pdb_id] = entry

    print("\nSequences (UniProt):")
    mismatches = []
    for acc, meta in UNIPROT_ENTRIES.items():
        print(f"  {acc}: {meta['label']}")
        entry = fetch_uniprot(acc)
        entry.update(meta)
        if entry.get("sequence") and not verify_gene(entry, meta.get("expect_gene")):
            mismatches.append(acc)
            entry["verified"] = False
        else:
            entry["verified"] = True
            if entry.get("sequence"):
                print(f"           verified {meta.get('expect_gene')}, "
                      f"{entry['length']} aa")
        manifest["uniprot"][acc] = entry

    if mismatches:
        print(f"\n  WARNING: {mismatches} did not match the expected gene "
              f"symbol. Fix the accessions before trusting any selectivity "
              f"comparison that uses them.")

    out = TARGET_DIR / "manifest.json"
    out.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"\nManifest -> {out}")

    missing = [k for k, v in manifest["pdb"].items() if not v["files"]]
    missing += [k for k, v in manifest["uniprot"].items() if not v.get("sequence")]
    if missing:
        print(f"\nWARNING: no data for {missing}. Downstream steps need these.")
        return 1
    print("\nAll source data present.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
