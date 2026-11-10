"""Fetch and verify every catalogue candidate; write the ones that hold up.

For each candidate this does four things and stops at the first failure:

  1. fetch the UniProt sequence and check the gene symbol matches
  2. fetch the PDB structure
  3. check a chain in that structure really is this protein
  4. derive a binding site from a partner bound to it

A candidate that fails any step is recorded in a "rejected" list with the
reason and is never offered to a user. Shipping a broken entry would be worse
than shipping fewer entries: the pipeline would run, produce numbers, and
those numbers would be about nothing.

Run: python scripts/build_catalog.py [--only KEY] [--refresh]
"""

from __future__ import annotations

import json
import re
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from plab.catalog import (  # noqa: E402
    CANDIDATES, CATALOG_DIR, VERIFIED_PATH, CatalogEntry, CANDIDATES_BY_KEY,
)
from plab.epitope import EpitopeError, derive_epitope, describe_surface  # noqa: E402

STRUCT_DIR = ROOT / "data" / "structures"
SEQ_DIR = ROOT / "data" / "sequences"
EPI_DIR = CATALOG_DIR / "epitopes"
USER_AGENT = "protein-lab/0.2 (educational in-silico design sandbox)"


def _get(url: str, timeout: int = 60) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8")


def fetch_sequence(acc: str, refresh: bool = False) -> tuple[str, str] | None:
    """Return (header, sequence) or None."""
    dest = SEQ_DIR / f"{acc}.fasta"
    if dest.exists() and dest.stat().st_size > 0 and not refresh:
        text = dest.read_text(encoding="utf-8")
    else:
        try:
            text = _get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta")
        except Exception as exc:
            print(f"      fetch failed: {exc}")
            return None
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(text, encoding="utf-8")
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    if not lines or not lines[0].startswith(">"):
        return None
    header = lines[0]
    seq = "".join(l for l in lines if not l.startswith(">"))
    return (header, seq) if seq else None


def fetch_structure(pdb_id: str, refresh: bool = False) -> Path | None:
    """Fetch a structure, trying the legacy PDB format then mmCIF.

    The fallback is not a nicety. Entries past 62 chains or 99,999 atoms
    cannot be written in the legacy format at all, and the PDB no longer
    generates .pdb files for a growing share of ordinary entries either. Two
    candidates in the second research pass -- 9EZA and 9XN4 -- were rejected
    as "could not fetch PDB", which reads like a missing structure when the
    structure is there and only the file format is absent.
    """
    for ext in ("pdb", "cif"):
        dest = STRUCT_DIR / f"{pdb_id}.{ext}"
        if dest.exists() and dest.stat().st_size > 0 and not refresh:
            return dest
        try:
            text = _get(f"https://files.rcsb.org/download/{pdb_id}.{ext}")
        except urllib.error.HTTPError as exc:
            if exc.code != 404:
                print(f"      PDB {pdb_id}.{ext}: HTTP {exc.code}")
            continue
        except Exception as exc:
            print(f"      PDB {pdb_id}.{ext}: {exc}")
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(text, encoding="utf-8")
        return dest
    print(f"      PDB {pdb_id}: neither .pdb nor .cif could be fetched")
    return None


def verify(entry: CatalogEntry, refresh: bool = False) -> tuple[dict | None, str]:
    """Returns (verified_record, reason_if_rejected)."""
    got = fetch_sequence(entry.uniprot, refresh)
    if not got:
        return None, f"could not fetch UniProt {entry.uniprot}"
    header, sequence = got

    # The gene symbol must match, or the accession points at a different
    # protein and every number downstream would be about the wrong thing.
    #
    # Matched on whole tokens. An earlier version also accepted the symbol as
    # a bare substring anywhere in the header, which is far too loose: "HA"
    # or "C5" occurs inside plenty of unrelated protein names, so a wrong
    # accession could have passed this gate on a coincidence.
    gene = entry.gene
    gene_tokens = set(re.findall(r"[A-Za-z0-9]+", header))
    gene_ok = (f"GN={gene}" in header
               or gene in gene_tokens
               or gene.upper() in {t.upper() for t in gene_tokens})
    if not gene_ok:
        return None, (f"gene mismatch: expected {gene}, header reads "
                      f"{header[:90]}")

    path = fetch_structure(entry.pdb_complex, refresh)
    if path is None:
        return None, f"could not fetch PDB {entry.pdb_complex}"

    try:
        result = derive_epitope(
            path, sequence, prefer_binder=entry.binder_preference)
    except EpitopeError as exc:
        return None, str(exc)
    except Exception as exc:
        return None, f"{type(exc).__name__}: {exc}"

    if not result.ok:
        return None, (f"interface too small: only {len(result.epitope)} "
                      f"residues contacted")

    if entry.pdb_apo:
        fetch_structure(entry.pdb_apo, refresh)

    # Off-target sequences, so the selectivity panel is ready at run time.
    missing_off = []
    for acc in entry.off_targets:
        if fetch_sequence(acc, refresh) is None:
            missing_off.append(acc)
    if missing_off:
        print(f"      note: off-target sequences unavailable: {missing_off}")

    EPI_DIR.mkdir(parents=True, exist_ok=True)
    (EPI_DIR / f"{entry.key}.json").write_text(
        json.dumps(result.to_dict(), indent=2), encoding="utf-8")

    surface = describe_surface(result.epitope_sequence)
    record = entry.to_dict()
    record.update({
        "sequence_length": len(sequence),
        "uniprot_header": header,
        "resolution": result.resolution,
        "structure_title": result.title,
        "target_chain": result.target_chain,
        "binder_chains": result.binder_chains,
        "binder_kind": result.binder_kind,
        "epitope_sequence": result.epitope_sequence,
        "paratope_sequence": result.paratope_sequence,
        "n_epitope_residues": len(result.epitope),
        "descriptors": result.descriptors,
        "surface": surface,
        "target_identity": result.target_identity,
        "target_coverage": result.target_coverage,
        "identity_method": result.identity_method,
        "warnings": result.warnings,
        "epitope_file": f"epitopes/{entry.key}.json",
    })
    return record, ""


def main() -> int:
    refresh = "--refresh" in sys.argv
    only = None
    if "--only" in sys.argv:
        only = sys.argv[sys.argv.index("--only") + 1]

    candidates = [c for c in CANDIDATES if only is None or c.key == only]
    CATALOG_DIR.mkdir(parents=True, exist_ok=True)

    verified: dict[str, dict] = {}
    rejected: list[dict] = []

    print(f"Verifying {len(candidates)} candidate targets.")
    print("Each must resolve to a real protein with something bound to it.\n")

    for i, entry in enumerate(candidates, 1):
        print(f"[{i:>2}/{len(candidates)}] {entry.key:<22} "
              f"{entry.gene:<8} {entry.pdb_complex}  {entry.disease[:34]}")
        t0 = time.time()
        record, reason = verify(entry, refresh)
        dt = time.time() - t0
        if record is None:
            print(f"      REJECTED: {reason}\n")
            rejected.append({"key": entry.key, "gene": entry.gene,
                             "disease": entry.disease,
                             "pdb": entry.pdb_complex, "reason": reason})
            continue
        s = record["surface"]
        print(f"      ok ({dt:.1f}s): {record['n_epitope_residues']} residues "
              f"on chain {record['target_chain']}, bound by "
              f"{record['binder_kind']} {','.join(record['binder_chains'])}")
        print(f"      surface {record['epitope_sequence']} -> wants a "
              f"{s['required_character']} binder")
        for w in record["warnings"]:
            print(f"      note: {w}")
        print()
        verified[entry.key] = record

    payload = {
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_candidates": len(candidates),
        "n_verified": len(verified),
        "targets": verified,
        "rejected": rejected,
    }
    # Preserve previously verified entries when running with --only.
    if only and VERIFIED_PATH.exists():
        old = json.loads(VERIFIED_PATH.read_text(encoding="utf-8"))
        merged = old.get("targets", {})
        merged.update(verified)
        payload["targets"] = merged
        payload["n_verified"] = len(merged)
        # Carry forward the other rejections, but drop any whose candidate
        # no longer exists. Renaming an entry -- itgb6_ipf became
        # itgav_fibrosis when it turned out the antibody binds the alpha
        # subunit -- otherwise leaves its rejection on the page forever,
        # describing a candidate nobody can find in the source.
        live = set(CANDIDATES_BY_KEY)
        payload["rejected"] = [r for r in old.get("rejected", [])
                               if r["key"] != only and r["key"] in live
                               and r["key"] not in merged] + rejected

    VERIFIED_PATH.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print("=" * 72)
    print(f"VERIFIED {payload['n_verified']} / {len(candidates)} candidates")
    print("=" * 72)
    by_area: dict[str, list[str]] = {}
    for rec in payload["targets"].values():
        by_area.setdefault(rec["area"], []).append(
            f"{rec['gene']} ({rec['disease']})")
    for area, items in sorted(by_area.items()):
        print(f"  {area}:")
        for it in sorted(items):
            print(f"    - {it}")

    if payload["rejected"]:
        print(f"\nREJECTED {len(payload['rejected'])}:")
        for r in payload["rejected"]:
            print(f"  {r['key']} ({r['gene']}, {r['pdb']}): {r['reason'][:100]}")
        print("\nRejections are kept in the catalogue file and shown in the UI.")
        print("A candidate that cannot be verified is never offered as a target.")

    print(f"\nWrote {VERIFIED_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
