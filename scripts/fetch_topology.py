"""Fetch where each catalogue protein actually sits: surface, soluble, inside.

The lineage screen in src/plab/lineage.py asks one question of every target:
which healthy cells carry it. That question assumes the protein is displayed
on a cell surface where a binder could meet it. For most of this catalogue
it is not.

Asking it anyway produces confident nonsense. Before this script existed the
screen reported alpha-synuclein at 103.2 nTPM on dendritic cells and listed
them as a population at risk -- but SNCA is cytoplasmic, and no antibody
reaches the cytoplasm, so the number is a measurement of what those cells
MAKE rather than of anything a binder could bind. TNF read 73.4 nTPM on
monocytes, which are the cells that produce it.

So the classification has to come from somewhere, and the one place it must
not come from is recall. Three of this catalogue's early rejections were PDB
ids remembered as the wrong protein; a remembered subcellular location is the
same mistake somewhere nothing downstream would catch it. This fetches the
annotation from UniProt instead.

What is stored is EVIDENCE, not a verdict: the keywords, the subcellular
location strings, and the transmembrane feature count, per accession. The
rule that turns those into "surface", "shed", "secreted" or "intracellular"
lives in lineage.classify_topology, where it can be read, tested and argued
with. A verdict frozen into a data file is a verdict nobody can check.

Two ways this API will quietly mislead you, both guarded here:

1. UniProt serves secondary accessions. Ask for an accession that has been
   merged into another entry and you get a perfectly valid record for a
   DIFFERENT protein, with no error. Every response is checked against the
   accession that was requested, and a mismatch is recorded by name rather
   than stored under the name we asked for.

2. An annotation that is absent for every protein means the field name is
   wrong, not that ninety-five proteins all lack it. The same trap as
   `t_RNA_skin` in fetch_lineage.py, which cost a wrong paragraph in a
   docstring before it was caught.

Run: python scripts/fetch_topology.py [ACCESSION ...]
"""

from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.request
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

TOPOLOGY_PATH = ROOT / "data" / "lineage" / "topology.json"

API = "https://rest.uniprot.org/uniprotkb"
USER_AGENT = "protein-lab/0.4 (educational in-silico design sandbox)"
SOURCE = "UniProtKB (uniprot.org), REST entry API"

#: Requested explicitly so a renamed field shows up as missing for every
#: entry rather than as an absent annotation.
FIELDS = ("accession,id,protein_name,cc_subcellular_location,"
          "ft_signal,ft_transmem,ft_topo_dom,ft_lipid,keyword")

#: Feature types that put part of the chain in the membrane.
MEMBRANE_FEATURES = ("Transmembrane", "Intramembrane")


def fetch(accession: str, timeout: int = 60) -> dict:
    """Topology evidence for exactly this accession.

    Raises on a mismatch rather than returning the wrong protein's
    annotation, because a silently substituted entry is worse than no entry:
    one is a gap and the other is a confident wrong answer.
    """
    url = f"{API}/{accession}.json?fields={FIELDS}"
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        doc = json.loads(resp.read().decode("utf-8"))

    got = doc.get("primaryAccession", "")
    if got.upper() != accession.upper():
        raise ValueError(
            f"asked for {accession}, served {got} -- a secondary accession "
            f"that now redirects to another entry")

    keywords = sorted({k.get("name", "") for k in doc.get("keywords", [])
                       if k.get("name")})

    locations: list[str] = []
    for c in doc.get("comments", []):
        if c.get("commentType") != "SUBCELLULAR LOCATION":
            continue
        for s in c.get("subcellularLocations", []):
            v = (s.get("location") or {}).get("value")
            if v and v not in locations:
                locations.append(v)

    features: dict[str, int] = {}
    for f in doc.get("features", []):
        t = f.get("type")
        if t:
            features[t] = features.get(t, 0) + 1

    return {
        "accession": got,
        "entry": doc.get("uniProtkbId", ""),
        "keywords": keywords,
        "locations": locations,
        "features": features,
        "membrane_spans": sum(features.get(t, 0) for t in MEMBRANE_FEATURES),
        "has_signal_peptide": features.get("Signal", 0) > 0,
    }


def catalogue() -> list[tuple[str, str]]:
    doc = json.loads((ROOT / "data" / "catalog" / "verified.json")
                     .read_text(encoding="utf-8"))
    pairs = {(t["gene"], t["uniprot"]) for t in doc["targets"].values()
             if t.get("gene") and t.get("uniprot")}
    return sorted(pairs)


def main() -> int:
    args = sys.argv[1:]
    pairs = ([(a, a) for a in args] if args else catalogue())
    print(f"{len(pairs)} accessions from UniProt\n")

    entries: dict[str, dict] = {}
    unresolved: dict[str, str] = {}
    t0 = time.time()

    for n, (gene, acc) in enumerate(pairs, 1):
        try:
            entries[gene.upper()] = fetch(acc)
        except (urllib.error.HTTPError, urllib.error.URLError,
                ValueError, json.JSONDecodeError) as exc:
            print(f"  {gene} ({acc}): {exc}")
            unresolved[gene.upper()] = f"{type(exc).__name__}: {exc}"
            continue
        if n % 20 == 0 or n == len(pairs):
            print(f"  {n}/{len(pairs)} ({time.time() - t0:.0f}s)")

    # An annotation absent for EVERY entry means the field name is wrong.
    # Stated here rather than left in the file for nobody to read.
    if entries:
        if not any(e["keywords"] for e in entries.values()):
            print("\n!! no entry returned any keyword, which means the "
                  "keyword field is not being requested correctly")
        if not any(e["locations"] for e in entries.values()):
            print("\n!! no entry returned a subcellular location, which "
                  "means that field is not being requested correctly")
        if not any(e["membrane_spans"] for e in entries.values()):
            print("\n!! no entry returned a transmembrane feature, which is "
                  "not credible for this catalogue -- check the field name")

    TOPOLOGY_PATH.parent.mkdir(parents=True, exist_ok=True)
    TOPOLOGY_PATH.write_text(json.dumps({
        "source": SOURCE,
        "fetched": date.today().isoformat(),
        "fields": FIELDS,
        "note": "Evidence only. The rule that turns this into a topology "
                "verdict lives in plab.lineage.classify_topology, so it can "
                "be read and tested rather than trusted.",
        "unresolved": unresolved,
        "entries": entries,
    }, indent=2, sort_keys=True), encoding="utf-8")

    print(f"\n{len(entries)} written, {len(unresolved)} unresolved, "
          f"in {time.time() - t0:.0f}s -> {TOPOLOGY_PATH.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
