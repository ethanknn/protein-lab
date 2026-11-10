"""Fetch the transcript screen behind src/plab/lineage.py.

Writes data/lineage/expression.json: for every catalogue gene, the Human
Protein Atlas nTPM in each blood cell type and screened tissue that
lineage.py asks about.

Two ways this API will quietly lie to you, both guarded here:

1. The search is fuzzy. Asking for CD33 returns fourteen genes -- SIGLEC6,
   SIGLEC5, ABCG2, FGFR1..4, JAG1, NCR1..3, SIGLEC15, SIGLEC8 -- in that
   order, with CD33 first only by luck. Taking the first row would store one
   protein's expression under another's name, which is the same class of
   mistake as the SELP symbol resolving to selenoprotein P instead of
   P-selectin. Every row is matched on the exact gene symbol.

2. An unrecognised column code is silently dropped. `t_RNA_skin` is not a
   column -- the real one is `t_RNA_skin_1` -- and asking for it returns a
   perfectly well-formed row with no skin in it and no error anywhere. A
   typo therefore reads downstream as "not expressed in skin". Every
   requested column is checked against what came back, and anything missing
   is recorded by name in the output rather than left to look like a zero.

Run: python scripts/fetch_lineage.py [GENE ...]
"""

from __future__ import annotations

import json
import sys
import time
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from plab.lineage import ALL_COLUMNS, EXPRESSION_PATH  # noqa: E402

API = "https://www.proteinatlas.org/api/search_download.php"
USER_AGENT = "protein-lab/0.4 (educational in-silico design sandbox)"
SOURCE = "Human Protein Atlas (proteinatlas.org), search_download API"
# Pinned: the column codes and cell-type vocabulary change between releases,
# and an unpinned fetch would silently mix two vocabularies in one file.
VERSION = "v24"


def fetch(gene: str, columns: tuple[str, ...], timeout: int = 60,
          ) -> tuple[dict[str, float], list[str], list[str]]:
    """nTPM per column for exactly this gene, the columns that vanished, and
    the gene's synonyms.

    The synonyms are not decoration. The antigen checklist a clinician fills
    in is read in CD numbers -- CD371, CD33, CD38 -- while this catalogue is
    keyed by HGNC symbol, and CLEC12A and CD371 are the same molecule. Asking
    somebody to recognise their own flow panel under the wrong vocabulary is
    a good way to have an antigen marked "not tested" when it was tested.
    """
    url = (f"{API}?search={urllib.parse.quote(gene)}&format=json"
           f"&columns=g,gs,{','.join(columns)}&compress=no")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        rows = json.loads(resp.read().decode("utf-8"))

    exact = [r for r in rows if (r.get("Gene") or "").upper() == gene.upper()]
    if not exact:
        return {}, list(columns), []
    row = exact[0]
    syn = row.get("Gene synonym") or []
    if isinstance(syn, str):
        syn = [s.strip() for s in syn.split(",") if s.strip()]

    # The API returns pretty labels ("Blood RNA - neutrophil [nTPM]"), not the
    # codes it was asked for, so match them back by the code's own tail.
    labels = {k: k for k in row if k != "Gene"}
    values: dict[str, float] = {}
    missing: list[str] = []
    for code in columns:
        tail = code.split("_RNA_", 1)[-1].replace("_", " ").lower()
        hit = next((k for k in labels
                    if k.lower().startswith(("blood rna - ", "tissue rna - "))
                    and k.lower().split(" - ", 1)[1]
                    .removesuffix(" [ntpm]").strip() == tail), None)
        if hit is None or row.get(hit) in (None, ""):
            missing.append(code)
            continue
        try:
            values[code] = float(row[hit])
        except (TypeError, ValueError):
            missing.append(code)
    return values, missing, [s for s in syn if s]


def catalogue_genes() -> list[str]:
    doc = json.loads((ROOT / "data" / "catalog" / "verified.json")
                     .read_text(encoding="utf-8"))
    genes = {t["gene"] for t in doc["targets"].values() if t.get("gene")}
    return sorted(genes)


def main() -> int:
    genes = sys.argv[1:] or catalogue_genes()
    print(f"{len(genes)} genes x {len(ALL_COLUMNS)} columns, {VERSION}\n")

    values: dict[str, dict[str, float]] = {}
    missing: dict[str, list[str]] = {}
    synonyms: dict[str, list[str]] = {}
    unresolved: list[str] = []
    t0 = time.time()

    for n, gene in enumerate(genes, 1):
        try:
            vals, gone, syn = fetch(gene, ALL_COLUMNS)
        except Exception as exc:
            print(f"  {gene}: fetch failed ({type(exc).__name__})")
            unresolved.append(gene)
            continue
        if not vals:
            print(f"  {gene}: no exact symbol match in the atlas")
            unresolved.append(gene)
            continue
        values[gene.upper()] = vals
        if gone:
            missing[gene.upper()] = gone
        if syn:
            synonyms[gene.upper()] = syn
        if n % 10 == 0 or n == len(genes):
            print(f"  {n}/{len(genes)} ({time.time() - t0:.0f}s)")

    # A column missing for EVERY gene is a bad column code, not biology. That
    # distinction is the whole reason `missing` is recorded, so say it here
    # rather than leaving it in the file for nobody to read.
    if values:
        always_gone = [c for c in ALL_COLUMNS
                       if all(c not in v for v in values.values())]
        if always_gone:
            print("\n!! these columns came back for NO gene, which means the "
                  "column code is wrong rather than the gene being silent:")
            for c in always_gone:
                print(f"     {c}")

    EXPRESSION_PATH.parent.mkdir(parents=True, exist_ok=True)
    EXPRESSION_PATH.write_text(json.dumps({
        "source": SOURCE,
        "version": VERSION,
        "fetched": date.today().isoformat(),
        "columns": list(ALL_COLUMNS),
        "unresolved": unresolved,
        "missing": missing,
        "synonyms": synonyms,
        "values": values,
    }, indent=2, sort_keys=True), encoding="utf-8")

    print(f"\n{len(values)} genes written, {len(unresolved)} unresolved, "
          f"in {time.time() - t0:.0f}s -> {EXPRESSION_PATH.name}")
    if unresolved:
        print("Unresolved: " + ", ".join(unresolved))
    return 0


if __name__ == "__main__":
    sys.exit(main())
