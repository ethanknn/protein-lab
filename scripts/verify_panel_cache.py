"""Prove the stored surrogate surfaces are what re-deriving them would give.

scripts/build_panels.py maps the epitope onto each panel member by local
alignment and writes the result into panels.json. Everything downstream then
re-derives exactly the same alignments, once per Evaluator, which is the
single most expensive thing this package does.

Using the stored surfaces instead is only safe if they still describe the
current binding site. Rebuilding the catalogue can move an epitope -- it has,
twice -- and a surrogate mapped onto the old one would be a confident wrong
answer rather than a visible failure.

So this does not simply stamp the panel file as current. It re-derives every
surface and compares, and records the epitope only for the panels that match.
A panel that disagrees is left without the marker, which makes the consumer
fall back to deriving the surfaces itself, and is reported here.

Run: python scripts/verify_panel_cache.py [KEY ...]
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from plab.panels import PANELS_PATH  # noqa: E402
from plab.scoring.selectivity import epitope_window, map_epitope_onto  # noqa: E402
from plab.targets import load_off_targets, load_target  # noqa: E402


def main() -> int:
    doc = json.loads(PANELS_PATH.read_text(encoding="utf-8"))
    panels = doc.get("panels", {})
    keys = sys.argv[1:] or sorted(panels)

    ok = stale = skipped = 0
    t0 = time.time()
    for key in keys:
        rec = panels.get(key)
        if not rec or not rec.get("members"):
            continue
        try:
            target = load_target(key)
        except Exception as exc:
            print(f"  {key}: cannot load ({type(exc).__name__})")
            skipped += 1
            continue

        got = epitope_window(target)
        if got is None:
            print(f"  {key}: epitope window unmappable; cache left unmarked")
            rec.pop("target_epitope", None)
            stale += 1
            continue
        window, positions = got

        by_acc = {o.uniprot: o for o in load_off_targets(target)}
        mismatches = []
        for m in rec["members"]:
            off = by_acc.get(m["uniprot"])
            stored = m.get("surrogate_surface")
            if off is None or not stored:
                mismatches.append(f"{m['gene']} missing")
                continue
            fresh = map_epitope_onto(window, positions, off.sequence)
            if fresh is None or fresh.surface != stored:
                mismatches.append(
                    f"{m['gene']} {stored} != "
                    f"{fresh.surface if fresh else 'unmappable'}")

        if mismatches:
            rec.pop("target_epitope", None)
            stale += 1
            print(f"  {key}: STALE -- {'; '.join(mismatches[:3])}")
        else:
            rec["target_epitope"] = target.epitope_sequence
            ok += 1

    doc["panels"] = panels
    PANELS_PATH.write_text(json.dumps(doc, indent=2), encoding="utf-8")
    print(f"\n{ok} panels verified and marked, {stale} stale, "
          f"{skipped} skipped, in {time.time() - t0:.0f}s")
    if stale:
        print("Stale panels fall back to deriving their surfaces at run time. "
              "Re-run scripts/build_panels.py --only <key> to fix one.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
