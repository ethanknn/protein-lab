"""Build a cross-reactivity panel for every catalogue target, from data.

The problem this solves
-----------------------
Selectivity -- would this binder also stick to something we need -- is the
axis that has actually killed therapeutic candidates. Until now only CLEC12A
had a panel of related proteins to score against, because I wrote that one by
hand. Forty-one targets reported "not assessed".

Writing forty-one more by hand is exactly the move that put a bacterial
encapsulin in this catalogue under the name of interleukin-5. So these are
derived instead:

  1. ask UniProt which domain families the target belongs to
  2. pull every reviewed HUMAN protein sharing one of those families
  3. align each one against the target's epitope WINDOW -- the binding site
     plus flanking context, not the whole protein
  4. keep the closest few that can actually have the epitope mapped onto them

Step 3 is the one that matters. Membership of a domain family says two
proteins are relatives; it does not say they resemble each other at the
particular patch of surface a binder would touch. EGFR and HER2 share a
family and most of a fold, but a binder against one face need not touch the
equivalent face of the other. Ranking on the epitope window asks the question
the score will actually be about.

What a panel does NOT mean
--------------------------
A panel is a list of proteins whose surfaces can be compared on composition.
It is not a cross-reactivity assay, and it says nothing about which cells
express these proteins. Those are the facts that decide whether an off-target
matters clinically, and they are measured on primary human samples.

A target with no human relatives -- a viral spike, a bacterial toxin, a
parasite protein -- gets an empty panel, and that is recorded as a distinct
state. "This protein has no human homologs" and "nobody looked" are different
claims and must not render the same way.

Run: python scripts/build_panels.py [--only KEY] [--max N]
"""

from __future__ import annotations

import json
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from plab.panels import (  # noqa: E402
    MAX_PANEL, PANELS_PATH, build_panel_for, save_sequences,
)
from plab.targets import available_targets, load_target  # noqa: E402


def build_panel(key: str, max_panel: int = MAX_PANEL) -> dict:
    """The catalogue-key form the CLI uses."""
    return build_panel_for(load_target(key), max_panel)


def main() -> int:
    only = None
    if "--only" in sys.argv:
        only = sys.argv[sys.argv.index("--only") + 1]
    max_panel = MAX_PANEL
    if "--max" in sys.argv:
        max_panel = int(sys.argv[sys.argv.index("--max") + 1])

    keys = sorted(available_targets())
    if only:
        keys = [k for k in keys if k == only]

    existing = {}
    if PANELS_PATH.exists():
        existing = json.loads(PANELS_PATH.read_text(encoding="utf-8")).get("panels", {})

    # --new: only targets that have no panel yet, for when the catalogue has
    # grown and rebuilding the rest would cost an hour for no change.
    if "--new" in sys.argv:
        keys = [k for k in keys if k not in existing]

    panels = dict(existing)
    print(f"Building cross-reactivity panels for {len(keys)} targets.\n")
    for i, key in enumerate(keys, 1):
        print(f"[{i:>2}/{len(keys)}] {key}")
        rec = build_panel(key, max_panel)
        save_sequences(rec)
        # The measured ceiling and its error bar survive a panel rebuild as
        # long as the members have not changed. They are what decides whether
        # a selectivity score is reported at all, and silently dropping them
        # falls back to the absolute cutoff -- the one that let VEGF-A report
        # 100% selective off a reachable margin of 0.029.
        prior = existing.get(key, {})
        same_members = ([m["uniprot"] for m in prior.get("members", [])]
                        == [m["uniprot"] for m in rec.get("members", [])])
        if same_members and prior.get("ceiling"):
            rec["ceiling"] = prior["ceiling"]
            rec["ceiling_noise"] = prior.get("ceiling_noise", {})
        elif prior.get("ceiling"):
            print("      panel membership changed; ceilings dropped, "
                  "re-run scripts/measure_ceilings.py")
        panels[key] = rec
        if rec["members"]:
            print(f"      {rec['state']}: "
                  + ", ".join(f"{m['gene']}({m['window_align_score']:.0f})"
                              for m in rec["members"]))
        else:
            print(f"      {rec['state']}: {rec['note'][:88]}")

    PANELS_PATH.parent.mkdir(parents=True, exist_ok=True)
    PANELS_PATH.write_text(json.dumps({
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "max_panel": max_panel,
        "panels": panels,
    }, indent=2), encoding="utf-8")

    by_state: dict[str, int] = {}
    for r in panels.values():
        by_state[r["state"]] = by_state.get(r["state"], 0) + 1
    print("\n" + "=" * 70)
    for state, n in sorted(by_state.items(), key=lambda kv: -kv[1]):
        print(f"  {state:<24}{n}")
    print(f"\nWrote {PANELS_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
