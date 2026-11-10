"""How big is each panel's reachable margin, and how big is its error bar?

The selectivity score is a fraction of the best margin a target's epitope and
panel allow, so that ceiling decides what 100% means. For CLEC12A the ceiling
is 0.350. For VEGF-A it is 0.029 -- and VEGF-A's optimised design was reported
as 100% selective, which is true and useless in the same breath.

`INDISTINGUISHABLE_MARGIN`, an absolute 0.02, was meant to catch exactly that
and missed it by nine thousandths. The reason it was the wrong kind of
threshold is measurable: drop one residue from one mapped surrogate surface --
a single alignment error, which at 25-40% identity is a certainty rather than
a hypothetical -- and recompute. The ceiling moves by 0.003 for alpha-synuclein
and by 0.077 for TROP2, a factor of twenty-five. No single number can serve as
a cutoff across that range.

So the cutoff becomes relative: a ceiling must stand two standard deviations
clear of zero. This script measures both numbers and writes them into
panels.json, per face length, because the jackknife is a few hundred sweeps
per target and has no business running inside a design job.

Run: python scripts/measure_ceilings.py [KEY ...]
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from plab.design.scaffolds import REGISTRY, get_scaffold  # noqa: E402
from plab.scoring.selectivity import (  # noqa: E402
    CEILING_NOISE_FACTOR, achievable_margin_range, build_paralog_epitopes,
    ceiling_noise,
)
from plab.targets import PANELS_PATH, load_off_targets, load_target  # noqa: E402


def face_lengths() -> list[int]:
    """Every face size a design run can ask about, so nothing is computed
    for a scaffold and then read for another."""
    out = set()
    for name in REGISTRY:
        n = len(get_scaffold(name).variable_positions)
        if n:
            out.add(n)
    return sorted(out)


def main() -> int:
    doc = json.loads(PANELS_PATH.read_text(encoding="utf-8"))
    panels = doc.get("panels", {})
    keys = sys.argv[1:] or sorted(panels)
    lengths = face_lengths()
    print(f"face lengths in use: {lengths}")
    print(f"{'target':<24}{'face':>5}{'ceiling':>9}{'noise':>8}{'ratio':>7}"
          f"  reading")
    print("-" * 70)

    for key in keys:
        if key not in panels:
            print(f"{key:<24}  not in panels.json")
            continue
        target = load_target(key)
        surrogates = build_paralog_epitopes(target, load_off_targets(target))
        if not surrogates:
            continue
        epi = target.epitope_sequence
        ceilings: dict[str, float] = {}
        noises: dict[str, float] = {}
        for n in lengths:
            t0 = time.time()
            hi = achievable_margin_range(epi, surrogates, n)[1]
            noise = ceiling_noise(epi, surrogates, n)
            ceilings[str(n)] = round(hi, 5)
            noises[str(n)] = round(noise, 5)
            ratio = hi / noise if noise else float("inf")
            reading = ("withheld: ceiling is within its own error bar"
                       if ratio <= CEILING_NOISE_FACTOR else "scored")
            print(f"{key:<24}{n:>5}{hi:>9.3f}{noise:>8.3f}{ratio:>7.1f}"
                  f"  {reading}  ({time.time() - t0:.0f}s)")
        panels[key]["ceiling"] = ceilings
        panels[key]["ceiling_noise"] = noises

    doc["panels"] = panels
    doc["ceiling_noise_factor"] = CEILING_NOISE_FACTOR
    PANELS_PATH.write_text(json.dumps(doc, indent=1), encoding="utf-8")
    print(f"\nwrote ceilings for {len(keys)} target(s) to {PANELS_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
