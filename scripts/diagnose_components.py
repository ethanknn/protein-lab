"""Health check on every scoring component at once.

Three bugs in this project have been the same bug: a component that could not
vary. Checking them one at a time only found each after it had already
distorted a run, so this checks all of them together and reports the two
failure modes that matter.

SATURATION -- a component stuck at one value supplies no gradient. It spends
weight in the composite and contributes nothing to the search, which is worse
than omitting it.

REDUNDANCY -- two components measuring the same underlying quantity double
count it. The composite then silently weights that quantity twice, and a
candidate is penalised or rewarded twice for one property. Pairwise
correlation across a realistic sample exposes it.

Run: python scripts/diagnose_components.py [--plm] [--n 400]
"""

from __future__ import annotations

import math
import random
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from plab.aminoacids import AA20  # noqa: E402
from plab.design.scaffolds import get_scaffold  # noqa: E402
from plab.scoring.composite import GROUPS, Evaluator  # noqa: E402
from plab.targets import load_target  # noqa: E402

REDUNDANCY_THRESHOLD = 0.80
SATURATION_THRESHOLD = 0.60   # fraction stuck at one end
LOW_VARIANCE_SD = 0.04        # in 0..1 units


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n < 3:
        return 0.0
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if dx == 0 or dy == 0:
        return 0.0
    return num / (dx * dy)


def group_of(name: str) -> str:
    for g, names in GROUPS.items():
        if name in names:
            return g
    return "other"


def build_sample(scaffold, n: int, seed: int = 31) -> list[str]:
    """A realistic spread: random library members plus mutants of a good
    candidate, so the sample covers both the search space and the region the
    optimiser actually ends up in."""
    rng = random.Random(seed)
    ss = scaffold.secondary_structure
    pool = {p: [a for a in AA20 if a != "C"
                and not (ss[p] == "H" and a in "PG")]
            for p in scaffold.variable_positions}

    seqs = []
    for _ in range(n // 2):
        assign = {p: rng.choice(pool[p]) for p in scaffold.variable_positions}
        seqs.append(scaffold.apply(scaffold.framework, assign))

    # Mutants around a good candidate from the reference run.
    good = "VDNKFNKEYYEAYKEINELPNLNQDQKDAFIWSLNDDPSQSANLLAEAKKLNDAQAPK"
    for _ in range(n - len(seqs)):
        chars = list(good)
        for p in scaffold.variable_positions:
            if rng.random() < 0.25:
                chars[p] = rng.choice(pool[p])
        seqs.append("".join(chars))
    return seqs


def main() -> int:
    use_plm = "--plm" in sys.argv
    n = 400
    if "--n" in sys.argv:
        n = int(sys.argv[sys.argv.index("--n") + 1])

    target = load_target()
    ev = Evaluator(target, use_plm=use_plm)
    sc = get_scaffold("affibody_Z")

    print(f"Sampling {n} sequences on {sc.name}"
          f"{' (with ESM-2)' if use_plm else ' (ESM-2 off; pass --plm to include it)'}\n")

    seqs = build_sample(sc, n)
    series: dict[str, list[float]] = {}
    weights: dict[str, float] = {}
    composites: list[float] = []

    for i, s in enumerate(seqs):
        card = ev.evaluate(s, sc)
        composites.append(card.score)
        for c in card.components:
            if c.weight <= 0:
                continue
            series.setdefault(c.name, []).append(c.normalised)
            weights[c.name] = c.weight
        if use_plm and (i + 1) % 50 == 0:
            print(f"  scored {i + 1}/{len(seqs)}")

    names = sorted(series, key=lambda k: (group_of(k), k))

    # --- saturation ------------------------------------------------------
    print(f"\n{'=' * 78}")
    print("SATURATION: can each component vary?")
    print("=" * 78)
    print(f"  {'component':<26} {'group':<14} {'w':>4} {'mean':>6} {'sd':>6} "
          f"{'@0':>5} {'@1':>5}  verdict")
    flagged_sat = []
    for nm in names:
        v = series[nm]
        sd = statistics.stdev(v) if len(v) > 1 else 0.0
        at0 = sum(1 for x in v if x <= 0.001) / len(v)
        at1 = sum(1 for x in v if x >= 0.999) / len(v)
        verdict = "ok"
        if sd < LOW_VARIANCE_SD:
            verdict = "FLAT -- no gradient"
        elif at0 > SATURATION_THRESHOLD:
            verdict = "pinned at floor"
        elif at1 > SATURATION_THRESHOLD:
            verdict = "pinned at ceiling"
        if verdict != "ok":
            flagged_sat.append((nm, verdict))
        print(f"  {nm:<26} {group_of(nm):<14} {weights[nm]:>4.1f} "
              f"{statistics.fmean(v):>6.3f} {sd:>6.3f} "
              f"{at0:>5.0%} {at1:>5.0%}  {verdict}")

    # --- redundancy ------------------------------------------------------
    print(f"\n{'=' * 78}")
    print("REDUNDANCY: are two components measuring the same thing?")
    print("=" * 78)
    pairs = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            r = pearson(series[a], series[b])
            pairs.append((abs(r), r, a, b))
    pairs.sort(reverse=True)

    print("  strongest correlations:")
    for absr, r, a, b in pairs[:10]:
        mark = "  <-- DOUBLE COUNTING" if absr >= REDUNDANCY_THRESHOLD else ""
        print(f"    {r:+.3f}  {a} <-> {b}{mark}")

    flagged_red = [(a, b, r) for absr, r, a, b in pairs
                   if absr >= REDUNDANCY_THRESHOLD]

    # --- influence -------------------------------------------------------
    print(f"\n{'=' * 78}")
    print("INFLUENCE: how much does each component actually steer the score?")
    print("=" * 78)
    print("  Correlation with the composite, and the spread it contributes.")
    print("  A high weight with near-zero influence means the weight is wasted.\n")
    infl = []
    for nm in names:
        r = pearson(series[nm], composites)
        sd = statistics.stdev(series[nm]) if len(series[nm]) > 1 else 0.0
        infl.append((weights[nm] * sd, r, nm))
    infl.sort(reverse=True)
    print(f"  {'component':<26} {'w*sd':>7} {'r(composite)':>13}")
    for contrib, r, nm in infl:
        print(f"  {nm:<26} {contrib:>7.4f} {r:>13.3f}")

    # --- verdict ---------------------------------------------------------
    print(f"\n{'=' * 78}")
    print("VERDICT")
    print("=" * 78)
    if not flagged_sat and not flagged_red:
        print("  No component is saturated and none is redundant.")
    for nm, verdict in flagged_sat:
        print(f"  SATURATED  {nm}: {verdict}")
    for a, b, r in flagged_red:
        print(f"  REDUNDANT  {a} and {b} correlate at r={r:+.3f}.")
        print(f"             The composite weights that quantity "
              f"{weights[a] + weights[b]:.1f} in total, across two components "
              f"that move together.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
