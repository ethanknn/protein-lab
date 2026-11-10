"""Find out why the paralog-selectivity score never moves.

It has read 50% in every run, with a kernel margin of exactly 0.075 each
time. A score that does not vary is not measuring anything, so before
changing the formula, establish what it is actually capable of reporting.

The test: sweep many different binding faces across the whole composition
range and see how much the margin can possibly vary. If the spread is tiny,
the metric has a structural ceiling and the optimiser has been blind on this
axis the whole time -- the same failure mode as the saturated language-model
score, in a different place.

Run: python scripts/diagnose_selectivity.py
"""

from __future__ import annotations

import itertools
import random
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from plab.aminoacids import AA20  # noqa: E402
from plab.design.scaffolds import get_scaffold  # noqa: E402
from plab.scoring.composite import Evaluator  # noqa: E402
from plab.scoring.selectivity import complementarity_kernel  # noqa: E402
from plab.targets import load_target  # noqa: E402


def rule(t: str) -> None:
    print(f"\n{'=' * 72}\n{t}\n{'=' * 72}")


def surface_profile(name: str, seq: str) -> str:
    cat = sum(1 for a in seq if a in "RK")
    acid = sum(1 for a in seq if a in "DE")
    arom = sum(1 for a in seq if a in "FWY")
    hyd = sum(1 for a in seq if a in "AVILMFWCY")
    return (f"  {name:<12} {seq:<16} cations {cat}  acidic {acid}  "
            f"aromatic {arom}  hydrophobic {hyd}")


def main() -> int:
    target = load_target()
    ev = Evaluator(target, use_plm=False)
    sc = get_scaffold("affibody_Z")

    epi = target.epitope_sequence
    paralogs = ev.paralog_epitopes

    rule("THE SURFACES BEING COMPARED")
    print(surface_profile("CLEC12A", epi) + "   <- the target")
    for gene, seq in paralogs.items():
        print(surface_profile(gene, seq))

    rule("CAN THE MARGIN MOVE AT ALL?")
    print("Sweeping binding faces across the composition range and recording")
    print("the margin between CLEC12A and the best-scoring paralog.\n")

    rng = random.Random(11)
    margins: list[tuple[float, str, float, float]] = []

    # Systematic sweep: vary acidic / basic / aromatic counts on a 13-mer face.
    for n_acid in range(0, 9):
        for n_arom in range(0, 7):
            for n_base in range(0, 5):
                n_other = 13 - n_acid - n_arom - n_base
                if n_other < 0:
                    continue
                face = ("D" * n_acid + "Y" * n_arom + "K" * n_base
                        + "S" * n_other)
                on = complementarity_kernel(face, epi)
                off = max((complementarity_kernel(face, p)
                           for p in paralogs.values()), default=0.0)
                margins.append((on - off, face, on, off))

    # Random faces too, in case the systematic sweep misses something.
    for _ in range(400):
        face = "".join(rng.choice(AA20) for _ in range(13))
        on = complementarity_kernel(face, epi)
        off = max((complementarity_kernel(face, p)
                   for p in paralogs.values()), default=0.0)
        margins.append((on - off, face, on, off))

    vals = [m[0] for m in margins]
    print(f"  faces tested     : {len(margins)}")
    print(f"  margin range     : {min(vals):+.4f} to {max(vals):+.4f}")
    print(f"  margin mean      : {statistics.fmean(vals):+.4f} "
          f"+/- {statistics.stdev(vals):.4f}")
    print(f"  distinct values  : {len(set(round(v, 4) for v in vals))}")

    best = max(margins, key=lambda m: m[0])
    print(f"\n  best possible    : {best[0]:+.4f}  face {best[1]}")
    print(f"                     on-target {best[2]:.3f}, paralog {best[3]:.3f}")

    lo, hi = ev.selectivity_scale(13)
    print(f"\n  scoring scale    : {lo:+.4f} (10th pct of random faces) to "
          f"{hi:+.4f} (reachable optimum)")
    span = len({round(100 * max(0.0, min(1.0, (v - lo) / (hi - lo))))
                for v in vals})
    print(f"  distinct scores  : {span} of 101 possible")
    if span < 10:
        print("  -> TOO FEW. The component cannot supply a useful gradient.")
    else:
        print("  -> the component can discriminate across the range.")

    rule("HOW THE SURFACES PULL APART")
    print("A face scored against the target and every paralog. Wide spread")
    print("means the metric can tell them apart; values bunched together mean")
    print("it cannot.\n")

    demo = "DDDDYYYKKSSSS"
    print(f"  example face: {demo}")
    print(f"    vs CLEC12A ({epi}): {complementarity_kernel(demo, epi):.4f}")
    for gene, seq in paralogs.items():
        print(f"    vs {gene:<8} ({seq}): "
              f"{complementarity_kernel(demo, seq):.4f}")
    spread = [complementarity_kernel(demo, p) for p in paralogs.values()]
    on_demo = complementarity_kernel(demo, epi)
    print(f"\n  on-target {on_demo:.3f}, closest paralog {max(spread):.3f}, "
          f"spread across paralogs {max(spread) - min(spread):.3f}")

    rule("WHAT THE OPTIMISER SEES")
    print("Scoring real candidates from the runs:\n")
    seqs = {
        "ESM run best":
            "VDNKFNKEYYEAYDEIKKLPNLNQEQQKAFIDSLWDDPSQSANLLAEAKKLNDAQAPK",
        "fast run best":
            "VDNKFNKEFNKAYREITQLPNLNDDQYDAFIESLYDDPSQSANLLAEAKKLNDAQAPK",
        "wild-type Z":
            sc.framework,
    }
    for label, seq in seqs.items():
        card = ev.evaluate(seq, sc)
        comp = card.component("paralog_selectivity")
        print(f"  {label:<16} margin {comp.raw:+.4f} -> "
              f"{comp.normalised * 100:5.1f}%")

    print("\nIf these are all nearly identical, the component supplies no")
    print("gradient and the search cannot optimise selectivity even in")
    print("principle. That is the thing to fix.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
