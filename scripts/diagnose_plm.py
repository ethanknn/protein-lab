"""Measure what ESM-2 actually thinks of each class of sequence here.

The composite run scored a designed candidate at 0% on plausibility. Before
treating that as a finding, check whether the scale is fair: the reference is
a 245-residue antibody scFv, and ESM-2 has seen vast numbers of antibodies,
so its perplexity there will be unusually low. Anchoring a 58-residue designed
affibody to that number may be comparing across protein families rather than
measuring anything about the design.

This prints the landscape so the anchoring choice can be made on data.

Run: python scripts/diagnose_plm.py
"""

from __future__ import annotations

import random
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from plab.aminoacids import AA20  # noqa: E402
from plab.design.scaffolds import get_scaffold  # noqa: E402
from plab.scoring.composite import Evaluator  # noqa: E402
from plab.scoring.plm import get_runner  # noqa: E402
from plab.targets import load_target  # noqa: E402

N_DECOYS = 24


def main() -> int:
    runner = get_runner()
    if not runner.load():
        print(f"ESM-2 unavailable: {runner.status}")
        return 1
    print(f"Model: {runner.status}\n")

    target = load_target()
    sc = get_scaffold("affibody_Z")
    ev = Evaluator(target, use_plm=False)

    rows: list[tuple[str, str, float]] = []

    def measure(label: str, seq: str, group: str) -> float:
        ppl = runner.pseudo_perplexity(seq)
        rows.append((group, label, ppl))
        print(f"  {label:<44} {ppl:7.3f}   ({len(seq)} aa)")
        return ppl

    print("NATURAL PROTEINS (the scale's lower end)")
    wt = measure("wild-type Z domain (natural parent)", sc.framework, "natural")
    ref = ev.reference
    if ref:
        measure("50C1 scFv (antibody, 245 aa)", ref.sequence, "natural")
    measure("CLEC12A full ectodomain region",
            target.full_sequence[60:200], "natural")

    print("\nRANDOM LIBRARY MEMBERS (the scale's upper end)")
    rng = random.Random(42)
    ss = sc.secondary_structure
    decoy_ppls = []
    for i in range(N_DECOYS):
        assign = {}
        for p in sc.variable_positions:
            pool = [a for a in AA20 if a != "C"
                    and not (ss[p] == "H" and a in "PG")]
            assign[p] = rng.choice(pool)
        seq = sc.apply(sc.framework, assign)
        ppl = runner.pseudo_perplexity(seq)
        decoy_ppls.append(ppl)
        if i < 6:
            print(f"  decoy {i + 1:<38} {ppl:7.3f}")
    d_mean = statistics.fmean(decoy_ppls)
    d_sd = statistics.stdev(decoy_ppls)
    print(f"  ... {N_DECOYS} decoys: mean {d_mean:.3f} +/- {d_sd:.3f}, "
          f"range {min(decoy_ppls):.2f}-{max(decoy_ppls):.2f}")

    print("\nDESIGNED CANDIDATES (from the two runs)")
    designed = {
        "ESM run best  (WDETTYRQRYEEW)":
            "VDNKFNKEWDEATTEIYRLPNLNQRQYEAFIESLWDDPSQSANLLAEAKKLNDAQAPK",
        "fast run best (FNKYRTQDDYDEY)":
            "VDNKFNKEFNKAYREITQLPNLNDDQYDAFIESLYDDPSQSANLLAEAKKLNDAQAPK",
    }
    for label, seq in designed.items():
        measure(label, seq, "designed")

    print("\nNEGATIVE CONTROLS")
    measure("poly-Glu face", sc.apply(sc.framework,
            {p: "E" for p in sc.variable_positions}), "control")
    measure("poly-Arg face", sc.apply(sc.framework,
            {p: "R" for p in sc.variable_positions}), "control")
    shuffled = list(sc.framework)
    random.Random(7).shuffle(shuffled)
    measure("shuffled Z domain (same composition)",
            "".join(shuffled), "control")

    # --- the verdict ------------------------------------------------------
    print(f"\n{'=' * 70}")
    print("WHAT THIS MEANS FOR THE SCALE")
    print("=" * 70)
    best_designed = min(r[2] for r in rows if r[0] == "designed")
    print(f"  wild-type Z domain (natural)   : {wt:.2f}")
    print(f"  random library members         : {d_mean:.2f} +/- {d_sd:.2f}")
    print(f"  best designed candidate        : {best_designed:.2f}")
    if ref:
        ref_ppl = next(r[2] for r in rows if "scFv" in r[1])
        print(f"  50C1 antibody scFv             : {ref_ppl:.2f}")
        print()
        print(f"  The antibody is {wt / ref_ppl:.1f}x more 'natural-looking' to "
              f"ESM-2 than")
        print(f"  the wild-type Z domain is -- and the Z domain is a real, "
              f"natural protein.")
        print(f"  Anchoring the affibody scale to the antibody is therefore "
              f"comparing")
        print(f"  across protein families, not measuring design quality.")
    print()
    print(f"  Fair anchors for this scaffold: good = {wt:.2f} (the natural "
          f"parent),")
    print(f"  poor = {d_mean:.2f} (a random library member). The designed "
          f"candidate at")
    print(f"  {best_designed:.2f} then sits at "
          f"{100 * max(0.0, min(1.0, (d_mean - best_designed) / (d_mean - wt))):.0f}% "
          f"on that scale.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
