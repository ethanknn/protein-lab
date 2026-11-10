"""Check every stage of the pipeline before trusting a full run.

Run: python scripts/smoke_test.py [--plm]
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from plab.design.scaffolds import get_scaffold, list_scaffolds  # noqa: E402
from plab.scoring.composite import Evaluator  # noqa: E402
from plab.targets import load_off_targets, load_target  # noqa: E402


def rule(title: str) -> None:
    print(f"\n{'=' * 72}\n{title}\n{'=' * 72}")


def main() -> int:
    use_plm = "--plm" in sys.argv

    rule("1. TARGET BRIEF")
    target = load_target()
    print(target.design_brief())
    print(f"\nWhy this target:\n  {target.rationale}")
    print(f"\nReality check:\n  {target.clinical_note}")

    rule("2. OFF-TARGET PANEL")
    for ot in load_off_targets():
        print(f"  {ot.gene:<8} ({ot.uniprot}, {len(ot.sequence)} aa)")
        print(f"           {ot.concern}")

    rule("3. SCAFFOLDS")
    for sc in list_scaffolds():
        print(sc.describe())
        print()

    rule("4. EVALUATOR STATUS")
    target_spec = target
    ev = Evaluator(target_spec, use_plm=use_plm)
    print(ev.status())

    ref = ev.reference
    if ref:
        rule("5. REFERENCE MOLECULE (the anchor for all percentages)")
        print(f"{ref.name}: {len(ref.sequence)} aa")
        print(f"  {ref.description}")
        print(f"\n  sequence: {ref.sequence}")
        print(f"\n  composite: {ref.score_100:.2f}/100")
        print("  by group:")
        for g, v in ref.group_scores.items():
            print(f"    {g:<16} {v * 100:6.2f}")
        print()
        print(ref.card.report())

    rule("6. SCORING THE UNMODIFIED SCAFFOLD (negative control)")
    sc = get_scaffold("affibody_Z")
    card = ev.evaluate(sc.framework, sc)
    print("The wild-type Z domain has no reason to bind CLEC12A. It should")
    print("score poorly on recognition. If it does not, the scorers are broken.")
    print()
    print(card.report())
    groups = ev.group_scores(card)
    print("\n  by group:")
    for g, v in groups.items():
        print(f"    {g:<16} {v * 100:6.2f}")
    pct = ev.percent_of_reference(card)
    rpct = ev.recognition_percent_of_reference(card)
    print(f"\n  composite as % of reference : "
          f"{'n/a' if pct is None else f'{pct:.1f}%'}")
    print(f"  recognition as % of reference: "
          f"{'n/a' if rpct is None else f'{rpct:.1f}%'}")

    rule("7. SANITY CONTROLS")
    controls = {
        "poly-alanine (should fail stability/recognition)": "A" * 58,
        "poly-glutamate (charge extreme)": "E" * 58,
        "poly-arginine (same charge as epitope: should repel)": "R" * 58,
        "wild-type Z domain": sc.framework,
    }
    rows = []
    for label, seq in controls.items():
        c = ev.evaluate(seq, sc)
        g = ev.group_scores(c)
        rows.append((label, c.score_100, g.get("recognition", 0) * 100,
                     g.get("selectivity", 0) * 100,
                     g.get("developability", 0) * 100))
    print(f"  {'control':<52} {'total':>6} {'recog':>6} {'selec':>6} {'devel':>6}")
    for label, tot, rec, sel, dev in rows:
        print(f"  {label:<52} {tot:6.2f} {rec:6.2f} {sel:6.2f} {dev:6.2f}")

    print("\nExpectations: poly-Arg should score worst on recognition (it has")
    print("the same charge as the epitope). Poly-Ala should score badly on")
    print("recognition too. If the ordering looks wrong, the scorers need work.")

    rule("8. CALIBRATION AGAINST RANDOM LIBRARY MEMBERS")
    print("A score means nothing without a null distribution. This draws 200")
    print("random members of the same 13-position library and scores them, then")
    print("checks where the real antibody falls in that distribution. If the")
    print("real binder does not stand out, the scoring function is weak and")
    print("every percentile it reports should be discounted.\n")
    cal = ev.calibrate(sc, n_decoys=200)
    print(cal.report())

    print("\n  decoy recognition percentiles:")
    for p, v in cal.recognition.percentiles.items():
        print(f"    {p:>3}th  {v:6.2f}")

    print("\n  where the controls fall:")
    for label, seq in controls.items():
        c = ev.evaluate(seq, sc)
        rec = ev.group_scores(c, exclude_self_referential=True).get(
            "recognition", 0.0) * 100
        interp = cal.interpret(c, rec)
        print(f"    {label[:46]:<46} recognition {rec:6.2f} "
              f"= {interp['recognition_percentile']:5.1f}th pct "
              f"(z {interp['recognition_z']:+.2f})")

    rule("9. STRUCTURE FILES PRESENT")
    for p in sorted((ROOT / "data" / "structures").glob("*.pdb")):
        print(f"  {p.name:<12} {p.stat().st_size:>10,} bytes")

    print("\nSmoke test complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
