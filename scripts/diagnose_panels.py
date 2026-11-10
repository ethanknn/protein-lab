"""Are the derived selectivity panels measuring anything?

A panel that produces a number is not the same as a panel that produces a
meaningful number. This project has already shipped a selectivity score
pinned at exactly 50% for every design, and the only reason that was caught
was running controls whose answers were known in advance. Forty-one new
panels built automatically deserve the same treatment before anyone reads a
number off them.

Four controls per target, each with an answer known before it runs:

  self          run the target's own sequence through the mapping used on
                every panel member; it must come back as exactly the epitope.
                Anything else means the machinery is broken.

  reference     the real binder's own contact residues. A molecule that
                demonstrably binds this target and not its relatives should
                not score badly. If it does, the panel or the kernel is
                wrong, not the antibody.

  random        random faces. These should sit low, and their spread is what
                tells you the metric can vary at all.

  anti-target   a face complementary to the closest RELATIVE instead of the
                target. This should score clearly worse than the reference.
                It is the sharpest control: it has the same composition
                budget and differs only in what it is aimed at.

The summary flags any target where the metric is saturated -- every face
scoring the same -- because that is the failure this codebase keeps hitting.

Run: python scripts/diagnose_panels.py [KEY ...]
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
from plab.scoring.selectivity import (  # noqa: E402
    CEILING_NOISE_FACTOR, INDISTINGUISHABLE_MARGIN,
    achievable_margin_range, build_paralog_epitopes, complementarity_kernel,
    margin_to_score,
)
from plab.targets import available_targets, load_off_targets, load_target  # noqa: E402

N_RANDOM = 200
FLAT_SD = 0.02


def best_face_for(epitope: str, length: int) -> str:
    """The composition this kernel likes most against a given surface."""
    best, best_v = "S" * length, -9.9
    for n_acid in range(length + 1):
        for n_arom in range(length - n_acid + 1):
            for n_base in range(length - n_acid - n_arom + 1):
                rest = length - n_acid - n_arom - n_base
                face = "D" * n_acid + "Y" * n_arom + "K" * n_base + "S" * rest
                v = complementarity_kernel(face, epitope)
                if v > best_v:
                    best, best_v = face, v
    return best


def report(key: str, scaffold) -> dict | None:
    t = load_target(key)
    panel = build_paralog_epitopes(t, load_off_targets(t))
    if not panel:
        return {"key": key, "gene": t.gene, "state": t.panel_state,
                "n_panel": 0}

    length = len(scaffold.variable_positions)
    epi = t.epitope_sequence
    lo, hi = achievable_margin_range(epi, panel, length)

    def margin(face: str) -> float:
        on = complementarity_kernel(face, epi)
        off = max(complementarity_kernel(face, p) for p in panel.values())
        return on - off

    def pct(face: str) -> float:
        # The scorer's own ramp, imported rather than reimplemented. This
        # control used to compute the percentage itself, which meant it was
        # checking a number nobody reads: it went on reporting the reference
        # at 49 while the scorer, on the same margin of -0.268, was free to
        # report anything at all.
        return 100.0 * margin_to_score(margin(face), lo, hi)

    # Self control. Run the target's OWN sequence through the same mapping
    # machinery the panel members go through: it must come back as exactly
    # the epitope. Subtracting a kernel from itself would always give zero
    # and prove nothing; this actually exercises the alignment, the window
    # and the numbering offset, which is where the mistakes live.
    from plab.scoring.selectivity import epitope_window, map_epitope_onto
    win = epitope_window(t)
    remapped = map_epitope_onto(*win, t.full_sequence) if win else None
    self_ok = remapped is not None and remapped.surface == epi

    # The real antibody's contact residues, whole. An earlier version of this
    # control cut them to the scaffold's 13 positions with [:length], which
    # for a 28-residue paratope keeps a contiguous prefix and throws away
    # whatever is past it -- for CD38 that discarded six of isatuximab's
    # nine tyrosines. Since the kernel reads composition, slicing the
    # molecule changed which molecule was being tested: adalimumab's face
    # scored 11, 46 or 31 depending only on how it was cut. The kernel is
    # length-invariant now, so the reference needs no cutting.
    ref = t.paratope_sequence or "S" * length
    rng = random.Random(11)
    alphabet = [a for a in AA20 if a != "C"]
    randoms = ["".join(rng.choice(alphabet) for _ in range(length))
               for _ in range(N_RANDOM)]
    rand_pcts = [pct(f) for f in randoms]

    worst_relative = max(panel, key=lambda g: complementarity_kernel(ref, panel[g]))
    anti = best_face_for(panel[worst_relative], length)
    face_len = str(length)
    noise = (t.panel_ceiling_noise or {}).get(face_len)
    stored_ceiling = (t.panel_ceiling or {}).get(face_len, hi)

    return {
        "key": key, "gene": t.gene, "state": t.panel_state,
        "n_panel": len(panel),
        "closest": worst_relative,
        "self_ok": self_ok,
        "reference": pct(ref),
        # The RAW margin, not just the percentage. A margin below zero
        # means the face prefers the relative to the target, and that
        # sign was invisible for every run this script has ever made:
        # nine of fifteen references were negative while reading 40-60.
        "ref_margin": margin(ref),
        "anti_margin": margin(anti),
        "random_mean": statistics.fmean(rand_pcts),
        "random_sd": statistics.pstdev(rand_pcts),
        "anti": pct(anti),
        "scale": (lo, hi),
        "ceiling": stored_ceiling,
        # No face can separate this target from its closest relative. The
        # score is withheld for these, so the ref/anti comparison below does
        # not apply to them and must not be reported as a failure.
        #
        # The scorer's own test, not a reimplementation of it: both the
        # absolute cutoff and the ceiling's measured error bar, against the
        # stored ceiling rather than the live scale, because
        # achievable_margin_range substitutes a usable value when the sweep
        # finds no positive margin at all.
        "indistinguishable": stored_ceiling <= max(
            INDISTINGUISHABLE_MARGIN, CEILING_NOISE_FACTOR * (noise or 0.0)),
        "ceiling_noise": noise,
    }


def main() -> int:
    keys = sys.argv[1:] or sorted(available_targets())
    scaffold = get_scaffold("affibody_Z")
    rows, skipped = [], []
    for k in keys:
        try:
            r = report(k, scaffold)
        except Exception as exc:
            print(f"  {k}: ERROR {type(exc).__name__}: {exc}")
            continue
        if r is None or not r["n_panel"]:
            skipped.append(r)
            continue
        rows.append(r)

    print(f"{'target':<22}{'n':>3}{'closest':>10}{'ref':>6}{'ref_m':>8}"
          f"{'rand':>6}{'anti':>6}{'ceiling':>9}  verdict")
    print("-" * 88)
    problems = []
    for r in sorted(rows, key=lambda x: x["key"]):
        flags = []
        if not r["self_ok"]:
            flags.append("SELF-MAPPING BROKEN")
        if r["indistinguishable"]:
            verdict = "indistinguishable - score withheld"
        else:
            if r["random_sd"] < FLAT_SD:
                flags.append("SATURATED")
            if r["reference"] <= r["anti"]:
                flags.append("ref below anti-target")
            # The sign, reported separately from the percentage. A negative
            # reference margin is not "the antibody scored low": it is this
            # metric claiming a licensed drug prefers the off-target, and it
            # is the single most useful thing this script can say about
            # whether the number above it means anything.
            if r["ref_margin"] < 0:
                flags.append("reference margin NEGATIVE")
            verdict = "; ".join(flags) if flags else "ok"
        if flags:
            problems.append((r["key"], verdict))
        print(f"{r['key']:<22}{r['n_panel']:>3}{r['closest'][:9]:>10}"
              f"{r['reference']:>6.0f}{r['ref_margin']:>+8.3f}"
              f"{r['random_mean']:>6.0f}{r['anti']:>6.0f}"
              f"{r['ceiling']:>+9.3f}  {verdict}")

    print("-" * 88)
    print(f"{len(rows)} targets with a panel, {len(skipped)} without")
    if rows:
        sr = [r for r in rows if not r["indistinguishable"]] or rows
        print(f"reference mean {statistics.fmean(r['reference'] for r in sr):.0f}, "
              f"anti-target mean {statistics.fmean(r['anti'] for r in sr):.0f}, "
              f"random mean {statistics.fmean(r['random_mean'] for r in sr):.0f}")
        scored_rows = [r for r in rows if not r["indistinguishable"]]
        beats = sum(1 for r in scored_rows if r["reference"] > r["anti"])
        print(f"reference beats the anti-target face in "
              f"{beats}/{len(scored_rows)} scored targets "
              f"({len(rows) - len(scored_rows)} withheld as indistinguishable)")
        pos = sum(1 for r in scored_rows if r["ref_margin"] > 0)
        print(f"the reference's own margin is positive -- it prefers the "
              f"target to its nearest relative -- in {pos}/{len(scored_rows)}")
    if problems:
        print(f"\n{len(problems)} target(s) need a look:")
        for k, v in problems:
            print(f"  {k}: {v}")
    states: dict[str, int] = {}
    for r in skipped:
        if r:
            states[r["state"]] = states.get(r["state"], 0) + 1
    if states:
        print("\nno panel, by reason:")
        for s, n in sorted(states.items(), key=lambda kv: -kv[1]):
            print(f"  {s:<24}{n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
