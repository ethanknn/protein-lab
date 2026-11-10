"""Run a full design cycle and build the dashboard.

    python run.py                          # default: affibody scaffold, ESM-2 on
    python run.py --no-plm                 # skip the language model (much faster)
    python run.py --generations 60 --pop 96
    python run.py --scaffold beta_hairpin_SS
    python run.py --decoys 400 --seed 11

Outputs land in out/:
    out/run_<scaffold>.json     machine-readable record of the whole search
    out/dashboard.html          the visual report -- open this in a browser
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from plab.design.loop import DesignLoop, LoopConfig  # noqa: E402
from plab.design.scaffolds import get_scaffold  # noqa: E402
from plab.report.dashboard import build_dashboard  # noqa: E402
from plab.scoring.composite import Evaluator  # noqa: E402
from plab.targets import load_target  # noqa: E402

STRUCT_DIR = ROOT / "data" / "structures"
OUT_DIR = ROOT / "out"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--target", default="CLEC12A")
    p.add_argument("--scaffold", default="affibody_Z")
    p.add_argument("--generations", type=int, default=30)
    p.add_argument("--pop", type=int, default=48)
    p.add_argument("--mutation-rate", type=float, default=0.18)
    p.add_argument("--patience", type=int, default=8)
    p.add_argument("--decoys", type=int, default=200)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--top", type=int, default=10)
    p.add_argument("--no-plm", action="store_true",
                   help="disable the ESM-2 language model scorer")
    p.add_argument("--target-score", type=float, default=None,
                   help="stop early at this composite score (0..1)")
    p.add_argument("--out", default=None, help="dashboard path")
    p.add_argument("--tag", default=None,
                   help="suffix for output filenames, so concurrent runs "
                        "do not overwrite each other's records")
    p.add_argument("--artifact", action="store_true",
                   help="emit the dashboard as a fragment for publishing as "
                        "an Artifact instead of a standalone file")
    return p.parse_args()


def rule(t: str) -> None:
    print(f"\n{'=' * 72}\n{t}\n{'=' * 72}")


def main() -> int:
    args = parse_args()
    t_start = time.time()

    rule("TARGET")
    target = load_target(args.target)
    print(target.design_brief())

    rule("SCAFFOLD")
    scaffold = get_scaffold(args.scaffold)
    print(scaffold.describe())

    rule("EVALUATOR")
    if not args.no_plm:
        print("Loading ESM-2. The first run downloads the weights (~2.5 GB) --")
        print("this is a one-time cost, cached in ~/.cache/torch/hub.\n")
    evaluator = Evaluator(target, use_plm=not args.no_plm,
                          struct_dir=STRUCT_DIR)
    print(evaluator.status())

    ref = evaluator.reference
    if ref is None:
        print("\nWARNING: no reference molecule could be built, so all")
        print("percentages relative to a real binder are unavailable. The")
        print("absolute scores still work but are much harder to interpret.")

    rule(f"CALIBRATION ({args.decoys} random library members)")
    print("Establishing the null distribution before optimising, so the")
    print("optimiser's output can be judged against chance.\n")
    calibration = evaluator.calibrate(scaffold, n_decoys=args.decoys,
                                      seed=args.seed + 1000)
    print(calibration.report())

    rule("DESIGN LOOP")
    config = LoopConfig(
        population_size=args.pop,
        generations=args.generations,
        mutation_rate=args.mutation_rate,
        patience=args.patience,
        target_score=args.target_score,
        seed=args.seed,
    )
    loop = DesignLoop(evaluator, scaffold, config)
    loop.run(verbose=True)

    rule("RESULT")
    bestI = loop.best()
    groups = evaluator.group_scores(bestI.card, exclude_self_referential=True)
    rec = groups.get("recognition", 0.0) * 100
    interp = calibration.interpret(bestI.card, rec)

    print(f"Best candidate (generation {bestI.generation}, via {bestI.origin}):")
    print(f"  full sequence : {bestI.sequence}")
    print(f"  designed face : {bestI.card.metadata.get('designed_residues','')}")
    print(f"  annotated     : {scaffold.annotate(bestI.sequence)}")
    print()
    print(f"  composite score        : {bestI.score_100:.2f} / 100")
    print(f"  recognition subscore   : {rec:.2f} / 100")
    print(f"  percentile vs random   : {interp['recognition_percentile']:.1f}th "
          f"(z = {interp['recognition_z']:+.2f})")
    pct = evaluator.recognition_percent_of_reference(bestI.card)
    if pct is not None:
        print(f"  vs the real antibody   : {pct:.1f}% of its recognition score")
    print(f"  scoring function power : {calibration.discriminating_power}")
    print()
    print(bestI.card.report())

    print("\n--- what this does and does not mean ---")
    print("This candidate scores well on criteria a laptop can compute. That")
    print("makes it a reasonable thing to synthesise first, and nothing more.")
    print("It has not been shown to bind CLEC12A. Establishing that needs")
    print("expression, a binding assay, then primary AML cells, then animal")
    print("work, then trials. Most candidates that look this good still fail.")

    rule("OUTPUTS")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    suffix = f"_{args.tag}" if args.tag else ""
    json_path = loop.save(OUT_DIR / f"run_{scaffold.name}{suffix}.json",
                          top_n=25)
    print(f"  run record : {json_path}")

    dash = Path(args.out) if args.out else (
        OUT_DIR / f"dashboard{suffix}.html")
    dash = build_dashboard(loop, evaluator, calibration, dash,
                           STRUCT_DIR, top_n=args.top,
                           artifact_mode=args.artifact)
    size_mb = dash.stat().st_size / 1e6
    print(f"  dashboard  : {dash}  ({size_mb:.1f} MB)")

    # In artifact mode also write the standalone file, so one expensive run
    # produces both the publishable fragment and something openable on disk.
    if args.artifact:
        from plab.report.dashboard import _wrap_local
        local = dash.with_name(dash.stem + "_standalone.html")
        local.write_text(
            _wrap_local(dash.read_text(encoding="utf-8")), encoding="utf-8")
        print(f"  standalone : {local}")
    print(f"\nTotal wall time: {time.time() - t_start:.1f}s")
    print(f"\nOpen the dashboard in a browser:\n  {dash}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
