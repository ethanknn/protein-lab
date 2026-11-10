"""Calibration: what does a score of 72 actually mean?

A raw score is uninterpretable on its own. This module supplies the two
things that give it meaning:

1. A DECOY PANEL -- random draws from the same scaffold library. A candidate
   is reported as a percentile against these. "Better than 98% of random
   library members" is a claim with content; "scored 72" is not. This is the
   same logic as an enrichment score in a display campaign.

2. DISCRIMINATING POWER of the scoring function itself. We know one molecule
   that genuinely engages this epitope: the 50C1 antibody. If the scoring
   function cannot separate that molecule from random decoys, then the
   function is weak and every percentile it produces should be discounted.
   That diagnostic is computed and reported rather than hidden, because a
   scoring function that cannot tell a real binder from noise is the single
   most important thing a user of this tool should know.

The honest expectation: composition-based scoring has genuinely limited
discriminating power for protein-protein recognition, because real
specificity is geometric. This module exists to quantify that limit instead
of pretending it away.
"""

from __future__ import annotations

import math
import random
import statistics
from dataclasses import dataclass, field

from ..design.scaffolds import Scaffold
from .base import ScoreCard


@dataclass
class DecoyStats:
    n: int
    mean: float
    stdev: float
    minimum: float
    maximum: float
    percentiles: dict[int, float] = field(default_factory=dict)
    values: list[float] = field(default_factory=list)

    def percentile_of(self, value: float) -> float:
        """Percentage of decoys this value beats."""
        if not self.values:
            return float("nan")
        below = sum(1 for v in self.values if v < value)
        return 100.0 * below / len(self.values)

    def z_of(self, value: float) -> float:
        if self.stdev <= 1e-9:
            return float("nan")
        return (value - self.mean) / self.stdev


@dataclass
class Calibration:
    """Everything needed to interpret a score for one scaffold + target."""

    scaffold_name: str
    n_decoys: int
    composite: DecoyStats
    recognition: DecoyStats
    reference_composite: float | None
    reference_recognition: float | None
    # z-score of the reference against the decoy distribution. This IS the
    # discriminating power of the scoring function.
    reference_z_composite: float | None
    reference_z_recognition: float | None
    notes: list[str] = field(default_factory=list)

    @property
    def discriminating_power(self) -> str:
        """Plain-language verdict on whether the score can be trusted."""
        z = self.reference_z_recognition
        if z is None or z != z:
            return "unknown (no reference available)"
        if z >= 3.0:
            return f"strong (reference sits {z:.1f} sd above random)"
        if z >= 2.0:
            return f"moderate (reference sits {z:.1f} sd above random)"
        if z >= 1.0:
            return f"weak (reference sits only {z:.1f} sd above random)"
        return (f"poor (reference sits {z:.1f} sd from random -- this scoring "
                f"function cannot reliably tell a real binder from noise)")

    @property
    def trust_discount(self) -> float:
        """How much to discount percentiles, 0..1. Used in the report."""
        z = self.reference_z_recognition
        if z is None or z != z:
            return 0.5
        return max(0.1, min(1.0, z / 3.0))

    def interpret(self, card: ScoreCard, recognition: float) -> dict:
        return {
            "composite_percentile": self.composite.percentile_of(card.score_100),
            "composite_z": self.composite.z_of(card.score_100),
            "recognition_percentile": self.recognition.percentile_of(recognition),
            "recognition_z": self.recognition.z_of(recognition),
            "decoy_mean_composite": self.composite.mean,
            "decoy_mean_recognition": self.recognition.mean,
            "discriminating_power": self.discriminating_power,
            "trust_discount": self.trust_discount,
        }

    def to_dict(self) -> dict:
        def pack(s: DecoyStats) -> dict:
            return {
                "n": s.n, "mean": round(s.mean, 4), "stdev": round(s.stdev, 4),
                "min": round(s.minimum, 4), "max": round(s.maximum, 4),
                "percentiles": {k: round(v, 4) for k, v in s.percentiles.items()},
                "values": [round(v, 4) for v in s.values],
            }

        return {
            "scaffold": self.scaffold_name,
            "n_decoys": self.n_decoys,
            "composite": pack(self.composite),
            "recognition": pack(self.recognition),
            "reference_composite": self.reference_composite,
            "reference_recognition": self.reference_recognition,
            "reference_z_composite": self.reference_z_composite,
            "reference_z_recognition": self.reference_z_recognition,
            "discriminating_power": self.discriminating_power,
            "trust_discount": round(self.trust_discount, 3),
            "notes": self.notes,
        }

    def report(self) -> str:
        lines = [
            f"Calibration for {self.scaffold_name} ({self.n_decoys} random decoys)",
            f"  decoy composite   : {self.composite.mean:.2f} "
            f"+/- {self.composite.stdev:.2f} "
            f"(range {self.composite.minimum:.1f}-{self.composite.maximum:.1f})",
            f"  decoy recognition : {self.recognition.mean:.2f} "
            f"+/- {self.recognition.stdev:.2f} "
            f"(range {self.recognition.minimum:.1f}-{self.recognition.maximum:.1f})",
        ]
        if self.reference_recognition is not None:
            lines.append(
                f"  reference binder  : composite "
                f"{self.reference_composite:.2f}, recognition "
                f"{self.reference_recognition:.2f}")
        lines.append(f"  DISCRIMINATING POWER: {self.discriminating_power}")
        for n in self.notes:
            lines.append(f"  note: {n}")
        return "\n".join(lines)


def _stats(values: list[float]) -> DecoyStats:
    if not values:
        return DecoyStats(0, 0.0, 0.0, 0.0, 0.0, {}, [])
    srt = sorted(values)

    def pct(p: int) -> float:
        if len(srt) == 1:
            return srt[0]
        idx = (p / 100.0) * (len(srt) - 1)
        lo, hi = int(math.floor(idx)), int(math.ceil(idx))
        if lo == hi:
            return srt[lo]
        return srt[lo] + (srt[hi] - srt[lo]) * (idx - lo)

    return DecoyStats(
        n=len(values),
        mean=statistics.fmean(values),
        stdev=statistics.stdev(values) if len(values) > 1 else 0.0,
        minimum=srt[0],
        maximum=srt[-1],
        percentiles={p: pct(p) for p in (5, 25, 50, 75, 95, 99)},
        values=values,
    )


def build_calibration(evaluator, scaffold: Scaffold, *, n_decoys: int = 200,
                      seed: int = 1234, verbose: bool = False) -> Calibration:
    """Score random library members to establish the null distribution.

    Decoys are random assignments at exactly the designed positions, i.e.
    genuine random members of the library this loop searches. That is the
    correct null: it asks whether the optimiser found something, not whether
    the scaffold happens to score well.
    """
    from ..aminoacids import AA20

    rng = random.Random(seed)
    excluded = set("C")  # keep the cysteine-free scaffold cysteine-free
    ss = scaffold.secondary_structure

    comp_vals: list[float] = []
    recog_vals: list[float] = []
    notes: list[str] = []

    for _ in range(n_decoys):
        assignment = {}
        for p in scaffold.variable_positions:
            pool = [a for a in AA20 if a not in excluded
                    and not (p < len(ss) and ss[p] == "H" and a in "PG")]
            assignment[p] = rng.choice(pool)
        seq = scaffold.apply(scaffold.framework, assignment)
        card = evaluator.evaluate(seq, scaffold)
        comp_vals.append(card.score_100)
        groups = evaluator.group_scores(card)
        recog_vals.append(groups.get("recognition", 0.0) * 100.0)

    composite = _stats(comp_vals)
    recognition = _stats(recog_vals)

    ref = evaluator.reference
    ref_comp = ref_recog = ref_z_comp = ref_z_recog = None
    if ref is not None:
        ref_comp = ref.score_100
        ref_recog = ref.group_scores.get("recognition", 0.0) * 100.0
        ref_z_comp = composite.z_of(ref_comp)
        ref_z_recog = recognition.z_of(ref_recog)
        notes.append(
            "The reference is an antibody scFv and the decoys are 58-residue "
            "affibodies, so the composite comparison is across formats. The "
            "recognition comparison is the meaningful one."
        )

    if recognition.stdev < 1.0:
        notes.append(
            "Decoy recognition scores barely vary, which means the recognition "
            "scorers are not discriminating between library members. Treat "
            "percentiles from this run as close to meaningless."
        )

    return Calibration(
        scaffold_name=scaffold.name,
        n_decoys=n_decoys,
        composite=composite,
        recognition=recognition,
        reference_composite=ref_comp,
        reference_recognition=ref_recog,
        reference_z_composite=ref_z_comp,
        reference_z_recognition=ref_z_recog,
        notes=notes,
    )
