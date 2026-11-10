"""Scoring framework.

Design principle: every number a candidate receives must be traceable to a
named, published criterion and must carry an explanation in plain language.
A score with no explanation is not allowed in this pipeline -- that is what
makes the final percentage interpretable rather than decorative.

Each ScoreComponent reports:
  raw        - the measured physical quantity, in its own units
  normalised - that quantity mapped to 0..1, where 1 is "as good as we ask for"
  weight     - how much it contributes to the composite
  detail     - what was measured and why it matters
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class ScoreComponent:
    name: str
    raw: float
    normalised: float
    weight: float
    detail: str
    units: str = ""
    # What this component genuinely cannot tell us. Surfaced in the report.
    caveat: str = ""
    # False when the component is constant across the library being searched.
    # It is still measured and still displayed -- it remains a real check on
    # the molecule -- but it is excluded from the ranking composite, because a
    # value that never varies cannot distinguish one candidate from another.
    # Leaving such components in was inflating every composite by a near
    # constant amount and compressing the useful range.
    ranked: bool = True
    note: str = ""
    # Non-empty when this component says the candidate does the OPPOSITE of
    # its job, as opposed to doing it badly. The two are different claims and
    # a weighted mean cannot tell them apart: it treats every axis as a
    # trade-off, so a candidate buys back a failed one with good scores
    # elsewhere. That is right for most of these -- aromatic content really
    # does trade against solubility -- and wrong for the few that have a
    # measured zero with a meaning.
    #
    # The one that has it today is paralog selectivity. A negative margin
    # means the face complements the nearest relative MORE than the target,
    # and no amount of manufacturability makes that the right molecule to
    # build. Measured: for TNF-alpha the search found 25 candidates with a
    # positive margin and ranked a -0.218 one above all of them, because the
    # 3.7 composite points selectivity can swing were outvoted by the other
    # nine axes. Raising its weight would have distorted the ranking among
    # the good candidates too; the problem is not that the axis is
    # underweighted, it is that disqualification is not a weight.
    disqualifying: str = ""

    def __post_init__(self) -> None:
        self.normalised = max(0.0, min(1.0, float(self.normalised)))

    @property
    def contribution(self) -> float:
        return self.normalised * self.weight if self.ranked else 0.0

    @property
    def effective_weight(self) -> float:
        return self.weight if self.ranked else 0.0

    def __str__(self) -> str:
        bar = "#" * int(round(self.normalised * 20)) + "." * (
            20 - int(round(self.normalised * 20)))
        unit = f" {self.units}" if self.units else ""
        return (f"  {self.name:<26} [{bar}] {self.normalised * 100:5.1f}%   "
                f"raw={self.raw:.3g}{unit}")


@dataclass
class ScoreCard:
    """The full assessment of one candidate sequence."""

    sequence: str
    components: list[ScoreComponent] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)

    def add(self, component: ScoreComponent) -> None:
        self.components.append(component)

    @property
    def total_weight(self) -> float:
        return sum(c.effective_weight for c in self.components)

    @property
    def unranked(self) -> list[ScoreComponent]:
        """Measured and reported, but constant across this library."""
        return [c for c in self.components if not c.ranked]

    @property
    def score(self) -> float:
        """Weighted composite in 0..1. Not a probability of anything."""
        tw = self.total_weight
        if tw <= 0:
            return 0.0
        return sum(c.contribution for c in self.components) / tw

    @property
    def score_100(self) -> float:
        return self.score * 100.0

    @property
    def disqualified(self) -> list[ScoreComponent]:
        """Components saying this candidate does the opposite of its job."""
        return [c for c in self.components if c.disqualifying]

    @property
    def rank_key(self) -> tuple[int, float]:
        """How candidates are ORDERED, which is not the same as how they score.

        The composite stays an honest weighted mean of everything measured --
        a counter-selective TNF-alpha design really does score 82.6 on the
        nine axes that are not selectivity, and hiding that would be its own
        kind of dishonesty. What it must not do is come first.

        So ordering is lexicographic: anything disqualified sorts below
        everything that is not, and within each group the composite decides.
        The composite still drives the search inside the disqualified group,
        which is why the counter-selective band is a band and not a flat
        zero -- it is the gradient a candidate climbs to get out.

        If nothing qualifies, every candidate is in one group and this
        degenerates to ranking by composite, exactly as before. That case is
        reported rather than hidden: it means the search never found a face
        that prefers the target to its nearest relative.
        """
        return (0 if self.disqualified else 1, self.score)

    def component(self, name: str) -> ScoreComponent | None:
        for c in self.components:
            if c.name == name:
                return c
        return None

    def weakest(self, n: int = 3) -> list[ScoreComponent]:
        """The components dragging this candidate down -- i.e. what to fix."""
        return sorted(self.components, key=lambda c: c.normalised)[:n]

    def strongest(self, n: int = 3) -> list[ScoreComponent]:
        return sorted(self.components, key=lambda c: -c.normalised)[:n]

    def to_dict(self) -> dict:
        return {
            "sequence": self.sequence,
            "length": len(self.sequence),
            "score": round(self.score, 6),
            "score_100": round(self.score_100, 2),
            "flags": list(self.flags),
            "metadata": dict(self.metadata),
            # Why this candidate cannot be ranked first, if it cannot. A list,
            # because disqualification is a property of specific components
            # and the reader is owed the name of the one responsible.
            "disqualifying": [
                {"name": c.name, "why": c.disqualifying}
                for c in self.disqualified
            ],
            "components": [
                {
                    "name": c.name,
                    "raw": None if c.raw != c.raw else round(float(c.raw), 6),
                    "normalised": round(c.normalised, 6),
                    "weight": c.weight,
                    "units": c.units,
                    "detail": c.detail,
                    "caveat": c.caveat,
                    "ranked": c.ranked,
                    "note": c.note,
                    "disqualifying": c.disqualifying,
                }
                for c in self.components
            ],
        }

    def report(self) -> str:
        ranked = [c for c in self.components if c.ranked]
        unranked = self.unranked
        lines = [
            f"Candidate ({len(self.sequence)} aa): {self.sequence}",
            f"Composite design score: {self.score_100:.1f} / 100"
            f"  (over {len(ranked)} ranking components)",
            "",
        ]
        for c in sorted(ranked, key=lambda c: -c.weight):
            lines.append(str(c))
        if unranked:
            lines.append("")
            lines.append("Measured but not ranked -- constant across this "
                         "library, so they cannot separate candidates:")
            for c in sorted(unranked, key=lambda c: c.name):
                lines.append(f"{str(c)}  [not ranked]")
        if self.flags:
            lines.append("")
            lines.append("Flags:")
            lines.extend(f"  ! {f}" for f in self.flags)
        return "\n".join(lines)


class Scorer:
    """Base class. Subclasses implement score() and declare a weight."""

    name: str = "unnamed"
    weight: float = 1.0

    def score(self, sequence: str, context: dict) -> ScoreComponent:
        raise NotImplementedError

    def available(self) -> bool:
        """False when an optional dependency or model is not present."""
        return True


# --- Normalisation helpers --------------------------------------------------


def window_score(value: float, lo: float, hi: float,
                 soft: float = 0.0) -> float:
    """1.0 inside [lo, hi], tapering to 0 over `soft` units outside it.

    Used wherever the criterion is 'should be in this range', which is most
    developability properties -- too little and too much are both bad.
    """
    if lo <= value <= hi:
        return 1.0
    if soft <= 0:
        return 0.0
    if value < lo:
        return max(0.0, 1.0 - (lo - value) / soft)
    return max(0.0, 1.0 - (value - hi) / soft)


def target_score(value: float, ideal: float, tolerance: float) -> float:
    """1.0 at `ideal`, falling linearly to 0 at `ideal +/- tolerance`."""
    if tolerance <= 0:
        return 1.0 if value == ideal else 0.0
    return max(0.0, 1.0 - abs(value - ideal) / tolerance)


def ramp_score(value: float, zero_at: float, one_at: float) -> float:
    """Linear ramp: 0 at `zero_at`, 1 at `one_at`. Works in either direction."""
    if zero_at == one_at:
        return 1.0 if value >= one_at else 0.0
    t = (value - zero_at) / (one_at - zero_at)
    return max(0.0, min(1.0, t))


def penalty_score(count: int, free: int = 0, each: float = 0.25) -> float:
    """1.0 while count <= free, then dropping by `each` per extra occurrence."""
    excess = max(0, count - free)
    return max(0.0, 1.0 - excess * each)
