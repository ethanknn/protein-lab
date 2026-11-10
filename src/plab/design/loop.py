"""The design loop: propose, score, select, repeat.

This is an evolutionary search over the scaffold's designed positions only.
The framework is never touched, so every individual in the population is a
molecule that plausibly folds.

Why evolutionary rather than gradient-based: the objective is a weighted sum
of thirteen non-differentiable criteria, several of which are discrete counts.
A population method handles that directly, keeps a diverse set of solutions
rather than one, and produces a legible history -- which matters here, because
being able to see why a candidate won is part of the point.

Honesty note on convergence: the loop stops when the composite score plateaus
or hits the configured target. That means "the search found nothing better",
not "the molecule works". Those are different claims and the report keeps
them separate.
"""

from __future__ import annotations

import json
import random
import statistics
import time
from dataclasses import dataclass, field
from pathlib import Path

from ..aminoacids import AA20
from ..scoring.base import ScoreCard
from .scaffolds import Scaffold

# Residues that break an alpha helix; never proposed at a helical position.
HELIX_BREAKERS = set("PG")
# Cysteine is excluded everywhere: a free thiol in a cysteine-free scaffold
# creates dimers and scrambling. Disulfide scaffolds freeze their own cysteines.
ALWAYS_EXCLUDED = set("C")


@dataclass
class LoopConfig:
    population_size: int = 48
    generations: int = 30
    mutation_rate: float = 0.18      # per designed position, per offspring
    crossover_rate: float = 0.5
    elite_count: int = 4
    tournament_size: int = 3
    # Stop early if the best score stops improving for this many generations.
    patience: int = 8
    min_improvement: float = 0.002
    # Stop early if this composite score (0..1) is reached.
    target_score: float | None = None
    # Fraction of mutations drawn from the brief-informed distribution rather
    # than uniformly. Explicit so the bias is visible, not hidden.
    informed_fraction: float = 0.65
    seed: int | None = 7

    def to_dict(self) -> dict:
        return dict(self.__dict__)


@dataclass
class Individual:
    sequence: str
    card: ScoreCard
    generation: int
    origin: str = "random"           # random | mutation | crossover | elite
    parents: tuple[str, ...] = ()

    @property
    def score(self) -> float:
        return self.card.score

    @property
    def score_100(self) -> float:
        return self.card.score_100

    @property
    def rank_key(self) -> tuple[int, float]:
        """Ordering, not scoring -- see ScoreCard.rank_key.

        Every comparison in this file goes through here rather than through
        .score, so selection pressure points the same way the final ranking
        does. Changing only the final sort would have left the search
        spending its whole budget in the disqualified region and then
        reporting the least bad thing it happened to pass on the way.
        """
        return self.card.rank_key


@dataclass
class Generation:
    index: int
    individuals: list[Individual]
    elapsed_s: float = 0.0
    n_evaluated: int = 0
    n_cached: int = 0

    @property
    def best(self) -> Individual:
        return max(self.individuals, key=lambda i: i.rank_key)

    @property
    def scores(self) -> list[float]:
        return [i.score for i in self.individuals]

    @property
    def mean(self) -> float:
        return statistics.fmean(self.scores) if self.individuals else 0.0

    @property
    def median(self) -> float:
        return statistics.median(self.scores) if self.individuals else 0.0

    @property
    def stdev(self) -> float:
        s = self.scores
        return statistics.stdev(s) if len(s) > 1 else 0.0

    @property
    def diversity(self) -> float:
        """Mean pairwise Hamming distance over designed positions, normalised.

        Watching this collapse is how you know the search has converged rather
        than merely stalled.
        """
        views = [i.card.metadata.get("designed_residues", "")
                 for i in self.individuals]
        views = [v for v in views if v]
        if len(views) < 2:
            return 0.0
        n = len(views[0])
        if n == 0:
            return 0.0
        total = pairs = 0
        for a in range(len(views)):
            for b in range(a + 1, len(views)):
                total += sum(1 for x, y in zip(views[a], views[b]) if x != y)
                pairs += 1
        return (total / pairs) / n if pairs else 0.0

    def summary(self) -> str:
        return (f"gen {self.index:>3}  best {self.best.score_100:6.2f}  "
                f"mean {self.mean * 100:6.2f}  sd {self.stdev * 100:5.2f}  "
                f"diversity {self.diversity:5.2f}  "
                f"({self.n_evaluated} new, {self.n_cached} cached, "
                f"{self.elapsed_s:.1f}s)")


class ProposalDistribution:
    """Which residues to try at which position.

    The informed component is derived from the target's own design brief --
    for an arginine-rich epitope it favours acidic and aromatic residues --
    rather than being a fixed table. The uniform component keeps the search
    from collapsing into whatever the brief suggests, so the loop can still
    discover something the brief did not anticipate.
    """

    def __init__(self, target, scaffold: Scaffold, rng: random.Random) -> None:
        self.rng = rng
        self.scaffold = scaffold
        self.informed = self._build_informed(target)

    @staticmethod
    def _build_informed(target) -> list[str]:
        """A weighted bag of residues, built from the epitope's chemistry."""
        bag: list[str] = []
        n_arg = len(target.arginine_positions)
        q = target.epitope_net_charge

        # Acidic residues to pair the epitope's cations.
        weight_acid = 6 if n_arg >= 2 else 3
        bag += ["D"] * weight_acid + ["E"] * weight_acid
        # Aromatics for cation-pi stacking; tyrosine weighted highest because
        # it adds a hydrogen-bonding hydroxyl, which is why real paratopes
        # are tyrosine-rich.
        bag += ["Y"] * 7 + ["W"] * 3 + ["F"] * 3
        # Polar residues for hydrogen bonding to the epitope's Ser/Thr/Gln.
        bag += ["S", "T", "N", "Q"] * 2
        # A few basic residues to answer the epitope's aspartates, but kept
        # scarce: excess positive charge drives nonspecific binding.
        if q > 0:
            bag += ["K", "R"]
        else:
            bag += ["K"] * 3 + ["R"] * 3
        # Small aliphatics for packing.
        bag += ["A", "L", "I", "V", "M"]
        return bag

    def allowed(self, position: int) -> list[str]:
        ss = self.scaffold.secondary_structure
        excluded = set(ALWAYS_EXCLUDED)
        if position < len(ss) and ss[position] == "H":
            excluded |= HELIX_BREAKERS
        return [a for a in AA20 if a not in excluded]

    def draw(self, position: int, informed_fraction: float,
             exclude: str | None = None) -> str:
        allowed = set(self.allowed(position))
        if exclude:
            allowed.discard(exclude)
        if not allowed:
            return exclude or "A"
        if self.rng.random() < informed_fraction:
            pool = [a for a in self.informed if a in allowed]
            if pool:
                return self.rng.choice(pool)
        return self.rng.choice(sorted(allowed))


class DesignLoop:
    """Runs the search and keeps the full history."""

    def __init__(self, evaluator, scaffold: Scaffold,
                 config: LoopConfig | None = None) -> None:
        self.evaluator = evaluator
        self.scaffold = scaffold
        self.config = config or LoopConfig()
        self.rng = random.Random(self.config.seed)
        self.proposals = ProposalDistribution(
            evaluator.target, scaffold, self.rng)
        self.history: list[Generation] = []
        self._cache: dict[str, ScoreCard] = {}
        self._n_evaluated = 0
        self._n_cached = 0
        self.stop_reason: str = "not started"

    # --- scoring with memoisation ----------------------------------------

    def _score(self, sequence: str) -> ScoreCard:
        hit = self._cache.get(sequence)
        if hit is not None:
            self._n_cached += 1
            return hit
        card = self.evaluator.evaluate(sequence, self.scaffold)
        self._cache[sequence] = card
        self._n_evaluated += 1
        return card

    # --- population operations -------------------------------------------

    def _random_sequence(self) -> str:
        assignment = {
            p: self.proposals.draw(p, self.config.informed_fraction)
            for p in self.scaffold.variable_positions
        }
        return self.scaffold.apply(self.scaffold.framework, assignment)

    def _mutate(self, sequence: str) -> str:
        chars = list(sequence)
        changed = False
        for p in self.scaffold.variable_positions:
            if self.rng.random() < self.config.mutation_rate:
                chars[p] = self.proposals.draw(
                    p, self.config.informed_fraction, exclude=chars[p])
                changed = True
        if not changed:  # guarantee at least one change
            p = self.rng.choice(self.scaffold.variable_positions)
            chars[p] = self.proposals.draw(
                p, self.config.informed_fraction, exclude=chars[p])
        return "".join(chars)

    def _crossover(self, a: str, b: str) -> str:
        """Uniform crossover over the designed positions only."""
        chars = list(a)
        for p in self.scaffold.variable_positions:
            if self.rng.random() < 0.5:
                chars[p] = b[p]
        return "".join(chars)

    def _tournament(self, pop: list[Individual]) -> Individual:
        k = min(self.config.tournament_size, len(pop))
        return max(self.rng.sample(pop, k), key=lambda i: i.rank_key)

    # --- the loop ---------------------------------------------------------

    def run(self, *, verbose: bool = True, on_generation=None) -> Generation:
        cfg = self.config
        t_start = time.time()

        # Generation 0: random draws from the proposal distribution.
        t0 = time.time()
        before_e, before_c = self._n_evaluated, self._n_cached
        seen: set[str] = set()
        pop: list[Individual] = []
        attempts = 0
        while len(pop) < cfg.population_size and attempts < cfg.population_size * 20:
            attempts += 1
            seq = self._random_sequence()
            if seq in seen:
                continue
            seen.add(seq)
            pop.append(Individual(seq, self._score(seq), 0, origin="random"))

        gen = Generation(0, pop, time.time() - t0,
                         self._n_evaluated - before_e, self._n_cached - before_c)
        self.history.append(gen)
        if verbose:
            print(gen.summary())
        if on_generation:
            on_generation(gen)

        best_key = gen.best.rank_key
        stagnant = 0

        for g in range(1, cfg.generations + 1):
            t0 = time.time()
            before_e, before_c = self._n_evaluated, self._n_cached

            ranked = sorted(pop, key=lambda i: i.rank_key,
                            reverse=True)
            elites = [
                Individual(i.sequence, i.card, g, origin="elite",
                           parents=(i.sequence,))
                for i in ranked[:cfg.elite_count]
            ]
            children: list[Individual] = list(elites)
            seen = {i.sequence for i in children}

            guard = 0
            while len(children) < cfg.population_size and guard < cfg.population_size * 30:
                guard += 1
                if self.rng.random() < cfg.crossover_rate and len(pop) >= 2:
                    pa, pb = self._tournament(pop), self._tournament(pop)
                    seq = self._crossover(pa.sequence, pb.sequence)
                    seq = self._mutate(seq) if self.rng.random() < 0.5 else seq
                    origin, parents = "crossover", (pa.sequence, pb.sequence)
                else:
                    pa = self._tournament(pop)
                    seq = self._mutate(pa.sequence)
                    origin, parents = "mutation", (pa.sequence,)

                if seq in seen:
                    continue
                seen.add(seq)
                children.append(
                    Individual(seq, self._score(seq), g, origin, parents))

            pop = children
            gen = Generation(g, pop, time.time() - t0,
                             self._n_evaluated - before_e,
                             self._n_cached - before_c)
            self.history.append(gen)
            if verbose:
                print(gen.summary())
            if on_generation:
                on_generation(gen)

            # --- stopping conditions ---
            # Measured against the rank key, not the composite. The first
            # generation to find a candidate that is not disqualified is the
            # most important progress a run can make, and its composite is
            # usually LOWER than the counter-selective one it displaces --
            # on TNF-alpha, 72.6 against 82.6. Compared as bare floats that
            # reads as four generations of stagnation and stops the search
            # exactly when it has started working.
            key = gen.best.rank_key
            if (key[0] > best_key[0]
                    or (key[0] == best_key[0]
                        and key[1] > best_key[1] + cfg.min_improvement)):
                best_key = key
                stagnant = 0
            else:
                stagnant += 1

            # A disqualified candidate never satisfies a target score, however
            # high its composite.
            if (cfg.target_score is not None and not gen.best.card.disqualified
                    and gen.best.score >= cfg.target_score):
                self.stop_reason = (
                    f"target composite score {cfg.target_score:.3f} reached at "
                    f"generation {g}")
                break
            if stagnant >= cfg.patience:
                self.stop_reason = (
                    f"no improvement above {cfg.min_improvement} for "
                    f"{cfg.patience} generations -- the search plateaued, which "
                    f"says the search is done, not that the molecule works")
                break
        else:
            self.stop_reason = f"completed all {cfg.generations} generations"

        self.total_elapsed = time.time() - t_start
        if verbose:
            print(f"\nStopped: {self.stop_reason}")
            print(f"Total: {self.total_elapsed:.1f}s, "
                  f"{self._n_evaluated} unique sequences scored, "
                  f"{self._n_cached} cache hits")
        return self.best_generation()

    # --- results ----------------------------------------------------------

    def best_generation(self) -> Generation:
        return max(self.history, key=lambda g: g.best.rank_key)

    def best(self) -> Individual:
        return max((i for g in self.history for i in g.individuals),
                   key=lambda i: i.rank_key)

    def top(self, n: int = 10) -> list[Individual]:
        """Best n distinct sequences across the whole run."""
        seen: dict[str, Individual] = {}
        for g in self.history:
            for i in g.individuals:
                prev = seen.get(i.sequence)
                if prev is None or i.rank_key > prev.rank_key:
                    seen[i.sequence] = i
        return sorted(seen.values(), key=lambda i: i.rank_key,
                      reverse=True)[:n]

    def trajectory(self) -> list[dict]:
        return [
            {
                "generation": g.index,
                "best": round(g.best.score_100, 4),
                "mean": round(g.mean * 100, 4),
                "median": round(g.median * 100, 4),
                "stdev": round(g.stdev * 100, 4),
                "diversity": round(g.diversity, 4),
                "elapsed_s": round(g.elapsed_s, 3),
                "n_evaluated": g.n_evaluated,
                "n_cached": g.n_cached,
                "best_sequence": g.best.sequence,
                "best_designed": g.best.card.metadata.get("designed_residues", ""),
                # "Best" is the rank key's best, so this line can step DOWN
                # at the generation where the search first finds a candidate
                # that is not disqualified -- on TNF-alpha a 72.6 displacing
                # an 82.6. That step is the run working, not regressing, and
                # the chart needs the flag to be able to say so.
                "best_disqualified": bool(g.best.card.disqualified),
            }
            for g in self.history
        ]

    def save(self, path: str | Path, top_n: int = 25) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        ref = self.evaluator.reference
        payload = {
            "config": self.config.to_dict(),
            "scaffold": {
                "name": self.scaffold.name,
                "description": self.scaffold.description,
                "framework": self.scaffold.framework,
                "secondary_structure": self.scaffold.secondary_structure,
                "variable_positions_1based": [
                    p + 1 for p in self.scaffold.variable_positions],
                "provenance": self.scaffold.provenance,
                "search_space": self.scaffold.search_space,
            },
            "stop_reason": self.stop_reason,
            "total_elapsed_s": round(getattr(self, "total_elapsed", 0.0), 2),
            "n_unique_scored": self._n_evaluated,
            "n_cache_hits": self._n_cached,
            "trajectory": self.trajectory(),
            "reference": None if ref is None else {
                "name": ref.name,
                "sequence": ref.sequence,
                "description": ref.description,
                "score_100": round(ref.score_100, 3),
                "group_scores": {k: round(v, 4) for k, v in ref.group_scores.items()},
                "card": ref.card.to_dict(),
            },
            "top": [
                {
                    **i.card.to_dict(),
                    "generation": i.generation,
                    "origin": i.origin,
                    "group_scores": {
                        k: round(v, 4)
                        for k, v in self.evaluator.group_scores(i.card).items()},
                    "percent_of_reference": (
                        None if self.evaluator.percent_of_reference(i.card) is None
                        else round(self.evaluator.percent_of_reference(i.card), 2)),
                    "recognition_percent_of_reference": (
                        None
                        if self.evaluator.recognition_percent_of_reference(i.card) is None
                        else round(
                            self.evaluator.recognition_percent_of_reference(i.card), 2)),
                }
                for i in self.top(top_n)
            ],
        }
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        return path
