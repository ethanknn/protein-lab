"""Developability scorers: can this molecule actually be made and used?

These are the properties that kill real candidates long before efficacy is
ever tested. A binder that aggregates in the vial, gets cleared in an hour,
or cannot be expressed is worthless no matter how well it docks.

All thresholds are from standard protein-therapeutic practice and are stated
explicitly so they can be argued with.
"""

from __future__ import annotations

from ..aminoacids import (
    AGGREGATION_PRONE, AROMATIC, HYDROPHOBIC, LIABILITY_MOTIFS, SHEET_PROPENSITY,
    aliphatic_index, aromaticity, extinction_coefficient, gravy,
    hydrophobic_moment, isoelectric_point, max_homopolymer_run,
    molecular_weight, net_charge, shannon_entropy,
)
from .base import ScoreComponent, Scorer, penalty_score, ramp_score, window_score


class SolubilityScorer(Scorer):
    """Net charge, pI and hydrophobicity -- will it stay in solution?

    Rationale: proteins are least soluble near their isoelectric point, where
    net charge is zero and there is no electrostatic repulsion between
    molecules. A formulation sitting at pH 7.4 therefore wants a pI well away
    from 7.4. Excess surface hydrophobicity drives aggregation.
    """

    name = "solubility"
    weight = 1.2

    def score(self, sequence: str, context: dict) -> ScoreComponent:
        q = net_charge(sequence, 7.4)
        pi = isoelectric_point(sequence)
        g = gravy(sequence)

        # Charge magnitude: want enough to repel, not so much it binds
        # everything electrostatically. |q| of 2-8 is the useful band.
        charge_term = window_score(abs(q), 2.0, 8.0, soft=4.0)
        # pI at least 1.5 units from formulation pH.
        pi_term = ramp_score(abs(pi - 7.4), 0.5, 2.0)
        # GRAVY below ~0 means net hydrophilic, which is what we want exposed.
        gravy_term = window_score(g, -2.0, 0.0, soft=0.8)

        norm = 0.4 * charge_term + 0.25 * pi_term + 0.35 * gravy_term

        return ScoreComponent(
            name=self.name,
            raw=g,
            normalised=norm,
            weight=self.weight,
            units="GRAVY",
            detail=(
                f"net charge {q:+.2f} at pH 7.4, theoretical pI {pi:.2f}, "
                f"GRAVY {g:+.2f}. Wants |charge| 2-8, pI >1.5 units from 7.4, "
                f"GRAVY below 0."
            ),
            caveat=(
                "Sequence-level solubility proxies. Real solubility depends on "
                "the folded surface and the formulation buffer."
            ),
        )


class StabilityScorer(Scorer):
    """Will the fold hold together and resist aggregation?"""

    name = "stability"
    weight = 1.0

    def score(self, sequence: str, context: dict) -> ScoreComponent:
        ai = aliphatic_index(sequence)
        ent = shannon_entropy(sequence)
        run = max_homopolymer_run(sequence)
        apr = self._aggregation_windows(sequence)

        # Ikai aliphatic index: higher is more thermostable. 70-110 is typical
        # for well-behaved globular proteins.
        ai_term = window_score(ai, 60.0, 110.0, soft=40.0)
        # Entropy: a 20-letter alphabet caps at 4.32 bits. Below ~3.3 bits the
        # sequence is low-complexity, which correlates with disorder.
        ent_term = ramp_score(ent, 2.6, 3.6)
        # Long single-residue runs are a classic expression and stickiness problem.
        run_term = penalty_score(run, free=3, each=0.3)
        # Aggregation-prone windows.
        apr_term = penalty_score(len(apr), free=0, each=0.35)

        norm = 0.3 * ai_term + 0.25 * ent_term + 0.15 * run_term + 0.3 * apr_term

        apr_text = (", ".join(f"{s}-{e}:{sequence[s:e]}" for s, e in apr[:3])
                    or "none")
        return ScoreComponent(
            name=self.name,
            raw=ai,
            normalised=norm,
            weight=self.weight,
            units="aliphatic index",
            detail=(
                f"aliphatic index {ai:.1f}, sequence entropy {ent:.2f} bits, "
                f"longest repeat run {run}, aggregation-prone windows: {apr_text}."
            ),
            caveat=(
                "Aggregation windows here are a transparent sliding-window "
                "heuristic, not a validated predictor like TANGO or AGGRESCAN."
            ),
        )

    @staticmethod
    def _aggregation_windows(seq: str, window: int = 6,
                             hydro_cut: float = 0.70,
                             beta_cut: float = 1.15) -> list[tuple[int, int]]:
        """Stretches that are both hydrophobic and beta-prone.

        That combination is the signature of an amyloid-like aggregation hot
        spot in short peptides. Labelled a heuristic on purpose.
        """
        hits: list[tuple[int, int]] = []
        if len(seq) < window:
            return hits
        for i in range(len(seq) - window + 1):
            w = seq[i:i + window]
            frac_agg = sum(1 for a in w if a in AGGREGATION_PRONE) / window
            beta = sum(SHEET_PROPENSITY[a] for a in w) / window
            if frac_agg >= hydro_cut and beta >= beta_cut:
                if hits and i <= hits[-1][1]:
                    hits[-1] = (hits[-1][0], i + window)  # merge overlap
                else:
                    hits.append((i, i + window))
        return hits


class ManufacturabilityScorer(Scorer):
    """Chemical liabilities that cause trouble during production and storage."""

    name = "manufacturability"
    weight = 0.9

    def score(self, sequence: str, context: dict) -> ScoreComponent:
        found: list[str] = []

        # Dipeptide liabilities: deamidation and isomerisation hotspots.
        for motif, desc in LIABILITY_MOTIFS.items():
            if len(motif) != 2:
                continue
            idx = sequence.find(motif)
            while idx != -1:
                found.append(f"{motif}@{idx + 1} ({desc})")
                idx = sequence.find(motif, idx + 1)

        # N-linked glycosylation sequons N-X-S/T where X is not proline.
        sequons = [
            i + 1 for i in range(len(sequence) - 2)
            if sequence[i] == "N" and sequence[i + 1] != "P"
            and sequence[i + 2] in "ST"
        ]

        # Unpaired cysteine: free thiols scramble disulfides and dimerise.
        n_cys = sequence.count("C")
        cys_unpaired = n_cys % 2 == 1

        # Oxidation-prone residues in excess.
        n_met = sequence.count("M")

        eps = extinction_coefficient(sequence)

        motif_term = penalty_score(len(found), free=1, each=0.18)
        sequon_term = penalty_score(len(sequons), free=0, each=0.3)
        cys_term = 0.4 if cys_unpaired else 1.0
        met_term = penalty_score(n_met, free=2, each=0.2)
        # Needs at least one Trp or Tyr to be quantifiable by A280 in a lab.
        quant_term = 1.0 if eps >= 1490 else 0.5

        norm = (0.3 * motif_term + 0.2 * sequon_term + 0.2 * cys_term
                + 0.15 * met_term + 0.15 * quant_term)

        bits = []
        if found:
            bits.append(f"{len(found)} liability motifs ({'; '.join(found[:3])})")
        if sequons:
            bits.append(f"N-glyc sequons at {sequons}")
        if cys_unpaired:
            bits.append(f"odd cysteine count ({n_cys}) leaves a free thiol")
        bits.append(f"A280 extinction {eps:.0f} M-1cm-1")

        comp = ScoreComponent(
            name=self.name,
            raw=float(len(found) + len(sequons)),
            normalised=norm,
            weight=self.weight,
            units="liabilities",
            detail="; ".join(bits) if bits else "no liabilities detected.",
            caveat=(
                "Motif-based flags. Whether a hotspot actually degrades depends "
                "on solvent exposure in the folded structure."
            ),
        )
        return comp


class ScaffoldFitScorer(Scorer):
    """Is the sequence consistent with the fold we claim it has?

    A designed helical binder whose sequence has beta-sheet propensity
    everywhere is not going to adopt the structure the design assumes. This
    checks the sequence against the scaffold's intended secondary structure.
    """

    name = "scaffold_consistency"
    weight = 1.0

    def score(self, sequence: str, context: dict) -> ScoreComponent:
        from ..aminoacids import HELIX_PROPENSITY

        scaffold = context.get("scaffold")
        ss = getattr(scaffold, "secondary_structure", None) if scaffold else None

        if not ss or len(ss) != len(sequence):
            # No structural claim to check: score the generic foldability proxy.
            mom = hydrophobic_moment(sequence)
            return ScoreComponent(
                name=self.name,
                raw=mom,
                normalised=ramp_score(mom, 0.05, 0.45),
                weight=self.weight,
                units="hydrophobic moment",
                detail=(
                    f"no scaffold secondary structure supplied; scored on "
                    f"amphipathicity instead (moment {mom:.3f}). Higher means "
                    f"hydrophobic and polar faces are better segregated."
                ),
                caveat="Not a structure prediction.",
            )

        agree = 0.0
        checked = 0
        for aa, code in zip(sequence, ss):
            if code == "H":
                agree += min(HELIX_PROPENSITY[aa] / 1.3, 1.0)
                checked += 1
            elif code == "E":
                agree += min(SHEET_PROPENSITY[aa] / 1.5, 1.0)
                checked += 1
            elif code in "LT":  # loop / turn: glycine and proline are assets
                agree += 1.0 if aa in "GPSDNT" else 0.55
                checked += 1

        norm = agree / checked if checked else 0.5

        # Helical faces should be amphipathic, so reward the moment too.
        if "H" in ss:
            mom = hydrophobic_moment(sequence)
            norm = 0.75 * norm + 0.25 * ramp_score(mom, 0.05, 0.40)

        return ScoreComponent(
            name=self.name,
            raw=norm,
            normalised=norm,
            weight=self.weight,
            units="agreement",
            detail=(
                f"sequence-vs-scaffold agreement over {checked} positions "
                f"using Chou-Fasman propensities; scaffold topology "
                f"{getattr(scaffold, 'name', '?')} ({ss})."
            ),
            caveat=(
                "Chou-Fasman propensities are a 1978 statistical model. They "
                "indicate consistency, not that the fold will form."
            ),
        )


class SizeScorer(Scorer):
    """Is the candidate the right size for the job?

    Too short and there is not enough surface to bind specifically; too long
    and it will not fold predictably, costs more to make, and penetrates
    tissue poorly.
    """

    name = "size_suitability"
    weight = 0.6

    def score(self, sequence: str, context: dict) -> ScoreComponent:
        n = len(sequence)
        mw = molecular_weight(sequence) / 1000.0

        # The acceptable band belongs to the scaffold: penalising an antibody
        # for being antibody-sized would make the reference comparison
        # meaningless rather than honest.
        scaffold = context.get("scaffold")
        lo, hi = getattr(scaffold, "size_band", (25.0, 65.0)) if scaffold \
            else (25.0, 65.0)
        norm = window_score(float(n), lo, hi, soft=max(0.3 * (hi - lo), 10.0))

        return ScoreComponent(
            name=self.name,
            raw=float(n),
            normalised=norm,
            weight=self.weight,
            units="residues",
            detail=(
                f"{n} residues, {mw:.2f} kDa against a target band of "
                f"{lo:.0f}-{hi:.0f} for the "
                f"{getattr(scaffold, 'name', 'default')} format. Small formats "
                f"fold autonomously and penetrate marrow better; large ones "
                f"bury more surface and last longer in circulation."
            ),
            caveat=(
                "The band is a format convention, not a measurement. Size "
                "trades off against half-life and tissue penetration."
            ),
        )
