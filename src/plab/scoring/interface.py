"""Interface complementarity: would this face actually grip the epitope?

The CLEC12A epitope read off PDB 8W9J is dominated by three arginines
(R201, R204, R232) sitting on the top face of the lectin domain. The real
antibody 50C1 answers that with acidic residues forming salt bridges and
tyrosines stacking against the guanidinium groups. Those two interactions --
salt bridge and cation-pi -- are the physical basis of this recognition
event, so they are what we score.

Each scorer below measures one physical mechanism, separately, so a candidate
can be diagnosed rather than just ranked.
"""

from __future__ import annotations

from ..aminoacids import (
    AROMATIC, HBOND_ACCEPTORS, HBOND_DONORS, HYDROPHOBIC, NEGATIVE_SIDECHAINS,
    POLAR_UNCHARGED, net_charge,
)
from .base import ScoreComponent, Scorer, ramp_score, target_score, window_score

# Physicochemical classes used for profile comparison. Deliberately coarse:
# we want "does this face have the right character", not sequence identity.
CLASSES = {
    "acidic": set("DE"),
    "basic": set("KRH"),
    "aromatic": set("FWY"),
    "aliphatic": set("AVILM"),
    "polar": set("STNQ"),
    "special": set("GPC"),
}


def class_profile(seq: str) -> dict[str, float]:
    if not seq:
        return {k: 0.0 for k in CLASSES}
    n = len(seq)
    return {k: sum(1 for a in seq if a in members) / n
            for k, members in CLASSES.items()}


def _interface_subsequence(sequence: str, context: dict) -> str:
    """The residues the scaffold presents to the target, if it declares any.

    Falls back to the whole sequence, which is the right behaviour for a short
    unstructured peptide where every residue is potentially in contact.
    """
    scaffold = context.get("scaffold")
    positions = getattr(scaffold, "interface_positions", None) if scaffold else None
    if not positions:
        return sequence
    return "".join(sequence[i] for i in positions if 0 <= i < len(sequence))


class SaltBridgeCapacityScorer(Scorer):
    """Can the binder pair the epitope's arginines?

    Three arginines need three acidic partners to be neutralised. Fewer leaves
    affinity on the table; many more adds nonspecific electrostatics and hurts
    selectivity, so this is a window, not a ramp.
    """

    name = "salt_bridge_capacity"
    # Reduced from 1.4. Measured across 500 sequences this correlates with
    # paralog_selectivity at r = +0.870: both are driven by how well the
    # face's acidic supply matches the epitope's cation demand, so at 1.4
    # each the composite was weighting that single quantity 2.8. Selectivity
    # keeps the higher weight because it measures the same matching AND the
    # margin over the paralog panel, which is strictly more information.
    weight = 0.9

    def score(self, sequence: str, context: dict) -> ScoreComponent:
        target = context["target"]
        face = _interface_subsequence(sequence, context)

        # Count all cations, not arginines only: a lysine-rich surface needs
        # acidic partners just as much, and several catalogue targets have no
        # arginine in the epitope at all.
        n_arg = len(target.cation_positions)
        n_acidic = sum(1 for a in face if a in "DE")
        # Epitope acidic residues want a basic partner on our side too.
        n_epi_acidic = sum(1 for r in target.epitope if r.aa in "DE")
        n_basic = sum(1 for a in face if a in "KR")

        if n_arg == 0 and n_epi_acidic == 0:
            # An uncharged surface: salt bridges are not the mechanism here,
            # so report that rather than scoring a criterion that does not
            # apply. Excess charge is still penalised, since it buys
            # nonspecific binding and nothing else.
            excess = max(0, n_acidic + n_basic - 3)
            return ScoreComponent(
                name=self.name, raw=float(n_acidic + n_basic),
                normalised=max(0.0, 1.0 - 0.15 * excess),
                weight=self.weight * 0.4, units="charged residues",
                detail=(f"the epitope carries no charged residues, so salt "
                        f"bridges are not the recognition mechanism here; "
                        f"scored only for excess charge ({n_acidic} acidic, "
                        f"{n_basic} basic on the face)."),
                caveat="Not applicable to this surface; down-weighted.",
            )

        # Want roughly one acidic residue per epitope arginine, up to ~2x.
        # The band is deliberately tight: a soft tail of 3 residues let the
        # wild-type scaffold (2 acidics, no affinity for this target) score
        # 77% here, which is not a distinction worth making.
        acidic_term = window_score(float(n_acidic), float(n_arg),
                                   float(n_arg * 2), soft=1.5)
        # A small number of basic residues to answer the epitope aspartates.
        basic_term = window_score(float(n_basic), max(1.0, n_epi_acidic - 1),
                                  float(n_epi_acidic + 1), soft=1.5)
        # Acid/base balance: an all-acidic face is a charge magnet, not a
        # specific binder. Reward having both, weighted toward acidic.
        if n_acidic + n_basic > 0:
            ratio = n_acidic / (n_acidic + n_basic)
            balance_term = target_score(ratio, 0.70, tolerance=0.45)
        else:
            balance_term = 0.0

        norm = 0.55 * acidic_term + 0.2 * basic_term + 0.25 * balance_term

        return ScoreComponent(
            name=self.name,
            raw=float(n_acidic),
            normalised=norm,
            weight=self.weight,
            units="acidic residues",
            detail=(
                f"binding face has {n_acidic} acidic and {n_basic} basic residues "
                f"for an epitope with {n_arg} arginines "
                f"({', '.join(target.arginine_positions)}) and "
                f"{n_epi_acidic} aspartates. Wants {n_arg}-{n_arg * 2} acidic."
            ),
            caveat=(
                "Counts capacity, not geometry. Whether a given Asp reaches a "
                "given Arg depends on the docked pose, which this does not model."
            ),
        )


class CationPiScorer(Scorer):
    """Aromatic residues for stacking against the arginine cluster.

    Cation-pi interactions between arginine guanidinium groups and aromatic
    rings are a well-characterised and often dominant contributor at
    protein interfaces. The reference antibody uses five tyrosines here,
    which is a strong hint that this is the mechanism that matters.
    """

    name = "cation_pi_potential"
    weight = 1.3

    def score(self, sequence: str, context: dict) -> ScoreComponent:
        target = context["target"]
        face = _interface_subsequence(sequence, context)
        if not face:
            return ScoreComponent(self.name, 0.0, 0.0, self.weight,
                                  "empty binding face")

        n_arom = sum(1 for a in face if a in AROMATIC)
        frac = n_arom / len(face)
        n_tyr = face.count("Y")
        ref_frac = target.paratope_aromatic_fraction
        n_cation = len(target.cation_positions)

        # Aim for the reference's aromatic density, which demonstrably works.
        frac_term = target_score(frac, ref_frac, tolerance=max(ref_frac, 0.25))
        # Tyrosine specifically: ring plus a hydroxyl that can hydrogen bond,
        # which is why it dominates real paratopes.
        tyr_term = ramp_score(float(n_tyr), 0.0, max(2.0, len(face) * 0.15))
        # One aromatic per epitope cation is the stacking budget. Where the
        # surface has no cations there is no cation-pi to score, so the
        # budget term falls back to the reference's aromatic density.
        budget_term = (ramp_score(float(n_arom), 0.0, float(n_cation))
                       if n_cation else frac_term)

        norm = 0.4 * frac_term + 0.3 * tyr_term + 0.3 * budget_term

        return ScoreComponent(
            name=self.name,
            raw=frac,
            normalised=norm,
            weight=self.weight,
            units="aromatic fraction",
            detail=(
                f"{n_arom} aromatics ({n_tyr} Tyr) across a {len(face)}-residue "
                f"face = {frac:.0%}; reference paratope is {ref_frac:.0%}. "
                f"Budget is one aromatic per epitope arginine "
                f"({len(target.arginine_positions)})."
            ),
            caveat=(
                "Stacking requires the ring and the guanidinium to be roughly "
                "parallel and within ~4-5 A. Composition alone cannot confirm that."
            ),
        )


class ChargeComplementarityScorer(Scorer):
    """Does the overall electrostatic field oppose the epitope's?

    Long-range Coulombic steering genuinely matters for association rate:
    opposite net charges increase the on-rate before any specific contact
    forms. This is the crudest but most reliable interface signal.
    """

    name = "charge_complementarity"
    weight = 1.1

    def score(self, sequence: str, context: dict) -> ScoreComponent:
        target = context["target"]
        face = _interface_subsequence(sequence, context)
        if not face:
            return ScoreComponent(self.name, 0.0, 0.0, self.weight, "empty face")

        # Charge DENSITY on the binding face, not net charge of the whole
        # molecule. An earlier version scored whole-molecule net charge and
        # gave full marks to anything mildly acidic, including the wild-type
        # scaffold, which binds a completely different protein. Density on the
        # contact face is the quantity that actually steers association.
        q_epi = target.epitope_net_charge
        epi_density = q_epi / max(len(target.epitope), 1)

        face_charge = sum(
            1.0 if a in "KR" else (0.1 if a == "H" else (-1.0 if a in "DE" else 0.0))
            for a in face
        )
        face_density = face_charge / len(face)

        ref_face = target.paratope_sequence
        ref_charge = sum(
            1.0 if a in "KR" else (0.1 if a == "H" else (-1.0 if a in "DE" else 0.0))
            for a in ref_face
        ) if ref_face else 0.0
        ref_density = ref_charge / max(len(ref_face), 1)

        if q_epi * face_charge >= 0:
            # Same sign or neutral: no electrostatic steering to speak of.
            norm = 0.15 if face_charge == 0 else 0.0
        else:
            # Opposite sign: how close is the magnitude to what the real
            # binder uses? Both too little and too much lose points -- excess
            # charge buys nonspecific electrostatics, not selectivity.
            norm = target_score(abs(face_density), abs(ref_density),
                                tolerance=max(abs(ref_density), 0.12))

        verdict = ("opposing (favourable)" if q_epi * face_charge < 0
                   else "same sign (repulsive)" if q_epi * face_charge > 0
                   else "no net charge on the face")

        return ScoreComponent(
            name=self.name,
            raw=face_density,
            normalised=norm,
            weight=self.weight,
            units="charge per face residue",
            detail=(
                f"face carries {face_charge:+.1f} over {len(face)} residues "
                f"({face_density:+.3f}/residue) against an epitope at "
                f"{epi_density:+.3f}/residue -- {verdict}. Reference paratope "
                f"runs {ref_density:+.3f}/residue."
            ),
            caveat=(
                "Charge density on the face is still an average; real steering "
                "depends on the dipole and where the charges actually sit."
            ),
        )


class HydrogenBondScorer(Scorer):
    """Donor/acceptor availability against the epitope's polar residues.

    The CLEC12A epitope carries Ser, Thr, Gln and Asp, all of which want
    hydrogen-bond partners. Interfaces that bury polar groups without
    satisfying them pay a desolvation penalty.
    """

    name = "hbond_complementarity"
    weight = 0.8

    def score(self, sequence: str, context: dict) -> ScoreComponent:
        target = context["target"]
        face = _interface_subsequence(sequence, context)
        if not face:
            return ScoreComponent(self.name, 0.0, 0.0, self.weight, "empty face")

        epi = target.epitope_sequence
        epi_donors = sum(1 for a in epi if a in HBOND_DONORS)
        epi_acceptors = sum(1 for a in epi if a in HBOND_ACCEPTORS)

        our_donors = sum(1 for a in face if a in HBOND_DONORS)
        our_acceptors = sum(1 for a in face if a in HBOND_ACCEPTORS)

        # Densities, matched to what the reference paratope achieves. Raw
        # counts were too easy to satisfy -- any polar protein maxed them out
        # regardless of whether the groups complemented anything.
        ref = target.paratope_sequence
        n_face = len(face)
        ref_n = max(len(ref), 1)
        ref_donor_d = sum(1 for a in ref if a in HBOND_DONORS) / ref_n
        ref_accept_d = sum(1 for a in ref if a in HBOND_ACCEPTORS) / ref_n

        donor_d = our_donors / n_face
        accept_d = our_acceptors / n_face

        d_term = target_score(accept_d, ref_accept_d,
                              tolerance=max(ref_accept_d, 0.2))
        a_term = target_score(donor_d, ref_donor_d,
                              tolerance=max(ref_donor_d, 0.2))
        # Polar groups buried without a partner cost desolvation energy, so a
        # face that is overwhelmingly polar is penalised, not rewarded.
        excess = max(0.0, (donor_d + accept_d) - 1.4)
        norm = max(0.0, 0.5 * d_term + 0.5 * a_term - excess)

        return ScoreComponent(
            name=self.name,
            raw=donor_d + accept_d,
            normalised=norm,
            weight=self.weight,
            units="polar groups per face residue",
            detail=(
                f"face offers {our_donors} donors / {our_acceptors} acceptors "
                f"over {n_face} residues ({donor_d:.2f}/{accept_d:.2f} per "
                f"residue) against an epitope with {epi_donors} donors / "
                f"{epi_acceptors} acceptors. Reference runs "
                f"{ref_donor_d:.2f}/{ref_accept_d:.2f}."
            ),
            caveat="Counts availability, not satisfied bonds in a docked pose.",
        )


class ParatopeMimicryScorer(Scorer):
    """How closely does the face resemble the character of a binder that works?

    This compares distributions over physicochemical classes, not sequence
    identity -- the point is to capture the strategy 50C1 uses (acidic,
    tyrosine-heavy, few aliphatics) without rewarding plagiarism of its
    sequence, which would be both trivial and patent-encumbered.
    """

    name = "paratope_character"
    weight = 1.0

    def score(self, sequence: str, context: dict) -> ScoreComponent:
        target = context["target"]
        face = _interface_subsequence(sequence, context)
        ref = target.paratope_sequence
        if not face or not ref:
            return ScoreComponent(self.name, 0.0, 0.0, self.weight,
                                  "no reference paratope available")

        p = class_profile(face)
        q = class_profile(ref)
        # Total variation distance between the two class distributions.
        tvd = 0.5 * sum(abs(p[k] - q[k]) for k in CLASSES)
        norm = 1.0 - tvd

        # Report the biggest divergence, which tells the designer what to fix.
        worst = max(CLASSES, key=lambda k: abs(p[k] - q[k]))
        direction = "too much" if p[worst] > q[worst] else "too little"

        return ScoreComponent(
            name=self.name,
            raw=tvd,
            normalised=norm,
            weight=self.weight,
            units="class distance",
            detail=(
                f"class profile distance {tvd:.3f} from the reference paratope. "
                f"Biggest gap: {direction} {worst} "
                f"({p[worst]:.0%} vs reference {q[worst]:.0%})."
            ),
            caveat=(
                "Resembling a working binder's chemistry is suggestive, not "
                "sufficient -- the geometry has to be right too."
            ),
        )


class FootprintScorer(Scorer):
    """Is the binding face big enough to cover the epitope, and no bigger?

    The epitope has a radius of gyration of ~8.7 A and 13 residues. A face
    much smaller cannot bury enough surface for useful affinity; much larger
    and it necessarily contacts neighbouring surface, which costs selectivity.
    """

    name = "footprint_match"
    # Low weight on purpose: for a fixed scaffold this is near-constant across
    # every candidate, so it carries almost no information within a run. It is
    # kept because it does discriminate between scaffolds.
    weight = 0.4

    def score(self, sequence: str, context: dict) -> ScoreComponent:
        target = context["target"]
        face = _interface_subsequence(sequence, context)
        n_epi = len(target.epitope)

        # Real interfaces pair roughly one binder residue per epitope residue.
        norm = window_score(float(len(face)), n_epi * 0.7, n_epi * 1.6,
                            soft=float(n_epi))

        return ScoreComponent(
            name=self.name,
            raw=float(len(face)),
            normalised=norm,
            weight=self.weight,
            units="face residues",
            detail=(
                f"{len(face)}-residue binding face against a {n_epi}-residue "
                f"epitope (Rg {target.epitope_radius:.1f} A). Wants "
                f"{int(n_epi * 0.7)}-{int(n_epi * 1.6)} residues."
            ),
            caveat=(
                "Residue counting is a stand-in for buried surface area, which "
                "needs a docked structure to compute properly."
            ),
        )
