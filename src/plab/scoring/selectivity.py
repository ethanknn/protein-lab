"""Selectivity: would this bind things we do not want it to bind?

This is the part of the problem that has actually killed leukaemia
candidates. CD33- and CD123-directed agents hit healthy haematopoietic
stem and progenitor cells and cause lasting marrow toxicity. CLEC12A was
chosen precisely because its expression spares stem cells.

What this module can and cannot do
----------------------------------
CAN: compare how well a candidate's chemistry complements the target
epitope versus the structurally equivalent surface on its closest human
relatives, and flag the charge/hydrophobicity signature associated with
promiscuous binding.

CANNOT: predict tissue distribution or stem-cell sparing. Whether a target
is present on a healthy cell is a fact about that cell's biology, measured
by flow cytometry on primary human samples. No sequence score reaches it.
That is stated here because conflating the two would be the single most
misleading thing this software could do.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..align import smith_waterman
from ..aminoacids import AA20, AROMATIC, HYDROPHOBIC, net_charge
from .base import ScoreComponent, Scorer, ramp_score, window_score


def _supply_match(supply: float, demand: float,
                  overshoot_penalty: float = 0.5) -> float:
    """How well does `supply` meet `demand`, penalising excess as well as shortfall?

    This is the heart of the selectivity fix. The original kernel used
    min(supply, demand) / demand, which saturates at 1.0 -- so a face with
    enough acidic residues for CLEC12A's three arginines ALSO scored a
    perfect 1.0 against CLEC9A's single arginine, and the difference between
    target and paralog collapsed. A sweep of 680 faces across the whole
    composition range found a maximum achievable margin of 0.075, meaning no
    design could ever score above 50%. The optimiser was not failing at
    selectivity; it was pinned at the ceiling of a metric that could not
    report higher.

    Penalising overshoot is what creates discrimination, and it is also the
    correct chemistry: an acidic residue with no cation to pair against on
    the target is not neutral, it is precisely the residue that goes and
    binds something else.
    """
    if demand <= 0:
        # Nothing to satisfy here, so any supply is unpaired and unhelpful.
        return max(0.0, 1.0 - supply * overshoot_penalty)
    ratio = supply / demand
    if ratio <= 1.0:
        return ratio
    return max(0.0, 1.0 - (ratio - 1.0) * overshoot_penalty)


def complementarity_kernel(face: str, epitope: str) -> float:
    """A single scalar: how well does `face` chemically complement `epitope`?

    Intentionally the same function for the real target and for a paralog, so
    the difference between the two is a fair comparison rather than two
    different yardsticks.

    The epitope's demand is scaled to the size of the face before it is used.
    Without that scaling the kernel reads surface *length* as if it were
    chemistry: writing the same epitope out twice -- identical residues,
    identical proportions, not one atom changed -- moved the result from
    0.475 to 0.725, and doubling the face moved it to 0.275. That is not a
    subtle bias. Epitopes in this catalogue run from 11 to 34 residues, and
    28 of 70 panel surrogates are a different length from the target surface
    they are compared against, so two fifths of every cross-reactivity
    comparison was partly a measurement of how many residues a curator had
    annotated. A face of 13 residues can only pair with about 13 epitope
    residues whatever the annotation says; scaling the demand to the face is
    what makes the ratio a statement about composition alone. Every other
    scorer in `interface.py` already works in densities for this reason.
    """
    if not face or not epitope:
        return 0.0

    # Pairing partners the face can actually reach, not everything the
    # epitope annotation happens to contain.
    reach = len(face) / len(epitope)

    n_arg = sum(1 for a in epitope if a == "R")
    n_lys = sum(1 for a in epitope if a == "K")
    n_cation = (n_arg + n_lys) * reach
    n_epi_acid = sum(1 for a in epitope if a in "DE") * reach
    n_epi_hydro = sum(1 for a in epitope if a in HYDROPHOBIC) * reach

    n_acid = sum(1 for a in face if a in "DE")
    n_base = sum(1 for a in face if a in "KR")
    n_arom = sum(1 for a in face if a in AROMATIC)
    n_hydro = sum(1 for a in face if a in HYDROPHOBIC)

    # Pair cations with acids, epitope acids with our bases, hydrophobic with
    # hydrophobic. Each term now peaks when supply meets demand and falls off
    # when it overshoots, which is what lets the target and a paralog of
    # different composition pull apart.
    salt = _supply_match(n_acid, n_cation)
    rev_salt = _supply_match(n_base, n_epi_acid)
    stack = _supply_match(n_arom, n_cation)
    # Hydrophobic burial is more forgiving of excess than charge is: a spare
    # leucine is far less promiscuous than a spare aspartate.
    apolar = _supply_match(n_hydro, n_epi_hydro, overshoot_penalty=0.25)

    # Electrostatic steering: reward opposing net charge on the side chains.
    q_epi = (n_cation - n_epi_acid)
    q_face = (n_base - n_acid)
    steer = 1.0 if q_epi * q_face < 0 else (0.5 if q_epi * q_face == 0 else 0.0)

    return 0.30 * salt + 0.15 * rev_salt + 0.30 * stack + 0.10 * apolar + 0.15 * steer


def margin_to_score(margin: float, lo: float, hi: float) -> float:
    """Turn a kernel margin into a 0..1 score, split at zero.

    Zero is a real boundary on this axis, not a point on a ramp: above it the
    face complements the target more than its nearest relative, below it the
    reverse. Scoring the whole range on one ramp anchored to random faces
    threw that sign away -- see COUNTER_SELECTIVE_BAND.

    Lives here, and is used by both the scorer and scripts/diagnose_panels.py,
    because a control that computes the percentage its own way is not checking
    the number anyone actually reads.
    """
    if margin < 0:
        return COUNTER_SELECTIVE_BAND * ramp_score(margin, min(lo, 0.0), 0.0)
    return COUNTER_SELECTIVE_BAND + (1.0 - COUNTER_SELECTIVE_BAND) * ramp_score(
        margin, 0.0, hi)


def achievable_margin_range(epitope: str, paralog_epitopes: dict[str, str],
                            face_length: int, *, n_samples: int = 800,
                            seed: int = 77) -> tuple[float, float]:
    """Measure what margins are actually reachable, to set the score's scale.

    Hard-coding the ramp bounds was how the old version ended up reporting
    50% for every design: the bounds assumed a range the metric could not
    reach. Sampling the real achievable range instead means the score always
    discriminates, and it stays correct when paralogs are added or removed.

    Two different samples set the two ends, deliberately:

    lo -- the 10th percentile of RANDOM faces, i.e. what a poor but realistic
          design looks like.
    hi -- the best margin found by a SYSTEMATIC sweep of compositions, i.e.
          the most selective this scaffold could possibly be.

    Anchoring the top to random faces instead was not enough: designs
    optimised for the other objectives already beat random on selectivity, so
    they all pinned at 100% and the gradient vanished again at the other end.
    Anchoring to the reachable optimum means 100% genuinely denotes "as
    selective as this epitope allows", and leaves the optimiser room to climb.

    Returns (lo, hi).
    """
    import random as _random

    if not epitope or not paralog_epitopes or face_length <= 0:
        return (0.0, 0.5)

    def margin_of(face: str) -> float:
        on = complementarity_kernel(face, epitope)
        off = max(complementarity_kernel(face, p)
                  for p in paralog_epitopes.values())
        return on - off

    # Lower anchor: random faces.
    rng = _random.Random(seed)
    alphabet = [a for a in AA20 if a != "C"]
    randoms = sorted(
        margin_of("".join(rng.choice(alphabet) for _ in range(face_length)))
        for _ in range(n_samples)
    )
    lo = randoms[int(0.10 * (len(randoms) - 1))]

    # Upper anchor: systematic sweep over the composition classes that drive
    # the kernel. Cheap -- this is arithmetic on counts, no model involved.
    #
    # The hydrophobic class is swept explicitly. Leaving it out, as an
    # earlier version did, made the ceiling too low: the kernel's apolar term
    # was unreachable, real antibody paratopes exceeded the modelled maximum,
    # and their scores clipped at exactly 100 -- the same silent saturation
    # this scorer has already been fixed for twice.
    #
    # The sweep enumerates compositions, so its cost is O(face_length^4). For
    # the 13-residue design face that is 2,380 faces and nothing to worry
    # about. The REFERENCE molecule's face is thirty or forty residues, and
    # 41^4/24 is a hundred and thirty-five thousand -- which is why scoring
    # one reference took twenty seconds and the test suite took twenty
    # minutes, measured rather than guessed.
    #
    # It does not need to scale with length. Since the kernel was fixed to
    # read densities rather than raw counts, its value depends on the face
    # only through the four class FRACTIONS, so a coarser grid on a long
    # face covers the same chemistry. The step is 1 for any face of 13 or
    # fewer residues, which is every scaffold in this project, so no stored
    # ceiling and no design score changes by a digit.
    step = max(1, face_length // 13)
    best = max(randoms)
    for n_acid in range(0, face_length + 1, step):
        for n_arom in range(0, face_length - n_acid + 1, step):
            for n_base in range(0, face_length - n_acid - n_arom + 1, step):
                left = face_length - n_acid - n_arom - n_base
                for n_apolar in range(0, left + 1, step):
                    face = ("D" * n_acid + "Y" * n_arom + "K" * n_base
                            + "L" * n_apolar + "S" * (left - n_apolar))
                    m = margin_of(face)
                    if m > best:
                        best = m

    hi = best
    if hi <= lo:
        hi = lo + 0.1
    return (lo, hi)


def ceiling_noise(epitope: str, paralog_epitopes: dict[str, str],
                  face_length: int) -> float:
    """How far the reachable ceiling moves if one mapped residue is wrong.

    The surrogate surfaces are read off local alignments at 25-40% identity,
    so at least one residue being misassigned is a certainty, not a worry.
    Dropping each mapped residue in turn and recomputing the ceiling gives
    the spread the ceiling has for reasons that are not chemistry.

    This exists because `INDISTINGUISHABLE_MARGIN` -- an absolute 0.02 --
    turned out to be the wrong KIND of threshold. VEGF-A's ceiling is 0.029,
    so it passed, and its optimised design was then reported as 100%
    selective. But VEGF-A's wobble is 0.027: the entire reachable margin is
    the size of the error bar on itself. A fixed cutoff cannot know that,
    because the wobble ranges from 0.003 to 0.077 across this catalogue -- a
    factor of twenty-five.

    Returned as a standard deviation, so the usual two-sigma convention
    applies: a ceiling is worth reading only if it is at least twice this.
    """
    import statistics

    his: list[float] = []
    for gene, surface in paralog_epitopes.items():
        for i in range(len(surface)):
            alt = dict(paralog_epitopes)
            alt[gene] = surface[:i] + surface[i + 1:]
            his.append(achievable_margin_range(
                epitope, alt, face_length, n_samples=16)[1])
    return statistics.pstdev(his) if len(his) > 1 else 0.0


class ParalogSelectivityScorer(Scorer):
    """Margin between complementarity to the target and to its closest relative.

    The epitope positions are mapped onto each panel member by local
    alignment, giving a surrogate surface for each. A candidate that scores
    equally well against the target and against a relative is not
    discriminating between them, whatever its affinity.

    The panel is built per target by scripts/build_panels.py from the human
    proteins sharing a domain family, ranked by similarity at the epitope
    window. For CLEC12A that independently returns CLEC9A as the nearest
    relative, which is what the hand-written panel said.
    """

    name = "paralog_selectivity"
    # Raised from 1.2: for this target selectivity is the clinically decisive
    # axis, and until the kernel was fixed this component could not vary, so
    # weighting it more heavily would have achieved nothing. Now that it
    # discriminates, the weight reflects what actually matters.
    weight = 1.4

    def score(self, sequence: str, context: dict) -> ScoreComponent:
        target = context["target"]
        face = context.get("face") or sequence
        surrogates: dict[str, str] = context.get("paralog_epitopes") or {}

        on_target = complementarity_kernel(face, target.epitope_sequence)

        if not surrogates:
            # Say WHY there is no panel. "This protein has no human relatives"
            # and "no panel was built" are opposite claims about how much is
            # known, and a single generic message makes them look identical.
            state = getattr(target, "panel_state", "none")
            why = {
                "no_human_relatives": (
                    "no other reviewed human protein shares a domain family "
                    "with this target, so cross-reactivity with a close human "
                    "relative is not the risk here. Binding to unrelated "
                    "human proteins still is, and is not assessed"),
                "no_mappable_relatives": (
                    "this target has human relatives, but none could have the "
                    "epitope mapped onto it well enough to compare -- they "
                    "diverge at this particular surface"),
                "unmappable": (
                    "the structure's residue numbering could not be reconciled "
                    "with the reference sequence, so the epitope cannot be "
                    "located in any relative"),
                "no_family": (
                    "UniProt lists no domain family for this protein, so there "
                    "is no principled set of relatives to compare against"),
            }.get(state, "no panel of related proteins was built for this target")
            return ScoreComponent(
                name=self.name,
                raw=on_target,
                normalised=0.5,
                weight=self.weight,
                units="margin",
                detail=(f"Selectivity was not measured: {why}. Scored neutral "
                        f"and excluded from the ranking rather than assumed "
                        f"good."),
                caveat="Neutral placeholder, not a measurement.",
            )

        ranked = sorted(
            ((name, complementarity_kernel(face, epi))
             for name, epi in surrogates.items()),
            key=lambda kv: -kv[1],
        )
        worst_name, worst_val = ranked[0]
        margin = on_target - worst_val

        # Scale taken from the measured achievable range rather than assumed.
        lo, hi = context.get("selectivity_scale", (-0.10, 0.25))

        # The ceiling has its own error bar, measured by dropping single
        # mapped residues (scripts/measure_ceilings.py). A ceiling no bigger
        # than that error bar cannot carry a percentage: VEGF-A's is 0.029
        # against a wobble of 0.027, and its design was reported 100%
        # selective off the back of it. The absolute 0.02 cutoff below was
        # the right idea with the wrong units -- the wobble spans 0.003 to
        # 0.077 across this catalogue -- so both tests apply.
        face_len = str(len(face))
        noise = (getattr(target, "panel_ceiling_noise", {}) or {}).get(face_len)
        floor = INDISTINGUISHABLE_MARGIN
        noise_limited = False
        if noise:
            noise_floor = CEILING_NOISE_FACTOR * noise
            if noise_floor > floor:
                floor, noise_limited = noise_floor, True

        # The stored ceiling, where there is one, rather than the live scale:
        # achievable_margin_range() substitutes lo + 0.1 when the sweep finds
        # no positive margin at all, which keeps the ramp usable but hides
        # the very case this test exists to catch. On a 10-residue face no
        # composition beats CD33's Siglecs -- the true ceiling is -0.012 --
        # and the substituted value would have waved it through.
        ceiling = (getattr(target, "panel_ceiling", {}) or {}).get(face_len)
        if ceiling is None:
            ceiling = hi

        if ceiling <= floor:
            # No face can separate this target from its closest relative --
            # the two present nearly the same surface. Reporting a score
            # here is not merely imprecise, it inverts the finding: the ramp
            # maps a best-possible margin of zero onto 100%, so the targets
            # where selectivity is IMPOSSIBLE came out looking perfectly
            # selective. Alpha-synuclein against beta-synuclein and IL-17A
            # against IL-17F both did exactly that.
            # Name the relative whose SURFACE most resembles the target's,
            # not whichever happened to win against this one face. What makes
            # the distinction impossible is a property of the two proteins,
            # so the explanation has to point at the protein responsible.
            def surface_identity(s: str) -> float:
                n = min(len(s), len(target.epitope_sequence))
                if not n:
                    return 0.0
                return sum(1 for a, b in zip(s, target.epitope_sequence)
                           if a == b) / n

            worst_name = max(surrogates, key=lambda g: surface_identity(surrogates[g]))
            twin = surrogates[worst_name]
            same = sum(1 for a, b in zip(twin, target.epitope_sequence)
                       if a == b)
            return ScoreComponent(
                name=self.name,
                raw=margin,
                normalised=0.5,
                weight=self.weight,
                ranked=False,
                units="kernel margin",
                # Two different reasons for the same silence, and saying the
                # wrong one is a claim about the proteins that the evidence
                # does not support. VEGF-A's nearest mapped surface is 4/17
                # identical -- nothing like a twin -- and the limit there is
                # the measurement, not the biology.
                detail=(
                    (f"This comparison is below its own resolution. The best "
                     f"margin any composition reaches against this panel is "
                     f"{ceiling:+.3f}, and that ceiling moves by {noise:.3f} "
                     f"on its own if a single one of the aligned residues in "
                     f"the panel is misassigned -- which at this level of "
                     f"identity is a certainty, not a worry. A percentage "
                     f"read off a ceiling that size would be a percentage of "
                     f"noise. The nearest mapped surface is {worst_name} at "
                     f"{same}/{len(target.epitope_sequence)} residues "
                     f"identical ({target.epitope_sequence} against {twin}), "
                     f"so what limits this is the measurement rather than any "
                     f"claim that the two proteins are alike."
                     if noise_limited else
                     f"No binding face can separate {target.gene} from this "
                     f"panel: the best margin any composition reaches is "
                     f"{ceiling:+.3f}, which is indistinguishable from zero. "
                     f"The most similar surface is {worst_name}, "
                     f"{same}/{len(target.epitope_sequence)} residues "
                     f"identical ({target.epitope_sequence} against {twin}). "
                     f"Whether the surfaces are genuinely too alike or this "
                     f"composition-based measure simply cannot see the "
                     f"difference, it has no signal to offer here, so nothing "
                     f"is scored.")),
                note=(f"reachable margin is within its own error bar; "
                      f"excluded from ranking" if noise_limited else
                      f"selectivity against {worst_name} is unreachable by "
                      f"composition; excluded from ranking"),
                caveat=(
                    "A real antibody can still discriminate here, because it "
                    "reads shape and backbone conformation that this "
                    "composition-based metric cannot see. What this says is "
                    "that THIS method has no signal to offer on the "
                    "distinction, not that the distinction cannot be made."),
            )

        norm = margin_to_score(margin, lo, hi)

        runners = ", ".join(f"{n} {v:.2f}" for n, v in ranked[1:3])
        direction = (
            f"This face prefers {worst_name} to {target.gene}: the margin is "
            f"negative, so it is not weakly selective but pointed at the wrong "
            f"protein of the two. Scored in the bottom quarter by "
            f"construction. "
            if margin < 0 else "")

        return ScoreComponent(
            name=self.name,
            raw=margin,
            normalised=norm,
            weight=self.weight,
            units="kernel margin",
            detail=(
                f"{direction}"
                f"complementarity {on_target:.3f} to {target.gene} against "
                f"{len(surrogates)} paralogs; closest competitor is "
                f"{worst_name} at {worst_val:.3f}, margin {margin:+.3f} "
                f"(0 to {hi:+.3f} is the reachable selective range; below 0 "
                f"is counter-selective). Next: {runners or 'none'}."
                + (f" The ceiling is {hi / noise:.1f} times its own error bar "
                   f"of {noise:.3f}"
                   + (", which is close enough to the line that the ranking "
                      "here deserves less weight than the number suggests."
                      if hi < 3.0 * noise else ".")
                   if noise else "")
            ),
            note=(f"counter-selective: prefers {worst_name}"
                  if margin < 0 else ""),
            # Not "low on this axis" but "aimed at the wrong protein", which
            # is a different claim and must not be tradeable against
            # manufacturability. See ScoreCard.rank_key.
            disqualifying=(
                f"complements {worst_name} more than {target.gene} "
                f"(margin {margin:+.3f})" if margin < 0 else ""),
            caveat=(
                "This number is scaled to what is REACHABLE for this target, "
                "so it is not comparable between targets: 90 here means 'near "
                "the best margin this epitope and this panel allow', and the "
                "best allowed margin is itself larger when the relatives are "
                "distant. A high score against a panel of distant relatives "
                "is an easier result than a moderate score against close "
                "ones. "
                "Compares composition at aligned surface positions, not docked "
                "geometry. A designed face can score ABOVE the reference "
                "antibody here without being more selective: this metric "
                "rewards matching a surface's composition, while an antibody's "
                "real selectivity comes from shape complementarity it cannot "
                "see. It also does NOT predict expression on healthy cells: "
                "whether an off-target protein sits on a cell that matters is "
                "measured by flow cytometry on primary human samples, and it "
                "is the reason CLEC12A was chosen over CD33 in the first "
                "place -- the two are not close relatives at all."
            ),
        )


class PolyspecificityScorer(Scorer):
    """Flags the chemistry associated with sticking to everything.

    In antibody engineering, patches of positive charge combined with exposed
    hydrophobicity predict nonspecific binding to membranes, heparin and
    off-target cells, and they correlate with fast clearance in vivo. This is
    one of the better-established sequence-level liabilities.
    """

    name = "polyspecificity_risk"
    weight = 1.0

    def score(self, sequence: str, context: dict) -> ScoreComponent:
        q = net_charge(sequence, 7.4)
        n = max(len(sequence), 1)
        frac_hydro = sum(1 for a in sequence if a in HYDROPHOBIC) / n
        frac_basic = sum(1 for a in sequence if a in "KR") / n

        patch = self._worst_basic_patch(sequence)

        # Strongly positive overall charge is the main risk factor.
        charge_term = 1.0 if q <= 1.0 else max(0.0, 1.0 - (q - 1.0) / 8.0)
        # Exposed hydrophobic fraction above ~0.45 gets sticky.
        hydro_term = window_score(frac_hydro, 0.15, 0.42, soft=0.18)
        # The product is what really drives promiscuity.
        combo = frac_basic * frac_hydro
        combo_term = 1.0 - ramp_score(combo, 0.03, 0.14)
        # Contiguous basic patches bind heparin and cell membranes.
        patch_term = 1.0 if patch <= 2 else max(0.0, 1.0 - (patch - 2) * 0.25)

        norm = (0.3 * charge_term + 0.25 * hydro_term
                + 0.25 * combo_term + 0.2 * patch_term)

        return ScoreComponent(
            name=self.name,
            raw=combo,
            normalised=norm,
            weight=self.weight,
            units="basic x hydrophobic",
            detail=(
                f"net charge {q:+.2f}, {frac_hydro:.0%} hydrophobic, "
                f"{frac_basic:.0%} basic, longest basic patch {patch}. "
                f"Risk rises with positive charge combined with exposed "
                f"hydrophobic surface."
            ),
            caveat=(
                "A developability heuristic from antibody engineering, applied "
                "to a sequence rather than a folded surface."
            ),
        )

    @staticmethod
    def _worst_basic_patch(seq: str) -> int:
        best = run = 0
        for a in seq:
            run = run + 1 if a in "KRH" else 0
            best = max(best, run)
        return best


# Gates on mapping an epitope onto another protein. A surrogate surface built
# from a weak alignment is not evidence about anything, and it would be scored
# with the same confidence as a good one.
MIN_MAP_ALIGN_SCORE = 20
MIN_MAPPED_POSITIONS = 5
WINDOW_PAD = 30
# Best-achievable margin at or below which the target and its closest
# relative are indistinguishable by composition. Not a tuning knob: at zero
# there is literally no face that complements one better than the other.
INDISTINGUISHABLE_MARGIN = 0.02
# A margin below zero means the face complements the nearest RELATIVE better
# than the target -- the design is not weakly selective, it is pointed the
# wrong way. Those designs used to be scored on a ramp whose floor was the
# 10th percentile of random faces, so merely beating a random face read as a
# majority-good score: the optimised TNF-alpha design reached a margin of
# -0.153 and was reported as 65.6% selective, and the EGFR design reached
# -0.011 and was reported as 70.7%. The sign of the measurement was being
# thrown away by its own scale.
#
# Zero is not an arbitrary point on this axis, so the whole region below it is
# compressed into the bottom quarter of the score and labelled. The band is
# kept rather than flattened to 0 because the optimiser needs a gradient to
# climb out of it, and a score pinned at one value is the failure this
# codebase has already hit five times.
COUNTER_SELECTIVE_BAND = 0.25
# A ceiling is worth reading only if it is at least this many standard
# deviations above zero, measured by ceiling_noise(). Two sigma is the
# ordinary convention for "distinguishable from nothing", not a figure chosen
# to make any particular target pass or fail -- and it does leave IL-1beta
# and CD38 sitting just the far side of the line, which is reported rather
# than smoothed over.
CEILING_NOISE_FACTOR = 2.0


def numbering_offset(target) -> tuple[int, float] | None:
    """Find the shift between structure numbering and reference numbering.

    Crystallographers number by the construct, which frequently differs from
    the UniProt sequence by a constant -- most often the length of a signal
    peptide that is cleaved off and so never appears in the structure. HER2
    is the clear case here: its epitope residues are numbered for the mature
    protein and sit 22 positions below their UniProt index.

    Assuming the two numberings agree, which is what this did at first,
    silently refuses to map the epitope of every target whose depositors
    numbered it the other way. Searching for the offset that makes the most
    epitope residues match their own amino acid recovers them, and it is
    self-checking: a wrong offset will not reproduce the right residue
    identities at a dozen separate positions.

    Returns (offset, fraction_matched) or None.
    """
    full = getattr(target, "full_sequence", "")
    if not full or not target.epitope:
        return None

    pairs = [(r.resseq - 1, r.aa) for r in target.epitope]
    span = max(len(full), max(i for i, _ in pairs) + 1)
    best, best_hits = 0, -1
    for off in range(-span, span + 1):
        hits = sum(1 for i, aa in pairs
                   if 0 <= i + off < len(full) and full[i + off] == aa)
        if hits > best_hits:
            best, best_hits = off, hits
        if hits == len(pairs):
            break
    return (best, best_hits / len(pairs)) if best_hits > 0 else None


def epitope_window(target) -> tuple[str, list[int]] | None:
    """The epitope plus flanking context, and where the epitope sits in it.

    Returns (window_sequence, positions_within_window), or None when the
    structure's residue numbering cannot be reconciled with the reference
    sequence -- in which case nothing downstream may map anything, because
    the positions would be pointing at the wrong residues.

    Shared deliberately with the panel builder in scripts/build_panels.py.
    If a panel were selected using one window and scored using another, the
    members would have been chosen for a similarity that is not the one being
    measured, which is the kind of mismatch that produces a confident number
    about nothing.
    """
    full = getattr(target, "full_sequence", "")
    if not full or not target.epitope:
        return None

    found = numbering_offset(target)
    if found is None:
        return None
    offset, matched = found
    if matched < 0.6:
        # The two numberings cannot be reconciled at all. Refusing here is
        # the point: a mapping built on positions that do not hold the right
        # residues would compare the wrong patch of surface.
        return None

    verified = [r.resseq - 1 + offset for r in target.epitope
                if 0 <= r.resseq - 1 + offset < len(full)
                and full[r.resseq - 1 + offset] == r.aa]
    if len(verified) < len(target.epitope) * 0.6:
        return None

    lo = max(0, min(verified) - WINDOW_PAD)
    hi = min(len(full), max(verified) + WINDOW_PAD + 1)
    return full[lo:hi], [i - lo for i in verified]


@dataclass
class MappedSurface:
    """One protein's surface at the positions equivalent to the epitope."""

    surface: str
    score: float        # raw Smith-Waterman score; grows with alignment length
    identity: float     # fraction identical over the aligned region
    coverage: float     # fraction of the epitope window that aligned at all


def map_epitope_onto(window: str, positions: list[int], sequence: str,
                     ) -> MappedSurface | None:
    """Read off the residues of `sequence` that align to the epitope.

    Returns None when the alignment is too weak or too gappy to mean
    anything. Positions that align to a gap are skipped rather than filled in
    with a guess.

    Note what `score` is and is not. Smith-Waterman score grows with the
    length of the alignment, so ranking relatives by it hands the top of
    every panel to whichever family member is longest. Ranking PD-L1's
    relatives that way returned MAG, nephrin and perlecan -- enormous
    multi-domain proteins that happen to contain an immunoglobulin fold --
    ahead of PD-L2, which is the actual nearest paralog. Rank on `identity`.
    """
    if not window or not sequence:
        return None
    aln = smith_waterman(window, sequence)
    if aln.score < MIN_MAP_ALIGN_SCORE:
        return None
    chars = [sequence[aln.mapping[p]] for p in positions if p in aln.mapping]
    if len(chars) < MIN_MAPPED_POSITIONS:
        return None
    paired = sum(1 for a, b in zip(aln.query_aligned, aln.target_aligned)
                 if a != "-" and b != "-")
    return MappedSurface("".join(chars), aln.score, aln.identity,
                         paired / len(window))


def build_paralog_epitopes(target, off_targets) -> dict[str, str]:
    """Map the epitope positions onto each paralog by local alignment.

    Returns {gene: surrogate_epitope_sequence}.

    scripts/build_panels.py has already done this and written the answer into
    panels.json, so the alignments here reproduce a stored result exactly --
    at a cost of up to a second per Evaluator, which made the test suite
    twenty minutes long and every design run a second slower than it needed
    to be.

    The cache is used only when scripts/verify_panel_cache.py has re-derived
    every surface and confirmed it matches, which it records by storing the
    epitope the surfaces belong to. Rebuilding the catalogue can move a
    binding site, and a surrogate mapped onto the old one would be a
    confident wrong answer rather than a visible failure; when the recorded
    epitope is absent or does not match the target's, the alignments run as
    before.
    """
    cached = getattr(target, "panel_surfaces", None)
    recorded = getattr(target, "panel_epitope", "")
    if cached and recorded and recorded == target.epitope_sequence:
        wanted = {ot.gene for ot in off_targets
                  if getattr(ot, "structural", True)}
        # Every structural off-target must be present, not merely every one
        # that happens to be cached. The first version of this compared the
        # filtered hit against its own intersection, which is true by
        # construction and therefore checked nothing: CLEC12A's panel
        # silently dropped from nine members to eight, because a curated
        # member has no stored surrogate, and its reachable ceiling moved
        # from 0.350 to 0.359 without a word.
        if wanted <= set(cached):
            return {g: cached[g] for g in wanted}

    got = epitope_window(target)
    if got is None:
        # Numbering does not line up; do not fabricate a mapping.
        return {}
    window, positions = got

    out: dict[str, str] = {}
    for ot in off_targets:
        # Only genuine structural relatives get the epitope mapped across.
        # Aligning onto a different fold family (CD33 is a Siglec) produces a
        # correspondence with no structural meaning, and a confident number
        # from a meaningless alignment is worse than no number.
        if not getattr(ot, "structural", True):
            continue
        mapped = map_epitope_onto(window, positions, ot.sequence)
        if mapped is not None:
            out[ot.gene] = mapped.surface
    return out
