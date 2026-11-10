"""The composite evaluator and the reference anchor.

What the headline number means
------------------------------
The composite is a weighted average of the components in this package, each
of which measures one named, published criterion. It is a *design triage*
score: it answers "is this worth making before the others", which is the
question a real campaign asks when it can only express twenty candidates.

It is NOT a probability of clinical benefit, and it is not an affinity
prediction. To stop the number floating free of reality, every run also
scores a real molecule -- a single-chain Fv built from the 50C1 antibody that
is crystallised against this exact epitope in PDB 8W9J -- through the
identical pipeline. Candidates are then reported both absolutely and as a
percentage of that reference.

Reaching 100% of the reference means "scores as well as a real anti-CLEC12A
antibody on the criteria a laptop can compute". It does not mean the molecule
binds, and it does not mean it helps anyone. The gap between those two
statements is the entire reason wet labs exist.
"""

from __future__ import annotations

import json

from dataclasses import dataclass, field
from pathlib import Path

from ..design.scaffolds import Scaffold, antibody_vh_reference
from ..targets import TargetSpec, load_off_targets
from .base import ScoreCard, Scorer
from .biophysics import (
    ManufacturabilityScorer, ScaffoldFitScorer, SizeScorer, SolubilityScorer,
    StabilityScorer,
)
from .interface import (
    CationPiScorer, ChargeComplementarityScorer, FootprintScorer,
    HydrogenBondScorer, ParatopeMimicryScorer, SaltBridgeCapacityScorer,
)
from .plm import PLMPlausibilityScorer
from .selectivity import (
    ParalogSelectivityScorer, PolyspecificityScorer, build_paralog_epitopes,
)

# Which group each component belongs to, for the grouped breakdown.
GROUPS = {
    "recognition": [
        "salt_bridge_capacity", "cation_pi_potential", "charge_complementarity",
        "hbond_complementarity", "paratope_character", "footprint_match",
    ],
    "selectivity": ["paralog_selectivity", "polyspecificity_risk"],
    "developability": [
        "solubility", "stability", "manufacturability", "size_suitability",
        "scaffold_consistency",
    ],
    "plausibility": ["plm_plausibility"],
}

SCFV_LINKER = "GGGGSGGGGSGGGGS"

# Random library members measured once per scaffold to set the upper end of
# the plausibility scale. Twelve is enough to place the mean; each one costs a
# few hundred milliseconds on a GPU.
PLM_ANCHOR_DECOYS = 12

# A component whose score varies less than this across the library is treated
# as constant: reported, but not used for ranking.
FLAT_COMPONENT_SD = 0.04

# Components that compare a candidate against the reference binder. Scoring
# the reference with these is circular -- it necessarily scores perfectly
# against itself -- so they are excluded when the reference's own group
# scores are computed. Leaving them in inflated the reference's recognition
# to exactly 100.0, which made every candidate percentage look worse than it
# was for a bookkeeping reason rather than a scientific one.
SELF_REFERENTIAL = {"paratope_character"}


@dataclass
class ReferenceAnchor:
    """A real molecule, scored by the same pipeline, to anchor the scale."""

    name: str
    sequence: str
    card: ScoreCard
    description: str
    group_scores: dict[str, float] = field(default_factory=dict)

    @property
    def score_100(self) -> float:
        return self.card.score_100

    @property
    def perplexity(self) -> float | None:
        c = self.card.component("plm_plausibility")
        if c is None or c.raw != c.raw:  # NaN check
            return None
        return float(c.raw)


# Built reference constructs, keyed by everything they depend on. Building
# one parses the whole structure -- seconds for a large assembly -- and the
# answer is a pure function of the entry, so it is computed once and kept.
# The key includes the epitope sequence and the contact list, which is what
# makes staleness impossible rather than merely checked: change the binding
# site and the key changes with it.
REFERENCE_CACHE = (Path(__file__).resolve().parents[3]
                   / "data" / "catalog" / "references.json")
_REF_MEM: dict[str, dict] | None = None


def _reference_key(target: TargetSpec) -> str:
    import hashlib
    raw = "|".join([
        target.key, target.complex_pdb or "", target.target_chain or "",
        ",".join(target.binder_chains or []), target.epitope_sequence,
        str(len(target.contacts or [])),
    ])
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:24]


def _reference_cache() -> dict[str, dict]:
    global _REF_MEM
    if _REF_MEM is None:
        try:
            _REF_MEM = json.loads(REFERENCE_CACHE.read_text(encoding="utf-8"))
        except Exception:
            _REF_MEM = {}
    return _REF_MEM


def _store_reference(key: str, value: tuple) -> None:
    cache = _reference_cache()
    cache[key] = {"sequence": value[0], "contacts": value[1],
                  "provenance": value[2], "format": value[3]}
    try:
        REFERENCE_CACHE.parent.mkdir(parents=True, exist_ok=True)
        REFERENCE_CACHE.write_text(json.dumps(cache, indent=1),
                                   encoding="utf-8")
    except Exception:
        pass  # a cache that cannot be written is a slow run, not a wrong one


def build_reference_scfv(target: TargetSpec,
                         struct_dir: Path) -> tuple[str, list[int], str, str] | None:
    """Build the molecule that is actually bound in the reference structure.

    Returns (sequence, contact_indices, provenance, format) or None.

    Two cases, because catalogue targets differ in what is bound to them:

    * An antibody is rebuilt as a single-chain Fv (VH + linker + VL). That is
      a real therapeutic format -- it is what sits on the front of a CAR --
      so comparing a designed miniprotein against it compares two things
      someone could actually build.
    * A protein partner (ACE2 on the SARS-CoV-2 spike, for instance) is taken
      as it is. It is already a real molecule that engages this surface.

    Either way the point is the same: anchor the scale to something that
    demonstrably binds, scored through the identical pipeline.
    """
    from ..pdbio import parse_structure

    cache_key = _reference_key(target)
    hit = _reference_cache().get(cache_key)
    if hit:
        return (hit["sequence"], list(hit["contacts"]),
                hit["provenance"], hit["format"])

    # Either format. An mmCIF-only entry -- 9EZA and 9XN4 among the
    # catalogue's own targets -- has no .pdb file, and looking only for one
    # left those targets with no reference molecule to anchor their scores
    # against.
    pdb_path = struct_dir / f"{target.complex_pdb}.pdb"
    if not pdb_path.exists():
        pdb_path = struct_dir / f"{target.complex_pdb}.cif"
    if not pdb_path.exists():
        return None

    struct = parse_structure(pdb_path)
    roles = target.chain_roles or {}
    touching = set(target.binder_chains or [])
    heavy = next((c for c, r in roles.items()
                  if r == "antibody_heavy" and c in touching), None)
    light = next((c for c, r in roles.items()
                  if r == "antibody_light" and c in touching), None)

    def contact_map(chain_residues: dict[str, list]) -> list[int]:
        """Indices of the crystallographic contact residues in the construct."""
        index: dict[tuple[str, int], int] = {}
        offset = 0
        for cid, residues in chain_residues.items():
            for i, r in enumerate(residues):
                index[(cid, r.resseq)] = offset + i
            offset += len(residues) + len(SCFV_LINKER)
        out, seen = [], set()
        for c in target.contacts:
            label = c.get("binder", "")
            if ":" not in label:
                continue
            cid, rest = label.split(":", 1)
            num = "".join(ch for ch in rest if ch.isdigit())
            if not num:
                continue
            key = (cid, int(num))
            if key in seen or key not in index:
                continue
            seen.add(key)
            out.append(index[key])
        return sorted(out)

    # --- antibody: rebuild as an scFv ---
    if heavy and light:
        vh_res = [r for r in struct.chain(heavy)
                  if r.is_standard and r.resseq <= 120]
        vl_res = [r for r in struct.chain(light)
                  if r.is_standard and r.resseq <= 110]
        if len(vh_res) >= 70 and len(vl_res) >= 70:
            vh = "".join(r.one_letter for r in vh_res)
            vl = "".join(r.one_letter for r in vl_res)
            seq = vh + SCFV_LINKER + vl
            idx = contact_map({heavy: vh_res, light: vl_res})
            prov = (
                f"scFv assembled from PDB {target.complex_pdb}: chain {heavy} "
                f"residues 1-120 (VH) + {SCFV_LINKER} linker + chain {light} "
                f"residues 1-110 (VL). {len(idx)} crystallographic contact "
                f"residues located in the construct."
            )
            result = (seq, idx, prov, "scFv")
            _store_reference(cache_key, result)
            return result

    # --- any other bound partner: take it as it is ---
    partner_chains = [c for c in struct.chain_ids
                      if c in touching and c != target.target_chain]
    residues = {c: [r for r in struct.chain(c) if r.is_standard]
                for c in partner_chains}
    residues = {c: r for c, r in residues.items() if len(r) >= 8}
    if not residues:
        return None

    seq = SCFV_LINKER.join("".join(r.one_letter for r in rs)
                           for rs in residues.values())
    idx = contact_map(residues)
    kind = "antibody fragment" if (heavy or light) else "protein partner"
    prov = (
        f"{target.reference_binder_name} taken from PDB "
        f"{target.complex_pdb}, chain(s) {','.join(residues)}"
        + (f" joined by a {SCFV_LINKER} linker" if len(residues) > 1 else "")
        + f". {len(idx)} crystallographic contact residues located."
    )
    result = (seq, idx, prov, kind)
    _store_reference(cache_key, result)
    return result


class Evaluator:
    """Assembles the scorers and produces ScoreCards."""

    def __init__(self, target: TargetSpec, *, use_plm: bool = True,
                 struct_dir: Path | None = None,
                 off_targets: list | None = None) -> None:
        self.target = target
        # Supplied directly for a target built in this session, whose panel
        # is deliberately not written into the curated catalogue file that
        # load_off_targets reads.
        self._given_off_targets = off_targets
        self.struct_dir = struct_dir or (
            Path(__file__).resolve().parents[3] / "data" / "structures")

        self.scorers: list[Scorer] = [
            SaltBridgeCapacityScorer(),
            CationPiScorer(),
            ChargeComplementarityScorer(),
            HydrogenBondScorer(),
            ParatopeMimicryScorer(),
            FootprintScorer(),
            ParalogSelectivityScorer(),
            PolyspecificityScorer(),
            SolubilityScorer(),
            StabilityScorer(),
            ManufacturabilityScorer(),
            SizeScorer(),
            ScaffoldFitScorer(),
        ]

        self.plm: PLMPlausibilityScorer | None = None
        if use_plm:
            self.plm = PLMPlausibilityScorer()
            self.scorers.append(self.plm)

        # Paralog surfaces, mapped once. Most catalogue targets have no
        # curated panel; the selectivity scorer then reports that it could
        # not assess selectivity rather than assuming the candidate is clean.
        self.off_targets = (self._given_off_targets
                            if self._given_off_targets is not None
                            else load_off_targets(target))
        self.paralog_epitopes = build_paralog_epitopes(target, self.off_targets)

        self._reference: ReferenceAnchor | None = None
        self._reference_built = False
        self._plm_anchors: dict[str, tuple[float, float] | None] = {}
        self._sel_scale: dict[int, tuple[float, float]] = {}
        # Components that do not vary across the library being searched.
        self._flat: dict[str, dict[str, float]] = {}
        self._active_flat: set[str] = set()
        self._measuring_variance = False

    # --- the reference ----------------------------------------------------

    @property
    def reference(self) -> ReferenceAnchor | None:
        if not self._reference_built:
            self._reference_built = True
            self._reference = self._build_reference()
        return self._reference

    def _build_reference(self) -> ReferenceAnchor | None:
        built = build_reference_scfv(self.target, self.struct_dir)
        if built is None:
            return None
        seq, para_idx, provenance, fmt = built

        scaffold = antibody_vh_reference(seq)
        scaffold.name = f"reference_{self.target.gene}"
        scaffold.interface_positions = para_idx
        # The band follows the reference's own format, so it is not penalised
        # for being the size that format is.
        scaffold.size_band = (len(seq) * 0.8, len(seq) * 1.25)
        scaffold.provenance = provenance

        card = self.evaluate(seq, scaffold, _is_reference=True)
        name = f"{self.target.reference_binder_name} ({fmt})"
        anchor = ReferenceAnchor(
            name=name,
            sequence=seq,
            card=card,
            description=(
                f"{self.target.reference_binder_name}, the molecule "
                f"crystallised against this exact surface in PDB "
                f"{self.target.complex_pdb}"
                + (f" at {self.target.source_resolution:.1f} A"
                   if self.target.source_resolution else "")
                + ". Scored by the identical pipeline so the scale is anchored "
                  "to something that demonstrably binds."
            ),
        )
        anchor.group_scores = self.group_scores(
            card, exclude_self_referential=True)
        return anchor

    def flat_components(self, scaffold: Scaffold | None,
                        n_samples: int = 150) -> dict[str, float]:
        """Components that cannot vary across this scaffold's library.

        On a fixed scaffold several criteria are constant by construction:
        every candidate is the same length, presents the same number of face
        residues, and shares the whole framework. Measured across 500
        sequences, size_suitability and footprint_match had a standard
        deviation of exactly zero, and four more sat near their ceiling --
        together about 27% of the composite's weight, contributing a nearly
        constant bonus that inflated every score and compressed the range.

        Such a component is still a real check on the molecule, so it is kept
        and displayed. It is only excluded from the RANKING, because a value
        that never varies cannot distinguish one candidate from another. On a
        scaffold where it does vary, it ranks again -- hence measuring rather
        than hard-coding a list.
        """
        if scaffold is None or not scaffold.variable_positions:
            return {}
        key = scaffold.name
        cached = self._flat.get(key)
        if cached is not None:
            return cached

        import random
        import statistics

        from ..aminoacids import AA20

        self._flat[key] = {}  # guard against re-entry
        self._measuring_variance = True
        try:
            rng = random.Random(4242)
            ss = scaffold.secondary_structure
            series: dict[str, list[float]] = {}
            for _ in range(n_samples):
                assign = {}
                for p in scaffold.variable_positions:
                    pool = [a for a in AA20 if a != "C"
                            and not (p < len(ss) and ss[p] == "H" and a in "PG")]
                    assign[p] = rng.choice(pool)
                card = self.evaluate(scaffold.apply(scaffold.framework, assign),
                                     scaffold)
                for c in card.components:
                    if c.weight > 0:
                        series.setdefault(c.name, []).append(c.normalised)
        finally:
            self._measuring_variance = False

        flat = {}
        for name, vals in series.items():
            sd = statistics.stdev(vals) if len(vals) > 1 else 0.0
            if sd < FLAT_COMPONENT_SD:
                flat[name] = sd

        self._flat[key] = flat
        self._active_flat = set(flat)
        return flat

    def selectivity_scale(self, face_length: int) -> tuple[float, float]:
        """The reachable range of paralog margins for a face of this length.

        Measured rather than assumed. The previous version hard-coded the
        ramp at -0.10 to +0.25 while the kernel could only reach +0.075, so
        every design scored exactly 50% no matter what it did. Sampling the
        real range keeps the component discriminating, and keeps it correct
        when paralogs are added to the panel.
        """
        if face_length <= 0:
            return (-0.10, 0.25)
        cached = self._sel_scale.get(face_length)
        if cached is not None:
            return cached
        from .selectivity import achievable_margin_range
        scale = achievable_margin_range(
            self.target.epitope_sequence, self.paralog_epitopes, face_length)
        self._sel_scale[face_length] = scale
        return scale

    # --- language-model scale --------------------------------------------

    def plm_anchors(self, scaffold: Scaffold | None) -> tuple[float, float] | None:
        """The plausibility scale for one scaffold, measured once and cached.

        good = the scaffold's natural parent sequence; poor = the mean of
        random library members on that same scaffold. Anchoring this way keeps
        the comparison inside one protein family. Anchoring it to the antibody
        reference instead rated the wild-type Z domain -- a real natural
        protein -- at 0%, and left every candidate saturated at zero, so the
        component carried weight while supplying no gradient.
        """
        if self.plm is None or scaffold is None:
            return None
        key = scaffold.name
        if key in self._plm_anchors:
            return self._plm_anchors[key]

        self._plm_anchors[key] = None  # guard against re-entry
        runner = self.plm.runner
        if not runner.load() or not scaffold.variable_positions:
            return None

        import random
        import statistics

        from ..aminoacids import AA20

        good = runner.pseudo_perplexity(scaffold.framework)
        if good is None:
            return None

        rng = random.Random(20240)
        ss = scaffold.secondary_structure
        vals: list[float] = []
        for _ in range(PLM_ANCHOR_DECOYS):
            assign = {}
            for p in scaffold.variable_positions:
                pool = [a for a in AA20 if a != "C"
                        and not (p < len(ss) and ss[p] == "H" and a in "PG")]
                assign[p] = rng.choice(pool)
            v = runner.pseudo_perplexity(
                scaffold.apply(scaffold.framework, assign))
            if v is not None:
                vals.append(v)

        if not vals:
            return None
        poor = statistics.fmean(vals)
        out = (good, poor) if poor > good else None
        self._plm_anchors[key] = out
        return out

    # --- evaluation -------------------------------------------------------

    def evaluate(self, sequence: str, scaffold: Scaffold | None = None,
                 *, _is_reference: bool = False) -> ScoreCard:
        face = scaffold.face(sequence) if scaffold else sequence
        if not face:
            face = sequence

        context = {
            "target": self.target,
            "scaffold": scaffold,
            "face": face,
            "paralog_epitopes": self.paralog_epitopes,
            "selectivity_scale": self.selectivity_scale(len(face)),
        }

        # Scaffold-anchored plausibility scale. Falls back to the reference
        # only when the scaffold cannot supply one (no designed positions).
        anchors = self.plm_anchors(scaffold)
        if anchors:
            context["plm_good"], context["plm_poor"] = anchors
        elif not _is_reference and self._reference is not None:
            ref_ppl = self._reference.perplexity
            if ref_ppl:
                context["reference_perplexity"] = ref_ppl

        card = ScoreCard(sequence=sequence)
        for scorer in self.scorers:
            try:
                card.add(scorer.score(sequence, context))
            except Exception as exc:
                card.flags.append(
                    f"{scorer.name} failed: {type(exc).__name__}: {exc}")

        # Drop constant components from the ranking (but keep reporting them).
        # Skipped while measuring variance, which would otherwise recurse.
        #
        # Skipped for the reference too, and not only for speed. "Constant
        # across this library" is a statement about the design library being
        # searched. The reference is one molecule; measuring variance by
        # mutating ITS positions invents a library that does not exist and
        # then annotates the reference's card with conclusions about it. The
        # comparison that matters, recognition_percent_of_reference, already
        # applies the DESIGN scaffold's flat set to both sides, so nothing
        # downstream wanted this. It cost 150 extra scorecards per target --
        # about ten seconds each, and the bulk of a twenty-minute test suite.
        if not self._measuring_variance and not _is_reference:
            flat = self.flat_components(scaffold)
            for c in card.components:
                if c.name in flat:
                    c.ranked = False
                    c.note = (
                        f"constant across this library (sd {flat[c.name]:.3f}) "
                        f"-- still measured and still a real check, but it "
                        f"cannot distinguish candidates, so it is excluded "
                        f"from the ranking score"
                    )

        card.metadata.update({
            "scaffold": getattr(scaffold, "name", None),
            "face": face,
            "designed_positions": (
                [p + 1 for p in scaffold.variable_positions] if scaffold else []),
            "designed_residues": (
                scaffold.variable_view(sequence) if scaffold else ""),
            "is_reference": _is_reference,
        })
        self._add_flags(card, sequence, scaffold)
        return card

    @staticmethod
    def _add_flags(card: ScoreCard, sequence: str,
                   scaffold: Scaffold | None) -> None:
        from ..aminoacids import max_homopolymer_run

        if sequence.count("C") % 2 == 1:
            card.flags.append(
                "odd number of cysteines: a free thiol will scramble disulfides")
        if max_homopolymer_run(sequence) >= 4:
            card.flags.append("contains a run of 4+ identical residues")
        for motif in ("DP", "NG", "DG"):
            if motif in sequence:
                card.flags.append(f"contains {motif}, a known degradation hotspot")
        if scaffold and scaffold.frozen_positions:
            for p in scaffold.frozen_positions:
                if p < len(sequence) and p < len(scaffold.framework):
                    if sequence[p] != scaffold.framework[p]:
                        card.flags.append(
                            f"frozen position {p + 1} was changed "
                            f"({scaffold.framework[p]} -> {sequence[p]}); "
                            f"the fold may not survive this")

    # --- grouped views ----------------------------------------------------

    @staticmethod
    def group_scores(card: ScoreCard,
                     exclude_self_referential: bool = False,
                     exclude_names: set[str] | None = None) -> dict[str, float]:
        """Weighted score per group, 0..1, over the ranking components only.

        Groups with no active component are omitted rather than reported as
        zero. `exclude_names` is used to hold both sides of a comparison to
        the same component set.
        """
        skip = set(SELF_REFERENTIAL) if exclude_self_referential else set()
        skip |= set(exclude_names or ())
        out: dict[str, float] = {}
        for group, names in GROUPS.items():
            comps = [c for c in card.components
                     if c.name in names and c.weight > 0 and c.ranked
                     and c.name not in skip]
            tw = sum(c.weight for c in comps)
            if tw > 0:
                out[group] = sum(c.normalised * c.weight for c in comps) / tw
        return out

    def percent_of_reference(self, card: ScoreCard) -> float | None:
        """Candidate composite as a percentage of the reference molecule's."""
        ref = self.reference
        if ref is None or ref.score_100 <= 0:
            return None
        return 100.0 * card.score_100 / ref.score_100

    def recognition_percent_of_reference(self, card: ScoreCard) -> float | None:
        """The same comparison restricted to the interface components.

        This is the more meaningful of the two numbers: it asks whether the
        candidate's binding face has the chemistry to engage this epitope as
        well as a real antibody's does, without the comparison being blurred
        by the two molecules' very different sizes and formats.
        """
        ref = self.reference
        if ref is None:
            return None
        # Both sides must use the same component set or the ratio is
        # meaningless: self-referential components dropped from both, and any
        # component that is constant for the design library dropped from the
        # reference too, even though it varies for an antibody.
        drop = self._active_flat
        mine = self.group_scores(card, exclude_self_referential=True,
                                 exclude_names=drop).get("recognition")
        theirs = self.group_scores(ref.card, exclude_self_referential=True,
                                   exclude_names=drop).get("recognition")
        if not mine or not theirs or theirs <= 0:
            return None
        return 100.0 * mine / theirs

    def calibrate(self, scaffold: Scaffold, *, n_decoys: int = 200,
                  seed: int = 1234):
        """Build the decoy null distribution for this scaffold."""
        from .calibration import build_calibration
        return build_calibration(self, scaffold, n_decoys=n_decoys, seed=seed)

    # --- diagnostics ------------------------------------------------------

    def status(self) -> str:
        lines = [
            f"Target: {self.target.gene} ({self.target.disease})",
            f"Scorers active: {sum(1 for s in self.scorers if s.available())}"
            f"/{len(self.scorers)}",
        ]
        if self.plm:
            lines.append(f"Language model: {self.plm.runner.status}")
        else:
            lines.append("Language model: disabled")
        if self.paralog_epitopes:
            for gene, epi in self.paralog_epitopes.items():
                lines.append(f"Paralog surface mapped: {gene} -> {epi}")
        else:
            lines.append(
                "Paralog surfaces: none mapped (selectivity scored neutral)")
        ref = self.reference
        if ref:
            lines.append(
                f"Reference anchor: {ref.name}, {len(ref.sequence)} aa, "
                f"composite {ref.score_100:.1f}/100")
        else:
            lines.append("Reference anchor: UNAVAILABLE (percentages suppressed)")
        return "\n".join(lines)
