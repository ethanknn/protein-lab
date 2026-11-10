"""Ranking catalogue targets against one specimen's antigen profile.

The input is what a flow cytometry report already says -- which surface
antigens the blast population carries -- plus which healthy lineages must
survive, plus any variants called in the target genes. The output is an
ordering over catalogue targets, with the reasoning attached.

This is target selection, not treatment selection. It says which protein a
lab should build a binder against first. It does not say whether a binder
would work, and it cannot: see the module docstring of scoring/selectivity.py
for what the design scores do and do not mean.

--- The ranking rule ------------------------------------------------------

Targets are ranked by HOW MUCH IS KNOWN about them, not by how clean they
look.

A naive ranking sorts by concerns found, ascending. That ranking puts a
target nobody has ever examined at the very top, because an unexamined target
has no concerns on record. It is the same error as reading nTPM 0.0 as
"absent" -- see lineage.py -- promoted from one cell of the table to the
ordering over the whole table, where it does more damage.

So a lineage positively CLEARED by cytometry counts in a target's favour, and
a lineage nobody has data for counts against it. Two targets with no recorded
concerns are not equivalent: the one with cytometry behind it is a better
first thing to build, and the other is mostly an absence of looking.

Disqualification is kept separate from ranking, as it is in
scoring/base.py: an antigen the blast does not carry is not a weak candidate
but a non-candidate, and a weighted score would let a strong result elsewhere
compensate for the one fact that makes the whole thing pointless.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .lineage import (CURATED_TISSUE_GENES, LINEAGE_BY_KEY, LineageCall,
                      SEVERITY_LABEL, TOPOLOGY_LABEL, assess_gene,
                      assess_tissues, decisive, is_human, load_expression,
                      load_topology, not_applicable, off_axis, phrase,
                      physiology_for)

# What a cytometry report says about one antigen on the blast population.
# "not_tested" is not a synonym for "negative" and the two must not collapse:
# one is a measurement, the other is its absence.
ANTIGEN_STATES = ("positive", "dim", "negative", "not_tested")

STATE_LABEL = {
    "positive": "positive on blasts",
    "dim": "dim on blasts",
    "negative": "negative on blasts",
    "not_tested": "not tested",
}


@dataclass(frozen=True)
class Variant:
    """A coding change called in a target gene, in reference numbering.

    `position` is 1-based against the UniProt sequence, which is what a
    sequencing report gives. Reconciling that with the crystal structure's
    own numbering is done in `check_variant`, never assumed.
    """

    gene: str
    position: int
    from_aa: str = ""
    to_aa: str = ""

    def label(self) -> str:
        if self.from_aa and self.to_aa:
            return f"{self.gene} p.{self.from_aa}{self.position}{self.to_aa}"
        return f"{self.gene} position {self.position}"


@dataclass
class Specimen:
    """One specimen's profile. Deliberately carries no identifiers.

    There is no name, date of birth or record number field, and nothing here
    is persisted between sessions. The tool reasons about a cell population,
    and giving it a patient's identity would add nothing to the reasoning
    while turning a design sandbox into a store of clinical records.
    """

    antigens: dict[str, str] = field(default_factory=dict)
    spared: tuple[str, ...] = ()
    variants: tuple[Variant, ...] = ()
    note: str = ""

    def state(self, gene: str) -> str:
        return self.antigens.get((gene or "").upper(), "not_tested")


# --- Variant check ----------------------------------------------------------


@dataclass
class VariantHit:
    variant: str
    hits_epitope: bool
    residue: str = ""
    distance: float | None = None
    detail: str = ""
    reference_mismatch: bool = False


def check_variant(target, variant: Variant) -> VariantHit:
    """Does this variant sit on a residue the reference binder touches?

    A binder is designed against the epitope of the REFERENCE sequence. If
    the specimen's copy of the gene differs at a residue in that epitope, the
    designed surface is complementary to a protein the specimen does not
    have. This is not hypothetical: epitope loss is an established route to
    relapse after CD19-directed therapy.

    The structure's residue numbering and the reference sequence's numbering
    routinely differ by a constant, so the offset is measured with
    numbering_offset() rather than assumed -- and if the variant's stated
    reference amino acid does not match the reference sequence at that
    position, that is reported instead of being worked around, because it
    means the two are not using the same numbering and no answer here would
    mean anything.
    """
    from .scoring.selectivity import numbering_offset

    full = getattr(target, "full_sequence", "") or ""
    found = numbering_offset(target)
    if found is None or not full:
        return VariantHit(
            variant=variant.label(), hits_epitope=False,
            detail="Cannot check: this target's structure numbering could "
                   "not be reconciled with its reference sequence, so a "
                   "position in one numbering cannot be located in the other.")

    offset, matched = found
    if matched < 0.6:
        return VariantHit(
            variant=variant.label(), hits_epitope=False,
            detail=f"Cannot check: structure and reference numbering agree "
                   f"on only {matched:.0%} of epitope residues, which is too "
                   f"little to place a variant position reliably.")

    i = variant.position - 1
    if variant.from_aa and (i < 0 or i >= len(full)
                            or full[i] != variant.from_aa):
        at = full[i] if 0 <= i < len(full) else "nothing"
        return VariantHit(
            variant=variant.label(), hits_epitope=False,
            reference_mismatch=True,
            detail=f"Position {variant.position} of {target.gene} is {at} in "
                   f"the reference sequence, not {variant.from_aa}. The "
                   f"variant and this catalogue entry are not using the same "
                   f"numbering, so no conclusion is drawn.")

    for r in target.epitope:
        if r.resseq + offset == variant.position:
            return VariantHit(
                variant=variant.label(), hits_epitope=True,
                residue=f"{r.aa}{variant.position}",
                distance=r.min_distance_to_binder,
                detail=(
                    f"This variant lands on {r.aa}{variant.position}, which "
                    f"is in the epitope, {r.min_distance_to_binder:.1f} A "
                    f"from the reference binder "
                    f"({target.reference_binder_name}). A binder designed "
                    f"against the reference surface is complementary to a "
                    f"residue this specimen does not have at that position, "
                    f"so the designed face may not fit. Whether it still "
                    f"binds depends on how much of the interface that one "
                    f"residue carries, which this cannot tell you."))

    near = min((abs(r.resseq + offset - variant.position)
                for r in target.epitope), default=None)
    return VariantHit(
        variant=variant.label(), hits_epitope=False,
        detail=(f"Not in the epitope"
                + (f"; the nearest epitope residue is {near} positions away"
                   if near is not None else "")
                + ". It may still change the protein's folding or its "
                  "surface at a distance, which is not assessed here."))


# --- Assessment -------------------------------------------------------------


@dataclass
class SpecimenCall:
    """One catalogue target weighed against one specimen."""

    key: str
    gene: str
    protein_name: str
    disease: str
    area: str
    antigen_state: str
    disqualifying: str = ""
    conflicts: list[LineageCall] = field(default_factory=list)
    cleared: list[LineageCall] = field(default_factory=list)
    unknowns: list[LineageCall] = field(default_factory=list)
    tissue_flags: list[LineageCall] = field(default_factory=list)
    # Non-haematopoietic liabilities established by tissue-level evidence,
    # split from tissue_flags because they are a different kind of claim: a
    # transcript flag says a gene is transcribed somewhere, while these say a
    # named organ was damaged in named patients. Kept apart from `conflicts`
    # too -- the clinician asked which blood populations must survive, and
    # never got asked about the heart.
    tissue_harms: list[LineageCall] = field(default_factory=list)
    tissue_spared: list[LineageCall] = field(default_factory=list)
    variant_hits: list[VariantHit] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    # Set only for targets the lineage question does not apply to -- the
    # viral and bacterial entries. See lineage.is_human.
    inapplicable: LineageCall | None = None
    # Set for targets that are not displayed on a cell surface at all: the
    # secreted and intracellular ones. Its own field for the same reason as
    # `inapplicable` -- a row saying the question does not apply must never
    # be counted among the lineages something positively cleared.
    off_axis: LineageCall | None = None
    topology: str = "unknown"
    # What the target normally does and what blocking it costs. Not a
    # lineage fact and not ranked on: the clinician asked which populations
    # must survive, not what the molecule is for. Carried because for a
    # target that is not on a cell, this is the only safety question left.
    physiology: object | None = None

    @property
    def disqualified(self) -> bool:
        return bool(self.disqualifying)

    @property
    def rank_key(self) -> tuple:
        """Lexicographic, mirroring ScoreCard.rank_key in scoring/base.py.

        Disqualified last, then by how much is positively known: more
        lineages cleared by cytometry is better, then fewer unknowns, then
        fewer soft concerns. There is no composite number here on purpose --
        averaging "three lineages cleared" against "two unknown" would let
        the second be bought off by the first, and they are not the same
        currency.

        `tissue_harms` sits between the two concern counts for a reason. An
        organ damaged in named patients is stronger evidence than a
        transcript flag and belongs above it; it is still weaker than a
        conflict with a population this clinician said must survive, because
        that constraint is the one they actually stated. So it outranks the
        screen and is outranked by the stated requirement.

        `tissue_spared` is deliberately NOT here. CLDN18.2's absence from the
        gastric stem cell zone is the best single fact in the tissue table,
        and it is still an answer about epithelium to a question about blood.
        It is carried for display and earns no rank credit.
        """
        return (
            0 if self.disqualified else 1,
            len(self.cleared),
            -len(self.unknowns),
            -len(self.conflicts),
            -len(self.tissue_harms),
            -len(self.tissue_flags),
        )

    def to_dict(self) -> dict:
        return {
            "key": self.key, "gene": self.gene,
            "protein_name": self.protein_name, "disease": self.disease,
            "area": self.area,
            "antigen_state": self.antigen_state,
            "antigen_label": STATE_LABEL.get(self.antigen_state,
                                             self.antigen_state),
            "disqualifying": self.disqualifying,
            "disqualified": self.disqualified,
            "conflicts": [c.to_dict() for c in self.conflicts],
            "cleared": [c.to_dict() for c in self.cleared],
            "unknowns": [c.to_dict() for c in self.unknowns],
            "tissue_flags": [c.to_dict() for c in self.tissue_flags],
            "tissue_harms": [c.to_dict() for c in self.tissue_harms],
            "tissue_spared": [c.to_dict() for c in self.tissue_spared],
            "inapplicable": (self.inapplicable.to_dict()
                             if self.inapplicable else None),
            "off_axis": (self.off_axis.to_dict() if self.off_axis else None),
            "topology": self.topology,
            "topology_label": TOPOLOGY_LABEL.get(self.topology,
                                                 self.topology),
            "physiology": ({
                "signal": self.physiology.signal,
                "cost": self.physiology.cost,
                "evidence": self.physiology.evidence,
                "source": self.physiology.source,
                "severity": self.physiology.severity,
                "severity_label": SEVERITY_LABEL.get(
                    self.physiology.severity, self.physiology.severity),
            } if self.physiology else None),
            "variant_hits": [{
                "variant": v.variant, "hits_epitope": v.hits_epitope,
                "residue": v.residue, "distance": v.distance,
                "detail": v.detail,
                "reference_mismatch": v.reference_mismatch,
            } for v in self.variant_hits],
            "notes": self.notes,
        }


def _assess_tissues_into(call, gene, expression) -> None:
    """Fill the three tissue buckets and the notes that go with them.

    Shared by both paths through assess_target. A secreted target still
    gets its tissue rows -- VEGFA's are the most informative thing on
    its card -- so this cannot live inside the branch that runs the
    lineage loop.
    """
    tissues = assess_tissues(gene, expression)
    curated_tissue = [c for c in tissues
                      if c.lineage.startswith("curated_tissue")]
    call.tissue_harms = [c for c in curated_tissue if c.harm == "observed"]
    call.tissue_spared = [c for c in curated_tissue if c.is_cleared]
    call.tissue_flags = [c for c in tissues
                         if not c.lineage.startswith("curated_tissue")
                         and c.level in ("moderate", "high")]

    # The sentence this feature was missing until the solid tumour antigens
    # went in.
    #
    # Everything above reasons about blood: the antigen panel comes from a
    # flow report on a marrow or blood specimen, and every population a
    # clinician can ask to spare is haematopoietic. For a blood cancer
    # antigen that is the right axis. For HER2 it is not. HER2's established
    # liability is cardiomyocytes, nectin-4's is skin, mesothelin's is the
    # lining of the pleura -- and a profile that comes back with no
    # haematopoietic conflicts has not cleared any of them. It has answered
    # a different question and said nothing about this one.
    #
    # Left unsaid, this is the worst failure mode available to this page,
    # because it fails QUIETLY and in the reassuring direction: an empty
    # conflict list reads as a clean target. The ranking already keeps these
    # antigens off the top by counting unknowns against them, but a rank is
    # not an explanation, and the clinician reading one card needs the
    # explanation.
    if call.tissue_harms:
        sites = ", ".join(c.label for c in call.tissue_harms)
        call.notes.append(
            f"{gene} has an established liability outside blood: {sites}. "
            f"Nothing in a flow panel on this specimen can speak to it, and "
            f"the sparing profile above does not try to -- every population "
            f"it offers is haematopoietic. Read the tissue rows as a separate "
            f"question that this specimen cannot answer, not as a milder "
            f"version of the rows above.")
    elif gene.upper() in CURATED_TISSUE_GENES:
        call.notes.append(
            f"{gene} has a curated non-haematopoietic finding, which the "
            f"sparing profile above cannot reach. The populations this tool "
            f"offers to spare are all blood populations.")
    elif call.area in ("Solid tumour", "Immuno-oncology") and not call.conflicts:
        # Only where the blood axis found nothing. KIT is why: it is a solid
        # tumour target whose liability genuinely IS a blood population, and
        # telling its card that "the liability for such targets usually sits
        # in tissue" while excluding it for carrying CD117 on the stem
        # compartment would have the page arguing with itself.
        article = "an" if call.area[0].lower() in "aeiou" else "a"
        call.notes.append(
            f"This is {article} {call.area.lower()} target, and the liability "
            f"for such targets usually sits in tissue rather than in blood. "
            f"Nothing is curated for {gene} outside blood, which means "
            f"nobody has recorded it here -- not that the tissue is clear. "
            f"The transcript rows below are flags to check, and they are the "
            f"weakest evidence on this page.")


def _check_variants_into(call, target, specimen) -> None:
    """Variants are a question about the epitope, not about topology,
    so every target gets this regardless of where the protein sits.
    """
    for v in specimen.variants:
        if v.gene.upper() != target.gene.upper():
            continue
        hit = check_variant(target, v)
        call.variant_hits.append(hit)
        if hit.hits_epitope:
            call.disqualifying = call.disqualifying or (
                f"{v.label()} changes a residue inside the epitope this "
                f"catalogue entry designs against, so a binder built from "
                f"the reference surface may not fit this specimen's protein.")


def assess_target(target, specimen: Specimen, expression=None,
                  topology=None) -> SpecimenCall:
    """Weigh one catalogue target against one specimen profile."""
    expression = expression if expression is not None else load_expression()
    topology = topology if topology is not None else load_topology()
    gene = target.gene
    state = specimen.state(gene)

    call = SpecimenCall(
        key=target.key, gene=gene, protein_name=target.protein_name,
        disease=target.disease, area=target.area, antigen_state=state)

    # Nothing to bind. Disqualifying rather than low-scoring: a binder cannot
    # engage an antigen the population does not display, and no amount of
    # selectivity elsewhere changes that.
    if state == "negative":
        call.disqualifying = (
            f"{gene} is reported negative on this specimen's blasts. There is "
            f"nothing for a binder to engage.")
    elif state == "not_tested":
        call.notes.append(
            f"{gene} was not in the antigen panel entered, so whether this "
            f"specimen carries it is unknown. Ranked on the sparing evidence "
            f"alone, which says nothing about whether the target is there.")
    elif state == "dim":
        call.notes.append(
            "Dim expression means fewer copies per cell for a binder to "
            "engage. Antigen density affects whether a given affinity is "
            "enough, and that threshold is measured in cell assays, not "
            "computed here.")

    if not is_human(target):
        # Its own field, not `cleared`. A non-human target has earned no
        # evidence credit -- nobody did any cytometry -- so counting it as a
        # cleared lineage would let it outrank a human target that somebody
        # actually examined, on the strength of a question not applying.
        call.inapplicable = not_applicable(target)
        call.notes.append(
            "Sparing healthy lineages is not the constraint for a non-human "
            "target. The constraint is the opposite one -- that the binder "
            "not engage a human protein by accident -- and that is what the "
            "paralog panel measures.")
        return call

    call.topology, why = topology.verdict(gene)
    call.physiology = physiology_for(gene)

    # The lineage loop runs only where its question has an answer.
    #
    # For a secreted or intracellular target it does not. The transcript
    # readings exist either way -- TNF reads 73.4 nTPM on monocytes, SNCA
    # 103.2 on dendritic cells -- and they measure which cells PRODUCE the
    # protein. Rendered through this loop they would come out as "carried on
    # monocytes", in the same words, the same pills and the same consequence
    # text as CD33 on progenitors. That is not a near-miss; it is a
    # different claim wearing the right answer's clothes, and the
    # consequence line would tell a clinician their monocytes are at risk
    # from a molecule that never touches a monocyte.
    #
    # So these targets get one honest row instead of eleven misleading ones,
    # exactly as the viral entries do.
    if call.topology in ("secreted", "intracellular"):
        call.off_axis = off_axis(gene, call.topology, why)
        call.notes.append(
            f"The populations this tool offers to spare are populations that "
            f"could carry a target on their surface. {gene} is "
            f"{TOPOLOGY_LABEL[call.topology]}, so that question was not put "
            f"to it and no absence of conflicts below should be read as one.")
        _assess_tissues_into(call, gene, expression)
        _check_variants_into(call, target, specimen)
        return call

    if call.topology == "shed":
        # Worth one line, because it changes what a binder has to do rather
        # than which cells are at risk. A circulating pool of the same
        # antigen is a sink: it binds the molecule before it reaches the
        # cell. Soluble mesothelin, soluble PD-L1 and soluble BCMA are all
        # measured as serum biomarkers for exactly that reason.
        call.notes.append(
            f"{gene} is shed as well as displayed, so a soluble pool of the "
            f"same antigen circulates. The lineage rows below still apply, "
            f"and a decoy that binds the molecule before it reaches a cell "
            f"is an additional problem they do not describe.")

    calls = {c.lineage: c for c in assess_gene(gene, expression)}
    for key in specimen.spared:
        c = calls.get(key)
        if c is None:
            continue
        if c.is_concern:
            call.conflicts.append(c)
            # Known harm to a population the clinician said must survive is
            # a reason to stop, not a penalty to trade off. Three conditions,
            # and the third is the one that took a correction to get right.
            #
            # Dense expression on a must-spare population disqualifies, but
            # only where nothing shows the population survives it. CD47 sits
            # on haematopoietic stem cells and blocking it still does not
            # increase phagocytosis of normal stem and progenitor cells --
            # the window there is a density difference, not an absence -- so
            # reading expression as destruction would rule out a target that
            # human and in vivo experience says is usable. "Nobody has
            # checked" is not that evidence and does not earn the exemption.
            if (c.level == "high" and decisive(c.evidence)
                    and c.harm != "not_observed"):
                lin = LINEAGE_BY_KEY.get(key)
                call.disqualifying = (
                    f"{gene} is carried at high level on "
                    f"{phrase(lin.label) if lin else key}, which this "
                    f"specimen's profile marks as a population that must "
                    f"survive. Evidence: {c.evidence}"
                    + ("; damage to this population has been observed in "
                       "patients" if c.harm == "observed" else
                       "; whether this population survives being engaged "
                       "has not been established")
                    + "." + (f" {lin.consequence}" if lin else ""))
            elif c.harm == "not_observed":
                call.notes.append(
                    f"{gene} is carried on "
                    f"{phrase(LINEAGE_BY_KEY[key].label) if key in LINEAGE_BY_KEY else key}, "
                    f"which this profile asks to spare, but that population "
                    f"is not destroyed by engaging it. Kept as a candidate "
                    f"and flagged rather than excluded.")
        elif c.is_cleared:
            call.cleared.append(c)
        else:
            call.unknowns.append(c)

    _assess_tissues_into(call, gene, expression)
    _check_variants_into(call, target, specimen)

    return call


def rank(targets, specimen: Specimen, expression=None,
         topology=None) -> list[SpecimenCall]:
    """Every target weighed, best first, disqualified ones kept and last.

    Disqualified entries are returned rather than filtered out. "CD33 was
    excluded because your profile says the marrow must recover" is the most
    useful sentence this produces, and dropping the row would delete it.
    """
    expression = expression if expression is not None else load_expression()
    topology = topology if topology is not None else load_topology()
    calls = [assess_target(t, specimen, expression, topology)
             for t in targets]
    calls.sort(key=lambda c: c.rank_key, reverse=True)
    return calls
