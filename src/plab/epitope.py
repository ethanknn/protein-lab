"""Derive a binding site from any protein-protein complex.

This is the generalised form of what scripts/analyze_epitope.py did for
CLEC12A. Given a structure containing a target and something bound to it, it
works out which chains are which, finds the largest interface, and returns
the contacted residues on both sides.

Why it can be general: nothing about the procedure is specific to one
protein. "The residues within 4.5 A of the bound partner" is the definition
of an epitope regardless of what the target is. The scorers downstream read
the epitope's chemistry out of this result rather than hard-coding any
particular surface, so pointing the pipeline at a new disease is a matter of
supplying a different complex.

What it cannot do: invent a binding site. If a structure has no second chain
in contact, there is no epitope to derive and this says so rather than
guessing. That honesty is what keeps the catalogue trustworthy.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from .pdbio import (
    Contact,
    Residue, Structure, find_contacts, interface_residues, parse_structure,
    radius_of_gyration,
)

CONTACT_CUTOFF = 4.5
# Low on purpose. Several real drug targets are short peptides rather than
# folded domains -- amyloid-beta in PDB 6CO3 is a 6-residue chain, and a
# 20-residue floor silently excluded it along with other peptide epitopes.
MIN_CHAIN_RESIDUES = 6
MIN_INTERFACE_RESIDUES = 5
# Above the floor but still thin. A typical antibody epitope covers 15-25
# residues; several real targets here are short peptides that cannot.
SMALL_INTERFACE_RESIDUES = 10
# Exact-k-mer identity is only ever used as a FAST PATH. It needs all k
# positions identical, so at sequence identity p the expected hit rate is
# p**k: a true ortholog at 90% identity scores 0.9**8 = 43% and a different
# viral strain of the same protein scores far lower. Using it as the decision
# rule rejected gp120 in 3NGB (75% identical, scored 13%) and mouse RANKL in
# 4GIQ (90% identical, scored 42%) -- both unambiguously the target protein.
IDENTITY_THRESHOLD = 0.5
SHORT_PEPTIDE_LENGTH = 25

# The real decision rule: aligned identity over most of the chain. The gap
# between a true ortholog and an unrelated chain is wide and these sit in it.
# Measured on this catalogue, chains that are NOT the target land at 10-22%
# identity (antibody chains against their antigen, a bacterial encapsulin
# against interleukin-5), while strain variants and mouse orthologs of the
# target land at 56-90%. A 45% floor is comfortably inside that gap.
ALIGNED_IDENTITY_THRESHOLD = 0.45
ALIGNED_COVERAGE_THRESHOLD = 0.70
# Below this, say so: the structure is a variant or ortholog of the human
# protein, not the human protein, and that is a real caveat about the design.
SAME_PROTEIN_IDENTITY = 0.90
# Contacts are read from side-chain positions, which are not reliably placed
# much beyond this.
RESOLUTION_WARN = 3.0
# Centroid distance within which two chains' contact patches count as the
# same site. Neighbouring residues in one contiguous epitope sit 5-9 A apart.
PATCH_LINKAGE = 8.0
# Two bound chains that are near-identical in sequence are two copies of one
# protein, so they are two molecules and must never be merged however closely
# they pack. Without this, gp120 in 3NGB picked up a second VRC01 heavy chain
# from the neighbouring copy and reported a 43-residue epitope for an
# antibody whose real footprint is about 25.
SAME_MOLECULE_MAX_IDENTITY = 0.8
# Residue-pair contacts two bound chains must share before they are called
# one molecule. A Fab's heavy and light chains bury a large interface against
# each other and are far above this; two separate Fabs that merely touch in
# the crystal are far below it.
MIN_PAIR_CONTACTS = 15

# Immunoglobulin J-segment signatures, used only to label antibody chains
# nicely in the report. Nothing depends on a binder being an antibody.
HEAVY_MOTIFS = ("WGQG", "WGRG", "WGAG", "WGKG", "WGPG", "WGMG")
LIGHT_MOTIFS = ("FGQG", "FGGG", "FGSG", "FGPG", "FGTG")


@dataclass
class EpitopeResidue:
    chain: str
    resseq: int
    icode: str
    resname: str
    aa: str
    label: str
    centroid: tuple[float, float, float]
    min_distance_to_binder: float


@dataclass
class EpitopeResult:
    """Everything derived from one complex."""

    pdb_id: str
    title: str
    resolution: float | None
    target_chain: str
    binder_chains: list[str]
    chain_roles: dict[str, str]
    epitope: list[EpitopeResidue]
    paratope_sequence: str
    epitope_sequence: str
    descriptors: dict
    contacts: list[dict]
    binder_kind: str          # "antibody" | "protein" | "unknown"
    target_identity: float    # identity of the matched chain to the reference
    target_coverage: float = 1.0
    identity_method: str = "exact k-mer"
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return len(self.epitope) >= MIN_INTERFACE_RESIDUES

    def to_dict(self) -> dict:
        d = asdict(self)
        d["epitope"] = [asdict(r) for r in self.epitope]
        return d

    def summary(self) -> str:
        return (
            f"{self.pdb_id}: {len(self.epitope)} epitope residues on chain "
            f"{self.target_chain} ({self.epitope_sequence}), bound by "
            f"{self.binder_kind} chains {','.join(self.binder_chains)}; "
            f"{self.descriptors.get('n_residue_contacts', 0)} contacts, "
            f"Rg {self.descriptors.get('epitope_radius_of_gyration', 0):.1f} A"
        )


class EpitopeError(RuntimeError):
    """Raised when no usable binding site can be derived."""


def kmer_identity(query: str, reference: str, k: int = 8) -> float:
    """Fraction of query k-mers found verbatim in reference.

    A fast path only. See the note on IDENTITY_THRESHOLD for why this must
    never be the thing that decides whether a chain is the target.
    """
    if len(query) < k or len(reference) < k:
        return 0.0
    ref = {reference[i:i + k] for i in range(len(reference) - k + 1)}
    hits = sum(1 for i in range(len(query) - k + 1) if query[i:i + k] in ref)
    return hits / (len(query) - k + 1)


@dataclass
class ChainMatch:
    """How well one chain matches the target, and by which measure."""

    identity: float
    coverage: float
    method: str          # "contained" | "exact k-mer" | "aligned"
    is_match: bool

    @property
    def same_protein(self) -> bool:
        """True when this is the target itself rather than a variant of it."""
        return self.identity >= SAME_PROTEIN_IDENTITY


def aligned_match(seq: str, target_sequence: str) -> ChainMatch:
    """Smith-Waterman identity and coverage of one chain against the target."""
    from .align import smith_waterman

    if not seq or not target_sequence:
        return ChainMatch(0.0, 0.0, "aligned", False)
    aln = smith_waterman(seq, target_sequence)
    paired = sum(1 for a, b in zip(aln.query_aligned, aln.target_aligned)
                 if a != "-" and b != "-")
    coverage = paired / len(seq)
    ok = (aln.identity >= ALIGNED_IDENTITY_THRESHOLD
          and coverage >= ALIGNED_COVERAGE_THRESHOLD)
    return ChainMatch(aln.identity, coverage, "aligned", ok)


def match_target(seq: str, target_sequence: str,
                 identity_threshold: float = IDENTITY_THRESHOLD) -> ChainMatch:
    """Is this chain the target protein, and how certain is that?

    Three tests, cheapest first:

      contained    a short peptide that the full-length protein contains --
                   exactly the relationship between amyloid-beta and APP, or
                   a circumsporozoite epitope peptide and PfCSP
      exact k-mer  an exact sequence match, which is the common case
      aligned      Smith-Waterman identity and coverage, which is what
                   catches orthologs, cleaved fragments and strain variants

    The third test exists because the first two cannot distinguish "a
    different protein" from "the same protein from a different species or
    virus strain". That distinction is the whole question here.
    """
    if not seq:
        return ChainMatch(0.0, 0.0, "contained", False)

    if len(seq) <= SHORT_PEPTIDE_LENGTH:
        if seq in target_sequence:
            return ChainMatch(1.0, 1.0, "contained", True)
        # Allow one construct difference at either end.
        core = seq[1:-1]
        if len(core) >= 4 and core in target_sequence:
            return ChainMatch(len(core) / len(seq), 1.0, "contained", True)
        # Fall through: a short peptide may still be a divergent variant.

    km = kmer_identity(seq, target_sequence)
    if km > identity_threshold:
        return ChainMatch(km, 1.0, "exact k-mer", True)

    return aligned_match(seq, target_sequence)


def matches_target(seq: str, target_sequence: str,
                   identity_threshold: float = IDENTITY_THRESHOLD) -> bool:
    """Boolean form of match_target, kept for callers that only need yes/no."""
    return match_target(seq, target_sequence, identity_threshold).is_match


# Sequence fingerprints of the proteins that bind a receptor from INSIDE the
# cell. These are not off-targets or design partners: a structure solved with
# one of them is showing how a receptor signals, and the interface it reveals
# is on the cytoplasmic face.
# arrestin one is its conserved polar core.
# Short cores rather than full motifs, because the constructs used in these
# structures are engineered: 9XN4's G-alpha carries GAGESGK**N**T, a
# dominant-negative P-loop mutation used to trap the nucleotide-free state,
# and the full wild-type motif misses it entirely.
#
# G-alpha and G-beta only. Arrestin is left out rather than guessed at: this
# project has a rule about writing down sequences from memory, and a motif
# that is nearly right would fail silently in exactly the way the thing it
# is meant to catch does.
INTRACELLULAR_MOTIFS = (
    "GAGESGK",     # G-alpha P-loop
    "DVGGQR",      # G-alpha switch II
    "TLRGH",       # G-beta WD40 repeat
)


def _looks_intracellular(struct: Structure, chains: list[str]) -> bool:
    """Does this bound partner look like a cytoplasmic signalling protein?

    Written after the glucagon receptor entry verified cleanly against a
    receptor-Gs complex and derived a 36-residue binding site from the G
    protein -- the inside of the membrane, where an injected binder cannot
    reach. Every automatic check passed: the chain was the right protein and
    something was certainly bound to it. Nothing asked which side it was on.
    """
    for cid in chains:
        seq = struct.sequence(cid) or ""
        if any(m in seq for m in INTRACELLULAR_MOTIFS):
            return True
    return False


def classify_chains(struct: Structure, target_sequence: str,
                    identity_threshold: float = IDENTITY_THRESHOLD,
                    ) -> tuple[dict[str, str], dict[str, ChainMatch]]:
    """Label every chain as the target, an antibody part, or another protein."""
    roles: dict[str, str] = {}
    matches: dict[str, ChainMatch] = {}
    for cid in struct.chain_ids:
        seq = struct.sequence(cid)
        if len(seq) < MIN_CHAIN_RESIDUES:
            roles[cid] = "too_short"
            continue
        m = match_target(seq, target_sequence, identity_threshold)
        matches[cid] = m
        if m.is_match:
            roles[cid] = "target"
        elif any(mo in seq for mo in HEAVY_MOTIFS):
            roles[cid] = "antibody_heavy"
        elif any(mo in seq for mo in LIGHT_MOTIFS):
            roles[cid] = "antibody_light"
        else:
            roles[cid] = "partner"
    return roles, matches


def partner_groups(roles: dict[str, str]) -> list[tuple[str, list[str]]]:
    """Candidate binders grouped by sequence role. Superseded; see below.

    Kept only because it states the problem clearly. Grouping antibody
    chains by their J-segment motif looked reasonable and was wrong in five
    of 43 catalogue entries: the motifs only occur in chains that include the
    J segment, so a Fab light chain from a non-human germline is labelled a
    generic partner, split off from its own heavy chain, and scored as a
    binder on its own. `contact_groups` decides it geometrically instead.
    """
    groups: list[tuple[str, list[str]]] = []
    antibody = [c for c, r in roles.items()
                if r in ("antibody_heavy", "antibody_light")]
    if antibody:
        groups.append(("antibody", sorted(antibody)))
    for c, r in roles.items():
        if r == "partner":
            groups.append(("protein", [c]))
    return groups


def _patch_distance(a: list[Residue], b: list[Residue]) -> float:
    """Closest approach between two sets of residues, centroid to centroid."""
    if not a or not b:
        return float("inf")
    pa = np.asarray([r.centroid for r in a], dtype=float)
    pb = np.asarray([r.centroid for r in b], dtype=float)
    return float(np.linalg.norm(pa[:, None, :] - pb[None, :, :], axis=-1).min())


def group_binders(patches: dict[str, list[Residue]],
                  pair_contacts: dict[tuple[str, str], int],
                  copies: set[tuple[str, str]] | None = None,
                  linkage: float = PATCH_LINKAGE,
                  min_pair_contacts: int = MIN_PAIR_CONTACTS,
                  ) -> list[list[str]]:
    """Group bound chains into binders: which chains are ONE molecule here.

    Two chains are merged when all of these hold:

      they are not two copies of the same protein, and
      their patches on the target are adjacent, so they cover one site, and
      they pack extensively against each other, so they are one molecule

    Both conditions are needed, and each was learned by getting it wrong.

    Sequence motifs alone split Fabs whose light chain lacked a J-segment
    motif. Chain-to-chain contact alone merged two DIFFERENT antibodies,
    because neighbouring Fabs in a crystal touch. Patch adjacency alone
    merged an antibody bound to one copy of the target with an antibody
    bound to the next copy over, since crystal neighbours land close by.
    Clustering the contacted target residues instead split interfaces that
    are genuinely elongated -- the SARS-CoV-2 RBD/ACE2 interface fell into
    two pieces and lost eight of its twenty residues.

    The two conditions together separate a heavy/light pair, which buries a
    large interface with itself, from a crystal neighbour, which grazes it.
    """
    ids = list(patches)
    parent = {c: c for c in ids}

    def root(c: str) -> str:
        while parent[c] != c:
            parent[c] = parent[parent[c]]
            c = parent[c]
        return c

    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            if root(a) == root(b):
                continue
            if copies and ((a, b) in copies or (b, a) in copies):
                continue
            n = pair_contacts.get((a, b), pair_contacts.get((b, a), 0))
            if (n >= min_pair_contacts
                    and _patch_distance(patches[a], patches[b]) <= linkage):
                parent[root(b)] = root(a)

    merged: dict[str, list[str]] = {}
    for c in ids:
        merged.setdefault(root(c), []).append(c)
    return [sorted(g) for g in merged.values()]


def derive_epitope(pdb_path: str | Path, target_sequence: str, *,
                   cutoff: float = CONTACT_CUTOFF,
                   prefer_chain: str | None = None,
                   prefer_binder: str | None = None) -> EpitopeResult:
    """Find the binding site on `target_sequence`'s chain in this structure.

    `prefer_binder` restricts which bound partner counts as the binder:
    "antibody" where a therapeutic antibody and a natural partner are both
    present and the antibody epitope is the one wanted, or an explicit chain
    id. Left unset, the largest single interface wins.

    Raises EpitopeError when the structure does not contain a target chain
    and a distinct partner in contact with it. Failing loudly matters here:
    a silently invented epitope would poison everything downstream.
    """
    pdb_path = Path(pdb_path)
    struct = parse_structure(pdb_path)
    roles, matches = classify_chains(struct, target_sequence)
    warnings: list[str] = []

    target_chains = [c for c, r in roles.items() if r == "target"]
    if prefer_chain and prefer_chain in target_chains:
        target_chains = [prefer_chain] + [c for c in target_chains
                                          if c != prefer_chain]
    if not target_chains:
        best_cid, best_match = max(
            matches.items(), key=lambda kv: kv[1].identity,
            default=(None, ChainMatch(0.0, 0.0, "aligned", False)))
        raise EpitopeError(
            f"{pdb_path.stem}: no chain matches the target sequence (best "
            f"match was chain {best_cid} at {best_match.identity:.0%} "
            f"identity over {best_match.coverage:.0%} of its length). Either "
            f"the accession and the structure disagree, or this entry does "
            f"not contain the target."
        )

    partner_chains = [c for c, r in roles.items()
                      if r in ("antibody_heavy", "antibody_light", "partner")]
    if not partner_chains:
        raise EpitopeError(
            f"{pdb_path.stem}: the target is present but nothing is bound to "
            f"it, so there is no binding site to derive. A complex with a "
            f"partner chain is required."
        )

    # Evaluate every copy of the target against every distinct binder and keep
    # the largest single interface. A crystal often holds several copies of
    # the complex and sometimes several different partners; the biggest
    # contact surface from ONE binder is the biological interface.
    def is_antibody(chains: list[str]) -> bool:
        return any(roles[c] in ("antibody_heavy", "antibody_light")
                   for c in chains)

    usable = partner_chains
    if prefer_binder == "antibody":
        usable = [c for c in partner_chains
                  if roles[c] in ("antibody_heavy", "antibody_light")]
        if not usable:
            raise EpitopeError(
                f"{pdb_path.stem}: an antibody binder was asked for but no "
                f"bound chain here looks like one.")
    elif prefer_binder:
        if prefer_binder not in partner_chains:
            raise EpitopeError(
                f"{pdb_path.stem}: chain {prefer_binder} was asked for as the "
                f"binder but is not a bound partner here.")
        usable = [prefer_binder]

    candidates = []
    for tc in target_chains:
        t_res = [r for r in struct.chain(tc) if r.is_standard]

        by_chain: dict[str, list[Contact]] = {}
        for pc in usable:
            cs = find_contacts(
                t_res, [r for r in struct.chain(pc) if r.is_standard], cutoff)
            if cs:
                by_chain[pc] = cs
        if not by_chain:
            continue

        patches = {pc: list({c.res_a.key: c.res_a for c in cs}.values())
                   for pc, cs in by_chain.items()}
        near = list(patches)
        pair_contacts = {}
        copies = set()
        for i, a in enumerate(near):
            ra = [r for r in struct.chain(a) if r.is_standard]
            sa = struct.sequence(a)
            for b in near[i + 1:]:
                rb = [r for r in struct.chain(b) if r.is_standard]
                pair_contacts[(a, b)] = len(find_contacts(ra, rb, cutoff))
                if kmer_identity(sa, struct.sequence(b)) >= SAME_MOLECULE_MAX_IDENTITY:
                    copies.add((a, b))
        groups = group_binders(patches, pair_contacts, copies)

        for chains in groups:
            contacts = [c for pc in chains for c in by_chain[pc]]
            epi = sorted({c.res_a.key: c.res_a for c in contacts}.values(),
                         key=lambda r: (r.chain, r.resseq, r.icode))
            para = sorted({c.res_b.key: c.res_b for c in contacts}.values(),
                          key=lambda r: (r.chain, r.resseq, r.icode))
            others = sorted({c for g in groups if g != chains for c in g})
            candidates.append((len(contacts), tc, contacts, epi, para,
                               "antibody" if is_antibody(chains) else "protein",
                               chains, others))

    # Identity outranks interface size. Picking the largest interface alone
    # reads the epitope off whichever chain the binder happens to touch, and
    # in a structure holding several paralogous subunits of one complex that
    # is not necessarily the target. 6Z6V holds all three C1q subunits; the
    # nanobody binds subunit B, and C1QC -- which matches its own chains at
    # 100% and subunit B at 51%, both above the aligned-identity threshold --
    # had its binding site read off subunit B, a different gene product.
    # Nothing downstream could have noticed: the sequence is real, the
    # contacts are real, and the epitope would simply have belonged to
    # another protein.
    #
    # So candidates are considered in bands of identity, best band first, and
    # interface size only decides within a band.
    # The split is SAME_PROTEIN_IDENTITY, the threshold this module already
    # uses for "this chain is the target rather than a relative of it". A
    # first attempt banded identity in fine steps and was too blunt: it
    # demoted chain Z of the PD-1 structure and chain D of the LAG-3 one,
    # which are 90% and 95% identical to their references because
    # crystallographic copies are partly disordered, not because they are
    # different proteins. Both lost their larger, better-ordered interface to
    # a second copy. Copies of the target compete on interface size exactly
    # as before; only a genuine relative is held back.
    def band(chain_id: str) -> int:
        m = matches.get(chain_id)
        return 1 if (m is not None and m.same_protein) else 0

    best = (max(candidates, key=lambda c: (band(c[1]), c[0]))
            if candidates else None)
    if best is not None and candidates:
        bands = {band(c[1]) for c in candidates}
        if len(bands) > 1:
            lost = max(c[0] for c in candidates if band(c[1]) < band(best[1]))
            if lost > best[0]:
                warnings.append(
                    f"a larger interface ({lost} contacts) exists on chain "
                    f"{max((c for c in candidates if band(c[1]) < band(best[1])), key=lambda c: c[0])[1]}, "
                    f"which is below "
                    f"{SAME_PROTEIN_IDENTITY:.0%} identical to the reference "
                    f"sequence. The site was taken "
                    f"from chain {best[1]} instead, which matches the target "
                    f"itself; the other chain is most likely a paralogous "
                    f"subunit in the same complex")
    if best is None:
        if prefer_binder:
            raise EpitopeError(
                f"{pdb_path.stem}: no bound partner matching "
                f"'{prefer_binder}' makes contact with the target.")
        raise EpitopeError(
            f"{pdb_path.stem}: target and partner chains are both present but "
            f"make no contacts within {cutoff} A. They may be separate copies "
            f"in the asymmetric unit rather than a bound complex."
        )

    (_, target_chain, contacts, epi, para, binder_kind, group_chains,
     set_aside) = best

    if set_aside:
        warnings.append(
            f"this structure holds more than one binder; the interface with "
            f"{binder_kind} chains {','.join(group_chains)} was used and "
            f"chains {','.join(set_aside)} were set aside so two separate "
            f"binding sites are not merged into one")

    # The fast paths decide membership cheaply but report on different
    # scales: an exact-k-mer fraction of 0.85 corresponds to roughly 98%
    # sequence identity, so warning on it directly called two entries
    # "orthologs" that differ from the reference by a couple of engineered
    # mutations. Re-measure the chain that actually won, by alignment, so the
    # number shown to a user means one thing.
    match = aligned_match(struct.sequence(target_chain), target_sequence)
    if match and not match.same_protein:
        warnings.append(
            f"chain {target_chain} is {match.identity:.0%} identical to the "
            f"reference sequence, so this structure is a variant or ortholog "
            f"rather than the protein itself; a binder designed against it "
            f"may not transfer")

    # A binding site on the wrong side of the membrane. The automatic checks
    # ask whether the chain is the right protein and whether something is
    # bound to it; neither asks whether the bound thing is somewhere a
    # therapeutic protein could ever reach. The glucagon receptor entry was
    # verified cleanly against 9XN4 and derived a 36-residue "binding site"
    # from the Gs heterotrimer -- the receptor's INTRACELLULAR face, which no
    # injected binder gets near. It was caught by eye, which is not a method,
    # so the sequences that give it away are listed here.
    if _looks_intracellular(struct, group_chains):
        warnings.append(
            f"the bound partner on chains {','.join(group_chains)} looks like "
            f"an intracellular signalling protein (a G protein or arrestin). "
            f"If so this binding site is on the inside of the membrane, where "
            f"an injected protein binder cannot reach it, and the structure "
            f"was solved to study signalling rather than to show a drug "
            f"site. Prefer an entry with an antibody or ligand bound to the "
            f"outside")

    if struct.resolution and struct.resolution > RESOLUTION_WARN:
        warnings.append(
            f"resolution is {struct.resolution:.2f} A; contacts are read from "
            f"side-chain positions, which are not reliably placed beyond "
            f"about {RESOLUTION_WARN:.1f} A, so this binding site is "
            f"approximate")

    if len(epi) < MIN_INTERFACE_RESIDUES:
        warnings.append(
            f"only {len(epi)} residues in the interface, which is small for a "
            f"protein-protein binding site; results will be weak")
    elif len(epi) < SMALL_INTERFACE_RESIDUES:
        warnings.append(
            f"the interface is {len(epi)} residues, which is small: a typical "
            f"antibody epitope covers 15-25. There is little surface here for "
            f"a designed binder to grip, and the scores have correspondingly "
            f"less to discriminate on")

    # Keep only the partner chains that actually touch the target.
    touching = sorted({c.res_b.chain for c in contacts})
    if set(touching) != set(group_chains):
        warnings.append(
            f"partner chains {sorted(set(group_chains) - set(touching))} "
            f"make no contact and were ignored")

    epitope = [
        EpitopeResidue(
            chain=r.chain, resseq=r.resseq, icode=r.icode.strip(),
            resname=r.resname, aa=r.one_letter, label=r.label,
            centroid=tuple(round(float(v), 3) for v in r.centroid),
            min_distance_to_binder=round(
                min(c.min_distance for c in contacts if c.res_a.key == r.key), 3),
        )
        for r in epi
    ]

    epi_seq = "".join(r.aa for r in epitope)
    para_seq = "".join(r.one_letter for r in para)

    n_salt = sum(
        1 for c in contacts
        if {c.res_a.one_letter, c.res_b.one_letter} & set("DE")
        and {c.res_a.one_letter, c.res_b.one_letter} & set("KR")
        and c.min_distance <= 4.0
    )
    n_hydrophobic = sum(
        1 for c in contacts
        if c.res_a.one_letter in "AVILMFWCY" and c.res_b.one_letter in "AVILMFWCY"
    )

    return EpitopeResult(
        pdb_id=struct.pdb_id,
        title=struct.title,
        resolution=struct.resolution,
        target_chain=target_chain,
        binder_chains=touching,
        chain_roles=roles,
        epitope=epitope,
        paratope_sequence=para_seq,
        epitope_sequence=epi_seq,
        descriptors={
            "n_residue_contacts": len(contacts),
            "n_salt_bridges": n_salt,
            "n_hydrophobic_contacts": n_hydrophobic,
            "epitope_radius_of_gyration": round(radius_of_gyration(epi), 3),
            "n_epitope_residues": len(epitope),
            "n_paratope_residues": len(para),
        },
        contacts=[
            {
                "target": c.res_a.label,
                "binder": c.res_b.label,
                "min_distance": round(c.min_distance, 3),
                "atom_pairs": c.n_atom_pairs,
            }
            for c in contacts
        ],
        binder_kind=binder_kind,
        target_identity=round(match.identity if match else 0.0, 4),
        target_coverage=round(match.coverage if match else 0.0, 4),
        identity_method=match.method if match else "unknown",
        warnings=warnings,
    )


def describe_surface(epitope_sequence: str) -> dict:
    """Plain-language character of a binding site, used to brief the designer.

    This is what makes the pipeline target-agnostic: instead of knowing that
    CLEC12A has three arginines, it reads whatever the surface happens to be
    and states the requirement that follows from it.
    """
    from .aminoacids import AROMATIC, CHARGED, HYDROPHOBIC, POLAR_UNCHARGED

    if not epitope_sequence:
        return {}
    n = len(epitope_sequence)
    cations = sum(1 for a in epitope_sequence if a in "KR")
    anions = sum(1 for a in epitope_sequence if a in "DE")
    arom = sum(1 for a in epitope_sequence if a in AROMATIC)
    hyd = sum(1 for a in epitope_sequence if a in HYDROPHOBIC)
    pol = sum(1 for a in epitope_sequence if a in POLAR_UNCHARGED)

    net = cations - anions
    if net >= 2:
        want = "acidic"
        why = (f"{cations} basic residues dominate the surface, so a binder "
               f"needs acidic groups to pair them and aromatics to stack "
               f"against them")
    elif net <= -2:
        want = "basic"
        why = (f"{anions} acidic residues dominate, so a binder needs basic "
               f"groups to pair them")
    elif hyd / n > 0.5:
        want = "hydrophobic"
        why = (f"{hyd} of {n} residues are hydrophobic, so burial of apolar "
               f"surface will drive binding more than electrostatics")
    else:
        want = "polar and balanced"
        why = (f"the surface is mixed ({pol} polar, {cations + anions} "
               f"charged), so hydrogen bonding and shape matter more than "
               f"net charge")

    return {
        "length": n,
        "cations": cations,
        "anions": anions,
        "net_formal_charge": net,
        "aromatic": arom,
        "hydrophobic": hyd,
        "polar": pol,
        "required_character": want,
        "rationale": why,
    }
