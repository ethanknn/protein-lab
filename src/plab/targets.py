"""Target definitions, loaded from the verified catalogue.

A TargetSpec is the design brief: the surface we aim at, the reference binder
we benchmark against, and the paralogs we want to avoid. Nothing here is
hand-written per target -- it is assembled from data that
scripts/build_catalog.py has already fetched and verified, so adding a disease
means adding a catalogue candidate, not editing this file.

The scorers read the epitope's chemistry out of a TargetSpec rather than
knowing anything about a particular protein. That is what lets one pipeline
serve acute myeloid leukaemia, COVID-19 and rheumatoid arthritis without
changing a line of scoring code.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

from .aminoacids import AROMATIC, CHARGED, HYDROPHOBIC, net_charge
from .catalog import CATALOG_DIR, load_verified

ROOT = Path(__file__).resolve().parent.parent.parent
DATA = ROOT / "data"
SEQ_DIRS = (DATA / "sequences", DATA / "targets")
STRUCT_DIR = DATA / "structures"


def _find_fasta(acc: str) -> Path | None:
    for d in SEQ_DIRS:
        p = d / f"{acc}.fasta"
        if p.exists() and p.stat().st_size > 0:
            return p
    return None


def _read_fasta(path: Path) -> tuple[str, str]:
    lines = [l.strip() for l in path.read_text(encoding="utf-8").splitlines()
             if l.strip()]
    header = lines[0] if lines and lines[0].startswith(">") else ""
    return header, "".join(l for l in lines if not l.startswith(">"))


def _gene_from_header(header: str, fallback: str) -> str:
    for tok in header.split():
        if tok.startswith("GN="):
            return tok[3:]
    return fallback


def organism_from_header(header: str) -> str:
    """The OS= field of a UniProt FASTA header, or "".

    UniProt writes the organism between OS= and the next tagged field, and it
    contains spaces ("Zaire ebolavirus (strain Mayinga-76)"), so it cannot be
    pulled out by splitting on whitespace the way GN= can.
    """
    if "OS=" not in header:
        return ""
    rest = header.split("OS=", 1)[1]
    for tag in (" OX=", " GN=", " PE=", " SV="):
        if tag in rest:
            rest = rest.split(tag, 1)[0]
    return rest.strip()


@dataclass(frozen=True)
class EpitopeResidue:
    chain: str
    resseq: int
    resname: str
    aa: str
    label: str
    centroid: tuple[float, float, float]
    min_distance_to_binder: float


@dataclass
class OffTarget:
    key: str
    gene: str
    label: str
    uniprot: str
    sequence: str
    concern: str
    structural: bool = True


@dataclass
class TargetSpec:
    """Everything the designer and scorers need to know about one target."""

    key: str
    gene: str
    protein_name: str
    disease: str
    area: str
    uniprot: str
    rationale: str
    clinical_note: str
    apo_pdb: str | None
    complex_pdb: str
    reference_binder_name: str
    binder_kind: str = "protein"

    epitope: list[EpitopeResidue] = field(default_factory=list)
    paratope_sequence: str = ""
    full_sequence: str = ""
    descriptors: dict = field(default_factory=dict)
    surface: dict = field(default_factory=dict)
    chain_roles: dict = field(default_factory=dict)
    target_chain: str = ""
    binder_chains: list[str] = field(default_factory=list)
    contacts: list[dict] = field(default_factory=list)
    source_resolution: float | None = None
    structure_title: str = ""
    # Source organism, from the UniProt FASTA header's OS= field. Carried
    # because "which healthy human cells carry this protein" is not a
    # question that can be asked about the Ebola glycoprotein, and a tool
    # that answered "unknown" there would be describing its own ignorance
    # rather than the target. See lineage.is_human.
    organism: str = ""
    off_target_accessions: list[str] = field(default_factory=list)
    # How the cross-reactivity panel was arrived at, and why it is the size
    # it is. "derived", "curated", "no_human_relatives", "unmappable",
    # "no_family", "no_mappable_relatives". These are NOT interchangeable:
    # an empty panel because a viral protein has no human relatives is a
    # fact about the target, while an empty panel because nobody built one
    # is a fact about this software, and they must never render alike.
    panel_state: str = "none"
    panel_note: str = ""
    panel_families: list[str] = field(default_factory=list)
    panel_considered: int = 0
    # Best reachable margin and its error bar, per face length, measured
    # by scripts/measure_ceilings.py. The ceiling decides what 100% means
    # and the error bar decides whether the ceiling means anything, so
    # both travel with the target rather than being re-derived per run.
    # Surrogate surfaces already mapped by scripts/build_panels.py, and
    # the epitope they were mapped onto. Re-deriving them means a
    # Smith-Waterman alignment per panel member per Evaluator, which is
    # the single most expensive thing this package does and produces the
    # identical answer every time.
    panel_surfaces: dict[str, str] = field(default_factory=dict)
    panel_epitope: str = ""
    panel_ceiling: dict[str, float] = field(default_factory=dict)
    panel_ceiling_noise: dict[str, float] = field(default_factory=dict)

    # --- derived properties of the surface we design against --------------

    @cached_property
    def epitope_sequence(self) -> str:
        return "".join(r.aa for r in self.epitope)

    @cached_property
    def epitope_net_charge(self) -> float:
        """Charge of the epitope side chains alone (termini excluded).

        Termini are excluded because an epitope is a patch on a folded
        protein, not a free peptide.
        """
        from .aminoacids import POSITIVE_SIDECHAINS, SIDECHAIN_PKA
        ph = 7.4
        q = 0.0
        for r in self.epitope:
            pka = SIDECHAIN_PKA.get(r.aa)
            if pka is None:
                continue
            if r.aa in POSITIVE_SIDECHAINS:
                q += 1.0 / (1.0 + 10 ** (ph - pka))
            else:
                q -= 1.0 / (1.0 + 10 ** (pka - ph))
        return q

    @cached_property
    def epitope_aromatic_fraction(self) -> float:
        if not self.epitope:
            return 0.0
        return sum(1 for r in self.epitope if r.aa in AROMATIC) / len(self.epitope)

    @cached_property
    def epitope_hydrophobic_fraction(self) -> float:
        if not self.epitope:
            return 0.0
        return sum(1 for r in self.epitope if r.aa in HYDROPHOBIC) / len(self.epitope)

    @cached_property
    def epitope_charged_fraction(self) -> float:
        if not self.epitope:
            return 0.0
        return sum(1 for r in self.epitope if r.aa in CHARGED) / len(self.epitope)

    @cached_property
    def epitope_radius(self) -> float:
        return float(self.descriptors.get("epitope_radius_of_gyration", 0.0))

    @cached_property
    def paratope_charge(self) -> float:
        return net_charge(self.paratope_sequence)

    @cached_property
    def paratope_aromatic_fraction(self) -> float:
        if not self.paratope_sequence:
            return 0.0
        return sum(1 for a in self.paratope_sequence if a in AROMATIC) / len(
            self.paratope_sequence)

    @cached_property
    def arginine_positions(self) -> list[str]:
        """Arginines in the epitope. For an arginine-rich surface these are
        the dominant recognition feature, engaged by acidic and aromatic
        partners. Empty for surfaces where charge is not the driver."""
        return [r.label for r in self.epitope if r.aa == "R"]

    @cached_property
    def cation_positions(self) -> list[str]:
        return [r.label for r in self.epitope if r.aa in "KR"]

    def design_brief(self) -> str:
        """A plain-language statement of what a good binder must do."""
        s = self.surface or {}
        want = s.get("required_character", "complementary")
        why = s.get("rationale", "")
        cations = self.cation_positions
        return (
            f"Target: {self.gene} ({self.protein_name}) in {self.disease}.\n"
            f"Surface: {len(self.epitope)} residues, sequence "
            f"{self.epitope_sequence}, radius of gyration "
            f"{self.epitope_radius:.1f} A.\n"
            f"Epitope side-chain charge at pH 7.4: "
            f"{self.epitope_net_charge:+.2f}"
            + (f" ({len(cations)} basic residues: {', '.join(cations[:6])})"
               if cations else "")
            + f".\n"
            f"Implied requirement: a compact, {want} binding face"
            + (f" -- {why}" if why else "") + ".\n"
            f"Reference binder ({self.reference_binder_name}) achieves this "
            f"with contact residues {self.paratope_sequence} "
            f"(charge {self.paratope_charge:+.2f}, "
            f"{self.paratope_aromatic_fraction:.0%} aromatic).\n"
            f"Experimental source: PDB {self.complex_pdb}"
            + (f" at {self.source_resolution:.1f} A." if self.source_resolution
               else ".")
        )


# --- Off-target panels ------------------------------------------------------

# Why each off-target matters. Keyed by accession; a target without an entry
# gets a generic note. Most catalogue targets have no curated panel at all,
# and the selectivity scorer reports that honestly rather than inventing one.
OFF_TARGET_CONCERNS = {
    "Q2HXU8": ("The immediate neighbour of CLEC12A in the NK gene complex and "
               "its closest paralog. If a binder cross-reacts with anything, "
               "this is the most likely thing."),
    "Q6UXN8": ("C-type lectin on dendritic cells. Cross-reactivity would hit "
               "antigen-presenting cells."),
    "Q9BXN2": ("Dectin-1, central to antifungal immunity. Engaging it would "
               "compromise host defence in already immunocompromised patients."),
    "Q8NC01": "C-type lectin in the same gene complex; cross-reactivity undesirable.",
    "Q12918": ("CD161 on NK and T cells -- the effector cells an immunotherapy "
               "depends on."),
    "P20138": ("CD33 is a cautionary comparison, not a structural relative: a "
               "Siglec, a different fold family. It sits on AML blasts but "
               "also on healthy haematopoietic stem cells, which is why "
               "CD33-directed agents carry lasting marrow toxicity."),
}

# Accessions whose fold family differs from the target, so mapping an epitope
# onto them by alignment would be structurally meaningless.
NON_STRUCTURAL = {"P20138"}


def _derived_concern(member: dict | None) -> str:
    """Say what is actually known about a derived panel member.

    "Related protein; cross-reactivity undesirable" is true of every entry
    and therefore tells a reader nothing. The numbers behind the selection
    do: how similar this protein is at the binding site, and how much of
    that site could be compared at all.
    """
    if not member:
        return "Related protein; cross-reactivity undesirable."
    ident = member.get("window_identity")
    cover = member.get("window_coverage")
    n = member.get("n_mapped")
    parts = ["Shares a domain family with the target"]
    if ident is not None:
        parts.append(f"and is {ident:.0%} identical across the binding site")
    if cover is not None and n:
        parts.append(f"({n} epitope positions mapped, {cover:.0%} of the "
                     f"surrounding window aligned)")
    return " ".join(parts) + "."


def load_off_targets(target: TargetSpec | None = None) -> list[OffTarget]:
    """Proteins a candidate should NOT engage.

    Cross-reactivity is assessed on sequence composition at aligned surface
    positions. That is a genuine triage signal and how real campaigns
    pre-screen paralogs, but it is not a substitute for a cross-reactivity
    panel measured on cells.
    """
    accessions = list(target.off_target_accessions) if target else []
    derived = {}
    if target:
        derived = {m["uniprot"]: m
                   for m in load_panels().get(target.key, {}).get("members", [])}
    out: list[OffTarget] = []
    for acc in accessions:
        path = _find_fasta(acc)
        if path is None:
            continue
        header, seq = _read_fasta(path)
        if not seq:
            continue
        gene = _gene_from_header(header, acc)
        out.append(OffTarget(
            key=f"{gene}_{acc}", gene=gene,
            label=header[1:110] if header else acc,
            uniprot=acc, sequence=seq,
            concern=OFF_TARGET_CONCERNS.get(acc) or _derived_concern(derived.get(acc)),
            structural=acc not in NON_STRUCTURAL,
        ))
    return out


# --- Loading a target -------------------------------------------------------

DEFAULT_TARGET = "clec12a_aml"


def available_targets() -> dict[str, dict]:
    return load_verified().get("targets", {})


def load_target(key: str = DEFAULT_TARGET) -> TargetSpec:
    """Build the design brief for any verified catalogue target."""
    catalog = available_targets()

    if key not in catalog:
        # Accept a gene symbol or disease word as a convenience.
        lowered = key.lower()
        for k, rec in catalog.items():
            if lowered in (rec["gene"].lower(), k.lower()) \
                    or lowered in rec["disease"].lower():
                key = k
                break
        else:
            raise KeyError(
                f"Unknown target {key!r}. Verified targets: "
                f"{sorted(catalog)}. Run scripts/build_catalog.py to add more."
            )

    rec = catalog[key]
    epi_path = CATALOG_DIR / rec["epitope_file"]
    if not epi_path.exists():
        raise FileNotFoundError(
            f"{epi_path} missing. Re-run scripts/build_catalog.py --only {key}.")
    epi_data = json.loads(epi_path.read_text(encoding="utf-8"))

    seq_path = _find_fasta(rec["uniprot"])
    full_sequence = _read_fasta(seq_path)[1] if seq_path else ""

    epitope = [
        EpitopeResidue(
            chain=r["chain"], resseq=r["resseq"], resname=r["resname"],
            aa=r["aa"], label=r["label"],
            centroid=tuple(r["centroid"]),
            min_distance_to_binder=r["min_distance_to_binder"],
        )
        for r in epi_data["epitope"]
    ]

    return TargetSpec(
        key=key,
        gene=rec["gene"],
        protein_name=rec["protein_name"],
        disease=rec["disease"],
        area=rec["area"],
        uniprot=rec["uniprot"],
        rationale=rec["rationale"],
        clinical_note=rec["clinical_note"],
        apo_pdb=rec.get("pdb_apo"),
        complex_pdb=rec["pdb_complex"],
        reference_binder_name=rec["known_binder"],
        binder_kind=rec.get("binder_kind", "protein"),
        epitope=epitope,
        paratope_sequence=rec.get("paratope_sequence", ""),
        full_sequence=full_sequence,
        descriptors=rec.get("descriptors", {}),
        surface=rec.get("surface", {}),
        chain_roles=epi_data.get("chain_roles", {}),
        target_chain=rec.get("target_chain", ""),
        binder_chains=rec.get("binder_chains", []),
        contacts=epi_data.get("contacts", []),
        source_resolution=rec.get("resolution"),
        structure_title=rec.get("structure_title", ""),
        organism=organism_from_header(rec.get("uniprot_header", "")),
        **_panel_fields(key, rec.get("off_targets", [])),
    )


PANELS_PATH = CATALOG_DIR / "panels.json"


_PANEL_CACHE: dict[str, dict] = {}


def load_panels() -> dict:
    """Cached because load_target is called in tight loops.

    The flat-component check alone evaluates 150 sequences per target, and
    re-reading and re-parsing the panel file on each of those turns a cheap
    lookup into the slowest thing in the scorer.
    """
    if not PANELS_PATH.exists():
        return {}
    stamp = str(PANELS_PATH.stat().st_mtime_ns)
    if _PANEL_CACHE.get("stamp") != stamp:
        _PANEL_CACHE.clear()
        _PANEL_CACHE["stamp"] = stamp
        _PANEL_CACHE["panels"] = json.loads(
            PANELS_PATH.read_text(encoding="utf-8")).get("panels", {})
    return _PANEL_CACHE["panels"]


def _panel_fields(key: str, curated: list[str]) -> dict:
    """Merge the hand-written panel with the one derived from UniProt.

    Curated members come first and keep their written-out concern, because
    they carry clinical knowledge a sequence database does not have: CD33 is
    in CLEC12A's panel not because it is a close relative -- it is a Siglec,
    a different fold -- but because CD33-directed drugs damage the marrow,
    which is the entire reason CLEC12A was picked over it.
    """
    panel = load_panels().get(key, {})
    derived = [m["uniprot"] for m in panel.get("members", [])]
    merged = list(dict.fromkeys(list(curated) + derived))

    if curated and derived:
        state = "curated+derived"
    elif curated:
        state = "curated"
    else:
        state = panel.get("state", "none")

    return {
        "off_target_accessions": merged,
        "panel_state": state,
        "panel_note": panel.get("note", ""),
        "panel_families": panel.get("families", []),
        "panel_considered": panel.get("n_considered", 0),
        "panel_surfaces": {m["gene"]: m["surrogate_surface"]
                           for m in panel.get("members", [])
                           if m.get("surrogate_surface")},
        "panel_epitope": panel.get("target_epitope", ""),
        "panel_ceiling": panel.get("ceiling", {}),
        "panel_ceiling_noise": panel.get("ceiling_noise", {}),
    }
