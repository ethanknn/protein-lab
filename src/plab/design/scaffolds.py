"""Binder scaffolds.

Design decision: we do not generate sequences from nothing. Real campaigns
start from a scaffold with a known fold and vary only the surface positions
that form the binding face. That keeps the fold trustworthy, shrinks the
search space from 20^58 to 20^13, and means a hit could actually be built.

The primary scaffold here is the affibody Z domain -- a 58-residue
three-helix bundle derived from the B domain of staphylococcal protein A.
It is used because:
  * the fold is small, cysteine-free, and folds without chaperones
  * affibody molecules are a genuinely validated scaffold class in humans
  * the 13 randomised surface positions below are the published library
    positions, on helices 1 and 2, which is where affibody binders bind

Framework residues are held fixed. Only the designated positions vary. That
is a deliberate constraint, not a limitation of the code.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# --- The Z domain -----------------------------------------------------------

# Z domain of staphylococcal protein A, 58 residues.
Z_DOMAIN_SEQ = "VDNKFNKEQQNAFYEILHLPNLNEEQRNAFIQSLKDDPSQSANLLAEAKKLNDAQAPK"

# Secondary structure: H helix, L loop/coil.
#                      1........10........20........30........40........50....58
Z_DOMAIN_SS = "LLLLLLLHHHHHHHHHHHLLLLLLHHHHHHHHHHHHLLLLLHHHHHHHHHHHHHHLLL"

# Published affibody library positions (1-based), on helices 1 and 2.
Z_DOMAIN_VARIABLE_1BASED = [9, 10, 11, 13, 14, 17, 18, 24, 25, 27, 28, 32, 35]


@dataclass
class Scaffold:
    """A fold plus the positions we are allowed to change."""

    name: str
    description: str
    framework: str                      # wild-type / parent sequence
    secondary_structure: str            # per-residue H / E / L / T
    variable_positions: list[int]       # 0-based indices we may mutate
    interface_positions: list[int]      # 0-based indices presented to target
    size_band: tuple[float, float] = (25.0, 65.0)
    pdb_template: str | None = None     # PDB id for 3D display
    pdb_chain: str | None = None
    provenance: str = ""
    # Positions we must never touch (prolines, glycines holding the fold).
    frozen_positions: list[int] = field(default_factory=list)

    def __post_init__(self) -> None:
        n = len(self.framework)
        if len(self.secondary_structure) != n:
            raise ValueError(
                f"{self.name}: secondary structure is "
                f"{len(self.secondary_structure)} long but framework is {n}"
            )
        bad = [p for p in self.variable_positions if not 0 <= p < n]
        if bad:
            raise ValueError(f"{self.name}: variable positions out of range: {bad}")

    @property
    def length(self) -> int:
        return len(self.framework)

    @property
    def n_variable(self) -> int:
        return len(self.variable_positions)

    @property
    def search_space(self) -> float:
        return 20.0 ** self.n_variable

    def face(self, sequence: str) -> str:
        """The residues this scaffold presents to the target."""
        return "".join(sequence[i] for i in self.interface_positions
                       if 0 <= i < len(sequence))

    def apply(self, base: str, assignment: dict[int, str]) -> str:
        """Return `base` with the given position -> residue substitutions."""
        chars = list(base)
        for pos, aa in assignment.items():
            if 0 <= pos < len(chars) and pos not in self.frozen_positions:
                chars[pos] = aa
        return "".join(chars)

    def variable_view(self, sequence: str) -> str:
        """The designed positions only, in order -- the 'designed face'."""
        return "".join(sequence[i] for i in self.variable_positions)

    def annotate(self, sequence: str) -> str:
        """Framework in lowercase, designed positions in uppercase."""
        var = set(self.variable_positions)
        return "".join(c.upper() if i in var else c.lower()
                       for i, c in enumerate(sequence))

    def describe(self) -> str:
        return (
            f"{self.name}: {self.description}\n"
            f"  length {self.length} aa, {self.n_variable} designed positions "
            f"(search space 20^{self.n_variable} = {self.search_space:.2e})\n"
            f"  designed positions (1-based): "
            f"{[p + 1 for p in self.variable_positions]}\n"
            f"  provenance: {self.provenance}"
        )


def z_domain_affibody() -> Scaffold:
    variable = [p - 1 for p in Z_DOMAIN_VARIABLE_1BASED]
    return Scaffold(
        name="affibody_Z",
        description=(
            "Three-helix bundle from the B domain of staphylococcal protein A. "
            "Cysteine-free, folds spontaneously, ~6.5 kDa, thermostable."
        ),
        framework=Z_DOMAIN_SEQ,
        secondary_structure=Z_DOMAIN_SS,
        variable_positions=variable,
        interface_positions=variable,  # the library face IS the binding face
        size_band=(50.0, 66.0),
        pdb_template="1Q2N",
        pdb_chain="A",
        provenance=(
            "Z domain sequence and the 13 randomised helix-1/helix-2 positions "
            "follow the published affibody library design (Nord et al., "
            "Nat Biotechnol 1997; Lofblom et al., FEBS Lett 2010)."
        ),
        frozen_positions=[19, 37, 56],  # P20, P38, P57 in 1-based numbering
    )


# --- A smaller, disulfide-constrained alternative ---------------------------

# A generic two-stranded beta hairpin with a stabilising disulfide. Shorter
# than the Z domain, so it buries less surface, but it diffuses better and is
# cheap to synthesise chemically. The framework here is a generic designed
# hairpin, not a specific published molecule -- flagged as such.
HAIRPIN_SEQ = "CTWEGNKLTCSGNKLTVNGKTLKGC"
HAIRPIN_SS = "LEEEEEELLLLEEEEEELLLEEEEL"
HAIRPIN_VARIABLE_1BASED = [3, 5, 7, 9, 13, 15, 17, 19, 22, 24]


def beta_hairpin() -> Scaffold:
    variable = [p - 1 for p in HAIRPIN_VARIABLE_1BASED]
    cys = [i for i, a in enumerate(HAIRPIN_SEQ) if a == "C"]
    return Scaffold(
        name="beta_hairpin_SS",
        description=(
            "Short disulfide-constrained beta hairpin. Small footprint, easy to "
            "synthesise, but buries less surface so affinity is harder to reach."
        ),
        framework=HAIRPIN_SEQ,
        secondary_structure=HAIRPIN_SS,
        variable_positions=variable,
        interface_positions=variable,
        size_band=(20.0, 35.0),
        pdb_template=None,
        provenance=(
            "Generic designed hairpin topology, not a specific published "
            "molecule. Included as a size contrast to the Z domain."
        ),
        frozen_positions=cys,  # the disulfide holds the fold together
    )


# --- The reference: the antibody that actually binds this epitope ----------


def antibody_vh_reference(sequence: str) -> Scaffold:
    """Wraps the 50C1 heavy-chain variable domain as a scaffold.

    Used only so the reference binder can be scored by the exact same
    machinery as the candidates. Its size band is widened to the antibody
    range, because penalising a Fab for being Fab-sized would make the
    comparison meaningless rather than honest.
    """
    n = len(sequence)
    return Scaffold(
        name="antibody_VH_50C1",
        description="Heavy-chain variable domain of the anti-CLEC12A antibody 50C1.",
        framework=sequence,
        secondary_structure="L" * n,  # no structural claim made
        variable_positions=[],
        interface_positions=[],       # filled in by the caller from real contacts
        size_band=(90.0, 140.0),
        pdb_template="8W9J",
        pdb_chain="H",
        provenance="Extracted from PDB 8W9J chain H.",
    )


REGISTRY = {
    "affibody_Z": z_domain_affibody,
    "beta_hairpin_SS": beta_hairpin,
}


def get_scaffold(name: str = "affibody_Z") -> Scaffold:
    if name not in REGISTRY:
        raise KeyError(f"Unknown scaffold {name!r}. Available: {list(REGISTRY)}")
    return REGISTRY[name]()


def list_scaffolds() -> list[Scaffold]:
    return [factory() for factory in REGISTRY.values()]
