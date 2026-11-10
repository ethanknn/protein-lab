"""Amino acid reference data.

Every table here comes from a published source, named in the comment above it.
Where a value is a heuristic of our own, it says so. Nothing in this file is
invented silently -- downstream scores are only as honest as these numbers.
"""

from __future__ import annotations

AA20 = "ACDEFGHIKLMNPQRSTVWY"

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    # Common modified / alternate residues mapped to their parent.
    "MSE": "M", "SEC": "C", "PYL": "K", "HSD": "H", "HSE": "H", "HSP": "H",
    "CSO": "C", "CME": "C", "PTR": "Y", "SEP": "S", "TPO": "T",
}

ONE_TO_THREE = {v: k for k, v in list(THREE_TO_ONE.items())[:20]}

FULL_NAME = {
    "A": "alanine", "C": "cysteine", "D": "aspartate", "E": "glutamate",
    "F": "phenylalanine", "G": "glycine", "H": "histidine", "I": "isoleucine",
    "K": "lysine", "L": "leucine", "M": "methionine", "N": "asparagine",
    "P": "proline", "Q": "glutamine", "R": "arginine", "S": "serine",
    "T": "threonine", "V": "valine", "W": "tryptophan", "Y": "tyrosine",
}

# Average residue mass in Da (monomer mass minus one water).
# Source: standard IUPAC average isotopic masses.
RESIDUE_MASS = {
    "A": 71.0788, "C": 103.1388, "D": 115.0886, "E": 129.1155, "F": 147.1766,
    "G": 57.0519, "H": 137.1411, "I": 113.1594, "K": 128.1741, "L": 113.1594,
    "M": 131.1926, "N": 114.1038, "P": 97.1167, "Q": 128.1307, "R": 156.1875,
    "S": 87.0782, "T": 101.1051, "V": 99.1326, "W": 186.2132, "Y": 163.1760,
}
WATER_MASS = 18.01528

# Kyte & Doolittle (1982) hydropathy index. J Mol Biol 157:105-132.
KD_HYDROPATHY = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}

# Eisenberg et al. (1984) consensus hydrophobicity scale. J Mol Biol 179:125-142.
EISENBERG = {
    "A": 0.62, "C": 0.29, "D": -0.90, "E": -0.74, "F": 1.19,
    "G": 0.48, "H": -0.40, "I": 1.38, "K": -1.50, "L": 1.06,
    "M": 0.64, "N": -0.78, "P": 0.12, "Q": -0.85, "R": -2.53,
    "S": -0.18, "T": -0.05, "V": 1.08, "W": 0.81, "Y": 0.26,
}

# Side-chain pKa values used for the Henderson-Hasselbalch charge model.
# Source: Nelson & Cox, Lehninger Principles of Biochemistry (standard values).
SIDECHAIN_PKA = {"D": 3.65, "E": 4.25, "H": 6.00, "C": 8.18, "Y": 10.07,
                 "K": 10.53, "R": 12.48}
POSITIVE_SIDECHAINS = {"H", "K", "R"}
NEGATIVE_SIDECHAINS = {"D", "E", "C", "Y"}
PKA_N_TERMINUS = 9.69
PKA_C_TERMINUS = 2.34

# Formal charge at physiological pH, used for quick electrostatic bookkeeping.
FORMAL_CHARGE = {aa: 0 for aa in AA20}
FORMAL_CHARGE.update({"D": -1, "E": -1, "K": 1, "R": 1, "H": 0.1})

# Chou & Fasman (1978) conformational parameters. Adv Enzymol 47:45-148.
HELIX_PROPENSITY = {
    "E": 1.51, "M": 1.45, "A": 1.42, "L": 1.21, "K": 1.16, "F": 1.13,
    "Q": 1.11, "W": 1.08, "I": 1.08, "V": 1.06, "D": 1.01, "H": 1.00,
    "R": 0.98, "T": 0.83, "S": 0.77, "C": 0.70, "Y": 0.69, "N": 0.67,
    "P": 0.57, "G": 0.57,
}
SHEET_PROPENSITY = {
    "V": 1.70, "I": 1.60, "Y": 1.47, "F": 1.38, "W": 1.37, "L": 1.30,
    "C": 1.19, "T": 1.19, "Q": 1.10, "M": 1.05, "R": 0.93, "N": 0.89,
    "H": 0.87, "A": 0.83, "S": 0.75, "G": 0.75, "K": 0.74, "P": 0.55,
    "D": 0.54, "E": 0.37,
}

# Residue side-chain volumes in cubic angstroms.
# Source: Zamyatnin (1972) Prog Biophys Mol Biol 24:107-123.
RESIDUE_VOLUME = {
    "A": 88.6, "C": 108.5, "D": 111.1, "E": 138.4, "F": 189.9,
    "G": 60.1, "H": 153.2, "I": 166.7, "K": 168.6, "L": 166.7,
    "M": 162.9, "N": 114.1, "P": 112.7, "Q": 143.8, "R": 173.4,
    "S": 89.0, "T": 116.1, "V": 140.0, "W": 227.8, "Y": 193.6,
}

# Maximum accessible surface area per residue (angstrom^2), theoretical.
# Source: Tien et al. (2013) PLoS ONE 8:e80635, Table 1 "theoretical".
MAX_ASA = {
    "A": 129.0, "C": 167.0, "D": 193.0, "E": 223.0, "F": 240.0,
    "G": 104.0, "H": 224.0, "I": 197.0, "K": 236.0, "L": 201.0,
    "M": 224.0, "N": 195.0, "P": 159.0, "Q": 225.0, "R": 274.0,
    "S": 155.0, "T": 172.0, "V": 174.0, "W": 285.0, "Y": 263.0,
}

# Molar extinction coefficients at 280 nm (M^-1 cm^-1).
# Source: Pace et al. (1995) Protein Sci 4:2411-2423.
EXT_COEFF_280 = {"W": 5500.0, "Y": 1490.0, "C": 125.0}  # C = cystine contribution

# Van der Waals radii in angstroms, by element. Bondi (1964) J Phys Chem 68:441.
VDW_RADII = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "H": 1.20,
             "P": 1.80, "SE": 1.90, "F": 1.47, "CL": 1.75}
DEFAULT_VDW = 1.70

# --- Groupings used by the mutation operators and the scorers ---------------

HYDROPHOBIC = set("AVILMFWCY")
AROMATIC = set("FWY")
ALIPHATIC = set("AVILM")
POLAR_UNCHARGED = set("STNQ")
CHARGED = set("DEKRH")
SMALL = set("AGSCP")
HBOND_DONORS = set("NQSTYWKRH")
HBOND_ACCEPTORS = set("NQSTYDE")

# Residues that most often drive beta-aggregation in short peptides. Used by the
# aggregation-propensity heuristic, which is clearly labelled as a heuristic.
AGGREGATION_PRONE = set("VILFWYM")

# Liabilities that would complicate real manufacture of a protein therapeutic.
# These are standard developability flags used in antibody engineering.
LIABILITY_MOTIFS = {
    "NG": "deamidation hotspot (Asn-Gly)",
    "NS": "deamidation risk (Asn-Ser)",
    "NT": "deamidation risk (Asn-Thr)",
    "DG": "isomerisation hotspot (Asp-Gly)",
    "DP": "acid-labile cleavage (Asp-Pro)",
    "DS": "isomerisation risk (Asp-Ser)",
    "M": "oxidation-prone methionine",
    "W": "oxidation-prone tryptophan",
}
NGLYC_NOTE = "N-X-S/T sequon: N-linked glycosylation site"


def is_valid_sequence(seq: str) -> bool:
    return bool(seq) and all(c in AA20 for c in seq)


def clean_sequence(seq: str) -> str:
    """Uppercase and strip anything that is not one of the standard 20."""
    return "".join(c for c in seq.upper() if c in AA20)


def molecular_weight(seq: str) -> float:
    """Average molecular weight of a linear peptide in Da."""
    if not seq:
        return 0.0
    return sum(RESIDUE_MASS[a] for a in seq) + WATER_MASS


def gravy(seq: str) -> float:
    """Grand average of hydropathy (Kyte-Doolittle mean)."""
    if not seq:
        return 0.0
    return sum(KD_HYDROPATHY[a] for a in seq) / len(seq)


def aliphatic_index(seq: str) -> float:
    """Ikai (1980) aliphatic index. J Biochem 88:1895-1898.

    Higher values correlate with greater thermostability.
    """
    if not seq:
        return 0.0
    n = len(seq)
    x_a = seq.count("A") / n
    x_v = seq.count("V") / n
    x_il = (seq.count("I") + seq.count("L")) / n
    return 100.0 * (x_a + 2.9 * x_v + 3.9 * x_il)


def aromaticity(seq: str) -> float:
    """Fraction of Phe + Trp + Tyr (Lobry & Gautier 1994)."""
    if not seq:
        return 0.0
    return sum(seq.count(a) for a in "FWY") / len(seq)


def net_charge(seq: str, ph: float = 7.4) -> float:
    """Net charge from the Henderson-Hasselbalch equation.

    Includes both termini. pH 7.4 is blood plasma, which is the condition a
    therapeutic candidate would actually encounter.
    """
    if not seq:
        return 0.0
    charge = 0.0
    # N-terminal amine: positive when protonated.
    charge += 1.0 / (1.0 + 10 ** (ph - PKA_N_TERMINUS))
    # C-terminal carboxyl: negative when deprotonated.
    charge -= 1.0 / (1.0 + 10 ** (PKA_C_TERMINUS - ph))
    for aa in seq:
        pka = SIDECHAIN_PKA.get(aa)
        if pka is None:
            continue
        if aa in POSITIVE_SIDECHAINS:
            charge += 1.0 / (1.0 + 10 ** (ph - pka))
        else:
            charge -= 1.0 / (1.0 + 10 ** (pka - ph))
    return charge


def isoelectric_point(seq: str) -> float:
    """Theoretical pI found by bisection on net_charge."""
    if not seq:
        return 7.0
    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2.0
        if net_charge(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2.0, 2)


def extinction_coefficient(seq: str, reduced: bool = True) -> float:
    """Molar extinction coefficient at 280 nm (M^-1 cm^-1).

    Matters practically: a candidate with no Trp or Tyr is very hard to
    quantify by absorbance in a real lab.
    """
    eps = seq.count("W") * EXT_COEFF_280["W"] + seq.count("Y") * EXT_COEFF_280["Y"]
    if not reduced:
        eps += (seq.count("C") // 2) * EXT_COEFF_280["C"]
    return eps


def hydrophobic_moment(seq: str, angle_deg: float = 100.0) -> float:
    """Eisenberg hydrophobic moment, normalised per residue.

    An alpha helix advances 100 degrees per residue, so a high moment means
    hydrophobic and polar faces are segregated -- the signature of an
    amphipathic helix, which is what a good helical binder face looks like.
    """
    import math

    if not seq:
        return 0.0
    rad = math.radians(angle_deg)
    sx = sum(EISENBERG[a] * math.cos(i * rad) for i, a in enumerate(seq))
    sy = sum(EISENBERG[a] * math.sin(i * rad) for i, a in enumerate(seq))
    return math.hypot(sx, sy) / len(seq)


def composition(seq: str) -> dict[str, float]:
    """Fractional composition over all 20 residues."""
    if not seq:
        return {a: 0.0 for a in AA20}
    n = len(seq)
    return {a: seq.count(a) / n for a in AA20}


def shannon_entropy(seq: str) -> float:
    """Sequence entropy in bits. Low values flag low-complexity repeats,
    which tend to be sticky, disordered, and hard to express."""
    import math

    if not seq:
        return 0.0
    n = len(seq)
    ent = 0.0
    for a in set(seq):
        p = seq.count(a) / n
        ent -= p * math.log2(p)
    return ent


def max_homopolymer_run(seq: str) -> int:
    """Longest run of a single repeated residue."""
    if not seq:
        return 0
    best = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, run)
    return best
