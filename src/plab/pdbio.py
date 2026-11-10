"""A small, strict PDB reader plus interface analysis.

Why hand-rolled: the pipeline needs exactly three things from a structure --
per-chain sequences, atom coordinates, and inter-chain contacts. A focused
parser that fails loudly on malformed input is easier to trust than a large
dependency, and the PDB fixed-column format is stable.

Column positions follow the official PDB format specification v3.30.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .aminoacids import DEFAULT_VDW, THREE_TO_ONE, VDW_RADII

# Solvent and common crystallisation additives we never want in the model.
IGNORED_HETERO = {
    "HOH", "DOD", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG", "ACT", "CL",
    "NA", "K", "MG", "ZN", "CA", "MN", "NI", "CO", "CD", "IOD", "BR",
    "TRS", "MES", "EPE", "DMS", "FMT", "NO3", "ACY", "PGE", "1PE", "P6G",
}


@dataclass(slots=True)
class Atom:
    serial: int
    name: str
    resname: str
    chain: str
    resseq: int
    icode: str
    x: float
    y: float
    z: float
    occupancy: float
    bfactor: float
    element: str
    hetero: bool

    @property
    def coord(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=np.float64)

    @property
    def vdw(self) -> float:
        return VDW_RADII.get(self.element.upper(), DEFAULT_VDW)

    @property
    def is_backbone(self) -> bool:
        return self.name in ("N", "CA", "C", "O")

    @property
    def is_sidechain(self) -> bool:
        return not self.is_backbone and self.element != "H"


@dataclass(slots=True)
class Residue:
    resname: str
    chain: str
    resseq: int
    icode: str
    atoms: list[Atom] = field(default_factory=list)

    @property
    def key(self) -> tuple[str, int, str]:
        return (self.chain, self.resseq, self.icode)

    @property
    def label(self) -> str:
        """Human-readable id, e.g. 'A:Lys244'."""
        one = self.one_letter
        return f"{self.chain}:{self.resname.capitalize()}{self.resseq}{self.icode.strip()}" \
            if one == "X" else f"{self.chain}:{one}{self.resseq}{self.icode.strip()}"

    @property
    def one_letter(self) -> str:
        return THREE_TO_ONE.get(self.resname.upper(), "X")

    @property
    def is_standard(self) -> bool:
        return self.resname.upper() in THREE_TO_ONE

    def atom(self, name: str) -> Atom | None:
        for a in self.atoms:
            if a.name == name:
                return a
        return None

    @property
    def ca(self) -> Atom | None:
        return self.atom("CA")

    def coords(self, sidechain_only: bool = False) -> np.ndarray:
        sel = [a for a in self.atoms if (a.is_sidechain if sidechain_only else True)]
        if not sel:
            sel = self.atoms
        return np.array([[a.x, a.y, a.z] for a in sel], dtype=np.float64)

    @property
    def centroid(self) -> np.ndarray:
        c = self.coords()
        return c.mean(axis=0) if len(c) else np.zeros(3)


@dataclass
class Structure:
    pdb_id: str
    residues: list[Residue] = field(default_factory=list)
    title: str = ""
    resolution: float | None = None

    @property
    def chain_ids(self) -> list[str]:
        seen: list[str] = []
        for r in self.residues:
            if r.chain not in seen:
                seen.append(r.chain)
        return seen

    def chain(self, chain_id: str) -> list[Residue]:
        return [r for r in self.residues if r.chain == chain_id]

    def sequence(self, chain_id: str) -> str:
        return "".join(r.one_letter for r in self.chain(chain_id) if r.is_standard)

    def atoms(self, chain_id: str | None = None) -> list[Atom]:
        rs = self.residues if chain_id is None else self.chain(chain_id)
        return [a for r in rs for a in r.atoms]

    def residue_by_number(self, chain_id: str, resseq: int) -> Residue | None:
        for r in self.chain(chain_id):
            if r.resseq == resseq:
                return r
        return None

    def summary(self) -> str:
        parts = [f"{self.pdb_id}: {self.title[:70]}"]
        if self.resolution:
            parts.append(f"  resolution: {self.resolution:.2f} A")
        for c in self.chain_ids:
            seq = self.sequence(c)
            parts.append(f"  chain {c}: {len(self.chain(c))} residues, {len(seq)} standard")
        return "\n".join(parts)


def parse_structure(path: str | Path, **kw) -> Structure:
    """Parse a structure file, PDB or mmCIF, chosen by extension.

    mmCIF is not an optional extra any more. The legacy PDB format cannot
    express an entry with more than 62 chains or 99,999 atoms, so the PDB
    stopped generating .pdb files for entries that exceed it -- and for a
    growing number of ordinary ones. Two of the targets added in the second
    research pass, the IL-31 receptor with nemolizumab (9EZA) and the
    glucagon receptor (9XN4), have no .pdb file at all: the download returns
    404 and the candidate was rejected for "could not fetch", which reads
    like a dead structure rather than a reader that cannot open it.

    That matters more for the design-your-own tab than for the catalogue,
    because a person typing in a gene will land on recent entries constantly
    and has no way to know that the reason is a file format.
    """
    path = Path(path)
    if path.suffix.lower() in (".cif", ".mmcif"):
        return parse_mmcif(path, **kw)
    return parse_pdb(path, **kw)


def parse_mmcif(path: str | Path, *, first_model_only: bool = True,
                keep_hetero: bool = False) -> Structure:
    """Parse the _atom_site loop of an mmCIF file.

    Deliberately narrow: this reads the same three things the PDB reader
    reads -- per-chain sequences, coordinates and the identifiers needed to
    group atoms into residues -- and ignores the rest of the dictionary. A
    full mmCIF parser is a large piece of software, and everything else this
    pipeline needs about an entry it gets from the RCSB API instead.

    Uses the auth_* identifiers rather than label_*, because those are the
    chain and residue numbers the rest of this project, the RCSB entity
    listing and every published paper all use.
    """
    path = Path(path)
    struct = Structure(pdb_id=path.stem.upper())
    by_key: dict[tuple[str, int, str], Residue] = {}
    seen_atom: set[tuple[str, int, str, str]] = set()

    fields: list[str] = []
    in_loop = False

    def col(name: str) -> int:
        try:
            return fields.index(name)
        except ValueError:
            return -1

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_struct.title"):
                struct.title = _cif_value(s, "_struct.title")
                continue
            if s.startswith("_refine.ls_d_res_high") or \
                    s.startswith("_em_3d_reconstruction.resolution"):
                v = _cif_value(s, s.split()[0])
                try:
                    struct.resolution = float(v)
                except (TypeError, ValueError):
                    pass
                continue
            if s == "loop_":
                fields, in_loop = [], False
                continue
            if s.startswith("_atom_site."):
                fields.append(s.split()[0][len("_atom_site."):])
                in_loop = True
                continue
            if not in_loop:
                continue
            if not s or s.startswith("#") or s.startswith("_") or s == "loop_":
                in_loop = False
                continue

            parts = _cif_split(s)
            if len(parts) < len(fields):
                continue
            idx = _CIF_COLS
            if idx.get("cache") is not fields:
                idx = {
                    "cache": fields,
                    "group": col("group_PDB"), "name": col("label_atom_id"),
                    "alt": col("label_alt_id"), "resname": col("label_comp_id"),
                    "chain": col("auth_asym_id"), "resseq": col("auth_seq_id"),
                    "icode": col("pdbx_PDB_ins_code"),
                    "x": col("Cartn_x"), "y": col("Cartn_y"), "z": col("Cartn_z"),
                    "occ": col("occupancy"), "b": col("B_iso_or_equiv"),
                    "element": col("type_symbol"),
                    "model": col("pdbx_PDB_model_num"),
                    "label_chain": col("label_asym_id"),
                    "label_seq": col("label_seq_id"),
                }
                _CIF_COLS.clear()
                _CIF_COLS.update(idx)

            if idx["model"] >= 0 and first_model_only and \
                    parts[idx["model"]] not in ("1", ".", "?"):
                continue

            hetero = parts[idx["group"]] == "HETATM" if idx["group"] >= 0 else False
            resname = parts[idx["resname"]].strip().upper()
            if hetero and resname in IGNORED_HETERO:
                continue
            if hetero and not keep_hetero and resname not in THREE_TO_ONE:
                continue

            chain_i = idx["chain"] if idx["chain"] >= 0 else idx["label_chain"]
            seq_i = idx["resseq"] if idx["resseq"] >= 0 else idx["label_seq"]
            chain = (parts[chain_i] or "_").strip() or "_"
            name = parts[idx["name"]].strip().strip('"')
            try:
                resseq = int(parts[seq_i])
                x = float(parts[idx["x"]])
                y = float(parts[idx["y"]])
                z = float(parts[idx["z"]])
            except (ValueError, IndexError):
                continue
            icode_raw = parts[idx["icode"]] if idx["icode"] >= 0 else "?"
            icode = " " if icode_raw in (".", "?") else icode_raw

            dedupe = (chain, resseq, icode, name)
            if dedupe in seen_atom:
                continue  # alternate location; keep the first
            seen_atom.add(dedupe)

            element = (parts[idx["element"]].strip()
                       if idx["element"] >= 0 else name[:1])
            atom = Atom(
                serial=0, name=name, resname=resname, chain=chain,
                resseq=resseq, icode=icode, x=x, y=y, z=z,
                occupancy=_cif_float(parts, idx["occ"], 1.0),
                bfactor=_cif_float(parts, idx["b"], 0.0),
                element=element or name[:1], hetero=hetero,
            )
            key = (chain, resseq, icode)
            res = by_key.get(key)
            if res is None:
                res = Residue(resname=resname, chain=chain, resseq=resseq,
                              icode=icode)
                by_key[key] = res
                struct.residues.append(res)
            res.atoms.append(atom)

    if not struct.residues:
        raise ValueError(f"No usable atoms parsed from {path}")
    return struct


def to_pdb_text(struct: Structure, chains: list[str] | None = None) -> str:
    """Write a parsed structure back out as PDB-format ATOM records.

    The viewer speaks PDB and the chain-subsetting code is written against
    it, so rather than teach both about mmCIF, an mmCIF-only entry is parsed
    and re-emitted here. Only the fields the viewer uses are written.

    mmCIF allows chain identifiers longer than the single column PDB gives
    them. Those are remapped to single characters, and because the epitope
    was derived from the same parse, a remapped entry's chain labels in the
    viewer will not match the ones named elsewhere on the page. Single
    character ids, which is nearly every entry, pass through untouched.
    """
    wanted = set(chains) if chains else None
    ids = sorted({r.chain for r in struct.residues})
    pool = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    remap: dict[str, str] = {}
    used: set[str] = {c for c in ids if len(c) == 1}
    for c in ids:
        if len(c) == 1:
            remap[c] = c
            continue
        pick = next((p for p in pool if p not in used), "X")
        used.add(pick)
        remap[c] = pick

    lines, serial = [], 1
    if struct.title:
        lines.append(f"TITLE     {struct.title[:70]}")
    if struct.resolution:
        lines.append(f"REMARK   2 RESOLUTION. {struct.resolution:.2f} ANGSTROMS.")
    for res in struct.residues:
        if wanted is not None and res.chain not in wanted:
            continue
        for a in res.atoms:
            rec = "HETATM" if a.hetero else "ATOM  "
            name = a.name if len(a.name) >= 4 else f" {a.name:<3}"
            lines.append(
                f"{rec}{serial:>5} {name:<4}{'':1}{a.resname:>3} "
                f"{remap.get(a.chain, 'X')}{a.resseq:>4}{a.icode or ' '}   "
                f"{a.x:>8.3f}{a.y:>8.3f}{a.z:>8.3f}"
                f"{a.occupancy:>6.2f}{a.bfactor:>6.2f}          "
                f"{a.element:>2}")
            serial += 1
    lines.append("END")
    return "\n".join(lines) + "\n"


_CIF_COLS: dict = {}


def _cif_float(parts: list[str], i: int, default: float) -> float:
    if i < 0 or i >= len(parts):
        return default
    try:
        return float(parts[i])
    except ValueError:
        return default


def _cif_value(line: str, key: str) -> str:
    v = line[len(key):].strip()
    if v.startswith("'") and v.endswith("'") and len(v) > 1:
        v = v[1:-1]
    elif v.startswith('"') and v.endswith('"') and len(v) > 1:
        v = v[1:-1]
    return v


def _cif_split(line: str) -> list[str]:
    """Split a loop row, honouring the single and double quoting mmCIF uses.

    Atom names like 'CA' and descriptions with spaces are quoted, so a plain
    split() silently shifts every column after them.
    """
    out, cur, quote = [], [], ""
    for ch in line:
        if quote:
            if ch == quote:
                quote = ""
                out.append("".join(cur))
                cur = []
            else:
                cur.append(ch)
        elif ch in "'\"":
            quote = ch
        elif ch.isspace():
            if cur:
                out.append("".join(cur))
                cur = []
        else:
            cur.append(ch)
    if cur:
        out.append("".join(cur))
    return out


def parse_pdb(path: str | Path, *, first_model_only: bool = True,
              keep_hetero: bool = False) -> Structure:
    """Parse a PDB file into a Structure.

    Takes the first altloc encountered per atom name, drops solvent, and by
    default reads only the first MODEL (NMR ensembles would otherwise stack
    dozens of copies on top of each other).
    """
    path = Path(path)
    struct = Structure(pdb_id=path.stem.upper())
    by_key: dict[tuple[str, int, str], Residue] = {}
    seen_atom: set[tuple[str, int, str, str]] = set()
    in_later_model = False

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            rec = line[:6]

            if rec == "TITLE ":
                struct.title = (struct.title + " " + line[10:80].strip()).strip()
                continue
            if rec == "REMARK" and " RESOLUTION. " in line:
                tok = line.split("RESOLUTION.")[-1].split()
                if tok and tok[0].replace(".", "", 1).isdigit():
                    struct.resolution = float(tok[0])
                continue
            if rec == "MODEL ":
                num = line[10:14].strip()
                in_later_model = first_model_only and num not in ("", "1")
                continue
            if rec == "ENDMDL":
                continue
            if rec not in ("ATOM  ", "HETATM"):
                continue
            if in_later_model:
                continue

            hetero = rec == "HETATM"
            resname = line[17:20].strip().upper()
            if hetero and resname in IGNORED_HETERO:
                continue
            if hetero and not keep_hetero and resname not in THREE_TO_ONE:
                continue

            name = line[12:16].strip()
            chain = line[21].strip() or "_"
            try:
                resseq = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue  # malformed coordinate line: skip rather than guess
            icode = line[26]

            dedupe = (chain, resseq, icode, name)
            if dedupe in seen_atom:
                continue  # alternate location; keep the first
            seen_atom.add(dedupe)

            element = line[76:78].strip() or name[:1]
            atom = Atom(
                serial=int(line[6:11]) if line[6:11].strip().isdigit() else 0,
                name=name, resname=resname, chain=chain, resseq=resseq,
                icode=icode, x=x, y=y, z=z,
                occupancy=float(line[54:60]) if line[54:60].strip() else 1.0,
                bfactor=float(line[60:66]) if line[60:66].strip() else 0.0,
                element=element, hetero=hetero,
            )

            key = (chain, resseq, icode)
            res = by_key.get(key)
            if res is None:
                res = Residue(resname=resname, chain=chain, resseq=resseq, icode=icode)
                by_key[key] = res
                struct.residues.append(res)
            res.atoms.append(atom)

    if not struct.residues:
        raise ValueError(f"No usable atoms parsed from {path}")
    return struct


# --- Interface analysis -----------------------------------------------------


@dataclass(slots=True)
class Contact:
    res_a: Residue
    res_b: Residue
    min_distance: float
    n_atom_pairs: int

    def describe(self) -> str:
        return (f"{self.res_a.label} -- {self.res_b.label}  "
                f"{self.min_distance:.2f} A  ({self.n_atom_pairs} pairs)")


def _atom_matrix(residues: list[Residue]) -> tuple[np.ndarray, np.ndarray]:
    """Return (coords, residue_index) for all heavy atoms."""
    coords, owner = [], []
    for i, r in enumerate(residues):
        for a in r.atoms:
            if a.element.upper() == "H":
                continue
            coords.append((a.x, a.y, a.z))
            owner.append(i)
    return np.array(coords, dtype=np.float64), np.array(owner, dtype=np.int32)


def find_contacts(residues_a: list[Residue], residues_b: list[Residue],
                  cutoff: float = 4.5) -> list[Contact]:
    """All residue pairs across two groups with any heavy-atom pair under cutoff.

    4.5 A is the conventional heavy-atom contact distance for protein-protein
    interfaces -- it captures hydrogen bonds, salt bridges and van der Waals
    packing without pulling in second-shell residues.
    """
    if not residues_a or not residues_b:
        return []
    ca, oa = _atom_matrix(residues_a)
    cb, ob = _atom_matrix(residues_b)
    if len(ca) == 0 or len(cb) == 0:
        return []

    # Pairwise distances. Interface sizes here are small enough that the dense
    # matrix is cheaper than building a spatial index.
    d = np.linalg.norm(ca[:, None, :] - cb[None, :, :], axis=-1)
    hits = np.argwhere(d <= cutoff)
    if len(hits) == 0:
        return []

    agg: dict[tuple[int, int], list[float]] = {}
    for ia, ib in hits:
        agg.setdefault((int(oa[ia]), int(ob[ib])), []).append(float(d[ia, ib]))

    out = [
        Contact(res_a=residues_a[i], res_b=residues_b[j],
                min_distance=min(ds), n_atom_pairs=len(ds))
        for (i, j), ds in agg.items()
    ]
    out.sort(key=lambda c: c.min_distance)
    return out


def interface_residues(residues_a: list[Residue], residues_b: list[Residue],
                       cutoff: float = 4.5) -> tuple[list[Residue], list[Residue]]:
    """The residues on each side that form the interface, in sequence order."""
    contacts = find_contacts(residues_a, residues_b, cutoff)
    seen_a: dict[tuple, Residue] = {}
    seen_b: dict[tuple, Residue] = {}
    for c in contacts:
        seen_a.setdefault(c.res_a.key, c.res_a)
        seen_b.setdefault(c.res_b.key, c.res_b)
    key = lambda r: (r.chain, r.resseq, r.icode)
    return sorted(seen_a.values(), key=key), sorted(seen_b.values(), key=key)


def center_of_mass(residues: list[Residue]) -> np.ndarray:
    coords = [[a.x, a.y, a.z] for r in residues for a in r.atoms]
    return np.array(coords).mean(axis=0) if coords else np.zeros(3)


def radius_of_gyration(residues: list[Residue]) -> float:
    coords = np.array([[a.x, a.y, a.z] for r in residues for a in r.atoms])
    if len(coords) == 0:
        return 0.0
    return float(np.sqrt(((coords - coords.mean(axis=0)) ** 2).sum(axis=1).mean()))
