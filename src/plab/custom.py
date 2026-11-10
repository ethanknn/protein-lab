"""Targets a person builds themselves, from any protein in the databases.

The catalogue is 71 targets somebody curated. This module is the same
pipeline pointed at whatever protein someone types in: resolve the name,
find a structure that holds it with something bound, read the binding site
off that structure, derive the cross-reactivity panel, measure what the
selectivity score can and cannot reach, and hand back a target the design
loop can run against.

The one thing it cannot do is invent a binding site. This project reads
epitopes off structures of a complex, so a protein solved only on its own,
or not solved at all, yields no target -- and that is reported as the fact it
is rather than as an error. Roughly nine in ten human proteins have no
structure with a bound partner, which is a statement about structural biology
and not about the protein's importance.

A target built here is NOT a catalogue target. Catalogue entries were
checked by a person, carry a written rationale and a clinical caveat, and
were kept or rejected deliberately. These carry none of that, and the
interface says so. The verification that CAN be automated -- the chain really
is this protein, something really is bound to it, the interface is big enough
to be a binding site -- is the same code the catalogue build runs, and its
warnings are passed through untouched.
"""

from __future__ import annotations

import json
import re
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

from .epitope import EpitopeError, derive_epitope, describe_surface
from .lookup import fetch_accession, rank_structures, resolve_symbol
from .panels import build_panel_for, save_sequences
from .targets import OffTarget, TargetSpec, _derived_concern

ROOT = Path(__file__).resolve().parent.parent.parent
STRUCT_DIR = ROOT / "data" / "structures"
SEQ_DIR = ROOT / "data" / "sequences"
USER_AGENT = "protein-lab/0.5 (educational in-silico design sandbox)"

# An accession, not a gene symbol: P01375, Q5QGZ9, A0A0B4J2F0.
ACCESSION = re.compile(r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|"
                       r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|"
                       r"^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]$")


def _get(url: str, timeout: int = 90) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8")


# --- Step one: what protein does this person mean? --------------------------


def search(query: str) -> dict:
    """Resolve a gene symbol or accession, and list usable structures.

    Returns every reviewed human protein matching the symbol rather than the
    first, because the first is not reliably the right one: `SELP` returns
    selenoprotein P ahead of P-selectin, and `MET` returns a transcription
    modulator alongside the kinase. Choosing between them is the person's
    job and they cannot do it if only one is shown.
    """
    q = query.strip()
    if not q:
        return {"query": q, "proteins": []}

    if ACCESSION.match(q.upper()):
        rec = fetch_accession(q.upper())
        hits = [rec] if rec else []
    else:
        hits = resolve_symbol(q)

    out = []
    for h in hits[:5]:
        out.append({
            "accession": h["accession"],
            "name": h["name"],
            "genes": h.get("genes") or [],
            "length": h.get("length"),
            "organism": h.get("organism"),
        })
    return {"query": q, "proteins": out}


def structures(accession: str) -> list[dict]:
    """Entries holding this protein with a partner bound, best first."""
    return rank_structures(accession, limit=40, top=8)


# --- Step two: fetch and read the binding site ------------------------------


def fetch_sequence(accession: str) -> str | None:
    dest = SEQ_DIR / f"{accession}.fasta"
    if dest.exists() and dest.stat().st_size > 0:
        text = dest.read_text(encoding="utf-8")
    else:
        try:
            text = _get(f"https://rest.uniprot.org/uniprotkb/{accession}.fasta")
        except Exception:
            return None
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(text, encoding="utf-8")
    lines = text.splitlines()
    return "".join(lines[1:]) if lines else None


def fetch_structure(pdb_id: str) -> Path | None:
    """Legacy PDB format where it exists, mmCIF otherwise.

    The fallback is what makes recent entries reachable at all: many have no
    .pdb file, and a person typing in a gene will land on them constantly.
    """
    for ext in ("pdb", "cif"):
        dest = STRUCT_DIR / f"{pdb_id}.{ext}"
        if dest.exists() and dest.stat().st_size > 0:
            return dest
        try:
            text = _get(f"https://files.rcsb.org/download/{pdb_id}.{ext}")
        except urllib.error.HTTPError:
            continue
        except Exception:
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(text, encoding="utf-8")
        return dest
    return None


@dataclass
class CustomTarget:
    """A target somebody built, plus everything known about how well."""
    key: str
    spec: TargetSpec
    record: dict
    panel: dict = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


def prepare(accession: str, pdb_id: str, *, disease: str = "",
            with_panel: bool = True, progress=None) -> CustomTarget:
    """Build a runnable target from an accession and a structure.

    Raises EpitopeError with the same wording the catalogue build uses when
    the structure cannot yield a binding site. Those messages were written to
    be read by somebody deciding what to do next, so they are not reworded
    here.
    """
    def say(stage: str, frac: float) -> None:
        if progress:
            progress(stage, frac)

    accession = accession.upper().strip()
    pdb_id = pdb_id.upper().strip()
    key = f"custom_{accession.lower()}_{pdb_id.lower()}"

    say(f"fetching the reference sequence for {accession}", 0.08)
    sequence = fetch_sequence(accession)
    if not sequence:
        raise EpitopeError(
            f"UniProt has no sequence for {accession}. Check the accession, "
            f"or search by gene symbol instead.")

    protein = fetch_accession(accession) or {}
    gene = (protein.get("genes") or [accession])[0]

    say(f"downloading structure {pdb_id}", 0.18)
    path = fetch_structure(pdb_id)
    if path is None:
        raise EpitopeError(
            f"{pdb_id} could not be downloaded in either PDB or mmCIF "
            f"format. The entry may be withdrawn, or the id may be wrong.")

    say(f"reading the binding site off {pdb_id}", 0.32)
    result = derive_epitope(str(path), sequence)
    if not result.ok:
        raise EpitopeError(
            f"the largest interface in {pdb_id} is only "
            f"{len(result.epitope)} residues, which is too small to treat as "
            f"a binding site. Try another structure of this protein.")

    record = {
        "key": key, "gene": gene,
        "protein_name": protein.get("name") or gene,
        "disease": disease or "User-defined target",
        "area": "Built by you",
        "uniprot": accession,
        "pdb_complex": pdb_id, "pdb_apo": None,
        "known_binder": result.binder_kind,
        "rationale": "",
        "clinical_note": "",
        "sequence_length": len(sequence),
        "resolution": result.resolution,
        "structure_title": result.title,
        "target_chain": result.target_chain,
        "binder_chains": result.binder_chains,
        "binder_kind": result.binder_kind,
        "epitope_sequence": result.epitope_sequence,
        "paratope_sequence": result.paratope_sequence,
        "n_epitope_residues": len(result.epitope),
        "descriptors": result.descriptors,
        "surface": describe_surface(result.epitope_sequence),
        "target_identity": result.target_identity,
        "target_coverage": result.target_coverage,
        "identity_method": result.identity_method,
        "warnings": result.warnings,
        "off_targets": [],
    }

    spec = TargetSpec(
        key=key, gene=gene,
        protein_name=record["protein_name"],
        disease=record["disease"], area=record["area"],
        uniprot=accession,
        rationale="", clinical_note="",
        apo_pdb=None, complex_pdb=pdb_id,
        reference_binder_name=f"the {result.binder_kind} bound in {pdb_id}",
        binder_kind=result.binder_kind,
        epitope=result.epitope,
        paratope_sequence=result.paratope_sequence,
        full_sequence=sequence,
        descriptors=result.descriptors,
        surface=record["surface"],
        chain_roles=result.chain_roles,
        target_chain=result.target_chain,
        binder_chains=result.binder_chains,
        contacts=result.contacts,
        source_resolution=result.resolution,
        structure_title=result.title,
        off_target_accessions=[],
        panel_state="not_built",
        panel_note="",
    )

    panel: dict = {}
    if with_panel:
        say("asking UniProt which human proteins share a domain family", 0.46)
        try:
            panel = build_panel_for(spec)
        except Exception as exc:
            spec.panel_state = "panel_failed"
            spec.panel_note = (
                f"the domain-family lookup did not complete ({type(exc).__name__}). "
                f"Selectivity is not assessed for this target.")
            panel = {}
        if panel:
            say("saving the panel members' sequences", 0.72)
            try:
                save_sequences(panel)
            except Exception:
                pass
            spec.panel_state = panel.get("state", "none")
            spec.panel_note = panel.get("note", "")
            spec.panel_families = panel.get("families", [])
            spec.panel_considered = panel.get("n_considered", 0)
            spec.off_target_accessions = [m["uniprot"]
                                          for m in panel.get("members", [])]
            record["off_targets"] = list(spec.off_target_accessions)

    if spec.off_target_accessions:
        say("measuring what selectivity is reachable against that panel", 0.84)
        _measure_ceilings(spec, panel)

    notes = list(result.warnings)
    return CustomTarget(key=key, spec=spec, record=record, panel=panel,
                        notes=notes)


def _measure_ceilings(spec: TargetSpec, panel: dict) -> None:
    """The best reachable margin and its error bar, for this target's panel.

    Without these the scorer falls back to an absolute cutoff, which is the
    one that let VEGF-A report 100% selective off a margin of 0.029. A custom
    target deserves the same guard as a catalogue one, so it is measured
    here even though the jackknife costs a few seconds.
    """
    from .design.scaffolds import REGISTRY, get_scaffold
    from .scoring.selectivity import (achievable_margin_range,
                                      build_paralog_epitopes, ceiling_noise)
    from .targets import load_off_targets

    surrogates = build_paralog_epitopes(spec, load_off_targets(spec))
    if not surrogates:
        return
    epi = spec.epitope_sequence
    ceilings, noises = {}, {}
    for name in REGISTRY:
        n = len(get_scaffold(name).variable_positions)
        if not n or str(n) in ceilings:
            continue
        ceilings[str(n)] = round(achievable_margin_range(epi, surrogates, n)[1], 5)
        noises[str(n)] = round(ceiling_noise(epi, surrogates, n), 5)
    spec.panel_ceiling = ceilings
    spec.panel_ceiling_noise = noises
    panel["ceiling"] = ceilings
    panel["ceiling_noise"] = noises


def off_targets_for(spec: TargetSpec, panel: dict) -> list[OffTarget]:
    """Panel members as OffTarget records, without touching panels.json.

    load_off_targets reads the stored panel file, which a custom target is
    deliberately not written into: the curated catalogue is a reviewed
    artefact and a session's experiments do not belong in it.
    """
    derived = {m["uniprot"]: m for m in panel.get("members", [])}
    out = []
    for acc in spec.off_target_accessions:
        path = SEQ_DIR / f"{acc}.fasta"
        if not path.exists():
            continue
        lines = path.read_text(encoding="utf-8").splitlines()
        if not lines:
            continue
        header, seq = lines[0], "".join(lines[1:])
        gene = derived.get(acc, {}).get("gene") or acc
        out.append(OffTarget(
            key=f"{gene}_{acc}", gene=gene, label=header[1:110],
            uniprot=acc, sequence=seq,
            concern=_derived_concern(derived.get(acc)),
            structural=True,
        ))
    return out
