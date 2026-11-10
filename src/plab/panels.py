"""Deriving a cross-reactivity panel for any target, catalogue or not.

Lifted out of scripts/build_panels.py when the web application needed to
build a panel for a protein a person had just typed in. The script is now a
CLI over this module, so the panel behind a custom target is built by exactly
the same code as the panel behind a catalogue one -- including every
threshold, every rejection and every honest "no panel, and here is why".

The reasoning behind the method is in that script's docstring; the short
version is that membership of a domain family says two proteins are
relatives, and only alignment at the epitope WINDOW says whether they
resemble each other at the patch a binder would touch.
"""

from __future__ import annotations

import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

from .catalog import CATALOG_DIR
from .epitope import kmer_identity
from .scoring.selectivity import epitope_window, map_epitope_onto

ROOT = Path(__file__).resolve().parent.parent.parent

SEQ_DIR = ROOT / "data" / "sequences"
PANELS_PATH = CATALOG_DIR / "panels.json"
USER_AGENT = "protein-lab/0.4 (educational in-silico design sandbox)"

MAX_PANEL = 8
# Above this, the "paralog" is an isoform or a near-duplicate entry of the
# target itself. Scoring selectivity against yourself is not a measurement.
MAX_SELF_IDENTITY = 0.95
# A family with more members than this is too generic to be informative --
# usually a domain shared across an entire superfamily.
MAX_FAMILY_SIZE = 400
# A relative whose alignment covers only a sliver of the epitope window is
# matching some other part of the protein, not this surface.
MIN_WINDOW_COVERAGE = 0.5
# Identity at the epitope window below which a "relative" is not one.
#
# This is the softest threshold in the project and it is worth saying why.
# There is no clean gap to put it in: measured across this catalogue, the
# genuine C-type lectin paralogs of CLEC12A -- CLEC9A, Dectin-1, CD161, the
# same proteins a human picked by hand -- sit at 26-33% identity at this
# surface, while the immunoglobulin-fold proteins that are NOT relatives of
# PD-L1 sit at 19-22%. Those ranges nearly touch. 25% is the conventional
# edge of the alignment twilight zone and it separates these two cases, but
# it is a judgement, not a discovered boundary, and a panel built at 26%
# identity is a far weaker instrument than one built at 70%. That is why
# every member carries its identity into the interface rather than being
# presented as simply "a relative".
MIN_WINDOW_IDENTITY = 0.25
# A "relative" whose mapped surface is identical to the target's own is not a
# relative, it is a duplicate database entry for the same protein. UniProt
# carries several of these for immunoglobulins: P0DOX4 is a full IgE heavy
# chain and P01854 is its constant region, so whole-sequence comparison puts
# them at 64% while the binding site is character-for-character the same.
# Left in a panel it is fatal rather than merely noisy: selectivity is the
# margin over the CLOSEST relative, so a member scoring identically to the
# target pins that margin at zero for every candidate forever.
MAX_SURFACE_IDENTITY = 0.95
# One distant relative is not a panel. Where the only member sits near the
# identity floor, the controls in scripts/diagnose_panels.py show the metric
# becoming unreliable -- for the malaria target it inverted, rating a face
# designed for the wrong protein above the real antibody.
MIN_MEMBERS_IF_DISTANT = 2
MODERATE_IDENTITY = 0.35
FAMILY_CACHE = ROOT / "data" / "families"


def _get(url: str, timeout: int = 120) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8")


def pfam_families(accession: str) -> list[str]:
    url = (f"https://rest.uniprot.org/uniprotkb/{accession}.json"
           f"?fields=xref_pfam")
    try:
        data = json.loads(_get(url))
    except Exception as exc:
        print(f"      family lookup failed: {exc}")
        return []
    return [x["id"] for x in data.get("uniProtKBCrossReferences", [])
            if x.get("database") == "Pfam"]


def parse_fasta(text: str) -> list[tuple[str, str, str]]:
    """-> [(accession, gene, sequence)]"""
    out, acc, gene, buf = [], None, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if acc:
                out.append((acc, gene or acc, "".join(buf)))
            parts = line[1:].split("|")
            acc = parts[1] if len(parts) > 2 else line[1:].split()[0]
            gene = acc
            for tok in line.split():
                if tok.startswith("GN="):
                    gene = tok[3:]
                    break
            buf = []
        else:
            buf.append(line.strip())
    if acc:
        out.append((acc, gene or acc, "".join(buf)))
    return out


def family_members(pfam: str) -> list[tuple[str, str, str]]:
    """Every reviewed human protein carrying this domain, with sequences."""
    cached = FAMILY_CACHE / f"{pfam}.fasta"
    if cached.exists() and cached.stat().st_size > 0:
        return parse_fasta(cached.read_text(encoding="utf-8"))
    query = f"(xref:pfam-{pfam}) AND (organism_id:9606) AND (reviewed:true)"
    url = ("https://rest.uniprot.org/uniprotkb/stream?query="
           + urllib.parse.quote(query) + "&format=fasta")
    try:
        text = _get(url)
    except Exception as exc:
        print(f"      {pfam}: fetch failed: {exc}")
        return []
    cached.parent.mkdir(parents=True, exist_ok=True)
    cached.write_text(text, encoding="utf-8")
    time.sleep(0.2)
    return parse_fasta(text)


def overall_identity(a: str, b: str) -> float:
    """Cheap "is this the same protein" test.

    Exact k-mer sharing is the wrong tool for recognising a DIFFERENT strain
    or species of a protein -- that mistake is documented at length in
    src/plab/epitope.py. It is the right tool here, where the question is the
    opposite one: an isoform or duplicate entry of the same human protein
    shares nearly every k-mer, and a genuine paralog does not come close.
    """
    return kmer_identity(a, b)


def build_panel_for(target, max_panel: int = MAX_PANEL) -> dict:
    key = target.key
    record: dict = {
        "key": key, "gene": target.gene, "uniprot": target.uniprot,
        # The binding site these surfaces were mapped onto. Stored so that a
        # reader can tell whether the cached surrogates still describe this
        # epitope: rebuilding the catalogue can move a binding site, and a
        # surrogate mapped onto the old one is a confident wrong answer.
        "target_epitope": target.epitope_sequence,
        "members": [], "families": [], "state": "", "note": "",
    }

    got = epitope_window(target)
    if got is None:
        record["state"] = "unmappable"
        record["note"] = (
            "The structure's residue numbering could not be reconciled with "
            "the reference sequence, so epitope positions cannot be located "
            "in any relative. No panel can be built without inventing the "
            "correspondence.")
        return record
    window, positions = got

    families = pfam_families(target.uniprot)
    record["families"] = families
    if not families:
        record["state"] = "no_family"
        record["note"] = (
            "UniProt lists no Pfam domain family for this protein, so there "
            "is no principled set of relatives to compare against.")
        return record

    seen: dict[str, tuple[str, str, str]] = {}
    oversized = []
    for pf in families:
        members = family_members(pf)
        if len(members) > MAX_FAMILY_SIZE:
            oversized.append(f"{pf} ({len(members)} human proteins)")
            continue
        for acc, gene, seq in members:
            if acc != target.uniprot and seq:
                seen.setdefault(acc, (acc, gene, seq))

    if not seen:
        if oversized:
            # Not the same claim as having no relatives. PD-1's only Pfam
            # domain is the immunoglobulin V-set, shared by several hundred
            # human proteins; that says nothing about which ones resemble
            # PD-1, and a panel drawn from it would be arbitrary. Reporting
            # that as "no close human relative exists" would be false.
            record["state"] = "family_too_generic"
            record["note"] = (
                f"The only domain family on record for this protein is shared "
                f"by too many human proteins to identify close relatives: "
                f"{', '.join(oversized)}. A panel drawn from it would be "
                f"arbitrary rather than informative, so none was built. This "
                f"is a limit of the method, not a finding about the target.")
        else:
            record["state"] = "no_human_relatives"
            record["note"] = (
                f"No other reviewed human protein shares this protein's "
                f"domain family ({', '.join(families)}). Cross-reactivity "
                f"with a close human relative is therefore not the risk here. "
                f"Off-target binding to unrelated human proteins remains a "
                f"real risk that this method cannot assess.")
        return record

    # Rank by similarity AT THE BINDING SITE, and by IDENTITY rather than raw
    # alignment score -- score grows with length, so scoring put perlecan and
    # nephrin at the head of PD-L1's panel ahead of PD-L2.
    scored, duplicates = [], []
    for acc, gene, seq in seen.values():
        mapped = map_epitope_onto(window, positions, seq)
        if mapped is None or mapped.coverage < MIN_WINDOW_COVERAGE:
            continue
        if mapped.identity < MIN_WINDOW_IDENTITY:
            continue
        same = sum(1 for a, b in zip(mapped.surface, target.epitope_sequence)
                   if a == b) / max(len(target.epitope_sequence), 1)
        if same >= MAX_SURFACE_IDENTITY:
            duplicates.append({"uniprot": acc, "gene": gene,
                               "surface_identity": round(same, 4)})
            continue
        scored.append({
            "surface_identity": round(same, 4),
            "uniprot": acc, "gene": gene,
            "window_identity": round(mapped.identity, 4),
            "window_coverage": round(mapped.coverage, 4),
            "window_align_score": round(mapped.score, 1),
            "surrogate_surface": mapped.surface,
            "n_mapped": len(mapped.surface),
        })
    scored.sort(key=lambda m: (-m["window_identity"], -m["window_coverage"]))

    # Only now pay for full-length comparison, and only on the handful that
    # could make the panel. Doing it for every family member first meant
    # aligning ~1200 residues against ~1200 for hundreds of proteins per
    # target: hours of arithmetic to discard almost all of it.
    members = []
    for m in scored:
        if len(members) >= max_panel:
            break
        if overall_identity(target.full_sequence,
                            seen[m["uniprot"]][2]) >= MAX_SELF_IDENTITY:
            continue  # an isoform of the target, not a different protein
        members.append(m)

    record["members"] = members
    record["n_considered"] = len(seen)
    record["duplicates_excluded"] = duplicates
    if members and len(members) < MIN_MEMBERS_IF_DISTANT             and members[0]["window_identity"] < MODERATE_IDENTITY:
        record["state"] = "panel_too_weak"
        record["note"] = (
            f"The only relative found, {members[0]['gene']}, is just "
            f"{members[0]['window_identity']:.0%} identical to this target at "
            f"the binding site. A single distant relative is not a panel: "
            f"tested against controls, the margin against one such protein "
            f"does not reliably separate a real binder from a face designed "
            f"for the wrong protein. No number is reported.")
        record["rejected_members"] = members
        record["members"] = []
    elif members:
        best = members[0]["window_identity"]
        record["closest_identity"] = best
        record["strength"] = ("close" if best >= 0.5
                              else "moderate" if best >= 0.35 else "distant")
        record["state"] = "derived"
        record["note"] = (
            f"the closest {len(members)} of {len(seen)} reviewed human "
            f"proteins sharing a domain family ({', '.join(families)}), "
            f"ranked by identity at the epitope window rather than across the "
            f"whole protein. The nearest is {best:.0%} identical at this "
            f"surface, which makes this a {record['strength']} comparison.")
    else:
        record["state"] = "no_close_relatives"
        record["note"] = (
            f"{len(seen)} human proteins share a domain family with this "
            f"target, but none is more than {MIN_WINDOW_IDENTITY:.0%} "
            f"identical to it at the binding site itself. They share a fold, "
            f"not this surface. Mapping the epitope onto them would produce "
            f"a surrogate surface too unreliable to score against, so no "
            f"number is reported.")
    return record


def save_sequences(record: dict) -> None:
    """Keep each panel member's sequence on disk so scoring can read it."""
    for m in record.get("members", []):
        dest = SEQ_DIR / f"{m['uniprot']}.fasta"
        if dest.exists() and dest.stat().st_size > 0:
            continue
        try:
            text = _get(f"https://rest.uniprot.org/uniprotkb/{m['uniprot']}.fasta")
        except Exception as exc:
            print(f"      could not save {m['uniprot']}: {exc}")
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(text, encoding="utf-8")


