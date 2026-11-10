"""Asking the databases what exists, instead of recalling it.

This module holds the two lookups the catalogue research depends on:

  resolve_symbol  gene symbol -> reviewed human UniProt accession
  rank_structures accession   -> PDB entries that really contain it, with
                                something bound to it, antibodies first

It was pulled out of scripts/rcsb_lookup.py and scripts/resolve_gene.py when
the web application needed the same two steps: a person typing a gene name
into the "design your own" tab is doing exactly what I do when adding a
catalogue entry, and the two must not drift apart. The scripts now import
from here.

Both exist for the same reason. Three of this catalogue's first rejections
were PDB ids recalled from memory that held the wrong protein -- 7KQ5 is a
bacterial encapsulin, not interleukin-5 -- and the symbol SELP resolves to
selenoprotein P before it resolves to P-selectin. Neither mistake is visible
downstream: the structure parses, the epitope derives, the panel builds, and
the target is simply not the one anybody meant.
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request

GRAPHQL = "https://data.rcsb.org/graphql"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
UNIPROT = "https://rest.uniprot.org/uniprotkb"
USER_AGENT = "protein-lab/0.3 (educational in-silico design sandbox)"
HUMAN = "9606"

ANTIBODY_WORDS = ("fab", "antibody", "immunoglobulin", "nanobody", "igg",
                  "heavy chain", "light chain", "scfv", "vhh", "fv ")

# A partner this small is a ligand or a tag, not a binding partner whose
# interface is worth designing against: 3BXM pairs glutamate carboxypeptidase
# with a three-residue dipeptide, which would derive an "epitope" of almost
# nothing.
MIN_PARTNER_RESIDUES = 20

# How many entries a structure search looks at. Results come back sorted by
# resolution, so a protein crystallised many times on its own buries its
# complexes behind them -- see rank_structures.
SCAN_LIMIT = 250
# GraphQL entries per call. Naming 250 at once invites truncation.
ENTITY_BATCH = 50


def _post(url: str, payload: dict, timeout: int = 60) -> dict:
    body = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url, data=body,
        headers={"Content-Type": "application/json", "User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def _get(url: str, timeout: int = 60) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def _get_text(url: str, timeout: int = 60) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8")


# --- UniProt ----------------------------------------------------------------


def resolve_symbol(symbol: str) -> list[dict]:
    """Reviewed human entries whose gene name is exactly this symbol.

    `gene_exact`, never free text: searching "TNF" as text returns hundreds of
    proteins with TNF somewhere in the description and the first is not the
    one anybody means.
    """
    q = (f"gene_exact:{symbol} AND organism_id:{HUMAN} AND reviewed:true")
    url = (f"{UNIPROT}/search?query={urllib.parse.quote(q)}&format=json&size=10"
           "&fields=accession,id,protein_name,gene_names,length")
    try:
        data = _get(url)
    except Exception:
        return []
    return [_uniprot_record(r) for r in data.get("results", [])]


def fetch_accession(accession: str) -> dict | None:
    """One entry by accession, for when the person typed P01375 directly."""
    try:
        data = _get(f"{UNIPROT}/{accession}.json"
                    "?fields=accession,id,protein_name,gene_names,length,organism_name")
    except Exception:
        return None
    rec = _uniprot_record(data)
    org = (data.get("organism") or {}).get("scientificName")
    rec["organism"] = org
    return rec


def _uniprot_record(r: dict) -> dict:
    desc = r.get("proteinDescription") or {}
    rec = (desc.get("recommendedName") or {}).get("fullName") or {}
    if not rec.get("value"):
        subs = desc.get("submissionNames") or []
        rec = (subs[0].get("fullName") if subs else {}) or {}
    genes = [g.get("geneName", {}).get("value") for g in r.get("genes", [])]
    return {
        "accession": r.get("primaryAccession"),
        "id": r.get("uniProtkbId"),
        "name": rec.get("value") or "",
        "genes": [g for g in genes if g],
        "length": (r.get("sequence") or {}).get("length"),
    }


def uniprot_header(accession: str) -> str:
    try:
        return _get_text(f"{UNIPROT}/{accession}.fasta").splitlines()[0]
    except Exception:
        return ""


# --- RCSB -------------------------------------------------------------------

ENTITY_QUERY = """
query Entries($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    struct { title }
    rcsb_entry_info { resolution_combined experimental_method }
    polymer_entities {
      rcsb_id
      rcsb_polymer_entity { pdbx_description }
      entity_poly { rcsb_sample_sequence_length }
      rcsb_polymer_entity_container_identifiers {
        auth_asym_ids
        reference_sequence_identifiers { database_accession database_name }
      }
      rcsb_entity_source_organism { scientific_name }
    }
  }
}
"""


def entity_map(pdb_ids: list[str]) -> list[dict]:
    """What each chain in these entries actually is."""
    if not pdb_ids:
        return []
    data = _post(GRAPHQL, {"query": ENTITY_QUERY,
                           "variables": {"ids": [p.upper() for p in pdb_ids]}})
    out = []
    for entry in (data.get("data", {}).get("entries") or []):
        if entry is None:
            continue
        info = entry.get("rcsb_entry_info") or {}
        res = (info.get("resolution_combined") or [None])[0]
        ents = []
        for pe in entry.get("polymer_entities") or []:
            ids = pe.get("rcsb_polymer_entity_container_identifiers") or {}
            refs = ids.get("reference_sequence_identifiers") or []
            orgs = pe.get("rcsb_entity_source_organism") or []
            ents.append({
                "entity": pe.get("rcsb_id"),
                "chains": ids.get("auth_asym_ids") or [],
                "description": ((pe.get("rcsb_polymer_entity") or {})
                                .get("pdbx_description") or ""),
                "length": ((pe.get("entity_poly") or {})
                           .get("rcsb_sample_sequence_length")),
                "accessions": [r.get("database_accession") for r in refs
                               if r.get("database_name") == "UniProt"],
                "organism": orgs[0].get("scientific_name") if orgs else None,
            })
        out.append({
            "pdb_id": entry.get("rcsb_id"),
            "title": (entry.get("struct") or {}).get("title", ""),
            "method": (info.get("experimental_method") or ""),
            "resolution": res,
            "entities": ents,
        })
    return out


def structures_for(accession: str, limit: int = 25) -> list[str]:
    """PDB ids whose polymer entities reference this UniProt accession."""
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": ("rcsb_polymer_entity_container_identifiers"
                              ".reference_sequence_identifiers"
                              ".database_accession"),
                "operator": "exact_match",
                "value": accession,
            },
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": limit},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined",
                      "direction": "asc"}],
        },
    }
    url = SEARCH + "?json=" + urllib.parse.quote(json.dumps(query))
    try:
        data = _get(url)
    except Exception:
        return []
    return [r["identifier"] for r in data.get("result_set", [])]


def looks_like_antibody(text: str) -> bool:
    low = text.lower()
    return any(w in low for w in ANTIBODY_WORDS)


def rank_structures(accession: str, limit: int = SCAN_LIMIT,
                    top: int = 8) -> list[dict]:
    """Entries usable as a design target, best first.

    Usable means the accession really is in the entry and something large
    enough to have an interface is bound to it. Antibody complexes rank
    first -- they are the cleanest two-body interfaces and the closest thing
    to a drug already bound -- then other partners, then by resolution.
    Huge assemblies rank last rather than being hidden, because they are slow
    to parse and rarely a clean interface, but they are sometimes all there
    is.

    The scan is deliberately deep. Entries come back sorted by resolution,
    and a protein that has been crystallised on its own a hundred times puts
    every one of those ahead of its complexes: at a limit of 40, superoxide
    dismutase reported NO usable structure when it has eleven, and PSMA
    reported a three-residue dipeptide when it has nanobody complexes. Both
    would have told somebody using the build-your-own tab that their protein
    could not be worked on, which is a far worse answer than a slow one.
    """
    ids = structures_for(accession, limit)
    if not ids:
        return []
    # Batched, because a GraphQL query naming 250 entries at once is a good
    # way to be rate-limited or truncated.
    records: list[dict] = []
    for i in range(0, len(ids), ENTITY_BATCH):
        try:
            records.extend(entity_map(ids[i:i + ENTITY_BATCH]))
        except Exception:
            continue
    ranked = []
    for rec in records:
        mine = [e for e in rec["entities"] if accession in e["accessions"]]
        others = [e for e in rec["entities"]
                  if accession not in e["accessions"]
                  and (e["length"] or 0) >= MIN_PARTNER_RESIDUES]
        if not mine or not others:
            continue
        ab = [e for e in others if looks_like_antibody(e["description"])]
        big = len(rec["entities"]) > 6
        ranked.append({
            "pdb_id": rec["pdb_id"],
            "title": rec["title"],
            "resolution": rec["resolution"],
            "method": rec["method"],
            "has_antibody": bool(ab),
            "large_assembly": big,
            "target_chains": mine[0]["chains"],
            "target_length": mine[0]["length"],
            "target_description": mine[0]["description"],
            "partners": [
                {"chains": o["chains"], "length": o["length"],
                 "description": o["description"],
                 "antibody": looks_like_antibody(o["description"])}
                for o in others[:4]
            ],
            "_sort": (not ab, big, rec["resolution"] or 99),
        })
    ranked.sort(key=lambda r: r["_sort"])
    for r in ranked:
        r.pop("_sort", None)
    return ranked[:top]
