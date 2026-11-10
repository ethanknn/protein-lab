"""The protein-lab web application.

    python server.py            # then open http://127.0.0.1:8800

Built on the standard library alone. That is deliberate: the scientific core
already avoids fragile dependencies, and a web framework would add install
risk for a server that needs to do very little -- serve one page, hand out
JSON, and run jobs on background threads.

Design runs are long (minutes with the language model enabled), so the API is
asynchronous: POST a job, poll its status, read the result when it finishes.
Progress is reported per generation so the page can draw the search as it
happens rather than freezing on a spinner.
"""

from __future__ import annotations

import json
import mimetypes
import sys
import threading
import traceback
import uuid
from datetime import datetime, timezone
from functools import lru_cache
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlparse

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from plab.catalog import load_verified  # noqa: E402
from plab.design.loop import DesignLoop, LoopConfig  # noqa: E402
from plab.design.scaffolds import get_scaffold, list_scaffolds  # noqa: E402
from plab.scoring.composite import Evaluator  # noqa: E402
from plab.targets import (  # noqa: E402
    available_targets, load_off_targets, load_target,
)

STRUCT_DIR = ROOT / "data" / "structures"
WEB_DIR = ROOT / "web"
HOST, PORT = "127.0.0.1", 8800

MAX_GENERATIONS = 60
MAX_POPULATION = 96
MAX_ACTIVE_JOBS = 3


# --- job registry -----------------------------------------------------------


class Job:
    def __init__(self, job_id: str, params: dict) -> None:
        self.id = job_id
        self.params = params
        self.state = "queued"          # queued | running | done | failed
        self.stage = "waiting to start"
        self.progress = 0.0
        self.error: str | None = None
        self.result: dict | None = None
        self.trajectory: list[dict] = []
        self.calibration: dict | None = None
        self.reference: dict | None = None
        self.started = datetime.now(timezone.utc).isoformat(timespec="seconds")
        self.lock = threading.Lock()

    def snapshot(self) -> dict:
        with self.lock:
            return {
                "id": self.id, "state": self.state, "stage": self.stage,
                "progress": round(self.progress, 4), "error": self.error,
                "params": self.params, "started": self.started,
                "trajectory": list(self.trajectory),
                "calibration": self.calibration,
                "reference": self.reference,
                "result": self.result,
            }

    def set(self, **kw) -> None:
        with self.lock:
            for k, v in kw.items():
                setattr(self, k, v)


JOBS: dict[str, Job] = {}
JOBS_LOCK = threading.Lock()


def active_job_count() -> int:
    with JOBS_LOCK:
        return sum(1 for j in JOBS.values() if j.state in ("queued", "running"))


# --- the worker -------------------------------------------------------------


# Targets somebody built in this session. Deliberately in memory and
# deliberately NOT written into data/catalog: that file is a reviewed
# artefact with a written rationale and a clinical caveat per entry, and an
# experiment somebody ran once is not one. They live as long as the server
# does, which the page says.
CUSTOM: dict[str, object] = {}
MAX_CUSTOM_TARGETS = 40
CUSTOM_LOCK = threading.Lock()


def resolve_target(key: str):
    """A catalogue target, or one built in this session."""
    with CUSTOM_LOCK:
        made = CUSTOM.get(key)
    if made is not None:
        return made.spec
    return load_target(key)


def custom_off_targets(key: str):
    with CUSTOM_LOCK:
        made = CUSTOM.get(key)
    if made is None:
        return None
    from plab.custom import off_targets_for
    return off_targets_for(made.spec, made.panel)


# --- specimen profiles ------------------------------------------------------

# Ranking walks the whole catalogue, and load_target re-reads an epitope file
# and a FASTA every time. Cached because the specimen form is re-submitted on
# every change to a tick-box, and the files cannot change under a running
# server without a restart anyway.
_cached_target = lru_cache(maxsize=256)(load_target)

MAX_UNTESTED = 12


def specimen_setup() -> dict:
    """Everything the form needs: antigens to tick, lineages to spare."""
    from plab.lineage import (EXPRESSION_PATH, LINEAGES, PRESETS,
                              load_expression)
    from plab.specimen import ANTIGEN_STATES

    exp = load_expression()
    cat = load_verified().get("targets", {})

    antigens = []
    for key, rec in sorted(cat.items(), key=lambda kv: kv[1]["gene"]):
        gene = rec["gene"]
        antigens.append({
            "key": key, "gene": gene,
            # CD numbers, so a clinician recognises their own flow panel.
            "aliases": exp.aliases(gene),
            "protein_name": rec["protein_name"],
            "disease": rec["disease"], "area": rec["area"],
            "has_expression": exp.has(gene),
        })

    return {
        "antigens": antigens,
        "lineages": [{
            "key": l.key, "label": l.label, "consequence": l.consequence,
            # Whether the transcript screen can say anything at all about
            # this population. False for the stem and progenitor
            # compartments, which is the most important thing on the form
            # and must not look like an oversight.
            "screened": bool(l.blood_columns or l.tissue_columns),
        } for l in LINEAGES],
        "presets": [{"key": k, "label": lab, "lineages": list(keys)}
                    for k, lab, keys in PRESETS],
        "antigen_states": list(ANTIGEN_STATES),
        "expression": {
            "source": exp.source, "version": exp.version,
            "fetched": exp.fetched, "n_genes": len(exp.values),
            "present": EXPRESSION_PATH.exists(),
        },
    }


def specimen_rank(body: dict) -> dict:
    """Rank the catalogue against one specimen profile.

    Partitioned rather than merged into one list. "CD33 is excluded because
    the marrow must recover" and "you did not test CLEC12A, and it would be
    the strongest candidate here" are different kinds of statement, and a
    single ordering would bury the second one among ninety entries the
    clinician never asked about.
    """
    from plab.lineage import load_expression
    from plab.specimen import ANTIGEN_STATES, Specimen, Variant, rank

    raw = body.get("antigens") or {}
    antigens = {str(g).upper(): s for g, s in raw.items()
                if s in ANTIGEN_STATES}
    spared = tuple(str(s) for s in (body.get("spared") or []))

    variants = []
    for v in (body.get("variants") or [])[:40]:
        try:
            pos = int(v.get("position"))
        except (TypeError, ValueError):
            continue
        variants.append(Variant(
            gene=str(v.get("gene") or "").upper(), position=pos,
            from_aa=str(v.get("from") or "").upper()[:1],
            to_aa=str(v.get("to") or "").upper()[:1]))

    spec = Specimen(antigens=antigens, spared=spared,
                    variants=tuple(variants))

    exp = load_expression()
    cat = load_verified().get("targets", {})
    areas = {str(a) for a in (body.get("areas") or [])}

    # Partitioned on the state, not on whether the gene appears in the form.
    # A box ticked "not tested" and a box left alone say the same thing, so
    # an explicit "not tested" belongs with the suggestions rather than
    # being ranked as a candidate on evidence nobody gathered.
    tested, untested = [], []
    for key, rec in cat.items():
        state = antigens.get(rec["gene"].upper())
        if state in ("positive", "dim", "negative"):
            tested.append(key)
        elif not areas or rec["area"] in areas:
            untested.append(key)

    def ranked(keys):
        return rank([_cached_target(k) for k in keys], spec, exp)

    calls = ranked(tested)
    candidates = [c.to_dict() for c in calls if not c.disqualified]
    excluded = [c.to_dict() for c in calls if c.disqualified]

    # Only the ones a panel could actually have found. An untested target
    # that would be disqualified anyway is not worth sending somebody back
    # to the flow lab for.
    suggest = [c.to_dict() for c in ranked(untested)
               if not c.disqualified][:MAX_UNTESTED]

    return {
        "n_antigens": len(antigens), "spared": list(spared),
        "candidates": candidates, "excluded": excluded,
        "untested": suggest,
        "n_untested_considered": len(untested),
        "variants": [v.label() for v in variants],
    }


def prepare_job(job: Job) -> None:
    """Build a target from an accession and a structure, as a job.

    A job rather than a request because the honest version of this takes
    time: a structure download, an epitope derivation, a domain-family
    lookup, every family member aligned at the epitope window, and then the
    jackknife that says whether the resulting selectivity score means
    anything. Thirty to ninety seconds, and the person should see which step
    it is on rather than a spinner.
    """
    from plab.custom import prepare
    from plab.epitope import EpitopeError
    p = job.params
    try:
        job.set(state="running", stage="resolving the protein", progress=0.02)

        def progress(stage: str, frac: float) -> None:
            job.set(stage=stage, progress=frac)

        made = prepare(p["accession"], p["pdb"], disease=p.get("disease", ""),
                       with_panel=bool(p.get("with_panel", True)),
                       progress=progress)
        with CUSTOM_LOCK:
            CUSTOM[made.key] = made
            # Bounded, because nothing here is ever freed otherwise and a
            # long session of exploring proteins would grow without limit.
            while len(CUSTOM) > MAX_CUSTOM_TARGETS:
                CUSTOM.pop(next(iter(CUSTOM)))
        job.set(state="done", stage="ready", progress=1.0,
                result={"key": made.key, "record": made.record,
                        "notes": made.notes,
                        "panel_state": made.spec.panel_state})
    except EpitopeError as exc:
        # The catalogue build's own wording, which was written to be read by
        # somebody deciding what to try next.
        job.set(state="failed", stage="could not build a target",
                error=str(exc))
    except Exception as exc:
        traceback.print_exc()
        job.set(state="failed", stage="could not build a target",
                error=f"{type(exc).__name__}: {exc}")


def run_job(job: Job) -> None:
    p = job.params
    try:
        job.set(state="running", stage="loading target", progress=0.02)
        target = resolve_target(p["target"])
        scaffold = get_scaffold(p.get("scaffold", "affibody_Z"))

        job.set(stage=("loading ESM-2 (first run downloads ~2.5 GB)"
                       if p.get("use_plm") else "preparing scorers"),
                progress=0.05)
        made_off = custom_off_targets(p["target"])
        evaluator = Evaluator(target, use_plm=bool(p.get("use_plm")),
                              struct_dir=STRUCT_DIR,
                              **({"off_targets": made_off}
                                 if made_off is not None else {}))

        job.set(stage="scoring the reference molecule", progress=0.10)
        ref = evaluator.reference
        if ref is not None:
            job.set(reference={
                "name": ref.name, "sequence": ref.sequence,
                "length": len(ref.sequence),
                "description": ref.description,
                "composite": round(ref.score_100, 2),
                "groups": {k: round(v * 100, 2)
                           for k, v in ref.group_scores.items()},
            })

        job.set(stage=f"calibrating against {p['decoys']} random library members",
                progress=0.16)
        calibration = evaluator.calibrate(scaffold, n_decoys=int(p["decoys"]))
        job.set(calibration=calibration.to_dict())

        gens = int(p["generations"])
        config = LoopConfig(
            population_size=int(p["population"]),
            generations=gens,
            patience=int(p.get("patience", 8)),
            seed=int(p.get("seed", 7)),
        )
        loop = DesignLoop(evaluator, scaffold, config)

        def on_generation(gen) -> None:
            frac = 0.20 + 0.75 * (gen.index / max(gens, 1))
            with job.lock:
                job.trajectory = loop.trajectory()
                job.stage = (f"generation {gen.index} of {gens} — "
                             f"best {gen.best.score_100:.1f}")
                job.progress = min(0.95, frac)

        job.set(stage="running the design loop", progress=0.20)
        loop.run(verbose=False, on_generation=on_generation)

        job.set(stage="assembling results", progress=0.96)
        job.set(result=build_result(loop, evaluator, calibration, target,
                                    scaffold),
                state="done", stage="finished", progress=1.0)

    except Exception as exc:
        job.set(state="failed",
                stage="failed",
                error=f"{type(exc).__name__}: {exc}",
                progress=1.0)
        traceback.print_exc()


def _panel_record(key: str) -> dict:
    from plab.targets import load_panels
    return load_panels().get(key, {})


def _member_identity(key: str, accession: str) -> float | None:
    for m in _panel_record(key).get("members", []):
        if m.get("uniprot") == accession:
            return m.get("window_identity")
    return None


@lru_cache(maxsize=256)
def _reference_check(key: str, face_length: int = 13) -> dict | None:
    """Score the real drug's own contact residues against the same panel.

    The single most useful calibration available: a molecule that
    demonstrably binds this target and not its relatives, measured by the
    metric that ranks the designs. When it comes out negative -- the metric
    claiming a licensed antibody prefers the off-target -- the ranking for
    that target is not to be trusted however confident the percentage looks,
    and the person reading the page is entitled to know that before they read
    the number.

    It was a line in a diagnostic script nobody but me ever runs. Nine of
    fifteen references were negative in there while the dashboard said
    nothing.
    """
    from plab.scoring.selectivity import (achievable_margin_range,
                                          build_paralog_epitopes,
                                          complementarity_kernel,
                                          margin_to_score)
    t = load_target(key)
    ref = t.paratope_sequence
    if not ref:
        return None
    panel = build_paralog_epitopes(t, load_off_targets(t))
    if not panel:
        return None

    epi = t.epitope_sequence
    on = complementarity_kernel(ref, epi)
    worst, worst_val = max(
        ((g, complementarity_kernel(ref, s)) for g, s in panel.items()),
        key=lambda kv: kv[1])
    margin = on - worst_val

    # The ceiling is precomputed; only the random floor is measured here,
    # and that is a few hundred count comparisons.
    hi = (t.panel_ceiling or {}).get(str(face_length))
    lo, swept_hi = achievable_margin_range(epi, panel, face_length)
    if hi is None:
        hi = swept_hi

    return {
        "binder": t.reference_binder_name,
        "margin": round(margin, 4),
        "score": round(100 * margin_to_score(margin, lo, hi), 1),
        "closest": worst,
        "counter_selective": margin < 0,
    }


def _panel_state(key: str) -> str:
    """Panel state for the catalogue listing, without loading a whole target."""
    from plab.targets import load_panels
    rec = available_targets().get(key, {})
    panel = load_panels().get(key, {})
    if rec.get("off_targets") and panel.get("members"):
        return "curated+derived"
    if rec.get("off_targets"):
        return "curated"
    return panel.get("state", "none")


def build_result(loop, evaluator, calibration, target, scaffold,
                 top_n: int = 12) -> dict:
    candidates = []
    for rank, ind in enumerate(loop.top(top_n), 1):
        groups = evaluator.group_scores(ind.card, exclude_self_referential=True)
        rec = groups.get("recognition", 0.0) * 100
        interp = calibration.interpret(ind.card, rec)
        candidates.append({
            "rank": rank,
            "sequence": ind.sequence,
            "designed": ind.card.metadata.get("designed_residues", ""),
            "annotated": scaffold.annotate(ind.sequence),
            "generation": ind.generation,
            "origin": ind.origin,
            "composite": round(ind.score_100, 2),
            "groups": {k: round(v * 100, 2) for k, v in groups.items()},
            "recognition": round(rec, 2),
            "recognition_percentile": round(interp["recognition_percentile"], 1),
            "recognition_z": round(interp["recognition_z"], 2),
            "percent_of_reference": (
                None if evaluator.recognition_percent_of_reference(ind.card) is None
                else round(evaluator.recognition_percent_of_reference(ind.card), 1)),
            "flags": ind.card.flags,
            # Why this candidate cannot be ranked first. Carried per
            # candidate rather than per run, because a table that silently
            # reorders rows is worse than one that explains the order.
            "disqualifying": [
                {"name": c.name, "why": c.disqualifying}
                for c in ind.card.disqualified
            ],
            "n_ranked": sum(1 for c in ind.card.components
                            if c.weight > 0 and c.ranked),
            "components": [
                {
                    "name": c.name,
                    "normalised": round(c.normalised * 100, 1),
                    "raw": None if c.raw != c.raw else round(float(c.raw), 4),
                    "units": c.units, "weight": c.weight,
                    "detail": c.detail, "caveat": c.caveat,
                    "ranked": c.ranked, "note": c.note,
                    "group": group_of(c.name),
                }
                for c in ind.card.components if c.weight > 0
            ],
        })

    return {
        "candidates": candidates,
        "stop_reason": loop.stop_reason,
        "elapsed_s": round(getattr(loop, "total_elapsed", 0.0), 1),
        "n_unique": loop._n_evaluated,
        "n_cached": loop._n_cached,
        "trajectory": loop.trajectory(),
        "scaffold": {
            "name": scaffold.name,
            "description": scaffold.description,
            "framework": scaffold.framework,
            "ss": scaffold.secondary_structure,
            "variable_1based": [p + 1 for p in scaffold.variable_positions],
            "provenance": scaffold.provenance,
            "search_space": scaffold.search_space,
            "pdb_template": scaffold.pdb_template,
        },
        "plm_status": (evaluator.plm.runner.status if evaluator.plm
                       else "disabled"),
        "paralogs": evaluator.paralog_epitopes,
        "gene": target.gene,
        # Did the search ever find a candidate that is not disqualified? If
        # not, ranking degenerates to the plain composite and the top row is
        # the least bad of a bad set -- which the page has to say out loud,
        # because the table looks identical either way.
        "any_qualified": any(
            not i.card.disqualified
            for g in loop.history for i in g.individuals),
        "n_disqualified_in_top": sum(
            1 for i in loop.top(top_n) if i.card.disqualified),
        # The reachable ceiling for THIS scaffold's face, and its error
        # bar. Both are needed to say why a selectivity score is absent:
        # a ceiling of zero means the proteins are twins, while a ceiling
        # no bigger than its own wobble means the measurement ran out of
        # resolution, and those are different things to tell someone.
        "panel_ceiling": (target.panel_ceiling or {}).get(
            str(len(scaffold.variable_positions))),
        "panel_ceiling_noise": (target.panel_ceiling_noise or {}).get(
            str(len(scaffold.variable_positions))),
        "panel_state": target.panel_state,
        "panel_note": target.panel_note,
        "off_targets": [
            {"gene": o.gene, "uniprot": o.uniprot, "concern": o.concern,
             "structural": o.structural}
            for o in evaluator.off_targets
        ],
    }


def group_of(name: str) -> str:
    from plab.scoring.composite import GROUPS
    for g, names in GROUPS.items():
        if name in names:
            return g
    return "other"


# --- target detail ----------------------------------------------------------


def custom_detail(made) -> dict:
    """The same payload shape a catalogue target returns.

    Same shape on purpose: a target somebody built should be examined with
    the same page, the same caveats and the same cross-reactivity panel as a
    curated one. What it does NOT get is a rationale or a clinical note,
    because nobody wrote them, and `custom` tells the page to say so rather
    than leave an empty space that reads like "no concerns".
    """
    t, rec, panel = made.spec, made.record, made.panel
    members = {m["uniprot"]: m for m in panel.get("members", [])}
    from plab.custom import off_targets_for
    return {
        "key": t.key, "gene": t.gene, "protein_name": t.protein_name,
        "disease": t.disease, "area": t.area, "uniprot": t.uniprot,
        "rationale": "", "clinical_note": "",
        "brief": t.design_brief(),
        "custom": True,
        "complex_pdb": t.complex_pdb, "apo_pdb": None,
        "resolution": t.source_resolution,
        "structure_title": t.structure_title,
        "known_binder": t.reference_binder_name,
        "binder_kind": t.binder_kind,
        "target_chain": t.target_chain, "binder_chains": t.binder_chains,
        "epitope_sequence": t.epitope_sequence,
        "paratope_sequence": t.paratope_sequence,
        "epitope": [
            {"label": r.label, "resseq": r.resseq, "aa": r.aa,
             "chain": r.chain, "dist": r.min_distance_to_binder}
            for r in t.epitope
        ],
        "epitope_charge": round(t.epitope_net_charge, 2),
        "epitope_radius": t.epitope_radius,
        "cations": t.cation_positions,
        "surface": t.surface,
        "descriptors": t.descriptors,
        "contacts": t.contacts,
        "off_target_accessions": t.off_target_accessions,
        "has_paralog_panel": bool(t.off_target_accessions),
        "panel_state": t.panel_state,
        "panel_note": t.panel_note,
        "panel_families": t.panel_families,
        "panel_considered": t.panel_considered,
        "panel_strength": panel.get("strength"),
        "panel_rejected": panel.get("rejected_members", []),
        "panel_closest_identity": panel.get("closest_identity"),
        "reference_check": None,
        "panel_ceiling": (t.panel_ceiling or {}).get("13"),
        "panel_ceiling_noise": (t.panel_ceiling_noise or {}).get("13"),
        "panel": [
            {"gene": o.gene, "uniprot": o.uniprot, "concern": o.concern,
             "structural": o.structural,
             "window_identity": members.get(o.uniprot, {}).get("window_identity")}
            for o in off_targets_for(t, panel)
        ],
        "sequence_length": rec.get("sequence_length"),
        "caveats": list(rec.get("warnings", [])),
        "target_identity": rec.get("target_identity"),
        "target_coverage": rec.get("target_coverage"),
        "identity_method": rec.get("identity_method"),
    }


def target_detail(key: str) -> dict:
    with CUSTOM_LOCK:
        made = CUSTOM.get(key)
    if made is not None:
        return custom_detail(made)
    t = load_target(key)
    rec = available_targets()[t.key]
    return {
        "key": t.key, "gene": t.gene, "protein_name": t.protein_name,
        "disease": t.disease, "area": t.area, "uniprot": t.uniprot,
        "rationale": t.rationale, "clinical_note": t.clinical_note,
        "brief": t.design_brief(),
        "complex_pdb": t.complex_pdb, "apo_pdb": t.apo_pdb,
        "resolution": t.source_resolution,
        "structure_title": t.structure_title,
        "known_binder": t.reference_binder_name,
        "binder_kind": rec.get("binder_kind", "protein"),
        "target_chain": t.target_chain, "binder_chains": t.binder_chains,
        "epitope_sequence": t.epitope_sequence,
        "paratope_sequence": t.paratope_sequence,
        "epitope": [
            {"label": r.label, "resseq": r.resseq, "aa": r.aa,
             "chain": r.chain, "dist": r.min_distance_to_binder}
            for r in t.epitope
        ],
        "epitope_charge": round(t.epitope_net_charge, 2),
        "epitope_radius": t.epitope_radius,
        "cations": t.cation_positions,
        "surface": t.surface,
        "descriptors": t.descriptors,
        "contacts": t.contacts,
        "off_target_accessions": t.off_target_accessions,
        "has_paralog_panel": bool(t.off_target_accessions),
        "panel_state": t.panel_state,
        "panel_note": t.panel_note,
        "panel_families": t.panel_families,
        "panel_considered": t.panel_considered,
        "panel_strength": _panel_record(key).get("strength"),
        "panel_rejected": _panel_record(key).get("rejected_members", []),
        "panel_closest_identity": _panel_record(key).get("closest_identity"),
        # How the real drug scores under the metric that ranks the
        # designs. A negative margin here means this measure claims a
        # licensed antibody prefers the off-target, which is a statement
        # about the measure and belongs on the page beside its output.
        "reference_check": _reference_check(key),
        "panel_ceiling": (t.panel_ceiling or {}).get("13"),
        "panel_ceiling_noise": (t.panel_ceiling_noise or {}).get("13"),
        "panel": [
            {"gene": o.gene, "uniprot": o.uniprot, "concern": o.concern,
             "structural": o.structural,
             # How identical this protein is to the target AT THE BINDING
             # SITE. Without it a panel of 30%-identical relatives and a
             # panel of 70%-identical ones look like equally strong evidence.
             "window_identity": _member_identity(key, o.uniprot)}
            for o in load_off_targets(t)
        ],
        "sequence_length": rec.get("sequence_length"),
        # Everything verification noticed but did not reject over: low
        # resolution, a structure that is a variant of the human protein, a
        # second binder set aside. These belong in front of the user, not in
        # a log file -- they bound how far the numbers can be trusted.
        "caveats": rec.get("warnings", []),
        "target_identity": rec.get("target_identity"),
        "target_coverage": rec.get("target_coverage"),
        "identity_method": rec.get("identity_method"),
    }


# --- HTTP -------------------------------------------------------------------


class Handler(BaseHTTPRequestHandler):
    server_version = "protein-lab/0.2"

    def log_message(self, fmt, *args):  # quieter console
        if "/api/job/" not in (args[0] if args else ""):
            sys.stderr.write(f"  {self.address_string()} {fmt % args}\n")

    # -- helpers --
    def _json(self, payload, status: int = 200) -> None:
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _text(self, body: str, ctype: str = "text/plain; charset=utf-8",
              status: int = 200, cache: bool = False) -> None:
        data = body.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control",
                         "public, max-age=3600" if cache else "no-store")
        self.end_headers()
        self.wfile.write(data)

    def _error(self, msg: str, status: int = 400) -> None:
        self._json({"error": msg}, status)

    # -- routes --
    def do_GET(self) -> None:
        path = unquote(urlparse(self.path).path)
        try:
            if path in ("/", "/index.html"):
                return self._serve_file(WEB_DIR / "index.html")
            if path.startswith("/static/"):
                return self._serve_file(WEB_DIR / "static" / path[len("/static/"):])
            if path == "/api/catalog":
                return self._json(self.catalog_payload())
            if path.startswith("/api/target/"):
                key = path[len("/api/target/"):]
                with CUSTOM_LOCK:
                    known = key in CUSTOM
                if key not in available_targets() and not known:
                    return self._error(f"unknown target {key!r}", 404)
                return self._json(target_detail(key))
            if path.startswith("/api/structure/"):
                return self._serve_structure(path[len("/api/structure/"):])
            if path.startswith("/api/job/"):
                job = JOBS.get(path[len("/api/job/"):])
                if job is None:
                    return self._error("no such job", 404)
                return self._json(job.snapshot())
            if path == "/api/custom/search":
                from plab.custom import search
                q = ""
                for part in urlparse(self.path).query.split("&"):
                    if part.startswith("q="):
                        q = unquote(part[2:])
                return self._json(search(q))
            if path == "/api/custom/structures":
                from plab.custom import structures
                acc = ""
                for part in urlparse(self.path).query.split("&"):
                    if part.startswith("accession="):
                        acc = unquote(part[len("accession="):])
                if not acc:
                    return self._error("accession required")
                return self._json({"accession": acc,
                                   "structures": structures(acc)})
            if path == "/api/custom/list":
                with CUSTOM_LOCK:
                    made = [{"key": k, "gene": v.record["gene"],
                             "protein_name": v.record["protein_name"],
                             "pdb": v.record["pdb_complex"],
                             "disease": v.record["disease"],
                             "epitope_sequence": v.record["epitope_sequence"],
                             "n_epitope_residues": v.record["n_epitope_residues"],
                             "panel_state": v.spec.panel_state}
                            for k, v in CUSTOM.items()]
                return self._json({"targets": made})
            if path == "/api/specimen/setup":
                return self._json(specimen_setup())
            if path == "/api/scaffolds":
                return self._json([
                    {"name": s.name, "description": s.description,
                     "length": s.length, "n_variable": s.n_variable,
                     "search_space": s.search_space,
                     "provenance": s.provenance}
                    for s in list_scaffolds()
                ])
            return self._error("not found", 404)
        except Exception as exc:
            traceback.print_exc()
            self._error(f"{type(exc).__name__}: {exc}", 500)

    def do_POST(self) -> None:
        path = unquote(urlparse(self.path).path)
        if path not in ("/api/run", "/api/custom/prepare",
                        "/api/specimen/rank"):
            return self._error("not found", 404)
        try:
            length = int(self.headers.get("Content-Length") or 0)
            body = json.loads(self.rfile.read(length) or b"{}")
        except Exception:
            return self._error("malformed JSON body")

        if path == "/api/specimen/rank":
            # Synchronous: this reads cached files and does no alignment, so
            # it returns in well under a second and does not need the job
            # machinery that design runs use.
            return self._json(specimen_rank(body))

        if path == "/api/custom/prepare":
            acc = (body.get("accession") or "").strip().upper()
            pdb = (body.get("pdb") or "").strip().upper()
            if not acc or not pdb:
                return self._error("accession and pdb are both required")
            if active_job_count() >= MAX_ACTIVE_JOBS:
                return self._error(
                    f"{MAX_ACTIVE_JOBS} jobs are already in flight", 429)
            params = {"accession": acc, "pdb": pdb,
                      "disease": (body.get("disease") or "").strip()[:120],
                      "with_panel": bool(body.get("with_panel", True))}
            job = Job(uuid.uuid4().hex[:12], params)
            with JOBS_LOCK:
                JOBS[job.id] = job
            threading.Thread(target=prepare_job, args=(job,),
                             daemon=True).start()
            return self._json({"job": job.id, "params": params}, 202)

        key = body.get("target")
        with CUSTOM_LOCK:
            is_custom = key in CUSTOM
        if key not in available_targets() and not is_custom:
            return self._error(f"unknown target {key!r}")
        if active_job_count() >= MAX_ACTIVE_JOBS:
            return self._error(
                f"{MAX_ACTIVE_JOBS} runs are already in flight; wait for one "
                f"to finish before starting another", 429)

        try:
            params = {
                "target": key,
                "scaffold": body.get("scaffold", "affibody_Z"),
                "generations": max(1, min(MAX_GENERATIONS,
                                          int(body.get("generations", 20)))),
                "population": max(8, min(MAX_POPULATION,
                                         int(body.get("population", 40)))),
                "decoys": max(40, min(400, int(body.get("decoys", 150)))),
                "seed": int(body.get("seed", 7)),
                "use_plm": bool(body.get("use_plm", False)),
            }
        except (TypeError, ValueError):
            return self._error("bad parameter types")

        job = Job(uuid.uuid4().hex[:12], params)
        with JOBS_LOCK:
            JOBS[job.id] = job
        threading.Thread(target=run_job, args=(job,), daemon=True).start()
        self._json({"job": job.id, "params": params}, 202)

    # -- static & data --
    def _serve_file(self, path: Path) -> None:
        path = path.resolve()
        if not str(path).startswith(str(WEB_DIR.resolve())) or not path.is_file():
            return self._error("not found", 404)
        ctype = mimetypes.guess_type(path.name)[0] or "application/octet-stream"
        if ctype.startswith("text/") or ctype in ("application/javascript",
                                                  "application/json"):
            ctype += "; charset=utf-8"
        # The app's own source is never cached. This is a local development
        # server whose files change while it runs, and an hour-long cache on
        # app.js means edits appear to do nothing -- which looks exactly like
        # a bug in the code you just changed. Structures and other downloaded
        # data are immutable and still cached.
        self._text(path.read_text(encoding="utf-8"), ctype, cache=False)

    def _serve_structure(self, name: str) -> None:
        """Serve a PDB file for the 3D viewer, optionally only some chains:
        /api/structure/8W9J?chains=C,H,L"""
        pdb_id = name.replace(".pdb", "").upper()
        if not pdb_id.isalnum() or len(pdb_id) > 8:
            return self._error("bad structure id")
        path = STRUCT_DIR / f"{pdb_id}.pdb"
        cif = STRUCT_DIR / f"{pdb_id}.cif"
        if not path.exists() and not cif.exists():
            return self._error(f"{pdb_id} not downloaded", 404)
        query = urlparse(self.path).query
        chains = None
        for part in query.split("&"):
            if part.startswith("chains="):
                chains = [c for c in unquote(part[7:]).split(",") if c]
        if not path.exists():
            # mmCIF-only entry. Parse it and re-emit PDB records rather than
            # teaching the viewer and the chain-subsetting code a second
            # format; see pdbio.to_pdb_text.
            from plab.pdbio import parse_structure, to_pdb_text
            text = to_pdb_text(parse_structure(cif), chains=chains)
        else:
            from plab.report.dashboard import extract_pdb_subset
            text = extract_pdb_subset(path, chains=chains)
        self._text(text, "chemical/x-pdb; charset=utf-8", cache=True)

    @staticmethod
    def catalog_payload() -> dict:
        cat = load_verified()
        targets = cat.get("targets", {})
        grouped: dict[str, list[dict]] = {}
        for key, rec in targets.items():
            grouped.setdefault(rec["area"], []).append({
                "key": key, "gene": rec["gene"], "disease": rec["disease"],
                "area": rec["area"], "protein_name": rec["protein_name"],
                "uniprot": rec["uniprot"], "pdb": rec["pdb_complex"],
                "known_binder": rec["known_binder"],
                "resolution": rec.get("resolution"),
                "epitope_sequence": rec.get("epitope_sequence", ""),
                "n_epitope_residues": rec.get("n_epitope_residues", 0),
                "surface": rec.get("surface", {}),
                "binder_kind": rec.get("binder_kind", "protein"),
                "has_paralog_panel": bool(rec.get("off_targets")),
                "panel_state": _panel_state(key),
            })
        for v in grouped.values():
            v.sort(key=lambda e: e["gene"])
        return {
            "built_at": cat.get("built_at"),
            "n_verified": len(targets),
            "n_candidates": cat.get("n_candidates"),
            "areas": dict(sorted(grouped.items())),
            "rejected": cat.get("rejected", []),
        }


def main() -> int:
    if not (WEB_DIR / "index.html").exists():
        print(f"Missing {WEB_DIR / 'index.html'}")
        return 1
    cat = load_verified()
    n = len(cat.get("targets", {}))
    if n == 0:
        print("No verified targets. Run: python scripts/build_catalog.py")
        return 1

    srv = ThreadingHTTPServer((HOST, PORT), Handler)
    print(f"protein-lab serving {n} verified targets")
    print(f"  http://{HOST}:{PORT}")
    print("  Ctrl-C to stop\n")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nstopping")
    finally:
        srv.server_close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
