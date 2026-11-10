"""Builds the self-contained HTML dashboard.

Everything is embedded -- structures, run data, charts -- so the output file
works offline, can be emailed, and will still open years from now. The only
external resource is the 3Dmol.js viewer from a CDN.

Design intent: the dashboard must make the limits of the method as visible as
the results. The caveat for every score component is in the page, the decoy
distribution is plotted next to the candidate so the reader can see how much
of the score is signal, and the discriminating-power diagnostic is given the
same visual weight as the headline number.
"""

from __future__ import annotations

import json
from pathlib import Path

VIEWER_CDN = "https://cdn.jsdelivr.net/npm/3dmol@2.4.2/build/3Dmol-min.js"


def extract_pdb_subset(path: Path, chains: list[str] | None = None,
                       max_resseq: int | None = None,
                       first_model_only: bool = True) -> str:
    """Pull a minimal PDB out of a large one, to keep the HTML small."""
    keep: list[str] = []
    in_later_model = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        rec = line[:6]
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
        chain = line[21]
        if chains and chain not in chains:
            continue
        if max_resseq is not None:
            try:
                if int(line[22:26]) > max_resseq:
                    continue
            except ValueError:
                continue
        # Drop solvent; it adds bulk and nothing visual.
        if line[17:20].strip() in ("HOH", "DOD", "WAT"):
            continue
        keep.append(line)
    keep.append("END")
    return "\n".join(keep)


def build_dashboard(loop, evaluator, calibration, out_path: Path,
                    struct_dir: Path, top_n: int = 10,
                    artifact_mode: bool = False) -> Path:
    """Write the dashboard.

    artifact_mode emits a bare fragment for publishing as an Artifact, which
    supplies its own document skeleton. The default writes a standalone file
    that opens from disk.
    """
    target = evaluator.target
    scaffold = loop.scaffold
    ref = evaluator.reference

    # --- gather the data ------------------------------------------------
    top = loop.top(top_n)
    candidates = []
    for rank, ind in enumerate(top, 1):
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
            "components": [
                {
                    "name": c.name,
                    "normalised": round(c.normalised * 100, 1),
                    "raw": None if c.raw != c.raw else round(float(c.raw), 4),
                    "units": c.units,
                    "weight": c.weight,
                    "detail": c.detail,
                    "caveat": c.caveat,
                    "group": _group_of(c.name),
                    "ranked": c.ranked,
                    "note": c.note,
                }
                for c in ind.card.components if c.weight > 0
            ],
            "n_ranked": sum(1 for c in ind.card.components
                            if c.weight > 0 and c.ranked),
        })

    epitope = [
        {"label": r.label, "resseq": r.resseq, "aa": r.aa, "chain": r.chain,
         "dist": r.min_distance_to_binder}
        for r in target.epitope
    ]

    payload = {
        "target": {
            "gene": target.gene,
            "protein": target.protein_name,
            "disease": target.disease,
            "uniprot": target.uniprot,
            "brief": target.design_brief(),
            "rationale": target.rationale,
            "clinical_note": target.clinical_note,
            "epitope": epitope,
            "epitope_sequence": target.epitope_sequence,
            "epitope_charge": round(target.epitope_net_charge, 2),
            "epitope_radius": target.epitope_radius,
            "arginines": target.arginine_positions,
            "apo_pdb": target.apo_pdb,
            "complex_pdb": target.complex_pdb,
            "resolution": target.source_resolution,
            "target_chain": target.target_chain,
            "paratope": target.paratope_sequence,
            "contacts": target.contacts,
        },
        "scaffold": {
            "name": scaffold.name,
            "description": scaffold.description,
            "framework": scaffold.framework,
            "ss": scaffold.secondary_structure,
            "variable_1based": [p + 1 for p in scaffold.variable_positions],
            "provenance": scaffold.provenance,
            "search_space": scaffold.search_space,
            "pdb_template": scaffold.pdb_template,
            "pdb_chain": scaffold.pdb_chain,
        },
        "reference": None if ref is None else {
            "name": ref.name,
            "sequence": ref.sequence,
            "length": len(ref.sequence),
            "description": ref.description,
            "composite": round(ref.score_100, 2),
            "groups": {k: round(v * 100, 2) for k, v in ref.group_scores.items()},
            "provenance": scaffold.provenance,
        },
        "calibration": calibration.to_dict(),
        "run": {
            "config": loop.config.to_dict(),
            "stop_reason": loop.stop_reason,
            "elapsed_s": round(getattr(loop, "total_elapsed", 0.0), 1),
            "n_unique": loop._n_evaluated,
            "n_cached": loop._n_cached,
            "trajectory": loop.trajectory(),
        },
        "candidates": candidates,
        "plm_status": (evaluator.plm.runner.status if evaluator.plm
                       else "disabled"),
        "paralogs": evaluator.paralog_epitopes,
    }

    # --- structures ------------------------------------------------------
    def _structure_text(pdb_id: str, chains: list[str]) -> str | None:
        """PDB records for the viewer, from whichever format exists.

        An mmCIF-only entry is parsed and re-emitted rather than skipped.
        Looking only for .pdb left the two cif-only catalogue targets with
        no structure in the static report at all.
        """
        legacy = struct_dir / f"{pdb_id}.pdb"
        if legacy.exists():
            return extract_pdb_subset(legacy, chains=chains)
        cif = struct_dir / f"{pdb_id}.cif"
        if cif.exists():
            from ..pdbio import parse_structure, to_pdb_text
            return to_pdb_text(parse_structure(cif), chains=chains)
        return None

    structures = {}
    complex_pdb = struct_dir / f"{target.complex_pdb}.pdb"
    if not complex_pdb.exists():
        complex_pdb = struct_dir / f"{target.complex_pdb}.cif"
    if complex_pdb.exists():
        heavy = next((c for c, r in (target.chain_roles or {}).items()
                      if r == "Fab_heavy"), "H")
        light = next((c for c, r in (target.chain_roles or {}).items()
                      if r == "Fab_light"), "L")
        structures["complex"] = _structure_text(
            target.complex_pdb,
            [target.target_chain or "C", heavy, light])
        structures["complex_chains"] = {
            "target": target.target_chain or "C",
            "heavy": heavy, "light": light,
        }
    if target.apo_pdb:
        apo = _structure_text(target.apo_pdb, ["A"])
        if apo:
            structures["apo"] = apo
    if scaffold.pdb_template:
        tmpl = _structure_text(scaffold.pdb_template,
                               [scaffold.pdb_chain or "A"])
        if tmpl:
            structures["scaffold"] = tmpl

    html = _render(payload, structures, artifact_mode=artifact_mode)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")
    return out_path


def _group_of(name: str) -> str:
    from ..scoring.composite import GROUPS
    for g, names in GROUPS.items():
        if name in names:
            return g
    return "other"


def _render(data: dict, structures: dict, artifact_mode: bool = False) -> str:
    data_json = json.dumps(data)
    struct_json = json.dumps(structures)
    fragment = _TEMPLATE.replace("/*__DATA__*/", data_json) \
                        .replace("/*__STRUCTURES__*/", struct_json) \
                        .replace("__VIEWER_CDN__", VIEWER_CDN)
    return fragment if artifact_mode else _wrap_local(fragment)


_TEMPLATE = r"""<title>CLEC12A Binder Design</title>
<meta name="claude:viewport-fit-sides" content="cover">
<script src="__VIEWER_CDN__"></script>
<style>
:root{
  color-scheme: light;
  --surface-1:#fcfcfb; --plane:#f9f9f7;
  --ink-1:#0b0b0b; --ink-2:#52514e; --ink-3:#898781;
  --grid:#e1e0d9; --axis:#c3c2b7; --border:rgba(11,11,11,0.10);
  --s1:#2a78d6; --s2:#eb6834; --s3:#1baf7a; --s4:#eda100;
  --seq-100:#cde2fb; --seq-250:#86b6ef; --seq-400:#3987e5; --seq-550:#1c5cab;
  --good:#0ca30c; --warning:#fab219; --serious:#ec835a; --critical:#d03b3b;
  --radius:10px;
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    color-scheme: dark;
    --surface-1:#1a1a19; --plane:#0d0d0d;
    --ink-1:#f0efec; --ink-2:#c3c2b7; --ink-3:#898781;
    --grid:#2c2c2a; --axis:#383835; --border:rgba(255,255,255,0.10);
    --s1:#3987e5; --s2:#d95926; --s3:#199e70; --s4:#c98500;
    --seq-100:#104281; --seq-250:#184f95; --seq-400:#2a78d6; --seq-550:#86b6ef;
  }
}
:root[data-theme="dark"]{
  color-scheme: dark;
  --surface-1:#1a1a19; --plane:#0d0d0d;
  --ink-1:#f0efec; --ink-2:#c3c2b7; --ink-3:#898781;
  --grid:#2c2c2a; --axis:#383835; --border:rgba(255,255,255,0.10);
  --s1:#3987e5; --s2:#d95926; --s3:#199e70; --s4:#c98500;
  --seq-100:#104281; --seq-250:#184f95; --seq-400:#2a78d6; --seq-550:#86b6ef;
}
*{box-sizing:border-box}
html,body{margin:0;padding:0}
body{
  background:var(--plane); color:var(--ink-1);
  font:15px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif;
  /* Side gutter clears a phone system bar at the screen edge. */
  padding-left:max(16px, env(safe-area-inset-left, 0px));
  padding-right:max(16px, env(safe-area-inset-right, 0px));
  padding-block:0 72px;
}
.wrap{max-width:1180px;margin:0 auto}
header{padding:40px 0 20px}
h1{font-size:clamp(24px,3.4vw,36px);line-height:1.15;margin:0 0 6px;letter-spacing:-.02em}
.sub{color:var(--ink-2);font-size:16px;max-width:70ch;margin:0}
h2{font-size:19px;margin:0 0 4px;letter-spacing:-.01em}
h3{font-size:15px;margin:0 0 8px;color:var(--ink-2);font-weight:600}
.secnote{color:var(--ink-3);font-size:13.5px;margin:0 0 18px;max-width:78ch}
section{margin:38px 0}
.card{
  background:var(--surface-1);border:1px solid var(--border);
  border-radius:var(--radius);padding:20px;
}
.grid{display:grid;gap:16px}
.g2{grid-template-columns:repeat(auto-fit,minmax(320px,1fr))}
.g3{grid-template-columns:repeat(auto-fit,minmax(220px,1fr))}
.g4{grid-template-columns:repeat(auto-fit,minmax(170px,1fr))}

/* --- the honesty banner: deliberately not dismissible --- */
.banner{
  border:1px solid var(--border);border-left:3px solid var(--warning);
  background:var(--surface-1);border-radius:var(--radius);
  padding:18px 20px;margin:24px 0 8px;
}
.banner h2{font-size:16px;margin:0 0 8px;display:flex;gap:8px;align-items:center}
.banner p{margin:0 0 10px;color:var(--ink-2);font-size:14px;max-width:88ch}
.banner p:last-child{margin-bottom:0}
.banner strong{color:var(--ink-1)}

.tile{background:var(--surface-1);border:1px solid var(--border);
  border-radius:var(--radius);padding:16px 18px}
.tile .lab{font-size:12px;letter-spacing:.04em;text-transform:uppercase;
  color:var(--ink-3);margin-bottom:6px}
.tile .val{font-size:30px;font-weight:650;line-height:1.1;letter-spacing:-.02em}
.tile .foot{font-size:12.5px;color:var(--ink-2);margin-top:6px;line-height:1.4}

.seq{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
  font-size:12.5px;word-break:break-all;line-height:1.7}
.seq .v{background:var(--seq-250);color:var(--ink-1);border-radius:3px;
  padding:1px 2px;font-weight:700}
:root[data-theme="dark"] .seq .v,
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]) .seq .v{color:#fff}}

table{width:100%;border-collapse:collapse;font-size:13.5px}
th,td{text-align:right;padding:8px 10px;border-bottom:1px solid var(--grid)}
th:first-child,td:first-child{text-align:left}
th{color:var(--ink-3);font-weight:600;font-size:12px;text-transform:uppercase;
  letter-spacing:.04em;white-space:nowrap}
td{font-variant-numeric:tabular-nums}
tbody tr:hover{background:var(--plane)}
tbody tr.sel{background:var(--plane);box-shadow:inset 2px 0 0 var(--s1)}
tbody tr{cursor:pointer}

.viewer{width:100%;height:380px;position:relative;border-radius:8px;
  overflow:hidden;background:var(--plane);border:1px solid var(--border)}
.vlegend{display:flex;flex-wrap:wrap;gap:12px;margin-top:10px;font-size:12.5px;
  color:var(--ink-2)}
.vlegend span{display:inline-flex;align-items:center;gap:6px}
.sw{width:11px;height:11px;border-radius:3px;display:inline-block;
  border:1px solid var(--border)}

.ctrls{display:flex;flex-wrap:wrap;gap:8px;margin-bottom:14px;align-items:center}
button,select{
  font:inherit;font-size:13px;color:var(--ink-1);background:var(--surface-1);
  border:1px solid var(--axis);border-radius:7px;padding:6px 11px;cursor:pointer;
}
button:hover,select:hover{background:var(--plane)}
button[aria-pressed="true"]{background:var(--s1);border-color:var(--s1);color:#fff}

.bars{display:flex;flex-direction:column;gap:9px}
.bar{display:grid;grid-template-columns:minmax(130px,190px) 1fr 58px;
  gap:12px;align-items:center;font-size:13px}
.bar .nm{color:var(--ink-2);white-space:nowrap;overflow:hidden;
  text-overflow:ellipsis}
.track{height:13px;background:var(--grid);border-radius:4px;overflow:hidden}
.fill{height:100%;border-radius:4px}
.bar .pc{text-align:right;font-variant-numeric:tabular-nums;color:var(--ink-1);
  font-weight:600}
.bar.hint{cursor:help}
.bar.unranked .nm,.bar.unranked .pc{color:var(--ink-3)}
.bar.unranked .track{opacity:.6}
.subhead{font-size:12px;letter-spacing:.04em;text-transform:uppercase;
  color:var(--ink-3);margin:20px 0 6px;padding-top:16px;
  border-top:1px solid var(--grid)}
.subnote{font-size:12.5px;color:var(--ink-3);margin:0 0 14px;max-width:82ch;
  line-height:1.5}

.flags{margin:14px 0 0;padding:0;list-style:none;font-size:13px}
.flags li{padding:6px 0 6px 22px;position:relative;color:var(--ink-2)}
.flags li::before{content:"!";position:absolute;left:6px;color:var(--serious);
  font-weight:700}

details{border-top:1px solid var(--grid);padding:12px 0 2px}
details summary{cursor:pointer;font-size:13.5px;color:var(--ink-2);
  font-weight:600;list-style:none}
details summary::-webkit-details-marker{display:none}
details summary::before{content:"+ ";color:var(--ink-3)}
details[open] summary::before{content:"- "}
details .body{padding:10px 0 4px;font-size:13.5px;color:var(--ink-2);
  max-width:84ch}
dl{margin:0;display:grid;grid-template-columns:auto 1fr;gap:4px 14px;
  font-size:13px}
dt{color:var(--ink-3)}
dd{margin:0;color:var(--ink-2)}
code{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:12.5px;
  background:var(--plane);padding:1px 5px;border-radius:4px}
a{color:var(--s1)}
.tip{position:fixed;pointer-events:none;background:var(--surface-1);
  border:1px solid var(--axis);border-radius:7px;padding:8px 11px;font-size:12.5px;
  color:var(--ink-1);box-shadow:0 6px 20px rgba(0,0,0,.14);max-width:300px;
  opacity:0;transition:opacity .1s;z-index:50;line-height:1.45}
svg{display:block;width:100%;height:auto;overflow:visible}
.axl{fill:var(--ink-3);font-size:11px}
.gl{stroke:var(--grid);stroke-width:1}
.bl{stroke:var(--axis);stroke-width:1}
.themer{position:fixed;z-index:60;
  top:calc(12px + env(safe-area-inset-top, 0px));
  right:calc(16px + env(safe-area-inset-right, 0px))}
@media (max-width:640px){
  .viewer{height:300px}
  .bar{grid-template-columns:110px 1fr 48px;gap:8px;font-size:12px}
  table{font-size:12px}
  th,td{padding:6px}
}
</style>
<div class="themer"><button id="tbtn" title="Toggle light / dark">theme</button></div>
<div class="wrap">

<header>
  <h1 id="hTitle">CLEC12A binder design</h1>
  <p class="sub" id="hSub"></p>
</header>

<div class="banner">
  <h2>What this number is, and what it is not</h2>
  <p id="bannerBody"></p>
  <p><strong>No score on this page predicts clinical benefit.</strong> A candidate
  scoring at or above the reference means it matches a real antibody on the
  criteria a desktop computer can evaluate &mdash; charge, aromatic content,
  solubility, sequence plausibility. It does not mean the molecule binds.
  Establishing that requires expressing the protein, measuring affinity, testing
  on primary patient cells, then animal models, then trials. That path takes
  years and a laboratory, and most candidates that look this good still fail.</p>
  <p id="bannerDisc"></p>
</div>

<section>
  <h2>Headline</h2>
  <p class="secnote">The best candidate from this run, measured three ways. The
  percentile is the most informative: it says how the design compares to random
  members of the same library, which is the question a real screening campaign
  asks.</p>
  <div class="grid g4" id="tiles"></div>
</section>

<section>
  <h2>The target, and what already binds it</h2>
  <p class="secnote" id="structNote"></p>
  <div class="grid g2">
    <div class="card">
      <h3>CLEC12A with the mapped epitope</h3>
      <div class="viewer" id="v1"></div>
      <div class="ctrls" style="margin:12px 0 0">
        <button id="b_surface" aria-pressed="false">Surface</button>
        <button id="b_ab" aria-pressed="true">Show antibody</button>
        <button id="b_spin1" aria-pressed="false">Spin</button>
      </div>
      <div class="vlegend">
        <span><i class="sw" style="background:#9aa0a6"></i>CLEC12A</span>
        <span><i class="sw" style="background:var(--s2)"></i>epitope (13 residues)</span>
        <span><i class="sw" style="background:var(--critical)"></i>the three arginines</span>
        <span><i class="sw" style="background:var(--s1)"></i>50C1 antibody</span>
      </div>
    </div>
    <div class="card">
      <h3>The designed candidate on its scaffold</h3>
      <div class="viewer" id="v2"></div>
      <div class="ctrls" style="margin:12px 0 0">
        <button id="b_surface2" aria-pressed="false">Surface</button>
        <button id="b_spin2" aria-pressed="false">Spin</button>
        <select id="candPick"></select>
      </div>
      <div class="vlegend">
        <span><i class="sw" style="background:#9aa0a6"></i>fixed framework</span>
        <span><i class="sw" style="background:var(--s3)"></i>designed positions</span>
        <span><i class="sw" style="background:var(--critical)"></i>acidic (D/E)</span>
        <span><i class="sw" style="background:var(--s4)"></i>aromatic (F/W/Y)</span>
      </div>
    </div>
  </div>
</section>

<section>
  <h2>The search</h2>
  <p class="secnote" id="runNote"></p>
  <div class="grid g2">
    <div class="card">
      <h3>Score over generations</h3>
      <div id="chConv"></div>
    </div>
    <div class="card">
      <h3>Population diversity</h3>
      <p class="secnote" style="margin:0 0 10px">Mean pairwise difference across
      the 13 designed positions. When this collapses the search has converged;
      a flat line near zero early means it converged prematurely.</p>
      <div id="chDiv"></div>
    </div>
  </div>
</section>

<section>
  <h2>Is this better than chance?</h2>
  <p class="secnote" id="calNote"></p>
  <div class="card">
    <div id="chHist"></div>
  </div>
</section>

<section>
  <h2>Score breakdown</h2>
  <p class="secnote">Every component measures one named criterion. Hover any bar
  for what was measured and what that measurement cannot tell you. Colour marks
  the group; each bar is labelled, so colour is never the only cue.</p>
  <div class="card">
    <div id="brkHead" style="margin-bottom:16px"></div>
    <div class="bars" id="brk"></div>
    <ul class="flags" id="brkFlags"></ul>
  </div>
</section>

<section>
  <h2>Candidates</h2>
  <p class="secnote">Click a row to load it into the 3D viewer and the breakdown
  above. <code>designed</code> shows only the 13 varied positions; the framework
  is identical in every candidate.</p>
  <div class="card" style="overflow-x:auto">
    <table id="tbl"></table>
  </div>
</section>

<section>
  <h2>Method, provenance and limits</h2>
  <div class="card" id="prov"></div>
</section>

</div>
<div class="tip" id="tip"></div>

<script>
const DATA = /*__DATA__*/;
const STRUCT = /*__STRUCTURES__*/;

/* ---------- theme ---------- */
const tbtn = document.getElementById('tbtn');
tbtn.onclick = () => {
  const cur = document.documentElement.getAttribute('data-theme');
  const dark = cur ? cur === 'dark'
    : matchMedia('(prefers-color-scheme: dark)').matches;
  document.documentElement.setAttribute('data-theme', dark ? 'light' : 'dark');
  redrawAll();
};
const cssv = n => getComputedStyle(document.documentElement)
  .getPropertyValue(n).trim();

/* ---------- tooltip ---------- */
const tip = document.getElementById('tip');
function showTip(html, x, y){
  tip.innerHTML = html; tip.style.opacity = 1;
  const r = tip.getBoundingClientRect();
  let left = x + 14, top = y + 14;
  if (left + r.width > innerWidth - 8) left = x - r.width - 14;
  if (top + r.height > innerHeight - 8) top = y - r.height - 14;
  tip.style.left = Math.max(8,left)+'px'; tip.style.top = Math.max(8,top)+'px';
}
const hideTip = () => tip.style.opacity = 0;

/* ---------- header & banner ---------- */
const T = DATA.target, S = DATA.scaffold, R = DATA.reference,
      C = DATA.calibration, RUN = DATA.run;
const best = DATA.candidates[0];

const scaffoldLabel = S.name.replace(/_/g, ' ');
const article = /^[aeiou]/i.test(scaffoldLabel) ? 'an' : 'a';
document.getElementById('hTitle').textContent =
  `Designing ${article} ${scaffoldLabel} binder against ${T.gene}`;
document.getElementById('hSub').textContent =
  `${T.protein} — a surface antigen on leukaemic blasts in ${T.disease}. ` +
  `Target surface read from PDB ${T.complex_pdb}` +
  (T.resolution ? ` at ${T.resolution} Å` : '') +
  `. ${RUN.n_unique.toLocaleString()} unique sequences scored in ${RUN.elapsed_s}s.`;

document.getElementById('bannerBody').innerHTML =
  `The headline score is a <strong>design-triage composite</strong>: a weighted ` +
  `average of ${best.components.length} computable criteria. It answers &ldquo;which ` +
  `of these should a lab build first&rdquo;, which is a real and useful question. ` +
  `It is anchored to <strong>${R ? R.name : 'a reference'}</strong>, a molecule ` +
  `crystallised against this exact epitope, scored through the identical pipeline.`;

document.getElementById('bannerDisc').innerHTML =
  `<strong>Discriminating power of this scoring function: ` +
  `${C.discriminating_power}.</strong> That figure is how far the real antibody ` +
  `sits above random library members. It is the honest ceiling on how much any ` +
  `percentile below should be trusted` +
  (C.notes && C.notes.length ? ` &mdash; ${C.notes.join(' ')}` : '') + `.`;

/* ---------- tiles ---------- */
const tiles = [
  {lab:'Best candidate, recognition', val:best.recognition.toFixed(1),
   foot:`out of 100 on the interface-complementarity components`},
  {lab:'Percentile vs random library', val:best.recognition_percentile.toFixed(0)+'th',
   foot:`beats ${best.recognition_percentile.toFixed(0)}% of ${C.n_decoys} random ` +
        `members of the same ${S.variable_1based.length}-position library ` +
        `(z = ${best.recognition_z >= 0 ? '+' : ''}${best.recognition_z})`},
  {lab:'Relative to the real antibody',
   val:(best.percent_of_reference !== null ? best.percent_of_reference.toFixed(0)+'%' : 'n/a'),
   foot:`of ${R ? R.name : 'reference'}'s recognition score — on computable ` +
        `criteria only, not on measured binding`},
  // Selectivity gets its own tile because it is the axis that actually sank
  // CD33- and CD123-directed agents in the clinic. Burying it inside the
  // composite would hide the most clinically decisive weakness.
  {lab:'Selectivity', val:(best.groups.selectivity ?? 0).toFixed(0),
   foot:`out of 100. The axis that matters most here and the hardest to ` +
        `compute — a low figure is a real warning, not a rounding detail`},
  {lab:'Search space covered',
   val:((RUN.n_unique / S.search_space) * 100).toExponential(1)+'%',
   foot:`${RUN.n_unique.toLocaleString()} of ${S.search_space.toExponential(1)} ` +
        `possible sequences — the loop is a guided search, not an exhaustive one`},
];
document.getElementById('tiles').innerHTML = tiles.map(t =>
  `<div class="tile"><div class="lab">${t.lab}</div>
   <div class="val">${t.val}</div><div class="foot">${t.foot}</div></div>`).join('');

document.getElementById('structNote').innerHTML =
  `Left: the real complex. The epitope is the ${T.epitope.length} CLEC12A residues ` +
  `within 4.5 Å of the 50C1 antibody — <code>${T.epitope_sequence}</code>, net ` +
  `side-chain charge ${T.epitope_charge >= 0 ? '+' : ''}${T.epitope_charge}, ` +
  `radius of gyration ${T.epitope_radius} Å. Three arginines ` +
  `(${T.arginines.join(', ')}) dominate it, which is why the design objective is ` +
  `an acidic, aromatic-rich face. Right: the candidate on the ` +
  (S.pdb_template ? `experimental ${S.pdb_template} fold` : 'scaffold') +
  `, with only the ${S.variable_1based.length} designed positions varied.`;

document.getElementById('runNote').textContent =
  `${RUN.config.population_size} individuals per generation, ` +
  `${RUN.trajectory.length - 1} generations. Stopped because: ${RUN.stop_reason}.`;

document.getElementById('calNote').innerHTML =
  `${C.n_decoys} random draws from the same library, scored identically. This is ` +
  `the null distribution. The real antibody is marked for scale — if it did not ` +
  `stand clear of the decoys, the scoring function would be worthless and so ` +
  `would every other number here.`;

/* ---------- charts (plain SVG, built here so they re-theme) ---------- */
const NS='http://www.w3.org/2000/svg';
const el=(n,a={})=>{const e=document.createElementNS(NS,n);
  for(const k in a)e.setAttribute(k,a[k]);return e;};

function lineChart(host, series, opts){
  host.innerHTML='';
  const W=640,H=260,m={t:16,r:16,b:34,l:44};
  const iw=W-m.l-m.r, ih=H-m.t-m.b;
  const svg=el('svg',{viewBox:`0 0 ${W} ${H}`,role:'img'});
  const xs=series[0].pts.map(p=>p[0]);
  const xmin=Math.min(...xs), xmax=Math.max(...xs);
  const allY=series.flatMap(s=>s.pts.map(p=>p[1]));
  let ymin=opts.ymin!==undefined?opts.ymin:Math.min(...allY);
  let ymax=opts.ymax!==undefined?opts.ymax:Math.max(...allY);
  const pad=(ymax-ymin)*0.12||1; ymin-=pad; ymax+=pad;
  if(opts.ymin!==undefined) ymin=opts.ymin;
  const X=v=>m.l+(xmax===xmin?iw/2:(v-xmin)/(xmax-xmin)*iw);
  const Y=v=>m.t+ih-(v-ymin)/(ymax-ymin)*ih;

  for(let i=0;i<=4;i++){
    const v=ymin+(ymax-ymin)*i/4, y=Y(v);
    svg.appendChild(el('line',{x1:m.l,x2:m.l+iw,y1:y,y2:y,class:'gl'}));
    const t=el('text',{x:m.l-8,y:y+4,class:'axl','text-anchor':'end'});
    t.textContent=v.toFixed(opts.dec??0); svg.appendChild(t);
  }
  svg.appendChild(el('line',{x1:m.l,x2:m.l+iw,y1:m.t+ih,y2:m.t+ih,class:'bl'}));
  const step=Math.max(1,Math.ceil(xs.length/8));
  xs.forEach((x,i)=>{ if(i%step&&i!==xs.length-1)return;
    const t=el('text',{x:X(x),y:H-12,class:'axl','text-anchor':'middle'});
    t.textContent=x; svg.appendChild(t);});
  const xl=el('text',{x:m.l+iw/2,y:H-0,class:'axl','text-anchor':'middle'});
  xl.textContent=opts.xlab||''; svg.appendChild(xl);

  series.forEach(s=>{
    const d=s.pts.map((p,i)=>(i?'L':'M')+X(p[0])+' '+Y(p[1])).join(' ');
    svg.appendChild(el('path',{d,fill:'none',stroke:cssv(s.color),
      'stroke-width':2,'stroke-linejoin':'round','stroke-linecap':'round'}));
    const last=s.pts[s.pts.length-1];
    svg.appendChild(el('circle',{cx:X(last[0]),cy:Y(last[1]),r:4,
      fill:cssv(s.color),stroke:cssv('--surface-1'),'stroke-width':2}));
    const lb=el('text',{x:X(last[0])-6,y:Y(last[1])-11,class:'axl',
      'text-anchor':'end',fill:cssv(s.color),'font-weight':600});
    lb.textContent=s.name; svg.appendChild(lb);
  });

  // hover crosshair
  const hl=el('line',{y1:m.t,y2:m.t+ih,class:'bl',opacity:0});
  svg.appendChild(hl);
  const hit=el('rect',{x:m.l,y:m.t,width:iw,height:ih,fill:'transparent'});
  svg.appendChild(hit);
  hit.addEventListener('mousemove',ev=>{
    const bb=svg.getBoundingClientRect();
    const px=(ev.clientX-bb.left)/bb.width*W;
    let bi=0,bd=1e9;
    xs.forEach((x,i)=>{const d=Math.abs(X(x)-px); if(d<bd){bd=d;bi=i;}});
    hl.setAttribute('x1',X(xs[bi])); hl.setAttribute('x2',X(xs[bi]));
    hl.setAttribute('opacity',1);
    const rows=series.map(s=>`<div style="display:flex;gap:8px;align-items:center">
      <i class="sw" style="background:${cssv(s.color)}"></i>${s.name}:
      <b>${s.pts[bi][1].toFixed(2)}</b></div>`).join('');
    showTip(`<b>${opts.xlab||'x'} ${xs[bi]}</b>${rows}`,ev.clientX,ev.clientY);
  });
  hit.addEventListener('mouseleave',()=>{hl.setAttribute('opacity',0);hideTip();});
  host.appendChild(svg);
}

function histogram(host){
  host.innerHTML='';
  const vals=C.recognition.values;
  const W=760,H=280,m={t:20,r:16,b:48,l:48};
  const iw=W-m.l-m.r, ih=H-m.t-m.b;
  const lo=Math.min(...vals, C.reference_recognition??1e9, best.recognition)-3;
  const hi=Math.max(...vals, C.reference_recognition??-1e9, best.recognition)+3;
  const nb=26, bw=(hi-lo)/nb;
  const bins=new Array(nb).fill(0);
  vals.forEach(v=>{let i=Math.floor((v-lo)/bw); i=Math.max(0,Math.min(nb-1,i));
    bins[i]++;});
  const ymax=Math.max(...bins)||1;
  const X=v=>m.l+(v-lo)/(hi-lo)*iw;
  const Y=c=>m.t+ih-c/ymax*ih;
  const svg=el('svg',{viewBox:`0 0 ${W} ${H}`,role:'img'});

  for(let i=0;i<=4;i++){
    const c=ymax*i/4,y=Y(c);
    svg.appendChild(el('line',{x1:m.l,x2:m.l+iw,y1:y,y2:y,class:'gl'}));
    const t=el('text',{x:m.l-8,y:y+4,class:'axl','text-anchor':'end'});
    t.textContent=Math.round(c); svg.appendChild(t);
  }
  bins.forEach((c,i)=>{
    if(!c)return;
    const x=X(lo+i*bw), w=Math.max(1,(iw/nb)-2), y=Y(c);
    const r=el('rect',{x,y,width:w,height:m.t+ih-y,rx:3,
      fill:cssv('--seq-250')});
    r.addEventListener('mousemove',ev=>showTip(
      `<b>${c}</b> of ${vals.length} random library members scored ` +
      `${(lo+i*bw).toFixed(1)}&ndash;${(lo+(i+1)*bw).toFixed(1)}`,
      ev.clientX,ev.clientY));
    r.addEventListener('mouseleave',hideTip);
    svg.appendChild(r);
  });
  svg.appendChild(el('line',{x1:m.l,x2:m.l+iw,y1:m.t+ih,y2:m.t+ih,class:'bl'}));
  for(let i=0;i<=5;i++){
    const v=lo+(hi-lo)*i/5;
    const t=el('text',{x:X(v),y:H-26,class:'axl','text-anchor':'middle'});
    t.textContent=v.toFixed(0); svg.appendChild(t);
  }
  const xl=el('text',{x:m.l+iw/2,y:H-6,class:'axl','text-anchor':'middle'});
  xl.textContent='recognition score (0–100)'; svg.appendChild(xl);

  const marks=[
    {v:C.recognition.mean,c:'--ink-3',n:'random mean'},
    {v:best.recognition,c:'--s3',n:'best candidate'},
  ];
  if(C.reference_recognition!=null)
    marks.push({v:C.reference_recognition,c:'--s2',n:'real antibody'});
  marks.forEach((k,i)=>{
    const x=X(k.v);
    svg.appendChild(el('line',{x1:x,x2:x,y1:m.t-6,y2:m.t+ih,
      stroke:cssv(k.c),'stroke-width':2,'stroke-dasharray':'5 4'}));
    const t=el('text',{x:x,y:m.t-10+(i%2?-13:0),class:'axl',
      'text-anchor':'middle',fill:cssv(k.c),'font-weight':700});
    t.textContent=k.n+' '+k.v.toFixed(1); svg.appendChild(t);
  });
  host.appendChild(svg);
}

/* ---------- breakdown ---------- */
const GCOL={recognition:'--s1',selectivity:'--s2',developability:'--s3',
            plausibility:'--s4',other:'--ink-3'};
function breakdown(c){
  document.getElementById('brkHead').innerHTML =
    `<div style="font-size:13px;color:var(--ink-2)">Showing candidate #${c.rank},
     found in generation ${c.generation} by ${c.origin}. Composite
     <b style="color:var(--ink-1)">${c.composite}</b>/100 over
     ${c.n_ranked ?? c.components.length} ranking components &middot;
     recognition <b style="color:var(--ink-1)">${c.recognition}</b>/100.</div>
     <div class="seq" style="margin-top:10px">${
       [...c.annotated].map(ch=> ch===ch.toUpperCase() && /[A-Z]/.test(ch)
         ? `<span class="v">${ch}</span>` : ch).join('')}</div>
     <div style="font-size:12px;color:var(--ink-3);margin-top:6px">
       Highlighted = the ${S.variable_1based.length} designed positions
       (${S.variable_1based.join(', ')}). Lowercase = fixed framework.</div>`;

  const order={recognition:0,selectivity:1,plausibility:2,developability:3,other:4};
  const sorter=(a,b)=>(order[a.group]-order[b.group])||(b.weight-a.weight);
  const ranked=c.components.filter(k=>k.ranked!==false).sort(sorter);
  const unranked=c.components.filter(k=>k.ranked===false).sort(sorter);

  const row=k=>`
    <div class="bar hint${k.ranked===false?' unranked':''}"
         data-detail="${encodeURIComponent(k.detail)}"
         data-caveat="${encodeURIComponent(k.caveat||'')}"
         data-note="${encodeURIComponent(k.note||'')}"
         data-nm="${k.name}" data-g="${k.group}" data-w="${k.weight}">
      <div class="nm">${k.name.replace(/_/g,' ')}</div>
      <div class="track"><div class="fill" style="width:${k.normalised}%;
        background:${k.ranked===false?cssv('--ink-3')
                     :cssv(GCOL[k.group]||'--ink-3')}"></div></div>
      <div class="pc">${k.normalised.toFixed(0)}%</div>
    </div>`;

  document.getElementById('brk').innerHTML =
    ranked.map(row).join('') +
    (unranked.length ? `
      <div class="subhead">Measured, but not used for ranking</div>
      <p class="subnote">Constant across every candidate on this scaffold
      &mdash; same length, same face size, same framework &mdash; so they
      cannot separate one design from another. They are still real checks on
      the molecule, and would rank again on a scaffold where they vary.
      Leaving them in the score added a near-constant bonus to every
      candidate and compressed the useful range by about a third.</p>`
      + unranked.map(row).join('') : '');
  document.querySelectorAll('#brk .bar').forEach(b=>{
    b.addEventListener('mousemove',ev=>showTip(
      `<b>${b.dataset.nm.replace(/_/g,' ')}</b> · ${b.dataset.g} · weight ${b.dataset.w}
       <div style="margin-top:6px">${decodeURIComponent(b.dataset.detail)}</div>` +
      (b.dataset.note ? `<div style="margin-top:6px;color:var(--ink-3)">
        <b>Not ranked:</b> ${decodeURIComponent(b.dataset.note)}</div>` : '') +
      (b.dataset.caveat ? `<div style="margin-top:6px;color:var(--ink-3)">
        <b>Cannot tell you:</b> ${decodeURIComponent(b.dataset.caveat)}</div>` : ''),
      ev.clientX,ev.clientY));
    b.addEventListener('mouseleave',hideTip);
  });
  document.getElementById('brkFlags').innerHTML =
    (c.flags||[]).map(f=>`<li>${f}</li>`).join('');
}

/* ---------- table ---------- */
function table(){
  const t=document.getElementById('tbl');
  t.innerHTML=`<thead><tr>
    <th>#</th><th>designed positions</th><th>gen</th><th>origin</th>
    <th>composite</th><th>recognition</th><th>pctile</th><th>z</th>
    <th>vs antibody</th></tr></thead><tbody>` +
    DATA.candidates.map(c=>`<tr data-r="${c.rank}" class="${c.rank===cur?'sel':''}">
      <td>${c.rank}</td>
      <td style="font-family:ui-monospace,Menlo,Consolas,monospace">${c.designed}</td>
      <td>${c.generation}</td><td style="text-align:left">${c.origin}</td>
      <td>${c.composite}</td><td>${c.recognition}</td>
      <td>${c.recognition_percentile}</td>
      <td>${c.recognition_z>=0?'+':''}${c.recognition_z}</td>
      <td>${c.percent_of_reference!==null?c.percent_of_reference+'%':'—'}</td>
    </tr>`).join('') + `</tbody>`;
  t.querySelectorAll('tbody tr').forEach(tr=>
    tr.onclick=()=>select(+tr.dataset.r));
}

/* ---------- provenance ---------- */
document.getElementById('prov').innerHTML=`
  <dl>
    <dt>Target</dt><dd>${T.protein} (${T.gene}), UniProt ${T.uniprot}</dd>
    <dt>Epitope source</dt><dd>PDB ${T.complex_pdb}${T.resolution?` at ${T.resolution} Å`:''},
      chain ${T.target_chain}; residues within 4.5 Å of the 50C1 Fab</dd>
    <dt>Apo structure</dt><dd>PDB ${T.apo_pdb}</dd>
    <dt>Scaffold</dt><dd>${S.name} — ${S.description}</dd>
    <dt>Scaffold source</dt><dd>${S.provenance}</dd>
    <dt>Reference binder</dt><dd>${R?R.description:'unavailable'}</dd>
    <dt>Language model</dt><dd>${DATA.plm_status}</dd>
    <dt>Paralog surfaces</dt><dd>${Object.entries(DATA.paralogs||{})
      .map(([k,v])=>`${k} → <code>${v}</code>`).join(', ')||'none mapped'}</dd>
  </dl>
  <details><summary>Why this target was chosen</summary>
    <div class="body">${T.rationale}</div></details>
  <details><summary>The clinical reality for this target class</summary>
    <div class="body">${T.clinical_note}</div></details>
  <details><summary>What the selectivity score does and does not cover</summary>
    <div class="body">The paralog-selectivity component compares how well a
    candidate's chemistry complements the CLEC12A epitope against the aligned
    surface on its closest human relative. It does <b>not</b> predict expression
    on healthy cells or sparing of haematopoietic stem cells. Whether a target
    appears on a healthy cell is a fact about that cell measured by flow
    cytometry on primary human samples; no sequence score reaches it. The reason
    CLEC12A was selected over CD33 and CD123 is published expression data, not
    anything computed here.</div></details>
  <details><summary>Full design brief</summary>
    <div class="body"><pre style="white-space:pre-wrap;margin:0;font-size:12.5px">${
      T.brief}</pre></div></details>
  <details><summary>Search configuration</summary>
    <div class="body"><pre style="white-space:pre-wrap;margin:0;font-size:12.5px">${
      JSON.stringify(RUN.config,null,2)}</pre></div></details>`;

/* ---------- 3D viewers ---------- */
let gv1=null, gv2=null, spin1=false, spin2=false;
const ACID=new Set(['ASP','GLU']), AROM=new Set(['PHE','TRP','TYR']);
const epiNums=T.epitope.map(e=>e.resseq);
const argNums=T.epitope.filter(e=>e.aa==='R').map(e=>e.resseq);

function styleV1(){
  if(!gv1) return;
  const ch=STRUCT.complex_chains||{target:'C',heavy:'H',light:'L'};
  gv1.setStyle({},{});
  gv1.setStyle({chain:ch.target},{cartoon:{color:'#9aa0a6',opacity:1}});
  gv1.setStyle({chain:ch.target,resi:epiNums},
    {cartoon:{color:cssv('--s2')},stick:{colorscheme:'default',radius:0.22}});
  gv1.setStyle({chain:ch.target,resi:argNums},
    {cartoon:{color:cssv('--critical')},stick:{color:cssv('--critical'),radius:0.26}});
  const showAb=document.getElementById('b_ab').getAttribute('aria-pressed')==='true';
  if(showAb){
    gv1.setStyle({chain:[ch.heavy,ch.light]},
      {cartoon:{color:cssv('--s1'),opacity:0.82}});
  }
  gv1.removeAllSurfaces();
  if(document.getElementById('b_surface').getAttribute('aria-pressed')==='true'){
    gv1.addSurface($3Dmol.SurfaceType.VDW,
      {opacity:0.72,color:'#9aa0a6'},{chain:ch.target});
    gv1.addSurface($3Dmol.SurfaceType.VDW,
      {opacity:0.92,color:cssv('--s2')},{chain:ch.target,resi:epiNums});
  }
  gv1.render();
}

function styleV2(cand){
  if(!gv2) return;
  const varPos=S.variable_1based;
  gv2.setStyle({},{cartoon:{color:'#9aa0a6'}});
  const three={A:'ALA',C:'CYS',D:'ASP',E:'GLU',F:'PHE',G:'GLY',H:'HIS',I:'ILE',
    K:'LYS',L:'LEU',M:'MET',N:'ASN',P:'PRO',Q:'GLN',R:'ARG',S:'SER',T:'THR',
    V:'VAL',W:'TRP',Y:'TYR'};
  varPos.forEach((p,i)=>{
    const aa=cand.sequence[p-1];
    const col = 'DE'.includes(aa) ? cssv('--critical')
              : 'FWY'.includes(aa) ? cssv('--s4')
              : cssv('--s3');
    gv2.setStyle({resi:[p]},
      {cartoon:{color:col},stick:{color:col,radius:0.3}});
    gv2.addLabel(aa+p,{position:{resi:p},backgroundColor:col,
      backgroundOpacity:0.85,fontColor:'#fff',fontSize:10,inFront:true});
  });
  gv2.removeAllSurfaces();
  if(document.getElementById('b_surface2').getAttribute('aria-pressed')==='true'){
    gv2.addSurface($3Dmol.SurfaceType.VDW,{opacity:0.7,color:'#9aa0a6'});
  }
  gv2.render();
}

function initViewers(){
  if(typeof $3Dmol==='undefined'){
    ['v1','v2'].forEach(id=>document.getElementById(id).innerHTML=
      `<div style="padding:24px;color:var(--ink-2);font-size:13.5px">
       3D viewer could not load (no network for the CDN). Every number on this
       page is unaffected; only the structure display is missing.</div>`);
    return;
  }
  if(STRUCT.complex){
    gv1=$3Dmol.createViewer(document.getElementById('v1'),
      {backgroundColor:cssv('--plane')});
    gv1.addModel(STRUCT.complex,'pdb');
    styleV1();
    const ch=STRUCT.complex_chains||{target:'C'};
    gv1.zoomTo({chain:ch.target,resi:epiNums}); gv1.zoom(0.55); gv1.render();
  }
  if(STRUCT.scaffold){
    gv2=$3Dmol.createViewer(document.getElementById('v2'),
      {backgroundColor:cssv('--plane')});
    gv2.addModel(STRUCT.scaffold,'pdb');
    styleV2(best); gv2.zoomTo(); gv2.render();
  } else {
    document.getElementById('v2').innerHTML=
      `<div style="padding:24px;color:var(--ink-2);font-size:13.5px">
       No experimental template for this scaffold, so the fold is not shown.
       The designed residues are listed in the breakdown below.</div>`;
  }
}

/* ---------- controls ---------- */
const toggle=(id,fn)=>{
  const b=document.getElementById(id);
  b.onclick=()=>{b.setAttribute('aria-pressed',
    b.getAttribute('aria-pressed')==='true'?'false':'true'); fn();};
};
toggle('b_surface',styleV1);
toggle('b_ab',styleV1);
toggle('b_surface2',()=>styleV2(DATA.candidates[cur-1]));
toggle('b_spin1',()=>{spin1=!spin1; if(gv1) gv1.spin(spin1?'y':false);});
toggle('b_spin2',()=>{spin2=!spin2; if(gv2) gv2.spin(spin2?'y':false);});

const pick=document.getElementById('candPick');
pick.innerHTML=DATA.candidates.map(c=>
  `<option value="${c.rank}">#${c.rank} — ${c.designed} (${c.recognition})</option>`
).join('');
pick.onchange=()=>select(+pick.value);

let cur=1;
function select(rank){
  cur=rank;
  const c=DATA.candidates[rank-1];
  breakdown(c); table(); pick.value=rank;
  if(gv2){ gv2.removeAllLabels(); styleV2(c); }
}

/* ---------- draw ---------- */
function redrawAll(){
  const tr=RUN.trajectory;
  lineChart(document.getElementById('chConv'),[
    {name:'best',color:'--s1',pts:tr.map(g=>[g.generation,g.best])},
    {name:'mean',color:'--s2',pts:tr.map(g=>[g.generation,g.mean])},
  ],{xlab:'generation',dec:0});
  lineChart(document.getElementById('chDiv'),[
    {name:'diversity',color:'--s3',pts:tr.map(g=>[g.generation,g.diversity])},
  ],{xlab:'generation',dec:2,ymin:0});
  histogram(document.getElementById('chHist'));
  breakdown(DATA.candidates[cur-1]); table();
  if(gv1){gv1.setBackgroundColor(cssv('--plane')); styleV1();}
  if(gv2){gv2.setBackgroundColor(cssv('--plane'));
    gv2.removeAllLabels(); styleV2(DATA.candidates[cur-1]);}
}
redrawAll();
initViewers();
addEventListener('resize',()=>{clearTimeout(window._rt);
  window._rt=setTimeout(redrawAll,180);});
</script>
"""

def _wrap_local(fragment: str) -> str:
    """Put a document skeleton back around the fragment.

    A published artifact is wrapped in its own skeleton at publish time, so
    _TEMPLATE is written as a bare fragment. A file opened straight from disk
    needs the skeleton supplied here. The split point is the end of the
    stylesheet: everything up to it belongs in <head>.
    """
    marker = "</style>"
    idx = fragment.find(marker)
    if idx == -1:
        head, body = "", fragment
    else:
        cut = idx + len(marker)
        head, body = fragment[:cut], fragment[cut:]
    return (
        '<!doctype html>\n<html lang="en">\n<head>\n'
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1, '
        'viewport-fit=cover">\n'
        f"{head}\n</head>\n<body>\n{body}\n</body>\n</html>\n"
    )
