/* Protein Lab — front end.
 *
 * Four views, routed off the hash: the catalogue, one target, a live run, and
 * its results. State is deliberately simple — one object, one render call per
 * view — because the interesting complexity belongs in the pipeline, not here.
 */

const S = {
  catalog: null,
  target: null,
  job: null,
  poll: null,
  selected: 1,
  viewers: {},
};

const $ = (id) => document.getElementById(id);
const cssv = (n) => getComputedStyle(document.documentElement)
  .getPropertyValue(n).trim();
const esc = (s) => String(s ?? '').replace(/[&<>"]/g,
  (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

/* ---------- tooltip ---------- */
const tip = $('tip');
function showTip(html, x, y) {
  tip.innerHTML = html; tip.style.opacity = 1;
  const r = tip.getBoundingClientRect();
  let left = x + 14, top = y + 14;
  if (left + r.width > innerWidth - 8) left = x - r.width - 14;
  if (top + r.height > innerHeight - 8) top = y - r.height - 14;
  tip.style.left = Math.max(8, left) + 'px';
  tip.style.top = Math.max(8, top) + 'px';
}
const hideTip = () => { tip.style.opacity = 0; };

/* ---------- theme ---------- */
$('themeBtn').onclick = () => {
  const cur = document.documentElement.getAttribute('data-theme');
  const dark = cur ? cur === 'dark'
    : matchMedia('(prefers-color-scheme: dark)').matches;
  document.documentElement.setAttribute('data-theme', dark ? 'light' : 'dark');
  try { localStorage.setItem('plab-theme', dark ? 'light' : 'dark'); } catch (e) {}
  rerender();
};
try {
  const saved = localStorage.getItem('plab-theme');
  if (saved) document.documentElement.setAttribute('data-theme', saved);
} catch (e) {}

/* ---------- api ---------- */
async function api(path, opts) {
  const r = await fetch(path, opts);
  const body = await r.json().catch(() => ({}));
  if (!r.ok) throw new Error(body.error || `${r.status} ${r.statusText}`);
  return body;
}

/* ---------- charts ---------- */
const NS = 'http://www.w3.org/2000/svg';
const el = (n, a = {}) => {
  const e = document.createElementNS(NS, n);
  for (const k in a) e.setAttribute(k, a[k]);
  return e;
};

function lineChart(host, series, opts = {}) {
  if (!host) return;
  host.innerHTML = '';
  const pts0 = series[0] && series[0].pts;
  if (!pts0 || pts0.length < 2) {
    host.innerHTML = `<p class="secnote" style="margin:0">Waiting for the
      first generations…</p>`;
    return;
  }
  const W = 640, H = 250, m = { t: 16, r: 18, b: 34, l: 44 };
  const iw = W - m.l - m.r, ih = H - m.t - m.b;
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img' });
  const xs = pts0.map((p) => p[0]);
  const xmin = Math.min(...xs), xmax = Math.max(...xs);
  const allY = series.flatMap((s) => s.pts.map((p) => p[1]));
  let ymin = opts.ymin ?? Math.min(...allY);
  let ymax = opts.ymax ?? Math.max(...allY);
  const pad = (ymax - ymin) * 0.12 || 1;
  ymin -= pad; ymax += pad;
  if (opts.ymin !== undefined) ymin = opts.ymin;
  const X = (v) => m.l + (xmax === xmin ? iw / 2 : (v - xmin) / (xmax - xmin) * iw);
  const Y = (v) => m.t + ih - (v - ymin) / (ymax - ymin) * ih;

  for (let i = 0; i <= 4; i++) {
    const v = ymin + (ymax - ymin) * i / 4, y = Y(v);
    svg.appendChild(el('line', { x1: m.l, x2: m.l + iw, y1: y, y2: y, class: 'gl' }));
    const t = el('text', { x: m.l - 8, y: y + 4, class: 'axl', 'text-anchor': 'end' });
    t.textContent = v.toFixed(opts.dec ?? 0);
    svg.appendChild(t);
  }
  svg.appendChild(el('line', { x1: m.l, x2: m.l + iw, y1: m.t + ih, y2: m.t + ih, class: 'bl' }));
  const step = Math.max(1, Math.ceil(xs.length / 8));
  xs.forEach((x, i) => {
    if (i % step && i !== xs.length - 1) return;
    const t = el('text', { x: X(x), y: H - 12, class: 'axl', 'text-anchor': 'middle' });
    t.textContent = x; svg.appendChild(t);
  });
  const xl = el('text', { x: m.l + iw / 2, y: H, class: 'axl', 'text-anchor': 'middle' });
  xl.textContent = opts.xlab || ''; svg.appendChild(xl);

  series.forEach((s) => {
    const d = s.pts.map((p, i) => (i ? 'L' : 'M') + X(p[0]) + ' ' + Y(p[1])).join(' ');
    svg.appendChild(el('path', {
      d, fill: 'none', stroke: cssv(s.color), 'stroke-width': 2,
      'stroke-linejoin': 'round', 'stroke-linecap': 'round',
    }));
    const last = s.pts[s.pts.length - 1];
    svg.appendChild(el('circle', {
      cx: X(last[0]), cy: Y(last[1]), r: 4, fill: cssv(s.color),
      stroke: cssv('--surface-1'), 'stroke-width': 2,
    }));
    const lb = el('text', {
      x: X(last[0]) - 6, y: Y(last[1]) - 11, class: 'axl',
      'text-anchor': 'end', fill: cssv(s.color), 'font-weight': 600,
    });
    lb.textContent = s.name; svg.appendChild(lb);
  });

  const hl = el('line', { y1: m.t, y2: m.t + ih, class: 'bl', opacity: 0 });
  svg.appendChild(hl);
  const hit = el('rect', { x: m.l, y: m.t, width: iw, height: ih, fill: 'transparent' });
  svg.appendChild(hit);
  hit.addEventListener('mousemove', (ev) => {
    const bb = svg.getBoundingClientRect();
    const px = (ev.clientX - bb.left) / bb.width * W;
    let bi = 0, bd = 1e9;
    xs.forEach((x, i) => { const d = Math.abs(X(x) - px); if (d < bd) { bd = d; bi = i; } });
    hl.setAttribute('x1', X(xs[bi])); hl.setAttribute('x2', X(xs[bi]));
    hl.setAttribute('opacity', 1);
    const rows = series.map((s) => `<div style="display:flex;gap:8px;align-items:center">
      <i class="sw" style="background:${cssv(s.color)}"></i>${s.name}:
      <b>${(s.pts[bi] ? s.pts[bi][1] : 0).toFixed(2)}</b></div>`).join('');
    showTip(`<b>${opts.xlab || 'x'} ${xs[bi]}</b>${rows}`, ev.clientX, ev.clientY);
  });
  hit.addEventListener('mouseleave', () => { hl.setAttribute('opacity', 0); hideTip(); });
  host.appendChild(svg);
}

function histogram(host, cal, best) {
  if (!host || !cal) return;
  host.innerHTML = '';
  const vals = cal.recognition.values || [];
  if (!vals.length) return;
  const W = 740, H = 270, m = { t: 22, r: 16, b: 46, l: 46 };
  const iw = W - m.l - m.r, ih = H - m.t - m.b;
  const marks = [{ v: cal.recognition.mean, c: '--ink-3', n: 'random mean' }];
  if (best != null) marks.push({ v: best, c: '--s3', n: 'best candidate' });
  if (cal.reference_recognition != null)
    marks.push({ v: cal.reference_recognition, c: '--s2', n: 'known binder' });
  const lo = Math.min(...vals, ...marks.map((k) => k.v)) - 3;
  const hi = Math.max(...vals, ...marks.map((k) => k.v)) + 3;
  const nb = 26, bw = (hi - lo) / nb, bins = new Array(nb).fill(0);
  vals.forEach((v) => {
    let i = Math.floor((v - lo) / bw);
    bins[Math.max(0, Math.min(nb - 1, i))]++;
  });
  const ymax = Math.max(...bins) || 1;
  const X = (v) => m.l + (v - lo) / (hi - lo) * iw;
  const Y = (c) => m.t + ih - c / ymax * ih;
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img' });

  for (let i = 0; i <= 4; i++) {
    const c = ymax * i / 4, y = Y(c);
    svg.appendChild(el('line', { x1: m.l, x2: m.l + iw, y1: y, y2: y, class: 'gl' }));
    const t = el('text', { x: m.l - 8, y: y + 4, class: 'axl', 'text-anchor': 'end' });
    t.textContent = Math.round(c); svg.appendChild(t);
  }
  bins.forEach((c, i) => {
    if (!c) return;
    const x = X(lo + i * bw), w = Math.max(1, iw / nb - 2), y = Y(c);
    const r = el('rect', { x, y, width: w, height: m.t + ih - y, rx: 3, fill: cssv('--seq-250') });
    r.addEventListener('mousemove', (ev) => showTip(
      `<b>${c}</b> of ${vals.length} random library members scored
       ${(lo + i * bw).toFixed(1)}&ndash;${(lo + (i + 1) * bw).toFixed(1)}`,
      ev.clientX, ev.clientY));
    r.addEventListener('mouseleave', hideTip);
    svg.appendChild(r);
  });
  svg.appendChild(el('line', { x1: m.l, x2: m.l + iw, y1: m.t + ih, y2: m.t + ih, class: 'bl' }));
  for (let i = 0; i <= 5; i++) {
    const v = lo + (hi - lo) * i / 5;
    const t = el('text', { x: X(v), y: H - 26, class: 'axl', 'text-anchor': 'middle' });
    t.textContent = v.toFixed(0); svg.appendChild(t);
  }
  const xl = el('text', { x: m.l + iw / 2, y: H - 6, class: 'axl', 'text-anchor': 'middle' });
  xl.textContent = 'recognition score (0–100)'; svg.appendChild(xl);
  marks.forEach((k, i) => {
    const x = X(k.v);
    svg.appendChild(el('line', {
      x1: x, x2: x, y1: m.t - 6, y2: m.t + ih, stroke: cssv(k.c),
      'stroke-width': 2, 'stroke-dasharray': '5 4',
    }));
    const t = el('text', {
      x, y: m.t - 10 + (i % 2 ? -13 : 0), class: 'axl', 'text-anchor': 'middle',
      fill: cssv(k.c), 'font-weight': 700,
    });
    t.textContent = `${k.n} ${k.v.toFixed(1)}`; svg.appendChild(t);
  });
  host.appendChild(svg);
}

/* ---------- 3D ---------- */
async function loadStructure(pdbId, chains) {
  const q = chains && chains.length ? `?chains=${chains.join(',')}` : '';
  const r = await fetch(`/api/structure/${pdbId}${q}`);
  if (!r.ok) throw new Error(`structure ${pdbId} unavailable`);
  return r.text();
}

function viewerUnavailable(host, msg) {
  host.innerHTML = `<div style="padding:22px;color:var(--ink-2);font-size:13px">
    ${esc(msg)}</div>`;
}

async function showComplex(hostId, t) {
  const host = $(hostId);
  if (!host) return;
  if (typeof $3Dmol === 'undefined')
    return viewerUnavailable(host, '3D viewer could not load (no network for the CDN). Every number on this page is unaffected.');
  try {
    const chains = [t.target_chain, ...(t.binder_chains || [])].filter(Boolean);
    const pdb = await loadStructure(t.complex_pdb, chains);
    const v = $3Dmol.createViewer(host, { backgroundColor: cssv('--plane') });
    v.addModel(pdb, 'pdb');
    S.viewers[hostId] = { v, t };
    styleComplex(hostId);
    const nums = t.epitope.map((e) => e.resseq);
    v.zoomTo({ chain: t.target_chain, resi: nums });
    v.zoom(0.6); v.render();
  } catch (e) {
    viewerUnavailable(host, `Could not render the structure: ${e.message}`);
  }
}

function styleComplex(hostId) {
  const entry = S.viewers[hostId];
  if (!entry) return;
  const { v, t } = entry;
  const nums = t.epitope.map((e) => e.resseq);
  const cationNums = t.epitope.filter((e) => 'KR'.includes(e.aa)).map((e) => e.resseq);
  const showBinder = $('b_binder') ? $('b_binder').getAttribute('aria-pressed') === 'true' : true;
  const surf = $('b_surface') ? $('b_surface').getAttribute('aria-pressed') === 'true' : false;

  v.setStyle({}, {});
  v.setStyle({ chain: t.target_chain }, { cartoon: { color: '#9aa0a6' } });
  v.setStyle({ chain: t.target_chain, resi: nums },
    { cartoon: { color: cssv('--s2') }, stick: { colorscheme: 'default', radius: 0.22 } });
  if (cationNums.length)
    v.setStyle({ chain: t.target_chain, resi: cationNums },
      { cartoon: { color: cssv('--critical') }, stick: { color: cssv('--critical'), radius: 0.26 } });
  if (showBinder && (t.binder_chains || []).length)
    v.setStyle({ chain: t.binder_chains }, { cartoon: { color: cssv('--s1'), opacity: 0.8 } });
  v.removeAllSurfaces();
  if (surf) {
    v.addSurface($3Dmol.SurfaceType.VDW, { opacity: 0.72, color: '#9aa0a6' },
      { chain: t.target_chain });
    v.addSurface($3Dmol.SurfaceType.VDW, { opacity: 0.92, color: cssv('--s2') },
      { chain: t.target_chain, resi: nums });
  }
  v.render();
}

async function showDesigned(hostId, scaffold, cand) {
  const host = $(hostId);
  if (!host) return;
  if (typeof $3Dmol === 'undefined' || !scaffold.pdb_template)
    return viewerUnavailable(host,
      'No experimental template for this scaffold, so the fold is not shown. The designed residues are listed below.');
  try {
    const pdb = await loadStructure(scaffold.pdb_template, ['A']);
    const v = $3Dmol.createViewer(host, { backgroundColor: cssv('--plane') });
    v.addModel(pdb, 'pdb');
    S.viewers[hostId] = { v, scaffold };
    styleDesigned(hostId, cand);
    v.zoomTo(); v.render();
  } catch (e) {
    viewerUnavailable(host, `Could not render the scaffold: ${e.message}`);
  }
}

function styleDesigned(hostId, cand) {
  const entry = S.viewers[hostId];
  if (!entry) return;
  const { v, scaffold } = entry;
  v.removeAllLabels();
  v.setStyle({}, { cartoon: { color: '#9aa0a6' } });
  scaffold.variable_1based.forEach((p) => {
    const aa = cand.sequence[p - 1];
    const col = 'DE'.includes(aa) ? cssv('--critical')
      : 'FWY'.includes(aa) ? cssv('--s4')
        : 'KR'.includes(aa) ? cssv('--s1') : cssv('--s3');
    v.setStyle({ resi: [p] }, { cartoon: { color: col }, stick: { color: col, radius: 0.3 } });
    v.addLabel(aa + p, {
      position: { resi: p }, backgroundColor: col, backgroundOpacity: 0.85,
      fontColor: '#fff', fontSize: 10, inFront: true,
    });
  });
  v.removeAllSurfaces();
  if ($('b_dsurface') && $('b_dsurface').getAttribute('aria-pressed') === 'true')
    v.addSurface($3Dmol.SurfaceType.VDW, { opacity: 0.7, color: '#9aa0a6' });
  v.render();
}

function toggle(id, fn) {
  const b = $(id);
  if (!b) return;
  b.onclick = () => {
    b.setAttribute('aria-pressed', b.getAttribute('aria-pressed') === 'true' ? 'false' : 'true');
    fn();
  };
}

/* ---------- banner ---------- */
const BANNER = `
<div class="banner">
  <h2>What this does, and what it cannot do</h2>
  <p>This runs a real <strong>design triage</strong> pipeline: it reads a
  binding site off an experimentally solved structure, generates candidate
  binder sequences on a validated protein scaffold, and scores them on
  published biophysical criteria plus a protein language model. It answers
  &ldquo;which of these should a lab build first&rdquo;, which is a real and
  useful question.</p>
  <p><strong>No score here predicts clinical benefit.</strong> A candidate
  matching the known binder means it matches on criteria a computer can
  evaluate &mdash; charge, aromatic content, solubility, sequence
  plausibility. It does not mean the molecule binds. Establishing that needs
  expression, binding assays, cell work, animal models and trials: years, in a
  laboratory. Most candidates that look this good still fail.</p>
  <p><strong>Not every illness can be worked on here.</strong> The method
  needs a protein target with a solved structure and something already bound
  to it. Diseases with no single protein driver have nothing to aim at. Every
  target offered below has been verified against real structural data; the
  ones that failed verification are listed with the reason.</p>
</div>`;

/* ---------- views ---------- */
function rerender() {
  const h = location.hash.slice(1);
  if (h.startsWith('target/')) {
    const rest = h.slice(7);
    return rest.endsWith('/run')
      ? renderTarget(rest.slice(0, -4), true)
      : renderTarget(rest);
  }
  if (h.startsWith('run/')) return renderRun(h.slice(4));
  if (h.startsWith('build')) return renderBuild();
  if (h.startsWith('specimen')) return renderSpecimen();
  return renderCatalog();
}

async function renderCatalog() {
  const v = $('view');
  if (!S.catalog) {
    v.innerHTML = `<p class="secnote"><span class="spin"></span>Loading the
      target catalogue…</p>`;
    S.catalog = await api('/api/catalog');
    renderRail();
  }
  const c = S.catalog;
  const q = ($('q').value || '').toLowerCase().trim();

  const match = (t) => !q || [t.gene, t.disease, t.protein_name, t.area,
    t.uniprot, t.pdb].join(' ').toLowerCase().includes(q);

  let areas = '';
  let shown = 0;
  for (const [area, list] of Object.entries(c.areas)) {
    const hits = list.filter(match);
    if (!hits.length) continue;
    shown += hits.length;
    areas += `<section id="area-${esc(area).replace(/\W+/g, '-')}">
      <h2>${esc(area)}</h2>
      <p class="secnote">${hits.length} verified target${hits.length > 1 ? 's' : ''}</p>
      <div class="grid g3">${hits.map(targetCard).join('')}</div>
    </section>`;
  }
  if (!shown) areas = `<p class="secnote">Nothing matches &ldquo;${esc(q)}&rdquo;.</p>`;

  v.innerHTML = `
    <h1>Design a protein binder against a disease target</h1>
    <p class="secnote">${c.n_verified} targets across
      ${Object.keys(c.areas).length} disease areas, each verified against an
      experimental structure. Pick one to see its binding site and run a
      design cycle — or
      <a href="#build">point the same pipeline at any protein you like</a>.</p>
    ${BANNER}
    ${areas}
    ${rejectedSection(c)}`;
  $('q').oninput = renderCatalog;
}

function targetCard(t) {
  const s = t.surface || {};
  return `<button class="tcard" onclick="location.hash='target/${t.key}'">
    <div class="gene">${esc(t.gene)}</div>
    <div class="dis">${esc(t.disease)}</div>
    <div class="meta">
      <span class="pill epi">${esc(t.epitope_sequence.slice(0, 14))}${t.epitope_sequence.length > 14 ? '…' : ''}</span>
      <span class="pill">${t.n_epitope_residues} residues</span>
      <span class="pill s1">PDB ${esc(t.pdb)}</span>
      ${s.required_character ? `<span class="pill">wants ${esc(s.required_character)}</span>` : ''}
      ${t.has_paralog_panel
      ? '<span class="pill ok">selectivity panel</span>'
      // "No close relatives" is a result; "no panel" is a missing result.
      // One pill for both would hide that difference.
      : t.panel_state === 'no_human_relatives'
        ? '<span class="pill">no human relatives</span>'
        : '<span class="pill warn">selectivity not assessed</span>'}
    </div>
  </button>`;
}

function rejectedSection(c) {
  if (!c.rejected || !c.rejected.length) return '';
  return `<section>
    <h2>Candidates that failed verification</h2>
    <p class="secnote">These were considered and rejected. They are shown
    because a tool that silently drops what it cannot do is harder to trust
    than one that says so. Worth knowing: an earlier build of this catalogue
    rejected five more than it does now, and every one of those turned out to
    be this pipeline's mistake rather than a limit of the data — three wrong
    PDB ids and a sequence-matching rule too brittle to recognise a protein
    from a different virus strain. A rejection is a claim about the evidence,
    so it has to be checked as carefully as a result.</p>
    <div class="card"><table><thead><tr>
      <th>Target</th><th>Gene</th><th>PDB</th><th>Reason</th></tr></thead>
      <tbody>${c.rejected.map((r) => `<tr style="cursor:default">
        <td>${esc(r.disease)}</td><td class="mono">${esc(r.gene)}</td>
        <td class="mono">${esc(r.pdb)}</td>
        <td style="text-align:left;color:var(--ink-3)">${esc(r.reason.split(':').slice(1).join(':').trim() || r.reason)}</td>
      </tr>`).join('')}</tbody></table></div>
  </section>`;
}

function renderRail() {
  const c = S.catalog;
  if (!c) return;
  $('railTag').textContent =
    `${c.n_verified} verified targets · ${Object.keys(c.areas).length} areas`;
  $('areaList').innerHTML = Object.entries(c.areas).map(([a, l]) =>
    `<div style="font-size:13px;padding:4px 0">
      <a href="#" onclick="document.getElementById('q').value='';location.hash='';
        setTimeout(()=>document.getElementById('area-${esc(a).replace(/\W+/g, '-')}')
        ?.scrollIntoView({behavior:'smooth'}),60);return false">
        ${esc(a)}</a>
      <span style="color:var(--ink-3)"> ${l.length}</span>
    </div>`).join('');
}

// Why an absent panel is absent. These are different claims about how much
// is known and must not render alike: a viral protein with no human
// relatives has been checked and found to have none, while an unreadable
// numbering scheme means nothing was checked at all.
const PANEL_ABSENT = {
  no_human_relatives: [
    'No close human relative exists',
    `No other reviewed human protein shares a domain family with this target,
     so cross-reactivity with a close relative is not the risk here. That is a
     finding, not a gap. Binding to unrelated human proteins is still a real
     risk, and this method cannot assess it.`],
  no_mappable_relatives: [
    'Relatives exist but cannot be compared',
    `Human proteins share a domain family with this target, but none could
     have the binding site mapped onto it closely enough to compare. They
     diverge at this particular surface — mildly reassuring, but not a
     measurement.`],
  panel_too_weak: [
    'Only one distant relative — not enough to score',
    `A single relative was found, and it is only distantly similar at this
     binding site. Tested against controls, a margin measured against one
     such protein does not reliably separate a real binder from a face
     designed for the wrong protein, so no number is reported.`],
  no_close_relatives: [
    'Relatives share a fold, not this surface',
    `Human proteins share a domain family with this target, but none is more
     than 25% identical to it at the binding site itself. A surrogate surface
     built from an alignment that weak is not reliable enough to score
     against, so no number is reported rather than a misleading one.`],
  family_too_generic: [
    'The only domain family is too broad to be useful',
    `This protein's sole domain family on record is shared by several hundred
     human proteins — an immunoglobulin fold, or similar. That identifies a
     shape, not a set of close relatives, and a panel drawn from it would be
     arbitrary. This is a limit of the method, not a finding about the
     target.`],
  unmappable: [
    'Binding site cannot be located in the reference sequence',
    `The structure's residue numbering could not be reconciled with the
     reference sequence for this protein, so the site cannot be found in any
     relative. Nothing was assessed.`],
  no_family: [
    'No domain family on record',
    `UniProt lists no Pfam domain family for this protein, so there is no
     principled set of relatives to compare it against.`],
  none: [
    'No panel was built',
    `Cross-reactivity was not assessed for this target. This is a gap in the
     evidence, never a finding about a candidate.`],
};

// Selectivity has no number — say which of two very different reasons.
// Either nothing was compared, or a comparison was made and found that no
// binder could tell these proteins apart. The second is a finding about the
// targets and is far more informative than silence.
function selectivityAbsentTile(r, best) {
  const comp = (best.components || [])
    .find((c) => c.name === 'paralog_selectivity');
  if (comp && !comp.ranked && (r.off_targets || []).length) {
    // Two reasons, two labels. A ceiling of zero means no binder could tell
    // these proteins apart — a finding about the proteins. A ceiling no
    // bigger than its own error bar means this measurement ran out of
    // resolution — a finding about the method. VEGF-A is the second, and
    // calling it "indistinguishable" would say its nearest relative looks
    // like it, when the mapped surfaces are 4/17 identical.
    const noise = r.panel_ceiling_noise, ceil = r.panel_ceiling;
    const belowResolution = noise != null && ceil != null
      && ceil > 0 && ceil <= 2 * noise;
    return tile('Selectivity',
      belowResolution ? 'below resolution' : 'indistinguishable',
      comp.detail);
  }
  return tile('Selectivity', 'not assessed',
    (PANEL_ABSENT[r.panel_state] || PANEL_ABSENT.none)[1]);
}

// Selectivity has a number. The number alone is not the whole finding: the
// margin it comes from has a sign, and below zero the candidate complements
// the nearest relative BETTER than the target. That is not weak selectivity,
// it is selectivity aimed at the wrong protein of the two, and the old ramp
// hid it — anchored to the 10th percentile of random faces, a TNF-alpha
// design with a margin of -0.153 was reported as 65.6% selective.
function selectivityTile(r, best) {
  const comp = (best.components || [])
    .find((c) => c.name === 'paralog_selectivity');
  const against = `margin over the closest of ${r.off_targets.length} related
    human proteins. Scaled to what this target allows, so it is not
    comparable with another target's number`;
  if (comp && comp.raw < 0) {
    return tile('Selectivity', best.groups.selectivity.toFixed(0),
      `<b style="color:var(--serious)">Counter-selective.</b> The margin is
       ${comp.raw.toFixed(3)} — below zero, so this face complements the
       nearest relative better than it complements ${esc(r.gene || 'the target')}.
       Scored in the bottom quarter for that reason, not merely low.`);
  }
  return tile('Selectivity', best.groups.selectivity.toFixed(0), against);
}

// The panel a candidate's selectivity is scored against. Shown in full,
// because "selectivity 90" means nothing without knowing what it was 90
// against — and because a panel is a list of sequences, not an assay.
function panelSection(t) {
  const members = t.panel || [];
  if (!members.length) {
    const [head, body] = PANEL_ABSENT[t.panel_state] || PANEL_ABSENT.none;
    const rej = (t.panel_rejected || [])[0];
    return `
      <section>
        <h2>Cross-reactivity panel</h2>
        <div class="card" style="border-left:3px solid var(--serious)">
          <h3>${esc(head)}</h3>
          <p style="font-size:13.5px;color:var(--ink-2);margin:6px 0 0">${body}</p>
          ${rej ? `<p class="subnote" style="margin-top:10px">The relative
            found was <span class="mono">${esc(rej.gene)}</span>, ${
            (rej.window_identity * 100).toFixed(0)}% identical at this site.
            It is named because it is still the most likely place for a
            binder to cross-react, even though the comparison is too weak to
            put a number on.</p>` : ''}
          ${t.panel_note && !rej ? `<p class="subnote" style="margin-top:10px">${esc(t.panel_note)}</p>` : ''}
        </div>
      </section>`;
  }
  const derived = (t.panel_state || '').includes('derived');
  const strength = t.panel_strength;
  const [weak, blurb] = {
    close: [false, `The nearest relative is highly similar at this surface, so
      telling them apart is a genuine and demanding test.`],
    moderate: [false, `The nearest relative is moderately similar at this
      surface. The comparison is informative, but the surrogate surface is a
      model rather than a structure.`],
    distant: [true, `The nearest relative is only distantly similar at this
      surface. These proteins share a fold more than they share this patch,
      so a high score here is a lower bar than it looks.`],
  }[strength] || [false, ''];
  return `
    <section>
      <h2>Cross-reactivity panel</h2>
      <p class="secnote">The proteins a binder for this target must
      <em>not</em> engage. A candidate is scored on the margin between how
      well it complements this binding site and how well it complements the
      equivalent surface on its closest relative here.</p>
      <div class="card">
        <table><thead><tr><th>Protein</th><th>UniProt</th>
          <th>Identity at this site</th>
          <th style="text-align:left">Why it is in the panel</th></tr></thead>
        <tbody>${members.map((m) => `<tr style="cursor:default">
          <td class="mono">${esc(m.gene)}</td>
          <td class="mono">${esc(m.uniprot)}</td>
          <td>${m.window_identity != null
            ? `${(m.window_identity * 100).toFixed(0)}%`
            : '<span style="color:var(--ink-3)">—</span>'}</td>
          <td style="text-align:left;color:var(--ink-3)">${esc(m.concern)}${
            m.structural ? '' : ' <b>Different fold family — compared for context, not scored.</b>'}</td>
        </tr>`).join('')}</tbody></table>
        ${blurb ? `<p class="subnote" style="margin-top:12px">
          <b style="color:var(--${weak ? 'serious' : 'ink-1'})">${
            esc(strength.replace(/^./, (c) => c.toUpperCase()))} comparison${
            t.panel_closest_identity != null
              ? ` — closest relative ${(t.panel_closest_identity * 100).toFixed(0)}% identical here`
              : ''}.</b> ${blurb}</p>` : ''}
        <p class="subnote" style="margin-top:8px">
          ${derived ? `Derived from ${esc(t.panel_note || 'a domain-family search')}`
                    : 'Hand-curated for this target.'}
        </p>
        ${referenceCheck(t)}
        <p class="subnote" style="margin-top:8px">
          <b style="color:var(--serious)">What this is not.</b> A panel
          compares surfaces on composition at aligned positions. It is not a
          cross-reactivity assay, and it says nothing about which cells
          express these proteins — the fact that actually decides whether an
          off-target matters. That is measured by flow cytometry on primary
          human samples.
        </p>
      </div>
    </section>`;
}

// Does this metric agree with a molecule known to work?
//
// The strongest calibration available for a selectivity panel: run the real
// drug's own contact residues through it. A positive margin means the measure
// and the molecule agree. A negative one means the measure claims a licensed
// antibody prefers the off-target — and whatever a design then scores against
// the same panel is a number from a measure that is wrong about the one case
// where the answer is already known.
//
// This lived in a diagnostic script only I ever ran. The dashboard showed a
// confident percentage for targets where the reference came out negative, and
// said nothing about it.
function referenceCheck(t) {
  const r = t.reference_check;
  if (!r || !r.binder) return '';
  const bad = r.counter_selective;
  return `
    <p class="subnote" style="margin-top:12px;padding-top:10px;
       border-top:1px solid var(--rule)">
      <b style="color:var(--${bad ? 'serious' : 'ink-1'})">Does this measure
      agree with a molecule that works?</b>
      ${esc(r.binder)}'s own contact residues score
      <span class="mono">${r.score.toFixed(0)}</span> against this panel,
      margin <span class="mono">${r.margin > 0 ? '+' : ''}${r.margin.toFixed(3)}</span>${
        bad
          ? `. That is <b>negative</b>: this measure says ${esc(r.binder)}
             complements <span class="mono">${esc(r.closest)}</span> better
             than it complements ${esc(t.gene)}. It is wrong about the one
             molecule here whose selectivity is established, so treat a
             design's score against this panel as weak evidence — the
             composition of a surface is evidently not what makes a binder
             selective for this target.`
          : `, so it prefers ${esc(t.gene)} to
             <span class="mono">${esc(r.closest)}</span> — the measure and the
             molecule agree on the case where the answer is already known.`}
    </p>`;
}

// What verification noticed about this structure without rejecting it. Shown
// on the target page rather than buried, because each item limits what the
// numbers from a run can mean -- a 3.5 A structure gives an approximate
// binding site however confident the score looks.
function caveatSection(t) {
  const items = t.caveats || [];
  const variant = t.target_identity != null && t.target_identity < 0.995;
  if (!items.length && !variant) return '';
  const idLine = variant
    ? `<li>The chain used is ${(t.target_identity * 100).toFixed(0)}% identical
       to the reference sequence for this protein, matched by
       ${esc(t.identity_method || 'alignment')}. Engineered constructs and
       strain or species variants are normal in crystallography, but a binder
       shaped for this copy may fit the human protein less well.</li>`
    : '';
  return `
    <section>
      <h2>Limits of this binding site</h2>
      <p class="secnote">Recorded by the catalogue build, which checked this
      structure against the sequence database before the target was offered.
      None of these was severe enough to reject the entry; all of them bound
      how far a score can be trusted.</p>
      <div class="card" style="border-left:3px solid var(--warning)">
        <ul style="font-size:13.5px;color:var(--ink-2);margin:4px 0;padding-left:20px;line-height:1.65">
          ${idLine}
          ${items.map(w => `<li>${esc(w.charAt(0).toUpperCase() + w.slice(1))}.</li>`).join('')}
        </ul>
      </div>
    </section>`;
}

async function renderTarget(key, auto = false) {
  const v = $('view');
  v.innerHTML = `<p class="secnote"><span class="spin"></span>Loading target…</p>`;
  let t;
  try { t = await api(`/api/target/${key}`); }
  catch (e) { v.innerHTML = `<p class="secnote">${esc(e.message)}</p>`; return; }
  S.target = t;
  const s = t.surface || {};

  v.innerHTML = `
    <div class="crumb">${t.custom
        ? '<a href="#build">Build your own</a>'
        : '<a href="#">All targets</a>'} <span>›</span>
      <span>${esc(t.area)}</span></div>
    <h1>${esc(t.gene)} — ${esc(t.disease)}</h1>
    ${t.custom ? `<p class="secnote"><span class="pill warn">Built by you,
      not reviewed</span> Derived automatically from
      <span class="mono">${esc(t.complex_pdb)}</span>. It is not in the
      catalogue and will be gone when the server restarts.</p>` : ''}
    <p class="secnote">${esc(t.protein_name)} · UniProt
      <code>${esc(t.uniprot)}</code> · structure <code>${esc(t.complex_pdb)}</code>
      ${t.resolution ? `at ${t.resolution} Å` : ''}</p>

    <div class="grid g2">
      <div class="card">
        <h3>The binding site, read from the structure</h3>
        <div class="viewer" id="v_complex"></div>
        <div class="ctrls">
          <button id="b_surface" aria-pressed="false">Surface</button>
          <button id="b_binder" aria-pressed="true">Show bound partner</button>
          <button id="b_spin" aria-pressed="false">Spin</button>
        </div>
        <div class="vlegend">
          <span><i class="sw" style="background:#9aa0a6"></i>${esc(t.gene)}</span>
          <span><i class="sw" style="background:var(--s2)"></i>binding site (${t.epitope.length})</span>
          ${t.cations.length ? '<span><i class="sw" style="background:var(--critical)"></i>basic residues</span>' : ''}
          <span><i class="sw" style="background:var(--s1)"></i>${esc(t.known_binder)}</span>
        </div>
      </div>

      <div class="card">
        <h3>What a binder has to do here</h3>
        <p style="font-size:13.5px;color:var(--ink-2);margin:6px 0 12px">
          ${esc(s.rationale || 'Surface character could not be summarised.')}</p>
        <dl>
          <dt>Binding site</dt><dd class="mono">${esc(t.epitope_sequence)}</dd>
          <dt>Residues</dt><dd>${t.epitope.length}, radius of gyration ${t.epitope_radius} Å</dd>
          <dt>Net charge</dt><dd>${t.epitope_charge >= 0 ? '+' : ''}${t.epitope_charge} at pH 7.4</dd>
          <dt>Composition</dt><dd>${s.cations ?? 0} basic, ${s.anions ?? 0} acidic,
            ${s.aromatic ?? 0} aromatic, ${s.hydrophobic ?? 0} hydrophobic</dd>
          <dt>Required</dt><dd><b style="color:var(--ink-1)">${esc(s.required_character || '—')}</b> binding face</dd>
          <dt>Known binder</dt><dd>${esc(t.known_binder)}, contacts
            <span class="mono">${esc(t.paratope_sequence)}</span></dd>
          <dt>Contacts</dt><dd>${(t.descriptors || {}).n_residue_contacts ?? '—'} residue pairs,
            ${(t.descriptors || {}).n_salt_bridges ?? 0} salt bridges</dd>
        </dl>
      </div>
    </div>

    ${panelSection(t)}

    ${t.custom ? `
    <section>
      <h2>Why this target, and the honest caveat</h2>
      <div class="card" style="border-left:3px solid var(--serious)">
        <h3>Nobody wrote one</h3>
        <p style="font-size:13.5px;color:var(--ink-2);margin:6px 0 0">
          You built this target, so it has no rationale and no clinical
          caveat. Every catalogue target carries both: why the protein is
          worth engaging, and what has already gone wrong for drugs against
          it. Those are judgements a person makes after reading the
          literature, and their absence here is a gap in the evidence rather
          than an absence of concerns.</p>
        <p style="font-size:13.5px;color:var(--ink-2);margin:8px 0 0">
          Everything else on this page — the binding site, the surface
          chemistry, the cross-reactivity panel, the scores — was produced by
          exactly the same code that produced them for the catalogue, and the
          automatic checks below flagged whatever they flagged.</p>
      </div>
    </section>` : `
    <section>
      <h2>Why this target, and the honest caveat</h2>
      <div class="grid g2">
        <div class="card"><h3>Rationale</h3>
          <p style="font-size:13.5px;color:var(--ink-2);margin:6px 0 0">${esc(t.rationale)}</p></div>
        <div class="card" style="border-left:3px solid var(--serious)">
          <h3>What has gone wrong with this target class</h3>
          <p style="font-size:13.5px;color:var(--ink-2);margin:6px 0 0">${esc(t.clinical_note)}</p></div>
      </div>
    </section>`}

    ${caveatSection(t)}

    <section>
      <h2>Run a design cycle</h2>
      <p class="secnote">An evolutionary search over the designed positions of
      a validated scaffold, scored every generation and calibrated against
      random members of the same library.</p>
      <div class="card">
        <div class="grid g4" style="margin-bottom:14px">
          <div><label for="f_gen">Generations</label>
            <input id="f_gen" type="number" min="4" max="60" value="20" style="width:100%"></div>
          <div><label for="f_pop">Population</label>
            <input id="f_pop" type="number" min="8" max="96" value="40" style="width:100%"></div>
          <div><label for="f_dec">Calibration decoys</label>
            <input id="f_dec" type="number" min="40" max="400" value="150" style="width:100%"></div>
          <div><label for="f_seed">Seed</label>
            <input id="f_seed" type="number" value="7" style="width:100%"></div>
        </div>
        <div class="row">
          <label style="display:flex;gap:8px;align-items:center;cursor:pointer">
            <input type="checkbox" id="f_plm" style="width:auto">
            Use ESM-2 protein language model (GPU; much slower, first run
            downloads ~2.5 GB)
          </label>
          <div class="spacer"></div>
          <button class="primary" id="goBtn">Run design cycle</button>
        </div>
        <p class="subnote" style="margin:12px 0 0">Without the language model
        a run takes a few seconds. With it, expect several minutes — every
        candidate needs one masked forward pass per residue.</p>
      </div>
    </section>

    <section>
      <h2>Provenance</h2>
      <div class="card">
        <dl>
          <dt>Structure</dt><dd>${esc(t.structure_title || t.complex_pdb)}</dd>
          <dt>Target chain</dt><dd>${esc(t.target_chain)}</dd>
          <dt>Partner chains</dt><dd>${esc((t.binder_chains || []).join(', ') || '—')}</dd>
          <dt>Sequence</dt><dd>${t.sequence_length || '—'} residues (UniProt ${esc(t.uniprot)})</dd>
        </dl>
        <details><summary>Full design brief</summary>
          <div class="body"><pre>${esc(t.brief)}</pre></div></details>
      </div>
    </section>`;

  showComplex('v_complex', t);
  toggle('b_surface', () => styleComplex('v_complex'));
  toggle('b_binder', () => styleComplex('v_complex'));
  toggle('b_spin', () => {
    const e = S.viewers['v_complex'];
    if (e) e.v.spin($('b_spin').getAttribute('aria-pressed') === 'true' ? 'y' : false);
  });
  $('goBtn').onclick = startRun;
  // #target/<key>/run — reached from the specimen tab, so picking a target
  // there and starting a design against it is one click rather than four.
  // The page still renders in full first: the run replaces the view a moment
  // later, and arriving straight at a progress bar with no sight of what is
  // being designed against would be the wrong trade.
  if (auto) startRun();
}

async function startRun() {
  const btn = $('goBtn');
  btn.disabled = true; btn.textContent = 'Starting…';
  try {
    const body = {
      target: S.target.key,
      generations: +$('f_gen').value,
      population: +$('f_pop').value,
      decoys: +$('f_dec').value,
      seed: +$('f_seed').value,
      use_plm: $('f_plm').checked,
    };
    const r = await api('/api/run', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body),
    });
    location.hash = `run/${r.job}`;
  } catch (e) {
    btn.disabled = false; btn.textContent = 'Run design cycle';
    alertInline(e.message);
  }
}

function alertInline(msg) {
  const p = document.createElement('p');
  p.className = 'subnote';
  p.style.color = 'var(--critical)';
  p.textContent = msg;
  $('goBtn').parentElement.appendChild(p);
}

async function renderRun(jobId) {
  if (S.poll) { clearInterval(S.poll); S.poll = null; }
  const v = $('view');
  v.innerHTML = `<p class="secnote"><span class="spin"></span>Starting…</p>`;
  let first = true;

  const tick = async () => {
    let job;
    try { job = await api(`/api/job/${jobId}`); }
    catch (e) {
      clearInterval(S.poll);
      v.innerHTML = `<p class="secnote">${esc(e.message)}</p>`;
      return;
    }
    S.job = job;
    if (first) { first = false; addRunToRail(job); }
    if (job.state === 'done') {
      clearInterval(S.poll); S.poll = null;
      renderResults(job);
    } else if (job.state === 'failed') {
      clearInterval(S.poll); S.poll = null;
      v.innerHTML = `
        <div class="crumb"><a href="#">All targets</a></div>
        <h1>The run failed</h1>
        <div class="card" style="border-left:3px solid var(--critical)">
          <pre>${esc(job.error || 'unknown error')}</pre>
        </div>`;
    } else {
      renderProgress(job);
    }
  };
  await tick();
  S.poll = setInterval(tick, 1200);
}

function renderProgress(job) {
  const v = $('view');
  const traj = job.trajectory || [];
  if (!$('progBar')) {
    v.innerHTML = `
      <div class="crumb"><a href="#">All targets</a> <span>›</span>
        <span>${esc(job.params.target)}</span></div>
      <h1>Designing…</h1>
      <p class="secnote">Population ${job.params.population},
        up to ${job.params.generations} generations,
        ${job.params.use_plm ? 'with' : 'without'} the language model.</p>
      <div class="card">
        <div class="prog"><div id="progBar" style="width:2%"></div></div>
        <div class="stage" id="stageTxt"></div>
      </div>
      <section><div class="card">
        <h3>Score over generations</h3>
        <div id="liveChart"></div>
      </div></section>
      <section id="calSec"></section>`;
  }
  $('progBar').style.width = `${Math.max(2, job.progress * 100)}%`;
  $('stageTxt').innerHTML = `<span class="spin"></span>${esc(job.stage)}`;
  if (traj.length > 1) {
    lineChart($('liveChart'), [
      { name: 'best', color: '--s1', pts: traj.map((g) => [g.generation, g.best]) },
      { name: 'mean', color: '--s2', pts: traj.map((g) => [g.generation, g.mean]) },
    ], { xlab: 'generation', dec: 0 });
  }
  if (job.calibration && !$('calBody')) {
    $('calSec').innerHTML = `
      <h2>Is this better than chance?</h2>
      <p class="secnote">Random members of the same library, scored
      identically — the null distribution the result is judged against.</p>
      <div class="card" id="calBody"><div id="histChart"></div>
        <p class="subnote" style="margin-top:12px">
        <b>Discriminating power: ${esc(job.calibration.discriminating_power)}.</b>
        That is how far the known binder sits above random members, and it is
        the honest ceiling on how much any percentile below should be
        trusted.</p></div>`;
    histogram($('histChart'), job.calibration, null);
  }
}

function addRunToRail(job) {
  const list = $('runList');
  if (list.dataset.empty !== 'no') { list.innerHTML = ''; list.dataset.empty = 'no'; }
  if ($(`run-${job.id}`)) return;
  const d = document.createElement('div');
  d.id = `run-${job.id}`;
  d.style.cssText = 'padding:4px 0;font-size:12.5px';
  d.innerHTML = `<a href="#run/${job.id}">${esc(job.params.target)}</a>`;
  list.prepend(d);
}

function renderResults(job) {
  const v = $('view');
  const r = job.result, cal = job.calibration, ref = job.reference;
  const best = r.candidates[0];
  S.selected = 1;

  v.innerHTML = `
    <div class="crumb"><a href="#">All targets</a> <span>›</span>
      <a href="#target/${esc(job.params.target)}">${esc(job.params.target)}</a>
      <span>›</span><span>result</span></div>
    <h1>Design result</h1>
    <p class="secnote">${r.n_unique.toLocaleString()} unique sequences scored in
      ${r.elapsed_s}s · ${esc(r.plm_status)} · stopped because
      ${esc(r.stop_reason)}</p>

    <div class="grid g4">
      ${tile('Recognition', best.recognition.toFixed(1),
        'out of 100 on the interface-complementarity components')}
      ${tile('Percentile vs random', `${best.recognition_percentile.toFixed(0)}th`,
        `beats ${best.recognition_percentile.toFixed(0)}% of ${cal.n_decoys}
         random library members (z = ${best.recognition_z >= 0 ? '+' : ''}${best.recognition_z})`)}
      ${tile('Vs the known binder',
        best.percent_of_reference != null ? `${best.percent_of_reference.toFixed(0)}%` : 'n/a',
        best.percent_of_reference != null && best.percent_of_reference > 100
          // Above 100% is a statement about the scoring function, not the
          // molecule. A real antibody binds this target; a sequence that
          // out-scores it has saturated the criteria being measured, which
          // means the remaining difference is in things not measured here.
          ? `higher than ${ref ? esc(ref.name) : 'the reference'} on the
             computable criteria. Read that as the scoring function running
             out of discriminating power, not as a better binder — the
             reference is known to work and this is untested.`
          : `of ${ref ? esc(ref.name) : 'the reference'}'s recognition score — on
             computable criteria only, not measured binding`)}
      ${best.groups.selectivity == null
      // Never render "not measured" as a number. With no paralog panel the
      // component is constant and correctly excluded from ranking, but
      // showing 0 would read as terrible selectivity rather than as an
      // absent measurement, which are opposite claims.
      ? selectivityAbsentTile(r, best)
      : selectivityTile(r, best)}
    </div>
    ${nothingQualifiedBanner(r)}

    <section>
      <h2>The designed protein</h2>
      <p class="secnote">The scaffold is a real fold; only the
        ${r.scaffold.variable_1based.length} designed positions vary. Click a
        row in the table below to load a different candidate.</p>
      <div class="grid g2">
        <div class="card">
          <h3>Candidate on the ${esc(r.scaffold.pdb_template || 'scaffold')} fold</h3>
          <div class="viewer sm" id="v_design"></div>
          <div class="ctrls">
            <button id="b_dsurface" aria-pressed="false">Surface</button>
            <button id="b_dspin" aria-pressed="false">Spin</button>
          </div>
          <div class="vlegend">
            <span><i class="sw" style="background:#9aa0a6"></i>fixed framework</span>
            <span><i class="sw" style="background:var(--critical)"></i>acidic</span>
            <span><i class="sw" style="background:var(--s4)"></i>aromatic</span>
            <span><i class="sw" style="background:var(--s1)"></i>basic</span>
            <span><i class="sw" style="background:var(--s3)"></i>other designed</span>
          </div>
        </div>
        <div class="card">
          <h3>Score breakdown</h3>
          <div id="brkHead" style="margin-bottom:14px"></div>
          <div class="bars" id="brk"></div>
          <ul class="flags" id="brkFlags"></ul>
        </div>
      </div>
    </section>

    <section>
      <h2>The search</h2>
      <div class="grid g2">
        <div class="card"><h3>Score over generations</h3><div id="convChart"></div>
          ${trajectoryNote(r)}</div>
        <div class="card"><h3>Population diversity</h3>
          <p class="subnote">Mean pairwise difference across the designed
          positions. When this collapses the search has converged.</p>
          <div id="divChart"></div></div>
      </div>
    </section>

    <section>
      <h2>Is this better than chance?</h2>
      <p class="secnote">${cal.n_decoys} random draws from the same library.
        The known binder is marked for scale — if it did not stand clear of
        the decoys, the scoring function would be worthless and so would every
        other number here.</p>
      <div class="card"><div id="histChart2"></div>
        <p class="subnote" style="margin-top:12px">
        <b>Discriminating power: ${esc(cal.discriminating_power)}.</b></p></div>
    </section>

    <section>
      <h2>Candidates</h2>
      <div class="card" style="overflow-x:auto"><table id="tbl"></table></div>
    </section>

    <section>
      <h2>Method and limits</h2>
      <div class="card">
        <dl>
          <dt>Scaffold</dt><dd>${esc(r.scaffold.name)} — ${esc(r.scaffold.description)}</dd>
          <dt>Search space</dt><dd>${r.scaffold.search_space.toExponential(1)} possible
            sequences; ${r.n_unique.toLocaleString()} explored
            (${((r.n_unique / r.scaffold.search_space) * 100).toExponential(1)}%)</dd>
          <dt>Language model</dt><dd>${esc(r.plm_status)}</dd>
          <dt>Reference</dt><dd>${ref ? esc(ref.description) : 'unavailable'}</dd>
        </dl>
        <details><summary>Scaffold provenance</summary>
          <div class="body">${esc(r.scaffold.provenance)}</div></details>
        ${r.off_targets.length ? `<details><summary>Off-target panel
          (${r.off_targets.length} proteins)</summary><div class="body">
          ${r.off_targets.map((o) => `<p style="margin:0 0 8px"><b>${esc(o.gene)}</b>
            (${esc(o.uniprot)})${o.structural ? '' : ' — not a structural relative'}:
            ${esc(o.concern)}</p>`).join('')}</div></details>`
      : `<details><summary>Why selectivity was not assessed</summary>
          <div class="body">No curated panel of related proteins exists for
          this target in the catalogue, so cross-reactivity could not be
          scored. The component reports a neutral placeholder rather than a
          measurement, and it is not evidence the candidate is selective.
          </div></details>`}
        <details><summary>What this result does not mean</summary>
          <div class="body">The candidate scores well on criteria a computer
          can evaluate. That makes it a reasonable thing to synthesise first,
          and nothing more. It has not been shown to bind. Establishing that
          needs expression, a binding assay, then primary cells, then animal
          work, then trials. Most candidates that look this good still fail.
          </div></details>
      </div>
    </section>`;

  const traj = r.trajectory;
  lineChart($('convChart'), [
    { name: 'best', color: '--s1', pts: traj.map((g) => [g.generation, g.best]) },
    { name: 'mean', color: '--s2', pts: traj.map((g) => [g.generation, g.mean]) },
  ], { xlab: 'generation', dec: 0 });
  lineChart($('divChart'), [
    { name: 'diversity', color: '--s3', pts: traj.map((g) => [g.generation, g.diversity]) },
  ], { xlab: 'generation', dec: 2, ymin: 0 });
  histogram($('histChart2'), cal, best.recognition);

  renderTable();
  selectCandidate(1);
  showDesigned('v_design', r.scaffold, best);
  toggle('b_dsurface', () => styleDesigned('v_design', S.job.result.candidates[S.selected - 1]));
  toggle('b_dspin', () => {
    const e = S.viewers['v_design'];
    if (e) e.v.spin($('b_dspin').getAttribute('aria-pressed') === 'true' ? 'y' : false);
  });
}

function tile(lab, val, foot) {
  // Non-numeric values ("not assessed") get smaller, muted type so they do
  // not read as a measurement sitting beside the real ones.
  const bare = /^[\d.+-]/.test(String(val));
  const style = bare ? '' :
    ' style="font-size:17px;color:var(--ink-3);font-weight:600"';
  return `<div class="tile"><div class="lab">${esc(lab)}</div>
    <div class="val"${style}>${esc(val)}</div>
    <div class="foot">${foot}</div></div>`;
}

// The one case a reader cannot detect from the table. When no candidate in
// the whole search complements the target more than its nearest relative,
// ranking has nothing to prefer and falls back to the plain composite — so
// the top row looks exactly as it always did, and is the least bad of a bad
// set. Saying nothing here would make an empty result indistinguishable from
// a good one.
// "Best" means the best RANKED candidate, and ranking puts anything
// counter-selective below everything else. So this line can step down at the
// generation where the search first finds a clean candidate — a 72.6
// displacing an 82.6 on TNF-alpha. Unexplained, that looks like the search
// going backwards; it is the opposite.
function trajectoryNote(r) {
  const t = r.trajectory || [];
  const switched = t.some((g, i) => i > 0
    && t[i - 1].best_disqualified && !g.best_disqualified);
  if (!switched) return '';
  const at = t.findIndex((g, i) => i > 0
    && t[i - 1].best_disqualified && !g.best_disqualified);
  const drop = t[at - 1].best - t[at].best;
  return `<p class="subnote" style="margin-top:10px">
    At generation ${t[at].generation} the best-ranked candidate changed kind.
    Up to there it was counter-selective — it complemented a relative more
    than ${esc(r.gene || 'the target')} — while scoring well on everything
    else. ${drop > 0.05
      ? `There the search first found one that was not, and the line steps
         <b>down</b> by ${drop.toFixed(1)} points: the clean candidate scores
         lower on the other axes, and is still the better molecule. The step
         is the run starting to work, not regressing.`
      : `There the search first found one that was not, at much the same
         composite — so the line barely moves even though what it is tracking
         changed.`}</p>`;
}

function nothingQualifiedBanner(r) {
  if (r.any_qualified !== false) return '';
  return `
    <div class="card" style="border-left:3px solid var(--serious);margin-top:14px">
      <h3>No candidate in this search was selective for ${esc(r.gene || 'the target')}</h3>
      <p style="font-size:13.5px;color:var(--ink-2);margin:6px 0 0">
        Every sequence the loop produced complements the nearest relative in
        the panel at least as well as it complements the target. The table
        below is therefore ordered by composite score alone, as it would have
        been with no selectivity measurement at all, and its top row is the
        least bad of a set in which nothing passed. Read it as a result about
        this scaffold against this epitope — thirteen positions on a
        58-residue helix bundle may simply not be enough surface to tell
        these two proteins apart — and not as a shortlist.
      </p>
    </div>`;
}

function renderTable() {
  const r = S.job.result;
  const t = $('tbl');
  t.innerHTML = `<thead><tr>
    <th>#</th><th>designed positions</th><th>gen</th><th>origin</th>
    <th>composite</th><th>recognition</th><th>pctile</th><th>z</th>
    <th>vs binder</th></tr></thead><tbody>` +
    r.candidates.map((c) => `<tr data-r="${c.rank}" class="${
      c.rank === S.selected ? 'sel' : ''}"${(c.disqualifying || []).length
        ? ` title="${esc((c.disqualifying || []).map((d) => d.why).join('; '))}"` : ''}>
      <td>${c.rank}</td><td class="mono">${esc(c.designed)}</td>
      <td>${c.generation}</td><td style="text-align:left">${esc(c.origin)}</td>
      <td>${c.composite}${(c.disqualifying || []).length
        ? ' <span title="ranked below every candidate that is not disqualified"'
          + ' style="color:var(--serious);font-weight:600">!</span>' : ''}</td>
      <td>${c.recognition}</td>
      <td>${c.recognition_percentile}</td>
      <td>${c.recognition_z >= 0 ? '+' : ''}${c.recognition_z}</td>
      <td>${c.percent_of_reference != null ? c.percent_of_reference + '%' : '—'}</td>
    </tr>`).join('') + '</tbody>';
  t.querySelectorAll('tbody tr').forEach((tr) => {
    tr.onclick = () => selectCandidate(+tr.dataset.r);
  });
}

const GCOL = {
  recognition: '--s1', selectivity: '--s2', developability: '--s3',
  plausibility: '--s4', other: '--ink-3',
};

function selectCandidate(rank) {
  S.selected = rank;
  const r = S.job.result;
  const c = r.candidates[rank - 1];
  const varPos = r.scaffold.variable_1based;

  $('brkHead').innerHTML = `
    <div style="font-size:12.5px;color:var(--ink-2)">Candidate #${c.rank},
      generation ${c.generation} via ${esc(c.origin)}. Composite
      <b style="color:var(--ink-1)">${c.composite}</b>/100 over ${c.n_ranked}
      ranking components.${(c.disqualifying || []).length
        ? ` <b style="color:var(--serious)">Disqualified:</b> ${
            esc(c.disqualifying.map((d) => d.why).join('; '))}. The composite
          is still the honest average of everything measured — this candidate
          genuinely does score that well on the other axes. It is ranked below
          every candidate that is not disqualified, because doing the opposite
          of the job is not a quantity that trades against solubility.` : ''}</div>
    <div class="seq" style="margin-top:9px">${[...c.annotated].map((ch) =>
      /[A-Z]/.test(ch) ? `<span class="v">${ch}</span>` : ch).join('')}</div>
    <div style="font-size:11.5px;color:var(--ink-3);margin-top:5px">
      Highlighted = the ${varPos.length} designed positions
      (${varPos.join(', ')}). Lowercase = fixed framework.</div>`;

  const order = { recognition: 0, selectivity: 1, plausibility: 2, developability: 3, other: 4 };
  const sorter = (a, b) => (order[a.group] - order[b.group]) || (b.weight - a.weight);
  const ranked = c.components.filter((k) => k.ranked !== false).sort(sorter);
  const unranked = c.components.filter((k) => k.ranked === false).sort(sorter);
  const row = (k) => `<div class="bar${k.ranked === false ? ' unranked' : ''}"
      data-detail="${encodeURIComponent(k.detail)}"
      data-caveat="${encodeURIComponent(k.caveat || '')}"
      data-note="${encodeURIComponent(k.note || '')}"
      data-nm="${esc(k.name)}" data-g="${esc(k.group)}" data-w="${k.weight}">
      <div class="nm">${esc(k.name.replace(/_/g, ' '))}</div>
      <div class="track"><div class="fill" style="width:${k.normalised}%;
        background:${k.ranked === false ? cssv('--ink-3') : cssv(GCOL[k.group] || '--ink-3')}"></div></div>
      <div class="pc">${k.normalised.toFixed(0)}%</div></div>`;

  $('brk').innerHTML = ranked.map(row).join('') +
    (unranked.length ? `<div class="subhead">Measured, not ranked</div>
      <p class="subnote">Constant across every candidate on this scaffold, so
      they cannot separate designs. Still real checks on the molecule; they
      would rank again on a scaffold where they vary.</p>`
      + unranked.map(row).join('') : '');

  document.querySelectorAll('#brk .bar').forEach((b) => {
    b.addEventListener('mousemove', (ev) => showTip(
      `<b>${esc(b.dataset.nm.replace(/_/g, ' '))}</b> · ${esc(b.dataset.g)} ·
       weight ${esc(b.dataset.w)}
       <div style="margin-top:6px">${esc(decodeURIComponent(b.dataset.detail))}</div>` +
      (b.dataset.note ? `<div style="margin-top:6px;color:var(--ink-3)">
        <b>Not ranked:</b> ${esc(decodeURIComponent(b.dataset.note))}</div>` : '') +
      (b.dataset.caveat ? `<div style="margin-top:6px;color:var(--ink-3)">
        <b>Cannot tell you:</b> ${esc(decodeURIComponent(b.dataset.caveat))}</div>` : ''),
      ev.clientX, ev.clientY));
    b.addEventListener('mouseleave', hideTip);
  });

  $('brkFlags').innerHTML = (c.flags || []).map((f) => `<li>${esc(f)}</li>`).join('');
  renderTable();
  if (S.viewers['v_design']) styleDesigned('v_design', c);
}


/* ---------- build your own target ---------- */
/*
 * The catalogue is 74 proteins somebody curated. This is the same pipeline
 * pointed at anything in UniProt: resolve the name, find a structure that
 * holds the protein with something bound to it, read the binding site off
 * that structure, derive the cross-reactivity panel, measure what
 * selectivity is reachable, and hand back a target the design loop runs
 * against exactly as it runs against a curated one.
 *
 * The one thing it cannot do is invent a binding site. This project reads
 * epitopes off structures of a complex. A protein solved only on its own, or
 * not solved at all, yields no target -- and most human proteins are in that
 * position. The page says so in those words, because "no results" would read
 * as a failure of the software or of the protein rather than as the state of
 * structural biology.
 */

const BUILD = { query: '', proteins: [], picked: null, structures: [],
                busy: false, error: '', job: null, made: [] };

async function renderBuild() {
  const v = $('view');
  v.innerHTML = `
    <div class="crumb"><a href="#catalog">All targets</a> &rsaquo; Build your own</div>
    <h1>Design against any protein</h1>
    <p class="secnote">Name a protein and this runs the whole pipeline on it:
      the binding site read off an experimental structure, a cross-reactivity
      panel derived from its domain family, and a design cycle scored the
      same way every catalogue target is scored.</p>

    <div class="card" style="border-left:3px solid var(--serious)">
      <h3>What this can and cannot do</h3>
      <p style="font-size:13.5px;color:var(--ink-2);margin:6px 0 0">
        It needs an <b>experimental structure of your protein with something
        bound to it</b>. The binding site is measured from that complex, never
        predicted — which is why the numbers here mean anything at all, and
        why a protein with no such structure cannot be used. Most human
        proteins have none. That is a fact about what has been solved, not
        about whether the protein matters.</p>
      <p style="font-size:13.5px;color:var(--ink-2);margin:8px 0 0">
        A target you build is <b>not reviewed</b>. Catalogue entries were
        checked by a person and carry a written rationale and a clinical
        caveat. These carry neither. The automatic checks are the same ones
        the catalogue build runs, and anything they flag is shown.</p>
    </div>

    <section>
      <h2>1 &middot; Which protein?</h2>
      <p class="secnote">A gene symbol (<span class="mono">TNF</span>,
        <span class="mono">KRAS</span>, <span class="mono">SOD1</span>) or a
        UniProt accession (<span class="mono">P01375</span>).</p>
      <div class="card">
        <div class="ctrls" style="gap:8px">
          <input type="search" id="bq" placeholder="gene symbol or accession"
                 value="${esc(BUILD.query)}" style="flex:1;min-width:220px"
                 aria-label="Protein to search for">
          <button id="bgo">Search</button>
        </div>
        <div id="bprot" style="margin-top:12px"></div>
      </div>
    </section>

    <section id="bstructsec"></section>
    <section id="bjobsec"></section>
    <section id="bmadesec"></section>`;

  $('bgo').onclick = buildSearch;
  $('bq').addEventListener('keydown', (e) => {
    if (e.key === 'Enter') buildSearch();
  });
  renderBuildProteins();
  renderBuildStructures();
  renderBuildJob();
  loadMade();
}

async function loadMade() {
  try {
    const r = await api('/api/custom/list');
    BUILD.made = r.targets || [];
  } catch (e) { BUILD.made = []; }
  const host = $('bmadesec');
  if (!host) return;
  if (!BUILD.made.length) { host.innerHTML = ''; return; }
  host.innerHTML = `
    <h2>Targets you have built</h2>
    <p class="secnote">These live as long as the server does. They are not
      written into the catalogue, which is a reviewed file.</p>
    <div class="grid g3">${BUILD.made.map((t) => `
      <a class="card tcard" href="#target/${esc(t.key)}">
        <div class="row"><b class="mono">${esc(t.gene)}</b>
          <span class="pill">${esc(t.pdb)}</span></div>
        <div style="font-size:12.5px;color:var(--ink-2);margin-top:4px">
          ${esc(t.protein_name)}</div>
        <div class="mono" style="font-size:11.5px;color:var(--ink-3);margin-top:6px">
          ${esc(t.epitope_sequence)}</div>
        <div style="font-size:11.5px;color:var(--ink-3);margin-top:6px">
          ${t.n_epitope_residues} contact residues</div>
      </a>`).join('')}</div>`;
}

async function buildSearch() {
  const q = ($('bq').value || '').trim();
  if (!q) return;
  BUILD.query = q; BUILD.busy = true; BUILD.error = '';
  BUILD.proteins = []; BUILD.picked = null; BUILD.structures = [];
  renderBuildProteins(); renderBuildStructures();
  try {
    const r = await api(`/api/custom/search?q=${encodeURIComponent(q)}`);
    BUILD.proteins = r.proteins || [];
  } catch (e) {
    BUILD.error = e.message;
  }
  BUILD.busy = false;
  renderBuildProteins();
}

function renderBuildProteins() {
  const host = $('bprot');
  if (!host) return;
  if (BUILD.busy) {
    host.innerHTML = `<p class="secnote"><span class="spin"></span>Asking
      UniProt…</p>`;
    return;
  }
  if (BUILD.error) {
    host.innerHTML = `<p class="secnote" style="color:var(--serious)">
      ${esc(BUILD.error)}</p>`;
    return;
  }
  if (!BUILD.query) { host.innerHTML = ''; return; }
  if (!BUILD.proteins.length) {
    host.innerHTML = `<p class="secnote">No reviewed human protein has the
      gene symbol &ldquo;${esc(BUILD.query)}&rdquo;. Reviewed means a curated
      Swiss-Prot entry; try the accession, or another name for the gene.</p>`;
    return;
  }
  // Every match, not the first. The first is not reliably the right one:
  // the symbol SELP returns selenoprotein P ahead of P-selectin.
  host.innerHTML = `
    ${BUILD.proteins.length > 1 ? `<p class="secnote">${BUILD.proteins.length}
      proteins answer to that name. Pick the one you mean.</p>` : ''}
    <table><thead><tr><th>Accession</th><th>Protein</th><th>Gene</th>
      <th>Length</th><th></th></tr></thead><tbody>
    ${BUILD.proteins.map((p) => `<tr>
      <td class="mono">${esc(p.accession)}</td>
      <td style="text-align:left">${esc(p.name)}</td>
      <td class="mono">${esc((p.genes || []).join(', '))}</td>
      <td>${p.length || '—'} aa</td>
      <td><button class="ghost" data-acc="${esc(p.accession)}">Structures</button></td>
    </tr>`).join('')}</tbody></table>`;
  host.querySelectorAll('button[data-acc]').forEach((b) => {
    b.onclick = () => pickProtein(b.dataset.acc);
  });
}

async function pickProtein(accession) {
  BUILD.picked = BUILD.proteins.find((p) => p.accession === accession) || null;
  BUILD.structures = []; BUILD.busy = true;
  renderBuildStructures();
  try {
    const r = await api(`/api/custom/structures?accession=${encodeURIComponent(accession)}`);
    BUILD.structures = r.structures || [];
  } catch (e) {
    BUILD.error = e.message;
  }
  BUILD.busy = false;
  renderBuildStructures();
}

function renderBuildStructures() {
  const host = $('bstructsec');
  if (!host) return;
  if (!BUILD.picked) { host.innerHTML = ''; return; }
  if (BUILD.busy) {
    host.innerHTML = `<h2>2 &middot; Which structure?</h2>
      <p class="secnote"><span class="spin"></span>Asking the Protein Data
      Bank which entries hold ${esc(BUILD.picked.accession)} with a partner
      bound…</p>`;
    return;
  }
  if (!BUILD.structures.length) {
    // The common case, and the one most worth explaining properly.
    host.innerHTML = `<h2>2 &middot; Which structure?</h2>
      <div class="card" style="border-left:3px solid var(--serious)">
        <h3>No usable structure for ${esc(BUILD.picked.name)}</h3>
        <p style="font-size:13.5px;color:var(--ink-2);margin:6px 0 0">
          The Protein Data Bank has no experimental entry holding this protein
          with another molecule bound to it. There may well be structures of
          it <em>alone</em> — but a binding site is read off a complex, and
          this pipeline will not invent one.</p>
        <p style="font-size:13.5px;color:var(--ink-2);margin:8px 0 0">
          This is the normal situation for most of the proteome, and it is the
          honest limit of the method rather than a verdict on the target. A
          protein in this position needs a structure solved, or a predicted
          complex that someone is prepared to stand behind.</p>
      </div>`;
    return;
  }
  host.innerHTML = `
    <h2>2 &middot; Which structure?</h2>
    <p class="secnote">Entries that hold
      <span class="mono">${esc(BUILD.picked.accession)}</span> with something
      bound to it. Antibody complexes come first: they are the cleanest
      two-body interfaces, and the bound molecule is usually a drug or a
      drug candidate.</p>
    <div class="card">
      <table><thead><tr><th>PDB</th><th>Resolution</th><th>Bound</th>
        <th style="text-align:left">Title</th><th></th></tr></thead><tbody>
      ${BUILD.structures.map((st) => `<tr>
        <td class="mono">${esc(st.pdb_id)}</td>
        <td>${st.resolution ? st.resolution.toFixed(2) + ' Å' : '—'}</td>
        <td>${st.has_antibody
          ? '<span class="pill">antibody</span>'
          : '<span class="pill">partner</span>'}${
          st.large_assembly ? ' <span class="pill">large</span>' : ''}</td>
        <td style="text-align:left;color:var(--ink-3)">${esc(st.title)}
          <div style="font-size:11.5px;margin-top:3px">target chain
            ${esc((st.target_chains || []).join(','))}, ${st.target_length} aa
            &middot; bound ${esc((st.partners || []).map((o) =>
              `${o.description} (${o.length} aa)`).join('; '))}</div></td>
        <td><button data-pdb="${esc(st.pdb_id)}">Build</button></td>
      </tr>`).join('')}</tbody></table>
      <p class="subnote" style="margin-top:10px">Resolution matters: below
        about 3 Å the contacts are well determined, above it the binding site
        is approximate and every number downstream inherits that.</p>
    </div>`;
  host.querySelectorAll('button[data-pdb]').forEach((b) => {
    b.onclick = () => startBuild(b.dataset.pdb);
  });
}

async function startBuild(pdb) {
  if (!BUILD.picked) return;
  try {
    const r = await api('/api/custom/prepare', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ accession: BUILD.picked.accession, pdb }),
    });
    BUILD.job = { id: r.job, state: 'queued', stage: 'queued', progress: 0 };
    renderBuildJob();
    pollBuild(r.job);
  } catch (e) {
    BUILD.job = { state: 'failed', error: e.message };
    renderBuildJob();
  }
}

async function pollBuild(id) {
  for (;;) {
    let j;
    try { j = await api(`/api/job/${id}`); } catch (e) { return; }
    BUILD.job = j;
    renderBuildJob();
    if (j.state === 'done' || j.state === 'failed') {
      if (j.state === 'done') loadMade();
      return;
    }
    await new Promise((r) => setTimeout(r, 700));
  }
}

function renderBuildJob() {
  const host = $('bjobsec');
  if (!host) return;
  const j = BUILD.job;
  if (!j) { host.innerHTML = ''; return; }
  if (j.state === 'failed') {
    // The pipeline's own wording, not a rephrasing of it: these messages
    // were written to tell someone what to try next.
    host.innerHTML = `
      <h2>3 &middot; Building</h2>
      <div class="card" style="border-left:3px solid var(--serious)">
        <h3>This structure did not yield a binding site</h3>
        <p style="font-size:13.5px;color:var(--ink-2);margin:6px 0 0">
          ${esc(j.error || 'unknown error')}</p>
        <p class="subnote" style="margin-top:10px">Nothing is wrong with your
          protein. Try another entry from the list — a different complex of
          the same protein often works.</p>
      </div>`;
    return;
  }
  if (j.state === 'done' && j.result) {
    const r = j.result;
    host.innerHTML = `
      <h2>3 &middot; Built</h2>
      <div class="card">
        <h3>${esc(r.record.gene)} &middot; ${esc(r.record.pdb_complex)}</h3>
        <p style="font-size:13.5px;color:var(--ink-2);margin:6px 0 0">
          ${r.record.n_epitope_residues} contact residues read off chain
          ${esc(r.record.target_chain)}, bound by
          ${esc((r.record.binder_chains || []).join(','))}.</p>
        <div class="mono" style="margin-top:8px;font-size:12.5px">
          ${esc(r.record.epitope_sequence)}</div>
        ${(r.notes || []).length ? `<ul class="subnote" style="margin-top:10px">
          ${r.notes.map((n) => `<li>${esc(n)}</li>`).join('')}</ul>` : ''}
        <div class="ctrls" style="margin-top:12px">
          <button class="primary" id="bopen">Open this target</button>
        </div>
      </div>`;
    const go = $('bopen');
    if (go) go.onclick = () => { location.hash = `#target/${r.key}`; };
    return;
  }
  host.innerHTML = `
    <h2>3 &middot; Building</h2>
    <div class="card">
      <div class="prog"><div style="width:${Math.max(2, (j.progress || 0) * 100)}%"
        class="progbar"></div></div>
      <div class="stage"><span class="spin"></span>${esc(j.stage || 'working')}</div>
      <p class="subnote" style="margin-top:10px">Downloading the structure,
        reading the interface, then asking UniProt which human proteins share
        a domain family with this one and aligning each at the binding site.
        The last step is what makes a selectivity number possible, and it is
        the slow one.</p>
    </div>`;
}

/* ================= specimen profile =================
 *
 * Target selection from an antigen panel. The whole reason this view is
 * separate from a design run is that it answers a different question — which
 * protein should a lab aim at for THIS specimen — and the honest answer leans
 * on what is known about healthy cells, which no design score can see.
 *
 * Two rules the rendering has to respect, because breaking either of them is
 * how a careful tool becomes a dangerous one:
 *
 *   1. "Unknown" is shown as loudly as "concern". A lineage nobody has data
 *      for is not a clean result, and if it renders faintly it reads as one.
 *   2. The evidence kind is on every single row, in words. "Absent by
 *      cytometry" and "no transcript detected" are not the same claim and
 *      must never look alike.
 */

const SP = {
  setup: null, area: 'Blood cancer', q: '',
  antigens: {},            // GENE -> positive | dim | negative | not_tested
  spared: new Set(['hsc', 'progenitor', 'erythroid', 'platelet']),
  variants: [],
  result: null, busy: false, error: '',
  open: {},                // target key -> expanded
};

const LEVEL_WORD = {
  high: 'carried, high', moderate: 'carried, moderate', low: 'carried, low',
  contested: 'sources disagree', absent: 'not carried', unknown: 'unknown',
};
const LEVEL_COLOR = {
  high: 'critical', moderate: 'serious', low: 'warning',
  contested: 'warning', absent: 'good', unknown: 'ink-3',
};
const STATE_SHORT = { positive: '+', dim: 'dim', negative: '−',
                      not_tested: '?' };

async function renderSpecimen() {
  const v = $('view');
  if (!SP.setup) {
    v.innerHTML = `<p class="secnote"><span class="spin"></span>Loading the
      antigen list…</p>`;
    try {
      SP.setup = await api('/api/specimen/setup');
    } catch (e) {
      v.innerHTML = `<p class="secnote">Could not load: ${esc(e.message)}</p>`;
      return;
    }
  }
  const s = SP.setup;
  const areas = [...new Set(s.antigens.map((a) => a.area))].sort();

  v.innerHTML = `
    <div class="crumb"><a href="#catalog">All targets</a> &rsaquo; Specimen profile</div>
    <h1>Which target, for this specimen</h1>
    <p class="secnote">Enter what the flow cytometry says the blasts carry and
      which healthy populations have to survive. This ranks the catalogue's
      targets against that profile and says why each one placed where it did.</p>

    <div class="card" style="border-left:3px solid var(--serious)">
      <h3>What this is, and what it is not</h3>
      <p style="font-size:13.5px;color:var(--ink-2);margin:6px 0 0">
        This is <b>target selection</b> &mdash; which protein a lab should
        build a binder against first. It is not treatment selection, it does
        not predict response, and it does not design a combination. Nothing
        here is a clinical recommendation.</p>
      <p style="font-size:13.5px;color:var(--ink-2);margin:8px 0 0">
        There are <b>no patient identifiers</b> in this form and nothing is
        stored. It reasons about a cell population.</p>
      <p style="font-size:13.5px;color:var(--ink-2);margin:8px 0 0">
        Targets rank well here for being <b>well characterised</b>, not for
        looking clean. A lineage nobody has examined counts against a target,
        because the alternative ranking &mdash; fewest concerns first &mdash;
        puts the least studied protein at the top of every list.</p>
    </div>

    <section>
      <h2>1 &middot; Antigen panel</h2>
      <p class="secnote">What the blasts carry. <b>Not tested</b> is not the
        same as <b>negative</b>, and the two are kept apart throughout.</p>
      <div class="card">
        <div class="ctrls" style="gap:8px;margin-top:0">
          <select id="sparea" aria-label="Disease area">
            ${areas.map((a) => `<option ${a === SP.area ? 'selected' : ''}
              >${esc(a)}</option>`).join('')}
          </select>
          <input type="search" id="spq" placeholder="filter antigens"
                 value="${esc(SP.q)}" style="flex:1;min-width:160px"
                 aria-label="Filter antigens">
          <span class="subnote" id="spcount"></span>
        </div>
        <div id="splist" style="margin-top:12px"></div>
      </div>
    </section>

    <section>
      <h2>2 &middot; What has to survive</h2>
      <p class="secnote">Each population here is one a binder could destroy
        along with the blasts. Tick the ones this specimen cannot afford to
        lose.</p>
      <div class="card">
        <div class="ctrls" style="margin-top:0">
          ${s.presets.map((p) => `<button class="ghost" data-preset="${esc(p.key)}"
            >${esc(p.label)}</button>`).join('')}
          <button class="ghost" data-preset="__none">Clear</button>
        </div>
        <div id="splin" style="margin-top:12px"></div>
      </div>
    </section>

    <section>
      <h2>3 &middot; Variants <span class="subnote">optional</span></h2>
      <p class="secnote">A coding change called in a target gene, numbered
        against the reference sequence. If it lands on a residue inside the
        binding site, a binder designed against the reference protein may not
        fit this specimen's copy &mdash; the mechanism behind epitope-loss
        relapse after CD19-directed therapy.</p>
      <div class="card">
        <div class="ctrls" style="margin-top:0">
          <input id="spvg" placeholder="gene" style="width:110px"
                 aria-label="Variant gene">
          <input id="spvf" placeholder="ref aa" style="width:78px" maxlength="1"
                 aria-label="Reference amino acid">
          <input id="spvp" placeholder="position" style="width:96px"
                 inputmode="numeric" aria-label="Position">
          <input id="spvt" placeholder="new aa" style="width:78px" maxlength="1"
                 aria-label="Variant amino acid">
          <button id="spvadd">Add variant</button>
        </div>
        <div id="spvlist" style="margin-top:10px"></div>
      </div>
    </section>

    <div class="ctrls" style="margin:18px 0">
      <button class="primary" id="sprank" ${SP.busy ? 'disabled' : ''}>
        ${SP.busy ? 'Ranking…' : 'Rank targets'}</button>
      <span class="subnote">${esc(s.expression.n_genes)} genes of transcript
        data, ${esc(s.expression.source)}, ${esc(s.expression.version)},
        fetched ${esc(s.expression.fetched)}.</span>
    </div>
    <div id="sperr"></div>
    <section id="spres"></section>`;

  $('sparea').onchange = (e) => { SP.area = e.target.value; spRenderList(); };
  $('spq').addEventListener('input', (e) => {
    SP.q = e.target.value; spRenderList();
  });
  v.querySelectorAll('[data-preset]').forEach((b) => {
    b.onclick = () => spPreset(b.dataset.preset);
  });
  $('spvadd').onclick = spAddVariant;
  $('sprank').onclick = spRank;

  spRenderList();
  spRenderLineages();
  spRenderVariants();
  spRenderResult();
}

function spRenderList() {
  const host = $('splist');
  if (!host) return;
  const q = SP.q.toLowerCase().trim();
  const rows = SP.setup.antigens.filter((a) => {
    if (q) {
      const hay = [a.gene, a.protein_name, a.disease, ...(a.aliases || [])]
        .join(' ').toLowerCase();
      return hay.includes(q);
    }
    return a.area === SP.area;
  });
  const n = Object.values(SP.antigens)
    .filter((x) => x && x !== 'not_tested').length;
  $('spcount').textContent = `${n} antigen${n === 1 ? '' : 's'} entered`;

  if (!rows.length) {
    host.innerHTML = `<p class="subnote">Nothing matches.</p>`;
    return;
  }
  host.innerHTML = rows.map((a) => {
    const cur = SP.antigens[a.gene.toUpperCase()] || '';
    const alias = (a.aliases || []).filter((x) => /^CD\d/i.test(x));
    return `
    <div style="display:flex;align-items:center;gap:10px;padding:6px 0;
                border-bottom:1px solid var(--border)">
      <div style="flex:1;min-width:0">
        <span class="mono" style="font-size:13.5px">${esc(a.gene)}</span>
        ${alias.length ? `<span class="subnote"> ${esc(alias.join(' / '))}</span>`
                       : ''}
        <div class="subnote" style="white-space:nowrap;overflow:hidden;
             text-overflow:ellipsis">${esc(a.disease)}</div>
      </div>
      ${!a.has_expression ? `<span class="pill warn" title="No human
        transcript data: this target is not a human protein."
        >non-human</span>` : ''}
      <div style="display:flex;gap:3px">
        ${SP.setup.antigen_states.map((st) => `
          <button class="ghost" data-ag="${esc(a.gene)}" data-st="${esc(st)}"
            title="${esc(st.replace('_', ' '))}"
            style="padding:3px 9px;font-size:12px;min-width:34px;
              ${cur === st ? `background:var(--s1);color:#fff;
                border-color:var(--s1)` : ''}"
            >${esc(STATE_SHORT[st] || st)}</button>`).join('')}
      </div>
    </div>`;
  }).join('');

  host.querySelectorAll('[data-ag]').forEach((b) => {
    b.onclick = () => {
      const g = b.dataset.ag.toUpperCase();
      SP.antigens[g] = SP.antigens[g] === b.dataset.st ? '' : b.dataset.st;
      spRenderList();
    };
  });
}

function spPreset(key) {
  if (key === '__none') { SP.spared = new Set(); spRenderLineages(); return; }
  const p = SP.setup.presets.find((x) => x.key === key);
  if (!p) return;
  // Additive: a patient can need their marrow back AND be neutropenic, and
  // making the presets exclusive would quietly drop one of the two.
  p.lineages.forEach((l) => SP.spared.add(l));
  spRenderLineages();
}

function spRenderLineages() {
  const host = $('splin');
  if (!host) return;
  host.innerHTML = SP.setup.lineages.map((l) => `
    <label style="display:flex;gap:10px;padding:7px 0;align-items:flex-start;
                  border-bottom:1px solid var(--border);cursor:pointer">
      <input type="checkbox" data-lin="${esc(l.key)}" style="margin-top:3px"
        ${SP.spared.has(l.key) ? 'checked' : ''}>
      <span style="flex:1">
        <b style="font-size:13.5px">${esc(l.label)}</b>
        ${!l.screened ? `<span class="pill warn" style="margin-left:6px"
          >no screen data</span>` : ''}
        <div class="subnote">${esc(l.consequence)}</div>
      </span>
    </label>`).join('')
    + `<p class="subnote" style="margin-top:10px"><b>no screen data</b> means
       the transcript atlas has nothing for that population at all &mdash; its
       bone marrow single-cell set contains no stem or progenitor cells. For
       those rows only cytometry recorded in this tool can say anything, and
       where there is none the answer is <i>unknown</i>, never <i>absent</i>.</p>`;

  host.querySelectorAll('[data-lin]').forEach((c) => {
    c.onchange = () => {
      if (c.checked) SP.spared.add(c.dataset.lin);
      else SP.spared.delete(c.dataset.lin);
    };
  });
}

function spAddVariant() {
  const gene = ($('spvg').value || '').trim().toUpperCase();
  const pos = parseInt(($('spvp').value || '').trim(), 10);
  if (!gene || !Number.isFinite(pos) || pos < 1) {
    SP.error = 'A variant needs a gene and a positive position.';
    spRenderResult(); return;
  }
  SP.variants.push({
    gene, position: pos,
    from: ($('spvf').value || '').trim().toUpperCase().slice(0, 1),
    to: ($('spvt').value || '').trim().toUpperCase().slice(0, 1),
  });
  ['spvg', 'spvf', 'spvp', 'spvt'].forEach((i) => { $(i).value = ''; });
  SP.error = '';
  spRenderVariants();
}

function spRenderVariants() {
  const host = $('spvlist');
  if (!host) return;
  if (!SP.variants.length) {
    host.innerHTML = `<p class="subnote">None entered.</p>`; return;
  }
  host.innerHTML = SP.variants.map((v, i) => `
    <span class="pill" style="margin:0 6px 6px 0">
      <span class="mono">${esc(v.gene)} p.${esc(v.from)}${esc(v.position)}${esc(v.to)}</span>
      <button class="ghost" data-rmv="${i}" style="padding:0 4px;font-size:14px;
        line-height:1">&times;</button></span>`).join('');
  host.querySelectorAll('[data-rmv]').forEach((b) => {
    b.onclick = () => {
      SP.variants.splice(parseInt(b.dataset.rmv, 10), 1); spRenderVariants();
    };
  });
}

async function spRank() {
  const antigens = {};
  Object.entries(SP.antigens).forEach(([g, s]) => { if (s) antigens[g] = s; });
  SP.busy = true; SP.error = ''; SP.result = null;
  $('sprank').disabled = true; $('sprank').textContent = 'Ranking…';
  try {
    SP.result = await api('/api/specimen/rank', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        antigens, spared: [...SP.spared], variants: SP.variants,
        areas: [SP.area],
      }),
    });
  } catch (e) {
    SP.error = e.message;
  }
  SP.busy = false;
  const b = $('sprank');
  if (b) { b.disabled = false; b.textContent = 'Rank targets'; }
  spRenderResult();
}

function levelPill(c) {
  const word = LEVEL_WORD[c.level] || c.level;
  const tone = LEVEL_COLOR[c.level] || 'ink-3';
  return `<span class="pill" style="border-color:color-mix(in srgb,
    var(--${tone}) 45%,transparent);color:var(--${tone})">${esc(word)}</span>`;
}

/* Severity of what is LOST by blocking a signal, which is a different
 * question from expression and from harm to a population. The two
 * well-tolerated entries matter as much as the boxed ones: a table where
 * everything is dangerous carries no information. */
const SEVERITY_TONE = {
  boxed_warning: 'critical',
  serious: 'critical',
  manageable: 'warning',
  well_tolerated: 'good',
  contested: 'warning',
};

/* Expression and harm are different questions, so they get different pills.
 * CD47 is carried on haematopoietic stem cells and blocking it still does
 * not destroy them; rendering that identically to CD33 on progenitors, which
 * does cause cytopenias, would collapse the one distinction this table was
 * expanded to make. */
const HARM_PILL = {
  observed: ['damage seen in patients', 'critical'],
  not_observed: ['population survives it', 'good'],
};

function lineageRows(calls, kind) {
  if (!calls.length) return '';
  return calls.map((c) => {
    const h = HARM_PILL[c.harm];
    return `
    <div style="padding:7px 0;border-top:1px solid var(--border)">
      <div style="display:flex;gap:8px;align-items:baseline;flex-wrap:wrap">
        <b style="font-size:13px">${esc(c.label)}</b>
        ${levelPill(c)}
        ${h ? `<span class="pill" style="border-color:color-mix(in srgb,
          var(--${h[1]}) 45%,transparent);color:var(--${h[1]})"
          >${esc(h[0])}</span>` : ''}
        <span class="subnote">${esc(c.evidence_label)}</span>
      </div>
      <div class="subnote" style="margin-top:3px">${esc(c.detail || c.note)}</div>
      ${c.source ? `<div class="subnote" style="margin-top:2px;opacity:.8">
        Source: ${esc(c.source)}</div>` : ''}
    </div>`;
  }).join('');
}

function specimenCard(c) {
  const open = !!SP.open[c.key];
  const nCon = c.conflicts.length, nUnk = c.unknowns.length,
        nCle = c.cleared.length;
  /* Tolerate a server that predates these fields rather than throwing on
   * `.length` of undefined -- the dev server is restarted by hand. */
  const tHarm = c.tissue_harms || [], tSpared = c.tissue_spared || [];
  return `
  <div class="card" style="${c.disqualified
      ? 'border-left:3px solid var(--critical)' : ''}">
    <div style="display:flex;gap:10px;align-items:baseline;flex-wrap:wrap">
      <b class="mono" style="font-size:15px">${esc(c.gene)}</b>
      <span class="subnote">${esc(c.protein_name)}</span>
      <span class="pill">${esc(c.antigen_label)}</span>
      <span style="flex:1"></span>
      ${nCle ? `<span class="pill ok">${nCle} cleared by cytometry</span>` : ''}
      ${nUnk ? `<span class="pill warn">${nUnk} unknown</span>` : ''}
      ${nCon ? `<span class="pill" style="border-color:color-mix(in srgb,
        var(--critical) 45%,transparent);color:var(--critical)"
        >${nCon} carried</span>` : ''}
      ${tHarm.length ? `<span class="pill" style="border-color:color-mix(in
        srgb,var(--critical) 45%,transparent);color:var(--critical)"
        >${tHarm.length} organ${tHarm.length > 1 ? 's' : ''} damaged in
        patients</span>` : ''}
      ${(c.physiology && SEVERITY_TONE[c.physiology.severity])
        ? `<span class="pill" style="border-color:color-mix(in srgb,
            var(--${SEVERITY_TONE[c.physiology.severity]}) 45%,transparent);
            color:var(--${SEVERITY_TONE[c.physiology.severity]})"
            >${esc(c.physiology.severity_label)}</span>` : ''}
      ${(c.topology && c.topology !== 'surface' && c.topology !== 'unknown')
        ? `<span class="pill warn">${esc(c.topology === 'shed'
            ? 'also shed as soluble' : c.topology === 'secreted'
            ? 'secreted, not on a cell' : 'intracellular')}</span>` : ''}
    </div>
    <div class="subnote" style="margin-top:4px">${esc(c.disease)}</div>

    ${c.disqualifying ? `<p style="margin:10px 0 0;font-size:13.5px;
      color:var(--ink-1)"><b>Excluded.</b> ${esc(c.disqualifying)}</p>` : ''}

    ${c.inapplicable ? `<p class="subnote" style="margin:10px 0 0">
      ${esc(c.inapplicable.detail)}</p>` : ''}

    ${c.off_axis ? `<p class="subnote" style="margin:10px 0 0;
      border-left:2px solid var(--border);padding-left:9px">
      ${esc(c.off_axis.detail)}</p>` : ''}

    ${c.physiology ? `<div style="margin:12px 0 0;padding:10px 11px;
      background:var(--surface-2,transparent);border:1px solid var(--border);
      border-radius:6px">
      <div style="display:flex;gap:8px;align-items:baseline;flex-wrap:wrap">
        <b style="font-size:13px">What is lost by blocking it</b>
        <span class="subnote">${esc(c.physiology.severity_label)}</span>
      </div>
      <div class="subnote" style="margin-top:5px"><i>Normally:</i>
        ${esc(c.physiology.signal)}</div>
      <div class="subnote" style="margin-top:5px"><i>Without it:</i>
        ${esc(c.physiology.cost)}</div>
      <div class="subnote" style="margin-top:5px;opacity:.8">Source:
        ${esc(c.physiology.source)}</div>
    </div>` : ''}

    ${(c.notes || []).map((n) => `<p class="subnote" style="margin:8px 0 0"
      >${esc(n)}</p>`).join('')}

    ${(c.variant_hits || []).map((v) => `
      <p style="margin:10px 0 0;font-size:13px;color:${v.hits_epitope
        ? 'var(--critical)' : 'var(--ink-2)'}">
        <span class="mono">${esc(v.variant)}</span> &mdash; ${esc(v.detail)}</p>`
      ).join('')}

    ${(nCon + nUnk + nCle + c.tissue_flags.length + tHarm.length
       + tSpared.length) ? `
      <div class="ctrls">
        <button class="ghost" data-sopen="${esc(c.key)}"
          >${open ? 'Hide' : 'Show'} the evidence</button>
        <a href="#target/${esc(c.key)}" style="font-size:13px">Open target</a>
        ${!c.disqualified ? `<a href="#target/${esc(c.key)}/run"
          class="mono" style="font-size:12.5px">design against it &rsaquo;</a>`
          : ''}
      </div>` : ''}

    ${open ? `
      <div style="margin-top:10px">
        ${nCon ? `<h3 style="margin:10px 0 0">Carried on a population you
          asked to spare</h3>${lineageRows(c.conflicts)}` : ''}
        ${nUnk ? `<h3 style="margin:14px 0 0">Not known for a population you
          asked to spare</h3>${lineageRows(c.unknowns)}` : ''}
        ${nCle ? `<h3 style="margin:14px 0 0">Positively cleared</h3>
          ${lineageRows(c.cleared)}` : ''}
        ${tHarm.length ? `<h3 style="margin:14px 0 0">Damaged outside the
          blood <span class="subnote">established in patients; nothing in a
          flow panel on this specimen can speak to it</span>
          </h3>${lineageRows(tHarm)}` : ''}
        ${tSpared.length ? `<h3 style="margin:14px 0 0">Spared outside the
          blood <span class="subnote">a tissue positively excluded, which is
          not an answer about blood</span>
          </h3>${lineageRows(tSpared)}` : ''}
        ${c.tissue_flags.length ? `<h3 style="margin:14px 0 0">Transcript
          flags outside the blood <span class="subnote">always screened,
          never asked for &mdash; the weakest evidence on this page</span>
          </h3>${lineageRows(c.tissue_flags)}` : ''}
      </div>` : ''}
  </div>`;
}

function spRenderResult() {
  const err = $('sperr');
  if (err) {
    err.innerHTML = SP.error
      ? `<div class="card" style="border-left:3px solid var(--critical)">
         <p style="margin:0;font-size:13.5px">${esc(SP.error)}</p></div>` : '';
  }
  const host = $('spres');
  if (!host) return;
  const r = SP.result;
  if (!r) { host.innerHTML = ''; return; }

  if (!r.n_antigens) {
    host.innerHTML = `
      <div class="card"><p style="margin:0;font-size:13.5px">No antigens
      entered, so there is nothing to rank. The suggestions below come from
      the sparing profile alone and say nothing about whether this specimen
      carries any of them.</p></div>
      ${r.untested.length ? `<h2 style="margin-top:20px">Worth putting on a
        panel</h2>${r.untested.map(specimenCard).join('')}` : ''}`;
    spWire(); return;
  }

  host.innerHTML = `
    <h2 style="margin-top:22px">Candidates
      <span class="subnote">${r.candidates.length}</span></h2>
    ${r.candidates.length
      ? `<p class="secnote">Ordered by how much is positively known about
         them, not by how clean they look.</p>`
         + r.candidates.map(specimenCard).join('')
      : `<div class="card" style="border-left:3px solid var(--serious)">
         <p style="margin:0;font-size:13.5px"><b>Nothing qualified.</b> Every
         antigen entered is either absent from the blasts or carried on a
         population this profile says must survive. That is a real answer
         about this specimen and the targets in this catalogue, not a
         failure of the search &mdash; and it is the situation in which a
         target nobody has built yet matters most.</p></div>`}

    ${r.excluded.length ? `
      <h2 style="margin-top:26px">Excluded
        <span class="subnote">${r.excluded.length}</span></h2>
      <p class="secnote">Kept rather than hidden. Why a target was ruled out
        is usually more useful than the ones that survived.</p>
      ${r.excluded.map(specimenCard).join('')}` : ''}

    ${r.untested.length ? `
      <h2 style="margin-top:26px">Worth adding to the panel
        <span class="subnote">${r.untested.length} of
        ${r.n_untested_considered} considered</span></h2>
      <p class="secnote">Catalogue targets in this disease area that the panel
        entered does not cover, and that the sparing profile does not rule
        out. Whether this specimen carries them is unknown until somebody
        looks.</p>
      ${r.untested.map(specimenCard).join('')}` : ''}`;
  spWire();
}

function spWire() {
  document.querySelectorAll('[data-sopen]').forEach((b) => {
    b.onclick = () => {
      SP.open[b.dataset.sopen] = !SP.open[b.dataset.sopen];
      spRenderResult();
    };
  });
}

/* ---------- boot ---------- */
addEventListener('hashchange', rerender);
addEventListener('resize', () => {
  clearTimeout(window._rt);
  window._rt = setTimeout(() => {
    if (S.job && S.job.state === 'done') renderResults(S.job);
  }, 200);
});
rerender();
