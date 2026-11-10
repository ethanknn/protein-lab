# protein-lab

An in-silico binder-design sandbox aimed at **CLEC12A (CLL-1)**, a surface
antigen on leukaemic blasts in acute myeloid leukaemia.

It generates candidate binder sequences on a real protein scaffold, scores them
against a real experimentally-determined epitope using published biophysical
criteria plus a protein language model, iterates with an evolutionary loop, and
produces a dashboard with 3D structures and the full reasoning behind every
number.

---

## Read this first

**This software cannot discover a cure, and no number it produces predicts
clinical benefit.**

What it does is *design triage*: the question "of these candidates, which should
a laboratory build first?" That is a real question that real campaigns ask, and
answering it well has value. It is not the same question as "does this molecule
help a patient", and conflating the two would be the most harmful thing this
code could do.

A candidate scoring at 95% of the reference means it matches a real anti-CLEC12A
antibody **on the criteria a desktop computer can evaluate** — charge
complementarity, aromatic content, solubility, manufacturability, sequence
plausibility. It does not mean the molecule binds. Establishing that requires:

1. expressing the protein (weeks)
2. measuring binding affinity — SPR or BLI (weeks)
3. testing on CLEC12A-positive cells, then primary patient blasts (months)
4. confirming healthy haematopoietic stem cells are spared (months)
5. animal models (1–2 years)
6. clinical trials (5–10 years)

Most candidates that look this good still fail. Every score in the dashboard
carries a `caveat` field stating what that specific measurement cannot tell you,
and those are displayed, not buried.

**On the target itself:** CLEC12A was chosen because published expression data
shows it on blasts in >90% of AML patients while being largely absent from
healthy haematopoietic stem cells — unlike CD33 and CD123, whose presence on
healthy stem cells causes lasting marrow toxicity. But CAR and antibody
therapies against CD33, CD123 *and* CLEC12A have **not** produced durable
remissions in most AML patients, and antigen escape is a documented resistance
mechanism. A good target is necessary and nowhere near sufficient.

---

## Quick start

```bash
pip install numpy torch fair-esm
python scripts/build_catalog.py    # fetch + verify every disease target
python server.py                   # then open http://127.0.0.1:8800
```

That gives you the web application: browse targets by disease area, inspect
any binding site in 3D, and run a design cycle with live progress.

For the command line against a single target:

```bash
python tests/test_pipeline.py      # 103 checks, including known-good controls
python run.py --no-plm --generations 18 --pop 36   # seconds
python run.py                                      # with ESM-2 on GPU, minutes
```

## The web application

`server.py` serves a single-page app on the standard library alone — no web
framework, because the server needs to do very little and the scientific core
already avoids fragile dependencies. Design runs happen on background threads
and the page polls, so a run that takes minutes shows its search converging
rather than freezing on a spinner.

| Route | What it does |
|---|---|
| `/` | catalogue, grouped by disease area, with search |
| `/api/catalog` | verified targets + the candidates that failed verification |
| `/api/target/<key>` | binding site, surface chemistry, design brief |
| `/api/structure/<pdb>?chains=` | PDB subset for the 3D viewer, converted from mmCIF when that is all there is |
| `POST /api/run` | start a design job |
| `/api/job/<id>` | progress, live trajectory, then results |
| `/api/custom/search?q=` | gene symbol or accession to reviewed human proteins |
| `/api/custom/structures?accession=` | entries holding it with a partner bound |
| `POST /api/custom/prepare` | build a target from an accession + a PDB id |
| `/api/custom/list` | targets built in this session |
| `/api/specimen/setup` | the antigen checklist, the lineages, the sparing presets |
| `POST /api/specimen/rank` | rank the catalogue against one antigen profile |

### What is in the catalogue

Ninety-five verified targets across 19 disease areas, listed by gene symbol:

| Area | Targets |
|---|---|
| Solid tumour | CEACAM5, CLDN18, EGFR, ERBB2, FOLH1, KIT, MET, MSLN, NECTIN4, TACSTD2, VEGFA |
| Autoimmune | FCGRT, IL17A, IL1B, IL23A, IL4R, IL6, IL6R, TNF, TNFSF13B |
| Immuno-oncology | CD274, CD47, CTLA4, HAVCR2, LAG3, NT5E, PDCD1, TIGIT, VSIR |
| Blood cancer | CD19, CD22, CD33, CD38, CLEC12A, GPRC5D, MS4A1, TNFRSF17 |
| Neurology | APP, AQP4, CALCRL, CHRNA1, MAPT, SNCA, SOD1, TREM2 |
| Rare disease | ACVR1, C3, C5, CFB, CFD, FGF23, KLKB1, MASP2 |
| Infectious disease | CSP, F, GP, HA, S, env, tcdB |
| Respiratory | IL11, IL33, IL5, IL5RA, ITGAV, TSLP |
| Metabolic | GCGR, GDF15, GLP1R, KLB, LEPR |
| Dermatology | IL13, IL31RA, IL36G, TNFRSF4 |
| Haematology | C1S, SELP, TFPI, VWF |
| Ophthalmology | ANGPT2, C1QB, IGF1R |
| Allergy | IGHE, SIGLEC8 |
| Bone | SOST, TNFSF11 |
| Cardiometabolic | F11, PCSK9 |
| Gastrointestinal | ITGB7, MADCAM1 |
| Musculoskeletal | MSTN, TGFB1 |
| Transplant | CD40LG, IL2RA |
| Pain | NGF |

Folate receptor alpha is the one deliberate failure; see below.

\* folate receptor alpha is the one deliberate failure; see below.

**Nothing enters the catalogue unverified.** `scripts/build_catalog.py` fetches
each candidate, checks the gene symbol matches the accession, checks a chain in
the structure really is that protein, and derives a binding site from whatever
is bound to it. A candidate that fails is listed in the UI with the reason,
because a tool that silently drops what it cannot do is harder to trust than
one that says so.

One candidate currently fails: **folate receptor alpha**. It is a validated
target with an approved drug against it and a 1.55 Å structure — but every
structure of it in the PDB is unbound, and this pipeline reads binding sites
off a bound partner. It is kept in the catalogue deliberately, as the clearest
statement of where the method stops.

#### The second research pass, and two bugs it exposed

The catalogue went from 42 targets in 11 areas to 74 in 19, filling the gaps
that mattered most: nothing in ophthalmology, metabolic disease,
dermatology, gastroenterology, transplant, pain or non-malignant haematology,
and two lonely entries in neurology.

**Nothing was typed from memory, at either step.** `scripts/resolve_gene.py`
turns a gene symbol into a reviewed human accession and
`scripts/rcsb_lookup.py pick` turns that accession into PDB entries that
really hold the protein with something bound to it. The first earned its keep
immediately: the symbol `SELP` resolves to **selenoprotein P** before it
resolves to P-selectin, and `MET` returns a transcription modulator alongside
the kinase. Either would have verified cleanly and been the wrong protein.
Both lookups now live in `src/plab/lookup.py`, because the build-your-own tab
runs the same two steps for whatever someone types in and a second copy would
be a second chance to get it wrong.

Thirty-two of thirty-two candidates verified. Two bugs had to be fixed first.

**The reader could not open modern structures.** Two candidates were rejected
as "could not fetch PDB" — the IL-31 receptor with nemolizumab (9EZA) and the
glucagon receptor (9XN4). Both structures exist and are fine. Neither has a
`.pdb` file: the legacy format cannot express an entry past 62 chains or
99,999 atoms, and the PDB has stopped generating it for a growing share of
ordinary entries too. The parser now reads mmCIF as well, and an mmCIF-only
entry is re-emitted as PDB records for the viewer. This matters far more for
the build-your-own tab than for the catalogue: someone typing in a gene lands
on recent entries constantly, and "could not fetch" reads like a dead
structure rather than a reader that cannot open it.

**A binding site was read off the wrong protein.** 6Z6V holds all three C1q
subunits. The nanobody binds subunit B. C1QC matches its own chains at 100%
and subunit B at 51% — both above the aligned-identity threshold, so both were
candidates, and the rule was *take the largest interface*. C1QC's binding site
was therefore read off a different gene product. Nothing downstream could have
noticed: the sequence is real, the contacts are real, and the epitope simply
belonged to another protein.

Identity now outranks interface size. The first attempt at that fix banded
identity in fine steps and was too blunt — it demoted chain Z of the PD-1
structure and chain D of the LAG-3 one, which are 90% and 95% identical to
their references because crystallographic copies are partly disordered, not
because they are different proteins. Both lost their larger, better-ordered
interface to a second copy, and two epitopes in the catalogue changed. The
split is now the threshold this module already had for "same protein", so
copies compete on interface size exactly as before and only a genuine relative
is held back. After the fix, **none of the original 42 epitopes changed** —
which is the check that mattered, and the reason to run it.

#### The rejections that were our fault

An earlier build rejected five more candidates than this one. Every one of
those turned out to be a defect here rather than a limit of the data, which is
why the rejection list now gets audited as carefully as the results:

| Candidate | Looked like | Actually was |
|---|---|---|
| HIV gp120 | "13% identity, not the target" | the right protein, a different clade |
| Influenza HA | "26% identity, not the target" | the right protein, 4FQI is H5N1 and the accession was H1N1 |
| RANKL | "42% identity, not the target" | the right protein, 4GIQ is the mouse ortholog |
| IL-6R | "0% identity" | PDB 4O9H holds interleukin-6, not its receptor |
| IL-5 | "0% identity" | PDB 7KQ5 is a bacterial encapsulin; the id was simply wrong |

The first three share one cause. Chain identity was measured as the fraction
of exact 8-mers shared with the reference, and an exact k-mer needs all k
positions identical — so at sequence identity *p* the expected hit rate is
*p*⁸. A true ortholog at 90% identity scores 0.9⁸ = 43%, below a 50%
threshold, and a different viral strain scores far lower still. The measure
was not slightly off; it was measuring the wrong quantity, and it said
"different protein" about proteins that were plainly the same.

Identity is now decided by Smith-Waterman alignment with a coverage
requirement, with the k-mer test kept only as a fast path. The threshold sits
in a wide empirical gap: on this catalogue, chains that are *not* the target
land at 10–22% identity, while orthologs and strain variants land at 56–90%.

The last two share a different cause: **a PDB id recalled from memory is a
guess.** `scripts/rcsb_lookup.py` now asks RCSB which accession each chain
maps to, and which entries hold a given protein with something bound; every
target added since was chosen from that output rather than from recall.

#### A binding site on the wrong side of the membrane

One entry verified cleanly and was still wrong. The glucagon receptor was
paired with 9XN4, a receptor–Gs complex, and the pipeline dutifully derived a
36-residue binding site — from the **G protein**, which grips the receptor's
*intracellular* face. Every automatic check passed. The chain was the right
protein, something was certainly bound to it, the interface was large and
well resolved. Nothing asked which side of the membrane it was on, and an
injected protein binder can never reach that surface.

The entry now uses 4LF3, an allosteric antagonist Fab on the extracellular
face, where its reachable selectivity margin is 0.443 against a ceiling noise
of 0.029 — a strong, usable target.

I caught it by eye, which is not a method, so `derive_epitope` now recognises
the sequence fingerprints of G proteins and warns. Getting that right took
the same care as everything else here: 9XN4's G-alpha carries
`GAGESGK`**`N`**`T`, a dominant-negative P-loop mutation used to trap the
nucleotide-free state, so the wild-type motif misses it completely. Arrestin
is deliberately *not* in the list — I do not know its conserved core well
enough to write it down, and a motif that is nearly right would fail silently
in exactly the way the thing it is meant to catch does.

#### The third pass: "no usable structure" was usually about us

The second pass left a list of targets I had written off — the Protein Data
Bank, I had concluded, held no structure of them with anything bound. Going
back over that list with a deeper search recovered most of it.

Structure searches return entries sorted by resolution. A protein that has
been crystallised a hundred times on its own therefore buries its complexes
behind all of them, and the scan stopped at forty.

| Target | What a shallow scan said | What is actually there |
|---|---|---|
| SOD1 (ALS) | no usable structure | **eleven**, the best a nanobody complex at 2.19 Å |
| PSMA (prostate cancer) | a three-residue dipeptide | five, including nanobody complexes |
| MET (lung cancer) | a 13-mer peptide bound to a ubiquitin ligase | seventeen, including onartuzumab and amivantamab |

Each of those is a disease target of real consequence, and each had been
recorded as impossible. Worse, the same limit sat behind the build-your-own
tab: typing `SOD1` told you your protein could not be worked on.

The scan now goes 250 entries deep, batched, which costs about a second per
protein. The regression test stubs the network and buries one usable complex
behind ninety-seven apo structures.

**What that leaves is a real answer.** After scanning 250–300 entries each,
these genuinely have nothing: transthyretin (300 scanned, every one either
alone or with a small-molecule stabiliser), folate receptor alpha, ANGPTL3,
CD52, DLL3, CCN2 and the PLA2 receptor. When the tool now says no usable
structure exists, it means it.

#### What the third pass added

Twenty-one candidates, filling the areas that were still thin and adding two
diseases whose *autoantigen* is known exactly:

- **Neurology** — SOD1 (ALS), aquaporin-4 (neuromyelitis optica), the
  acetylcholine receptor (myasthenia gravis). The last two are unusual: the
  bound molecule in each structure is a **patient's own pathogenic
  antibody**, so the binding site is literally the site that causes the
  disease.
- **Rare disease** — plasma kallikrein with lanadelumab's own Fab, ALK2
  (fibrodysplasia ossificans progressiva), FGF23, factor B, MASP-2.
- **Haematology** — von Willebrand factor, C1s.
- **Autoimmune** — FcRn with nipocalimab's Fab, the mechanism that strips
  pathogenic antibodies out of circulation whatever they target.
- **Solid tumour / immuno-oncology** — PSMA, MET, CD73, VISTA.
- **Respiratory / fibrosis** — integrin β6, IL-5Rα, IL-11, TGF-β1.
- **Allergy and bone** — Siglec-8 with lirentelimab's Fab, sclerostin.

One of the twenty-one did not survive, and usefully. The integrin entry was
written against **β6** and rejected: *no bound partner matching 'antibody'
makes contact with the target*. The 17E6 Fab in that structure binds the
**αV** subunit, not β6 — the same mistake as the C1q entry, caught this time
by the guard rather than by eye. Re-pointed at αV it gives a clean
22-residue site with no warnings at all.

The panels the new targets produced are the best check the method has had.
FcRn lands on the HLA class I molecules — it *is* an MHC class I fold, which
no part of this pipeline was told. Plasma kallikrein lands on factor XI at
67%, and both are now catalogue entries in their own right. Siglec-8 lands on
CD33 at 64%, so two catalogue targets are each other's cross-reactivity risk.
Aquaporin-4 lands on the other aquaporins, the acetylcholine receptor on its
own sibling subunits, factor B on C2, sclerostin on its single paralog
SOSTDC1.

The ceiling guard then withheld five of the fifteen new scores, including
TGF-β1 — whose nearest relative TGF-β3 is 67% identical at the binding site,
so no composition can separate them.

The catalogue now holds **95 verified targets across 19 disease areas**, with
**53** carrying a cross-reactivity panel and a measured selectivity ceiling.

### Design against any protein

The catalogue is 95 proteins somebody curated. The **Build your own** tab is
the same pipeline pointed at anything in UniProt: resolve the name, list the
structures that hold it with something bound, read the binding site off the
one you choose, derive the cross-reactivity panel from its domain family,
measure what selectivity is reachable, and run the design loop — all scored by
exactly the same code that scores a catalogue target.

Typing `KRAS` and picking 9E3S takes about a minute and produces a target
whose panel is the RAS superfamily, ranked by identity *at the binding site*:
RRAS2 at 69%, RRAS at 66%, MRAS at 64%. The selectivity score is then withheld
with the reason given — the KRAS patch is `EYDPIEDEY` and RRAS2's is
`DYDPIEDEF`, seven of nine identical, and no composition can separate them.
That is the honest answer and the machinery arrives at it unaided.

Three things the tab is careful about:

- **It cannot invent a binding site.** It needs an experimental structure of
  the protein with something bound. Most human proteins have none — typing
  `SOD1` returns sixty structures and no usable one, because every entry is
  the enzyme alone. The page says that is a fact about what has been solved
  rather than a verdict on the protein, because "no results" would read as
  the software failing.
- **A built target is not a catalogue target.** Catalogue entries were checked
  by a person and carry a written rationale and a clinical caveat. A built one
  carries neither, and its page says so in those words rather than leaving an
  empty space that reads as "no concerns". Everything automatic — the chain
  really is this protein, something really is bound to it, the interface is
  big enough — is the same code the catalogue build runs, and its warnings are
  passed through untouched.
- **It is not written into the catalogue.** `data/catalog/` is a reviewed
  file. Targets built in a session live in memory and disappear when the
  server stops, which the page also says.

### Which target, for this specimen

The selectivity scorer answers one question: could a binder confuse this
target with a structurally *related* protein? Its own caveat names what it
cannot do — whether an off-target protein sits on a cell that matters is
measured by cytometry, and it is why CLEC12A was chosen over CD33, two
proteins that share no domain family at all. No paralog panel will ever put
them in the same room.

The **Specimen profile** tab is that second axis. Enter what the flow
cytometry says the blasts carry and which healthy populations have to
survive, and it ranks the catalogue against that profile. For a specimen
whose marrow must recover, CD33 is excluded — carried on committed
progenitors, which is the expression behind gemtuzumab's myelosuppression —
while CLEC12A is not, and the tool names the paper that settles it. That is
target selection, not treatment selection, and the page says so.

Two tables sit behind it: a lineage table for blood populations, and a tissue
table for everything else. Which one carries the argument depends entirely on
the target — for a blood cancer antigen it is the first, and for a solid
tumour antigen it is the second, as the page has to say out loud.

#### Transcript abundance may raise a concern. It may never clear one.

This is the rule the whole feature is built on, and it came out of checking
rather than caution. GPRC5D is a myeloma target; in the clinic 75% of
patients on talquetamab get nail changes and 63% get oral toxicity, from
on-target expression in the nail matrix, hair cortex and the filiform
papillae of the tongue. Ask the Human Protein Atlas and the answer is split:
**26.4 nTPM in skin**, which is a real flag, against **0.0 in tongue**.

The screen sees part of the toxicity and is blind to the rest, which is more
instructive than a clean miss. A column reading zero cannot be read as
clearance even for a protein the same screen flags elsewhere — there is no
nail matrix in a bulk tongue sample. So there is deliberately **no nTPM value
that returns "absent"**; a population with no transcript and no cytometry is
reported *unknown*, and only cytometry or observed human toxicity can
establish absence.

An earlier draft of this claimed the screen missed GPRC5D entirely. That was
wrong for an instructive reason: the probe behind it asked for `t_RNA_skin`,
which is not a column code, and that API answers an unknown column by
silently leaving it out — a typo therefore reads downstream as "not expressed
in skin". `scripts/fetch_lineage.py` now checks every requested column
against what came back, and a column missing for *every* gene is reported as
a bad column code rather than as biology. The same guard catches the other
trap: that search is fuzzy, and asking it for `CD33` returns fourteen genes
with CD33 first only by luck, so every row is matched on the exact symbol.

#### Ranking by how much is known, not by how clean it looks

A naive ranking sorts by concerns found, ascending — which puts the protein
nobody has ever examined at the top of every list, because an unexamined
target has no concerns on record. That is the nTPM mistake promoted from one
cell of the table to the ordering over the whole table, where it does more
damage.

So a lineage positively **cleared by cytometry counts in a target's favour**,
and a lineage nobody has data for counts against it. Two targets with no
recorded concerns are not equivalent: one has evidence behind it and the
other is mostly an absence of looking. Disqualification stays separate from
ranking, as it is in the design scores — an antigen the blasts do not carry
is a non-candidate, not a weak one, and no amount of selectivity elsewhere
compensates.

#### Expression is not harm, and CD47 is why

The first version of this model disqualified a target whenever it was densely
carried on a population the clinician asked to spare. That is wrong, and
CD47 is the counterexample that forced the fix: it sits on haematopoietic
stem cells — it is the "do not eat me" signal on essentially every cell — and
blocking it still produces **no increase in phagocytosis of normal stem or
progenitor cells**, because the window there is a density difference rather
than an absence. Reading "expressed here" as "destroys this" would have ruled
out one of the few targets where human and in vivo experience says otherwise.

So every curated cell now carries a **harm** verdict alongside its level:
*damage seen in patients*, *population survives it*, or unknown. Expression
on a must-spare lineage disqualifies unless something shows the population
survives it, and "nobody has checked" does not earn that exemption. On one
CD47 card all three stem, erythroid and platelet rows read *carried, high*
and are separated only by what happened to patients — the stem row is
exempted with a citation, the erythroid row is the on-target anaemia, and the
platelet row is excluded precisely because nothing reassures.

#### The target that was perfectly selective and still failed

TIM-3 is the most instructive entry in the table. It discriminates leukaemic
stem cells from normal haematopoietic stem cells **inside the same
CD34+CD38− gate** — as clean a therapeutic window as this axis can describe,
and the tool duly ranks it as a candidate where CD33 is excluded.

The phase 3 STIMULUS-MDS2 trial of sabatolimab plus azacitidine then missed
its primary overall-survival endpoint, and the programme was discontinued in
2024.

That pairing is recorded on the TIM-3 stem-cell row on purpose. Everything
this software measures can come out right — the structure is real, the
epitope is measured, the selectivity margin is positive, the healthy-cell
evidence is cytometric and clean — and the drug can still not help anybody.
Design triage narrows which molecule is worth building. It does not make the
molecule work.

Three further things the tab is careful about:

- **Stem cells and committed progenitors are separate rows.** They were one
  bucket in the first draft, which would have been wrong: CD33 is largely
  absent from primitive CD34+CD38− cells but present on committed
  progenitors, and collapsing the two would have made the contrast with
  CLEC12A — absent from both — disappear into a single misleading cell.
  Neither population appears in the atlas at all, whose bone marrow
  single-cell set is 3,484 cells with no stem or progenitor cluster in it, so
  both rows carry a visible *no screen data* marker.
- **Non-blood tissue is always screened, never asked for.** A blood-only
  model scores GPRC5D as beautifully restricted and says nothing about the
  nails of three patients in four. The clinician does not know to ask.
- **Seven catalogue targets are viral or bacterial**, and the atlas correctly
  returns nothing for them. That renders as "no human cell carries this
  protein, so there is no lineage to spare" rather than as eleven unknowns —
  the same distinction the paralog panels draw with `no_human_relatives`.

There are no patient identifiers in the form and nothing is stored.

#### Blood is not the only axis, and for most of the catalogue it is the wrong one

This feature started as a blood story, because the antigens that motivated it
were blood cancer antigens. Extending it to the twelve solid tumour and nine
immuno-oncology targets moved the weight off the lineage table entirely.

HER2's established liability is cardiomyocytes. Nectin-4's is skin.
Mesothelin's is the lining of the pleura, the pericardium and the peritoneum.
PSMA's is the salivary acini. **None of those is a blood population**, and no
flow panel on a marrow specimen can address any of them. So one of these
targets can come back from the sparing check with no conflicts at all and be
the most dangerous thing on the list.

That is the worst failure available to this page, because it fails *quietly*
and in the reassuring direction: an empty conflict list reads as a clean
target. The ranking already keeps these antigens off the top — unknowns count
against them, so HER2 cannot outrank CLEC12A on a marrow profile — but a rank
is not an explanation, and a clinician reading one card needs the
explanation. So the card now says it:

> ERBB2 has an established liability outside blood: heart muscle. Nothing in
> a flow panel on this specimen can speak to it, and the sparing profile
> above does not try to — every population it offers is haematopoietic. Read
> the tissue rows as a separate question that this specimen cannot answer,
> not as a milder version of the rows above.

And when nothing is curated for a solid tumour target at all, it says *that*
instead, rather than letting silence pass for a result: "Nothing is curated
for MET outside blood, which means nobody has recorded it here — not that the
tissue is clear."

The tissue table is now keyed by `(gene, site)` rather than by gene. One cell
per gene could hold one note, and these antigens do not have one liability
each: PSMA reaches salivary gland, kidney *and* duodenum, and EGFR's rash and
its hypomagnesaemia are different organs, different mechanisms and different
severities. Collapsing them would repeat the first draft's single
stem/progenitor bucket in a new place.

#### Three ways a transcript screen misleads, not one

The GPRC5D false negative (above) was the first. The solid tumour antigens
produced two more, and the second is worse.

**The flag can be on the wrong isoform.** CLDN18 has two splice variants:
CLDN18.1 in lung alveolar epithelium, CLDN18.2 in differentiated gastric
mucosa. Zolbetuximab binds the first extracellular domain of 18.2 and
explicitly does not bind 18.1. The screen reads the *gene*. So it was
reporting a loud 237.8 nTPM flag in lung — the one isoform the drug is
specific against — and saying nothing whatever about the stomach, because
there was no stomach column. There is now, and it reads **684.2 nTPM**, which
is where the emesis in 76% of patients comes from. A screen can be precise,
well-formed, and pointed at the wrong molecule.

Adding the column fixes the blind spot. The isoform blindness is not fixable
at transcript level and remains true of every automated row on the page.

**The flag can be right about everything and informative about nothing.**
VEGFA is detected in all seventeen screened tissues — thyroid 255.6, heart
202.6, liver 202.0, kidney 105.1. All true. The deeper mismatch is that VEGFA
is *secreted*: "which healthy cells carry it" is the wrong question for a
soluble ligand, and the right one — what normal physiology depends on the
signal — is not an expression table at all. Roughly two dozen catalogue
targets are secreted and this axis mis-frames every one of them. That is
recorded as a limit, not papered over.

Five columns were added, each for a named reason rather than for coverage:
stomach for CLDN18, duodenum and pancreas because PSMA, CEACAM5 and
mesothelin all reach the upper GI tract, urinary bladder for nectin-4, breast
for HER2 and KIT. The duodenum column paid for itself immediately: **PSMA
reads 250.9 nTPM there**, its highest value in any screened tissue, above
salivary gland at 49.8. The twelve-column screen could not see it, because
the panel stopped at colon and oesophagus. The lesson is about the screen
rather than about PSMA — a tissue panel is a list somebody chose, and a
gene's loudest site can sit just outside it and read as silence.

HER2 is the counter-case worth keeping in view: heart muscle reads 53.9 nTPM,
in the top band, with no curation needed to raise the flag. The screen is not
useless and it is not sufficient, and **you cannot tell which case you are in
from the number alone.**

#### What no expression table can see

The most honest finding in this pass is an absence.

Mirvetuximab soravtansine's dose-limiting toxicity is corneal — keratopathy
or blurred vision in 50% of patients — and FRα is *absent from human corneal
tissue*. The ocular injury appears across ADCs irrespective of where their
targets are expressed. Sacituzumab govitecan's dose-limiting neutropenia and
diarrhoea are SN-38 effects, the same ones irinotecan produces with no
antibody at all, and the cynomolgus studies found no histopathological change
in the normal tissues that do carry TROP2. Telisotuzumab vedotin's peripheral
neuropathy in 45% of patients is MMAE.

In each case this axis returns a clean answer and is *right about the thing it
measures*. The healthy-cell question is one input to deciding what to build.
It was never going to be the one that decides whether a molecule is
tolerable.

TROP2 makes the point sharpest. It draws more transcript flags than anything
else in the catalogue — oesophagus 1245, salivary gland 585, skin 555,
urinary bladder 219, kidney 186, breast 127, pancreas 107 — which by the
logic of the rest of this page is the worst profile here. And its
dose-limiting toxicities come from the payload.

#### The Fc decides, not the antigen

The checkpoint receptors exposed a flaw in the `harm` field. PD-L1 sits on
monocytes and dendritic cells. Avelumab, a wild-type IgG1, kills
PD-L1-positive dendritic cells by NK-mediated ADCC. Atezolizumab and
durvalumab, which carry an engineered Fc, do not lyse the same cells in the
same assay. **Same target, same population, opposite answers, decided by the
format.**

So `harm` is not a property of the antigen and the cell. It is a property of
the antigen, the cell and the molecule, and this table records the verdict for
the format that was actually tested. Every entry whose answer depends on that
now says so. PD-1 on T cells reads "population survives it" *only* because
approved anti-PD-1 antibodies were deliberately built on Fc backbones that do
not recruit killing — an ADC or a T-cell engager against PD-1 would be
depleting the cells it needs, and nothing on this page would tell you.

The one thing this field must never do is let a blocking antibody's safety
record be read as clearance for a cytotoxic payload aimed at the same antigen.

Two evidence tiers were added for the same reason — to stop a finding wearing
a provenance it does not have. `ihc` because the tissue rows are read by
immunohistochemistry, and filing nail-matrix or gastric-mucosa staining under
"flow cytometry on primary human cells" would be a lie about how a thing is
known. `animal` because TROP2's "the tissue survives" verdict rests on
cynomolgus histopathology, which is real in vivo evidence and is not a human
measurement — so it may soften a flag and may never clear a human population
or exclude a target. Both thresholds now read the evidence ranking instead of
hard-coding a pair of strings.

#### KIT, the exception that proves the rule

Every other solid tumour antigen has its liability in epithelium. KIT has it
in the stem compartment, and it is the worst case in the table.

CD117 is a defining marker of haematopoietic stem and progenitor cells, and
the evidence here is unusually direct: **briquilimab is an anti-CD117 antibody
developed specifically as transplant conditioning**, because blocking stem
cell factor depletes the stem compartment. A drug whose entire purpose is to
empty the marrow niche settles the question of whether engaging this antigen
reaches that population. So a GIST target is now excluded by a
marrow-recovery profile on exactly the same terms as CD33.

Barzolvolimab's neutropenia and hair depigmentation are the same expression
seen from the clinic — and the depigmentation is the reminder that KIT's other
normal homes, mast cells, melanocytes, germ cells and the interstitial cells
of Cajal, are populations **no row of this table represents at all.**

#### The stem-compartment argument, in epithelium

CLDN18.2 is the reason to take the whole lineage/stem distinction seriously
outside blood. In normal stomach it sits at the apical side of tight
junctions in differentiated epithelial cells and is **absent from the gastric
stem cell zone** — a direct echo of the CLEC12A argument in a different
organ. That is why the target is tolerable at all: the mucosa regenerates, so
the emesis is severe in cycle 1, under 10% by cycle 2, and self-limiting.
Sparing the stem compartment is what makes an on-target toxicity recoverable,
in epithelium exactly as in marrow.

It earns no rank credit, deliberately. It is the best single fact in the
tissue table and it is still an answer about epithelium to a question about
blood.

#### What the curated table does and does not cover

Forty-three lineage cells across sixteen antigens, with the stem-cell row
filled for eleven, plus fourteen tissue cells across ten antigens. Against a
grid of 11 lineages × 95 targets = 1,045 lineage cells, that is a small
fraction, and everything else honestly reads *unknown*. The table is
deliberately small: an entry invented to fill a gap would be worse than the
gap, which at least renders as "not assessed".

Two gaps are left open on purpose, because looking is what found them:

- **GPRC5D has no stem-cell row.** Nothing found tested it on normal
  haematopoietic stem cells. The stem-cell exclusion that is easy to
  misattribute here — "restricted to plasma cells, absent from HSC" —
  belongs to **BCMA**, a different myeloma target. What *is* supported for
  GPRC5D is a ratio: 500- to 1000-fold more mRNA in marrow plasma cells than
  in blood B cells. That is a large difference, not an absence, and it is
  recorded as *low* rather than *absent*.
- **CD22 is not quite "B-lineage restricted".** CD22 appears *before* CD19 in
  development, on CD34+CD19−CD22+ cells. Those are B-committed rather than
  multipotent, so the stem compartment is still spared — but a CD34+
  progenitor population does carry the antigen, and a sparing argument
  resting on the tidy summary is resting on something imprecise.

- **CEACAM5 claims no granulocyte row.** The granulocyte expression usually
  attached to "CEA" belongs to **CEACAM8 and CEACAM3**, which are
  granulocyte-restricted, and to CEACAM6; CEACAM5 itself is mostly
  epithelial. Nor does the apical confinement that would protect a
  circulating binder belong to it — that is **CEACAM7**. The evidence names
  its relatives, not it, so the cell stays unknown. This is also the one
  catalogue target where this page and the paralog panel are both live at
  once: the CEACAMs really are close relatives, so confusability is a genuine
  risk here in a way it never was for CLEC12A against CD33.
- **CTLA-4's harm verdict is left unknown**, for a specific reason rather than
  for want of looking. In mice, anti-CTLA-4 depends on Fc-receptor-bearing
  macrophages to deplete intratumoral Tregs. In human tumours, a study looking
  for the same depletion did not find it, and the inhibitory receptor
  FcγRIIB appears to blunt it. Two species disagree, and flattening that into
  either answer would invent a result.
- **PSMA's kidney row records a concern with harm unknown.** Proximal tubule
  uptake is dose-limiting for radioligands; frank nephrotoxicity is described
  as a potential rather than an observed outcome. Overstating it would be
  inventing a result in the direction that feels safer.

One entry cuts against the tool's own favourite example. **CLEC12A is carried
on granulocytes and monocytes**, so a profile that needs mature myeloid cells
excludes it. Its selectivity is against the *stem* compartment, not against
myeloid cells generally, and stating that plainly matters more than keeping
the headline example clean.

One entry reports a source weaker than the rest. The **VISTA** rows rest on
flow cytometry disclosed across one family of anti-VISTA antibody patents —
near-identical disclosures, not independent peer-reviewed work. The direction
is corroborated by the atlas's own single-cell classification and by the first
anti-VISTA antibody stopping for cytokine release, which is what engaging a
dense myeloid antigen looks like. It is reported with the provenance stated
rather than omitted, because a neutropenic patient's profile needs to see it.

#### The question was wrong for half the catalogue, and it answered anyway

The sections above are about *which* cells carry a target. This one is about
a prior question nobody had asked: whether the target is on a cell at all.

The lineage screen reported **alpha-synuclein at 103.2 nTPM on dendritic
cells** and listed them among the populations at risk, in the same words, the
same pills and the same consequence text as CD33 on progenitors. SNCA is a
cytosolic protein. No antibody reaches the cytosol. The number was a
measurement of what those cells *make*, rendered as a claim about what a
binder would *meet* — and the consequence line told a clinician their
dendritic cells were threatened by a molecule that never touches one.

TNF read 73.4 nTPM on monocytes, which are the cells that produce it. IL-6,
VEGFA, PCSK9, C5 and about thirty more sat in the same position.

So the fix is a topology axis, and the one place it must not come from is
recall — three of this catalogue's early rejections were PDB ids remembered
as the wrong protein, and a remembered subcellular location is the same
mistake where nothing downstream would catch it. `scripts/fetch_topology.py`
pulls the annotation from UniProt for all 95 accessions. It stores
**evidence, not a verdict**: keywords, subcellular locations and
transmembrane feature counts. The rule that turns those into a classification
lives in `lineage.classify_topology`, where it can be read, tested and argued
with. A verdict frozen into a data file is a verdict nobody can check.

Four answers, because the honest handling differs for each:

| | what it means | what the page does |
|---|---|---|
| **surface** (44) | displayed where a binder meets it | the lineage question works as asked |
| **shed** (18) | surface *and* a soluble pool | lineage rows kept, plus a decoy warning |
| **secreted** (32) | never on a cell | one honest row instead of eleven wrong ones |
| **intracellular** (1) | inside, out of reach | same, plus the reachability problem |

A secreted or intracellular target now gets a single row saying the question
does not apply — exactly as the viral entries already did — rather than
eleven rows answering a question nobody asked.

**The surface test is structural, and that matters.** It asks for a
membrane-spanning feature or a GPI anchor, not for the `Cell membrane`
keyword. UniProt gives that keyword to proteins which merely *associate* with
a membrane, and the first version of this rule used it — which put tau, a
microtubule-binding protein of the cytosol, in the same bucket as CLEC12A.
Insertion and association are different things and only one of them puts a
protein where a binder can reach it.

Two overrides remain, both for the same reason, and the reason is more
interesting than the correction. UniProt records a genuine minor secreted
pool for **SNCA** and **MAPT**, so the automated verdict is not a bug in the
data — it is dominant biology outvoted by a true footnote. And that minor
pool *is* the drug rationale: anti-tau and anti-synuclein antibodies are
aimed at the small extracellular fraction that propagates between neurons,
precisely because the large intracellular pool is beyond reach. That is a
hard brief, and the entry should say so rather than let "secreted" imply an
accessible target.

**TNF is the case that would have been got wrong by lumping the cytokines
together.** It is made as a transmembrane protein, so monocytes really do
display it, and anti-TNF antibodies really do bind that form — infliximab
and adalimumab bind both soluble and membrane TNF, exert complement-dependent
cytotoxicity against transmembrane-TNF-bearing cells where etanercept is
considerably weaker, and induce outside-to-inside signalling through it.
Those authors suggest this may explain why the antibodies work in Crohn's
disease and etanercept does not; a Leiden review says plainly that it is
still unclear why. The hedging is theirs and is kept. The point for this page
is that TNF's lineage rows are *correct* and must survive, and only a
measured topology distinguishes it from IL-6 three rows away.

The shed class earns one extra line, because it changes what a binder has to
do rather than which cells are at risk: **a circulating pool of the same
antigen is a sink.** It binds the molecule before it reaches a cell. Soluble
mesothelin, soluble PD-L1 and soluble TREM2 are all measured as serum
biomarkers for exactly that reason.

Two limits stated rather than fixed. An off-axis target accumulates no
unknowns, so it is not penalised in the ranking the way an unexamined surface
target is — consistent with how the viral entries have always been treated,
and still a position the ranking does not really earn. And the classifier
only sees what UniProt keywords: BCMA sheds a soluble form that is a known
sink for BCMA bispecifics, and it comes back `surface` here because the entry
carries no `Secreted` keyword.

#### Finishing the sentence: what is lost when the signal is blocked

The topology axis left a half-sentence on every secreted card. It said that
binding such a target *neutralises a signal rather than depleting a
population, so the risk moves to whatever normal physiology depends on that
signal* — and then stopped, because no expression table can reach that.

Thirty entries now finish it, covering twenty-nine of the thirty-two
off-axis targets. Each records what the molecule normally does, what is lost
without it, the source, and how bad: `boxed_warning` (6), `serious` (8),
`manageable` (5), `well_tolerated` (7), `contested` (1) or `unknown` (3).

**The strongest evidence here is not a trial.** For several of these the
human knockout already exists, and it ran first:

- People born with **IL-17 pathway defects** develop chronic mucocutaneous
  candidiasis. Secukinumab and ixekizumab cause oral candidiasis at
  1.7–4.0% against 0.3% on placebo. The inborn error predicted the adverse
  event.
- People with **congenital terminal complement deficiency** get recurrent
  neisserial infection. Eculizumab carries a boxed warning for fatal
  meningococcal disease, a REMS, and — per CDC — a risk up to **2,000
  times** that of an otherwise healthy person. Vaccination reduces it and
  does not remove it.
- People with **congenital factor XI deficiency** largely do not bleed. That
  observation is the entire thesis of the factor XI drug class.

**NGF is the purest case in the table.** Tanezumab blocks nerve growth
factor, nociception falls, and patients load damaged joints they would
otherwise protect. Rapidly progressive osteoarthritis followed; of those
reaching the severe form, around 60% went on to total joint replacement, and
the FDA found that stopping the drug *does not appear to prevent further
damage*. An advisory committee voted 19–1 against. The drug worked exactly
as designed, and the loss of the signal was the harm.

**IL-6 shows a second-order version of the same thing.** Blocking it causes
lower gastrointestinal perforation at roughly 1–2 per 1,000 patient-years,
concentrated in patients with prior diverticulitis. But it also abolishes the
CRP rise and blunts the pain, so those patients present with mild
non-specific symptoms instead of an acute abdomen — and the laboratory
marker a clinician would reach for has been removed by the drug itself. Part
of the harm is that the warning signs are gone.

**Two controls, and they matter as much as the boxed entries.** A table where
everything is dangerous carries no information. PCSK9 homozygous
loss-of-function carriers live with LDL-C below 15 mg/dL showing no
neurocognitive deficit, diabetes, cataract or stroke; regulators required no
special monitoring for the class. And factor XI has a twist that earns it the
best row here: abelacimab produced a large reduction in major bleeding
against rivaroxaban, and the asundexian phase 3 **OCEANIC-AF trial was
stopped early for superior safety and failed non-inferior efficacy.**
Separating a signal from haemostasis may mean separating it from the clotting
you wanted to stop. The safety and the failure may be one fact seen twice —
the same shape as TIM-3 elsewhere in this README.

**One entry is recorded as contested rather than resolved.** Romosozumab
carries a boxed warning for myocardial infarction, stroke and cardiovascular
death. The imbalance appeared in ARCH against alendronate (2.5% vs 1.9% MACE
in year one) and **not** in the placebo-controlled FRAME trial, which leaves
open whether ARCH measured alendronate protecting the comparator arm rather
than romosozumab harming its own. Recorded as boxed because it is, and as
contested because it is that too.

#### The guard that paid for itself twice

Four complement entries carry boxed warnings for encapsulated bacteria — C3,
factor B, factor D and C5 — so the tempting move was to extend the class
warning across the family. C1s, C1q and MASP2 were held back instead, on the
grounds that they sit in the classical and lectin pathways and leave the
alternative pathway running: the CEACAM5 granulocyte mistake in a new
costume, where the evidence names relatives rather than the thing itself.

Researching them properly turned up something better than a correction.
**Those three have a different answer, not a weaker version of the same
one.** Blocking early classical complement points at **autoimmunity**, not at
encapsulated-organism sepsis:

- Hereditary **C1q** deficiency is among the strongest monogenic causes of
  lupus known — one review puts SLE prevalence in C1q-deficient individuals
  near 90%. C1q opsonises apoptotic cells and immune complexes so they are
  cleared before they can be presented as self-antigen; without it, they
  aren't.
- **Sutimlimab**'s own label carries encapsulated-bacteria infection under
  *Warnings and Precautions* — not a boxed warning, unlike the alternative
  and terminal agents — and then warns that, based on its mechanism, the
  drug **may increase the risk of developing autoimmune disease such as
  SLE**. A regulator naming the inborn error's phenotype as a predicted drug
  risk is the clearest possible endorsement of the reasoning this whole
  table rests on.

Had the class warning been extended across the family in the first pass, the
distinction would have been erased before it was found.

**MASP2 got the table's one `contested` row**, which exists for the same
reason `LEVELS` has "contested": it is a real answer, not a hedge. MASP-2
deficiency occurs in about 4% of Caucasians and up to 18% of some African
populations, healthy homozygotes for the classic variant are documented, and
the penetrance has been openly doubted — while the disease database still
lists autoimmunity and recurrent severe infection. A deficiency carried by
one person in twenty-five without obvious consequence is the strongest
argument that a pathway is redundant, and it is not the only published
answer. "Contested" is also distinct from "unknown": one means the evidence
disagrees with itself, the other means nobody looked.

#### Two more entries that cut against the shape of the table

**TGF-beta is the sharpest.** It restrains growth in normal epithelium — a
tumour suppressor — and shields established tumours from immune attack. Same
molecule, opposite jobs. In a phase 1 study of fresolimumab in 29 patients,
four developed keratoacanthoma-type lesions and one of those a
well-differentiated squamous cell carcinoma; they resolved after treatment
stopped. **An anti-cancer drug that causes cancer** is the clearest statement
available of what it means to block a signal with more than one job. (Cardiac
lesions are raised as a concern for chronic TGF-beta blockade in secondary
sources and were *not* verified here, so the entry names them as unconfirmed
rather than recording them.)

**IL-23 against IL-17 is the useful comparison**, and the table now holds
both. IL-23 sustains the Th17 response that produces the IL-17 defending
mucosa against Candida — yet the selective p19 inhibitors show no increase
in opportunistic infection, TB reactivation, oral candidiasis or IBD, while
IL-17A blockade causes candidiasis at several times placebo. Blocking the
upstream driver of a pathway is evidently not the same as blocking its
effector.

A third guard fired while writing these. The reassuring dataset for IL-33 is
**astegolimab, which blocks ST2 — the receptor, not the ligand**. Recording
the receptor's clean safety record under the ligand's name would have been
the same borrowing, so the IL-33 row is graded `unknown` and says why.

#### The pattern in the controls

Seven entries are graded `well_tolerated`, and four of them — **factor XI,
myostatin, tau and alpha-synuclein** — share something the table did not set
out to show: *the drug also didn't work*.

- Factor XI: superior safety, failed non-inferior efficacy (OCEANIC-AF).
- Myostatin: bimagrumab, domagrozumab, LY2495655 and ACE-083 all added muscle
  volume with little or no gain in strength; Pfizer ended the Duchenne
  programme after a 166-child trial showed no effect on mass or strength.
- Tau: gosuranemab missed in TANGO, tilavonemab was well tolerated in 453
  patients and its authors concluded further investigation was not warranted.
- Alpha-synuclein: cinpanemab's SPARK trial was stopped early for futility;
  prasinezumab's PASADENA missed its primary endpoint with a favourable
  safety profile.

None of these failed on toxicity. **A signal you can remove without
consequence may be a signal that wasn't doing much** — which is an
uncomfortable thing for a tool that ranks targets by how safely they can be
engaged, and belongs in the README for that reason rather than despite it.

#### The three left blank, and why that is the point

IL11, IL13 and SOD1 have no entry. Searching turned up nothing solid for
them, and **"my search found nothing" is a fact about the search, not about
the world.** Writing "no safety consequence is established" on the strength
of it would be precisely the nTPM-0.0-means-absent error this package is
built around, committed one layer up from where it was first caught.

So they stay empty, and the card says that half of the question is open
rather than going quiet. None of this touches the ranking: the clinician
asked which populations must survive, not what the molecule is for.

### Scope, honestly

The method needs a protein target with a solved structure and ideally
something already bound to it. That rules out diseases with no single protein
driver — most psychiatric illness, most polygenic conditions — and any target
with no structure. What remains is a few thousand proteins, which is a large
and useful space but is not "any illness".

Adding a target means adding a candidate to `src/plab/catalog.py` and
re-running the build. The scorers never change: they read the binding site's
chemistry out of the structure rather than knowing anything about a particular
protein. HER2's epitope is acidic (net −3.14) and the brief asks for a basic
binding face; TIGIT's is acidic too; CLEC12A's is mixed at +1.00 with three
arginines on it, and the brief says so and asks for shape and hydrogen
bonding rather than net charge. Every one of those is derived from the
structure, not written.

**Three gaps worth naming:**

- Selectivity is scored for 53 of 95 targets. For the rest the component is
  *not assessed* and says which of six reasons applies — these are not the
  same claim and must not render alike. See below.
- Six targets come from structures above 3.0 Å, where side-chain positions
  are not reliably placed and a 4.5 Å contact test is therefore approximate.
  CLEC12A, the target this project started from, is one of them at 3.5 Å.
- Several targets are short peptides rather than folded domains — PfCSP at
  eight interface residues, amyloid-β at six. There is very little surface
  there for a binder to grip, and correspondingly little for the scores to
  discriminate on.

All three are shown on the target page under *Limits of this binding site*,
alongside anything else the catalogue build noticed without rejecting over.

---

## How it works

### 1. The target surface is read, not invented

`scripts/analyze_epitope.py` parses **PDB 8W9J** — human CLEC12A bound to the
inhibitory antibody **50C1** at 3.5 Å — and takes every CLEC12A residue within
4.5 Å of the antibody as the epitope:

```
D198 S199 T200 R201 G202 M203 R204 V205 D206 R232 L233 Y234 Q236
```

Three arginines (R201, R204, R232) dominate that surface. The real antibody
answers them with acidic residues forming salt bridges and five tyrosines
stacking against the guanidinium groups (cation-π). That mechanism — read off a
molecule that demonstrably engages this epitope — becomes the design objective.

### 2. Candidates are built on a validated scaffold

The primary scaffold is the **affibody Z domain**: a 58-residue three-helix
bundle from staphylococcal protein A. Cysteine-free, folds without chaperones,
and affibody molecules are a genuinely validated scaffold class.

Only the **13 published affibody library positions** on helices 1 and 2 are
varied (9, 10, 11, 13, 14, 17, 18, 24, 25, 27, 28, 32, 35). The framework never
changes, so every candidate is a molecule that plausibly folds. Search space:
20¹³ ≈ 8.2 × 10¹⁶.

### 3. Scoring is decomposed and every part is explainable

Thirteen components in four groups. Each reports a raw physical quantity, a
normalised 0–1 score, a weight, an explanation, and a caveat.

| Group | Components |
|---|---|
| **recognition** | salt-bridge capacity, cation-π potential, charge complementarity, H-bond complementarity, paratope character, footprint match |
| **selectivity** | paralog selectivity (vs CLEC12B, CLEC9A, CLEC7A, CLEC1A, KLRB1), polyspecificity risk |
| **developability** | solubility, stability, manufacturability, size, scaffold consistency |
| **plausibility** | ESM-2 masked-marginal pseudo-perplexity |

Sources are named in the code — Kyte-Doolittle hydropathy, Eisenberg
consensus scale, Chou-Fasman propensities, Ikai aliphatic index, Tien et al.
maximum ASA, Pace et al. extinction coefficients.

The composite is a weighted mean of these, which is a *trade-off* model, and
one component escapes it: a design whose paralog-selectivity margin is
negative is **disqualified** rather than merely marked down, and sorts below
every candidate that is not. Its composite is still reported unchanged —
[disqualification is not a weight](#what-a-weighted-mean-cannot-say-disqualification-is-not-a-weight)
explains why that distinction had to exist and what it cost.

### 4. The score is anchored to a real molecule

A **50C1 scFv** is assembled from PDB 8W9J (VH + GGGGS×3 linker + VL) and scored
through the identical pipeline. An scFv is a real therapeutic format — it is
what sits on the front of a CAR — so the comparison is between two things
somebody could actually build.

Components that compare a candidate *to* the reference are excluded when the
reference scores itself, since that is circular.

### 5. The score is calibrated against chance

200+ random draws from the same 13-position library are scored to build a null
distribution. A candidate is reported as a **percentile against that null**,
which is the same logic as an enrichment score in a display campaign.

Critically, the module also measures **the scoring function's own discriminating
power**: how many standard deviations the real antibody sits above the random
decoys. If a scoring function cannot tell a real binder from noise, every
percentile it produces is worthless — so that figure is computed and displayed
rather than hidden.

Observed on the reference run:

| Control | Recognition | Percentile | Expected |
|---|---|---|---|
| 50C1 antibody (real binder) | 97.6 | **z = +3.5 σ** | should stand out ✓ |
| Wild-type Z domain | 50.8 | 54.5th | should be average ✓ |
| poly-Glu / poly-Ala / poly-Arg | 13–18 | 0th (z ≈ −2.9) | should be rejected ✓ |

The wild-type scaffold landing at exactly the 54th percentile is the important
one: it binds IgG Fc, not CLEC12A, and the scorer correctly finds it no better
than chance.

### Checking all the components at once

After the third saturation bug, finding them one at a time stopped being good
enough. `scripts/diagnose_components.py` scores a few hundred sequences and
reports the two failure modes that matter, for every component together:

- **Saturation** — a component stuck at one value supplies no gradient. It
  spends weight in the composite and contributes nothing to the search.
- **Redundancy** — two components measuring the same underlying quantity
  double-count it, so the composite silently weights that quantity twice.

Run it after changing any scorer. It found two real problems immediately.

**Four components were constant by construction.** On a fixed scaffold every
candidate is the same length with the same face size and the same framework,
so `size_suitability` and `footprint_match` had a standard deviation of
**exactly zero**, with `scaffold_consistency` and `polyspecificity_risk` close
behind. Together about **27% of the composite's weight** was a near-constant
bonus added to every candidate.

They are not wrong — they are real checks on the molecule — they just cannot
*rank*. So they are now measured, displayed, and excluded from the ranking
score, detected by measuring variance on the decoy panel rather than by a
hard-coded list (on a scaffold where they vary, they rank again).

| | with them | without them |
|---|---|---|
| ranking components | 13 | **9** |
| random library sd | 7.43 | **10.38** |
| random library range | 35.8 pts | **51.1 pts** |

Excluding them widened the usable range by about **43%** and lowered every
absolute score, which is the honest direction.

**And a redundancy I had introduced myself.** Fixing selectivity made
`paralog_selectivity` correlate with `salt_bridge_capacity` at **r = +0.870**
— both are driven by how well the face's acidic supply matches the epitope's
cation demand. At weight 1.4 each, the composite was weighting one quantity
2.8. `salt_bridge_capacity` dropped to 0.9; selectivity keeps the higher
weight because it measures the same matching *plus* the margin over the
paralog panel, which is strictly more information.

**What the diagnostic also showed:** `charge_complementarity`, which looked
like the weakest component at 71%, is one of the *best* — the highest variance
of any scorer (sd 0.320) against a sample mean of 0.330. It reads low because
it is demanding, which is what a useful scorer looks like. A low score is not
the same as a bad component, and the only way to tell them apart is to measure
the distribution.

### Three calibration bugs this project hit, and what they teach

All three were caught by controls, not by reading the code. All three would
have produced confident, meaningless numbers. They are the same bug wearing
different clothes: **a scoring component that cannot vary is not a
measurement, and it is worse than having no component at all** — it spends
weight in the composite while contributing no gradient to the search.

If you take one thing from this repository, take the habit that caught them:
score sequences whose answer you already know (a real binder, the unmodified
scaffold, poly-X controls), and check that the numbers come out in the order
reality says they should.

**1. Scorers that saturate tell you nothing.** The first version rated the
wild-type Z domain at 82% of a real anti-CLEC12A antibody. The Z domain binds
IgG Fc; it has no affinity for CLEC12A at all. Three recognition components
were awarding full marks for being vaguely acidic. The fix was to score charge
and hydrogen bonding as *densities on the binding face* matched to what the
reference achieves, rather than as raw counts with a generous tolerance.

**2. A metric anchored to the wrong family rates a real protein as fake.**
The ESM-2 plausibility scale was anchored to the 50C1 antibody scFv. Measured
perplexities on `esm2_t33_650M`:

| Sequence | Pseudo-perplexity |
|---|---|
| 50C1 scFv (antibody) | 2.91 |
| CLEC12A ectodomain | 5.12 |
| **Wild-type Z domain** (natural protein) | **10.08** |
| Random library members | 16.74 ± 1.38 |
| Shuffled Z domain | 20.76 |

ESM-2 has seen enormous numbers of antibodies and comparatively few
three-helix bundles, so anchoring an affibody to an antibody compares across
protein families. On that scale the natural parent protein scored **0%** — a
self-evident failure.

Worse, because *every* candidate saturated at zero, the component carried its
full weight while supplying **no gradient**. The search was blind on that axis
and drifted to a candidate at perplexity 17.67 — less natural-looking than a
random library member. A run with ESM switched off entirely found a better
sequence (13.89). **A saturated scorer is worse than no scorer.**

The fix anchors the scale to the scaffold's own family, measured at runtime:
`good` = the natural parent sequence, `poor` = the mean of random library
members on the same scaffold. The question becomes "closer to the natural
parent, or to random?" — which is answerable within one family.

`tests/test_pipeline.py --plm` now guards both: a natural protein must never be
rated implausible, and random members must not all saturate at one end.

**A caveat that remains.** Perplexity rewards predictability. A poly-arginine
face scores 10.45, close to the natural protein's 10.08, despite being a
useless binder. Only the complexity and recognition components catch that, so
the plausibility score is never meaningful alone.

**And a known limit of the fixed scale.** The top anchor is the natural parent
sequence, so anything at least that plausible scores 100%. In practice the
optimiser reaches perplexities below the parent's (8.75 in the reference run),
and 19 of 25 finalists sit at the ceiling. That is deliberate rather than a
regression: *more natural-looking than a natural protein* is not a meaningful
objective, and foldability genuinely saturates. The consequence is that among
finalists this component acts as a **filter** (reject the implausible) and not
as a **ranker**. Read it that way.

**3. A metric with a mathematical ceiling.** `paralog_selectivity` read 50% in
every single run. That constancy was the tell. A sweep of 680 binding faces
spanning the entire composition range found the **maximum reachable margin was
+0.075**, against a ramp that awarded full marks at +0.25 — so **50% was a
ceiling no design could exceed.** The optimiser had not been failing at
selectivity; it was already pinned at the top of a metric that could not
report higher.

The cause was again saturation: `min(supply, demand) / demand` caps at 1.0, so
a face with enough acidic residues for CLEC12A's three arginines also had
enough for CLEC9A's one. Both scored 1.0 and the difference collapsed to a
small steering term.

The fix penalises overshoot as well as shortfall, which is also the correct
chemistry: *an acidic residue with no cation to pair against on the target is
not neutral — it is precisely the residue that goes and binds something else.*

| | before | after |
|---|---|---|
| margin range | −0.55 to +0.075 | −0.80 to +0.35 |
| distinct values (680 faces) | 41 | 152 |
| real antibody 50C1 | 61% | **85%** |
| wild-type Z domain | 44% | **29%** |
| designed candidate | 50% (ceiling) | 90% |

Two further changes came with it. The off-target panel grew from 2 proteins to
7, now including **CLEC12B** — the immediate neighbour of CLEC12A in the NK
gene complex and by far the most likely cross-reaction. And **CD33 was removed
from the alignment**: it is a Siglec, a different fold family, so mapping the
CLEC12A epitope onto it produced a correspondence with no structural meaning.
It stays in the panel as documented context. Every accession is now verified
against its expected gene symbol at download time, because a typo would
silently substitute a different protein and the comparison would still look
fine.

The scoring scale is no longer hard-coded. It is measured per target from the
reachable margin range: the floor is the 10th percentile of random faces, the
ceiling the best a systematic composition sweep can reach. That keeps the
component discriminating and stays correct when paralogs are added.

**A caveat that remains.** A designed face can now score *above* the reference
antibody here (90% vs 75%) without being more selective. The metric rewards
matching a surface's composition; an antibody's real selectivity comes from
shape complementarity this cannot see.

### Cross-reactivity panels, and what an absent one means

Selectivity — would this binder also stick to something the patient needs —
is the axis that has actually killed candidates. Scoring it needs a panel of
proteins the binder must *not* engage.

These are derived, not written. `scripts/build_panels.py` asks UniProt which
domain families a target belongs to, pulls every reviewed human protein in
those families, and ranks them by identity **at the epitope window** rather
than across the whole protein. Family membership says two proteins are
relatives; it does not say they resemble each other at the patch a binder
would touch, which is the only thing the score is about.

The results are checkable against biology. CD33's panel is the Siglecs,
EGFR's is ERBB2/3/4, VEGF-A's is VEGF-B/C/D and PlGF, IL-17A's is IL-17F,
BAFF's is APRIL, α-synuclein's is β- and γ-synuclein. For CLEC12A the
derived panel independently returns CLEC9A as the nearest relative, which is
what the hand-written panel said.

**An absent panel records which of six things happened**, because they are
very different claims:

| State | Means | n |
|---|---|---|
| derived / curated | a panel exists and selectivity is scored | 15 |
| `no_human_relatives` | no reviewed human protein shares a domain family — a *finding*, not a gap | 13 |
| `no_close_relatives` | relatives exist but none above 25% identity at this surface | 6 |
| `panel_too_weak` | one distant relative, which the controls show is unreliable | 3 |
| `family_too_generic` | the only family on record is shared by hundreds of proteins | 3 |
| `unmappable` / `no_family` | the site cannot be located, or no family is annotated | 2 |

"This protein has no human relatives" and "nobody looked" are opposite
statements about how much is known. Collapsing them into one grey *not
assessed* would be the single most misleading thing this interface could do.

#### The controls, and the inversion they caught

`scripts/diagnose_panels.py` runs four controls per target whose answers are
known in advance: the target's own sequence must map back to exactly its own
epitope; the real binder's contact residues should score well; random faces
should not; and an **anti-target face** — one optimised to complement the
closest *relative* — should score clearly worse than the real binder. Across
the 37 scored targets the reference averages 28, random 19 and the anti-target
8, and the reference beats the anti-target in 33 of them. The control also
reports the raw margin's **sign** now, which it never did — nine of fifteen
references were negative while the table above them read 40 to 60.

Running them caught two real defects and one inversion.

A UniProt entry for the full IgE heavy chain, P0DOX4, entered IgE's panel as
a "relative". Its mapped surface was character-for-character the target's own
epitope. That is not noise but fatal: selectivity is the margin over the
*closest* relative, so one member scoring identically to the target pins that
margin at zero for every candidate forever. Panel members are now rejected on
surface identity, which is the thing that matters, rather than on
whole-sequence identity, which put these two at 64% and let it through.

The inversion was worse. For α-synuclein against β-synuclein, and IL-17A
against IL-17F, the best margin reachable by **any** composition is exactly
0.000 — the surfaces are 9/11 and 8/11 identical. The ramp mapped that onto
**100%**. The targets where selectivity is *impossible* were scoring as
perfectly selective. Those cases now report that no face can separate the
pair, and the component is excluded from ranking rather than scored. This is
a real result, not a shortcoming: anti-IL-17A antibodies genuinely cross-react
with IL-17F, which is why bimekizumab was designed to hit both.

A score that cannot vary has now appeared in this project five times. It has
been a floor, a ceiling, a constant, a saturated maximum — and this time a
maximum that meant its own opposite. The section below adds a sixth: a score
that varied when nothing had varied at all.

#### Following the two flagged targets, and what they were really about

The controls flagged TNF-α (the reference scoring below the anti-target face)
and CD38. Both flags were true, and neither was about those targets. They
were two defects in the metric itself, plus a bug in the control that found
them.

**The kernel was reading surface length as if it were chemistry.** Write the
same epitope out twice — same residues, same proportions, not one atom
changed — and the kernel moved 0.475 → 0.725 → 0.508 across one, two and
three copies; doubling the face instead moved it to 0.275. It compared the
face's raw residue counts against the epitope's, so the answer depended on
the ratio of two lengths. Epitopes in this catalogue run from 11 to 34
residues, and **28 of the 70 panel surrogates are a different length from the
surface they are compared against**, so two fifths of every cross-reactivity
comparison was partly a measurement of how many residues a curator had
annotated. The epitope's demand is now scaled to the size of the face, which
is what a 13-residue binder can actually reach. Every other scorer in
`interface.py` had been working in densities all along.

**The control was testing the wrong molecule.** It cut each reference
paratope to the scaffold's 13 positions with `[:13]` — a contiguous prefix.
For CD38 that discarded six of isatuximab's nine tyrosines. Since the kernel
reads composition, slicing the molecule changed which molecule was being
tested: adalimumab's face scored 11, 46 or 31 depending only on where it was
cut. With a length-invariant kernel the reference needs no cutting at all.

**The ramp was throwing away the sign.** This is the one that mattered. The
margin's zero is a real boundary — above it a face complements the target
more than its nearest relative, below it the reverse — and the score ran from
the 10th percentile of *random* faces up to the ceiling. When random faces
are catastrophically anti-selective, merely being less bad than random reads
as a majority-good score:

| target | optimised design's margin | was reported as |
|---|---|---|
| TNF-α | **−0.153** | 65.6% selective |
| EGFR | **−0.011** | 70.7% selective |
| HER2 | **−0.030** | 54.2% selective |

A design that prefers the off-target is not weakly selective; it is pointed
at the wrong protein of the two. The whole region below zero is now
compressed into the bottom quarter of the score and labelled
**counter-selective** on the tile. The band is kept rather than flattened to
zero because the optimiser needs a gradient to climb out of it, and a score
pinned at one value is the failure this codebase has already hit five times.

#### A fourth inversion, and a threshold in the wrong units

`INDISTINGUISHABLE_MARGIN` — an absolute 0.02 — was built to catch exactly
the α-synuclein case. VEGF-A's ceiling is **0.029**, so it passed, and the
VEGF-A design was then reported as **100% selective**: a perfect score for
3% of the margin CLEC12A can reach.

The cutoff was the right idea in the wrong units. Drop one residue from one
mapped surrogate — a single alignment error, which at 25–40% identity is a
certainty rather than a hypothetical — and recompute the ceiling. It moves by
0.003 for α-synuclein and by 0.077 for TROP2, a factor of twenty-five. No
fixed number can span that. `scripts/measure_ceilings.py` now measures each
ceiling and its own error bar, per face length, and stores both in
`panels.json`; a ceiling must stand two standard deviations clear of zero to
carry a percentage. VEGF-A, at 1.1, no longer does. CLEC12A, at 17.8, does
comfortably.

The stored ceiling is used rather than the live one, because
`achievable_margin_range` substitutes a usable value when the sweep finds no
positive margin at all — on a 10-residue face no composition beats CD33's
Siglecs, the true ceiling is **−0.012**, and the substituted value would have
waved it straight through.

#### The number that is now on the page: does this measure agree with a drug?

After all of that, adalimumab's margin against TNF-α is still **−0.268**, and
isatuximab's against CD38 is **−0.087**. Across the 37 targets that carry a
score — the catalogue has since grown to 95 — the reference molecule's own
contact residues prefer the target to its nearest relative in **11**.
Twenty-six of thirty-seven do not, and that proportion has now held while the
catalogue more than doubled, which makes it a property of the measure rather
than a quirk of the first dozen targets.

The fix for that is not to tune the kernel until the antibodies pass — that
is fitting four free parameters to eleven points. It is to put the fact where
the person reading the score can see it. Every target page with a panel now
carries the reference check: *isatuximab's own contact residues score 2
against this panel, margin −0.087. That is negative: this measure says
isatuximab complements BST1 better than it complements CD38.* A measure that
is wrong about the one molecule whose selectivity is established on that
target has not earned confidence in its ranking of untested ones.

Why these two in particular: CD38's only human relative, BST1, presents a
surface of almost exactly the same composition — 5 basic/5 acidic/1
aromatic/10 apolar against 6/5/2/10 — at 29% sequence identity. A
composition-only measure has nothing to work with there, however different
the two sequences are. For TNF-α the epitope is strongly acidic, and the
kernel therefore wants a basic, aromatic face; adalimumab's selectivity comes
from CDR shape complementarity that no composition term can see.

The score is also **scaled per target** — 100 means "as selective as this
epitope and this panel allow" — so it is not comparable between targets, and
a high score against distant relatives is an easier result than a moderate
one against close relatives. Both the interface and the component caveat say
so.

### What a weighted mean cannot say: disqualification is not a weight

Having made counter-selectivity visible, the obvious next question was
whether the composite should let a counter-selective design score 82.6 and
come first. It should not, and the problem was worse than presentation.

The composite is a weighted mean over nine ranking components. A mean is a
*compensatory* aggregation: it assumes every axis trades against every other,
so nine good scores outvote one bad one. That is right for most of these —
aromatic content really does trade against solubility — and wrong for an axis
whose zero means *the molecule does the opposite of its job*.

Measured over the candidates each search actually produced:

| target | ranked first | its margin | positive-margin candidates found | best such |
|---|---|---|---|---|
| TNF-α | 82.6 | **−0.218** | 25 | 72.6 |
| CD22 | 80.7 | **−0.325** | 167 | 75.8 |
| HER2 | 78.5 | **−0.030** | 71 | 75.7 |
| IL-1β | 79.9 | **−0.008** | 133 | 77.3 |

The search had found candidates that complement the target more than its
nearest relative, and the composite ranked one aimed at the *wrong protein*
above all of them. This was not a display problem: the composite is what the
loop optimises, so it had been spending its budget in the wrong region and
then recommending the result.

**Raising selectivity's weight would not fix it, and that is the useful
part.** The axis can swing about 3.7 composite points between a
counter-selective design and a decent one; the other nine can differ by more.
Any weight large enough to veto would also distort the ranking among the good
candidates, where selectivity should be one consideration among several. The
problem is not that disqualification is underweighted. It is that
disqualification is not a weight.

So a component may now be **disqualifying**, and ordering is lexicographic —
`(not disqualified, composite)` — while the composite itself is left exactly
as it was:

- **The number does not change.** A counter-selective TNF-α design really
  does score 82.6 on the nine axes that are not selectivity. Suppressing that
  would be its own dishonesty, and it is what tells you the candidate is
  one axis away from being useful.
- **The order changes.** Nothing disqualified outranks anything that is not.
  The dashboard marks those rows and the breakdown says which component
  disqualified the candidate and why.
- **The search changes too.** Every comparison in `loop.py` goes through the
  rank key — tournament selection, elitism, the stagnation counter, the
  early stop. Fixing only the final sort would have left the search
  optimising toward counter-selective designs and reporting the least bad
  thing it passed on the way.
- **The counter-selective band earns its keep.** Scoring the whole negative
  region as a flat zero would strand the search with no gradient to climb.
  The 0–25% band is the slope out.

Across all 42 targets every returned candidate now has a positive margin, and
no disqualified candidate appears anywhere in a top-12. The cost is real and
small: on the three targets where the winner changed, the composite fell by
5.8, 5.3 and 1.0 points — and selectivity rose from 19 to 50, 11 to 67, and
21 to 92.

Two details that only show up once it is running. The stagnation counter had
to move to the rank key as well, because the generation that first finds a
clean candidate usually has a *lower* composite than the counter-selective
one it displaces — compared as bare floats that reads as the search
stalling, and stops it exactly when it starts working. And the
score-over-generations chart can step down at that generation for the same
reason, so the page explains the step rather than leaving it looking like a
regression.

If nothing qualifies, ordering degenerates to the plain composite — the old
behaviour — and the run says so in a banner, because a table of candidates
none of which is selective looks identical to a good one.

### A fourth bug, of a different kind: telling one binder from two

The three calibration bugs above were all *a number that could not vary*.
This one is worth recording because it is a different failure and it took
four attempts to get right.

When a structure holds a target bound by more than one thing — two Fabs at
two epitopes, a Fab and the protein's own partner subunit, an antibody bound
to the next copy over in the crystal — the pipeline has to decide which
chains constitute *one binder*. Merge too eagerly and you brief the designer
on a site spanning two faces that no single molecule can occupy. Split too
eagerly and you hand it half an epitope.

| Rule | Why it looked right | What it broke |
|---|---|---|
| Antibody chains share a J-segment motif | Fabs do have those motifs | A light chain from a different germline has none, so the Fab was split and scored as half a binder — five entries |
| Chains that pack against each other are one molecule | True of a heavy/light pair | Neighbouring Fabs in a crystal also pack; two bezlotoxumab Fabs merged into a 40-residue "epitope" across two sites |
| Cluster the contacted residues into patches | Asks the question directly | Genuinely elongated interfaces split; SARS-CoV-2 RBD/ACE2 lost eight of its twenty residues |
| Patches adjacent **and** chains packed **and** not two copies of one sequence | — | holds on all 42 |

Each rule was correct about the case that motivated it and wrong about a case
it had never seen. What made the difference was not a cleverer rule but
checking the output of every entry in the catalogue after each change —
epitope sizes and binder chains, all forty-two, every time. Two of the four
rules looked perfect on the target the work started from.

The invariant that finally held is stated as three conditions because each
one was paid for: `src/plab/epitope.py:group_binders`.

### 6. The loop

An evolutionary search over the designed positions only — tournament selection,
uniform crossover, helix-aware mutation (never proline or glycine at a helical
position), elitism, and memoisation. It stops on a score plateau or a target
score.

A plateau means *the search found nothing better*, not *the molecule works*.
The report keeps those separate.

---

### Making it fast, by measuring rather than guessing

The suite reached eighteen minutes as the catalogue grew. I had a theory
about why, wrote it down, and it was wrong — twice. A profiler settled it in
one run:

```
86 of 86 seconds in achievable_margin_range
2,066,540 calls to complementarity_kernel   (for four targets)
```

The composition sweep is O(face_length⁴). For the 13-residue design face
that is 2,380 compositions and nothing to worry about. But the sweep also
runs for the **reference molecule**, whose face is thirty or forty residues
— and 41⁴/24 is a hundred and thirty-five thousand. One reference took
twenty seconds; ninety-five of them took the afternoon.

It does not need to scale with length. Since the kernel was fixed to read
densities rather than raw counts, its value depends on the face only through
the four class *fractions*, so a coarser grid on a long face covers the same
chemistry. The step is 1 for any face of 13 residues or fewer, which is
every scaffold here, so **all 106 stored ceilings reproduce to the digit**.

The second cost was real waste: `build_panels.py` maps the epitope onto each
panel member and stores the answer, and then everything downstream re-derived
the identical alignments once per Evaluator. They are read from the panel
file now — but only after `scripts/verify_panel_cache.py` has re-derived every
surface and confirmed it matches, which it records by storing the epitope the
surfaces belong to. All 53 panels verified, 0 stale, in 39 seconds. A
mismatched or unmarked panel falls back to deriving its surfaces.

**18 minutes → 1 minute 36 seconds, with five more tests than before.** A
design run went from several seconds to under two.

The interesting part is the bug this caused. My first version of the cache
guard compared the filtered result against its own intersection — true by
construction, so it checked nothing. CLEC12A's panel quietly fell from nine
members to eight, because one curated member has no stored surrogate, and its
reachable ceiling moved from 0.350 to 0.359 without an error or a warning. I
caught it only because I re-derived all 106 stored ceilings and compared.
`test_panel_cache` now asserts the cached path and the derived path return
the identical dictionary, and that a moved binding site invalidates the
cache. A cache whose correctness is not tested is a way of being confidently
wrong faster.

## Layout

```
data/structures/     PDB and mmCIF files, one per catalogue entry
data/catalog/        verified.json, panels.json, derived epitope JSON
data/sequences/      UniProt FASTAs for targets and panel members
data/families/       cached Pfam family FASTAs
data/lineage/        expression.json: the transcript screen, with provenance
                     topology.json: UniProt location evidence, no verdicts
scripts/             fetch_data, analyze_epitope, resolve_gene, rcsb_lookup,
                     build_catalog, build_panels, measure_ceilings,
                     verify_panel_cache, fetch_lineage, fetch_topology,
                     diagnose_*,
                     smoke_test
src/plab/
  aminoacids.py      property tables, each with its published source
  pdbio.py           PDB and mmCIF readers + interface analysis
  align.py           Smith-Waterman with a property-derived matrix
  lookup.py          gene symbol -> accession -> real PDB entries
  epitope.py         chain matching, binder grouping, binding-site derivation
  targets.py         the design brief, loaded from data/
  panels.py          cross-reactivity panels from domain families
  custom.py          targets built from any protein, for the web app
  lineage.py         which healthy cells carry a target, and how well that
                     is known -- evidence grading, never a bare number.
                     Two curated tables: blood lineages, and tissue by site.
                     Also the topology axis: whether the target is on a cell
                     surface at all, and so whether the question applies;
                     and what is lost when a signal is blocked, for the
                     targets where that is the only safety question left
  specimen.py        ranking the catalogue against one antigen profile
  design/            scaffolds, the evolutionary loop
  scoring/           base, biophysics, interface, selectivity, plm, calibration, composite
  report/            dashboard generator
out/                 run JSON + dashboard.html
run.py               CLI entry point
```

## Extending to another target

1. Find a structure, don't recall one:
   `python scripts/rcsb_lookup.py pick <UNIPROT_ACCESSION>` lists the entries
   that really contain that protein with something bound to it, antibody
   complexes first.
2. Add a `CatalogEntry` to `src/plab/catalog.py` with the accession and the
   PDB id from step 1. Where a structure holds two different binders, set
   `binder_preference` to say which one the entry is about.
3. `python scripts/build_catalog.py --only <key>` to verify and add it.
4. `python scripts/diagnose_catalog.py <key>` if it is rejected — it prints
   every chain's identity to the reference by both measures, which is usually
   enough to see whether the problem is the data or the entry.
5. `python scripts/build_panels.py --only <key>` to derive its
   cross-reactivity panel, then `python scripts/measure_ceilings.py <key>` to
   measure how big a margin that panel makes reachable and how big the error
   bar on it is. Without the second step the target falls back to an absolute
   cutoff, which is the one that let VEGF-A report 100%.
6. `python scripts/diagnose_panels.py <key>` to check the panel discriminates
   before trusting a selectivity number from it. Read the `ref_m` column, not
   just the percentage: a negative reference margin means this measure
   disagrees with a molecule known to work on that target.
7. `python scripts/verify_panel_cache.py <key>` so the mapped surfaces can be
   read from the panel file instead of re-aligned on every run. It re-derives
   each one and compares before marking it, so a panel that no longer matches
   its target is left alone rather than trusted.

The scorers are target-agnostic — they read the epitope's chemistry from the
brief rather than hard-coding CLEC12A's — so nothing below the catalogue
needs touching.

## Data sources

Structures from the [RCSB PDB](https://www.rcsb.org/); sequences from
[UniProt](https://www.uniprot.org/), which also supplies the subcellular
location and membrane-feature annotation behind the topology axis
(`scripts/fetch_topology.py` -> `data/lineage/topology.json`). Transcript
abundance per blood cell type and tissue from the [Human Protein
Atlas](https://www.proteinatlas.org/) (`scripts/fetch_lineage.py` ->
`data/lineage/expression.json`). ESM-2 weights from
[facebookresearch/esm](https://github.com/facebookresearch/esm).

All downloaded data is public. Nothing in `data/` is generated or guessed; if a
file is missing the pipeline fails loudly rather than substituting invented
values.
