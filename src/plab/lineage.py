"""Which healthy cells carry the target, and how well that is actually known.

The selectivity scorer in scoring/selectivity.py answers one question: could a
binder confuse this target with a *structurally related* protein? It does that
by mapping the epitope onto paralogs and comparing composition. Its own caveat
says what it cannot do:

    "It also does NOT predict expression on healthy cells: whether an
    off-target protein sits on a cell that matters is measured by flow
    cytometry on primary human samples, and it is the reason CLEC12A was
    chosen over CD33 in the first place -- the two are not close relatives
    at all."

That is this module. CD33 and CLEC12A share no domain family, so no paralog
panel will ever put them in the same room; the entire clinical argument for
preferring CLEC12A is about which healthy cells carry it. The two axes are
orthogonal and neither substitutes for the other.

--- The rule that shapes everything here ---------------------------------

Transcript abundance may RAISE a concern. It may never CLEAR one.

This is not caution for its own sake; it is the result of checking. GPRC5D is
a myeloma target, and in the clinic 75% of patients on talquetamab get nail
changes and 63% get oral toxicity, from on-target expression in the nail
matrix, hair cortex and the filiform papillae of the tongue. Ask the Human
Protein Atlas for GPRC5D and the answer is split: 26.4 nTPM in skin, which is
a real flag, against 0.0 in tongue and 0.2 in bone marrow. The screen sees
part of the toxicity and is blind to the rest of it.

That split is the useful finding, more than a clean miss would have been. A
column reading zero cannot be taken as clearance even for a protein the same
screen flags elsewhere, because the tissue the atlas samples and the structure
the toxicity lands in are not the same thing -- there is no nail matrix in a
bulk tongue sample.

So a lineage with no detected transcript and no cytometry behind it is
reported UNKNOWN, never ABSENT. Only flow cytometry on primary human cells, or
observed human toxicity, can establish absence. A tool that let nTPM < 1 mean
"safe" would hand back a confident all-clear on the oral toxicity in two
thirds of patients, and the whole point of this package is to not do that.

(The first draft of this file claimed the screen missed GPRC5D entirely. That
was wrong, and wrong for an instructive reason: the probe behind it asked for
`t_RNA_skin`, which is not a column code, and the atlas answers an unknown
column by silently leaving it out. The guard in scripts/fetch_lineage.py
exists because of this, and it is what produced the corrected number.)

--- Blood is not the only axis, and for most targets it is not the one ------

This module began as a blood story because the first antigens in it were
blood cancer antigens, where the lineage table carries the whole argument.
Extending it to the solid tumour and immuno-oncology entries moved the weight
to the other table. HER2's established liability is cardiomyocytes,
nectin-4's is skin, mesothelin's is the lining of the pleura, PSMA's is the
salivary acini; none of those is a blood population, and a flow panel on a
marrow specimen cannot address any of them. So a target can come back with no
haematopoietic conflicts at all and be the most dangerous thing on the list.

That is the quiet failure, and it fails in the reassuring direction, so
specimen.assess_target now says it on the card rather than leaving an empty
conflict list to be read as a clean result.

The same extension turned up two more ways this screen misleads, to set
beside the GPRC5D false negative below. CLDN18 has two splice variants:
CLDN18.1 in lung, CLDN18.2 in gastric mucosa, and zolbetuximab binds only
the second. The screen reads the GENE, so it was reporting a loud 237.8 nTPM
flag in the lung isoform the drug is specific against while saying nothing
about the stomach -- there was no stomach column until this work added one,
and it reads 684.2. A screen can be precise, well-formed and pointed at the
wrong molecule. And VEGFA is detected in all seventeen screened tissues,
which is true, on-target and nearly useless: a secreted ligand makes "which
healthy cells carry it" the wrong question in the first place.

--- What no expression table can see --------------------------------------

Worth stating because the most honest entry in this file is an absence.
Mirvetuximab soravtansine's dose-limiting toxicity is corneal, and FRalpha is
absent from human corneal tissue -- the ocular injury appears across ADCs
irrespective of where their targets are expressed. Sacituzumab govitecan's
neutropenia and diarrhoea are SN-38 effects, the same ones irinotecan
produces with no antibody at all. Telisotuzumab vedotin's neuropathy in 45%
of patients is MMAE.

In each case the payload carries toxicity of its own, and this axis would
return a clean answer and be right about the thing it measures. The healthy-
cell question is one input to deciding what to build. It was never going to
be the one that decides whether a molecule is tolerable.

--- Why stem cells and progenitors are separate --------------------------

They were one bucket in the first draft of this file, which would have been
wrong. CD33 is largely absent from the most primitive CD34+CD38- stem cells
but present on committed multilineage progenitors, with intensity falling as
cells mature, and it is the progenitor expression that drives gemtuzumab's
myelosuppression. Collapsing the two would make the contrast with CLEC12A --
which is absent from both -- disappear into a single misleading cell.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
EXPRESSION_PATH = ROOT / "data" / "lineage" / "expression.json"

# --- Evidence ---------------------------------------------------------------

# What kind of observation stands behind a level, strongest first. These are
# not interchangeable and must never render alike: "nobody looked" and
# "cytometry says it is not there" are opposite claims about how much is
# known, and the same mistake in the paralog panels (panel_state) is the
# reason that field exists.
# "animal" was added for the solid tumour antigens and sits deliberately
# below human cytometry. TROP2 is why: sacituzumab govitecan's dose-limiting
# neutropenia and diarrhoea are SN-38 effects, and the cynomolgus tox studies
# -- in a species whose TROP2 distribution resembles ours -- found no
# histopathological change in the normal tissues that carry it. That is real
# in vivo evidence about a population surviving, and it is not a human
# measurement, so it must not be able to CLEAR a human population or to
# exclude a target. Both of those thresholds now read this ranking rather
# than hard-coding a pair of strings, so the ordering is stated once.
# "ihc" was added at the same time and for a plainer reason: the solid tumour
# rows are about tissue, and tissue is read by immunohistochemistry, not by
# cytometry. Filing nail-matrix or gastric-mucosa staining under "flow
# cytometry on primary human cells" would attach a provenance the finding
# does not have, in a field whose whole job is to say how a thing is known.
# It sits alongside cytometry rather than below it -- both measure protein on
# primary human material -- so both can settle a question.
EVIDENCE_RANK = {"clinical": 5, "flow": 4, "ihc": 3, "animal": 2,
                 "rna": 1, "none": 0}

#: Evidence at or above this rank can establish an absence, or exclude a
#: target outright. Below it, a finding can only raise or soften a flag.
EVIDENCE_DECISIVE = EVIDENCE_RANK["ihc"]

EVIDENCE_LABEL = {
    "clinical": "observed human toxicity",
    "flow": "flow cytometry on primary human cells",
    "ihc": "immunohistochemistry on human tissue",
    "animal": "in vivo animal study",
    "rna": "transcript screen only",
    "none": "not assessed",
}


def decisive(evidence: str) -> bool:
    """Whether this kind of observation can settle a question, not just flag."""
    return EVIDENCE_RANK.get(evidence, 0) >= EVIDENCE_DECISIVE


def phrase(label: str) -> str:
    """A lineage label as it reads mid-sentence.

    `label.lower()` was doing this and turning "T cells" into "t cells" and
    "NK cells" into "nk cells". Lowercase the first letter only where the
    second character shows the word is ordinary prose rather than an
    initialism: "Monocytes and macrophages" lowercases, "T cells" and "NK
    cells" do not.
    """
    if len(label) > 1 and (label[1].isupper() or label[1] == " "):
        return label
    return label[:1].lower() + label[1:]

# Levels. "contested" is a real answer, not a hedge -- for some antigens the
# published cytometry genuinely disagrees, and flattening that to a number
# would invent a consensus that does not exist.
LEVELS = ("absent", "low", "moderate", "high", "contested", "unknown")

# Only these can be reached from transcript data. Note the omission of
# "absent": see the module docstring.
RNA_BANDS = ((50.0, "high"), (10.0, "moderate"), (1.0, "low"))


def band_for_ntpm(ntpm: float) -> str:
    """Transcript abundance -> level, with no route to "absent".

    The thresholds follow the atlas's own convention that nTPM >= 1 counts as
    detected; the bands above that are coarse on purpose, because the map from
    transcript to surface protein density is not linear and for some proteins
    not even monotonic. Three bands can be defended. Ten could not.
    """
    for cut, level in RNA_BANDS:
        if ntpm >= cut:
            return level
    return "unknown"


# --- Lineages ---------------------------------------------------------------


@dataclass(frozen=True)
class Lineage:
    """One population a clinician may need to keep alive.

    `blood_columns` and `tissue_columns` are Human Protein Atlas column codes.
    An empty pair means the atlas has nothing for this population, which is a
    fact about the data rather than about the biology -- the atlas's bone
    marrow single-cell set is 3,484 cells annotated as T-cells, B-cells,
    macrophages, erythroid and plasma cells, with no stem or progenitor
    population in it at all.
    """

    key: str
    label: str
    consequence: str
    blood_columns: tuple[str, ...] = ()
    tissue_columns: tuple[str, ...] = ()
    group: str = "haematopoietic"


LINEAGES: tuple[Lineage, ...] = (
    Lineage(
        "hsc", "Haematopoietic stem cells",
        "Lasting marrow failure. The population that has to survive for blood "
        "counts to recover at all, and the reason a target can be lethal "
        "despite perfect tumour selectivity.",
        # Nothing. The atlas cannot see this population, and it is the single
        # most consequential row in the table.
    ),
    Lineage(
        "progenitor", "Committed progenitors",
        "Prolonged cytopenias while the marrow reconstitutes from upstream "
        "stem cells. Recoverable, unlike loss of the stem compartment.",
    ),
    Lineage(
        "neutrophil", "Neutrophils and granulocytes",
        "Bacterial and fungal infection, in patients who are frequently "
        "already neutropenic from induction chemotherapy.",
        blood_columns=("blood_RNA_neutrophil", "blood_RNA_eosinophil",
                       "blood_RNA_basophil"),
    ),
    Lineage(
        "monocyte", "Monocytes and macrophages",
        "Impaired clearance of pathogens and of debris; monocytes are also "
        "the cells that execute antibody-dependent phagocytosis.",
        blood_columns=("blood_RNA_classical_monocyte",
                       "blood_RNA_intermediate_monocyte",
                       "blood_RNA_non-classical_monocyte"),
    ),
    Lineage(
        "dendritic", "Dendritic cells",
        "Loss of antigen presentation, which blunts any adaptive response "
        "the therapy is trying to recruit.",
        blood_columns=("blood_RNA_myeloid_DC", "blood_RNA_plasmacytoid_DC"),
    ),
    Lineage(
        "b_cell", "B cells",
        "Loss of new antibody responses and hypogammaglobulinaemia, usually "
        "manageable with immunoglobulin replacement.",
        blood_columns=("blood_RNA_naive_B-cell", "blood_RNA_memory_B-cell"),
    ),
    Lineage(
        "plasma_cell", "Plasma cells",
        "Loss of existing antibody titres, including vaccine protection "
        "built up over a lifetime.",
        # Bulk bone marrow was the first draft's proxy here and it is not
        # one. Marrow is overwhelmingly not plasma cells, so a myeloid
        # antigen's marrow signal would be read as plasma-cell expression --
        # CD33 reads 30.6 nTPM in marrow and belongs to none of it. Better
        # to have no screen for this population and say so than to have a
        # number that means something else.
    ),
    Lineage(
        "t_cell", "T cells",
        "Immunosuppression, and for any T-cell-engaging format it is "
        "self-defeating: the therapy destroys the cells it recruits.",
        blood_columns=("blood_RNA_naive_CD4_T-cell",
                       "blood_RNA_naive_CD8_T-cell",
                       "blood_RNA_memory_CD4_T-cell",
                       "blood_RNA_memory_CD8_T-cell",
                       "blood_RNA_T-reg", "blood_RNA_gdT-cell",
                       "blood_RNA_MAIT_T-cell"),
    ),
    Lineage(
        "nk_cell", "NK cells",
        "Removes the effector arm that antibody-dependent killing depends "
        "on, so the binder disables its own mechanism.",
        blood_columns=("blood_RNA_NK-cell",),
    ),
    Lineage(
        "erythroid", "Erythroid cells",
        "Transfusion-dependent anaemia. Anti-CD47 agents demonstrate this as "
        "an on-target effect rather than a theoretical one.",
    ),
    Lineage(
        "platelet", "Platelets and megakaryocytes",
        "Bleeding risk, compounded in patients already thrombocytopenic.",
    ),
)

LINEAGE_BY_KEY = {l.key: l for l in LINEAGES}

# Sparing presets a clinician actually reasons in. The first is deliberately
# the default: for any patient who needs their own marrow back, the stem and
# progenitor compartments are the ones that cannot be replaced.
PRESETS: tuple[tuple[str, str, tuple[str, ...]], ...] = (
    ("marrow_recovery", "Marrow must recover",
     ("hsc", "progenitor", "erythroid", "platelet")),
    ("effector_intact", "Keep the effector cells",
     ("t_cell", "nk_cell", "monocyte")),
    ("infection_risk", "Patient is neutropenic / infected",
     ("neutrophil", "monocyte", "dendritic")),
    ("humoral_intact", "Preserve antibody immunity",
     ("b_cell", "plasma_cell")),
)

# --- Non-haematopoietic screen ---------------------------------------------

# Always run, never opt-in, because the clinician does not know to ask. A
# blood-only sparing model scores GPRC5D as beautifully restricted -- plasma
# cells and nothing else in blood -- and says nothing about the nails of
# three patients in four. Being wrong in that direction is the failure this
# screen exists to make less likely, and on GPRC5D it half succeeds: skin
# comes back at 26.4 nTPM, which is the flag that would make somebody look,
# while tongue reads 0.0. Raising half a concern is worth a great deal more
# than raising none, and is still not clearance.
TISSUE_SCREEN: tuple[tuple[str, str], ...] = (
    ("t_RNA_skin_1", "skin"),
    ("t_RNA_tongue", "tongue"),
    ("t_RNA_salivary_gland", "salivary gland"),
    ("t_RNA_esophagus", "oesophagus"),
    ("t_RNA_stomach_1", "stomach"),
    ("t_RNA_duodenum", "duodenum"),
    ("t_RNA_pancreas", "pancreas"),
    ("t_RNA_heart_muscle", "heart muscle"),
    ("t_RNA_lung", "lung"),
    ("t_RNA_liver", "liver"),
    ("t_RNA_kidney", "kidney"),
    ("t_RNA_urinary_bladder", "urinary bladder"),
    ("t_RNA_colon", "colon"),
    ("t_RNA_breast", "breast"),
    ("t_RNA_cerebral_cortex", "cerebral cortex"),
    ("t_RNA_thyroid_gland", "thyroid gland"),
    ("t_RNA_testis", "testis"),
)
# The five added after the solid tumour antigens were curated, each for a
# named reason rather than for coverage: stomach because CLDN18.2 lives in
# gastric mucosa and the screen was reporting only the lung isoform it is not
# aimed at (see CURATED_TISSUE["CLDN18", "stomach"]); duodenum and pancreas
# because PSMA, CEACAM5 and mesothelin all reach the upper GI tract, which
# colon and oesophagus alone do not cover; urinary bladder for NECTIN4; and
# breast for HER2 and KIT. Bone marrow is still deliberately absent -- see
# the plasma_cell lineage above for why a bulk marrow column is worse than
# no column.

ALL_COLUMNS: tuple[str, ...] = tuple(dict.fromkeys(
    [c for l in LINEAGES for c in l.blood_columns]
    + [c for l in LINEAGES for c in l.tissue_columns]
    + [c for c, _ in TISSUE_SCREEN]
))


# --- Curated cytometry ------------------------------------------------------


# Whether engaging this antigen is actually known to damage this population.
# Separate from `level` because CD47 proves they are different questions: it
# is carried on haematopoietic stem cells, and blocking it still does not
# increase phagocytosis of normal stem and progenitor cells, because the
# therapeutic window is a density difference rather than an absence. A model
# that read "expressed here" as "destroys this" would rule out one of the few
# targets where human experience says otherwise.
#
#   observed      human experience shows damage to this population
#   not_observed  human or in vivo experience shows the population survives
#   unknown       nobody has established either way -- NOT a free pass
#
# The checkpoint antigens then showed that `harm` is not a property of the
# antigen and the cell alone. PD-L1 sits on monocytes and dendritic cells.
# Avelumab, a wild-type IgG1, kills PD-L1-positive dendritic cells by ADCC;
# atezolizumab and durvalumab, which carry an engineered Fc, do not lyse the
# same cells in the same assay. Same target, same population, opposite
# answers, decided by the format rather than by the biology.
#
# This table records the verdict for the format that was actually tested, and
# every entry whose answer depends on that says so in its note. An antibody
# chosen to be silent, an ADC carrying a microtubule poison and a T-cell
# engager are three different questions about one row, and the one thing this
# field must never do is let a blocking antibody's safety record be read as
# clearance for a cytotoxic payload aimed at the same antigen.
HARM = ("observed", "not_observed", "unknown")


@dataclass(frozen=True)
class Curated:
    level: str
    evidence: str          # "flow" or "clinical"
    source: str
    note: str = ""
    harm: str = "unknown"


# Facts established by cytometry on primary human cells, or by observed human
# toxicity. Every entry carries its source, because three of this catalogue's
# first rejections were PDB ids recalled from memory that held the wrong
# protein, and a recalled expression pattern is the same mistake in a place
# where nothing downstream would catch it.
#
# This table is deliberately small. It covers the antigens where the question
# has a settled answer that matters clinically, and stops there; an entry
# invented to fill a gap would be worse than the gap, which at least renders
# as "not assessed".
CURATED: dict[tuple[str, str], Curated] = {
    # The case the whole feature exists for. CLL-1 sits on the CD34+CD38-
    # compartment in 77 of 89 AML patients and is completely absent from
    # CD34+CD38- cells in normal (n=11) and regenerating (n=6) marrow.
    ("CLEC12A", "hsc"): Curated(
        "absent", "flow",
        "van Rhenen et al., Blood 2007;110(7):2659-2666 "
        "(doi:10.1182/blood-2007-03-083048)",
        "Absent on CD34+CD38- cells in normal and regenerating marrow, while "
        "present on the same compartment in most AML. This is the entire "
        "argument for preferring CLEC12A to CD33, and no paralog panel can "
        "see it, because the two proteins are not relatives."),
    ("CLEC12A", "progenitor"): Curated(
        "low", "flow",
        "van Rhenen et al., Blood 2007;110(7):2659-2666",
        "Myeloid-committed cells acquire it as they mature, so the sparing "
        "argument is strongest for the stem compartment specifically."),

    # The cautionary comparison, and more nuanced than the usual summary.
    ("CD33", "hsc"): Curated(
        "contested", "flow",
        "Reviewed literature disagrees; see e.g. Haematologica 2012 and "
        "subsequent reviews of CD33 ontogeny",
        "Reported absent from pluripotent stem cells by several sources, but "
        "the published cytometry does not agree, and estimates of density on "
        "CD34+/CD33+ marrow cells vary widely. Recorded as contested rather "
        "than resolved in either direction."),
    ("CD33", "progenitor"): Curated(
        "high", "clinical",
        "Gemtuzumab ozogamicin label and trial experience; CD33-deleted "
        "donor graft programmes (trem-cel) exist specifically to avoid it",
        "Present on early multilineage progenitors with intensity falling as "
        "cells mature. This is the expression behind gemtuzumab's "
        "myelosuppression and persistent cytopenias, and it is why a "
        "CRISPR-edited CD33-negative graft is being developed to work "
        "around it.", harm="observed"),
    ("CD33", "monocyte"): Curated(
        "high", "flow", "Established myeloid ontogeny; CD33 is a Siglec of "
        "the myeloid lineage",
        "A myeloid-lineage antigen, carried through monocyte maturation.",
        harm="observed"),
    ("CD33", "neutrophil"): Curated(
        "high", "flow", "Established myeloid ontogeny",
        "Carried on granulocytes, at falling intensity as they mature.",
        harm="observed"),

    # --- The target whose selectivity was perfect and whose drug failed ---
    #
    # TIM-3 discriminates leukaemic stem cells from normal haematopoietic
    # stem cells WITHIN the same CD34+CD38- fraction, which is as clean a
    # therapeutic window as this axis can describe. The phase 3 trial still
    # missed its primary endpoint and the programme was discontinued. That
    # pairing is the most useful thing in this table: everything this tool
    # measures can come out right and the drug can still not help anybody.
    ("HAVCR2", "hsc"): Curated(
        "absent", "flow",
        "Kikushige et al., Cell Stem Cell 2010; reviewed in the sabatolimab "
        "programme rationale",
        "Highly expressed on the large majority of leukaemic stem cells and "
        "not on normal haematopoietic stem cells, discriminating the two "
        "inside the same CD34+CD38- gate. Note what followed: the phase 3 "
        "STIMULUS-MDS2 trial of sabatolimab plus azacitidine did not meet "
        "its primary overall-survival endpoint and the programme was "
        "discontinued in 2024. A target can satisfy every selectivity "
        "argument on this page and still not produce a drug that helps "
        "patients.", harm="not_observed"),
    ("HAVCR2", "t_cell"): Curated(
        "high", "flow",
        "TIM-3 is an immune checkpoint receptor on T cells",
        "TIM-3 is a checkpoint receptor, so T cells carry it by definition. "
        "For a T-cell-engaging format this is self-defeating; for a "
        "checkpoint-blocking one it is the intended mechanism. Which of "
        "those applies depends on the molecule, not on the target."),
    ("HAVCR2", "nk_cell"): Curated(
        "high", "flow", "TIM-3 is expressed on NK cells and myeloid cells",
        "Carried on NK cells, the effector arm of antibody-dependent "
        "killing."),

    # --- CD38: the gate that defines the stem compartment ----------------
    ("CD38", "hsc"): Curated(
        "absent", "flow",
        "Standard CD34+CD38- immunophenotyping of the primitive "
        "haematopoietic compartment",
        "Absent by the operational definition of the compartment: the most "
        "primitive human stem cells are identified AS the CD38-negative "
        "fraction of CD34+ cells. That is partly circular and worth saying "
        "out loud -- it is a statement about how the population is defined "
        "as much as about the protein -- but it does mean a CD38-directed "
        "binder does not engage the cells that reconstitute marrow."),
    ("CD38", "progenitor"): Curated(
        "high", "flow",
        "CD34+CD38+ is the committed progenitor compartment; extensive "
        "daratumumab clinical experience",
        "Committed progenitors are CD38-positive by the same definition. "
        "Against that, anti-CD38 therapy is given for years in myeloma "
        "without marrow failure, so the expression is real and the "
        "consequence for the progenitor compartment is not what the "
        "expression alone would predict.", harm="not_observed"),
    ("CD38", "plasma_cell"): Curated(
        "high", "flow", "The target population in multiple myeloma",
        "The intended target."),
    ("CD38", "nk_cell"): Curated(
        "high", "clinical", "Daratumumab clinical experience",
        "NK cells carry CD38 and are depleted by anti-CD38 therapy, which "
        "removes part of the effector arm the antibody itself relies on.",
        harm="observed"),
    ("CD38", "erythroid"): Curated(
        "low", "clinical",
        "Daratumumab interference with pre-transfusion testing; mitigated "
        "by dithiothreitol treatment of reagent red cells",
        "Red cells carry CD38 at low level -- enough that daratumumab in "
        "patient plasma binds reagent red cells and produces false-positive "
        "antibody screens, which blood banks work around with DTT. ABO and "
        "Rh(D) typing are unaffected and no clinically significant "
        "haemolysis has been attributed to it, so this is an interference "
        "with laboratory testing rather than destruction of the lineage.",
        harm="not_observed"),

    # --- CD47: expressed on the population it does not destroy -----------
    ("CD47", "hsc"): Curated(
        "high", "flow",
        "Jaiswal et al., Cell 2009; subsequent AML stem-cell work",
        "Carried on haematopoietic stem cells -- CD47 is the 'do not eat "
        "me' signal on essentially all cells. The window is quantitative: "
        "AML stem cells carry more of it than normal marrow stem cells, and "
        "blocking it produced no increase in phagocytosis of normal stem or "
        "progenitor cells. This is the clearest case in this table where "
        "expression on a population does not mean destruction of it.",
        harm="not_observed"),
    ("CD47", "platelet"): Curated(
        "high", "flow",
        "CD47 is expressed on platelets and erythrocytes as a self marker",
        "Platelets carry it. Thrombocytopenia has been seen with agents in "
        "this class, so unlike the stem-cell row this one is not reassured "
        "by experience."),

    # Ubiquitous, and the erythroid row is the clinically decisive one.
    ("CD47", "erythroid"): Curated(
        "high", "clinical",
        "Magrolimab phase 1 experience in higher-risk MDS; CD47 loss shortens "
        "red cell half-life in preclinical models",
        "On-target anaemia, not a theoretical risk. Managed by a low priming "
        "dose that clears older red cells first, and accompanied by near-"
        "complete loss of red-cell CD47 without the same loss on white "
        "cells. Programmes targeting SIRPa instead exist to avoid it.",
        harm="observed"),

    # The case that proves the transcript screen cannot clear a target.
    ("GPRC5D", "plasma_cell"): Curated(
        "high", "flow",
        "Talquetamab development programme; GPRC5D is the myeloma target",
        "The intended target population.", harm="observed"),
    ("GPRC5D", "b_cell"): Curated(
        "low", "flow",
        "Immunohistochemistry and expression profiling in the GPRC5D "
        "myeloma literature",
        "Marrow plasma cells carry 500- to 1000-fold more GPRC5D mRNA than "
        "peripheral blood B cells, and surface expression is markedly "
        "higher on plasma cells than on B, T, NK or NKT cells. That is a "
        "large difference, not an absence, and it is recorded as low rather "
        "than absent because the published work measures a ratio. "
        "Deliberately absent from this table: a stem-cell row. Nothing "
        "found tested GPRC5D on normal haematopoietic stem cells -- the "
        "stem-cell exclusion that is easy to misattribute here belongs to "
        "BCMA -- so that cell stays unknown."),

    # Plasma-cell restriction within blood, which is why anti-CD20 does not
    # drop existing antibody titres the way a plasma-cell-directed agent does.
    ("MS4A1", "plasma_cell"): Curated(
        "absent", "flow",
        "Established B-cell ontogeny; CD20 is lost on terminal plasma-cell "
        "differentiation",
        "Lost as B cells become plasma cells, so existing titres largely "
        "survive anti-CD20 therapy even though new responses do not."),
    ("MS4A1", "hsc"): Curated(
        "absent", "flow",
        "Established B-cell ontogeny: CD20 appears at the pre-B stage",
        "B-lineage restricted and acquired well after commitment."),
    ("CD19", "hsc"): Curated(
        "absent", "flow",
        "Established B-cell ontogeny: CD19 appears at the pro-B stage",
        "B-lineage restricted, which is why CD19-directed therapy causes "
        "recoverable B-cell aplasia rather than marrow failure."),
    ("CD19", "b_cell"): Curated(
        "high", "clinical",
        "CD19-directed CAR-T experience; B-cell aplasia is the expected "
        "on-target effect",
        "Carried on every cell committed to the B lineage, which is exactly "
        "why CD19-directed therapy produces B-cell aplasia that can persist "
        "for years. Expected and managed, not a surprise.", harm="observed"),
    ("CD19", "plasma_cell"): Curated(
        "low", "clinical",
        "Persistence of antibody titres during CD19-directed B-cell aplasia",
        "A CD19-negative long-lived plasma cell population survives CD19 "
        "aplasia and keeps producing antibody, including alloantibody. The "
        "evidence here is indirect -- it is inferred from titres persisting "
        "rather than from cytometry on the cells -- so it is recorded as "
        "low rather than absent."),
    ("CD22", "hsc"): Curated(
        "absent", "flow",
        "Established B-cell ontogeny; CD22 is restricted to the B lineage "
        "in normal tissue",
        "Absent from the multipotent stem compartment. See the progenitor "
        "row, which is less simple than this one."),
    ("CD22", "progenitor"): Curated(
        "low", "flow",
        "CD22 precedes CD19 in B-cell development; CD34+CD19-CD22+ "
        "progenitors are described",
        "The tidy summary 'B-lineage restricted' hides something here: CD22 "
        "appears BEFORE CD19 during development, on CD34+CD19-CD22+ cells. "
        "Those are B-committed rather than multipotent, so the stem "
        "compartment is still spared, but a CD34+ progenitor population "
        "does carry this antigen and a sparing argument that leans on "
        "'restricted to B cells' is leaning on something imprecise."),
    ("CD22", "b_cell"): Curated(
        "high", "clinical",
        "CD22-directed CAR-T experience; normal tissue expression is "
        "restricted to the B lineage",
        "The target population, and the source of B-cell aplasia on the "
        "same terms as CD19.", harm="observed"),
    ("MS4A1", "b_cell"): Curated(
        "high", "clinical",
        "Rituximab and the wider anti-CD20 clinical record",
        "The target population. Depletion is expected, recoverable, and "
        "spares existing antibody titres because plasma cells have lost "
        "CD20 by then.", harm="observed"),
    ("TNFRSF17", "hsc"): Curated(
        "absent", "flow",
        "BCMA expression profiling: restricted to mature B lymphocytes, "
        "plasmablasts and plasma cells, not detected on naive B cells or "
        "haematopoietic stem cells",
        "Absent from the stem compartment, which is why BCMA-directed "
        "therapy spares marrow reconstitution. This is the stem-cell "
        "exclusion that is sometimes misattributed to GPRC5D."),
    ("TNFRSF17", "b_cell"): Curated(
        "low", "flow",
        "BCMA is acquired on late-stage B cells and plasma cells, and is "
        "not detected on naive B cells",
        "Largely spared on naive and memory B cells."),
    ("TNFRSF17", "plasma_cell"): Curated(
        "high", "clinical",
        "BCMA-directed bispecific and CAR-T experience in myeloma",
        "The intended target, and the reason these agents cause "
        "hypogammaglobulinaemia that needs immunoglobulin replacement.",
        harm="observed"),
    ("CLEC12A", "monocyte"): Curated(
        "high", "flow",
        "van Rhenen et al., Blood 2007; CLEC12A is a myeloid-lineage "
        "C-type lectin",
        "Carried on monocytes. The CLEC12A sparing argument is about the "
        "stem compartment specifically; mature myeloid cells are not "
        "spared, and any profile that needs them should say so.",
        harm="observed"),
    ("CLEC12A", "neutrophil"): Curated(
        "high", "flow",
        "van Rhenen et al., Blood 2007; CLEC12A is carried on granulocytes",
        "Carried on granulocytes. Worth stating plainly, because CLEC12A "
        "is usually introduced as the selective option and that selectivity "
        "is against stem cells, not against mature myeloid cells.",
        harm="observed"),

    # --- KIT: a solid tumour target whose blood liability is the stem
    # --- compartment, proven by a drug that depletes it on purpose --------
    #
    # Every other entry in the solid tumour group has its liability in
    # epithelium, where a haematopoietic sparing profile has nothing to say.
    # KIT is the exception and it is the worst case in the table: CD117 is a
    # defining marker of the stem and progenitor compartment, so a KIT-
    # directed binder engages exactly the population that has to survive.
    ("KIT", "hsc"): Curated(
        "high", "clinical",
        "Briquilimab (JSP191/AMG 191) is an anti-CD117 antibody developed AS "
        "a conditioning agent: NCI drug dictionary; phase 1 in older adults "
        "with AML/MDS, ASH 2023 abstract (NCT04429191)",
        "CD117 is carried on haematopoietic stem cells, and the evidence for "
        "it is unusually direct: an anti-CD117 antibody is in trials "
        "specifically BECAUSE blocking stem cell factor depletes the stem "
        "compartment, which is what transplant conditioning is for. A drug "
        "whose entire purpose is to empty the marrow niche settles the "
        "question of whether this antigen's engagement reaches that "
        "population. For a GIST or mastocytosis programme that is the "
        "liability; for a conditioning programme it is the mechanism.",
        harm="observed"),
    ("KIT", "progenitor"): Curated(
        "high", "flow",
        "CD117 immunophenotyping of multipotent and myeloid progenitors; "
        "barzolvolimab neutropenia (see this catalogue's KIT entry)",
        "Carried through the multipotent and myeloid progenitor stages. "
        "Barzolvolimab's neutropenia and hair depigmentation are on-target "
        "effects of the same expression -- and the depigmentation is the "
        "reminder that this antigen's other normal homes, mast cells, "
        "melanocytes, germ cells and the interstitial cells of Cajal, are "
        "populations no row of this table represents at all.",
        harm="observed"),

    # --- The checkpoint receptors, where expression IS the mechanism ------
    #
    # For these the question this table asks needs restating. T cells carry
    # PD-1 by definition, and that is not a toxicity: blocking it is the
    # entire point. What matters is whether engaging the antigen DESTROYS the
    # population, and for a checkpoint inhibitor the answer was engineered in
    # on purpose. See the HARM comment above.
    ("PDCD1", "t_cell"): Curated(
        "high", "flow",
        "PD-1 is a T-cell inhibitory receptor; anti-PD-1 antibodies are "
        "built on Fc backbones chosen not to engage activating FcgammaR, "
        "precisely to avoid depleting the CD8 effectors they bind",
        "T cells carry PD-1, and for checkpoint blockade that is the "
        "mechanism rather than a liability -- the recorded harm here is "
        "'population survives' only because approved anti-PD-1 antibodies "
        "were deliberately given an Fc that does not recruit killing. This "
        "row is a fact about those molecules, not about the antigen. An ADC "
        "or a T-cell engager aimed at PD-1 would be depleting the cells it "
        "needs, and nothing on this page would tell you so.",
        harm="not_observed"),
    ("CTLA4", "t_cell"): Curated(
        "moderate", "flow",
        "CTLA-4 is constitutive on regulatory T cells and induced on "
        "activated T cells; Simpson et al., J Exp Med 2013 (Fc-dependent "
        "intratumoral Treg depletion in mice) against Sharma et al. 2018, "
        "'Anti-CTLA-4 immunotherapy does not deplete FOXP3+ regulatory T "
        "cells in human cancers'",
        "Not the whole T-cell compartment: CTLA-4 is constitutive on "
        "regulatory T cells and induced on activation, so 'moderate' is the "
        "honest level for the lineage as a whole. The harm verdict is left "
        "unknown for a specific reason rather than for want of looking. In "
        "mice, anti-CTLA-4 depends on Fc-receptor-bearing macrophages to "
        "deplete intratumoral Tregs; in human tumours a study looking for "
        "the same depletion did not find it, and the inhibitory receptor "
        "FcgammaRIIB appears to blunt it. Two species disagree, and "
        "flattening that into either answer would invent a result.",
        harm="unknown"),
    ("CD274", "monocyte"): Curated(
        "high", "flow",
        "PD-L1 is constitutive or induced on myeloid cells; Fc comparison in "
        "Front/PMC11835948 (2025): avelumab, a wild-type IgG1, lyses "
        "PD-L1-positive monocyte-derived dendritic cells via NK-cell ADCC, "
        "while Fc-mutated durvalumab and atezolizumab do not",
        "Myeloid cells carry PD-L1. Whether that matters is the clearest "
        "format dependence in this table: two approved anti-PD-L1 "
        "antibodies give opposite answers on the same cells in the same "
        "assay, separated by an engineered Fc. Recorded as 'survives' for "
        "the Fc-silenced molecules, which is what that row means and all it "
        "means.",
        harm="not_observed"),
    ("CD274", "dendritic"): Curated(
        "high", "flow",
        "Same Fc comparison on monocyte-derived dendritic cells from healthy "
        "donors (PMC11835948, 2025); the assay is in vitro chromium release",
        "Dendritic cells from healthy donors carry PD-L1 by cytometry, and "
        "they are the cells the Fc comparison was actually run on. Losing "
        "them would blunt the antigen presentation a checkpoint inhibitor is "
        "trying to recruit, which makes this the one row where the format "
        "choice and the mechanism point the same way.",
        harm="not_observed"),
    ("NT5E", "b_cell"): Curated(
        "high", "flow",
        "Reviewed human cytometry: CD73 is expressed on most B cells and "
        "certainly on all mature naive B cells",
        "The main home of CD73 in normal human blood. B-cell CD73 is low at "
        "birth and reaches adult levels in the first six months of life, so "
        "this is a developmentally acquired pattern rather than a fixed one."),
    ("NT5E", "t_cell"): Curated(
        "low", "flow",
        "Reviewed human cytometry: CD73 is rare on T cells, found mainly on "
        "naive CD8 cells and most abundantly on CD4+CD25+ regulatory T cells",
        "Low across the T compartment as a whole, and concentrated in two "
        "subsets a lineage-level row cannot separate. One detail cuts "
        "against reading this as settled: human conventional T cells LOSE "
        "surface CD73 on activation, so what a resting panel measures is not "
        "what a binder would meet in an inflamed tumour. Mouse expression is "
        "wider -- most T cells, NK cells, peritoneal macrophages -- and a "
        "mouse model of this target is therefore not reading the same "
        "biology."),
    ("VSIR", "monocyte"): Curated(
        "high", "flow",
        "Flow cytometry on healthy human whole blood and PBMC disclosed in "
        "the anti-VISTA antibody patent family (e.g. US 11,009,509): "
        "highest VISTA on monocytes, then neutrophils, low on CD4 and CD8 "
        "T cells",
        "Monocytes carry the most VISTA of any population in healthy blood. "
        "Note the source: these are near-identical disclosures across one "
        "family of antibody patents, not independent peer-reviewed "
        "cytometry, so the ranking between monocytes and neutrophils should "
        "be treated as one laboratory's result. The direction is corroborated "
        "by the atlas's own single-cell classification of VSIR as neutrophil- "
        "and monocyte-enhanced, and by the first anti-VISTA antibody's "
        "severe cytokine release, which is what engaging a dense myeloid "
        "antigen looks like."),
    ("VSIR", "neutrophil"): Curated(
        "high", "flow",
        "Same anti-VISTA patent-family flow data; HPA single-cell "
        "classification of VSIR as neutrophil-enhanced",
        "Second only to monocytes in healthy blood, on the same single "
        "source. Reported here rather than omitted because a dense myeloid "
        "antigen is exactly what a neutropenic patient's profile needs to "
        "see, and because the programme that engaged it stopped for cytokine "
        "release."),
}

# Non-haematopoietic facts, keyed by (gene, site). Carried separately from the
# lineage table because they are not about blood at all.
#
# Keyed by SITE as well as gene since the solid tumour antigens were added. A
# single cell per gene could hold only one note, and these antigens do not
# have one liability each: PSMA reaches salivary gland and kidney and
# duodenum, and EGFR's rash and its hypomagnesaemia are different organs with
# different mechanisms and different severities. Collapsing them would be the
# same mistake as the first draft's single stem/progenitor bucket.
#
# `site` is a free label, not a TISSUE_SCREEN code. Some of the most important
# entries here name tissue the screen has no column for at all -- there is no
# mesothelium column, and no nail matrix -- and assess_tissues says so on the
# row rather than letting the absence pass as agreement.
#
# --- Why this table matters more for solid tumour antigens ----------------
#
# For a blood cancer antigen the lineage table carries the argument and this
# one is a backstop. For a solid tumour antigen it is the reverse. HER2's
# liability is cardiomyocytes; nectin-4's is skin; mesothelin's is the lining
# of the pleura. A flow panel on a marrow specimen cannot speak to any of
# them, so a sparing profile that comes back with no haematopoietic conflicts
# has not cleared these targets -- it has asked a different question. The
# specimen assessment says this out loud; see specimen.assess_target.
CURATED_TISSUE: dict[tuple[str, str], Curated] = {
    ("GPRC5D", "keratinised structures"): Curated(
        "high", "clinical",
        "Talquetamab dermatologic reviews (nail changes in ~75-79% of "
        "patients, onychomadesis ~57%) and IMWG immunotherapy registry "
        "(nail 75.5%, oral 62.9%, dysgeusia 60.4%, skin 39.7%)",
        "Expressed in hard keratinised structures -- the upper nail matrix "
        "and keratogenous zone, the cortical cells of the hair shaft, and "
        "the filiform papillae of the tongue. The transcript screen catches "
        "part of this and not the rest: skin reads 26.4 nTPM, while tongue "
        "reads 0.0. Grade 3+ events stayed at or below 4%, so this is a "
        "quality-of-life toxicity rather than a dose-limiting one, but it "
        "affects most patients treated.",
        harm="observed"),

    # --- HER2: the case where the transcript screen actually works --------
    ("ERBB2", "heart muscle"): Curated(
        "high", "clinical",
        "Trastuzumab cardiotoxicity literature; CHF in ~2.5% across eight "
        "RCTs, ~0.92% severe in anthracycline-naive early breast cancer; "
        "cardiac ERBB2 knockout mice develop progressive ventricular "
        "dilation over months",
        "Cardiomyocytes carry HER2 and use it: neuregulin-1 signalling "
        "through HER2 is a survival pathway, and blocking it is the leading "
        "explanation for the cardiomyopathy. Worth recording as the one "
        "place the always-on screen does its job unaided -- heart muscle "
        "reads 53.9 nTPM, in the top band, with no curation needed to raise "
        "the flag. Set against GPRC5D, which the same screen half-misses, "
        "that is the useful comparison: the screen is not useless and it is "
        "not sufficient, and you cannot tell which case you are in from the "
        "number alone. Mechanism still debated, and the risk is strongly "
        "modified by prior anthracycline exposure.",
        harm="observed"),

    # --- EGFR: two organs, two mechanisms, two severities -----------------
    ("EGFR", "skin"): Curated(
        "high", "clinical",
        "Cetuximab and panitumumab labels and reviews: acneiform rash in "
        "80-95% of patients on monotherapy",
        "EGFR in adult skin concentrates in the basal keratinocytes, the "
        "most proliferative epithelial population, which is where it "
        "regulates normal skin growth. The rash follows in the large "
        "majority of patients. The precise mechanism -- keratinocyte "
        "alteration, follicular proliferation and maturation -- is not "
        "settled, so what is established here is the association and the "
        "expression, not the causal chain.",
        harm="observed"),
    ("EGFR", "kidney"): Curated(
        "moderate", "clinical",
        "Cetuximab label (hypomagnesaemia in >10% of patients); EGF-"
        "dependent activation of the TRPM6 magnesium channel in the distal "
        "convoluted tubule",
        "A different organ and a cleaner mechanism than the rash: EGF "
        "normally stimulates TRPM6, the magnesium channel of the distal "
        "convoluted tubule, and blocking the receptor leaves it "
        "under-stimulated. The result is hypomagnesaemia rather than cell "
        "death -- a transport function lost, with the cells intact. Kept as "
        "a separate row because 'EGFR is expressed in kidney' and 'EGFR is "
        "expressed in skin' lead to different problems, and one note could "
        "not carry both.",
        harm="observed"),

    # --- PSMA: and the column that was missing until this was written -----
    ("FOLH1", "salivary gland"): Curated(
        "high", "clinical",
        "Lu-177-PSMA-617 radioligand experience: grade 1-2 dry mouth in 39% "
        "(VISION) and 60% (TheraP); PSMA localises to salivary acinar cells",
        "PSMA is carried by the acinar cells of the parotid and "
        "submandibular glands, and radioligand uptake there is visible on "
        "every PSMA PET scan. Radiation to the acini can be irreversible in "
        "some patients. This row is also where the modality matters most in "
        "this table: the toxicity is dose-deposited radiation, so a cold "
        "antibody against the same epitope would not reproduce it, and "
        "nothing in the level or the evidence column tells you that.",
        harm="observed"),
    ("FOLH1", "duodenum"): Curated(
        "high", "ihc",
        "Human Protein Atlas transcript screen, duodenum column added "
        "2026-11-09; PSMA expression in gastrointestinal mucosa is "
        "described in the radioligand literature",
        "Duodenum reads 250.9 nTPM -- the highest value for this gene in "
        "any screened tissue, above salivary gland at 49.8 and kidney at "
        "12.1. The screen could not see it at all until the column was "
        "added, because the earlier twelve covered colon and oesophagus and "
        "stopped. The lesson is about the screen rather than about PSMA: a "
        "tissue panel is a list somebody chose, and a gene's loudest site "
        "can sit just outside it and read as silence. No bowel toxicity "
        "matching this signal is established, which is why the harm verdict "
        "is left open.",
        harm="unknown"),
    ("FOLH1", "kidney"): Curated(
        "high", "clinical",
        "PSMA radioligand dosimetry: kidney uptake is dose-limiting "
        "alongside xerostomia, with potential for nephrotoxicity",
        "Proximal tubule uptake is one of the two dose-limiting sites for "
        "PSMA radioligands. Recorded as a concern with harm left unknown: "
        "the uptake and the dose constraint are established, while frank "
        "nephrotoxicity is described as a potential rather than an observed "
        "outcome, and overstating it here would be inventing a result in the "
        "direction that feels safer."),

    # --- Mesothelin: a tissue with no column in the screen ---------------
    ("MSLN", "mesothelium (pleura, pericardium, peritoneum)"): Curated(
        "high", "clinical",
        "Phase 1 mesothelin CAR-T: pleuritis, pericarditis and peritonitis "
        "as on-target off-tumour events; SS1P immunotoxin pleuritis; severe "
        "pulmonary toxicity in a high-dose cohort (2023)",
        "In normal adults mesothelin is essentially confined to mesothelium "
        "-- which sounds like excellent selectivity until you notice that "
        "the mesothelium lines the lung, the heart and the abdomen. "
        "Engaging it inflames those linings, and the events are named after "
        "them: pleuritis, pericarditis, peritonitis. Severity tracks CAR "
        "affinity and dose rather than the presence or absence of the "
        "antigen, and affinity-tuned CARs were more selective for "
        "mesothelin-high tumour. No column of the tissue screen samples "
        "mesothelium, so this entire liability is invisible to the "
        "automated half of this page.",
        harm="observed"),

    # --- CLDN18: the screen's loudest flag is the isoform the drug avoids -
    ("CLDN18", "stomach"): Curated(
        "high", "clinical",
        "GLOW (zolbetuximab + CAPOX): nausea 68.5% vs 50.2% placebo, "
        "vomiting 66.1% vs 30.9%, grade 3+ vomiting 12.2%; pooled "
        "SPOTLIGHT + GLOW (n=533): nausea 76.0%, vomiting 66.8%",
        "The mirror image of GPRC5D, and the more dangerous failure of the "
        "two. CLDN18 has two splice variants: CLDN18.1 in lung alveolar "
        "epithelium, CLDN18.2 in differentiated gastric mucosa. Zolbetuximab "
        "binds the first extracellular domain of 18.2 and does not bind "
        "18.1. The transcript screen reads the GENE, so before stomach was "
        "added it reported 237.8 nTPM in lung -- a loud flag on the one "
        "isoform the drug is specific against -- and nothing at all about "
        "the stomach, which now reads 684.2 and is where the emesis comes "
        "from. A screen can be precise, well-formed and pointed at the wrong "
        "molecule. Adding the column fixes the blind spot; the isoform "
        "blindness is not fixable at transcript level and remains true of "
        "every row of the automated screen.",
        harm="observed"),
    ("CLDN18", "gastric stem cell zone"): Curated(
        "absent", "ihc",
        "Immunohistochemistry in normal gastric tissue: CLDN18.2 sits at the "
        "apical side of tight junctions in differentiated epithelial cells "
        "and is absent from the gastric stem cell zone",
        "The reason this target is tolerable at all, and a direct echo of "
        "the CLEC12A argument in a different organ. CLDN18.2 is confined to "
        "differentiated cells and absent from the gastric stem cell zone, so "
        "the mucosa can regenerate: the emesis is severe in cycle 1, under "
        "10% by cycle 2, and self-limiting. Sparing the stem compartment is "
        "what makes an on-target toxicity recoverable, in epithelium exactly "
        "as in marrow. Its apical localisation behind tight junctions also "
        "limits what a circulating binder can reach, which is a protection "
        "no expression number describes."),

    # --- Nectin-4: expression in skin, and a boxed warning ---------------
    ("NECTIN4", "skin"): Curated(
        "high", "clinical",
        "Enfortumab vedotin label (boxed warning for severe and fatal "
        "cutaneous reactions including SJS and TEN); EV-301: skin reactions "
        "the most common grade 3+ event of special interest, 14.5% vs 0.7% "
        "with chemotherapy; GLP human tissue cross-reactivity study confirms "
        "binding to skin structures",
        "Nectin-4 is carried by epidermal keratinocytes and skin "
        "appendages, the cross-reactivity study confirms the drug binds "
        "them, and the transcript screen agrees loudly at 110.2 nTPM in "
        "skin. So the expression, the binding and the toxicity all line up "
        "-- and the sources still stop short of showing that normal-skin "
        "nectin-4 causes the SJS and TEN specifically, as opposed to the "
        "milder reactions. Recorded as the strongest kind of flag this table "
        "can raise, with that gap stated rather than closed by assumption."),

    # --- TROP2: the screen's widest flag, and a toxicity from elsewhere ---
    ("TACSTD2", "epithelium (broad)"): Curated(
        "high", "animal",
        "Cynomolgus toxicology for sacituzumab govitecan: neutropenia and "
        "diarrhoea at higher doses, with evidence for significant "
        "histopathological change in TROP2-expressing normal tissues "
        "lacking; phase 1 dose-limiting toxicity was neutropenia",
        "TROP2 draws more transcript flags than anything else in this "
        "catalogue: oesophagus 1245, salivary gland 585, skin 555, urinary "
        "bladder 219, kidney 186, breast 127, pancreas 107. By the logic of "
        "the rest of this page that is the worst profile here. And the "
        "dose-limiting toxicities -- neutropenia and diarrhoea -- are SN-38 "
        "effects, the same ones irinotecan produces without any antibody at "
        "all, while the monkey studies found no histopathological change in "
        "the normal tissues that carry the antigen. Recorded as 'survives' "
        "on animal evidence, which by design cannot clear a human population "
        "and does not here. The real lesson is that this page's axis can be "
        "maximally alarmed and still be looking in the wrong place, because "
        "the payload has toxicity of its own that no expression table "
        "contains.",
        harm="not_observed"),

    # --- CEACAM5: the one target where both axes are live ----------------
    ("CEACAM5", "colon"): Curated(
        "high", "ihc",
        "Western blotting of normal human colon lysates shows small amounts "
        "of CEACAM1, CEACAM5 and CEACAM6 against strong CEACAM5 in 15 of 16 "
        "colon cancer lysates; CEACAM7, not CEACAM5, is the one confined to "
        "the apical surface of differentiated colorectal epithelium",
        "Normal colon carries CEACAM5, at low level against a very high "
        "transcript signal (920.3 nTPM), and the window is a quantity "
        "difference rather than an absence. Two things to keep straight. "
        "First, the apical confinement that protects a circulating binder "
        "belongs to CEACAM7, not to CEACAM5, so that reassurance does not "
        "transfer. Second, the granulocyte expression often attached to "
        "'CEA' belongs to CEACAM8 and CEACAM3, which are granulocyte-"
        "restricted, and to CEACAM6; CEACAM5 itself is mostly epithelial. "
        "That is why this table has no neutrophil row for CEACAM5 -- the "
        "evidence names its relatives, not it. It also makes CEACAM5 the "
        "one catalogue target where this page and the paralog panel are "
        "both live at once: the family really are close relatives, so "
        "confusability is a genuine risk here in a way it never was for "
        "CLEC12A against CD33.",
        harm="unknown"),

    # --- VEGFA: a secreted ligand, where the question is mis-shaped ------
    ("VEGFA", "every screened tissue (secreted ligand)"): Curated(
        "high", "clinical",
        "Bevacizumab class effects: hypertension, proteinuria, impaired "
        "wound healing and haemorrhage, from blocking a normal "
        "physiological signal rather than from destroying a cell population",
        "VEGFA reads in the detected band in every one of the seventeen "
        "screened tissues -- thyroid 255.6, heart 202.6, liver 202.0, "
        "urinary bladder 115.2, kidney 105.1 -- which is true and nearly "
        "useless. This is the third way this screen fails, after the false "
        "negative (GPRC5D) and the flag on the wrong isoform (CLDN18): it "
        "can be right about everything and informative about nothing. The "
        "deeper mismatch is that VEGFA is secreted. 'Which healthy cells "
        "carry it' is the wrong question for a soluble ligand; the right one "
        "is what normal physiology depends on the signal, and the answer -- "
        "vascular tone, glomerular integrity, wound healing -- is not an "
        "expression table. The catalogue holds roughly two dozen secreted "
        "targets and this axis mis-frames all of them. One row cannot fix "
        "that, and pretending the silence was an assessment would be worse.",
        harm="observed"),
}

#: Genes with at least one curated non-haematopoietic row. A target in here
#: has an established liability outside blood, which is the thing a flow
#: panel on a marrow specimen cannot reach.
CURATED_TISSUE_GENES: frozenset[str] = frozenset(g for g, _ in CURATED_TISSUE)


# --- Assessment -------------------------------------------------------------


@dataclass
class LineageCall:
    """What is known about one antigen on one population."""

    lineage: str
    label: str
    level: str
    evidence: str
    source: str = ""
    note: str = ""
    ntpm: float | None = None
    detail: str = ""
    harm: str = "unknown"

    @property
    def is_concern(self) -> bool:
        return self.level in ("low", "moderate", "high", "contested")

    @property
    def is_cleared(self) -> bool:
        """True only when something actually looked and found nothing.

        Deliberately not `level == "absent"` alone: absence reached any way
        other than by cytometry or clinical observation is not reachable in
        the first place (see band_for_ntpm), so this is belt and braces
        against a future caller inventing one.
        """
        return self.level == "absent" and decisive(self.evidence)

    def to_dict(self) -> dict:
        return {
            "lineage": self.lineage, "label": self.label, "level": self.level,
            "evidence": self.evidence,
            "evidence_label": EVIDENCE_LABEL.get(self.evidence, self.evidence),
            "source": self.source, "note": self.note, "ntpm": self.ntpm,
            "detail": self.detail, "concern": self.is_concern,
            "cleared": self.is_cleared, "harm": self.harm,
        }


@dataclass
class ExpressionData:
    """The cached transcript screen, with its provenance.

    `missing` lists columns the atlas did not return. It is kept rather than
    dropped because an unrecognised column code is silently ignored by that
    API -- asking for `t_RNA_skin` instead of `t_RNA_skin_1` returns a row
    with no skin in it and no error at all -- so a column that simply
    vanished would otherwise read as "not expressed anywhere".
    """

    values: dict[str, dict[str, float]] = field(default_factory=dict)
    missing: dict[str, list[str]] = field(default_factory=dict)
    synonyms: dict[str, list[str]] = field(default_factory=dict)
    fetched: str = ""
    source: str = ""
    version: str = ""

    def for_gene(self, gene: str) -> dict[str, float]:
        return self.values.get(gene.upper(), {})

    def has(self, gene: str) -> bool:
        return gene.upper() in self.values

    def aliases(self, gene: str) -> list[str]:
        """Other names for this antigen, CD numbers included."""
        return self.synonyms.get(gene.upper(), [])


def load_expression(path: Path | None = None) -> ExpressionData:
    p = path or EXPRESSION_PATH
    if not p.exists():
        return ExpressionData()
    doc = json.loads(p.read_text(encoding="utf-8"))
    return ExpressionData(
        values={k.upper(): v for k, v in (doc.get("values") or {}).items()},
        missing=doc.get("missing") or {},
        synonyms={k.upper(): v
                  for k, v in (doc.get("synonyms") or {}).items()},
        fetched=doc.get("fetched", ""),
        source=doc.get("source", ""),
        version=doc.get("version", ""),
    )


def _rna_level(gene: str, lineage: Lineage,
               data: ExpressionData) -> tuple[str, float | None, str]:
    """Highest transcript level across the columns making up this lineage.

    Max, not mean: a lineage is at risk if ANY of its subsets carries the
    antigen, and averaging CD19 across seven T-cell subsets would dilute a
    real signal in one of them into nothing.
    """
    cols = lineage.blood_columns + lineage.tissue_columns
    if not cols:
        return "unknown", None, ""
    vals = data.for_gene(gene)
    seen = [(c, vals[c]) for c in cols if c in vals]
    if not seen:
        return "unknown", None, ""
    col, top = max(seen, key=lambda kv: kv[1])
    pretty = col.split("_RNA_", 1)[-1].replace("_", " ")
    return band_for_ntpm(top), top, pretty


# --- Where the protein sits, and whether this axis applies at all ----------

#: Raw UniProt evidence, fetched by scripts/fetch_topology.py.
TOPOLOGY_PATH = ROOT / "data" / "lineage" / "topology.json"

TOPOLOGIES = ("surface", "shed", "secreted", "intracellular", "unknown")

#: Phrased to read as a sentence fragment, because these are substituted
#: into prose on the card: "VEGFA is <label>, so ...".
TOPOLOGY_LABEL = {
    "surface": "displayed on the cell surface",
    "shed": "displayed on the cell surface and also released in soluble form",
    "secreted": "a secreted protein rather than one displayed on a cell",
    "intracellular": "an intracellular protein, where an antibody does not go",
    "unknown": "of a topology this page has not established",
}


def classify_topology(entry: dict) -> str:
    """UniProt evidence -> where a binder would actually meet this protein.

    The lineage table asks which healthy cells CARRY a target. That question
    silently assumes the protein is displayed somewhere a binder can reach,
    and for most of this catalogue it is not. Before this existed the screen
    reported alpha-synuclein at 103.2 nTPM on dendritic cells and listed them
    as a population at risk, when SNCA is cytoplasmic and the number measures
    what those cells make. TNF read "high" on monocytes, which are the cells
    that produce it.

    The test for a surface protein here is STRUCTURAL -- a membrane-spanning
    feature, or a GPI anchor -- and deliberately not the "Cell membrane"
    keyword. That keyword is also given to proteins that merely associate
    with the membrane: using it pulled tau, a microtubule-binding protein of
    the cytosol, into the same bucket as CLEC12A.

    Four answers, because the honest handling differs for each:

      surface        the lineage question works as asked
      shed           it works, AND a soluble pool exists that acts as a sink
      secreted       "which cells carry it" is the wrong question
      intracellular  wrong question, and mostly out of reach besides
    """
    kw = set(entry.get("keywords") or ())
    surface = bool(entry.get("membrane_spans", 0)) or "GPI-anchor" in kw
    soluble = "Secreted" in kw
    if surface and soluble:
        return "shed"
    if surface:
        return "surface"
    if soluble:
        return "secreted"
    return "intracellular"


# Where the annotation is right and still misleading. UniProt records a real
# minor secreted pool for both of these, so the automated verdict is not a
# bug in the data -- it is the dominant biology being outvoted by a true
# footnote. Each override says why, as the curated expression tables do.
#
# These two are also the case where the minor pool IS the drug rationale,
# which is worth more than the correction: anti-tau and anti-synuclein
# antibodies are aimed at the small extracellular fraction that propagates
# between neurons, precisely because the large intracellular pool is beyond
# reach. That is a hard brief, and the entry should say so rather than let
# "secreted" imply an accessible target.
CURATED_TOPOLOGY: dict[str, tuple[str, str]] = {
    "SNCA": ("intracellular",
             "UniProt lists Cytoplasm, Nucleus, Synapse and axon first and "
             "Secreted among them, and the secreted fraction is real but "
             "small. The aggregates that matter are intraneuronal. An "
             "antibody can only engage the extracellular seeds that "
             "propagate between neurons, which is the rationale for the "
             "programme and also why it is a hard one."),
    "MAPT": ("intracellular",
             "Tau is a microtubule-associated protein of the cytosol and "
             "axon. UniProt also records Cell membrane and Secreted, but it "
             "spans no membrane and has no GPI anchor; the membrane entry "
             "is association, not insertion. Same rationale and same "
             "difficulty as SNCA."),
}


@dataclass
class TopologyData:
    entries: dict[str, dict] = field(default_factory=dict)
    unresolved: dict[str, str] = field(default_factory=dict)
    fetched: str = ""
    source: str = ""

    def verdict(self, gene: str) -> tuple[str, str]:
        """(topology, why). Curated overrides win, as everywhere else here."""
        g = (gene or "").upper()
        cur = CURATED_TOPOLOGY.get(g)
        if cur is not None:
            return cur
        entry = self.entries.get(g)
        if entry is None:
            return "unknown", (
                "No UniProt topology was fetched for this gene, so whether a "
                "binder could meet it on a cell surface is not established "
                "here. Reported as unknown rather than assumed.")
        verdict = classify_topology(entry)
        kw = [k for k in (entry.get("keywords") or ())
              if k in ("Secreted", "Transmembrane", "GPI-anchor",
                       "Cell membrane", "Cytoplasm", "Nucleus")]
        return verdict, (
            f"From UniProt {entry.get('accession', '')}: "
            f"{entry.get('membrane_spans', 0)} membrane-spanning feature(s), "
            f"keywords {', '.join(kw) or 'none of the topology set'}.")


def load_topology(path: Path | None = None) -> TopologyData:
    p = path or TOPOLOGY_PATH
    if not p.exists():
        return TopologyData()
    doc = json.loads(p.read_text(encoding="utf-8"))
    return TopologyData(
        entries={k.upper(): v for k, v in (doc.get("entries") or {}).items()},
        unresolved=doc.get("unresolved") or {},
        fetched=doc.get("fetched", ""),
        source=doc.get("source", ""),
    )


# --- What is lost when the signal is blocked -------------------------------

# The topology axis above finishes a sentence it could not previously end.
# For a secreted target the card said: binding it neutralises a signal rather
# than depleting a population, so the risk moves to whatever normal
# physiology depends on that signal -- and then stopped. This is the other
# half.
#
# It is a different kind of claim from everything else in this module. The
# lineage and tissue tables answer "what does this protein sit on"; this one
# answers "what is this protein FOR, and what happens to a person without
# it". Expression data cannot reach it at any resolution.
#
# The strongest evidence here is not a trial. For several of these the human
# knockout already exists: people with inherited IL-17 pathway defects get
# chronic mucocutaneous candidiasis, people with congenital complement
# deficiency get recurrent neisserial infection, and people with congenital
# factor XI deficiency largely do not bleed. A drug that blocks the signal is
# running the same experiment on purpose, and the inborn error says in
# advance what to expect.
#
# Severity is about the consequence, not about the regulatory paperwork,
# except where regulators boxing it IS the fact worth carrying:
#
#   boxed_warning   regulators put this specific consequence in a box
#   serious         documented serious harm, not boxed
#   manageable      real, common, routinely managed
#   well_tolerated  blocked in people and little is lost -- the controls
#   contested       looked at, and the published answer disagrees with itself
#   unknown         nobody has established it
#
# "contested" is here for the same reason LEVELS has it: it is a real answer
# rather than a hedge. MASP2 deficiency occurs in about 4% of Caucasians and
# healthy homozygotes are documented, while the disease database still lists
# autoimmunity and recurrent infection. Collapsing that into either
# "well_tolerated" or "unknown" would invent a consensus. It is also distinct
# from "unknown", which means nobody looked.
SEVERITIES = ("boxed_warning", "serious", "manageable", "well_tolerated",
              "contested", "unknown")

SEVERITY_LABEL = {
    "boxed_warning": "boxed warning",
    "serious": "serious harm documented",
    "manageable": "common but managed",
    "well_tolerated": "little lost by blocking it",
    "contested": "published evidence disagrees",
    "unknown": "not established",
}


@dataclass(frozen=True)
class Physiology:
    """What a signal does, and what blocking it costs."""

    signal: str
    cost: str
    evidence: str
    source: str
    severity: str = "unknown"


# Thirty entries, covering twenty-nine of the thirty-two off-axis targets
# plus TNF, which is shed rather than secreted and needs this row anyway.
#
# Three are still absent -- IL11, IL13 and SOD1 -- and the reason is worth
# more than the entries would have been. Searching turned up nothing solid
# for them, and "my search found nothing" is a fact about the search, not
# about the world. Writing "no safety consequence is established" on the
# strength of it would be exactly the nTPM-0.0-means-absent error this
# package is built around, committed one layer up. They stay empty, and the
# card says that half of the question is open.
#
# The classical and lectin pathway entries (C1QB, C1S, MASP2) were held back
# in the first pass precisely because the four alternative/terminal pathway
# entries all carry boxed warnings for encapsulated bacteria, and extending
# a class warning across a family is the CEACAM5 granulocyte mistake. The
# research vindicated the caution and then paid for it: those three turn out
# to have a DIFFERENT answer, not a weaker version of the same one. Blocking
# early classical complement points at autoimmunity rather than at
# encapsulated-organism sepsis, because that is what the inborn errors show,
# and the sutimlimab label says so in as many words.
PHYSIOLOGY: dict[str, Physiology] = {
    # --- The complement block, where the inborn error came first ---------
    "C5": Physiology(
        signal="Terminal complement. C5 cleavage builds the membrane attack "
               "complex, which is how the immune system lyses Neisseria.",
        cost="Invasive meningococcal disease. Eculizumab and ravulizumab "
             "carry a boxed warning for life-threatening and fatal "
             "meningococcal infection and are available only through a "
             "restricted REMS programme, with vaccination required at least "
             "two weeks before the first dose. Vaccination reduces the risk "
             "and does not remove it: CDC data suggest the vaccines give "
             "incomplete protection in people on eculizumab, and CDC puts "
             "the risk at up to 2,000 times that of an otherwise healthy "
             "person. Congenital terminal complement deficiency produces the "
             "same vulnerability, which is why this was predictable.",
        evidence="clinical",
        source="Soliris and Ultomiris labels (boxed warning, REMS); CDC "
               "clinical guidance on complement inhibitors",
        severity="boxed_warning"),
    "C3": Physiology(
        signal="The convergence point of all three complement pathways, and "
               "the opsonin that marks encapsulated bacteria for clearance.",
        cost="Serious infection by encapsulated organisms -- Streptococcus "
             "pneumoniae, Neisseria meningitidis, Haemophilus influenzae "
             "type b. Pegcetacoplan carries this as its boxed warning, and "
             "starting it during an unresolved such infection is "
             "contraindicated. Blocking at C3 removes opsonisation as well "
             "as the membrane attack complex, so the exposure is wider than "
             "at C5, not narrower.",
        evidence="clinical",
        source="Empaveli (pegcetacoplan) label, boxed warning",
        severity="boxed_warning"),
    "CFB": Physiology(
        signal="Factor B drives the alternative pathway amplification loop, "
               "which is most of the complement response to a bacterial "
               "surface.",
        cost="Serious infection by encapsulated bacteria. Iptacopan carries "
             "the boxed warning and a restricted REMS. The FDA review calls "
             "this a known risk of complement inhibitors acting proximally "
             "in the cascade -- so it is a property of the position in the "
             "pathway rather than of the molecule, and the label is explicit "
             "that patients remain at increased risk even after developing "
             "antibodies to vaccination.",
        evidence="clinical",
        source="Fabhalta (iptacopan) label and FDA risk review",
        severity="boxed_warning"),
    "CFD": Physiology(
        signal="Factor D is the rate-limiting enzyme of the alternative "
               "pathway.",
        cost="Serious infection by encapsulated bacteria; danicopan carries "
             "the same boxed warning and its own REMS, and is contraindicated "
             "in an unresolved such infection.",
        evidence="clinical",
        source="Voydeya (danicopan) label, boxed warning and REMS",
        severity="boxed_warning"),

    # --- Signals whose loss shows up as infection ------------------------
    "TNF": Physiology(
        signal="TNF holds the tuberculous granuloma together. It is not "
               "merely inflammatory: the structure that walls off "
               "Mycobacterium tuberculosis depends on it.",
        cost="Reactivation of latent tuberculosis, often disseminated. "
             "Reported relative risk ranges from about 1.6 to 25 times the "
             "general population, with the monoclonal antibodies at the high "
             "end and etanercept lower -- the same split as the "
             "transmembrane-TNF binding. Latent TB screening before and "
             "during therapy is standard and does not abrogate the risk: "
             "tuberculin and IGRA testing can be unreliable in these "
             "patients, and miliary TB after a negative screen is reported.",
        evidence="clinical",
        source="Anti-TNF labels (boxed warning for serious infections); CDC "
               "guidance on LTBI screening before and during TNFI therapy",
        severity="boxed_warning"),
    "IL17A": Physiology(
        signal="IL-17A is central to mucosal antifungal defence.",
        cost="Mucocutaneous candidiasis, mostly oral and mostly mild. Pooled "
             "trial rates run about 1.7-4.0% against 0.3% on placebo, with "
             "real-world cohorts higher (2.68% secukinumab, 6.6% "
             "ixekizumab). The mechanism is not inferred from the drugs: "
             "people born with IL-17 pathway defects develop chronic "
             "mucocutaneous candidiasis, so the inborn error predicted the "
             "adverse event. Separately, IL-17 blockade is associated with "
             "exacerbation and new onset of inflammatory bowel disease, "
             "which is why a history of IBD is a caution -- that signal "
             "comes from pharmacovigilance and case data and should be read "
             "as a signal rather than a rate.",
        evidence="clinical",
        source="Pooled anti-IL-17 psoriasis trial reviews; Spanish "
               "real-world cohort; FAERS pharmacovigilance analysis",
        severity="manageable"),

    # --- Signals whose loss shows up as tissue failure -------------------
    "IL6": Physiology(
        signal="IL-6 drives the acute phase response and supports repair of "
               "intestinal mucosa.",
        cost="Lower gastrointestinal perforation, reported at roughly 1-2 "
             "per 1,000 patient-years and concentrated in patients with "
             "prior diverticulitis. The second half matters more than the "
             "first: blocking IL-6 also abolishes the CRP rise and blunts "
             "the pain, so these patients present with mild, non-specific "
             "symptoms rather than an acute abdomen, and the laboratory "
             "marker a clinician would reach for has been removed by the "
             "drug. The harm is partly that the warning signs are gone.",
        evidence="clinical",
        source="Tocilizumab real-world cohorts and case series on lower "
               "intestinal perforation; CRP suppression is a labelled effect",
        severity="serious"),
    "VEGFA": Physiology(
        signal="VEGF maintains vascular homeostasis, glomerular filtration "
               "integrity, blood pressure regulation and wound healing. It "
               "is a normal physiological signal, not a tumour product.",
        cost="Hypertension, proteinuria, impaired wound healing, haemorrhage "
             "and gastrointestinal perforation (0.3-3.2% on bevacizumab "
             "labelling). Grade 3+ proteinuria runs about 2.2%, roughly five "
             "times chemotherapy alone, and the mechanism is explicitly "
             "compared to pre-eclampsia, where VEGF is low. Surgery needs a "
             "28-day window either side because healing is impaired for as "
             "long as the antibody persists.",
        evidence="clinical",
        source="Bevacizumab labelling and toxicity reviews",
        severity="serious"),
    "NGF": Physiology(
        signal="Nerve growth factor carries nociception from the joint. Pain "
               "is the signal that stops a person loading a damaged joint.",
        cost="Rapidly progressive osteoarthritis. An FDA advisory committee "
             "voted 19-1 against tanezumab, with RPOA the central objection "
             "and the agency finding the proposed REMS insufficient to "
             "mitigate it. Adjudicated RPOA ran about 1.3% in one pooled "
             "analysis and over 6% at the 5 mg dose in a phase 3 report; of "
             "those reaching the more severe RPOA2, around 60% went on to "
             "total joint replacement. The FDA also found that stopping the "
             "drug does not appear to prevent further joint damage. This is "
             "the clearest case in the table of a drug working exactly as "
             "designed and the loss of the signal being the harm.",
        evidence="clinical",
        source="FDA advisory committee briefing and review, tanezumab (2021)",
        severity="serious"),
    "SOST": Physiology(
        signal="Sclerostin restrains bone formation; blocking it builds "
               "bone. What else the Wnt pathway does in the vasculature is "
               "the open question.",
        cost="Romosozumab carries a boxed warning for myocardial infarction, "
             "stroke and cardiovascular death, and is not to be started "
             "within a year of either event. The evidence behind the box is "
             "genuinely contested and this entry should not pretend "
             "otherwise: the imbalance appeared in ARCH against alendronate "
             "(2.5% vs 1.9% MACE in year one) and NOT in the "
             "placebo-controlled FRAME trial, which raises the possibility "
             "that what ARCH measured was alendronate protecting the "
             "comparator arm rather than romosozumab harming its own. Later "
             "retrospective work in patients with prior cardiovascular "
             "disease found higher risk. Recorded as boxed because it is, "
             "and as contested because it is that too.",
        evidence="clinical",
        source="Evenity label (boxed warning); ARCH and FRAME trial results; "
               "subsequent retrospective cohorts",
        severity="boxed_warning"),

    # --- The controls. A table where everything is dangerous teaches
    # --- nothing, and these two are the reason to believe the rest.
    "PCSK9": Physiology(
        signal="PCSK9 degrades the LDL receptor. Blocking it leaves more "
               "receptor on hepatocytes and clears more LDL.",
        cost="Very little that has been established. People carrying "
             "homozygous loss-of-function variants live with LDL-C below 15 "
             "mg/dL and generally show no neurocognitive deficit, diabetes, "
             "cataract or stroke -- the lifelong human knockout is healthy, "
             "which is why this target was pursued. Regulators required no "
             "special safety monitoring for the class. The honest caveat is "
             "duration: randomised data beyond ten years at very low LDL-C "
             "are limited, and one pooled alirocumab analysis noted a "
             "cataract signal.",
        evidence="clinical",
        source="PCSK9 loss-of-function genetic studies; pooled alirocumab "
               "safety analyses; class labelling",
        severity="well_tolerated"),
    "F11": Physiology(
        signal="Factor XI amplifies thrombin generation in pathological "
               "clotting but contributes little to the haemostasis that "
               "stops bleeding from a wound.",
        cost="Strikingly little, and that is the entire thesis of the drug "
             "class. People with congenital factor XI deficiency have little "
             "or no spontaneous bleeding, and abelacimab produced a large "
             "reduction in major bleeding against rivaroxaban in "
             "AZALEA-TIMI 71. The catch belongs on the same row: the "
             "asundexian phase 3 OCEANIC-AF trial was stopped early because "
             "it showed superior safety and failed to show non-inferior "
             "efficacy. Separating a signal from haemostasis may mean "
             "separating it from the clotting you wanted to stop, and the "
             "safety and the failure may be one fact seen twice.",
        evidence="clinical",
        source="Congenital FXI deficiency literature; AZALEA-TIMI 71; "
               "OCEANIC-AF early termination",
        severity="well_tolerated"),

    # --- An entry that exists to say nobody knows -------------------------
    "IL5": Physiology(
        signal="IL-5 drives eosinophil maturation and survival. Eosinophils "
               "are prominent in the response to helminths.",
        cost="Not established. The concern is theoretical and the sponsors "
             "say so: completed phase 3 trials reported no apparent effect "
             "on helminth infection incidence, while a WHO VigiBase "
             "disproportionality analysis reported parasitic infections more "
             "often with benralizumab than comparators (0.23% vs 0.07% for "
             "mepolizumab). Spontaneous reports cannot establish causation, "
             "and trials excluded untreated parasitic infection, so the "
             "population most at risk was never enrolled. Recorded as "
             "unknown rather than well tolerated for that reason.",
        evidence="animal",
        source="Benralizumab and mepolizumab trial protocols; WHO VigiBase "
               "disproportionality analysis",
        severity="unknown"),

    # --- Early classical complement: a different answer, not a weaker one --
    "C1QB": Physiology(
        signal="C1q starts the classical pathway and, as much to the point, "
               "opsonises apoptotic cells and immune complexes so they are "
               "cleared before they can be presented as self-antigen.",
        cost="Autoimmunity rather than sepsis. Hereditary C1q deficiency is "
             "among the strongest monogenic causes of systemic lupus "
             "erythematosus known -- one review puts the SLE prevalence in "
             "C1q-deficient individuals at about 90%, ahead of C4 and "
             "C1r/C1s -- and a fraction of those patients are also "
             "susceptible to encapsulated bacteria. Penetrance and "
             "expressivity are undetermined because the cases are so few. "
             "Note what this row is: the consequence of losing the signal, "
             "established by the inborn error. Whether a C1q-blocking drug "
             "reproduces it over a treatment course is a separate question "
             "and is not established.",
        evidence="clinical",
        source="Reviews of hereditary C1q deficiency and monogenic SLE",
        severity="serious"),
    "C1S": Physiology(
        signal="C1s is the protease of the C1 complex, so blocking it stops "
               "the classical pathway while leaving the alternative and "
               "lectin pathways running.",
        cost="Two things, and the second is the one the inborn errors "
             "predict. Sutimlimab's label carries serious infection by "
             "encapsulated bacteria under Warnings and Precautions -- NOT a "
             "boxed warning, unlike the alternative and terminal pathway "
             "agents -- and records that life-threatening and fatal "
             "infections have occurred in vaccinated and unvaccinated "
             "patients alike. It then warns that, based on its mechanism of "
             "action, the drug may increase the risk of developing "
             "autoimmune disease such as SLE. A label naming the inborn "
             "error's phenotype as a predicted drug risk is the clearest "
             "endorsement in this table of the reasoning the whole table "
             "rests on.",
        evidence="clinical",
        source="Enjaymo (sutimlimab) label, Warnings and Precautions",
        severity="serious"),
    "MASP2": Physiology(
        signal="MASP-2 is the protease of the lectin pathway, activated by "
               "mannose-binding lectin on microbial carbohydrate.",
        cost="Genuinely disputed, which is why this is the table's one "
             "contested row. MASP-2 deficiency is common -- around 4% of "
             "Caucasians and up to 18% in some African populations -- and "
             "most deficient people are asymptomatic. Early reports linked "
             "the homozygous D120G state to autoimmunity and recurrent "
             "pneumococcal infection; later reports described healthy adults "
             "homozygous for the same variant and cast doubt on the "
             "penetrance, while the disease database still lists "
             "autoimmunity, recurrent severe infection and chronic "
             "inflammation. A deficiency carried by one person in twenty-"
             "five without obvious consequence is the strongest argument "
             "that this pathway is redundant, and it is not the only "
             "published answer.",
        evidence="clinical",
        source="MASP2 deficiency prevalence and penetrance literature; "
               "OMIM 613791",
        severity="contested"),
    "KLKB1": Physiology(
        signal="Prekallikrein starts the contact activation pathway. Like "
               "factor XI downstream of it, it contributes to the "
               "amplification of clotting rather than to the haemostasis "
               "that seals a wound.",
        cost="Very little, which is why the contact pathway attracts drug "
             "developers. Congenital prekallikrein (Fletcher factor) "
             "deficiency typically presents as an isolated prolonged aPTT "
             "found incidentally, with no bleeding tendency. A systematic "
             "review of 53 reported patients found 47.1% asymptomatic and "
             "bleeding episodes in 7.5%, so bleeding is rare rather than "
             "impossible, and the authors could not tie the cardiovascular "
             "comorbidities seen in a third of cases to the deficiency "
             "itself. The same caution that belongs on factor XI belongs "
             "here: a pathway you can remove without bleeding may be a "
             "pathway that was not doing much.",
        evidence="clinical",
        source="Systematic review of prekallikrein deficiency case reports "
               "(53 patients); MedlinePlus Genetics",
        severity="well_tolerated"),

    # --- Cytokines -------------------------------------------------------
    "IL1B": Physiology(
        signal="IL-1 beta is an apex inflammatory cytokine and a principal "
               "driver of the fever and acute-phase response to infection.",
        cost="Fatal infection and sepsis. In CANTOS, deaths attributed to "
             "infection or sepsis were rare but more frequent on canakinumab "
             "than placebo, and this is the reason regulators on both sides "
             "of the Atlantic declined the secondary cardiovascular "
             "prevention indication despite the trial meeting its "
             "cardiovascular endpoint. Opportunistic infection was not seen. "
             "The shape of the decision is the lesson: the drug worked, and "
             "the cost of removing the signal was judged to outweigh it.",
        evidence="clinical",
        source="CANTOS trial safety data; regulatory decisions on "
               "canakinumab for cardiovascular prevention",
        severity="serious"),
    "IL23A": Physiology(
        signal="IL-23 sustains the Th17 response, upstream of the IL-17 "
               "that defends mucosal surfaces against Candida.",
        cost="Less than blocking IL-17 directly, which is the interesting "
             "part. Reviews of the selective p19 inhibitors -- guselkumab, "
             "tildrakizumab, risankizumab -- report no increase in "
             "opportunistic infection, tuberculosis reactivation, oral "
             "candidiasis or inflammatory bowel disease, and their labelled "
             "precautions are limited to infection, hypersensitivity and TB. "
             "Blocking the upstream driver of a pathway is evidently not the "
             "same as blocking its effector, and this table holds both so "
             "the comparison is visible. Latent TB screening is still "
             "advised before starting.",
        evidence="clinical",
        source="Narrative reviews of selective IL-23 p19 inhibitor safety; "
               "class labelling",
        severity="well_tolerated"),
    "IL36G": Physiology(
        signal="IL-36 gamma drives neutrophilic skin inflammation; the "
               "disease it is blocked for, generalised pustular psoriasis, "
               "is caused by unopposed IL-36 signalling.",
        cost="Infection, concentrated in the acute treatment window. In the "
             "Effisayil 1 trial, infections occurred in 6 of 35 patients "
             "(17%) on spesolimab in the first week against 6% on placebo, "
             "and in 24 of 51 (47%) by week 12 among anyone who received it; "
             "antidrug antibodies appeared in 23 of 50 (46%), with systemic "
             "drug reactions. Later and pooled data are more reassuring, "
             "with severe, serious or opportunistic infection low across six "
             "trials. Guidance is to withhold during significant active "
             "infection including tuberculosis.",
        evidence="clinical",
        source="Effisayil 1 and Effisayil 2 trials; pooled spesolimab "
               "safety analysis",
        severity="manageable"),
    "IL33": Physiology(
        signal="IL-33 is an alarmin released by damaged epithelium, "
               "signalling through ST2 to start type 2 inflammation.",
        cost="Not established for IL-33 itself. The reassuring dataset here "
             "is astegolimab, and astegolimab blocks ST2, the receptor -- "
             "across five trials its adverse and serious adverse event rates "
             "matched placebo, with no increase in infection or major "
             "adverse cardiac events and no treatment-related deaths. That "
             "is evidence about the receptor arm of the pathway, and this "
             "row is about the ligand. They are not interchangeable, the "
             "COPD trials are still running, and recording the receptor's "
             "safety record under the ligand's name would be the kind of "
             "borrowing this package keeps catching.",
        evidence="clinical",
        source="Narrative review of five astegolimab (anti-ST2) trials",
        severity="unknown"),
    "TSLP": Physiology(
        signal="TSLP is an epithelial alarmin sitting at the top of the "
               "asthma inflammatory cascade, upstream of both type 2 and "
               "non-type-2 pathways.",
        cost="Mild and mostly unremarkable: arthralgia at 3.8% and "
             "pharyngitis at 4.1% are the commonest labelled reactions, and "
             "bronchitis ran 5% against 6% on placebo in NAVIGATOR. One "
             "signal is worth carrying rather than smoothing away: a "
             "numerical imbalance in serious cardiac adverse events against "
             "placebo appeared in a long-term study, with no causal "
             "relationship established. Live attenuated vaccines are to be "
             "avoided and treatment stopped during serious infection.",
        evidence="clinical",
        source="Tezepelumab EU label; NAVIGATOR trial",
        severity="manageable"),

    # --- Metabolic, musculoskeletal, vascular -----------------------------
    "FGF23": Physiology(
        signal="FGF23 is the hormone that tells the kidney to excrete "
               "phosphate. Blocking it raises serum phosphate, which is the "
               "therapeutic point in X-linked hypophosphataemia and the "
               "hazard everywhere else.",
        cost="Hyperphosphataemia, and ectopic mineralisation if it is not "
             "watched. Labelling requires fasting phosphate monitoring and "
             "keeping the level in the lower part of the normal range "
             "specifically to limit mineral deposition; concurrent oral "
             "phosphate and active vitamin D are contraindicated and must "
             "stop at least a week before starting. Guidance calls for renal "
             "ultrasound for nephrocalcinosis at baseline, six-monthly in "
             "the first year and annually after. Ectopic mineralisation "
             "appeared in animals at supraphysiological phosphate. This is a "
             "managed risk rather than an accepted one, and the management "
             "is the entry.",
        evidence="clinical",
        source="Burosumab (Crysvita) labelling; NICE-based monitoring "
               "guidance; animal toxicology",
        severity="manageable"),
    "GDF15": Physiology(
        signal="GDF15 acts on GFRAL in the area postrema and nucleus "
               "tractus solitarius -- the brainstem nausea and appetite "
               "circuitry -- as a stress signal that suppresses eating.",
        cost="Not established, and the short-term data point the opposite "
             "way from the usual worry. In the phase 2 ponsegromab trial in "
             "cancer cachexia, patients gained weight and activity, and "
             "nausea and vomiting were LOWER on drug than placebo (4% vs 16% "
             "and 5% vs 13%). That is twelve weeks in advanced cancer, which "
             "is not a safety dataset for a signal that normally reports "
             "cellular stress to the brainstem, so the honest grade is "
             "unknown rather than well tolerated.",
        evidence="clinical",
        source="Phase 2 ponsegromab trial (NCT05546476); GDF15-GFRAL "
               "pathway reviews",
        severity="unknown"),
    "MSTN": Physiology(
        signal="Myostatin restrains skeletal muscle growth. Removing the "
               "brake adds muscle.",
        cost="Strikingly little harm, and strikingly little benefit. "
             "Bimagrumab, domagrozumab, LY2495655 and ACE-083 all increased "
             "muscle volume with little or no gain in strength or function; "
             "a meta-analysis in sarcopenia described bimagrumab as safe "
             "with minimal effects on strength, and Pfizer ended the "
             "domagrozumab Duchenne programme after a 166-child trial showed "
             "no effect on muscle mass or strength. Mass and strength do not "
             "track together. This row belongs beside factor XI and the "
             "neurodegeneration antibodies: safe to block, and the blocking "
             "did not deliver.",
        evidence="clinical",
        source="Bimagrumab sarcopenia trials and meta-analysis; "
               "domagrozumab Duchenne phase 2/3",
        severity="well_tolerated"),
    "VWF": Physiology(
        signal="Von Willebrand factor tethers platelets to damaged vessel "
               "wall. Caplacizumab blocks the A1 domain-GPIb interaction, "
               "which is the same interaction that stops bleeding.",
        cost="Bleeding, as the mechanism requires. In HERCULES bleeding-"
             "related adverse events were substantially more frequent than "
             "placebo, mostly epistaxis and gingival bleeding; serious "
             "bleeding occurred in 8 patients (11%) against 1 (1%). "
             "Published all-grade rates disagree across reports of the same "
             "trial -- 45.6% vs 23.3% in one, 65% vs 48% in another -- so "
             "the direction is certain and the magnitude is not, and that "
             "disagreement is recorded rather than averaged away.",
        evidence="clinical",
        source="HERCULES trial reports (figures differ between summaries)",
        severity="serious"),
    "ANGPT2": Physiology(
        signal="Angiopoietin-2 destabilises vessels and sensitises them to "
               "VEGF; blocking both is the design of faricimab.",
        cost="Low in trials and not zero in practice. Pivotal trials "
             "reported intraocular inflammation at rates similar across "
             "arms with no retinal vasculitis, while the label records that "
             "retinal vasculitis and/or retinal vascular occlusion, usually "
             "with intraocular inflammation, have been reported. "
             "Post-marketing case series put faricimab-related inflammation "
             "near 0.62% and describe eyes progressing from mild "
             "inflammation to irreversible vision loss. A rate that is low "
             "across a population and catastrophic in the individual is a "
             "shape no single number conveys.",
        evidence="clinical",
        source="Faricimab pivotal trial reviews; Vabysmo labelling; "
               "post-marketing intraocular inflammation case series",
        severity="manageable"),

    # --- The aggregate targets: well tolerated, and they did not work -----
    "MAPT": Physiology(
        signal="Tau binds and stabilises axonal microtubules. The "
               "therapeutic target is not that pool but the extracellular "
               "seeds thought to propagate pathology between neurons.",
        cost="Nothing demonstrable -- including benefit. Gosuranemab missed "
             "in the phase 2 TANGO trial in early Alzheimer's after an "
             "earlier failure in progressive supranuclear palsy. Tilavonemab "
             "was generally well tolerated in 453 patients and did not "
             "demonstrate efficacy, its authors concluding that further "
             "investigation in early Alzheimer's is not warranted. "
             "Semorinemab slowed one of two primary endpoints and not the "
             "other. These antibodies did not fail on toxicity. They failed "
             "on efficacy, which for this table is the more useful fact: "
             "reaching only the extracellular fraction of an intracellular "
             "protein may mean reaching too little of it to matter.",
        evidence="clinical",
        source="TANGO (gosuranemab); tilavonemab phase 2 in early "
               "Alzheimer's; semorinemab phase 2",
        severity="well_tolerated"),
    "SNCA": Physiology(
        signal="Alpha-synuclein is a presynaptic protein of the neuronal "
               "cytosol. As with tau, the accessible target is the small "
               "extracellular pool, not the aggregates inside the neuron.",
        cost="The same shape as tau, and the same disappointment. "
             "Cinpanemab's SPARK trial was stopped early for lack of "
             "efficacy and Biogen discontinued the programme. Prasinezumab's "
             "PASADENA trial missed its primary endpoint on MDS-UPDRS total "
             "score at week 52 while the sponsors reported a favourable "
             "safety profile and a signal on Part III; a later exploratory "
             "analysis suggested more effect in faster-progressing patients, "
             "with the usual small-subgroup caveats. Safety was not the "
             "problem in either programme.",
        evidence="clinical",
        source="SPARK (cinpanemab, stopped for futility); PASADENA "
               "(prasinezumab, primary endpoint missed)",
        severity="well_tolerated"),

    # --- Fibrosis and immune evasion --------------------------------------
    "TGFB1": Physiology(
        signal="TGF-beta restrains growth and promotes differentiation in "
               "normal epithelium. It is a tumour suppressor there, and an "
               "immune-evasion signal in established tumours -- the same "
               "molecule doing opposite jobs depending on where it is.",
        cost="Skin tumours, among other things. In a phase 1 study of "
             "fresolimumab in 29 patients with melanoma or renal cell "
             "carcinoma, four developed keratoacanthoma-type lesions and one "
             "of those also had a well-differentiated squamous cell "
             "carcinoma; the lesions resolved spontaneously after treatment "
             "ended. Anaemia and gingival, nasal and subconjunctival "
             "bleeding were also reported, and repeated dosing in animals "
             "produced epithelial hyperplasia of gingiva, bladder and nasal "
             "turbinate. An anti-cancer drug that causes cancer is the "
             "sharpest statement in this table of what it means to block a "
             "signal with more than one job. Cardiac lesions are raised as a "
             "concern for chronic TGF-beta blockade in secondary sources and "
             "were NOT verified here, so they are named as unconfirmed "
             "rather than recorded.",
        evidence="clinical",
        source="Fresolimumab phase 1 in melanoma and renal cell carcinoma; "
               "animal repeat-dose findings",
        severity="serious"),
}


def physiology_for(gene: str) -> Physiology | None:
    return PHYSIOLOGY.get((gene or "").upper())


def off_axis(gene: str, verdict: str, why: str) -> LineageCall:
    """The single honest row for a target the lineage question mis-fits.

    Mirrors not_applicable() for viral targets: one row that says the
    question does not apply, rather than eleven rows answering a question
    nobody asked. Returning the per-lineage transcript readings here would
    be worse than silence, because they measure which cells MAKE the
    protein and would render identically to which cells display it.
    """
    if verdict == "intracellular":
        detail = (
            f"{gene} is an intracellular protein. Which healthy cells carry "
            f"it is the wrong question twice over: every cell that expresses "
            f"it holds it inside, where an antibody does not go, so there is "
            f"no surface pool for a binder to engage and no population that "
            f"a binder would deplete. Any brief against it is really a brief "
            f"against the small extracellular fraction -- aggregates or "
            f"seeds -- and how much of the target that fraction represents "
            f"is the question that decides the programme. This page cannot "
            f"answer it. {why}")
    else:
        detail = (
            f"{gene} is secreted, not displayed on a cell. 'Which healthy "
            f"cells carry it' therefore has no answer of the kind this page "
            f"is built to give: the transcript readings for this gene "
            f"measure which cells PRODUCE it, which is a different claim "
            f"and would be read as the same one. Binding it neutralises a "
            f"signal rather than depleting a population, so the risk moves "
            f"to whatever normal physiology depends on that signal -- which "
            f"is not an expression table. "
            + (f"What that signal does, and what is lost without it, is "
               f"recorded separately for {gene}. "
               if physiology_for(gene) else
               f"What {gene} normally does is not recorded here, so that "
               f"half of the question is open rather than answered. ")
            + why)
    # level/evidence are deliberately inert. not_applicable() can say
    # "absent" because no human cell carries a viral protein; here the
    # protein is emphatically present, just not where this axis looks. Any
    # other pair would make this row clearable or countable as a concern in
    # some future caller, and it is neither.
    return LineageCall(
        lineage="off_axis", label="Not a cell-surface target",
        level="unknown", evidence="none",
        source="UniProt subcellular location and feature annotation",
        detail=detail)


HUMAN = "Homo sapiens"


def is_human(target) -> bool:
    """Whether the lineage question applies to this target at all.

    Seven catalogue entries are viral or bacterial -- SARS-CoV-2 spike, HIV
    envelope, influenza haemagglutinin, RSV F, the Ebola glycoprotein,
    Plasmodium circumsporozoite protein, C. difficile toxin B -- and the
    Human Protein Atlas returns nothing for any of them, correctly. Letting
    that empty result flow through as "expression unknown on every lineage"
    would be the familiar mistake in a new place: it states ignorance where
    the truth is that no human cell carries the protein, so there is no
    healthy population to spare. The same distinction the paralog panels
    draw with panel_state "no_human_relatives".
    """
    org = (getattr(target, "organism", "") or "").strip()
    return not org or org.startswith(HUMAN)


def not_applicable(target) -> LineageCall:
    """The single honest row for a non-human target."""
    org = getattr(target, "organism", "") or "an organism other than human"
    return LineageCall(
        lineage="not_applicable", label="Not applicable",
        level="absent", evidence="clinical",
        source="UniProt source organism",
        detail=(f"{target.gene} is a protein of {org}, not of the patient. "
                f"Healthy human cells do not carry it, so there is no "
                f"lineage to spare and cross-reactivity with human proteins "
                f"-- which is a real risk -- is the separate question the "
                f"paralog panel addresses."))


def assess_gene(gene: str, data: ExpressionData | None = None,
                ) -> list[LineageCall]:
    """Everything known about one antigen across every lineage.

    Curated cytometry wins over the transcript screen wherever both exist,
    because they are not the same kind of claim: one measures surface protein
    on sorted primary cells and the other measures mRNA in a bulk population.
    Where only transcript exists, a detected signal becomes a graded concern
    and a non-detection becomes "unknown" -- see the module docstring for the
    GPRC5D case that settles why.
    """
    data = data if data is not None else load_expression()
    gene_u = (gene or "").upper()
    calls: list[LineageCall] = []

    for lin in LINEAGES:
        cur = CURATED.get((gene_u, lin.key))
        rna_level, ntpm, where = _rna_level(gene_u, lin, data)

        if cur is not None:
            detail = cur.note
            if ntpm is not None:
                agrees = (cur.level == rna_level
                          or (cur.level in ("absent", "low")
                              and rna_level == "unknown"))
                detail += (
                    f" Transcript screen: {ntpm:.1f} nTPM in {where}"
                    + ("" if agrees else
                       f", which reads as {rna_level} and does not agree -- "
                       f"the cytometry is the stronger evidence and is what "
                       f"is reported") + ".")
            calls.append(LineageCall(
                lineage=lin.key, label=lin.label, level=cur.level,
                evidence=cur.evidence, source=cur.source, note=cur.note,
                ntpm=ntpm, detail=detail.strip(), harm=cur.harm))
            continue

        if rna_level == "unknown":
            if not (lin.blood_columns or lin.tissue_columns):
                detail = (
                    f"Not assessed. The Human Protein Atlas has no data for "
                    f"this population at all -- its bone marrow single-cell "
                    f"set contains no stem or progenitor cells -- and no "
                    f"cytometry for {gene_u} is recorded here. This is a gap "
                    f"in what is known, not a clean result.")
            elif ntpm is not None:
                detail = (
                    f"Transcript not detected ({ntpm:.1f} nTPM in {where}). "
                    f"Reported as unknown rather than absent: GPRC5D reads "
                    f"0.0 nTPM in tongue and still causes oral toxicity in "
                    f"63% of patients, so a non-detection here does not "
                    f"clear this population. Cytometry would.")
            else:
                detail = ("No transcript data returned for this population, "
                          "and no cytometry recorded.")
            calls.append(LineageCall(
                lineage=lin.key, label=lin.label, level="unknown",
                evidence="none", ntpm=ntpm, detail=detail))
            continue

        calls.append(LineageCall(
            lineage=lin.key, label=lin.label, level=rna_level, evidence="rna",
            source=data.source or "Human Protein Atlas", ntpm=ntpm,
            detail=(f"Transcript screen only: {ntpm:.1f} nTPM in {where}. "
                    f"This says the gene is transcribed in this population, "
                    f"not how much protein reaches the surface, and surface "
                    f"density is what a binder reads. Treat as a flag to "
                    f"check by cytometry.")))

    return calls


def assess_tissues(gene: str, data: ExpressionData | None = None,
                   ) -> list[LineageCall]:
    """The always-on non-haematopoietic screen."""
    data = data if data is not None else load_expression()
    gene_u = (gene or "").upper()
    vals = data.for_gene(gene_u)
    out: list[LineageCall] = []

    column_for = {label: col for col, label in TISSUE_SCREEN}
    curated_sites: set[str] = set()

    for (g, site), cur in CURATED_TISSUE.items():
        if g != gene_u:
            continue
        curated_sites.add(site)
        # Three cases, and the third is the one worth the code: the screen
        # agrees, the screen disagrees, or the screen has no column for this
        # tissue at all. Mesothelium and the nail matrix are the latter, and
        # letting that read like agreement is the failure this whole module
        # is shaped around.
        col = column_for.get(site)
        ntpm = vals.get(col) if col else None
        if ntpm is None:
            screen = (f" The tissue screen has no column for this site, so "
                      f"it is not corroborated or contradicted by the "
                      f"automated half of this page -- it is outside what "
                      f"that screen samples.")
        else:
            rna = band_for_ntpm(ntpm)
            agrees = (cur.level == rna
                      or (cur.level in ("absent", "low") and rna == "unknown"))
            screen = (f" Transcript screen: {ntpm:.1f} nTPM"
                      + ("" if agrees else
                         f", which reads as {rna} and does not agree -- the "
                         f"tissue-level evidence is the stronger of the two "
                         f"and is what is reported") + ".")
        out.append(LineageCall(
            lineage=f"curated_tissue:{site}", label=site,
            level=cur.level, evidence=cur.evidence, source=cur.source,
            note=cur.note, detail=(cur.note + screen).strip(),
            ntpm=ntpm, harm=cur.harm))

    for col, label in TISSUE_SCREEN:
        if col not in vals or label in curated_sites:
            continue
        ntpm = vals[col]
        level = band_for_ntpm(ntpm)
        if level == "unknown":
            continue
        out.append(LineageCall(
            lineage=f"tissue:{label}", label=label, level=level,
            evidence="rna", source=data.source or "Human Protein Atlas",
            ntpm=ntpm,
            detail=(f"{ntpm:.1f} nTPM. A transcript flag in a tissue the "
                    f"therapy is not aimed at; whether protein reaches the "
                    f"surface there is a separate question.")))
    # Curated rows first, strongest evidence first within them, and only then
    # the screen loudest first. Sorting everything by nTPM together would put
    # PSMA's duodenum transcript above PSMA's salivary-gland xerostomia purely
    # because 250.9 is a bigger number than 49.8, which inverts the evidence
    # ordering this module spends most of its length defending.
    out.sort(key=lambda c: (0 if c.lineage.startswith("curated_tissue")
                            else 1,
                            -EVIDENCE_RANK.get(c.evidence, 0),
                            -(c.ntpm or 0.0)))
    return out
