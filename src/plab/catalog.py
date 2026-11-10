"""The disease/target catalogue.

Scope, stated honestly
----------------------
This pipeline needs three things: a protein implicated in a disease, an
experimentally solved structure of it, and ideally something already bound to
that structure so the binding site can be read rather than guessed.

That rules a lot out. Diseases with no single protein driver -- most
psychiatric illness, most metabolic syndrome, polygenic conditions -- have
nothing here to aim at. Targets with no solved structure are equally out of
reach. What remains is still a few thousand proteins spanning cancer,
infection, autoimmunity, cardiovascular and neurodegenerative disease, which
is a large and genuinely useful space, but it is not "any illness".

The catalogue below is a list of CANDIDATES. Nothing here is trusted until
scripts/build_catalog.py has fetched it, matched the structure to the
accession, and derived a binding site. Entries that fail verification are
dropped with a recorded reason and never offered to a user. That is the
difference between a catalogue and a wish list.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
CATALOG_DIR = ROOT / "data" / "catalog"
VERIFIED_PATH = CATALOG_DIR / "verified.json"


@dataclass
class CatalogEntry:
    key: str                 # stable id, used in URLs
    disease: str
    area: str                # disease area, for grouping in the UI
    gene: str
    protein_name: str
    uniprot: str
    pdb_complex: str         # structure containing target + bound partner
    known_binder: str        # what is bound in that structure
    rationale: str           # why this target, in plain language
    clinical_note: str       # the honest caveat about this target
    pdb_apo: str | None = None
    # Which bound partner counts as the binder, where a structure holds more
    # than one. "antibody", or an explicit chain id. Unset means "the largest
    # single interface".
    binder_preference: str | None = None
    # Other proteins a binder should avoid. Filled per target where a
    # meaningful paralog panel exists.
    off_targets: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


# --- Candidates -------------------------------------------------------------
# Every PDB id here is a candidate for verification, not an assertion. The
# build script checks each one actually contains this protein with a partner
# bound to it.

CANDIDATES: list[CatalogEntry] = [
    # ---------------- Blood cancers ----------------
    CatalogEntry(
        key="clec12a_aml", disease="Acute myeloid leukaemia", area="Blood cancer",
        gene="CLEC12A", protein_name="C-type lectin domain family 12 member A (CLL-1)",
        uniprot="Q5QGZ9", pdb_complex="8W9J", pdb_apo="8W8T",
        known_binder="50C1 inhibitory antibody",
        rationale=(
            "Expressed by leukaemic blasts in over 90% of AML patients and, "
            "unlike CD33 and CD123, largely absent from healthy haematopoietic "
            "stem cells. That expression gap is why it is considered the "
            "cleanest myeloid target: engaging it should spare the stem-cell "
            "compartment that regenerates the blood."),
        clinical_note=(
            "CAR and antibody therapies against CD33, CD123 and CLEC12A have "
            "NOT produced durable remissions in most AML patients, and antigen "
            "escape is a documented resistance mechanism. A good target is "
            "necessary and nowhere near sufficient."),
        off_targets=["Q2HXU8", "Q6UXN8", "Q9BXN2", "Q8NC01", "Q12918", "P20138"],
    ),
    CatalogEntry(
        key="cd20_lymphoma", disease="B-cell non-Hodgkin lymphoma", area="Blood cancer",
        gene="MS4A1", protein_name="B-lymphocyte antigen CD20",
        uniprot="P11836", pdb_complex="2OSL",
        known_binder="rituximab Fab",
        rationale=(
            "The target of rituximab, the first monoclonal antibody approved "
            "for cancer. Present on B cells from the pre-B stage through to "
            "memory B cells but absent from stem cells and plasma cells, so "
            "the B-cell compartment can regenerate after treatment."),
        clinical_note=(
            "CD20 is a well-validated target, which makes it a good sanity "
            "check for this pipeline rather than an open problem. Resistance "
            "through antigen loss is well documented."),
    ),
    CatalogEntry(
        key="cd19_ball", disease="B-cell acute lymphoblastic leukaemia",
        area="Blood cancer",
        gene="CD19", protein_name="B-lymphocyte antigen CD19",
        uniprot="P15391", pdb_complex="6AL5",
        known_binder="FMC63 antibody fragment",
        rationale=(
            "The target of the first approved CAR-T therapies. CD19-directed "
            "CAR T cells achieve complete remission in up to 90% of patients "
            "with relapsed B-ALL, making this the clearest proof that "
            "engineered binders against a surface antigen can work."),
        clinical_note=(
            "Durable responses are far rarer in CLL than in B-ALL, and CD19 "
            "loss is a common escape route. The same antigen behaves very "
            "differently in different diseases."),
    ),

    # ---------------- Solid tumours ----------------
    CatalogEntry(
        key="her2_breast", disease="HER2-positive breast cancer", area="Solid tumour",
        gene="ERBB2", protein_name="Receptor tyrosine-protein kinase erbB-2 (HER2)",
        uniprot="P04626", pdb_complex="1N8Z",
        known_binder="trastuzumab Fab",
        rationale=(
            "Amplified in roughly one in five breast cancers and historically "
            "a marker of poor prognosis. Trastuzumab binding to domain IV of "
            "the extracellular region changed that outcome, so the binding "
            "site here is one of the best-validated in oncology."),
        clinical_note=(
            "Cardiotoxicity is a real and dose-limiting effect, because HER2 "
            "is also expressed on cardiomyocytes. A binder that cannot "
            "distinguish tumour from heart tissue is not useful."),
    ),
    CatalogEntry(
        key="egfr_crc", disease="Colorectal and head-and-neck cancer",
        area="Solid tumour",
        gene="EGFR", protein_name="Epidermal growth factor receptor",
        uniprot="P00533", pdb_complex="1YY9",
        known_binder="cetuximab Fab",
        rationale=(
            "A growth-factor receptor driving proliferation in several "
            "epithelial cancers. Cetuximab blocks ligand binding to domain "
            "III, which both prevents activation and flags the cell for "
            "immune clearance."),
        clinical_note=(
            "Tumours carrying downstream KRAS mutations do not respond, "
            "because the signal no longer depends on the receptor. Blocking "
            "the right protein in the wrong patient does nothing."),
    ),
    CatalogEntry(
        key="vegfa_angio", disease="Metastatic cancer and macular degeneration",
        area="Solid tumour",
        gene="VEGFA", protein_name="Vascular endothelial growth factor A",
        uniprot="P15692", pdb_complex="1BJ1",
        known_binder="bevacizumab Fab",
        rationale=(
            "Tumours beyond a few millimetres must recruit a blood supply. "
            "VEGF-A is the principal signal they use, so sequestering it "
            "starves the tumour rather than attacking it directly."),
        clinical_note=(
            "Benefits in solid tumours are usually measured in months, and "
            "blocking a normal physiological signal brings hypertension and "
            "impaired wound healing."),
    ),

    # ---------------- Immuno-oncology ----------------
    CatalogEntry(
        key="pd1_checkpoint", disease="Checkpoint blockade across many cancers",
        area="Immuno-oncology",
        gene="PDCD1", protein_name="Programmed cell death protein 1 (PD-1)",
        uniprot="Q15116", pdb_complex="5GGS",
        known_binder="pembrolizumab Fab",
        rationale=(
            "Tumours exploit an immune brake. Blocking PD-1 releases T cells "
            "that are already present but switched off, which is why it works "
            "across many tumour types rather than one."),
        clinical_note=(
            "Only a minority of patients respond, and releasing an immune "
            "brake causes autoimmune toxicity in a substantial fraction."),
    ),
    CatalogEntry(
        key="pdl1_checkpoint", disease="Checkpoint blockade across many cancers",
        area="Immuno-oncology",
        gene="CD274", protein_name="Programmed cell death 1 ligand 1 (PD-L1)",
        uniprot="Q9NZQ7", pdb_complex="5XXY",
        known_binder="atezolizumab Fab",
        rationale=(
            "The ligand side of the same brake. Targeting the ligand rather "
            "than the receptor leaves one PD-1 interaction intact, which may "
            "shift the toxicity profile."),
        clinical_note=(
            "Expression is heterogeneous and dynamic, so patient selection by "
            "PD-L1 staining is unreliable."),
    ),
    CatalogEntry(
        key="ctla4_checkpoint", disease="Melanoma and other cancers",
        area="Immuno-oncology",
        gene="CTLA4", protein_name="Cytotoxic T-lymphocyte protein 4",
        uniprot="P16410", pdb_complex="5XJ3",
        known_binder="ipilimumab Fab",
        rationale=(
            "A second, earlier immune checkpoint acting during T-cell "
            "priming rather than in the tumour itself."),
        clinical_note=(
            "The most toxic of the approved checkpoint targets; severe "
            "immune-related adverse events are common."),
    ),

    # ---------------- Infectious disease ----------------
    CatalogEntry(
        key="sars2_rbd", disease="COVID-19", area="Infectious disease",
        gene="S", protein_name="SARS-CoV-2 spike receptor-binding domain",
        uniprot="P0DTC2", pdb_complex="6M0J",
        known_binder="human ACE2 receptor",
        rationale=(
            "The receptor-binding domain is how the virus attaches to human "
            "cells through ACE2. A binder covering that interface blocks entry "
            "directly, which is the mechanism of every neutralising antibody "
            "developed for this virus."),
        clinical_note=(
            "This surface mutates rapidly under immune pressure. Most "
            "monoclonals developed early in the pandemic were rendered "
            "useless by later variants -- a binder for a moving target has a "
            "short shelf life."),
    ),
    CatalogEntry(
        key="hiv_gp120", disease="HIV-1 infection", area="Infectious disease",
        gene="env", protein_name="HIV-1 envelope glycoprotein gp120",
        # Strain 93TH057, which is the isolate actually crystallised in 3NGB.
        # This entry first failed verification against P04578 (strain HXB2):
        # the two are 75% identical, the same protein from a different clade,
        # but exact-k-mer matching scored that at 13% and rejected it.
        uniprot="Q0ED31", pdb_complex="3NGB",
        known_binder="VRC01 broadly neutralising antibody",
        rationale=(
            "The CD4-binding site on gp120 is functionally conserved because "
            "the virus needs it to enter cells, so it cannot mutate freely. "
            "VRC01 exploits exactly that constraint."),
        clinical_note=(
            "The site is shielded by a dense glycan coat and sits in a "
            "recessed pocket, which is why broadly neutralising antibodies "
            "against it took decades to find and are hard to elicit. The "
            "structure is one clinical isolate; HIV-1 diversity between "
            "clades is larger than between many separate viral species."),
    ),
    CatalogEntry(
        key="flu_ha", disease="Influenza A", area="Infectious disease",
        gene="HA", protein_name="Influenza A haemagglutinin (H5N1)",
        # 4FQI is an H5N1 haemagglutinin, not the H1N1 one this entry was
        # first paired with. RCSB maps both its chains to Q6DQ33, so that is
        # the sequence the structure should be checked against.
        uniprot="Q6DQ33", pdb_complex="4FQI",
        known_binder="CR9114 broadly neutralising antibody",
        rationale=(
            "The haemagglutinin stem is far more conserved than the head, so "
            "a stem binder can cover many strains at once. This is the basis "
            "of the universal flu vaccine effort, and CR9114 binds that stem "
            "across both influenza A groups."),
        clinical_note=(
            "The stem is poorly immunogenic and partly occluded, which is why "
            "natural infection mostly produces strain-specific head "
            "antibodies instead. Haemagglutinin is also cleaved into HA1 and "
            "HA2, so a designed binder must be checked against the fragment "
            "it actually covers."),
    ),
    CatalogEntry(
        key="rsv_f", disease="Respiratory syncytial virus", area="Infectious disease",
        gene="F", protein_name="RSV fusion glycoprotein F",
        uniprot="P03420", pdb_complex="3IXT",
        known_binder="motavizumab Fab",
        rationale=(
            "RSV is a leading cause of infant hospitalisation. Antigenic site "
            "II on the F protein is conserved between subtypes and is the "
            "target of the prophylactic antibodies used in at-risk infants."),
        clinical_note=(
            "The F protein changes shape dramatically during fusion; a binder "
            "must engage the right conformation to neutralise."),
    ),

    # ---------------- Autoimmune and inflammatory ----------------
    CatalogEntry(
        key="tnf_ra", disease="Rheumatoid arthritis and Crohn's disease",
        area="Autoimmune",
        gene="TNF", protein_name="Tumour necrosis factor alpha",
        uniprot="P01375", pdb_complex="3WD5",
        known_binder="adalimumab Fab",
        rationale=(
            "A central inflammatory cytokine. Neutralising it transformed "
            "treatment of rheumatoid arthritis and inflammatory bowel "
            "disease, making this one of the most commercially and clinically "
            "validated binding sites in medicine."),
        clinical_note=(
            "Suppressing TNF suppresses defence against tuberculosis and "
            "other intracellular pathogens. The toxicity is the mechanism."),
    ),
    CatalogEntry(
        key="il17a_psoriasis", disease="Psoriasis and psoriatic arthritis",
        area="Autoimmune",
        gene="IL17A", protein_name="Interleukin-17A",
        uniprot="Q16552", pdb_complex="5HI3",
        known_binder="CAT-2000 Fab",
        # 5HI3 also contains a synthetic 16-residue macrocyclic peptide
        # inhibitor bound in a different pocket. Both are real binders; the
        # antibody epitope is the one this entry is about.
        binder_preference="antibody",
        rationale=(
            "The effector cytokine of the Th17 axis driving psoriatic skin "
            "inflammation. Blocking it produces some of the highest skin "
            "clearance rates of any biologic."),
        clinical_note=(
            "IL-17 is needed for mucosal antifungal defence, so candidiasis "
            "is a recognised consequence."),
    ),
    CatalogEntry(
        key="il6r_ra", disease="Rheumatoid arthritis", area="Autoimmune",
        gene="IL6R", protein_name="Interleukin-6 receptor subunit alpha",
        # This entry first pointed at 4O9H, which holds interleukin-6 itself
        # rather than its receptor. RCSB reports 7DC8 as an IL-6R/Fab complex.
        uniprot="P08887", pdb_complex="7DC8",
        known_binder="anti-IL-6R Fab",
        rationale=(
            "Blocking the receptor rather than the cytokine shuts down both "
            "classical and trans-signalling by IL-6, a driver of inflammation "
            "in rheumatoid arthritis and of cytokine release syndrome."),
        clinical_note=(
            "IL-6 blockade masks the acute-phase response, so infections can "
            "progress without the usual warning signs. The antibody in this "
            "structure is a research construct rather than tocilizumab, so "
            "the site here is an IL-6R epitope, not the approved drug's."),
    ),
    CatalogEntry(
        key="il6_cytokine", disease="Rheumatoid arthritis and Castleman disease",
        area="Autoimmune",
        gene="IL6", protein_name="Interleukin-6",
        # Found by checking what 4O9H actually contains after this pipeline
        # rejected it as an IL-6 receptor entry. It is a perfectly good
        # target; it was simply filed under the wrong protein.
        uniprot="P05231", pdb_complex="4O9H",
        known_binder="camelid Fab fragment",
        rationale=(
            "The cytokine side of the same axis as IL-6R. Neutralising IL-6 "
            "itself is the mechanism of siltuximab in Castleman disease, and "
            "the comparison with the receptor entry shows how differently a "
            "ligand and its receptor present as design problems."),
        clinical_note=(
            "Same infection risk as receptor blockade. Trials of IL-6 "
            "blockade in diseases beyond its approved indications have "
            "largely disappointed."),
    ),

    # ---------------- Cardiometabolic ----------------
    CatalogEntry(
        key="pcsk9_cholesterol", disease="Familial hypercholesterolaemia",
        area="Cardiometabolic",
        gene="PCSK9", protein_name="Proprotein convertase subtilisin/kexin type 9",
        uniprot="Q8NBP7", pdb_complex="3BPS",
        known_binder="inhibitory Fab",
        rationale=(
            "PCSK9 degrades the LDL receptor. Blocking it leaves more "
            "receptor on the liver surface to clear cholesterol. Humans with "
            "natural loss-of-function variants have low LDL and no apparent "
            "harm, which is about as strong a safety argument as genetics "
            "ever provides for a target."),
        clinical_note=(
            "Efficacy is well established; the barrier has been cost and the "
            "need for injection rather than biology."),
    ),

    # ---------------- Neurodegeneration ----------------
    CatalogEntry(
        key="abeta_alzheimers", disease="Alzheimer's disease", area="Neurology",
        gene="APP", protein_name="Amyloid-beta peptide (from APP)",
        uniprot="P05067", pdb_complex="6CO3",
        known_binder="aducanumab Fab",
        rationale=(
            "Aggregated amyloid-beta is the defining pathology of Alzheimer's "
            "disease, and antibodies that clear plaque demonstrably reduce it "
            "on imaging."),
        clinical_note=(
            "Clearing plaque has produced at best modest cognitive benefit, "
            "and the approvals in this class have been genuinely contested. "
            "This is the clearest example in the catalogue of hitting a "
            "target successfully without reliably helping patients."),
    ),

    # ---------------- Bone and other ----------------
    CatalogEntry(
        key="rankl_osteoporosis", disease="Osteoporosis", area="Bone",
        gene="TNFSF11", protein_name="Receptor activator of NF-kB ligand (RANKL)",
        # First paired with 4GIQ, which is the MOUSE protein bound to mouse
        # RANK. 3URF is the human ligand bound to osteoprotegerin, its own
        # natural decoy receptor -- which is the interaction denosumab was
        # designed to imitate, so it is also the better site to brief on.
        uniprot="O14788", pdb_complex="3URF",
        known_binder="osteoprotegerin (natural decoy receptor)",
        rationale=(
            "RANKL drives osteoclast formation and therefore bone "
            "resorption. Sequestering it shifts the balance back toward bone "
            "formation. The body already does this with osteoprotegerin, so "
            "the interface here has a natural proof of concept."),
        clinical_note=(
            "Stopping treatment causes rapid rebound bone loss, so this is a "
            "commitment rather than a course."),
    ),
    CatalogEntry(
        key="il5_asthma", disease="Severe eosinophilic asthma", area="Respiratory",
        gene="IL5", protein_name="Interleukin-5",
        # 7KQ5, the id first used here, is a bacterial encapsulin and has
        # nothing to do with interleukin-5. 9GVN was found by asking RCSB
        # which entries contain P05113 with something bound to them.
        uniprot="P05113", pdb_complex="9GVN",
        known_binder="depemokimab Fab",
        rationale=(
            "IL-5 controls eosinophil maturation and survival. In the "
            "eosinophilic asthma subtype, removing it reduces exacerbations "
            "substantially."),
        clinical_note=(
            "Only helps the eosinophil-high subtype; in other asthma "
            "phenotypes it does nothing, so the biomarker matters as much as "
            "the drug."),
    ),

    # ---------------- Added in the second catalogue pass ----------------
    # Every structure below was chosen by asking RCSB which entries contain
    # the accession with a partner bound, not from memory. The three wrong
    # PDB ids this catalogue started with are the reason for that rule.

    CatalogEntry(
        key="bcma_myeloma", disease="Multiple myeloma", area="Blood cancer",
        gene="TNFRSF17", protein_name="B-cell maturation antigen (BCMA)",
        uniprot="Q02223", pdb_complex="4ZFO",
        known_binder="J22.9-xi Fab",
        rationale=(
            "Restricted almost entirely to plasma cells, which is why it "
            "became the target of the first approved myeloma CAR-T and "
            "T-cell engagers. The extracellular domain is tiny -- about 50 "
            "residues -- so almost all of it is binding surface."),
        clinical_note=(
            "BCMA is shed from the cell surface, and soluble BCMA acts as a "
            "decoy that soaks up the drug. Relapse after an initial response "
            "is the norm rather than the exception."),
    ),
    CatalogEntry(
        key="cd38_myeloma", disease="Multiple myeloma", area="Blood cancer",
        gene="CD38", protein_name="ADP-ribosyl cyclase/cyclic ADP-ribose hydrolase 1",
        uniprot="P28907", pdb_complex="4CMH",
        known_binder="isatuximab (SAR650984) Fab",
        rationale=(
            "Highly expressed on malignant plasma cells. Daratumumab and "
            "isatuximab against this antigen changed the standard of care in "
            "myeloma, and this structure is the isatuximab epitope itself."),
        clinical_note=(
            "CD38 is also on red cells, so the drugs interfere with blood "
            "cross-matching, and on NK cells, which the treatment depletes."),
    ),
    CatalogEntry(
        key="cd33_aml", disease="Acute myeloid leukaemia", area="Blood cancer",
        gene="CD33", protein_name="Myeloid cell surface antigen CD33 (Siglec-3)",
        uniprot="P20138", pdb_complex="9VL2",
        known_binder="10C8 Fab",
        rationale=(
            "The other myeloid target, and the instructive contrast with "
            "CLEC12A: gemtuzumab ozogamicin against CD33 is approved, but "
            "CD33 is also on healthy haematopoietic progenitors, so hitting "
            "it costs the stem-cell compartment."),
        clinical_note=(
            "Veno-occlusive liver disease and prolonged cytopenias follow "
            "from that lack of selectivity. This entry exists so a designed "
            "CLEC12A binder can be scored against the surface it must NOT "
            "engage."),
    ),
    CatalogEntry(
        key="cd22_ball", disease="B-cell acute lymphoblastic leukaemia",
        area="Blood cancer",
        gene="CD22", protein_name="B-cell receptor CD22",
        uniprot="P20273", pdb_complex="5VL3",
        known_binder="epratuzumab Fab",
        rationale=(
            "A second B-cell antigen, which matters because CD19 loss is the "
            "main escape route from CD19 CAR-T. Inotuzumab ozogamicin "
            "against CD22 works in exactly those relapsed patients."),
        clinical_note=(
            "CD22 is rapidly internalised, which suits a drug conjugate but "
            "works against a binder that needs to stay on the surface."),
    ),

    CatalogEntry(
        key="trop2_tnbc", disease="Triple-negative breast and urothelial cancer",
        area="Solid tumour",
        gene="TACSTD2", protein_name="Trophoblast cell-surface antigen 2 (TROP2)",
        uniprot="P09758", pdb_complex="9PI9",
        known_binder="sacituzumab Fab",
        rationale=(
            "Over-expressed across many epithelial cancers including "
            "triple-negative breast cancer, where few targeted options "
            "exist. This is the sacituzumab govitecan epitope."),
        clinical_note=(
            "TROP2 is present on normal epithelium too, so neutropenia and "
            "diarrhoea are dose-limiting. The therapeutic window comes from "
            "expression level, not from presence or absence."),
    ),
    CatalogEntry(
        key="nectin4_urothelial", disease="Urothelial carcinoma",
        area="Solid tumour",
        gene="NECTIN4", protein_name="Nectin-4",
        uniprot="Q96NY8", pdb_complex="9KKJ",
        known_binder="9MW2821 Fab",
        rationale=(
            "An adhesion molecule re-expressed in bladder cancer and largely "
            "absent from adult normal tissue, which is a rare and valuable "
            "property. Enfortumab vedotin against it is approved."),
        clinical_note=(
            "Skin toxicity and hyperglycaemia are common, and the residual "
            "normal expression in skin is thought to explain the former."),
    ),
    CatalogEntry(
        key="msln_mesothelioma", disease="Mesothelioma and pancreatic cancer",
        area="Solid tumour",
        gene="MSLN", protein_name="Mesothelin (C-terminal fragment)",
        uniprot="Q13421", pdb_complex="9P4C",
        known_binder="RO4 Fab",
        rationale=(
            "Mesothelin is almost absent from vital tissue and strongly "
            "expressed in mesothelioma and pancreatic adenocarcinoma, which "
            "makes it one of the most attractive solid-tumour antigens on "
            "paper."),
        clinical_note=(
            "Repeated clinical failure despite that profile. The target here "
            "is a 17-residue C-terminal peptide, so the binding surface is "
            "very small -- a genuinely hard design problem."),
    ),

    CatalogEntry(
        key="lag3_checkpoint", disease="Melanoma", area="Immuno-oncology",
        gene="LAG3", protein_name="Lymphocyte activation gene 3 protein",
        uniprot="P18627", pdb_complex="7TZH",
        known_binder="scFv F7",
        rationale=(
            "The third checkpoint to reach approval, as relatlimab combined "
            "with anti-PD-1 in melanoma. It works where PD-1 blockade alone "
            "does not, which is the argument for a third brake."),
        clinical_note=(
            "The benefit over PD-1 blockade alone is real but modest, and "
            "combining checkpoints compounds autoimmune toxicity."),
    ),
    CatalogEntry(
        key="tigit_checkpoint", disease="Lung cancer", area="Immuno-oncology",
        gene="TIGIT", protein_name="T-cell immunoreceptor with Ig and ITIM domains",
        uniprot="Q495A1", pdb_complex="8VTD",
        known_binder="vibostolimab Fab",
        rationale=(
            "A checkpoint that looked as promising as PD-1 preclinically, "
            "with a clean 1.2 A structure of a clinical antibody bound to "
            "it, which makes it an unusually well-defined design target."),
        clinical_note=(
            "Phase III results have been largely negative. This is the "
            "catalogue's clearest live example of a well-validated binding "
            "site that has not translated into patient benefit."),
    ),
    CatalogEntry(
        key="cd47_macrophage", disease="Myeloid malignancies and lymphoma",
        area="Immuno-oncology",
        gene="CD47", protein_name="Leukocyte surface antigen CD47",
        uniprot="Q08722", pdb_complex="5TZU",
        known_binder="B6H12.2 Fab",
        rationale=(
            "CD47 is the 'do not eat me' signal tumours use to evade "
            "macrophages. Blocking it engages innate rather than adaptive "
            "immunity, a different mechanism from the T-cell checkpoints."),
        clinical_note=(
            "CD47 is on every red blood cell, so anaemia is intrinsic to the "
            "mechanism, and the lead programme in this class was halted. "
            "Selectivity is the entire problem here."),
    ),

    CatalogEntry(
        key="csp_malaria", disease="Malaria", area="Infectious disease",
        gene="CSP", protein_name="Plasmodium falciparum circumsporozoite protein",
        uniprot="P19597", pdb_complex="8F9U",
        known_binder="Ky15.7 Fab",
        rationale=(
            "Malaria kills several hundred thousand people a year, most of "
            "them children. The circumsporozoite protein is the target of "
            "both approved vaccines and of monoclonal antibodies that "
            "prevent infection outright in trials."),
        clinical_note=(
            "Protection from the vaccines is partial and wanes within "
            "months. The epitope here is a short junctional peptide, so the "
            "available binding surface is small."),
    ),
    CatalogEntry(
        key="ebov_gp", disease="Ebola virus disease", area="Infectious disease",
        gene="GP", protein_name="Zaire ebolavirus glycoprotein GP2 (stem)",
        uniprot="Q05320", pdb_complex="9TGC",
        known_binder="B10 anti-stem antibody",
        rationale=(
            "Antibodies against the Ebola glycoprotein are the first "
            "treatment ever shown to reduce mortality in an acute viral "
            "haemorrhagic fever. The GP2 stem is the conserved part."),
        clinical_note=(
            "Approved antibodies cover Zaire ebolavirus and do not "
            "neutralise Sudan or Bundibugyo species, so an outbreak of the "
            "wrong species has no treatment."),
    ),
    CatalogEntry(
        key="tcdb_cdiff", disease="Clostridioides difficile infection",
        area="Infectious disease",
        gene="tcdB", protein_name="C. difficile toxin B (CROP domain)",
        uniprot="P18177", pdb_complex="4NP4",
        known_binder="bezlotoxumab Fab",
        rationale=(
            "The disease is caused by a toxin, not by the bacterium "
            "invading, so neutralising the toxin treats the illness without "
            "further disturbing the gut flora that antibiotics already "
            "damaged. Bezlotoxumab does exactly this."),
        clinical_note=(
            "It reduces recurrence rather than curing the acute episode, and "
            "showed excess heart failure in patients who already had it."),
    ),

    CatalogEntry(
        key="il23_psoriasis", disease="Psoriasis and Crohn's disease",
        area="Autoimmune",
        gene="IL23A", protein_name="Interleukin-23 subunit p19",
        uniprot="Q9NPF7", pdb_complex="3D85",
        known_binder="7G10 Fab",
        # 3D85 holds p19 bound BOTH by a neutralising Fab and by its own p40
        # partner subunit. Without this hint the larger of the two interfaces
        # wins and the brief describes the wrong site.
        binder_preference="antibody",
        rationale=(
            "Blocking the p19 subunit hits IL-23 alone, where blocking p40 "
            "hits IL-12 as well. Risankizumab and guselkumab exploit that "
            "distinction and clear psoriatic skin better than the older "
            "p40 antibodies."),
        clinical_note=(
            "A good illustration that subunit choice is the whole design "
            "question for a shared-subunit cytokine."),
    ),
    CatalogEntry(
        key="il4r_eczema", disease="Atopic dermatitis and asthma",
        area="Autoimmune",
        gene="IL4R", protein_name="Interleukin-4 receptor subunit alpha",
        uniprot="P24394", pdb_complex="6WGL",
        known_binder="dupilumab Fab",
        rationale=(
            "IL-4Ra is shared by IL-4 and IL-13, so one binder shuts down "
            "both arms of type-2 inflammation. Dupilumab against this exact "
            "epitope works across eczema, asthma and nasal polyposis."),
        clinical_note=(
            "Conjunctivitis is common and unexplained, and a transient rise "
            "in blood eosinophils follows from the mechanism."),
    ),
    CatalogEntry(
        key="baff_lupus", disease="Systemic lupus erythematosus",
        area="Autoimmune",
        gene="TNFSF13B", protein_name="B-cell activating factor (BAFF/BLyS)",
        uniprot="Q9Y275", pdb_complex="5Y9J",
        known_binder="belimumab Fab",
        rationale=(
            "Belimumab was the first new lupus drug in over fifty years, and "
            "this is its epitope. BAFF controls B-cell survival, so "
            "sequestering it removes the autoreactive cells gradually."),
        clinical_note=(
            "The effect size is small and takes months to appear. Lupus is a "
            "reminder that a validated target can still make a marginal drug."),
    ),
    CatalogEntry(
        key="il1b_inflammation",
        disease="Periodic fever syndromes and atherosclerosis",
        area="Autoimmune",
        gene="IL1B", protein_name="Interleukin-1 beta",
        uniprot="P01584", pdb_complex="5BVP",
        known_binder="canakinumab Fab",
        rationale=(
            "Canakinumab against IL-1b is curative in the hereditary "
            "periodic fever syndromes, and in CANTOS it reduced "
            "cardiovascular events -- the first direct evidence that "
            "inflammation itself drives atherosclerosis."),
        clinical_note=(
            "CANTOS also showed increased fatal infection. It was never "
            "approved for cardiovascular prevention, on cost and risk."),
    ),

    CatalogEntry(
        key="c5_pnh", disease="Paroxysmal nocturnal haemoglobinuria",
        area="Rare disease",
        gene="C5", protein_name="Complement C5 (beta chain)",
        uniprot="P01031", pdb_complex="5B71",
        known_binder="crovalimab (SKY59) Fab",
        rationale=(
            "Blocking complement C5 stops the terminal lytic pathway that "
            "destroys red cells in PNH. Eculizumab against this protein "
            "turned a fatal disease into a managed one, which is about as "
            "complete a success as protein therapeutics offer."),
        clinical_note=(
            "Patients become vulnerable to meningococcal infection for life "
            "and need vaccination and standby antibiotics. The cost is also "
            "among the highest of any medicine."),
    ),
    CatalogEntry(
        key="tslp_asthma", disease="Severe asthma", area="Respiratory",
        gene="TSLP", protein_name="Thymic stromal lymphopoietin",
        uniprot="Q969D9", pdb_complex="5J13",
        known_binder="tezepelumab Fab",
        rationale=(
            "TSLP sits upstream of the whole type-2 cascade, so blocking it "
            "helps severe asthma regardless of eosinophil count -- unlike "
            "IL-5 blockade, which needs the right subtype."),
        clinical_note=(
            "Acting upstream of several pathways means broader "
            "immunosuppression, and the long-term consequences are not yet "
            "known for this one."),
    ),
    CatalogEntry(
        key="ige_allergy", disease="Allergic asthma and chronic urticaria",
        area="Allergy",
        gene="IGHE", protein_name="Immunoglobulin E Fc region",
        uniprot="P01854", pdb_complex="4J4P",
        known_binder="anti-IgE Fab",
        rationale=(
            "IgE is the antibody that triggers allergic reactions. "
            "Omalizumab prevents it docking onto mast cells, and the same "
            "mechanism now underpins food-allergy treatment."),
        clinical_note=(
            "The binder must engage free IgE and not IgE already bound to "
            "mast cells, or it would trigger the very reaction it prevents. "
            "That constraint is a selectivity problem, not an affinity one."),
    ),
    CatalogEntry(
        key="synuclein_parkinsons", disease="Parkinson's disease",
        area="Neurology",
        gene="SNCA", protein_name="Alpha-synuclein (peptide)",
        uniprot="P37840", pdb_complex="8JJV",
        known_binder="nanobody",
        rationale=(
            "Aggregated alpha-synuclein is the defining pathology of "
            "Parkinson's disease, and clearing it is the most direct "
            "disease-modifying strategy anyone has for the condition."),
        clinical_note=(
            "Antibodies against alpha-synuclein have so far failed to slow "
            "progression in trials. Like amyloid-beta, this is a target that "
            "can be hit without the patient benefiting."),
    ),
    CatalogEntry(
        key="folr1_ovarian", disease="Ovarian cancer", area="Solid tumour",
        gene="FOLR1", protein_name="Folate receptor alpha",
        # Deliberately kept although it cannot pass: every structure of this
        # protein in the PDB is unbound, so there is no interface to read.
        # It is the clearest example of the method's actual limit -- a real,
        # drugged, well-structured target that this pipeline cannot brief.
        uniprot="P15328", pdb_complex="4KM6",
        known_binder="none in the PDB",
        rationale=(
            "Mirvetuximab soravtansine against folate receptor alpha is "
            "approved for platinum-resistant ovarian cancer, so the target "
            "is validated in patients."),
        clinical_note=(
            "Included to show a limit rather than a result. The protein is "
            "solved at 1.55 A but only on its own, and this pipeline reads "
            "binding sites off a bound partner. No partner, no brief."),
    ),

    # ================= Second research pass =================
    # Every accession below came from scripts/resolve_gene.py and every PDB id
    # from scripts/rcsb_lookup.py pick. Nothing here was typed from memory.
    # The resolver earned its keep immediately: the symbol SELP returns
    # selenoprotein P before P-selectin, and the symbol MET returns a
    # transcription modulator alongside the kinase.

    # ---------------- Haematology ----------------
    CatalogEntry(
        key="tfpi_haemophilia", disease="Haemophilia A and B", area="Haematology",
        gene="TFPI", protein_name="Tissue factor pathway inhibitor",
        uniprot="P10646", pdb_complex="5NMV",
        known_binder="2F22 Fab",
        rationale=(
            "Haemophilia is missing a clotting factor, and the obvious fix is "
            "to replace it. Blocking TFPI instead takes the brake off the "
            "clotting cascade that remains, which works whichever factor is "
            "missing and works in patients who have developed antibodies "
            "against replacement factor. Concizumab does this in the clinic."),
        clinical_note=(
            "Taking a brake off coagulation is as dangerous as it sounds. The "
            "concizumab programme was paused for non-fatal thrombotic events "
            "and resumed at a lower dose with added monitoring. The "
            "therapeutic window here is narrow and is the whole problem."),
    ),
    CatalogEntry(
        key="selp_sicklecell", disease="Sickle cell disease", area="Haematology",
        gene="SELP", protein_name="P-selectin",
        uniprot="P16109", pdb_complex="1G1S",
        known_binder="PSGL-1 peptide",
        rationale=(
            "Vaso-occlusive crisis is cells sticking to vessel walls, and "
            "P-selectin is the adhesion molecule that starts it. Crizanlizumab "
            "blocks it to reduce crisis frequency. The structure here is "
            "P-selectin holding its natural ligand PSGL-1, so the binding site "
            "is the one a blocker has to cover."),
        clinical_note=(
            "Crizanlizumab's confirmatory trial missed its primary endpoint "
            "and the EU authorisation was withdrawn in 2023. Reducing adhesion "
            "measurably is not the same as reducing crises, which is the "
            "distinction this whole project is about."),
    ),

    # ---------------- Cardiometabolic ----------------
    CatalogEntry(
        key="f11_thrombosis", disease="Thrombosis", area="Cardiometabolic",
        gene="F11", protein_name="Coagulation factor XI",
        uniprot="P03951", pdb_complex="6AOD",
        known_binder="anti-FXIa Fab",
        rationale=(
            "Every anticoagulant in wide use trades clot prevention against "
            "bleeding. Factor XI is the exception worth chasing: people born "
            "deficient in it are protected from thrombosis but bleed little, "
            "which suggests the two can be separated. Abelacimab and "
            "milvexian are testing exactly that."),
        clinical_note=(
            "The human genetics are suggestive, not conclusive, and bleeding "
            "less is not the same as not bleeding. Phase 3 outcome data are "
            "what decides this, and they are not all in."),
    ),

    # ---------------- Metabolic ----------------
    CatalogEntry(
        key="glp1r_diabetes", disease="Type 2 diabetes and obesity",
        area="Metabolic",
        gene="GLP1R", protein_name="Glucagon-like peptide 1 receptor",
        uniprot="P43220", pdb_complex="5E94",
        known_binder="anti-GLP-1R Fab",
        rationale=(
            "The receptor behind semaglutide and tirzepatide, and the most "
            "consequential metabolic target of the decade. The structure here "
            "is the extracellular domain with a Fab bound, which is the part "
            "a protein binder could reach."),
        clinical_note=(
            "The approved drugs are peptide AGONISTS that must switch this "
            "receptor on. A designed binder that merely sticks to the "
            "extracellular domain would at best block it, which is the "
            "opposite effect. Binding and agonism are different problems and "
            "this pipeline scores only the first."),
    ),
    CatalogEntry(
        key="gcgr_diabetes", disease="Type 2 diabetes", area="Metabolic",
        gene="GCGR", protein_name="Glucagon receptor",
        uniprot="P47871", pdb_complex="4LF3",
        known_binder="mAb1 allosteric Fab",
        binder_preference="antibody",
        rationale=(
            "Glucagon raises blood glucose; blocking its receptor lowers it, "
            "and does so by a mechanism independent of insulin. Antagonist "
            "antibodies have reached the clinic, and dual GLP-1/glucagon "
            "agonists are now a major class. The structure here is an "
            "allosteric antagonist Fab bound to the receptor's extracellular "
            "face -- a surface a protein binder could actually reach."),
        clinical_note=(
            "Glucagon-receptor blockade raises liver enzymes and LDL "
            "cholesterol and causes alpha-cell hyperplasia. Several "
            "programmes stopped for this. A target can be real and still be "
            "the wrong thing to block."),
    ),
    CatalogEntry(
        key="lepr_lipodystrophy", disease="Leptin receptor deficiency",
        area="Metabolic",
        gene="LEPR", protein_name="Leptin receptor",
        uniprot="P48357", pdb_complex="3V6O",
        known_binder="9F8 Fab",
        rationale=(
            "Congenital leptin-receptor deficiency causes severe childhood "
            "obesity that no amount of diet addresses, and an agonist at this "
            "receptor is one of the few genuinely causal treatments in "
            "obesity medicine. A rare disease where the mechanism is known "
            "exactly."),
        clinical_note=(
            "An agonist has to assemble a signalling complex, not just stick "
            "to the receptor. Mibavademab is in trials; most leptin-pathway "
            "drugs do nothing for common obesity, where leptin is already "
            "high and the problem is resistance to it."),
    ),
    CatalogEntry(
        key="klb_mash", disease="Metabolic dysfunction-associated steatohepatitis",
        area="Metabolic",
        gene="KLB", protein_name="Beta-klotho",
        uniprot="Q86Z14", pdb_complex="5WI9",
        known_binder="39F7 agonist Fab",
        rationale=(
            "Beta-klotho is the co-receptor that makes FGF21 signalling "
            "tissue-specific, and agonist antibodies against it reproduce the "
            "metabolic benefits of FGF21 with a far longer half-life. The "
            "structure is one of those agonist Fabs bound to it."),
        clinical_note=(
            "MASH trials have repeatedly shown improved liver fat and "
            "biomarkers without the fibrosis improvement that regulators "
            "require. Moving a marker is not treating the disease."),
    ),
    CatalogEntry(
        key="gdf15_cachexia", disease="Cancer cachexia", area="Metabolic",
        gene="GDF15", protein_name="Growth/differentiation factor 15",
        uniprot="Q99988", pdb_complex="5VZ4",
        known_binder="GFRAL receptor ectodomain",
        rationale=(
            "Cachexia, the wasting that kills many cancer patients before the "
            "tumour does, is driven in part by GDF15 acting on GFRAL in the "
            "brainstem. Ponsegromab blocks it and produced weight gain in "
            "trials. The structure here is GDF15 held by its own receptor, so "
            "the interface is the one that matters."),
        clinical_note=(
            "Weight gained is not survival or function gained, and cachexia "
            "has multiple drivers. GDF15 also rises in pregnancy and "
            "hyperemesis, so the biology is not solely pathological."),
    ),

    # ---------------- Ophthalmology ----------------
    CatalogEntry(
        key="angpt2_dme", disease="Diabetic macular oedema and wet AMD",
        area="Ophthalmology",
        gene="ANGPT2", protein_name="Angiopoietin-2",
        uniprot="O15123", pdb_complex="4ZFG",
        known_binder="5A12 Fab",
        rationale=(
            "Anti-VEGF injections preserve sight but need repeating for life. "
            "Faricimab blocks angiopoietin-2 as well as VEGF-A, stabilising "
            "vessels by a second mechanism and stretching the interval "
            "between injections. The structure is an anti-ANG2 Fab bound to "
            "the fibrinogen-like domain."),
        clinical_note=(
            "The gain is dosing interval, not better vision: faricimab was "
            "non-inferior to aflibercept on acuity. A real benefit to "
            "patients, and not the cure the mechanism might suggest."),
    ),
    CatalogEntry(
        key="c1q_geographic", disease="Geographic atrophy", area="Ophthalmology",
        gene="C1QB", protein_name="Complement C1q subcomponent subunit B",
        uniprot="P02746", pdb_complex="6Z6V",
        known_binder="nanobody C1qNb75",
        binder_preference="antibody",
        rationale=(
            "C1q starts the classical complement pathway and also tags "
            "synapses for pruning, which is why it is pursued in both retinal "
            "degeneration and neurodegeneration. Blocking it is meant to stop "
            "the photoreceptor loss that complement drives in geographic "
            "atrophy. C1q is a six-armed assembly of three different "
            "subunits; this entry is subunit B because that is the one the "
            "nanobody in this structure actually contacts."),
        clinical_note=(
            "The two approved geographic-atrophy drugs slow lesion growth "
            "without demonstrating functional vision benefit at one and two "
            "years. Anti-C1q antibodies remain investigational. Slowing a "
            "lesion on an image is not preserving sight."),
    ),
    CatalogEntry(
        key="igf1r_ted", disease="Thyroid eye disease", area="Ophthalmology",
        gene="IGF1R", protein_name="Insulin-like growth factor 1 receptor",
        uniprot="P08069", pdb_complex="5U8R",
        known_binder="Fv 24-60",
        rationale=(
            "Teprotumumab against IGF1R was the first drug to reverse "
            "proptosis in thyroid eye disease rather than merely suppress "
            "inflammation, an unusually clean example of blocking a receptor "
            "and reversing an anatomical change."),
        clinical_note=(
            "Permanent hearing loss occurs in a meaningful minority and is "
            "the reason the label carries a warning. The receptor is "
            "expressed far beyond the orbit, which is exactly the "
            "cross-reactivity problem this software tries to measure."),
    ),

    # ---------------- Rare disease / complement ----------------
    CatalogEntry(
        key="c3_complement", disease="Paroxysmal nocturnal haemoglobinuria",
        area="Rare disease",
        gene="C3", protein_name="Complement C3",
        uniprot="P01024", pdb_complex="2WY7",
        known_binder="Staphylococcus aureus Sbi-IV domain",
        rationale=(
            "C3 sits upstream of C5, so blocking it catches the extravascular "
            "haemolysis that C5 inhibitors leave behind; pegcetacoplan works "
            "this way. The bound partner here is not a drug but a bacterial "
            "immune-evasion protein that evolved to grip C3d, which makes it "
            "a template for what a designed binder has to do."),
        clinical_note=(
            "Blocking complement at C3 removes more of the defence against "
            "encapsulated bacteria than blocking at C5. Vaccination and "
            "prophylaxis are mandatory, and breakthrough infection is the "
            "known risk of going further upstream."),
    ),
    CatalogEntry(
        key="cfd_pnh", disease="Paroxysmal nocturnal haemoglobinuria",
        area="Rare disease",
        gene="CFD", protein_name="Complement factor D",
        uniprot="P00746", pdb_complex="4D9R",
        known_binder="anti-factor D Fab",
        rationale=(
            "Factor D is the rate-limiting enzyme of the alternative "
            "complement pathway and circulates at very low concentration, so "
            "a little inhibition goes a long way. Danicopan and iptacopan "
            "target this pathway orally; the structure here is an inhibitory "
            "Fab."),
        clinical_note=(
            "Lampalizumab, an anti-factor D Fab, failed two large phase 3 "
            "trials in geographic atrophy despite good target engagement. "
            "Engaging the right protein is not the same as changing the "
            "disease."),
    ),

    # ---------------- Gastrointestinal ----------------
    CatalogEntry(
        key="itgb7_ibd", disease="Ulcerative colitis and Crohn disease",
        area="Gastrointestinal",
        gene="ITGB7", protein_name="Integrin beta-7",
        uniprot="P26010", pdb_complex="3V4V",
        known_binder="Act-1 Fab",
        binder_preference="antibody",
        rationale=(
            "Vedolizumab blocks the alpha-4-beta-7 integrin that homes "
            "lymphocytes specifically to the gut, so the immunosuppression "
            "stays where the disease is. The structure is the alpha-4-beta-7 "
            "headpiece with Act-1, the parent antibody of vedolizumab."),
        clinical_note=(
            "Gut selectivity is why vedolizumab avoids the brain infection "
            "that sank natalizumab, which blocks alpha-4 more broadly. It is "
            "also why it works slowly and not for everyone."),
    ),
    CatalogEntry(
        key="madcam1_ibd", disease="Inflammatory bowel disease",
        area="Gastrointestinal",
        gene="MADCAM1", protein_name="Mucosal addressin cell adhesion molecule 1",
        uniprot="Q13477", pdb_complex="4HCR",
        known_binder="PF-547659 Fab",
        rationale=(
            "The other half of the gut-homing interaction: MAdCAM-1 is the "
            "address on gut blood vessels that alpha-4-beta-7 reads. Blocking "
            "the address rather than the receptor is an alternative route to "
            "the same selectivity."),
        clinical_note=(
            "Ontamalimab reached phase 3 and was discontinued for portfolio "
            "reasons rather than clear failure, which leaves the approach "
            "unresolved rather than disproved."),
    ),

    # ---------------- Dermatology ----------------
    CatalogEntry(
        key="il13_dermatitis", disease="Atopic dermatitis and asthma",
        area="Dermatology",
        gene="IL13", protein_name="Interleukin-13",
        uniprot="P35225", pdb_complex="3L5X",
        known_binder="H2L6 Fab",
        rationale=(
            "Tralokinumab and lebrikizumab target IL-13 directly rather than "
            "the shared IL-4 receptor, which narrows the blockade to one "
            "cytokine. A clean test of whether narrower is better."),
        clinical_note=(
            "Dupilumab, which blocks the receptor and so both IL-4 and "
            "IL-13, generally performs better. Narrower is not automatically "
            "safer or more effective."),
    ),
    CatalogEntry(
        key="il31ra_pruritus", disease="Prurigo nodularis and atopic dermatitis",
        area="Dermatology",
        gene="IL31RA", protein_name="Interleukin-31 receptor subunit alpha",
        uniprot="Q8NI17", pdb_complex="9EZA",
        known_binder="nemolizumab scFv",
        rationale=(
            "IL-31 is the itch cytokine, signalling directly to sensory "
            "neurons. Nemolizumab blocks its receptor and relieves itch "
            "faster than drugs aimed at the inflammation beneath it, which is "
            "why it matters for prurigo nodularis."),
        clinical_note=(
            "Relieving itch is not clearing disease, and nemolizumab carries "
            "a signal of new-onset or worsening asthma in atopic dermatitis "
            "trials."),
    ),
    CatalogEntry(
        key="il36g_gpp", disease="Generalised pustular psoriasis",
        area="Dermatology",
        gene="IL36G", protein_name="Interleukin-36 gamma",
        uniprot="Q9NZH8", pdb_complex="9I7X",
        known_binder="IL-36 receptor and accessory protein",
        rationale=(
            "Generalised pustular psoriasis is driven by unopposed IL-36 "
            "signalling, often from loss-of-function mutations in its natural "
            "antagonist. Spesolimab blocks the receptor and clears flares "
            "within a week, one of the clearest genotype-to-drug stories in "
            "dermatology."),
        clinical_note=(
            "This is a rare, severe disease and the evidence base is small. "
            "The structure here is the cytokine held by its receptor, so the "
            "interface is the signalling one rather than a drug's."),
    ),
    CatalogEntry(
        key="ox40_dermatitis", disease="Atopic dermatitis", area="Dermatology",
        gene="TNFRSF4", protein_name="OX40 (TNF receptor superfamily member 4)",
        uniprot="P43489", pdb_complex="6OKM",
        known_binder="Fab 3C8",
        rationale=(
            "OX40 is a costimulatory receptor on activated T cells. Depleting "
            "or blocking OX40-positive cells aims at the pathogenic T cells "
            "themselves rather than the cytokines they release, which is why "
            "rocatinlimab and amlitelimab are being tested for durable "
            "remission off drug."),
        clinical_note=(
            "Depleting activated T cells carries infection and malignancy "
            "questions that cytokine blockade does not, and the durability "
            "claim is not yet established."),
    ),

    # ---------------- Respiratory ----------------
    CatalogEntry(
        key="il33_copd", disease="COPD and severe asthma", area="Respiratory",
        gene="IL33", protein_name="Interleukin-33",
        uniprot="O95760", pdb_complex="9WWH",
        known_binder="tozorakimab Fab",
        rationale=(
            "IL-33 is an alarmin released by damaged epithelium, upstream of "
            "the type-2 cytokines that existing asthma biologics block. COPD "
            "has almost no targeted therapy, and IL-33 is the most advanced "
            "attempt: tozorakimab and itepekimab are both in late trials."),
        clinical_note=(
            "Itepekimab worked in former smokers and not in current smokers "
            "across its two phase 3 trials. COPD has defeated many "
            "mechanistically sound drugs."),
    ),

    # ---------------- Transplant and autoimmunity ----------------
    CatalogEntry(
        key="cd40lg_lupus", disease="Lupus and transplant rejection",
        area="Transplant",
        gene="CD40LG", protein_name="CD40 ligand (CD154)",
        uniprot="P29965", pdb_complex="6W9G",
        known_binder="humanised 5c8 Fab",
        rationale=(
            "Blocking the CD40 to CD40L costimulation signal interrupts "
            "T-cell help to B cells, which is central to lupus and to graft "
            "rejection. This exact structure is the engineered 5c8 Fab, "
            "redesigned to avoid the clotting problem described below."),
        clinical_note=(
            "The original anti-CD40L antibodies caused thromboembolism by "
            "crosslinking platelet Fc receptors, because CD40L is on "
            "platelets. The field returned only with Fc-silenced and fragment "
            "formats. A cautionary tale about where a target is expressed."),
    ),
    CatalogEntry(
        key="il2ra_transplant", disease="Transplant rejection and autoimmunity",
        area="Transplant",
        gene="IL2RA", protein_name="Interleukin-2 receptor subunit alpha (CD25)",
        uniprot="P01589", pdb_complex="6YIO",
        known_binder="RG6292 Fab",
        rationale=(
            "Basiliximab against CD25 is standard induction therapy in kidney "
            "transplantation. The same receptor is now being targeted the "
            "other way round in oncology, to deplete regulatory T cells "
            "without disturbing IL-2 signalling on effector cells."),
        clinical_note=(
            "Opposite goals on one target: tolerance in transplant, depletion "
            "in cancer. Whether a binder helps or harms depends entirely on "
            "which cells carry the receptor in that patient, which no "
            "sequence score can tell you."),
    ),

    # ---------------- Pain and neurology ----------------
    CatalogEntry(
        key="ngf_pain", disease="Chronic osteoarthritis pain", area="Pain",
        gene="NGF", protein_name="Beta-nerve growth factor",
        uniprot="P01138", pdb_complex="4EDW",
        known_binder="tanezumab Fab",
        rationale=(
            "Nerve growth factor sensitises pain neurons, and blocking it "
            "relieved osteoarthritis pain better than opioids in trials: a "
            "non-addictive analgesic mechanism, which is rare and badly "
            "needed."),
        clinical_note=(
            "Tanezumab was rejected by the FDA. Rapidly progressive "
            "osteoarthritis in a minority of patients led to joint "
            "replacements, and better pain relief turned out to come with "
            "faster joint destruction. Efficacy on the intended axis is not "
            "the same as benefit."),
    ),
    CatalogEntry(
        key="calcrl_migraine", disease="Migraine", area="Neurology",
        gene="CALCRL", protein_name="Calcitonin gene-related peptide receptor",
        uniprot="Q16602", pdb_complex="6UMG",
        known_binder="erenumab Fab",
        binder_preference="antibody",
        rationale=(
            "Erenumab is the only approved antibody against a GPCR, and this "
            "is its structure bound to the CGRP receptor's extracellular "
            "domain with RAMP1. Migraine prevention went from repurposed "
            "drugs to a mechanism-designed class on the strength of it."),
        clinical_note=(
            "Roughly half of patients reach a 50% reduction in migraine days, "
            "which is a real result and far from a cure. Constipation can be "
            "severe, and the long-term effect of blocking a vasodilatory "
            "pathway is still being watched."),
    ),
    CatalogEntry(
        key="tau_alzheimers", disease="Alzheimer disease", area="Neurology",
        gene="MAPT", protein_name="Microtubule-associated protein tau",
        uniprot="P10636", pdb_complex="5MO3",
        known_binder="DC8E8 Fab",
        rationale=(
            "Tau pathology tracks cognitive decline far better than amyloid "
            "does, which is why it is the second great Alzheimer target. "
            "This structure is a Fab bound to the 14-residue tau epitope it "
            "recognises."),
        clinical_note=(
            "Every anti-tau antibody tested in Alzheimer disease so far, "
            "gosuranemab, tilavonemab, semorinemab and zagotenemab among "
            "them, has failed to slow decline. The epitope here is a short "
            "peptide, so the binding site is small and the brief "
            "correspondingly thin."),
    ),
    CatalogEntry(
        key="trem2_alzheimers", disease="Alzheimer disease", area="Neurology",
        gene="TREM2", protein_name="Triggering receptor expressed on myeloid cells 2",
        uniprot="Q9NZC2", pdb_complex="8T51",
        known_binder="Fab 3.10C2",
        rationale=(
            "Rare TREM2 variants roughly triple Alzheimer risk, making it one "
            "of the few targets with direct human genetic support. Agonist "
            "antibodies aim to switch microglia into a protective state "
            "rather than block anything."),
        clinical_note=(
            "Genetic support says the protein matters, not that a drug "
            "against it helps. Microglial activation can be protective or "
            "damaging depending on stage, and the structure here covers only "
            "a 20-residue stalk peptide."),
    ),

    # ---------------- Musculoskeletal ----------------
    CatalogEntry(
        key="mstn_muscle", disease="Muscle wasting and sarcopenia",
        area="Musculoskeletal",
        gene="MSTN", protein_name="Myostatin (GDF-8)",
        uniprot="O14793", pdb_complex="5F3B",
        known_binder="RK35 chimeric antibody",
        rationale=(
            "Myostatin limits muscle growth; animals and people lacking it "
            "are visibly more muscular. Blocking it is the obvious approach "
            "to muscular dystrophy, sarcopenia, and now the muscle loss that "
            "comes with GLP-1 weight loss."),
        clinical_note=(
            "Myostatin blockade reliably adds muscle MASS and has repeatedly "
            "failed to add FUNCTION: trials in muscular dystrophy and "
            "sarcopenia missed their functional endpoints. Bigger muscles are "
            "not stronger patients."),
    ),

    # ---------------- Solid tumours ----------------
    CatalogEntry(
        key="cldn18_gastric", disease="Gastric and pancreatic cancer",
        area="Solid tumour",
        gene="CLDN18", protein_name="Claudin-18 isoform 2",
        uniprot="P56856", pdb_complex="9V2U",
        known_binder="osemitamab (TST001) Fab",
        rationale=(
            "Claudin-18.2 is a tight-junction protein normally buried between "
            "stomach cells and exposed only when tumours disrupt the "
            "architecture. Zolbetuximab's approval made it the first new "
            "gastric cancer target in years."),
        clinical_note=(
            "Nausea and vomiting are frequent and dose-limiting, because the "
            "normal stomach expresses the target too. Survival gains are "
            "measured in months. The structure is 3.8 A, so the binding site "
            "is approximate."),
    ),
    CatalogEntry(
        key="ceacam5_nsclc", disease="Non-small-cell lung cancer",
        area="Solid tumour",
        gene="CEACAM5",
        protein_name="Carcinoembryonic antigen-related cell adhesion molecule 5",
        uniprot="P06731", pdb_complex="8BW0",
        known_binder="tusamitamab Fab",
        rationale=(
            "CEA is the classic tumour marker, and antibody-drug conjugates "
            "finally made it a target rather than a blood test. The A3-B3 "
            "domain in this structure is the part tusamitamab grips."),
        clinical_note=(
            "Tusamitamab ravtansine was discontinued after its phase 3 "
            "interim analysis. CEACAM5 is also on normal gut and lung "
            "epithelium and has close relatives across the CEACAM family, "
            "which is precisely the cross-reactivity problem this software "
            "measures."),
    ),
    CatalogEntry(
        key="kit_gist", disease="GIST and systemic mastocytosis",
        area="Solid tumour",
        gene="KIT", protein_name="Mast/stem cell growth factor receptor Kit",
        uniprot="P10721", pdb_complex="4K94",
        known_binder="Fab19",
        rationale=(
            "KIT is where targeted cancer therapy began, with imatinib in "
            "GIST. Antibodies against its D4-D5 membrane-proximal region "
            "block the receptor dimerisation that activating mutations drive, "
            "a different mechanism from the kinase inhibitors that resistance "
            "mutations defeat."),
        clinical_note=(
            "KIT is required by haematopoietic stem cells, mast cells and "
            "melanocytes. Barzolvolimab causes hair depigmentation and "
            "neutropenia for exactly that reason: an on-target effect, not an "
            "off-target one."),
    ),

    # ---------------- Blood cancer ----------------
    CatalogEntry(
        key="gprc5d_myeloma", disease="Relapsed multiple myeloma",
        area="Blood cancer",
        gene="GPRC5D",
        protein_name="G protein-coupled receptor family C group 5 member D",
        uniprot="Q9NZD1", pdb_complex="9IMA",
        known_binder="talquetamab Fab",
        binder_preference="antibody",
        rationale=(
            "The myeloma target that is not BCMA, which matters because "
            "patients relapse after BCMA therapy. Talquetamab, a GPRC5D-CD3 "
            "bispecific, works after BCMA failure, and this structure is "
            "talquetamab's own Fab bound to the receptor."),
        clinical_note=(
            "GPRC5D is also expressed in hard keratinised tissue, so skin, "
            "nail and taste toxicity are common and sometimes "
            "treatment-limiting. Responses are deep but relapse is the rule."),
    ),

    # ---------------- Immuno-oncology ----------------
    CatalogEntry(
        key="tim3_mds", disease="Myelodysplastic syndrome and AML",
        area="Immuno-oncology",
        gene="HAVCR2", protein_name="TIM-3 (hepatitis A virus cellular receptor 2)",
        uniprot="Q8TDQ0", pdb_complex="7KQL",
        known_binder="Tim3.18 Fab",
        rationale=(
            "TIM-3 is a checkpoint on exhausted T cells and also on the "
            "leukaemic stem cells themselves, so blocking it could act on "
            "both sides at once. That dual expression is unusual and is why "
            "it was pursued in myeloid disease rather than solid tumours."),
        clinical_note=(
            "Sabatolimab's phase 3 programme in higher-risk MDS was "
            "discontinued for lack of efficacy. TIM-3 is also on regulatory "
            "and NK cells, so the net immunological effect is hard to predict "
            "from the target alone."),
    ),

    # ================= Third research pass =================
    # Every one of these had previously come back as "no usable structure".
    # Most of them had one all along: structure searches return entries
    # sorted by resolution, and a protein crystallised a hundred times on its
    # own buries its complexes behind all of them. Superoxide dismutase has
    # eleven usable complexes and reported none; PSMA reported a
    # three-residue dipeptide. The scan depth, not the data, was the limit.

    # ---------------- Neurology ----------------
    CatalogEntry(
        key="sod1_als", disease="Amyotrophic lateral sclerosis", area="Neurology",
        gene="SOD1", protein_name="Superoxide dismutase [Cu-Zn]",
        uniprot="P00441", pdb_complex="7NXX",
        known_binder="nanobody 2",
        rationale=(
            "SOD1 mutations were the first genetic cause of ALS found, and "
            "tofersen -- which lowers SOD1 protein -- is the first approved "
            "drug for a genetic form of the disease. The misfolded protein "
            "aggregates, so a binder that recognises it is a route to both "
            "imaging and clearance."),
        clinical_note=(
            "SOD1 mutations explain only about 2% of ALS. Tofersen lowers "
            "neurofilament convincingly and its effect on function is "
            "modest and debated. This is a target where the biology is "
            "certain and the benefit is not."),
    ),
    CatalogEntry(
        key="aqp4_nmo", disease="Neuromyelitis optica", area="Neurology",
        gene="AQP4", protein_name="Aquaporin-4",
        uniprot="P55087", pdb_complex="8V91",
        known_binder="pathogenic autoantibody rAB-58",
        binder_preference="antibody",
        rationale=(
            "Neuromyelitis optica is one of the few diseases where the "
            "pathogenic antibody is known exactly: patient IgG against "
            "aquaporin-4 destroys astrocytes. This structure is such an "
            "autoantibody from a patient, bound to its epitope -- so the "
            "binding site here is literally the site that causes the "
            "disease."),
        clinical_note=(
            "Unusually, the goal here is the opposite of the usual one. A "
            "binder that merely sticks to this epitope could block the "
            "autoantibody, or could mimic it and do the same damage. Which "
            "of those happens is decided by effector function and "
            "complement recruitment, neither of which this pipeline "
            "models at all."),
    ),
    CatalogEntry(
        key="chrna1_myasthenia", disease="Myasthenia gravis", area="Neurology",
        gene="CHRNA1", protein_name="Acetylcholine receptor subunit alpha",
        uniprot="P02708", pdb_complex="9DMS",
        known_binder="Fab6",
        binder_preference="antibody",
        rationale=(
            "The other classic antibody-mediated disease with a known "
            "autoantigen. Patient antibodies against the muscle "
            "acetylcholine receptor block neuromuscular transmission, and "
            "this is the human receptor with a Fab bound at 1.92 A -- one of "
            "the best-resolved autoimmune targets available."),
        clinical_note=(
            "The receptor is a five-subunit channel and this entry covers "
            "the alpha subunit only. Blocking the acetylcholine site would "
            "cause the disease rather than treat it; useful binders here "
            "would have to compete with patient antibodies without blocking "
            "the neurotransmitter, which is a far subtler requirement than a "
            "binding score can express."),
    ),

    # ---------------- Autoimmune ----------------
    CatalogEntry(
        key="fcrn_myasthenia", disease="Myasthenia gravis and ITP",
        area="Autoimmune",
        gene="FCGRT", protein_name="IgG receptor FcRn large subunit p51",
        uniprot="P55899", pdb_complex="9MI6",
        known_binder="nipocalimab Fab",
        binder_preference="antibody",
        rationale=(
            "FcRn is what gives IgG its three-week half-life by rescuing it "
            "from degradation. Blocking it strips pathogenic autoantibodies "
            "out of the circulation within days -- a mechanism that works "
            "for any antibody-driven disease regardless of what the antibody "
            "targets. Efgartigimod and nipocalimab are both approved, and "
            "this is nipocalimab's own structure."),
        clinical_note=(
            "It lowers ALL IgG, not just the harmful kind, so infection risk "
            "and loss of vaccine protection are the trade. The effect is "
            "also temporary: antibody levels recover between cycles."),
    ),

    # ---------------- Haematology ----------------
    CatalogEntry(
        key="vwf_thrombosis", disease="Thrombotic thrombocytopenic purpura",
        area="Haematology",
        gene="VWF", protein_name="von Willebrand factor A3 domain",
        uniprot="P04275", pdb_complex="2ADF",
        known_binder="82D6A3 antithrombotic antibody",
        rationale=(
            "In acquired TTP, von Willebrand factor multimers go uncleaved "
            "and catch platelets in small vessels. Caplacizumab, a nanobody "
            "against the A1 domain, blocks that and was the first "
            "nanobody-derived drug approved anywhere. This structure is an "
            "antithrombotic antibody on the collagen-binding A3 domain."),
        clinical_note=(
            "Bleeding is the obvious and observed risk of blocking the "
            "protein that makes platelets stick. Caplacizumab shortens time "
            "to platelet recovery; its effect on mortality is less clear."),
    ),
    CatalogEntry(
        key="c1s_cad", disease="Cold agglutinin disease", area="Haematology",
        gene="C1S", protein_name="Complement C1s subcomponent",
        uniprot="P09871", pdb_complex="9X1H",
        known_binder="RAY121 Fab",
        binder_preference="antibody",
        rationale=(
            "C1s is the protease that carries out the classical complement "
            "pathway's first cut. Sutimlimab blocks it and stops the "
            "antibody-driven haemolysis of cold agglutinin disease while "
            "leaving the alternative pathway -- and so most antibacterial "
            "defence -- intact. A good example of blocking one pathway "
            "rather than the whole system."),
        clinical_note=(
            "Classical-pathway blockade still raises the risk of "
            "encapsulated bacterial infection, and vaccination is required. "
            "Cold agglutinin disease is rare, so the evidence base is small."),
    ),

    # ---------------- Rare disease ----------------
    CatalogEntry(
        key="klkb1_hae", disease="Hereditary angioedema", area="Rare disease",
        gene="KLKB1", protein_name="Plasma kallikrein",
        uniprot="P03952", pdb_complex="4OGY",
        known_binder="DX-2930 (lanadelumab) Fab",
        rationale=(
            "Hereditary angioedema is unchecked bradykinin production, and "
            "plasma kallikrein is the enzyme that makes it. Lanadelumab "
            "blocks the enzyme and prevents attacks; this structure is "
            "lanadelumab's own Fab bound to it. One of the cleanest "
            "mechanism-to-drug stories in rare disease."),
        clinical_note=(
            "Prophylaxis, not rescue: it prevents attacks rather than "
            "treating one in progress, and an airway attack remains an "
            "emergency. Inhibiting a protease of the contact pathway also "
            "has theoretical bleeding implications."),
    ),
    CatalogEntry(
        key="acvr1_fop", disease="Fibrodysplasia ossificans progressiva",
        area="Rare disease",
        gene="ACVR1", protein_name="Activin receptor type-1 (ALK2)",
        uniprot="Q04771", pdb_complex="7YRU",
        known_binder="anti-ALK2 Fab",
        rationale=(
            "FOP turns muscle into bone after injury, and the cause is a "
            "single recurrent mutation in ALK2 that makes it respond to "
            "activin A as though it were a bone-forming signal. One "
            "mutation, one receptor, one gained function -- which is as "
            "clean as human genetics gets."),
        clinical_note=(
            "Garetosmab, which blocks activin A rather than the receptor, "
            "reduced new bone but was placed on clinical hold after deaths "
            "in its trial. The receptor is required for normal bone and "
            "iron metabolism, so blocking it is not free."),
    ),
    CatalogEntry(
        key="fgf23_xlh", disease="X-linked hypophosphataemia", area="Rare disease",
        gene="FGF23", protein_name="Fibroblast growth factor 23",
        uniprot="Q9GZV9", pdb_complex="7YSH",
        known_binder="FGFR1c and alpha-klotho",
        rationale=(
            "Excess FGF23 makes the kidney waste phosphate, which softens "
            "growing bone. Burosumab blocks it and corrects the rickets of "
            "X-linked hypophosphataemia -- a childhood disease with a "
            "mechanism-matched antibody. The structure is FGF23 held by the "
            "receptor and co-receptor it has to engage."),
        clinical_note=(
            "The bound partner here is the physiological receptor complex, "
            "not a drug, so the interface is the signalling one. "
            "Over-blocking causes hyperphosphataemia and ectopic "
            "calcification, which is why dosing is titrated on blood "
            "phosphate."),
    ),
    CatalogEntry(
        key="cfb_complement", disease="Complement-mediated kidney disease",
        area="Rare disease",
        gene="CFB", protein_name="Complement factor B",
        uniprot="P00751", pdb_complex="3HRZ",
        known_binder="cobra venom factor",
        rationale=(
            "Factor B is the catalytic half of the alternative pathway's C3 "
            "convertase, and iptacopan -- which inhibits the pathway -- is "
            "now approved for PNH and for IgA nephropathy and C3 "
            "glomerulopathy. The bound partner here is cobra venom factor, a "
            "C3b analogue from snake venom that forms the same complex: an "
            "evolved C3b mimic, which makes it a template for what a "
            "designed binder has to do."),
        clinical_note=(
            "Blocking the alternative pathway raises the risk of "
            "encapsulated bacterial infection and requires vaccination. The "
            "interface here is the one factor B makes with C3b, not a drug "
            "site anybody has validated."),
    ),
    CatalogEntry(
        key="masp2_tma", disease="Thrombotic microangiopathy", area="Rare disease",
        gene="MASP2", protein_name="Mannan-binding lectin serine protease 2",
        uniprot="O00187", pdb_complex="7PQN",
        known_binder="ecotin",
        rationale=(
            "MASP-2 is the effector protease of the lectin complement "
            "pathway, the arm that activates on carbohydrate patterns "
            "without needing an antibody. Narsoplimab blocks it for "
            "transplant-associated thrombotic microangiopathy, where "
            "complement damages small vessels."),
        clinical_note=(
            "Narsoplimab's path to approval has been long and contested. "
            "The bound partner here is ecotin, a bacterial protease "
            "inhibitor, so the interface is the active-site groove rather "
            "than a therapeutic epitope."),
    ),

    # ---------------- Solid tumours ----------------
    CatalogEntry(
        key="psma_prostate", disease="Metastatic prostate cancer",
        area="Solid tumour",
        gene="FOLH1", protein_name="Prostate-specific membrane antigen (PSMA)",
        uniprot="Q04609", pdb_complex="9HVL",
        known_binder="nanobody 37",
        rationale=(
            "PSMA is the basis of both modern prostate cancer imaging and "
            "lutetium-177 radioligand therapy, which improves survival in "
            "metastatic disease. A small, stable binder is exactly what "
            "radioligand therapy needs -- it has to reach the tumour and "
            "clear quickly -- and this structure is a nanobody bound to it."),
        clinical_note=(
            "PSMA is also on kidney, salivary gland and small bowel, which "
            "is why dry mouth is the characteristic toxicity of PSMA "
            "radioligands. The target is excellent and not exclusive."),
    ),
    CatalogEntry(
        key="met_nsclc", disease="MET-altered lung cancer", area="Solid tumour",
        gene="MET", protein_name="Hepatocyte growth factor receptor",
        uniprot="P08581", pdb_complex="4K3J",
        known_binder="onartuzumab Fab",
        binder_preference="antibody",
        rationale=(
            "MET amplification and exon-14 skipping drive a slice of lung "
            "cancer and are a common escape route from EGFR inhibitors. "
            "Amivantamab, a bispecific that engages both EGFR and MET, is "
            "approved on exactly that logic."),
        clinical_note=(
            "This structure is onartuzumab, which engaged MET well and "
            "failed its phase 3 trial -- the antibody worked and the "
            "patients did not benefit. A useful reminder that a good "
            "structure and good target engagement decide very little."),
    ),

    # ---------------- Immuno-oncology ----------------
    CatalogEntry(
        key="cd73_io", disease="Solid tumours", area="Immuno-oncology",
        gene="NT5E", protein_name="Ecto-5'-nucleotidase (CD73)",
        uniprot="P21589", pdb_complex="7BBJ",
        known_binder="mAb19",
        rationale=(
            "Tumours suppress the immune cells around them by flooding the "
            "tissue with adenosine, and CD73 is the enzyme that makes it. "
            "Blocking it is a way to relieve immunosuppression that does not "
            "go through the PD-1 axis, which matters because most patients "
            "do not respond to that axis."),
        clinical_note=(
            "Adenosine-pathway drugs have so far produced modest results as "
            "single agents. CD73 is also on normal lymphocytes and "
            "endothelium, and adenosine signalling is protective in several "
            "tissues."),
    ),
    CatalogEntry(
        key="vista_io", disease="Solid tumours", area="Immuno-oncology",
        gene="VSIR", protein_name="V-domain Ig suppressor of T-cell activation",
        uniprot="Q9H7M9", pdb_complex="6MVL",
        known_binder="pH-selective Fab",
        rationale=(
            "VISTA suppresses T cells and is most active in the acidic "
            "interior of a tumour. The antibody in this structure binds "
            "selectively at low pH, so it engages inside the tumour and not "
            "in healthy tissue -- targeting by local chemistry rather than "
            "by expression, which is an unusual and clever idea."),
        clinical_note=(
            "The first anti-VISTA antibody caused severe cytokine release "
            "and the programme stopped; pH-selective designs are the "
            "response to that. None is approved."),
    ),

    # ---------------- Respiratory ----------------
    CatalogEntry(
        key="itgav_fibrosis", disease="Fibrosis and solid tumours",
        area="Respiratory",
        gene="ITGAV", protein_name="Integrin alpha-V",
        uniprot="P06756", pdb_complex="9CZD",
        known_binder="17E6 Fab",
        binder_preference="antibody",
        rationale=(
            "The alpha-v integrins activate latent TGF-beta directly in "
            "injured epithelium, which makes them the local switch for "
            "fibrosis rather than a systemic one. Blocking them aims to stop "
            "fibrosis where it starts without suppressing TGF-beta "
            "everywhere. The entry is the alpha-V subunit because that is "
            "the one 17E6 in this structure actually binds -- paired with "
            "beta-6 it had no antibody contact at all."),
        clinical_note=(
            "A trial of an alpha-v-beta-6 inhibitor in pulmonary fibrosis "
            "was stopped early for worsening outcomes, and abituzumab "
            "against alpha-v missed its endpoint in colorectal cancer. "
            "Alpha-v pairs with five different beta subunits across many "
            "tissues, so blocking the shared subunit is far less selective "
            "than blocking one pair."),
    ),
    CatalogEntry(
        key="il5ra_asthma", disease="Severe eosinophilic asthma",
        area="Respiratory",
        gene="IL5RA", protein_name="Interleukin-5 receptor subunit alpha",
        uniprot="Q01344", pdb_complex="3QT2",
        known_binder="interleukin-5",
        rationale=(
            "The receptor half of the IL-5 axis. Benralizumab binds it and "
            "additionally recruits NK cells to kill eosinophils outright, "
            "which depletes them almost completely rather than merely "
            "neutralising the cytokine. This structure shows IL-5 held in "
            "the site a blocker has to cover."),
        clinical_note=(
            "Near-complete eosinophil depletion is more than cytokine "
            "blockade achieves, and whether that translates into better "
            "outcomes than anti-IL-5 antibodies is still argued. Eosinophils "
            "have roles in defence that chronic depletion may matter for."),
    ),
    CatalogEntry(
        key="il11_fibrosis", disease="Fibrotic disease", area="Respiratory",
        gene="IL11", protein_name="Interleukin-11",
        uniprot="P20809", pdb_complex="8QY4",
        known_binder="IL-11 receptor alpha and gp130",
        rationale=(
            "IL-11 was thought to be anti-inflammatory for twenty years; the "
            "work of the last decade indicates it is a central driver of "
            "fibrosis across heart, lung, liver and kidney, and blocking it "
            "reverses fibrosis in animals. A target whose sign in the "
            "literature flipped."),
        clinical_note=(
            "No anti-IL-11 drug is approved and the human evidence is thin. "
            "The structure is the signalling complex at 3.1 A, so the site "
            "is the receptor interface and is approximate."),
    ),

    # ---------------- Allergy ----------------
    CatalogEntry(
        key="siglec8_eoe", disease="Eosinophilic gastritis", area="Allergy",
        gene="SIGLEC8", protein_name="Sialic acid-binding Ig-like lectin 8",
        uniprot="Q9NYZ4", pdb_complex="7QUH",
        known_binder="lirentelimab (AK002) Fab",
        rationale=(
            "Siglec-8 sits almost exclusively on eosinophils and mast cells, "
            "the two cell types that drive allergic tissue damage, and "
            "engaging it tells them to die or stand down. That expression "
            "pattern is about as selective as a human cell-surface target "
            "gets, and this structure is lirentelimab's own Fab."),
        clinical_note=(
            "Lirentelimab depleted eosinophils convincingly and missed its "
            "symptom endpoints in eosinophilic gastritis and oesophagitis. "
            "Removing the cells and relieving the disease turned out not to "
            "be the same thing."),
    ),

    # ---------------- Bone ----------------
    CatalogEntry(
        key="sost_osteoporosis", disease="Osteoporosis", area="Bone",
        gene="SOST", protein_name="Sclerostin",
        uniprot="Q9BQB4", pdb_complex="6L6R",
        known_binder="LRP6 E1E2",
        rationale=(
            "Sclerostin is the brake on bone formation; people who lack it "
            "develop dense, strong bone. Romosozumab releases that brake and "
            "is the only approved drug that both builds bone and reduces "
            "resorption. The structure is sclerostin held by LRP6, the "
            "receptor it blocks."),
        clinical_note=(
            "Romosozumab carries a boxed warning for cardiovascular events, "
            "which is why it is reserved for high fracture risk. The "
            "structure is 3.8 A, so this binding site is approximate."),
    ),

    # ---------------- Fibrosis / oncology ----------------
    CatalogEntry(
        key="tgfb1_fibrosis", disease="Fibrosis and immune evasion",
        area="Musculoskeletal",
        gene="TGFB1", protein_name="Transforming growth factor beta-1",
        uniprot="P01137", pdb_complex="8UDZ",
        known_binder="LTBP-49247 Fab",
        rationale=(
            "TGF-beta drives fibrosis in every organ and shields tumours "
            "from immune attack. The interesting antibodies do not "
            "neutralise the free cytokine but lock the LATENT complex shut, "
            "preventing activation where it happens -- which is what this "
            "structure shows."),
        clinical_note=(
            "Broad TGF-beta blockade caused heart-valve lesions in animals "
            "and stopped several programmes. The cytokine is essential for "
            "immune tolerance and wound healing, so the therapeutic window "
            "is the entire difficulty."),
    ),
]


CANDIDATES_BY_KEY = {e.key: e for e in CANDIDATES}


# --- Loading the verified catalogue -----------------------------------------


def load_verified() -> dict:
    """The catalogue of targets that passed verification.

    Returns {"targets": {...}, "rejected": [...], "built_at": ...}.
    """
    if not VERIFIED_PATH.exists():
        return {"targets": {}, "rejected": [], "built_at": None}
    return json.loads(VERIFIED_PATH.read_text(encoding="utf-8"))


def verified_keys() -> list[str]:
    return sorted(load_verified().get("targets", {}))


def areas() -> dict[str, list[dict]]:
    """Verified targets grouped by disease area, for the catalogue UI."""
    cat = load_verified().get("targets", {})
    out: dict[str, list[dict]] = {}
    for key, entry in cat.items():
        out.setdefault(entry["area"], []).append(entry)
    for v in out.values():
        v.sort(key=lambda e: e["disease"])
    return dict(sorted(out.items()))
