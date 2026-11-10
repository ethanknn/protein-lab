"""Self-contained tests. No pytest needed:  python tests/test_pipeline.py

These guard the things that would silently corrupt results rather than crash:
incomplete property tables, a parser that drops chains, a scoring function that
stops discriminating, a loop that does not actually improve anything.
"""

from __future__ import annotations

import random
import copy
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from plab import aminoacids as aa  # noqa: E402
from plab.align import smith_waterman  # noqa: E402
from plab.design.loop import DesignLoop, LoopConfig  # noqa: E402
from plab.design.scaffolds import get_scaffold, list_scaffolds  # noqa: E402
from plab.pdbio import find_contacts, parse_pdb  # noqa: E402
from plab.scoring.composite import Evaluator  # noqa: E402
from plab.targets import load_off_targets, load_target  # noqa: E402

STRUCT = ROOT / "data" / "structures"
SEQ = ROOT / "data" / "sequences"

_passed, _failed = 0, []


def check(name: str, cond: bool, detail: str = "") -> None:
    global _passed
    if cond:
        _passed += 1
        print(f"  pass  {name}")
    else:
        _failed.append(name)
        print(f"  FAIL  {name}  {detail}")


def section(t: str) -> None:
    print(f"\n{t}")


# --- property tables --------------------------------------------------------

def test_tables() -> None:
    section("amino acid tables")
    tables = {
        "RESIDUE_MASS": aa.RESIDUE_MASS, "KD_HYDROPATHY": aa.KD_HYDROPATHY,
        "EISENBERG": aa.EISENBERG, "HELIX_PROPENSITY": aa.HELIX_PROPENSITY,
        "SHEET_PROPENSITY": aa.SHEET_PROPENSITY, "RESIDUE_VOLUME": aa.RESIDUE_VOLUME,
        "MAX_ASA": aa.MAX_ASA, "FULL_NAME": aa.FULL_NAME,
    }
    for name, tbl in tables.items():
        missing = set(aa.AA20) - set(tbl)
        check(f"{name} covers all 20", not missing, f"missing {sorted(missing)}")

    # Glycine is the lightest residue, tryptophan the heaviest.
    check("mass ordering G < W",
          aa.RESIDUE_MASS["G"] < aa.RESIDUE_MASS["W"])
    # Alanine monomer is 89.09 Da; residue mass + water should recover it.
    check("MW of 'A' ~ 89.09", abs(aa.molecular_weight("A") - 89.09) < 0.02,
          f"got {aa.molecular_weight('A'):.3f}")
    # A known tripeptide: Gly-Gly-Gly = 189.17 Da
    check("MW of 'GGG' ~ 189.17", abs(aa.molecular_weight("GGG") - 189.17) < 0.02,
          f"got {aa.molecular_weight('GGG'):.3f}")


def test_charge() -> None:
    section("charge model")
    check("poly-K is positive at pH 7.4", aa.net_charge("K" * 10, 7.4) > 8)
    check("poly-E is negative at pH 7.4", aa.net_charge("E" * 10, 7.4) < -8)
    check("poly-A is near neutral", abs(aa.net_charge("A" * 10, 7.4)) < 0.5)
    check("charge decreases as pH rises",
          aa.net_charge("KKKEEE", 3.0) > aa.net_charge("KKKEEE", 11.0))
    check("pI of poly-K is basic", aa.isoelectric_point("K" * 10) > 10)
    check("pI of poly-E is acidic", aa.isoelectric_point("E" * 10) < 5)
    # At its own pI a sequence should carry ~zero net charge.
    seq = "VDNKFNKEQQNAFYEILHLPNLNEEQRNAFIQSLKDDPSQSANLLAEAKKLNDAQAPK"
    pi = aa.isoelectric_point(seq)
    check("net charge ~0 at the computed pI", abs(aa.net_charge(seq, pi)) < 0.05,
          f"charge {aa.net_charge(seq, pi):.4f} at pI {pi}")


def test_descriptors() -> None:
    section("sequence descriptors")
    check("GRAVY: poly-I hydrophobic", aa.gravy("I" * 10) > 4)
    check("GRAVY: poly-R hydrophilic", aa.gravy("R" * 10) < -4)
    check("entropy of homopolymer is 0", aa.shannon_entropy("A" * 20) == 0)
    check("entropy of all-20 is max",
          abs(aa.shannon_entropy(aa.AA20) - 4.321928) < 1e-4)
    check("homopolymer run detected", aa.max_homopolymer_run("ABAAAAB") == 4)
    check("aromaticity counts F/W/Y", abs(aa.aromaticity("FWYA") - 0.75) < 1e-9)
    # An ideal amphipathic helix should beat a uniformly hydrophobic one.
    amphi = aa.hydrophobic_moment("LKKLLKLLKKLLKLL")
    uniform = aa.hydrophobic_moment("LLLLLLLLLLLLLLL")
    check("amphipathic helix has higher moment", amphi > uniform,
          f"{amphi:.3f} vs {uniform:.3f}")


# --- parsing ----------------------------------------------------------------

def test_parser() -> None:
    section("PDB parser")
    p = STRUCT / "8W9J.pdb"
    if not p.exists():
        check("8W9J present", False, "run scripts/fetch_data.py")
        return
    s = parse_pdb(p)
    check("8W9J has 6 chains", len(s.chain_ids) == 6, f"got {s.chain_ids}")
    check("resolution parsed", s.resolution == 3.5, f"got {s.resolution}")
    check("CLEC12A chain C is 155 residues", len(s.chain("C")) == 155,
          f"got {len(s.chain('C'))}")
    check("sequence has no X", "X" not in s.sequence("C"))
    check("every residue has a CA",
          all(r.ca is not None for r in s.chain("C")))
    # Contacts must be symmetric in count.
    c1 = find_contacts(s.chain("C"), s.chain("H"))
    c2 = find_contacts(s.chain("H"), s.chain("C"))
    check("contact detection is symmetric", len(c1) == len(c2),
          f"{len(c1)} vs {len(c2)}")
    check("all contacts within cutoff",
          all(c.min_distance <= 4.5 for c in c1))
    # A chain cannot contact itself across a 4.5 A gap in this sense.
    check("distinct chains found contacts", len(c1) > 0)


def test_align() -> None:
    section("alignment")
    seq = "VDNKFNKEQQNAFYEILHLPNLNEEQRNAFIQ"
    a = smith_waterman(seq, seq)
    check("self-alignment is 100% identity", a.identity == 1.0,
          f"got {a.identity}")
    check("self-alignment spans the sequence", a.length == len(seq),
          f"got {a.length}")
    b = smith_waterman("AAAAAAAA", "WWWWWWWW")
    check("dissimilar sequences score low", b.score < a.score)
    c = smith_waterman("KLMNPQ", "XXXKLMNPQXXX".replace("X", "G"))
    check("substring is located", c.identity >= 0.9, f"got {c.identity}")


# --- scaffolds --------------------------------------------------------------

def test_scaffolds() -> None:
    section("scaffolds")
    for sc in list_scaffolds():
        check(f"{sc.name}: ss length matches framework",
              len(sc.secondary_structure) == len(sc.framework))
        check(f"{sc.name}: variable positions in range",
              all(0 <= p < sc.length for p in sc.variable_positions))
        check(f"{sc.name}: face length matches variable count",
              len(sc.face(sc.framework)) == len(sc.interface_positions))

    z = get_scaffold("affibody_Z")
    check("Z domain is 58 residues", z.length == 58, f"got {z.length}")
    check("Z domain has 13 designed positions", z.n_variable == 13)
    check("Z domain is cysteine-free", "C" not in z.framework)
    # Frozen prolines must survive an apply() that tries to change them.
    mutated = z.apply(z.framework, {p: "W" for p in range(z.length)})
    for p in z.frozen_positions:
        check(f"frozen position {p + 1} preserved",
              mutated[p] == z.framework[p],
              f"{z.framework[p]} -> {mutated[p]}")
    check("apply changed non-frozen positions", mutated != z.framework)
    # annotate must not alter the sequence, only its case.
    check("annotate preserves sequence",
          z.annotate(z.framework).upper() == z.framework.upper())


# --- scoring ----------------------------------------------------------------

def test_scoring() -> None:
    section("scoring")
    target = load_target()
    check("epitope has 13 residues", len(target.epitope) == 13,
          f"got {len(target.epitope)}")
    check("epitope has 3 arginines", len(target.arginine_positions) == 3)
    check("epitope net charge is positive", target.epitope_net_charge > 0)
    check("off-target panel loaded", len(load_off_targets(target)) >= 2)

    ev = Evaluator(target, use_plm=False)
    sc = get_scaffold("affibody_Z")

    check("reference molecule built", ev.reference is not None)
    check("paralog surfaces mapped", len(ev.paralog_epitopes) >= 1,
          f"got {ev.paralog_epitopes}")

    def recog(seq: str) -> float:
        card = ev.evaluate(seq, sc)
        return ev.group_scores(card, exclude_self_referential=True).get(
            "recognition", 0.0) * 100

    # The epitope is arginine-rich, so an all-arginine face must do badly and
    # an acidic aromatic face must do well. If this ordering ever inverts the
    # scoring function has broken.
    poly_r = recog(sc.apply(sc.framework, {p: "R" for p in sc.variable_positions}))
    poly_e = recog(sc.apply(sc.framework, {p: "E" for p in sc.variable_positions}))
    designed = recog("VDNKFNKEFNKAYREITQLPNLNDDQYDAFIESLYDDPSQSANLLAEAKKLNDAQAPK")
    wild = recog(sc.framework)

    check("designed face beats poly-Arg", designed > poly_r,
          f"{designed:.1f} vs {poly_r:.1f}")
    check("designed face beats poly-Glu", designed > poly_e,
          f"{designed:.1f} vs {poly_e:.1f}")
    check("designed face beats wild-type", designed > wild,
          f"{designed:.1f} vs {wild:.1f}")
    check("poly-Arg scores poorly in absolute terms", poly_r < 40,
          f"got {poly_r:.1f}")

    # Every component must carry an explanation. A silent number is a bug.
    card = ev.evaluate(sc.framework, sc)
    check("every component has a detail string",
          all(c.detail for c in card.components))
    check("all normalised scores in [0,1]",
          all(0.0 <= c.normalised <= 1.0 for c in card.components))
    check("composite in [0,100]", 0 <= card.score_100 <= 100)


def test_selectivity() -> None:
    """Regression guard for the third saturating-component bug.

    The paralog kernel used min(supply, demand)/demand, which caps at 1.0.
    Once a face had enough acidic residues for CLEC12A's three arginines it
    also had enough for every paralog, so both sides scored 1.0 and the
    margin collapsed. A sweep of 680 faces found a maximum reachable margin
    of 0.075 against a ramp that awarded full marks at 0.25 -- so 50% was a
    mathematical ceiling and no design could ever beat it.

    The rule this enforces: the component must produce a real gradient, and
    must saturate at neither end.
    """
    import random as _random

    section("paralog selectivity")
    from plab.scoring.selectivity import (complementarity_kernel,
                                          margin_to_score)

    target = load_target()
    ev = Evaluator(target, use_plm=False)
    sc = get_scaffold("affibody_Z")

    paralogs = ev.paralog_epitopes
    check("several paralog surfaces mapped", len(paralogs) >= 4,
          f"got {sorted(paralogs)}")
    check("CLEC12B is in the panel", "CLEC12B" in paralogs,
          "the nearest paralog is the most important comparison")
    check("CD33 is NOT aligned as a structural homolog",
          "CD33" not in paralogs,
          "CD33 is a Siglec, a different fold family")

    lo, hi = ev.selectivity_scale(13)
    check("scale is ordered", lo < hi, f"lo {lo:.3f} hi {hi:.3f}")
    check("scale spans a usable range", hi - lo > 0.2,
          f"span {hi - lo:.3f}")

    # The kernel must separate the target from the paralogs for a face built
    # to the target's specific demands.
    epi = target.epitope_sequence
    tuned = "DDDYYYKSSSSSS"  # 3 acidic for 3 arginines, 3 aromatics to stack
    on = complementarity_kernel(tuned, epi)
    offs = {g: complementarity_kernel(tuned, p) for g, p in paralogs.items()}
    check("a target-tuned face beats every paralog",
          all(on > v for v in offs.values()),
          f"on {on:.3f} vs {({k: round(v, 3) for k, v in offs.items()})}")

    # Excess acidic residues must cost selectivity, which is the whole point
    # of penalising overshoot: unpaired charge goes and binds something else.
    lean = complementarity_kernel("DDDYYYSSSSSSS", epi)
    excess = complementarity_kernel("DDDDDDDDYYYSS", epi)
    check("excess acidic residues score worse on target", excess < lean,
          f"excess {excess:.3f} vs matched {lean:.3f}")

    # Gradient check: random faces must span many distinct scores.
    rng = _random.Random(3)
    scores = []
    for _ in range(120):
        face = "".join(rng.choice(aa.AA20) for _ in range(13))
        m = complementarity_kernel(face, epi) - max(
            complementarity_kernel(face, p) for p in paralogs.values())
        scores.append(round(100 * margin_to_score(m, lo, hi), 0))
    distinct = len(set(scores))
    check("component produces a real gradient", distinct >= 15,
          f"only {distinct} distinct scores across 120 random faces")
    check("does not saturate at the floor",
          sum(1 for s in scores if s == 0) < len(scores) * 0.6,
          f"{sum(1 for s in scores if s == 0)}/{len(scores)} at zero")
    check("does not saturate at the ceiling",
          sum(1 for s in scores if s == 100) < len(scores) * 0.6,
          f"{sum(1 for s in scores if s == 100)}/{len(scores)} at 100")

    # End to end: a designed candidate must beat the wild-type scaffold.
    designed = ev.evaluate(
        "VDNKFNKEYYEAYDEIKKLPNLNQEQQKAFIDSLWDDPSQSANLLAEAKKLNDAQAPK", sc
    ).component("paralog_selectivity")
    wild = ev.evaluate(sc.framework, sc).component("paralog_selectivity")
    check("designed candidate beats wild-type on selectivity",
          designed.normalised > wild.normalised,
          f"{designed.normalised:.3f} vs {wild.normalised:.3f}")
    check("designed candidate is not pinned at the ceiling",
          designed.normalised < 1.0,
          "no headroom left for the optimiser to climb")


def test_chain_matching() -> None:
    """A protein must be recognised as itself across species and strains.

    Five catalogue candidates were rejected because this was decided on
    exact 8-mer matching, which needs all eight positions identical and so
    scores a true ortholog at p**8. These checks pin the failure directly:
    the brittle measure must be seen to break, and the measure that replaced
    it must be seen to hold, or the bug comes back the next time someone
    "simplifies" the matching.
    """
    section("chain matching across variants")
    from plab.epitope import (
        ALIGNED_IDENTITY_THRESHOLD, aligned_match, kmer_identity, match_target,
    )

    real = load_target("clec12a_aml")
    seq = real.epitope_sequence
    # A real protein sequence long enough for both measures to apply.
    native = ("MSEEVTYADLQFQNSSEMEKIPEIGKFGEKAPPAPSHVWRPAALFLTLLCLLLLIGLGVLASMFHVTLKI"
              "EMKKMNKLQNISEELQRNISLQLMSNMNISNKIRNLSTTLQTIATKLCRELYSKEQEHKCKPCPRRWIWH"
              "KDSCYFLSDDVQTWQESKMACAAQNASLLKINNKNALEFIKSQSRSYDYWLGLSPEEDSTRGMRVDRLYQ")
    variant = "".join(
        ("A" if i % 11 == 0 and c not in "C" else c) for i, c in enumerate(native))
    identity = sum(1 for a, b in zip(native, variant) if a == b) / len(native)

    km = kmer_identity(variant, native)
    check("a 9-in-10 identical variant breaks exact k-mer matching",
          identity > 0.9 and km < 0.5,
          f"identity {identity:.2f} but k-mer {km:.2f}")
    check("alignment recognises the same variant as the protein",
          match_target(variant, native).is_match,
          f"aligned identity {aligned_match(variant, native).identity:.2f}")
    check("aligned identity is close to the true identity",
          abs(aligned_match(variant, native).identity - identity) < 0.08)

    unrelated = ("QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQ"
                 "KFQGRVTITADESTSTAYMELSSLRSEDTAVYYCARDRGDAFDIWGQGTMVTVSS")
    um = match_target(unrelated, native)
    check("an unrelated protein is still rejected", not um.is_match,
          f"antibody chain matched at {um.identity:.2f}")
    check("the gap between ortholog and stranger is wide",
          um.identity + 0.15 < ALIGNED_IDENTITY_THRESHOLD < identity,
          f"stranger {um.identity:.2f}, threshold "
          f"{ALIGNED_IDENTITY_THRESHOLD}, ortholog {identity:.2f}")

    peptide = seq
    check("a short peptide is matched by containment",
          match_target(peptide, native).method == "contained")


def test_binder_grouping() -> None:
    """Two antibodies bound to one antigen must not become one binder.

    Every rule tried here failed on a real structure before this one: motifs
    split Fabs, chain packing merged neighbouring Fabs, and clustering the
    epitope split the elongated RBD/ACE2 interface. The checks below are the
    three failures, written down.
    """
    section("telling one binder from two")
    from plab.epitope import MIN_PAIR_CONTACTS, group_binders

    class R:
        def __init__(self, xyz):
            self.centroid = xyz

    near_a = [R((0.0, 0.0, 0.0)), R((4.0, 0.0, 0.0))]
    near_b = [R((9.0, 0.0, 0.0)), R((13.0, 0.0, 0.0))]
    far = [R((60.0, 0.0, 0.0))]

    heavy_light = {"H": near_a, "L": near_b}
    pair = {("H", "L"): MIN_PAIR_CONTACTS + 20}
    check("a heavy/light pair is one binder",
          group_binders(heavy_light, pair) == [["H", "L"]])

    check("chains that only graze each other stay separate",
          sorted(group_binders(heavy_light, {("H", "L"): 2})) == [["H"], ["L"]])

    check("adjacent patches alone do not merge chains",
          sorted(group_binders(heavy_light, {})) == [["H"], ["L"]])

    distant = {"H": near_a, "X": far}
    check("a second antibody at a distant epitope stays separate",
          sorted(group_binders(distant, {("H", "X"): 400})) == [["H"], ["X"]])

    copies = {"H": near_a, "H2": near_b}
    check("two copies of one chain are never one molecule",
          sorted(group_binders(copies, {("H", "H2"): 400},
                               {("H", "H2")})) == [["H"], ["H2"]])


def test_numbering_offset() -> None:
    """Structure numbering and reference numbering are often not the same.

    Assuming they agree silently refused to map the binding site of 16 of 42
    targets -- everything whose depositors numbered by the mature protein
    rather than the precursor. The offsets recovered here are checkable
    against biology: amyloid-beta is residues 672-713 of APP, mature TNF
    starts at 77, HER2 carries a 22-residue signal peptide.
    """
    section("structure vs reference numbering")
    from plab.scoring.selectivity import epitope_window, numbering_offset

    class FakeRes:
        def __init__(self, resseq, aa):
            self.resseq, self.aa = resseq, aa

    class FakeTarget:
        def __init__(self, full, residues):
            self.full_sequence, self.epitope = full, residues

    full = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEENFKALVLIAFAQYLQQCPF"
    # A structure numbered from residue 1 of the MATURE protein, i.e. 18
    # positions along the precursor.
    shift = 18
    residues = [FakeRes(i + 1 - shift, full[i]) for i in range(20, 32)]
    got = numbering_offset(FakeTarget(full, residues))
    check("a constant numbering shift is recovered", got is not None and got[0] == shift,
          f"got {got}")
    check("and every epitope residue then matches", got is not None and got[1] == 1.0)

    scrambled = [FakeRes(i * 7 + 3, "W") for i in range(12)]
    w = epitope_window(FakeTarget(full, scrambled))
    check("numbering that cannot be reconciled is refused, not guessed",
          w is None)

    known = {"abeta_alzheimers": 671, "tnf_ra": 76, "her2_breast": 22}
    for key, expected in known.items():
        try:
            got = numbering_offset(load_target(key))
        except Exception:
            continue
        check(f"{key} offset is {expected}", got is not None and got[0] == expected,
              f"got {got}")


def test_kernel_scale_and_sign() -> None:
    """Two things the kernel must not read: surface length, and nothing.

    Both were found by controls whose answer was known in advance, and both
    had been shipping.

    LENGTH. Writing the same epitope out twice changes no chemistry at all --
    same residues, same proportions, same everything per residue. The kernel
    moved 0.475 -> 0.725 -> 0.508 across one, two and three copies, and
    0.475 -> 0.275 when the face was doubled instead. Epitopes here run from
    11 to 34 residues and 28 of 70 panel surrogates are a different length
    from the surface they are compared against, so two fifths of the
    cross-reactivity comparisons were partly measuring annotation length.

    SIGN. The margin's zero is a real boundary: above it the face prefers the
    target, below it the face prefers the relative. The score was a ramp from
    the 10th percentile of random faces up to the ceiling, so a design that
    preferred the off-target still read well above half -- the optimised
    TNF-alpha design reached -0.153 and was reported as 65.6% selective.
    """
    section("kernel scale and sign")
    from plab.scoring.selectivity import (COUNTER_SELECTIVE_BAND,
                                          complementarity_kernel,
                                          margin_to_score)

    face, epi = "DDYYKKSSSSSSS", "PQEKGQEAPYDYDFAES"
    base = complementarity_kernel(face, epi)
    for n in (2, 3, 5):
        check(f"epitope written out {n} times scores the same",
              abs(complementarity_kernel(face, epi * n) - base) < 1e-9,
              f"{complementarity_kernel(face, epi * n):.3f} vs {base:.3f}")
        check(f"face written out {n} times scores the same",
              abs(complementarity_kernel(face * n, epi) - base) < 1e-9,
              f"{complementarity_kernel(face * n, epi):.3f} vs {base:.3f}")

    # A real pair from the catalogue: the TNF-alpha epitope against the
    # surrogate mapped onto TNFSF13, which is 14 residues to the target's 17.
    # Padding the shorter one to equal length must not move the comparison.
    short = "DVEPDAYRALNLSP"
    check("a surrogate of unequal length is compared on chemistry alone",
          abs(complementarity_kernel(face, short)
              - complementarity_kernel(face, short * 2)) < 1e-9,
          "length still leaking into the paralog comparison")

    lo, hi = -0.50, 0.35
    check("a negative margin never leaves the bottom quarter",
          all(margin_to_score(m, lo, hi) <= COUNTER_SELECTIVE_BAND + 1e-9
              for m in (-0.5, -0.3, -0.153, -0.01, -1e-6)),
          "a design that prefers the off-target must not read as selective")
    check("zero margin sits exactly at the band edge",
          abs(margin_to_score(0.0, lo, hi) - COUNTER_SELECTIVE_BAND) < 1e-9,
          f"got {margin_to_score(0.0, lo, hi):.3f}")
    check("the reachable ceiling is 100%",
          abs(margin_to_score(hi, lo, hi) - 1.0) < 1e-9,
          f"got {margin_to_score(hi, lo, hi):.3f}")
    ladder = [margin_to_score(m / 100, lo, hi) for m in range(-50, 36)]
    check("the ramp is monotonic across the sign change",
          all(b >= a - 1e-12 for a, b in zip(ladder, ladder[1:])),
          "a split ramp must not introduce a step backwards")
    # The regression that motivated the band, with its real numbers.
    check("the TNF-alpha design's -0.153 no longer reads as 65.6%",
          margin_to_score(-0.153, -0.50, 0.086) < 0.25,
          f"got {margin_to_score(-0.153, -0.50, 0.086) * 100:.1f}%")


def test_ceiling_noise() -> None:
    """A ceiling smaller than its own error bar cannot carry a percentage.

    `INDISTINGUISHABLE_MARGIN` was an absolute 0.02, and VEGF-A's ceiling is
    0.029 -- so it passed, and the VEGF-A design was reported as 100%
    selective off a reachable margin of 0.029 against CLEC12A's 0.350. The
    cutoff was the right idea in the wrong units: drop one residue from one
    mapped surrogate, which at 25-40% identity is a certainty rather than a
    hypothetical, and the ceiling moves by 0.003 for alpha-synuclein and
    0.077 for TROP2. Nothing fixed can span that.
    """
    section("ceiling noise")
    from plab.scoring.selectivity import (CEILING_NOISE_FACTOR,
                                          achievable_margin_range,
                                          build_paralog_epitopes,
                                          ceiling_noise)
    from plab.targets import load_off_targets, load_target

    t = load_target("clec12a_aml")
    panel = build_paralog_epitopes(t, load_off_targets(t))
    epi = t.epitope_sequence
    hi = achievable_margin_range(epi, panel, 13)[1]
    noise = ceiling_noise(epi, panel, 13)
    check("the ceiling has a measurable error bar", noise > 0,
          "a jackknife that returns zero is not measuring anything")
    check("CLEC12A's ceiling stands clear of its error bar",
          hi > CEILING_NOISE_FACTOR * noise,
          f"ceiling {hi:.3f}, noise {noise:.3f}")

    # A panel whose surrogate IS the target surface: the ceiling is zero, so
    # the error bar must be too, and the pair must not slip through.
    twin = ceiling_noise(epi, {"TWIN": epi}, 13)
    twin_hi = achievable_margin_range(epi, {"TWIN": epi}, 13)[1]
    check("an identical surrogate leaves no reachable margin",
          twin_hi <= CEILING_NOISE_FACTOR * twin + 0.02,
          f"ceiling {twin_hi:.3f}, noise {twin:.3f}")

    # Every target carrying a stored ceiling must carry its error bar too,
    # or the scorer silently falls back to the absolute cutoff that missed
    # VEGF-A in the first place.
    from plab.targets import available_targets
    missing = []
    for key in sorted(available_targets()):
        spec = load_target(key)
        if spec.panel_ceiling and not spec.panel_ceiling_noise:
            missing.append(key)
    check("a stored ceiling always comes with its error bar", not missing,
          f"missing for {missing}")


def test_panels() -> None:
    """A panel that produces a number is not the same as a useful panel.

    These are the controls, not the feature: a surrogate surface identical to
    the target's makes the margin structurally zero, and a metric that cannot
    vary across faces is the bug this project has already shipped three
    times in other components.
    """
    section("cross-reactivity panels")
    import statistics

    from plab.aminoacids import AA20
    from plab.scoring.selectivity import (
        build_paralog_epitopes, complementarity_kernel,
    )
    from plab.targets import available_targets, load_panels

    panels = load_panels()
    if not panels:
        check("panels built", False, "run scripts/build_panels.py")
        return

    check("a panel was attempted for every verified target",
          set(panels) >= set(available_targets()),
          f"missing {sorted(set(available_targets()) - set(panels))[:4]}")
    check("every panel records why it is the size it is",
          all(p.get("state") and (p.get("members") or p.get("note"))
              for p in panels.values()))

    sc = get_scaffold("affibody_Z")
    length = len(sc.variable_positions)
    rng = random.Random(3)
    alphabet = [a for a in AA20 if a != "C"]
    faces = ["".join(rng.choice(alphabet) for _ in range(length))
             for _ in range(60)]

    with_panel, identical, flat, degenerate = [], [], [], []
    for key in sorted(available_targets()):
        t = load_target(key)
        surrogates = build_paralog_epitopes(t, load_off_targets(t))
        if not surrogates:
            continue
        with_panel.append(key)
        if any(s == t.epitope_sequence for s in surrogates.values()):
            identical.append(key)
        margins = []
        for f in faces:
            on = complementarity_kernel(f, t.epitope_sequence)
            off = max(complementarity_kernel(f, s) for s in surrogates.values())
            margins.append(on - off)
        if statistics.pstdev(margins) < 1e-9:
            flat.append(key)
        if all(abs(m) < 1e-9 for m in margins):
            degenerate.append(key)

    check("a useful number of targets have a panel", len(with_panel) >= 12,
          f"only {len(with_panel)}")
    # The invariant that actually matters is not how many panels exist, but
    # that every target without one records which of several very different
    # reasons applies.
    reasons = {load_target(k).panel_state for k in available_targets()
               if not load_target(k).off_target_accessions}
    check("every target without a panel records why", "none" not in reasons,
          f"states: {sorted(reasons)}")
    check("no surrogate surface is identical to the target's",
          not identical, f"identical for {identical}")
    check("the margin varies across faces for every panel",
          not flat, f"constant for {flat}")
    check("no panel makes the margin identically zero",
          not degenerate, f"degenerate for {degenerate}")


def test_catalog() -> None:
    """Every verified target must load and score, not just the first one.

    The pipeline was built around CLEC12A and then generalised. This checks
    the generalisation really holds: each catalogue entry has to produce a
    real epitope, a surface description, a reference molecule and a working
    evaluator, with no leftover assumptions about arginines or antibodies.
    """
    section("catalogue generalisation")
    from plab.catalog import load_verified
    from plab.targets import available_targets

    cat = load_verified()
    targets = available_targets()
    check("catalogue has been built", len(targets) >= 5,
          f"only {len(targets)} verified; run scripts/build_catalog.py")
    check("several disease areas covered",
          len({r["area"] for r in targets.values()}) >= 4,
          f"areas: {sorted({r['area'] for r in targets.values()})}")
    check("rejected candidates are recorded with reasons",
          all(r.get("reason") for r in cat.get("rejected", [])))

    sc = get_scaffold("affibody_Z")
    problems = []
    for key in sorted(targets):
        try:
            t = load_target(key)
            if not t.epitope:
                problems.append(f"{key}: empty epitope")
                continue
            if not t.surface.get("required_character"):
                problems.append(f"{key}: no surface description")
            ev = Evaluator(t, use_plm=False)
            card = ev.evaluate(sc.framework, sc)
            if not (0 <= card.score_100 <= 100):
                problems.append(f"{key}: score out of range {card.score_100}")
            if any(c.normalised != c.normalised for c in card.components):
                problems.append(f"{key}: NaN in a component")
            if ev.reference is None:
                problems.append(f"{key}: no reference molecule")
        except Exception as exc:
            problems.append(f"{key}: {type(exc).__name__}: {exc}")

    check("every verified target loads and scores", not problems,
          "; ".join(problems[:3]))

    # A target with no panel AT ALL must say it was not measured, and must
    # say why, rather than scoring well by default.
    bare = [k for k in targets if not load_target(k).off_target_accessions]
    if bare:
        t = load_target(bare[0])
        ev = Evaluator(t, use_plm=False)
        comp = ev.evaluate(sc.framework, sc).component("paralog_selectivity")
        check("missing selectivity panel is reported, not assumed good",
              comp is not None and "not measured" in comp.detail.lower(),
              f"detail was: {comp.detail[:70] if comp else 'none'}")
        check("and the reason is specific, not generic",
              comp is not None and len(comp.detail) > 90,
              f"detail was: {comp.detail[:70] if comp else 'none'}")

    # A surface with no cations must not crash the cation-driven scorers.
    uncharged = [k for k in targets
                 if not load_target(k).cation_positions]
    if uncharged:
        t = load_target(uncharged[0])
        ev = Evaluator(t, use_plm=False)
        card = ev.evaluate(sc.framework, sc)
        check(f"uncharged surface ({t.gene}) scores without error",
              0 <= card.score_100 <= 100)
    else:
        print("  skip  no uncharged epitope in the catalogue to test")


def test_flat_components() -> None:
    """Constant components must not carry weight in the ranking.

    Measured across 500 sequences, size_suitability and footprint_match had a
    standard deviation of exactly zero -- on a fixed scaffold every candidate
    is the same length with the same face size. Four components together
    carried about 27% of the composite's weight while contributing a near
    constant bonus, which inflated every score and compressed the usable
    range. They are still measured and displayed; they are just not ranked.
    """
    import random as _random
    import statistics as _stats

    section("constant-component detection")
    ev = Evaluator(load_target(), use_plm=False)
    sc = get_scaffold("affibody_Z")

    flat = ev.flat_components(sc)
    check("constant components detected", len(flat) >= 2, f"got {sorted(flat)}")
    check("size_suitability is detected as constant", "size_suitability" in flat)
    check("footprint_match is detected as constant", "footprint_match" in flat)
    check("a discriminating component is NOT excluded",
          "paralog_selectivity" not in flat and "charge_complementarity" not in flat,
          f"wrongly excluded from {sorted(flat)}")

    card = ev.evaluate(sc.framework, sc)
    unranked = {c.name for c in card.unranked}
    check("card marks them unranked", unranked == set(flat),
          f"{sorted(unranked)} vs {sorted(flat)}")
    check("unranked components still carry a detail string",
          all(c.detail for c in card.unranked))
    check("unranked components explain why",
          all("not" in c.note.lower() for c in card.unranked))
    check("unranked contribute nothing to the score",
          all(c.contribution == 0.0 for c in card.unranked))
    check("ranking weight excludes them",
          abs(card.total_weight
              - sum(c.weight for c in card.components if c.ranked)) < 1e-9)

    # The point of the exclusion: a wider, more useful dynamic range.
    rng = _random.Random(21)
    ss = sc.secondary_structure
    ranked_scores, all_scores = [], []
    for _ in range(80):
        assign = {}
        for p in sc.variable_positions:
            pool = [a for a in aa.AA20 if a != "C"
                    and not (ss[p] == "H" and a in "PG")]
            assign[p] = rng.choice(pool)
        c = ev.evaluate(sc.apply(sc.framework, assign), sc)
        ranked_scores.append(c.score_100)
        tw = sum(x.weight for x in c.components if x.weight > 0)
        all_scores.append(
            100 * sum(x.normalised * x.weight
                      for x in c.components if x.weight > 0) / tw)
    check("excluding constants widens the score range",
          _stats.stdev(ranked_scores) > _stats.stdev(all_scores),
          f"sd {_stats.stdev(ranked_scores):.2f} vs "
          f"{_stats.stdev(all_scores):.2f} with them included")


def test_calibration() -> None:
    section("calibration")
    ev = Evaluator(load_target(), use_plm=False)
    sc = get_scaffold("affibody_Z")
    cal = ev.calibrate(sc, n_decoys=120, seed=99)

    check("decoy distribution has spread", cal.recognition.stdev > 2.0,
          f"sd {cal.recognition.stdev:.2f}")
    # The whole method rests on this: a real binder must stand out from noise.
    check("reference stands above random decoys",
          cal.reference_z_recognition is not None
          and cal.reference_z_recognition > 2.0,
          f"z {cal.reference_z_recognition}")
    check("discriminating power is reported",
          "strong" in cal.discriminating_power
          or "moderate" in cal.discriminating_power,
          cal.discriminating_power)
    check("percentile of the reference is near the top",
          cal.recognition.percentile_of(cal.reference_recognition) > 95)


def test_loop() -> None:
    section("design loop")
    ev = Evaluator(load_target(), use_plm=False)
    sc = get_scaffold("affibody_Z")
    loop = DesignLoop(ev, sc, LoopConfig(population_size=20, generations=8,
                                         patience=99, seed=3))
    loop.run(verbose=False)

    first = loop.history[0].best.score
    last = loop.best().score
    check("search improves on its starting point", last > first,
          f"{first:.4f} -> {last:.4f}")
    check("history recorded for every generation", len(loop.history) == 9,
          f"got {len(loop.history)}")
    # Monotone in the RANK KEY, which is what the loop maximises. The
    # composite alone is not monotone and must not be asserted to be: the
    # generation that first finds a candidate which is not disqualified
    # usually has a LOWER composite than the counter-selective one it
    # displaces, and that is the run working.
    check("best is monotone across history",
          all(loop.history[i].best.rank_key <= loop.best().rank_key
              for i in range(len(loop.history))))
    check("framework preserved in every candidate",
          all(all(ind.sequence[i] == sc.framework[i]
                  for i in range(sc.length)
                  if i not in sc.variable_positions)
              for g in loop.history for ind in g.individuals))
    check("no cysteines introduced",
          all("C" not in ind.sequence
              for g in loop.history for ind in g.individuals))
    check("no helix breakers at helical designed positions",
          all(ind.sequence[p] not in "PG"
              for g in loop.history for ind in g.individuals
              for p in sc.variable_positions
              if sc.secondary_structure[p] == "H"))
    check("top() returns distinct sequences",
          len({i.sequence for i in loop.top(10)}) == len(loop.top(10)))
    check("reproducible with the same seed",
          DesignLoop(ev, sc, LoopConfig(population_size=20, generations=8,
                                        patience=99, seed=3)
                     ).run(verbose=False).best.sequence == loop.best_generation().best.sequence)


def test_mmcif() -> None:
    """The same structure must read the same from either file format.

    Not an optional extra. The legacy PDB format cannot hold an entry past
    62 chains or 99,999 atoms, and the PDB no longer generates .pdb files
    for a growing share of ordinary entries either. Two targets in the
    second research pass -- the IL-31 receptor with nemolizumab (9EZA) and
    the glucagon receptor (9XN4) -- have no .pdb file at all, and were
    rejected as "could not fetch PDB", which reads like a missing structure
    rather than a reader that cannot open one. It matters more for the
    build-your-own tab, where a person typing in a gene lands on recent
    entries constantly and has no way to know the reason is a file format.
    """
    section("mmCIF reader")
    from plab.pdbio import parse_mmcif, parse_structure, to_pdb_text

    both = [p.stem for p in STRUCT.glob("*.pdb")
            if (STRUCT / f"{p.stem}.cif").exists()]
    if not both:
        check("a structure in both formats is available", False,
              "run scripts/fetch_data.py")
    for pdb_id in sorted(both)[:2]:
        a = parse_structure(STRUCT / f"{pdb_id}.pdb")
        b = parse_structure(STRUCT / f"{pdb_id}.cif")
        ids_a = sorted({r.chain for r in a.residues})
        chains_a = {c: a.sequence(c) for c in ids_a}
        ids_b = sorted({r.chain for r in b.residues})
        chains_b = {c: b.sequence(c) for c in ids_b}
        check(f"{pdb_id}: both formats give the same chains",
              set(chains_a) == set(chains_b),
              f"{sorted(chains_a)} vs {sorted(chains_b)}")
        same = [c for c in chains_a if chains_a[c] == chains_b.get(c)]
        check(f"{pdb_id}: every chain sequence matches",
              len(same) == len(chains_a),
              f"{len(same)}/{len(chains_a)} identical")

    cif = STRUCT / "9EZA.cif"
    if cif.exists():
        s = parse_mmcif(cif)
        check("a cif-only entry parses", len(s.residues) > 500,
              f"{len(s.residues)} residues")
        check("and carries its resolution",
              s.resolution is not None and 1.0 < s.resolution < 5.0,
              f"got {s.resolution}")
        # The viewer speaks PDB, so a cif-only entry is re-emitted. Round
        # trip it rather than trusting the writer.
        import tempfile
        from pathlib import Path as _P
        with tempfile.TemporaryDirectory() as tmp:
            out = _P(tmp) / "round.pdb"
            out.write_text(to_pdb_text(s), encoding="utf-8")
            back = parse_structure(out)
        check("re-emitted as PDB it reads back the same",
              len(back.residues) == len(s.residues),
              f"{len(back.residues)} vs {len(s.residues)}")


def test_paralog_chain_selection() -> None:
    """A binding site must be read off the target, not off its neighbour.

    6Z6V holds all three C1q subunits. The nanobody binds subunit B, and
    C1QC matches its own chains at 100% and subunit B at 51% -- both above
    the aligned-identity threshold, so both were candidates. Picking the
    largest interface read C1QC's binding site off a different gene product,
    and nothing downstream could have noticed: the sequence is real, the
    contacts are real, and the epitope simply belonged to another protein.

    The fix must not go too far the other way. A first attempt banded
    identity finely and demoted chain Z of the PD-1 structure and chain D of
    the LAG-3 one, which are 90% and 95% identical because crystallographic
    copies are partly disordered, not because they are different proteins.
    Both lost their larger, better-ordered interface to a second copy.
    """
    section("paralog chain selection")
    from plab.epitope import SAME_PROTEIN_IDENTITY, derive_epitope

    p = STRUCT / "6Z6V.pdb"
    seq_c = SEQ / "P02747.fasta"
    if not p.exists() or not seq_c.exists():
        check("6Z6V and C1QC are available", False,
              "run scripts/build_catalog.py")
        return
    ref = "".join(seq_c.read_text(encoding="utf-8").splitlines()[1:])
    # The catalogue entry's own configuration: the nanobody is the binder of
    # interest, and it sits on subunit B. Without that preference the search
    # settles on a C1q subunit contact instead and the paralog question never
    # arises.
    r = derive_epitope(str(p), ref, prefer_binder="antibody")
    check("the site is read off a chain that IS the target",
          r.target_identity >= SAME_PROTEIN_IDENTITY,
          f"chain {r.target_chain} at {r.target_identity:.0%} identity")
    check("and the larger interface on the relative is reported, not hidden",
          any("identical to the reference" in w for w in r.warnings),
          f"warnings: {r.warnings}")

    # The copies case: both PD-1 chains are the same protein, so the larger
    # interface must still win.
    from plab.targets import available_targets
    cat = available_targets()
    for key in ("pd1_checkpoint", "lag3_checkpoint"):
        rec = cat.get(key)
        if not rec:
            continue
        check(f"{key} still uses its largest interface",
              rec["n_epitope_residues"] >= 15,
              f"{rec['n_epitope_residues']} residues on chain "
              f"{rec.get('target_chain')}")


def test_custom_targets() -> None:
    """The build-your-own path, without the network.

    What is checked here is the part that decides what a person is allowed
    to build: an accession has to look like an accession, and a structure is
    only offered if something substantial is bound to the protein. 3BXM
    pairs glutamate carboxypeptidase with a three-residue dipeptide, which
    would derive an "epitope" of almost nothing.
    """
    section("custom targets")
    from plab.custom import ACCESSION
    from plab.lookup import MIN_PARTNER_RESIDUES, looks_like_antibody

    for good in ("P01375", "Q5QGZ9", "O15123", "A6NIZ1"):
        check(f"{good} is recognised as an accession",
              bool(ACCESSION.match(good)))
    for bad in ("TNF", "KRAS", "SOD1", "CLEC12A", ""):
        check(f"{bad or '(empty)'} is not mistaken for one",
              not ACCESSION.match(bad),
              "gene symbols must go to the symbol lookup instead")

    check("Fab, nanobody and scFv all read as antibodies",
          all(looks_like_antibody(t) for t in
              ("tanezumab Fab heavy chain", "Nanobody C1q75",
               "Nemolizumab scFv", "IgG-binding protein")))
    check("a receptor does not", not looks_like_antibody("Glucagon receptor"))
    check("a tripeptide ligand is too small to be a binding partner",
          MIN_PARTNER_RESIDUES > 3,
          "or 3BXM's N-acetyl-aspartyl-glutamate would become an epitope")

    # Scan depth, with the network stubbed out. Entries come back sorted by
    # resolution, so a protein crystallised many times on its own buries its
    # complexes behind all of them. At the old limit of 40, superoxide
    # dismutase reported NO usable structure when it has eleven, and a
    # person typing SOD1 into the build-your-own tab would have been told
    # their protein could not be worked on.
    from plab import lookup as _lk
    buried_at = 97
    fake_ids = [f"X{i:03d}" for i in range(140)]

    def fake_structures_for(acc, limit=0):
        return fake_ids[:limit]

    def fake_entity_map(ids):
        out = []
        for i in ids:
            n = int(i[1:])
            ents = [{"chains": ["A"], "length": 150, "description": "target",
                     "accessions": ["P00441"], "organism": None,
                     "entity": i}]
            if n == buried_at:
                ents.append({"chains": ["B"], "length": 130,
                             "description": "nanobody 2", "accessions": [],
                             "organism": None, "entity": i})
            out.append({"pdb_id": i, "title": "t", "method": "X-ray",
                        "resolution": 2.0, "entities": ents})
        return out

    real_sf, real_em = _lk.structures_for, _lk.entity_map
    _lk.structures_for, _lk.entity_map = fake_structures_for, fake_entity_map
    try:
        found = _lk.rank_structures("P00441")
        shallow = _lk.rank_structures("P00441", limit=40)
    finally:
        _lk.structures_for, _lk.entity_map = real_sf, real_em
    check("a complex buried behind 97 apo structures is still found",
          len(found) == 1 and found[0]["pdb_id"] == f"X{buried_at:03d}",
          f"got {[r['pdb_id'] for r in found]}")
    check("and a shallow scan is what used to miss it",
          not shallow,
          "the regression this depth exists to prevent")


def test_panel_cache() -> None:
    """The cached surrogate surfaces must equal the derived ones, exactly.

    scripts/build_panels.py maps the epitope onto each panel member and
    stores the result; everything downstream used to redo those alignments
    once per Evaluator, which was most of a twenty-minute test suite and a
    second on every design run. Reading the stored answer instead is only
    safe while it is the same answer.

    This is here because the first version of that shortcut checked nothing.
    Its guard compared the filtered result against its own intersection --
    true by construction -- so CLEC12A's panel silently fell from nine
    members to eight, and its reachable ceiling moved from 0.350 to 0.359
    with no error and no warning. A cache whose correctness is not tested is
    a way of being confidently wrong faster.
    """
    section("panel cache")
    from plab.scoring.selectivity import build_paralog_epitopes
    from plab.targets import available_targets

    keys = [k for k in sorted(available_targets())
            if load_target(k).panel_surfaces][:12]
    check("some targets carry a verified cache", bool(keys),
          "run scripts/verify_panel_cache.py")

    mismatched, unmarked = [], []
    for key in keys:
        t = load_target(key)
        offs = load_off_targets(t)
        cached = build_paralog_epitopes(t, offs)
        # The same call with the cache taken away.
        bare = copy.copy(t)
        bare.panel_surfaces = {}
        bare.panel_epitope = ""
        derived = build_paralog_epitopes(bare, offs)
        if cached != derived:
            mismatched.append(
                f"{key}: cached {sorted(cached)} vs derived {sorted(derived)}")
        if not t.panel_epitope:
            unmarked.append(key)

    check("the cache returns exactly what deriving returns", not mismatched,
          "; ".join(mismatched[:2]))
    check("every panel member survives the cache, not just the stored ones",
          not mismatched,
          "a missing curated member is the failure that got through before")
    check("an unverified panel is not used from cache", not unmarked,
          f"unmarked: {unmarked[:3]}")

    # And the guard itself: a target whose recorded epitope no longer matches
    # must fall back rather than answer from a cache about a different site.
    t = load_target(keys[0])
    moved = copy.copy(t)
    moved.panel_epitope = "THISISNOTTHEEPITOPE"
    check("a moved binding site invalidates the cache",
          build_paralog_epitopes(moved, load_off_targets(t))
          == build_paralog_epitopes(t, load_off_targets(t)),
          "falling back must still give the right answer, just slower")


def test_disqualification() -> None:
    """Can a candidate aimed at the wrong protein still come first?

    It could, and it did. A weighted mean treats every axis as a trade-off,
    so nine good scores outvote one that says the molecule does the opposite
    of its job. On TNF-alpha the search found 25 candidates whose face
    complements TNF more than its nearest relative, and ranked a -0.218 one
    above all of them: 82.6 composite against 72.6, bought with
    manufacturability and solubility.

    Raising selectivity's weight would not have fixed it. The axis can only
    swing about 3.7 composite points between a counter-selective design and a
    decent one, and the other nine can differ by more than that. The problem
    is not that disqualification is underweighted; it is that
    disqualification is not a weight at all.
    """
    section("disqualification")
    from plab.scoring.base import ScoreCard, ScoreComponent

    def card(score: float, dq: str = "") -> ScoreCard:
        c = ScoreCard("X")
        c.add(ScoreComponent("a", 0.0, score, 1.0, "", disqualifying=dq))
        return c

    good, bad = card(0.60), card(0.95, "prefers the off-target")
    check("a disqualified candidate ranks below a lower-scoring clean one",
          good.rank_key > bad.rank_key,
          f"{good.rank_key} vs {bad.rank_key}")
    check("the composite itself is untouched",
          abs(bad.score - 0.95) < 1e-9,
          "the number must stay an honest mean of what was measured")
    check("disqualification names the component responsible",
          [c.name for c in bad.disqualified] == ["a"])
    # The fallback. If nothing qualifies, ranking must not collapse or throw;
    # it degenerates to the composite, which is the old behaviour, and the
    # run reports that no candidate passed.
    allbad = sorted([card(0.3, "x"), card(0.8, "x"), card(0.5, "x")],
                    key=lambda c: c.rank_key, reverse=True)
    check("with nothing qualified, order falls back to the composite",
          [round(c.score, 1) for c in allbad] == [0.8, 0.5, 0.3],
          f"got {[round(c.score, 1) for c in allbad]}")

    # End to end on the target that exposed it. The invariant is not a
    # number but a guarantee: nothing the dashboard offers may be
    # counter-selective while a clean candidate exists.
    ev = Evaluator(load_target("tnf_ra"), use_plm=False)
    sc = get_scaffold("affibody_Z")
    loop = DesignLoop(ev, sc, LoopConfig(population_size=40, generations=14,
                                         patience=8, seed=7))
    loop.run(verbose=False)
    everyone = {i.sequence: i for g in loop.history for i in g.individuals}
    clean = [i for i in everyone.values() if not i.card.disqualified]
    check("the search finds candidates that are not disqualified",
          len(clean) > 0, "nothing to rank, so nothing to check")
    if clean:
        comp = loop.best().card.component("paralog_selectivity")
        check("the returned best is not counter-selective",
              not loop.best().card.disqualified and comp.raw > 0,
              f"margin {comp.raw:+.3f}")
        check("no disqualified candidate appears in the top 12",
              not any(i.card.disqualified for i in loop.top(12)),
              f"{sum(1 for i in loop.top(12) if i.card.disqualified)} of 12")
        # And the thing that used to win is still in the population, still
        # scoring well on everything else -- it is ordered differently, not
        # suppressed.
        by_composite = max(everyone.values(), key=lambda i: i.card.score)
        check("the old winner is still measured and still reported",
              by_composite.card.score >= loop.best().card.score
              or not by_composite.card.disqualified,
              "disqualification must reorder, never delete")


def test_plm_scale() -> None:
    """Regression guard for a bug that shipped once.

    The plausibility scale was anchored to the 50C1 antibody (perplexity
    2.91). Measured against that, the wild-type Z domain -- a real natural
    protein at 10.08 -- scored 0%, and so did every candidate. A saturated
    component is worse than no component: it carries weight while supplying
    no gradient, and the search drifted to a sequence less natural-looking
    than a random library member.

    The rule this enforces: a real natural protein must never be rated
    implausible, and random sequences must rate worse than natural ones.
    """
    section("language model scale (slow: loads ESM-2)")
    import random as _random

    from plab.scoring.plm import get_runner

    runner = get_runner()
    if not runner.load():
        print(f"  skip  ESM-2 unavailable ({runner.status})")
        return

    ev = Evaluator(load_target(), use_plm=True)
    sc = get_scaffold("affibody_Z")

    anchors = ev.plm_anchors(sc)
    check("scaffold anchors measured", anchors is not None)
    if not anchors:
        return
    good, poor = anchors
    check("natural parent beats random library members", good < poor,
          f"parent {good:.2f} vs random {poor:.2f}")

    card = ev.evaluate(sc.framework, sc)
    comp = card.component("plm_plausibility")
    check("component is active", comp is not None and comp.weight > 0)
    check("the natural parent is NOT rated implausible",
          comp is not None and comp.normalised > 0.8,
          f"got {comp.normalised:.3f} -- this is the bug that shipped")

    # A shuffled framework keeps composition but destroys the fold signal, so
    # it must score worse than the real sequence.
    chars = list(sc.framework)
    _random.Random(7).shuffle(chars)
    shuffled = ev.evaluate("".join(chars), sc).component("plm_plausibility")
    check("shuffled sequence scores worse than the real one",
          shuffled is not None and comp is not None
          and shuffled.normalised < comp.normalised,
          f"shuffled {shuffled.normalised:.3f} vs real {comp.normalised:.3f}")

    # The scale must leave room to discriminate: if every candidate saturates
    # at either end, the component contributes no gradient to the search.
    rng = _random.Random(5)
    ss = sc.secondary_structure
    scores = []
    for _ in range(6):
        assign = {}
        for p in sc.variable_positions:
            pool = [a for a in aa.AA20 if a != "C"
                    and not (ss[p] == "H" and a in "PG")]
            assign[p] = rng.choice(pool)
        c = ev.evaluate(sc.apply(sc.framework, assign), sc)
        scores.append(c.component("plm_plausibility").normalised)
    check("random members do not all saturate at one end",
          len({round(s, 2) for s in scores}) > 1,
          f"got {[round(s, 3) for s in scores]}")


def test_lineage() -> None:
    """The evidence asymmetry, which is the whole basis of the feature.

    Every check here has an answer settled before the code runs, and each one
    guards a failure that would look like an ordinary list of results rather
    than like a bug.
    """
    print("\nlineage evidence")
    from plab.lineage import (CURATED, assess_gene, assess_tissues,
                              band_for_ntpm, is_human, load_expression)

    exp = load_expression()
    check("expression data present", len(exp.values) > 50,
          f"only {len(exp.values)} genes")

    # The rule: transcript abundance can raise a concern, never clear one.
    # There is deliberately no nTPM value that returns "absent".
    check("no nTPM yields absent",
          all(band_for_ntpm(v) != "absent"
              for v in (0.0, 0.1, 0.9, 1.0, 9.9, 10.0, 49.9, 50.0, 5000.0)))
    check("below detection yields unknown", band_for_ntpm(0.4) == "unknown")
    check("detected yields a graded level",
          (band_for_ntpm(1.0), band_for_ntpm(10.0), band_for_ntpm(50.0))
          == ("low", "moderate", "high"))

    # GPRC5D is the case that settles the rule: the screen flags it in skin
    # and reads 0.0 in tongue, while the clinic sees oral toxicity in 63% of
    # patients. A lineage with no transcript must therefore not be cleared.
    calls = {c.lineage: c for c in assess_gene("GPRC5D", exp)}
    neut = calls.get("neutrophil")
    check("GPRC5D neutrophil not cleared",
          neut is not None and not neut.is_cleared and neut.level == "unknown",
          f"got {neut.level if neut else None}")

    flags = {c.label: c for c in assess_tissues("GPRC5D", exp)}
    check("GPRC5D flagged outside the blood",
          "skin" in flags and flags["skin"].level in ("moderate", "high"),
          f"skin={flags['skin'].level if 'skin' in flags else 'absent'}")
    check("GPRC5D carries its curated keratinised-tissue entry",
          "keratinised structures" in flags)

    # Only cytometry or observed toxicity may clear a population, and the
    # curated table is the only route to one.
    for (gene, lin), cur in CURATED.items():
        if cur.level == "absent":
            check(f"{gene}/{lin} absence has evidence behind it",
                  cur.evidence in ("flow", "clinical") and bool(cur.source),
                  f"evidence={cur.evidence!r} source={cur.source!r}")

    cleared = [c for c in assess_gene("CLEC12A", exp) if c.is_cleared]
    check("CLEC12A cleared on stem cells",
          any(c.lineage == "hsc" for c in cleared),
          f"cleared={[c.lineage for c in cleared]}")

    from plab.lineage import HARM
    check("every curated harm value is a known one",
          all(c.harm in HARM for c in CURATED.values()),
          str([k for k, c in CURATED.items() if c.harm not in HARM]))

    # GPRC5D has no stem-cell row, deliberately. Nothing found tested it on
    # normal haematopoietic stem cells -- the stem-cell exclusion that is
    # easy to misattribute here belongs to BCMA -- so the cell must stay
    # unknown rather than borrow a neighbour's evidence.
    check("GPRC5D has no stem-cell claim",
          ("GPRC5D", "hsc") not in CURATED)
    check("BCMA does have one", ("TNFRSF17", "hsc") in CURATED)
    gp = {c.lineage: c for c in assess_gene("GPRC5D", exp)}
    check("GPRC5D stem cells read unknown", gp["hsc"].level == "unknown",
          gp["hsc"].level)

    # A viral target has no human lineage to spare, and that is a statement
    # about the target rather than about our ignorance.
    check("human target reads as human", is_human(load_target("clec12a_aml")))
    check("viral target does not", not is_human(load_target("sars2_rbd")))

    # --- The solid tumour axis -------------------------------------------
    #
    # Blood is not where these antigens do their damage, so these checks are
    # about the tissue table and the two new ways the screen was caught
    # misleading.
    from plab.lineage import (CURATED_TISSUE, EVIDENCE_RANK, LineageCall,
                              TISSUE_SCREEN, decisive, phrase)

    # The exact column code, which is the trap that already bit this file
    # once: `t_RNA_skin` silently returned a row with no skin in it. A
    # misspelt stomach column would read as "CLDN18 is not in the stomach",
    # which is the single worst sentence this module could produce.
    check("stomach is screened under its real column code",
          ("t_RNA_stomach_1", "stomach") in TISSUE_SCREEN)
    check("the stomach column actually returns data",
          exp.for_gene("CLDN18").get("t_RNA_stomach_1", 0) > 100,
          str(exp.for_gene("CLDN18").get("t_RNA_stomach_1")))

    # CLDN18: the screen's loud flag is the lung isoform zolbetuximab does
    # not bind, and the stomach -- where the emesis comes from -- reads
    # higher still. Both must be true for the entry to mean what it says.
    cl = exp.for_gene("CLDN18")
    check("CLDN18 reads higher in stomach than in lung",
          cl["t_RNA_stomach_1"] > cl["t_RNA_lung"] > 100,
          f"stomach={cl['t_RNA_stomach_1']} lung={cl['t_RNA_lung']}")

    ct = {c.label: c for c in assess_tissues("CLDN18", exp)}
    check("the curated stomach row supersedes the transcript row",
          ct["stomach"].evidence == "clinical",
          ct["stomach"].evidence)
    check("and there is only one stomach row",
          sum(1 for c in assess_tissues("CLDN18", exp)
              if c.label == "stomach") == 1)
    check("the gastric stem cell zone is positively excluded",
          ct["gastric stem cell zone"].is_cleared
          and ct["gastric stem cell zone"].evidence == "ihc",
          ct["gastric stem cell zone"].evidence)

    # A tissue the screen has no column for must say so rather than letting
    # the silence read as agreement. There is no mesothelium column.
    ms = {c.label: c for c in assess_tissues("MSLN", exp)}
    meso = next(c for c in ms.values()
                if c.lineage.startswith("curated_tissue"))
    check("mesothelium has no screen column and says so",
          meso.ntpm is None and "no column" in meso.detail,
          meso.detail[-90:])

    # The duodenum column was added for PSMA and immediately became its
    # loudest tissue. If this stops being true the column is miscoded.
    fh = exp.for_gene("FOLH1")
    check("PSMA reads loudest in duodenum",
          fh["t_RNA_duodenum"] > fh["t_RNA_salivary_gland"] > 10,
          f"duodenum={fh['t_RNA_duodenum']} salivary={fh['t_RNA_salivary_gland']}")

    # Animal evidence may soften a flag and may never clear a human
    # population. TROP2 is the entry that needed the tier; this asserts the
    # threshold rather than that one entry.
    check("animal evidence is not decisive", not decisive("animal"))
    check("human tissue and cytometry are",
          decisive("ihc") and decisive("flow") and decisive("clinical"))
    check("transcript is not", not decisive("rna"))
    check("an absence on animal evidence cannot clear",
          not LineageCall("x", "x", "absent", "animal").is_cleared)
    check("the same absence on cytometry can",
          LineageCall("x", "x", "absent", "flow").is_cleared)
    check("clinical outranks every other kind",
          EVIDENCE_RANK["clinical"] == max(EVIDENCE_RANK.values()))

    for (gene, site), cur in CURATED_TISSUE.items():
        check(f"{gene}/{site} is sourced",
              bool(cur.source) and cur.harm in HARM,
              f"source={cur.source!r} harm={cur.harm!r}")
        if cur.level == "absent":
            check(f"{gene}/{site} absence has decisive evidence behind it",
                  decisive(cur.evidence), cur.evidence)

    # KIT is the one solid tumour target whose liability IS a blood
    # population, and the evidence is a drug that depletes it on purpose.
    check("KIT is carried on stem cells, with harm observed",
          CURATED[("KIT", "hsc")].level == "high"
          and CURATED[("KIT", "hsc")].harm == "observed")

    # The GPRC5D/BCMA discipline applied again: the granulocyte expression
    # often attached to "CEA" belongs to CEACAM8 and CEACAM3, so CEACAM5
    # must not borrow it.
    check("CEACAM5 claims no granulocyte row",
          ("CEACAM5", "neutrophil") not in CURATED)
    ce = {c.lineage: c for c in assess_gene("CEACAM5", exp)}
    check("so its neutrophil cell reads unknown",
          ce["neutrophil"].level == "unknown", ce["neutrophil"].level)

    # Checkpoint receptors: expression on T cells is the mechanism, and the
    # "survives" verdict belongs to the Fc that was actually tested.
    check("PD-1 on T cells is not read as destruction",
          CURATED[("PDCD1", "t_cell")].harm == "not_observed")
    check("CTLA-4's harm is left unknown, not guessed",
          CURATED[("CTLA4", "t_cell")].harm == "unknown")

    check("labels read correctly mid-sentence",
          (phrase("T cells"), phrase("NK cells"),
           phrase("Monocytes and macrophages"))
          == ("T cells", "NK cells", "monocytes and macrophages"))

    # --- Topology: whether the lineage question applies at all -----------
    #
    # The failure these guard was live and silent. The screen reported
    # alpha-synuclein at 103.2 nTPM on dendritic cells, in the same words
    # and pills as CD33 on progenitors, for a protein that sits in the
    # cytosol where no antibody reaches.
    from plab.lineage import (CURATED_TOPOLOGY, TOPOLOGIES, classify_topology,
                              load_topology, off_axis)

    topo = load_topology()
    check("topology data present", len(topo.entries) > 50,
          f"only {len(topo.entries)} entries")
    check("no accession was silently substituted", not topo.unresolved,
          str(topo.unresolved))

    # The rule, on synthetic evidence, so these assert the rule rather than
    # today's annotation.
    def ev(spans=0, **kw):
        return {"membrane_spans": spans, "keywords": kw.get("keywords", [])}

    check("a membrane span means surface",
          classify_topology(ev(1, keywords=["Transmembrane"])) == "surface")
    check("a GPI anchor means surface too",
          classify_topology(ev(0, keywords=["GPI-anchor"])) == "surface")
    check("span plus secreted means shed",
          classify_topology(ev(1, keywords=["Transmembrane", "Secreted"]))
          == "shed")
    check("secreted alone means secreted",
          classify_topology(ev(0, keywords=["Secreted"])) == "secreted")
    check("neither means intracellular",
          classify_topology(ev(0, keywords=["Cytoplasm"])) == "intracellular")

    # The trap that made the first version of this rule wrong. "Cell
    # membrane" is given to proteins that merely associate with a membrane,
    # and using it as the surface test pulled tau in beside CLEC12A.
    check("membrane ASSOCIATION is not membrane insertion",
          classify_topology(ev(0, keywords=["Cell membrane", "Cytoplasm"]))
          == "intracellular")

    # Controls on the real annotation. Each has an answer settled before the
    # code runs, so a future UniProt change that breaks one fails loudly
    # instead of quietly reclassifying a target.
    for gene, want in (("CLEC12A", "surface"), ("CD33", "surface"),
                       ("CEACAM5", "surface"), ("ERBB2", "surface"),
                       ("TNF", "shed"), ("MSLN", "shed"), ("CD274", "shed"),
                       ("IL6", "secreted"), ("VEGFA", "secreted"),
                       ("PCSK9", "secreted"), ("C5", "secreted"),
                       ("SOD1", "intracellular"), ("SNCA", "intracellular"),
                       ("MAPT", "intracellular")):
        got, why = topo.verdict(gene)
        check(f"{gene} is {want}", got == want, f"got {got}: {why}")

    check("every verdict is a known one",
          all(topo.verdict(g)[0] in TOPOLOGIES for g in topo.entries))
    check("every verdict carries its reasoning",
          all(len(topo.verdict(g)[1]) > 20 for g in topo.entries))
    check("each curated override says why",
          all(v in TOPOLOGIES and len(why) > 40
              for v, why in CURATED_TOPOLOGY.values()),
          str(CURATED_TOPOLOGY))

    # The off-axis row must be inert in every bucket. It is a statement that
    # the question does not apply, and if it could ever read as "positively
    # cleared" it would earn rank credit for not having been asked.
    row = off_axis("IL6", "secreted", "because")
    check("the off-axis row clears nothing", not row.is_cleared)
    check("and is not a concern either", not row.is_concern)

    # --- What is lost when the signal is blocked -------------------------
    #
    # The other half of the sentence the topology axis left unfinished.
    from plab.lineage import (PHYSIOLOGY, SEVERITIES, SEVERITY_LABEL,
                              load_topology as _lt, physiology_for)

    check("every physiology entry is well formed",
          all(p.severity in SEVERITIES and p.evidence in EVIDENCE_RANK
              and len(p.signal) > 30 and len(p.cost) > 60 and bool(p.source)
              for p in PHYSIOLOGY.values()),
          str([g for g, p in PHYSIOLOGY.items()
               if p.severity not in SEVERITIES or not p.source]))
    check("every severity has a label",
          all(s_ in SEVERITY_LABEL for s_ in SEVERITIES))

    # A table where everything is dangerous carries no information. The
    # controls are what make the boxed entries mean something.
    tolerated = [g for g, p in PHYSIOLOGY.items()
                 if p.severity == "well_tolerated"]
    boxed = [g for g, p in PHYSIOLOGY.items()
             if p.severity == "boxed_warning"]
    check("the table records targets that are safely blocked",
          len(tolerated) >= 2, str(tolerated))
    check("as well as ones regulators boxed", len(boxed) >= 2, str(boxed))
    check("PCSK9 is one of the controls",
          physiology_for("PCSK9").severity == "well_tolerated")

    # The misattribution guard, same discipline as CEACAM5's missing
    # granulocyte row, and now with the answer it was protecting. Blocking
    # the alternative or terminal pathway carries a boxed warning for
    # encapsulated bacteria; blocking early classical complement points at
    # autoimmunity instead, because that is what the inborn errors show.
    # Had the class warning been extended across the family in the first
    # pass, this distinction would have been erased before it was found.
    for gene in ("C3", "C5", "CFB", "CFD"):
        check(f"{gene} carries the encapsulated-bacteria box",
              physiology_for(gene).severity == "boxed_warning")
    for gene in ("C1S", "C1QB", "MASP2"):
        check(f"{gene} did not inherit that box",
              physiology_for(gene).severity != "boxed_warning",
              physiology_for(gene).severity)
    check("the early classical entries point at autoimmunity",
          all("lupus" in physiology_for(g).cost.lower()
              or "autoimmun" in physiology_for(g).cost.lower()
              for g in ("C1QB", "C1S", "MASP2")))

    # "contested" has to be reachable and distinct from "unknown": one means
    # the published answer disagrees with itself, the other means nobody
    # looked, and collapsing them is the mistake this package is built on.
    check("contested is a severity in its own right",
          "contested" in SEVERITIES and "unknown" in SEVERITIES)
    check("and something uses it",
          any(p.severity == "contested" for p in PHYSIOLOGY.values()))

    # Lookup is case-insensitive and absent genes return None rather than
    # raising -- most of the catalogue has no entry and that is the point.
    check("lookup is case-insensitive",
          physiology_for("tnf") is physiology_for("TNF"))
    check("an absent gene is None, not an error",
          physiology_for("CLEC12A") is None)

    # Coverage, asserted as a floor rather than an exact number so adding an
    # entry does not break the test that an entry was added.
    _topo = _lt()
    off = [g for g in _topo.entries
           if _topo.verdict(g)[0] in ("secreted", "intracellular")]
    covered = [g for g in off if g in PHYSIOLOGY]
    check("most off-axis targets now carry what blocking them costs",
          len(covered) >= len(off) - 5,
          f"{len(covered)} of {len(off)}")

    # Several entries are "safe to block, and the drug did not work".
    # Factor XI, myostatin, tau and alpha-synuclein all land there, and the
    # pattern is worth protecting: a signal that can be removed without
    # consequence may be a signal that was not doing much.
    quiet_failures = [g for g in ("F11", "MSTN", "MAPT", "SNCA")
                      if physiology_for(g)
                      and physiology_for(g).severity == "well_tolerated"]
    check("the safe-but-ineffective group is recorded as such",
          len(quiet_failures) == 4, str(quiet_failures))


def test_specimen() -> None:
    """Specimen ranking, including the ordering rule that defines it."""
    print("\nspecimen ranking")
    import copy

    from plab.lineage import load_expression
    from plab.scoring.selectivity import numbering_offset
    from plab.specimen import Specimen, Variant, assess_target, rank

    exp = load_expression()
    clec = load_target("clec12a_aml")
    cd33 = load_target("cd33_aml")
    marrow = ("hsc", "progenitor", "erythroid", "platelet")

    # The motivating case. CD33 and CLEC12A share no domain family, so no
    # paralog panel can separate them; only expression on healthy cells can,
    # and getting this backwards is the clinical error the feature exists
    # to prevent.
    spec = Specimen(antigens={"CLEC12A": "positive", "CD33": "positive"},
                    spared=marrow)
    a = assess_target(clec, spec, exp)
    b = assess_target(cd33, spec, exp)
    check("CD33 disqualified when marrow must recover", b.disqualified,
          b.disqualifying)
    check("CLEC12A not disqualified", not a.disqualified, a.disqualifying)
    check("CLEC12A outranks CD33", a.rank_key > b.rank_key)

    # The ordering rule: being examined counts in a target's favour, so an
    # unexamined target must not float to the top for having no findings.
    #
    # The comparator is found rather than named. An earlier version of this
    # test used CD38 as the unexamined case and broke the moment CD38's
    # stem-cell row was curated -- which is the test asserting a fact about
    # today's table instead of the rule it exists to defend.
    from plab.lineage import CURATED
    from plab.catalog import load_verified
    blank = next(
        (k for k, rec in sorted(load_verified()["targets"].items())
         if (rec["gene"].upper(), "hsc") not in CURATED
         and rec["area"] in ("Blood cancer", "Immuno-oncology")), None)
    check("a target with no stem-cell cytometry exists to compare against",
          blank is not None)

    examined = assess_target(clec, Specimen(
        antigens={"CLEC12A": "positive"}, spared=("hsc",)), exp)
    other = load_target(blank)
    unexamined = assess_target(other, Specimen(
        antigens={other.gene.upper(): "positive"}, spared=("hsc",)), exp)
    check("examined target has a cleared lineage", len(examined.cleared) == 1)
    check("unexamined target has an unknown one",
          len(unexamined.unknowns) == 1,
          f"{other.gene}: {[c.lineage for c in unexamined.unknowns]}")
    check("cleared outranks unknown",
          examined.rank_key > unexamined.rank_key)

    # An antigen the blasts do not carry is a non-candidate, not a weak one.
    gone = assess_target(clec, Specimen(antigens={"CLEC12A": "negative"}), exp)
    check("antigen-negative disqualifies", gone.disqualified)
    check("not_tested does not disqualify",
          not assess_target(
              clec, Specimen(antigens={"CLEC12A": "not_tested"}),
              exp).disqualified)

    # Variants. The first lands on a real epitope residue; the second states
    # a reference amino acid that is not there, which means the two are not
    # using the same numbering and no answer would mean anything.
    off, _ = numbering_offset(clec)
    r = clec.epitope[0]
    hit = assess_target(clec, Specimen(
        antigens={"CLEC12A": "positive"},
        variants=(Variant("CLEC12A", r.resseq + off, r.aa, "A"),)), exp)
    check("epitope variant flagged", hit.variant_hits[0].hits_epitope)
    check("epitope variant disqualifies", hit.disqualified)

    wrong = "W" if r.aa != "W" else "A"
    bad = assess_target(clec, Specimen(
        antigens={"CLEC12A": "positive"},
        variants=(Variant("CLEC12A", r.resseq + off, wrong, "A"),)), exp)
    check("numbering mismatch refused",
          bad.variant_hits[0].reference_mismatch)
    check("numbering mismatch draws no conclusion", not bad.disqualified)

    # Disqualified entries are kept and sorted last, never dropped: the
    # reason a target was excluded is the most useful thing produced here.
    ordered = rank([clec, cd33], spec, exp)
    check("ranking keeps disqualified entries", len(ordered) == 2)
    check("disqualified sorts last", ordered[-1].gene == "CD33")

    # Expression on a spared population is not the same as destroying it.
    # CD47 sits on haematopoietic stem cells and blocking it still does not
    # increase phagocytosis of normal stem and progenitor cells, so the
    # stem-cell row must not exclude it. The platelet row, where nothing
    # reassures, still does -- which is what keeps this from being a
    # blanket exemption.
    cd47 = assess_target(load_target("cd47_macrophage"), Specimen(
        antigens={"CD47": "positive"}, spared=("hsc",)), exp)
    stem = [c for c in cd47.conflicts if c.lineage == "hsc"]
    check("CD47 is recorded as carried on stem cells",
          len(stem) == 1 and stem[0].level == "high")
    check("CD47 not excluded on the stem-cell row", not cd47.disqualified,
          cd47.disqualifying)
    check("the exemption is explained, not silent",
          any("not destroyed" in n for n in cd47.notes), str(cd47.notes))
    plat = assess_target(load_target("cd47_macrophage"), Specimen(
        antigens={"CD47": "positive"}, spared=("platelet",)), exp)
    check("CD47 still excluded where nothing reassures", plat.disqualified)

    # TIM-3 is absent from normal stem cells, which is the cleanest window
    # this axis can describe -- and its phase 3 trial still failed.
    tim3 = assess_target(load_target("tim3_mds"), Specimen(
        antigens={"HAVCR2": "positive"}, spared=marrow), exp)
    check("TIM-3 cleared on stem cells",
          any(c.lineage == "hsc" for c in tim3.cleared))
    check("TIM-3 not excluded for marrow recovery", not tim3.disqualified)

    # CLEC12A spares the stem compartment, not mature myeloid cells. The
    # tool has to say so rather than let it read as selective everywhere.
    neutro = assess_target(clec, Specimen(
        antigens={"CLEC12A": "positive"},
        spared=("neutrophil", "monocyte")), exp)
    check("CLEC12A excluded when mature myeloid cells must survive",
          neutro.disqualified, "CLEC12A is on granulocytes and monocytes")

    # A non-human target is neither cleared nor unknown on human lineages.
    viral = assess_target(load_target("sars2_rbd"), Specimen(
        antigens={"S": "positive"}, spared=marrow), exp)
    check("viral target marked not applicable",
          viral.inapplicable is not None)
    check("viral target earns no cleared-lineage credit",
          not viral.cleared and not viral.unknowns)

    # --- Solid tumour antigens, where the sparing axis is the wrong one ---

    # KIT is the exception that proves the rule: a GIST target whose
    # haematopoietic liability is the stem compartment itself, established by
    # an anti-CD117 antibody developed AS transplant conditioning. It must be
    # excluded by a marrow-recovery profile exactly as CD33 is.
    kit = assess_target(load_target("kit_gist"), Specimen(
        antigens={"KIT": "positive"}, spared=marrow), exp)
    check("KIT excluded when marrow must recover", kit.disqualified,
          kit.disqualifying)
    # ...and its card must not then claim the liability is usually elsewhere.
    check("KIT is not told its liability sits in tissue",
          not any("sits in tissue" in n for n in kit.notes), str(kit.notes))

    # HER2's liability is cardiac, and no flow panel on this specimen can
    # reach it. The silent failure would be an empty conflict list reading as
    # a clean target, so the card has to say which question went unanswered.
    her2 = assess_target(load_target("her2_breast"), Specimen(
        antigens={"ERBB2": "positive"}, spared=marrow), exp)
    check("HER2 carries no haematopoietic conflict", not her2.conflicts)
    check("HER2 is not thereby clean",
          len(her2.tissue_harms) == 1
          and her2.tissue_harms[0].label == "heart muscle",
          str([c.label for c in her2.tissue_harms]))
    check("and the card says the sparing profile cannot reach it",
          any("outside blood" in n and "flow panel" in n for n in her2.notes),
          str(her2.notes))
    check("HER2's cardiac row is clinical, not a transcript flag",
          her2.tissue_harms[0].evidence == "clinical")

    # An absence of haematopoietic findings must not outrank an antigen
    # somebody actually examined on the populations in question.
    check("CLEC12A still outranks HER2 on a marrow profile",
          a.rank_key > her2.rank_key,
          f"clec={a.rank_key} her2={her2.rank_key}")

    # Observed organ damage outranks a transcript flag, and is outranked by a
    # conflict with a population the clinician named. Both halves matter: the
    # first stops a clinical finding being a last-place tiebreak, the second
    # stops it overriding the constraint actually stated.
    msln = assess_target(load_target("msln_mesothelioma"), Specimen(
        antigens={"MSLN": "positive"}, spared=marrow), exp)
    trop2 = assess_target(load_target("trop2_tnbc"), Specimen(
        antigens={"TACSTD2": "positive"}, spared=marrow), exp)
    check("mesothelin has an observed organ harm", len(msln.tissue_harms) == 1)
    check("TROP2 has none, on far more transcript flags",
          not trop2.tissue_harms and len(trop2.tissue_flags) > 5,
          f"flags={len(trop2.tissue_flags)}")
    check("one organ damaged in patients costs more than nine weak flags",
          trop2.rank_key > msln.rank_key,
          f"trop2={trop2.rank_key} msln={msln.rank_key}")

    # A target where neither table has an entry must not look examined. This
    # is the ranking rule from the top of the file applied to tissue.
    met = assess_target(load_target("met_nsclc"), Specimen(
        antigens={"MET": "positive"}, spared=marrow), exp)
    check("an uncurated solid tumour target says nobody looked",
          any("not that the tissue is clear" in n for n in met.notes),
          str(met.notes))

    # CLDN18.2's exclusion from the gastric stem cell zone is the best fact
    # in the tissue table and still earns no rank credit, because it is an
    # answer about epithelium to a question about blood.
    cldn = assess_target(load_target("cldn18_gastric"), Specimen(
        antigens={"CLDN18": "positive"}, spared=marrow), exp)
    check("the gastric stem-cell exclusion is carried", len(cldn.tissue_spared) == 1,
          str([c.label for c in cldn.tissue_spared]))
    stripped = copy.copy(cldn)
    stripped.tissue_spared = []
    check("and it changes no ranking", stripped.rank_key == cldn.rank_key)

    # --- Targets that are not on a cell surface at all -------------------
    #
    # The headline fix. Before this, SNCA's transcript reading on dendritic
    # cells came through the lineage loop as "carried on dendritic cells",
    # with the consequence line warning about lost antigen presentation --
    # for a cytosolic protein a binder never meets.
    myeloid = ("monocyte", "dendritic", "neutrophil")
    snca = assess_target(load_target("synuclein_parkinsons"), Specimen(
        antigens={"SNCA": "positive"}, spared=myeloid), exp)
    check("alpha-synuclein is read as intracellular",
          snca.topology == "intracellular", snca.topology)
    check("and no longer names populations at risk",
          not snca.conflicts and not snca.unknowns,
          f"{len(snca.conflicts)} conflicts, {len(snca.unknowns)} unknowns")
    check("it says the question does not apply instead",
          snca.off_axis is not None
          and "does not go" in snca.off_axis.detail)
    check("the off-axis row earns no cleared credit", not snca.cleared)

    il6 = assess_target(load_target("il6_cytokine"), Specimen(
        antigens={"IL6": "positive"}, spared=myeloid), exp)
    check("a secreted cytokine is off this axis too",
          il6.topology == "secreted" and il6.off_axis is not None)
    check("and the card says the transcript measures production",
          "PRODUCE" in il6.off_axis.detail)

    # TNF is the case that would have been got wrong by lumping cytokines
    # together: it is made as a transmembrane protein, monocytes really do
    # display it, and anti-TNF antibodies really do bind that form. Its
    # lineage rows must survive.
    tnf = assess_target(load_target("tnf_ra"), Specimen(
        antigens={"TNF": "positive"}, spared=myeloid), exp)
    check("TNF is shed, not simply secreted", tnf.topology == "shed")
    check("so its lineage rows are kept", tnf.conflicts and not tnf.off_axis,
          f"{len(tnf.conflicts)} conflicts")
    check("and the soluble pool is called out as a sink",
          any("soluble pool" in n for n in tnf.notes), str(tnf.notes))

    # A target being off-axis must not cost it its tissue rows: VEGFA's are
    # the most informative thing on its card.
    vegfa = assess_target(load_target("vegfa_angio"), Specimen(
        antigens={"VEGFA": "positive"}, spared=myeloid), exp)
    check("an off-axis target keeps its tissue rows",
          bool(vegfa.tissue_harms) and vegfa.off_axis is not None,
          f"{len(vegfa.tissue_harms)} tissue harms")

    # The off-axis card used to stop at "the risk moves to whatever normal
    # physiology depends on that signal" and say nothing further. Now it
    # either carries the answer or says that half is open.
    c5 = assess_target(load_target("c5_pnh"), Specimen(
        antigens={"C5": "positive"}, spared=myeloid), exp)
    check("a secreted target carries what blocking it costs",
          c5.physiology is not None
          and c5.physiology.severity == "boxed_warning")
    check("and the off-axis row points at it",
          "recorded separately" in c5.off_axis.detail)

    # IL11, IL13 and SOD1 are deliberately left without an entry: searching
    # turned up nothing solid, and "the search found nothing" is a fact
    # about the search rather than about the world. Writing an entry on the
    # strength of it would be nTPM-0.0-means-absent one layer up.
    il11 = assess_target(load_target("il11_fibrosis"), Specimen(
        antigens={"IL11": "positive"}, spared=myeloid), exp)
    check("a target with no entry says that half is open",
          il11.physiology is None
          and "is not recorded here" in il11.off_axis.detail,
          il11.off_axis.detail[-120:])

    # TNF is shed, so it is NOT off-axis -- and it still needs the
    # physiology row, because tuberculosis reactivation is the thing a
    # clinician most needs to see about this target and no lineage row
    # would ever show it.
    check("a shed target gets the physiology row too",
          tnf.physiology is not None and tnf.off_axis is None)
    check("and it is the tuberculosis one",
          "tuberculosis" in tnf.physiology.cost)

    # It must not touch the ranking: the clinician asked which populations
    # must survive, not what the molecule is for.
    before = c5.rank_key
    c5.physiology = None
    check("physiology earns no rank credit either way",
          c5.rank_key == before)


def main() -> int:
    slow = "--plm" in sys.argv
    tests = [test_tables, test_charge, test_descriptors, test_parser,
             test_align, test_scaffolds, test_scoring, test_selectivity,
             test_chain_matching, test_binder_grouping,
             test_numbering_offset, test_kernel_scale_and_sign,
             test_ceiling_noise, test_panels,
             test_catalog, test_flat_components, test_calibration,
             test_loop, test_panel_cache, test_disqualification,
             test_mmcif,
             test_paralog_chain_selection, test_custom_targets,
             test_lineage, test_specimen]
    if slow:
        tests.append(test_plm_scale)
    else:
        print("(pass --plm to also run the ESM-2 scale tests)")

    for fn in tests:
        try:
            fn()
        except Exception:
            _failed.append(fn.__name__)
            print(f"  ERROR in {fn.__name__}:")
            traceback.print_exc()

    print(f"\n{'=' * 60}")
    print(f"{_passed} passed, {len(_failed)} failed")
    if _failed:
        for f in _failed:
            print(f"  - {f}")
        return 1
    print("All tests passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
