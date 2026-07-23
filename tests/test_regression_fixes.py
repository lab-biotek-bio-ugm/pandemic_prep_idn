"""Regression tests for defects found by the post-fix adversarial review.

These pin behaviour that the first round of fixes got WRONG — either by being so
strict it would break working pipelines, or by claiming to close a leak it did not
actually close.
"""
import json

import pytest

from conftest import load_module  # noqa: E402

iupac = load_module("scripts/tools/iupac_utils.py", "iupac_utils")
es = load_module("scripts/export_supplementary.py", "export_supplementary")
m08 = load_module("scripts/m08_degenerate_rescue.py", "m08_degenerate_rescue")


# --- B1: the triplex export leak the first fix did NOT close -----------------

def test_uninformative_dir_name_alone_cannot_identify_targets():
    es.ACTIVE_TARGETS = ["influenza-a", "influenza-b"]
    # The orchestrator writes every triplex here; the name carries no virus token.
    assert es._viruses_in("multiplex-triplex-all-3") == []


def test_triplex_leak_closed_via_result_content(tmp_path):
    """A Flu+Measles triplex under the orchestrator's literal directory name must
    still be excluded when only Influenza A/B are active."""
    es.ACTIVE_TARGETS = ["influenza-a-mp-indo-only", "influenza-b-mp-indo-only"]
    f = tmp_path / "triplex_candidates.json"
    f.write_text(json.dumps([{"triplex_id": "FluA_R9_FluB_R14_Measles_R1", "score": 88}]))
    assert es.is_multiplex_active("multiplex-triplex-all-3", f) is False


def test_flu_only_triplex_content_is_kept(tmp_path):
    es.ACTIVE_TARGETS = ["influenza-a", "influenza-b"]
    f = tmp_path / "multiplex_candidates.json"
    f.write_text(json.dumps([{"pair_id": "influenza-a_R1_influenza-b_R2"}]))
    assert es.is_multiplex_active("multiplex-triplex-all-3", f) is True


# --- B2: empty dye must not read as "separated" ------------------------------

@pytest.mark.parametrize("dyes", [("FAM", ""), ("", ""), ("FAM", None), ()])
def test_missing_dye_never_passes(dyes):
    _, ok = iupac.check_dye_separation(*dyes)
    assert ok is False


def test_valid_dyes_still_pass():
    assert iupac.check_dye_separation("FAM", "HEX", "Cy5")[1] is True


# --- A2: M08 hairpin must not be looser than the PrimerLab engine ------------

def test_hairpin_default_matches_engine_template():
    """The engine templates set hairpin_dg_min: -2.0. Defaulting the rescue step to
    the doc's -3.0 would make it LOOSER (more negative = more permissive), letting
    the degenerate track accept a hairpin the non-degenerate track rejects."""
    import argparse
    parser = argparse.ArgumentParser()
    # Mirror the script's own default rather than re-parsing the whole CLI.
    src = open("scripts/m08_degenerate_rescue.py").read()
    assert '"--hairpin-dg-min", type=float, default=-2.0' in src


def test_thermo_fail_sentinel_not_shown_as_real_dg():
    """A primer3 failure must not print as '-99.00 kcal/mol', which a bench
    scientist would read as catastrophic secondary structure."""
    assert "kcal/mol" not in m08._fmt_dg(m08._THERMO_FAIL_DG)
    assert "NA" in m08._fmt_dg(m08._THERMO_FAIL_DG)
    # Real values still render normally.
    assert m08._fmt_dg(-3.25) == "-3.25 kcal/mol"


# --- C: entropy rescaling consistency ---------------------------------------

def test_entropy_normalisers_agree_across_modules():
    """Every module that normalises M07 entropy must use log2(5), otherwise gappy
    columns produce values >1.0 (plot clipping, out-of-range GFF3 attributes).

    The consumers are checked by source inspection so this test does not require
    matplotlib (m07_plot_entropy) or a BLAST toolchain to be installed.
    """
    import math
    m07 = load_module("scripts/m07_bridge_builder.py", "m07_bridge_builder")
    assert m07.MAX_ENTROPY_BITS == pytest.approx(math.log2(5))

    plot_src = open("scripts/m07_plot_entropy.py").read()
    assert "H_MAX = math.log2(5)" in plot_src, "m07_plot_entropy still uses a 4-symbol H_MAX"

    m09_src = open("scripts/m09_blast_verification.py").read()
    assert "val / math.log2(5)" in m09_src, "m09 entropy normalisation still divides by 2.0"

    sweep_src = open("scripts/tools/m07_sensitivity_sweep.py").read()
    assert "math.log2(5)" in sweep_src, "sensitivity sweep still lacks entropy normalisation"


def test_entropy_thresholds_are_not_rescaled_by_the_normaliser():
    """The M07 entropy thresholds are compared against RAW bits, while
    MAX_ENTROPY_BITS only normalises the display score. A gap-free column has
    identical raw entropy under the 4- and 5-symbol alphabets, so multiplying the
    thresholds by log2(5)/log2(4) — as an earlier revision did — did not "rescale"
    anything, it just loosened every criterion by 16%."""
    import math
    import yaml
    cfg = yaml.safe_load(open("config/config.yaml"))
    block = cfg["bridge"]["scoring"]["multipliers"]["conserved_block"]["entropy_threshold"]
    assert cfg["bridge"]["entropy_threshold"] == pytest.approx(0.020, abs=1e-4)
    assert block == pytest.approx(0.10, abs=1e-4)
    assert cfg["primerlab"]["design"]["entropy_threshold"] == pytest.approx(0.50, abs=1e-4)

    # Pin the premise: a column with no gaps is unaffected by the alphabet change,
    # so no threshold applied to raw bits should have moved.
    def h(maf):
        return -((1 - maf) * math.log2(1 - maf) + maf * math.log2(maf))
    assert h(0.10) == pytest.approx(0.4690, abs=1e-3)

    # The comparison really is against raw bits, not a normalised value.
    src = open("scripts/m07_bridge_builder.py").read()
    assert "if mean_e < entropy_threshold:" in src
    assert "mean_e / MAX_ENTROPY_BITS" in src  # normaliser used only for the score


# --- M11 threshold unification ------------------------------------------------

def test_multiplex_thresholds_are_shared_not_duplicated():
    """Pairwise and triplex must apply the SAME ΔTm / cross-dimer limits. They
    previously held separate literals and had already drifted (3.5 vs 4.0°C), so
    an identical pair was PASS in one report and WARNING in the other."""
    assert iupac.MAX_DELTA_TM_C == 4.0
    pair_src = open("scripts/m11_multiplex_checker.py").read()
    tri_src = open("scripts/m11_triplex_checker.py").read()
    for src, name in ((pair_src, "pairwise"), (tri_src, "triplex")):
        assert "MAX_DELTA_TM_C" in src, f"{name} does not use the shared ΔTm threshold"
        assert "delta_tm <= 3.5" not in src, f"{name} still hardcodes 3.5"
        assert "delta_tm <= 4.0" not in src, f"{name} still hardcodes 4.0"


def test_compatibility_score_formula_is_shared():
    """The two checkers used different weights (pairwise ΔG×3.5/ΔTm×4.0, triplex
    ×4.0/×3.0), so scores on a shared 0-100 scale were not comparable between the
    two reports. Both must now call the one implementation."""
    for path, name in (("scripts/m11_multiplex_checker.py", "pairwise"),
                       ("scripts/m11_triplex_checker.py", "triplex")):
        src = open(path).read()
        assert "compatibility_score(" in src, f"{name} does not use the shared score"
        assert "* 3.5" not in src and "* 4.0" not in src and "* 3.0" not in src, \
            f"{name} still computes the score with its own weights"

    # ΔG is weighted above ΔTm: at 39°C strand invasion is recombinase-driven, so
    # dimer formation dominates over inter-primer Tm matching.
    assert iupac.SCORE_WEIGHT_CROSS_DIMER > iupac.SCORE_WEIGHT_DELTA_TM
    assert iupac.compatibility_score(0.0, 0.0) == 100.0
    assert iupac.compatibility_score(-5.0, 2.0) == 74.0
    assert iupac.compatibility_score(-1.0, 0.5, dye_ok=False) == 0.0
    # Positive ΔG (no structure) must not inflate the score above 100.
    assert iupac.compatibility_score(3.0, 0.0) == 100.0


def test_score_is_never_reported_as_a_percentage():
    """A 0-100 heuristic with no experimental calibration must not be rendered as
    '%', which reads as a probability of the multiplex working."""
    for path in ("scripts/m11_multiplex_checker.py", "scripts/m11_triplex_checker.py"):
        assert "['score']}%" not in open(path).read(), f"{path} still prints score as %"


# --- 3'-end anchoring ---------------------------------------------------------

def test_three_prime_anchor_detects_designed_complementarity():
    """Only an annealed 3' terminus can be extended, and those products amplify as
    primer noise — the factor the TwistAmp manual names as limiting sensitivity."""
    core = "ACGTGACTTAGCCTGAATCGGATTCACGGCATT"
    rc8 = iupac.reverse_complement(core[-8:])
    partner = "GCATTAGCTCAGATTGCAAGTCC" + rc8
    assert iupac.three_prime_anchor_length(core, partner) >= 8
    assert iupac.three_prime_anchor_length(core, partner) > iupac.MAX_THREE_PRIME_ANCHOR_BP


def test_three_prime_anchor_is_symmetric_and_safe_on_empty():
    a, b = "ACGTACGTACGTAAAA", "TTTTACGTACGTACGT"
    assert iupac.three_prime_anchor_length(a, b) == iupac.three_prime_anchor_length(b, a)
    assert iupac.three_prime_anchor_length("", b) == 0
    assert iupac.three_prime_anchor_length(a, None) == 0


def test_three_prime_dg_is_reported_but_not_thresholded():
    """A ΔG over a short terminal window saturates and cannot separate designed
    complementarity from background, so no PASS/WARNING decision may depend on it.
    """
    for path in ("scripts/m11_multiplex_checker.py", "scripts/m11_triplex_checker.py"):
        src = open(path).read()
        assert "max_anchor_3p <= MAX_THREE_PRIME_ANCHOR_BP" in src
        assert "min_dg_3p >=" not in src, f"{path} gates on the saturating ΔG metric"

    # Demonstrate the saturation the design note claims: doubling the length of a
    # perfect 3' complement does not move a 5 nt window.
    core = "ACGTGACTTAGCCTGAATCGGATTCACGGCGCC"
    p5 = "ATCAGTCAGTCAGTCA" + iupac.reverse_complement(core[-5:])
    p10 = "ATCAGTCAGTCAGTCA" + iupac.reverse_complement(core[-10:])
    dg5 = iupac.evaluate_three_prime_dimer(core, p5, 39.0, window=5)
    dg10 = iupac.evaluate_three_prime_dimer(core, p10, 39.0, window=5)
    assert abs(dg10 - dg5) < 1.0, "5 nt window unexpectedly resolved the difference"
    # The anchor length does resolve it.
    assert (iupac.three_prime_anchor_length(core, p10)
            > iupac.three_prime_anchor_length(core, p5))


def test_thf_cleavage_window_claim_removed_from_docs():
    """The documented "mismatches within 2 bp of the THF site are prohibited" rule
    had no basis. The TwistAmp manual §3.1.1 states the opposite — a mismatch on
    one of the dT-label positions, which sit 0-2 bases from the THF, is an
    explicitly sanctioned design fallback."""
    doc = open("docs/pipeline-architecture/m10-insilico-validation.md").read()
    # The rule must not be stated as active guidance...
    assert "**Cleavage Window**" not in doc
    assert "within **2 bp**" not in doc
    # ...and the retraction must be explicit, since a silent deletion would leave
    # anyone who read the old page still believing it.
    assert "there is no 2 bp cleavage window" in doc
    assert "5–9 base pair mismatches" in doc, "the contradicting source is not cited"


# --- M10: a probe must not be judged with primer 3' kinetics ------------------

def _m10():
    return load_module("scripts/m10_insilico_pcr.py", "m10_insilico_pcr")


def test_probe_3prime_mismatches_no_longer_reject_detection():
    """An RAA probe carries a 3' blocking group so the polymerase cannot extend
    it. Running it through assess_raa_binding() applied primer extension kinetics
    to an oligo that cannot be extended, rejecting probes for mismatches at a
    position where the chemistry gives them no effect — understating coverage for
    exactly the sequences with 3'-proximal probe variation."""
    m10 = _m10()
    probe = "GAATTTCAGAGGCTATAGCGATCTCAGGTCAATCGATAGATCGCTA"
    flip = {"A": "G", "G": "A", "C": "T", "T": "C"}

    def mutate(positions):
        s = list(probe)
        for i in positions:
            s[i] = flip[s[i]]
        return "ACGTACGTAC" + "".join(s) + "TTGCAATTGC"

    near_3p = mutate([len(probe) - 3, len(probe) - 6])
    # The old criterion rejected this; the mismatch load is only 2.
    assert m10.assess_raa_binding(probe, near_3p, critical_3prime_bp=8,
                                  max_3prime_mm=1, max_internal_mm=4)[0] is False
    ok, n_mm = m10.assess_probe_binding(probe, near_3p, max_mismatches=5)
    assert ok is True and n_mm == 2


def test_probe_still_rejected_above_the_documented_tolerance():
    """The TwistAmp manual's 5-9 range is a tolerance, not a licence: a probe
    carrying more mismatches than the threshold must still fail."""
    m10 = _m10()
    probe = "GAATTTCAGAGGCTATAGCGATCTCAGGTCAATCGATAGATCGCTA"
    flip = {"A": "G", "G": "A", "C": "T", "T": "C"}
    s = list(probe)
    for i in (3, 9, 15, 21, 28, 35, 42):
        s[i] = flip[s[i]]
    amp = "ACGTACGTAC" + "".join(s) + "TTGCAATTGC"
    ok, n_mm = m10.assess_probe_binding(probe, amp, max_mismatches=5)
    assert ok is False and n_mm == 7


def test_absent_probe_does_not_block_detection():
    """A candidate with no probe means detection == amplification, not failure."""
    m10 = _m10()
    assert m10.assess_probe_binding("", "ACGTACGT", max_mismatches=5) == (True, 0)


def test_probe_threshold_is_wired_end_to_end():
    """A threshold that exists only as a CLI default is not configurable."""
    import yaml
    cfg = yaml.safe_load(open("config/config.yaml"))
    assert cfg["insilico_pcr"]["raa_probe_max_mm"] == 5
    smk = open("rules/m10_insilico_pcr.smk").read()
    assert "raa_probe_max_mm=config" in smk
    assert "--raa-probe-max-mm {params.raa_probe_max_mm}" in smk
    src = open("scripts/m10_insilico_pcr.py").read()
    # No probe path may still call the primer routine.
    assert "assess_raa_binding(\n                            probe_deg" not in src


# --- Synthesis specification for RAA probes -----------------------------------

def test_probe_spec_is_not_a_taqman_layout():
    """Reports described every probe as "5'-FAM / 3'-BHQ1 / Internal-THF" — a
    hydrolysis-probe layout. An RAA exo probe carries both dyes on INTERNAL dT
    residues flanking the abasic site, with a blocking group at the 3' end. An
    oligo house given the old string would have synthesised the wrong molecule."""
    spec = iupac.probe_synthesis_spec(
        {"probe_seq": "ACGT", "probe_type": "exo"}, "FAM")
    assert "Internal [FAM-dT]" in spec["modification"]
    assert "blocked, non-extendable" in spec["modification"]
    assert "5'-FAM" not in spec["modification"]
    for path in ("scripts/m11_multiplex_checker.py", "scripts/m11_triplex_checker.py"):
        assert "Internal-THF" not in open(path).read(), f"{path} still emits the old label"


def test_probe_spec_uses_primerlab_markup_and_swaps_the_dye():
    labeled = "GAATTTCAGAGG[FAM-dT]C[THF]A[BHQ1-dT]CGATAGATCGCTA[C3-spacer]"
    spec = iupac.probe_synthesis_spec(
        {"probe_labeled_seq": labeled, "probe_thf_index": 30,
         "probe_type": "exo", "probe_compliant": True}, "HEX")
    assert "[HEX-dT]" in spec["order_sequence"]
    assert "[FAM-dT]" not in spec["order_sequence"]
    assert spec["thf_index"] == 30
    assert spec["notes"] == ""


def test_probe_spec_refuses_to_look_orderable_without_labelling():
    """A bare sequence with no abasic-site position is not an orderable probe;
    saying so is the difference between a warning and a wasted synthesis order."""
    spec = iupac.probe_synthesis_spec({"probe_seq": "ACGTACGT"}, "FAM")
    assert "do not order" in spec["notes"]
    assert "position unknown" in spec["notes"]


def test_probe_spec_reports_deliberate_label_mismatches():
    spec = iupac.probe_synthesis_spec(
        {"probe_labeled_seq": "AC[FAM-dT][THF][BHQ1-dT]GT", "probe_thf_index": 30,
         "probe_label_mismatches": 2, "probe_compliant": False}, "FAM")
    assert "non-T base" in spec["notes"]


def test_probe_labelling_is_forwarded_by_m09():
    """PrimerLab exports the fields, but M11 only sees what M09 puts in
    final_candidates.json."""
    src = open("scripts/m09_blast_verification.py").read()
    for field in ("probe_labeled_seq", "probe_thf_index", "probe_type",
                  "probe_label_mismatches", "probe_compliant"):
        assert f'"{field}":' in src, f"M09 does not forward {field}"


# --- Primer 3' anchor region --------------------------------------------------

def test_raa_3prime_window_matches_the_published_anchor_region():
    """Higgins et al. 2022 (J Mol Diagn 24:1207-1216) define the RPA primer anchor
    region as pre-ante-penultimate to 3'-terminal — 4 nt — and found lone
    mismatches outside the terminal position not to significantly affect
    time-to-positivity. The previous value of 8 had no source and charged
    positions 5-8 against the strict 3' quota, understating coverage."""
    import yaml
    cfg = yaml.safe_load(open("config/config.yaml"))
    assert cfg["insilico_pcr"]["raa_3prime_bp"] == 4
    raw = open("config/config.yaml").read()
    assert "jmoldx.2022.08.005" in raw, "the threshold is not traceable to its source"


# --- M05/M07: selection-significance thresholds -------------------------------

def test_meme_threshold_is_read_from_config_not_hardcoded():
    """Every config carried hyphy.p_value_threshold and nothing read it; M07
    counted MEME danger sites against a literal 0.05. A user tuning the config
    value got no effect at all."""
    src = open("scripts/m07_bridge_builder.py").read()
    assert 'config.get("hyphy", {}).get("p_value_threshold"' in src
    assert "meme_pvals if p < meme_pvalue_threshold" in src
    assert "meme_pvals if p < 0.05" not in src, "the hardcoded literal is still there"


def test_meme_threshold_follows_the_hyphy_recommendation():
    """HyPhy recommends p <= 0.1 for MEME rather than the conventional 0.05,
    because MEME is conservative on small alignments. FUBAR is a posterior
    probability, not a p-value, and keeps its own >= 0.9 convention."""
    import yaml
    m07 = load_module("scripts/m07_bridge_builder.py", "m07_bridge_builder")
    assert m07.MEME_PVALUE_THRESHOLD == 0.1
    for name in ("config.yaml", "config_influenza_a_mp.yaml", "config_measles.yaml"):
        cfg = yaml.safe_load(open(f"config/{name}"))
        assert cfg["hyphy"]["p_value_threshold"] == 0.1, name
        # FUBAR keeps a posterior-probability criterion, not a p-value.
        tiers = cfg.get("bridge", {}).get("tiers", {})
        if tiers.get("strict"):
            assert tiers["strict"].get("min_purifying_prob", 0.9) >= 0.9, name


def test_dead_bridge_p_value_threshold_removed():
    """bridge.p_value_threshold (0.01) was read by nothing and contradicted the
    hyphy value (0.1). Two unused numbers disagreeing about the same concept is
    worse than one, because either looks authoritative to a reader."""
    import glob
    import yaml
    for path in glob.glob("config/config*.yaml"):
        cfg = yaml.safe_load(open(path))
        assert "p_value_threshold" not in cfg.get("bridge", {}), path
