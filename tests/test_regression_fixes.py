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


def test_config_entropy_thresholds_rescaled():
    """Raw-bits thresholds in config must be rescaled to the 5-symbol alphabet,
    otherwise a 2% gap column silently breaks conserved blocks."""
    import yaml
    cfg = yaml.safe_load(open("config/config.yaml"))
    block = cfg["bridge"]["scoring"]["multipliers"]["conserved_block"]["entropy_threshold"]
    assert block == pytest.approx(0.116, abs=1e-3)
    assert cfg["bridge"]["entropy_threshold"] == pytest.approx(0.023, abs=1e-3)


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
