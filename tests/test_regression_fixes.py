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
