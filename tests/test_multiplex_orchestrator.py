"""Tests for dye/label resolution in the multiplex orchestrator.

Guards the wiring that makes the M11 dye-channel check meaningful: dyes must come
from each target's own config chain (target config -> primerlab.template ->
parameters.probe.labels.fluorophore), and an unresolvable dye must degrade to a
sentinel that the channel check reports as UNKNOWN — never to a silent FAM/HEX
default that would fabricate a "separated" verdict.
"""
import textwrap

import pytest

from conftest import load_module  # noqa: E402

batch = load_module("scripts/run_multiplex_batch.py", "run_multiplex_batch")
iupac = load_module("scripts/tools/iupac_utils.py", "iupac_utils")


def _write_target(tmp_path, name, dye=None, template_ref=True, build="build-x"):
    tpl = tmp_path / f"tpl_{name}.yaml"
    if dye is not None:
        tpl.write_text(textwrap.dedent(f"""
            parameters:
              probe:
                labels:
                  fluorophore: "{dye}"
        """))
    cfg = tmp_path / f"cfg_{name}.yaml"
    body = f'build_name: "{build}"\n'
    if template_ref:
        body += f'primerlab:\n  template: "{tpl}"\n'
    cfg.write_text(body)
    return cfg


def test_resolves_dye_from_primerlab_template(tmp_path):
    cfg = _write_target(tmp_path, "a", dye="HEX", build="flu-b")
    build, dye, label = batch.parse_target_meta(str(cfg))
    assert build == "flu-b" and dye == "HEX" and label == "flu-b"


def test_explicit_dye_overrides_template(tmp_path):
    cfg = _write_target(tmp_path, "b", dye="HEX", build="t")
    cfg.write_text(cfg.read_text() + '  dye: "Cy5"\n')
    _, dye, _ = batch.parse_target_meta(str(cfg))
    assert dye == "Cy5"


def test_missing_template_yields_unresolved_sentinel(tmp_path):
    cfg = _write_target(tmp_path, "c", template_ref=False, build="t")
    _, dye, _ = batch.parse_target_meta(str(cfg))
    assert dye == batch.UNRESOLVED_DYE


def test_unresolved_dye_is_flagged_not_silently_passed(tmp_path):
    # The whole point: an unknown dye must NOT read as a successful separation.
    status, ok = iupac.check_dye_separation("FAM", batch.UNRESOLVED_DYE)
    assert not ok and "UNKNOWN" in status


def test_target_label_is_used_when_provided(tmp_path):
    cfg = _write_target(tmp_path, "d", dye="FAM", build="b1")
    cfg.write_text(cfg.read_text() + 'target_label: "Influenza A MP"\n')
    _, _, label = batch.parse_target_meta(str(cfg))
    assert label == "Influenza A MP"


@pytest.mark.parametrize("cfg_path,expected_dye", [
    ("config/config_influenza_a_mp.yaml", "FAM"),
    ("config/config_influenza_b_mp.yaml", "HEX"),
    ("config/config_measles.yaml", "Cy5"),
])
def test_real_repo_configs_resolve_expected_dyes(cfg_path, expected_dye):
    _, dye, _ = batch.parse_target_meta(cfg_path)
    assert dye == expected_dye


def test_real_flu_a_b_measles_are_spectrally_separated():
    dyes = [batch.parse_target_meta(p)[1] for p in (
        "config/config_influenza_a_mp.yaml",
        "config/config_influenza_b_mp.yaml",
        "config/config_measles.yaml",
    )]
    _, ok = iupac.check_dye_separation(*dyes)
    assert ok


def test_export_filter_accepts_build_names_as_active_targets():
    es = load_module("scripts/export_supplementary.py", "export_supplementary")
    es.ACTIVE_TARGETS = ["influenza-a-mp-indo-only", "influenza-b-mp-indo-only"]
    assert es.is_multiplex_active(
        "multiplex-influenza-a-mp-indo-only_x_influenza-b-mp-indo-only") is True
    assert es.is_multiplex_active(
        "multiplex-triplex-influenza-a_influenza-b_measles") is False
