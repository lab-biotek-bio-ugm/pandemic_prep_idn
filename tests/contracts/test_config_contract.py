"""Cross-file agreement checks for config/config*.yaml.

The bug class this exists for is threshold drift between the seven parallel
pathogen configs, and dead keys that look like configuration but are read by
nothing. Both are invisible at run time: Snakefile line 6 declares
`configfile: "config/config.yaml"`, and Snakemake deep-merges any --configfile
on top of it rather than replacing it, so an omitted key silently takes the
influenza-A-MP value. Every check here therefore loads each file STANDALONE.

Only pyyaml is needed, so this module runs in the cheap CI job that installs no
scientific dependencies. Keep it that way.
"""
import re

import pytest
import yaml

from conftest import CONFIG_PATHS, REPO, flatten, load_config

pytestmark = pytest.mark.contract

CONTRACT = yaml.safe_load(
    (REPO / "tests" / "contracts" / "config_contract.yaml").read_text(encoding="utf-8")
)
_GAPS = CONTRACT["known_gaps"]
_SIMPLE_TYPES = {"str": str, "list": list, "dict": dict, "bool": bool}

_RATCHET_HINT = (
    "\n\nIf you have just fixed this, delete the entry from known_gaps in "
    "tests/contracts/config_contract.yaml. The lists there are only allowed to "
    "shrink."
)


def _known(section, name):
    entry = _GAPS.get(section, {}).get(name)
    return list(entry["keys"]) if entry else []


def _check_value(value, spec, path, config_name):
    wanted = spec["type"]
    if wanted == "int":
        # `type(v) is int` rather than isinstance: bool subclasses int, and this
        # also catches the YAML hazard where a value written on the line below
        # its key folds into a scalar and silently becomes a string.
        assert type(value) is int, (
            f"{config_name}: {path} is {type(value).__name__} {value!r}, expected int"
        )
    elif wanted == "float":
        assert type(value) in (int, float), (
            f"{config_name}: {path} is {type(value).__name__}, expected a number"
        )
    else:
        assert isinstance(value, _SIMPLE_TYPES[wanted]), (
            f"{config_name}: {path} is {type(value).__name__}, expected {wanted}"
        )

    if "min" in spec:
        assert value >= spec["min"], f"{config_name}: {path}={value} is below {spec['min']}"
    if "max" in spec:
        assert value <= spec["max"], f"{config_name}: {path}={value} is above {spec['max']}"
    if "one_of" in spec:
        assert value in spec["one_of"], (
            f"{config_name}: {path}={value!r} is not one of {spec['one_of']}"
        )
    if spec.get("non_empty"):
        assert str(value).strip(), f"{config_name}: {path} is empty"
    if spec.get("is_repo_path"):
        assert (REPO / str(value)).exists(), (
            f"{config_name}: {path} points at a file that does not exist: {value}"
        )


@pytest.mark.parametrize("cfg", CONFIG_PATHS, ids=lambda p: p.name)
def test_required_keys_are_present(cfg):
    flat = flatten(load_config(cfg))
    missing = sorted(p for p in CONTRACT["required"] if p not in flat)
    expected = sorted(_known("missing_required_keys", cfg.name))

    assert missing == expected, (
        f"{cfg.name}: the set of missing required keys is not the recorded one.\n"
        f"  missing now: {missing}\n"
        f"  recorded   : {expected}\n"
        "A key listed as required is read without a guard somewhere in rules/ or "
        "scripts/. It does not crash today only because config/config.yaml is "
        "merged underneath every --configfile, so the influenza value is used "
        "instead of yours." + _RATCHET_HINT
    )


@pytest.mark.parametrize("cfg", CONFIG_PATHS, ids=lambda p: p.name)
def test_required_keys_are_well_typed(cfg):
    flat = flatten(load_config(cfg))
    for path, spec in CONTRACT["required"].items():
        if path in flat:
            _check_value(flat[path], spec, path, cfg.name)


@pytest.mark.parametrize("key", CONTRACT["identical_across_configs"])
def test_pipeline_wide_thresholds_have_not_drifted(key):
    values = {p.name: flatten(load_config(p)).get(key, "<MISSING>") for p in CONFIG_PATHS}
    drifted = len({repr(v) for v in values.values()}) > 1
    known = key in _GAPS.get("drifted_thresholds", {})

    if known:
        assert drifted, (
            f"{key} no longer drifts, but is still listed under "
            f"known_gaps.drifted_thresholds." + _RATCHET_HINT
        )
        return

    assert not drifted, (
        f"{key} has drifted across the pathogen configs:\n"
        + "\n".join(f"  {name}: {value!r}" for name, value in sorted(values.items()))
        + "\nThis key is pipeline-wide policy. If it genuinely has to vary per "
        "pathogen, move it to may_vary in tests/contracts/config_contract.yaml "
        "with a reason."
    )


@pytest.mark.parametrize("cfg", CONFIG_PATHS, ids=lambda p: p.name)
def test_no_keys_that_nothing_reads(cfg):
    top = load_config(cfg)
    flat = flatten(top)
    present = sorted(
        bad
        for bad in CONTRACT["forbidden"]
        if bad in top or any(k == bad or k.startswith(bad + ".") for k in flat)
    )
    expected = sorted(_known("forbidden_keys_present", cfg.name))

    detail = "\n".join(
        f"  {bad}: "
        + (f"use {CONTRACT['forbidden'][bad]['alias_of']!r} instead. " if "alias_of" in CONTRACT["forbidden"][bad] else "")
        + " ".join(CONTRACT["forbidden"][bad]["reason"].split())
        for bad in present
        if bad not in expected
    )
    assert present == expected, (
        f"{cfg.name}: the set of keys read by nothing is not the recorded one.\n"
        f"  present now: {present}\n"
        f"  recorded   : {expected}\n" + (detail and "\n" + detail) + _RATCHET_HINT
    )


@pytest.mark.parametrize("cfg", CONFIG_PATHS, ids=lambda p: p.name)
def test_no_duplicate_keys(cfg):
    # load_config uses a loader that raises on a repeated mapping key. Plain
    # PyYAML keeps the last occurrence silently, which in these files would be
    # an invisible threshold override.
    load_config(cfg)


@pytest.mark.parametrize("cfg", CONFIG_PATHS, ids=lambda p: p.name)
def test_no_number_folded_into_a_string(cfg):
    """Guard for the "value written on the line below its key" YAML hazard.

    config/config_influenza_b_mp.yaml writes hyphy.threads in that shape today
    and happens to parse as int 8, because every continuation line starts with
    '#'. Delete one '#' and the value silently becomes the string
    "8 is sufficient ...", which reaches HyPhy as a thread count. This sweep
    covers every key rather than only the ones enumerated under `required`.
    """
    folded = [
        f"  {path} = {value!r}"
        for path, value in flatten(load_config(cfg)).items()
        if isinstance(value, str) and re.fullmatch(r"\s*-?\d+(\.\d+)?\s+\S.*", value, re.S)
    ]
    assert not folded, (
        f"{cfg.name}: a number has folded together with prose into one string:\n"
        + "\n".join(folded)
        + "\nPut the value on the same line as its key."
    )


def test_every_config_declares_a_distinct_build_name():
    """Two configs sharing a build_name would write into the same results/
    directory, and the second run would silently consume the first one's
    artefacts."""
    names = {}
    for cfg in CONFIG_PATHS:
        names.setdefault(load_config(cfg).get("build_name"), []).append(cfg.name)
    clashes = {name: files for name, files in names.items() if len(files) > 1}
    assert not clashes, f"configs share a build_name and so share results/: {clashes}"
