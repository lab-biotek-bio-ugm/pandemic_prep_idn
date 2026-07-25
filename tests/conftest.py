"""Shared pytest fixtures / path setup for the pipeline unit tests.

Makes `scripts/` and `scripts/tools/` importable so the core scientific
functions can be tested in isolation from Snakemake.
"""
import importlib.util
import os
import sys
from pathlib import Path

import yaml

_REPO = Path(__file__).resolve().parent.parent
_SCRIPTS = _REPO / "scripts"

# Public alias. `_REPO` is kept because existing tests import it.
REPO = _REPO

for p in (_SCRIPTS, _SCRIPTS / "tools"):
    sp = str(p)
    if sp not in sys.path:
        sys.path.insert(0, sp)

# Pin the working directory at conftest import time, i.e. before collection.
# Dozens of tests use bare relative paths -- open("config/config.yaml"),
# glob.glob("scripts/*.py"), subprocess.run([..., "scripts/m07_bridge_builder.py"])
# -- and seven test modules call load_module() at module scope, which runs during
# collection. A session-scoped autouse fixture would run too late for those, so
# this has to stay a module-level statement.
os.chdir(_REPO)

# The pathogen configs. glob("config*.yaml") excludes config/primerlab_template*.yaml,
# which are PrimerLab engine templates rather than pipeline configs. The project
# currently targets influenza (Orthomyxoviridae); a new pathogen config dropped
# in here is picked up automatically and must satisfy the contract.
CONFIG_PATHS = sorted((_REPO / "config").glob("config*.yaml"))


def load_module(rel_path: str, name: str):
    """Load a script by file path (handles scripts that aren't importable by name)."""
    spec = importlib.util.spec_from_file_location(name, _REPO / rel_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class _NoDuplicateLoader(yaml.SafeLoader):
    """SafeLoader that rejects a repeated mapping key instead of keeping the last.

    PyYAML silently keeps the final occurrence. In a 700-line config maintained
    by search-and-replace across seven files, a duplicated key is a threshold
    override with no diff signal at all.
    """


def _no_duplicate_mapping(loader, node, deep=False):
    seen, mapping = set(), {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in seen:
            raise yaml.constructor.ConstructorError(
                None, None, f"duplicate key {key!r}", key_node.start_mark
            )
        seen.add(key)
        mapping[key] = loader.construct_object(value_node, deep=deep)
    return mapping


_NoDuplicateLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _no_duplicate_mapping
)


def load_config(path):
    """Load one config file STANDALONE, with duplicate keys rejected.

    Loading each file alone is the entire point of the contract suite. Snakefile
    line 6 declares `configfile: "config/config.yaml"`, and Snakemake deep-merges
    any `--configfile` on top of that rather than replacing it. A per-pathogen
    config that omits a key therefore does not fail; it silently inherits the
    influenza-A-MP value. No dry run can see that hole, because after the merge
    there is no hole. This can.
    """
    return yaml.load(Path(path).read_text(encoding="utf-8"), _NoDuplicateLoader)


def merged_config(path):
    """Reproduce the config Snakemake actually hands the rules for `path`.

    Mirrors snakemake.utils.update_config: a recursive dict merge in which the
    chosen config wins over config/config.yaml, leaf by leaf. Verified to agree
    with a real `snakemake --configfile <path>` run for every shipped config.

    Use this when asserting on what the pipeline WILL do; use load_config when
    asserting on what a file DECLARES.
    """
    base = load_config(_REPO / "config" / "config.yaml")
    if Path(path).resolve() == (_REPO / "config" / "config.yaml").resolve():
        return base
    return _deep_merge(base, load_config(path))


def _deep_merge(base: dict, override: dict) -> dict:
    out = dict(base)
    for key, value in override.items():
        if isinstance(out.get(key), dict) and isinstance(value, dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = value
    return out


def flatten(mapping, prefix=""):
    """Dotted leaf path -> value. Lists are leaves.

    An empty dict is also treated as a leaf. Recursing into it without emitting
    anything would make `multiplex.dyes: {}` invisible to key-parity checks --
    which is exactly the shape a block takes when it is stubbed out in one
    config and populated in another.
    """
    out = {}
    for key, value in (mapping or {}).items():
        path = f"{prefix}.{key}" if prefix else str(key)
        if isinstance(value, dict) and value:
            out.update(flatten(value, path))
        else:
            out[path] = value
    return out
