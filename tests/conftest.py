"""Shared pytest fixtures / path setup for the pipeline unit tests.

Makes `scripts/` and `scripts/tools/` importable so the core scientific
functions can be tested in isolation from Snakemake.
"""
import importlib.util
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parent.parent
_SCRIPTS = _REPO / "scripts"

for p in (_SCRIPTS, _SCRIPTS / "tools"):
    sp = str(p)
    if sp not in sys.path:
        sys.path.insert(0, sp)


def load_module(rel_path: str, name: str):
    """Load a script by file path (handles scripts that aren't importable by name)."""
    spec = importlib.util.spec_from_file_location(name, _REPO / rel_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod
