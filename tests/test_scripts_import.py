"""Import-smoke for pipeline scripts that have no other coverage.

Byte-compilation (CI's compileall step) catches syntax errors, and ruff F821
catches undefined names, but neither runs a script's module body -- a broken
top-level import or a statement that fails at load time slips through both and
only surfaces on the cluster. Importing each script here closes that gap cheaply.

Excluded, with reason:
  - scripts/prune_tree.py, scripts/m02_filter_wrapper.py: read the `snakemake`
    global that Snakemake injects into script:-directive scripts, so they cannot
    be imported standalone. ruff F821 (with builtins=["snakemake"]) and the
    pyflakes sweep in test_scripts_execute.py vet them instead.
  - scripts/m07_plot_entropy.py: sys.exit(1) at import if matplotlib is absent,
    and matplotlib is deliberately not a CI dependency.
"""
import importlib.util

import pytest

from conftest import REPO

pytestmark = pytest.mark.unit

# Scripts with no execution or import coverage elsewhere, that load cleanly with
# only the CI dependency set (pandas, Biopython, pyyaml).
IMPORTABLE = [
    "scripts/early_filter.py",
    "scripts/enrich_metadata.py",
    "scripts/sanitize_metadata.py",
    "scripts/extract_genes.py",
    "scripts/hyphy_subsample.py",
    "scripts/run_hyphy_gard_safe.py",
    "scripts/m08_payload_generator.py",
    "scripts/m10_compare_runs.py",
    "scripts/build_custom_blast_db.py",
    "scripts/tools/design_template.py",
    "scripts/tools/fetch_flu_references.py",
    "scripts/tools/fetch_h3n2_ref.py",
    "scripts/tools/fetch_universal_refs.py",
]


@pytest.mark.parametrize("rel_path", IMPORTABLE)
def test_script_imports_cleanly(rel_path):
    spec = importlib.util.spec_from_file_location(rel_path, REPO / rel_path)
    module = importlib.util.module_from_spec(spec)
    # A main() guarded by `if __name__ == "__main__"` does not run on import, so
    # this exercises the module body without triggering the CLI.
    spec.loader.exec_module(module)
