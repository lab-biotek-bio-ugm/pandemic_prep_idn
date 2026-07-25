"""End-to-end tests for scripts/m09_blast_verification.py.

The script is 75 KB and had no execution coverage at all. Its defining bug was
class 4 in the changelog: `primerlab blast` exits 1 to mean "this oligo set is
NOT specific" -- a result, not an error -- and the old code reported those 28
routine verdicts as "[FAIL] BLAST verification failed", making them
indistinguishable from a genuine failure (missing binary, unreadable database).

These tests drive the real script as a subprocess against a fake `primerlab`
shim whose exit code and written artefact are controlled independently, so the
separation the bug was about is asserted directly. No blastn, no network.
"""
import json
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from conftest import REPO
from fixtures import install_shim, path_with

pytestmark = pytest.mark.smoke

M09 = REPO / "scripts" / "m09_blast_verification.py"

# A primerlab shim. It parses --output from argv, optionally writes a
# blast_result.json verdict there, prints a line, then exits with $SHIM_EXIT.
# Behaviour is chosen through the environment so one shim body serves every case.
_SHIM = textwrap.dedent(
    """
    import json, os, sys
    from pathlib import Path

    argv = sys.argv[1:]
    out = None
    for i, a in enumerate(argv):
        if a == "--output":
            out = argv[i + 1]

    if os.environ.get("SHIM_ARTEFACT") == "1" and out is not None:
        Path(out).mkdir(parents=True, exist_ok=True)
        verdict = {
            "is_specific": os.environ["SHIM_SPECIFIC"] == "1",
            "combined_score": float(os.environ.get("SHIM_SCORE", "0")),
            "hits": [],
        }
        (Path(out) / "blast_result.json").write_text(json.dumps(verdict))

    print("fake primerlab blast: done")
    sys.exit(int(os.environ.get("SHIM_EXIT", "0")))
    """
)


def _write_alternatives(path, count=2):
    """Minimal alternatives.json: a candidate needs a forward and a reverse
    primer to be processed rather than skipped."""
    candidates = [
        {
            "final_rank": i + 1,
            "primers": {
                "forward": {"sequence": "ACGTACGTACGTACGTACGT"},
                "reverse": {"sequence": "TGCATGCATGCATGCATGCA"},
                "probe": {"sequence": "GATTACAGATTACAGATTAC"},
            },
        }
        for i in range(count)
    ]
    Path(path).write_text(json.dumps(candidates), encoding="utf-8")


def _run_m09(tmp_path, *, specific, exit_code, artefact, candidates=2):
    alternatives = tmp_path / "alternatives.json"
    _write_alternatives(alternatives, candidates)
    outdir = tmp_path / "m09"
    bin_dir = install_shim(tmp_path / "bin", "primerlab", _SHIM)

    env = {
        "PATH": path_with(bin_dir),
        # Keep the benchmark instrumentation from writing results/<build>/.
        "PANDEMIC_BENCHMARK_OFF": "1",
        "SHIM_EXIT": str(exit_code),
        "SHIM_SPECIFIC": "1" if specific else "0",
        "SHIM_ARTEFACT": "1" if artefact else "0",
        "SHIM_SCORE": "0.99" if specific else "0.10",
    }
    proc = subprocess.run(
        [sys.executable, str(M09),
         "--alternatives", str(alternatives),
         "--database", "custom_local_db",
         "--mode", "local",
         "--outdir", str(outdir),
         "--num-candidates", str(candidates)],
        capture_output=True, text=True, env=env, cwd=str(REPO),
    )
    return proc, outdir


def test_specific_verdict_is_a_success(tmp_path):
    proc, outdir = _run_m09(tmp_path, specific=True, exit_code=0, artefact=True)
    assert proc.returncode == 0, proc.stderr
    assert "[SUCCESS]" in proc.stdout
    assert "[FAIL]" not in proc.stdout
    assert (outdir / "verification_complete.txt").exists()


def test_nonspecific_verdict_is_a_result_not_a_failure(tmp_path):
    # The heart of bug class 4: primerlab exits 1, but it wrote a complete
    # verdict saying "not specific". That is a finished verification, not an
    # error, so the run must succeed and must not report [FAIL].
    proc, outdir = _run_m09(tmp_path, specific=False, exit_code=1, artefact=True)
    assert proc.returncode == 0, proc.stderr
    assert "[NOT SPECIFIC]" in proc.stdout
    assert "[FAIL]" not in proc.stdout
    assert (outdir / "verification_complete.txt").exists()


def test_missing_artefact_is_a_genuine_failure(tmp_path):
    # primerlab exited non-zero AND left no readable artefact: specificity is
    # unknown. With every candidate failing and none skipped, the script must
    # abort non-zero so downstream modules never treat unverified primers as
    # specific.
    proc, _ = _run_m09(tmp_path, specific=False, exit_code=1, artefact=False)
    assert proc.returncode == 1
    assert "[FAIL]" in proc.stdout
    assert "FATAL" in proc.stdout


def test_missing_binary_is_a_fatal_failure(tmp_path):
    # No shim on PATH at all -> the subprocess call raises -> genuine failure.
    alternatives = tmp_path / "alternatives.json"
    _write_alternatives(alternatives, 1)
    proc = subprocess.run(
        [sys.executable, str(M09),
         "--alternatives", str(alternatives),
         "--database", "custom_local_db",
         "--mode", "local",
         "--outdir", str(tmp_path / "m09"),
         "--num-candidates", "1"],
        capture_output=True, text=True,
        env={"PATH": "/nonexistent", "PANDEMIC_BENCHMARK_OFF": "1"},
        cwd=str(REPO),
    )
    assert proc.returncode == 1
    assert "FATAL" in proc.stdout


def test_empty_alternatives_completes_without_running_blast(tmp_path):
    alternatives = tmp_path / "alternatives.json"
    alternatives.write_text("[]", encoding="utf-8")
    outdir = tmp_path / "m09"
    proc = subprocess.run(
        [sys.executable, str(M09),
         "--alternatives", str(alternatives),
         "--database", "custom_local_db",
         "--mode", "local",
         "--outdir", str(outdir)],
        capture_output=True, text=True,
        env={"PATH": "/nonexistent", "PANDEMIC_BENCHMARK_OFF": "1"},
        cwd=str(REPO),
    )
    assert proc.returncode == 0, proc.stderr
    assert (outdir / "verification_complete.txt").exists()
