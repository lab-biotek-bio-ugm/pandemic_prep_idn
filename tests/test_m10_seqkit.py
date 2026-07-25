"""Tests for the seqkit amplicon wrapper in scripts/m10_insilico_pcr.py.

The wrapper is fail-closed on purpose: a seqkit crash or a missing binary must
raise, so it is never confused with a legitimate empty (no-amplicon) result --
an in-silico PCR that silently reports "no coverage" because the tool fell over
would reject good primers. These tests drive that wrapper against a fake seqkit
shim whose exit code and output are controlled, so both branches are exercised
without the real binary.
"""
import subprocess
import textwrap
from pathlib import Path

import pytest

from conftest import REPO, load_module
from fixtures import install_shim, path_with

pytestmark = pytest.mark.smoke

m10 = load_module("scripts/m10_insilico_pcr.py", "m10_insilico_pcr")

# seqkit is invoked as `seqkit amplicon -F <fwd> -R <rev> -m N -j T --bed <fasta>`
# and its --bed output is parsed as tab-separated columns, seq id in column 0 and
# amplicon sequence in column 6. The shim emits that shape, or nothing, or fails,
# according to $SEQKIT_MODE.
_SEQKIT_SHIM = textwrap.dedent(
    """
    import os, sys
    mode = os.environ.get("SEQKIT_MODE", "hit")
    if mode == "crash":
        sys.stderr.write("fake seqkit: boom\\n")
        sys.exit(1)
    if mode == "hit":
        # chrom start end name score strand amplicon
        sys.stdout.write("seqA\\t0\\t120\\tamp\\t0\\t+\\tACGTACGTACGTACGTACGT\\n")
        sys.stdout.write("seqB\\t0\\t118\\tamp\\t0\\t+\\tACGTACGTACGTACGTAC\\n")
    # mode == "empty": print nothing, exit 0
    sys.exit(0)
    """
)

FWD = "ACGTACGTACGT"   # non-degenerate, so the -F/-R branch is taken
REV = "TGCATGCATGCA"


def _with_seqkit(monkeypatch, tmp_path, mode):
    bin_dir = install_shim(tmp_path / "bin", "seqkit", _SEQKIT_SHIM)
    monkeypatch.setenv("PATH", path_with(bin_dir))
    monkeypatch.setenv("SEQKIT_MODE", mode)
    fasta = tmp_path / "seqs.fasta"
    fasta.write_text(">seqA\nACGT\n>seqB\nACGT\n", encoding="utf-8")
    return str(fasta)


def test_amplicons_are_parsed_from_bed_output(monkeypatch, tmp_path):
    fasta = _with_seqkit(monkeypatch, tmp_path, "hit")
    matches = m10.run_seqkit_amplicon(FWD, REV, max_mismatch=2, fasta_path=fasta)
    assert set(matches) == {"seqA", "seqB"}
    assert matches["seqA"] == "ACGTACGTACGTACGTACGT"


def test_no_amplicon_is_an_empty_dict_not_an_error(monkeypatch, tmp_path):
    fasta = _with_seqkit(monkeypatch, tmp_path, "empty")
    matches = m10.run_seqkit_amplicon(FWD, REV, max_mismatch=2, fasta_path=fasta)
    assert matches == {}


def test_a_seqkit_crash_raises_rather_than_reporting_no_coverage(monkeypatch, tmp_path):
    fasta = _with_seqkit(monkeypatch, tmp_path, "crash")
    # Fail-closed: the wrapper must raise, not return {} (which downstream would
    # read as "this primer amplifies nothing").
    with pytest.raises((m10.ToolError, subprocess.CalledProcessError, RuntimeError)):
        m10.run_seqkit_amplicon(FWD, REV, max_mismatch=2, fasta_path=fasta)


def test_a_missing_seqkit_binary_raises(monkeypatch, tmp_path):
    monkeypatch.setenv("PATH", "/nonexistent")
    fasta = tmp_path / "seqs.fasta"
    fasta.write_text(">seqA\nACGT\n", encoding="utf-8")
    with pytest.raises((m10.ToolError, FileNotFoundError, RuntimeError)):
        m10.run_seqkit_amplicon(FWD, REV, max_mismatch=2, fasta_path=str(fasta))
