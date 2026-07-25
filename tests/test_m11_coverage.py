"""Producer/consumer contract for the M10 -> M11 coverage hand-off.

M10 writes a per-candidate RAA detection coverage into a summary TSV, and the
multiplex and triplex checkers read it back to decide what to order. The two
sides agreeing on the column name is load-bearing: the changelog records that
this number "never reached the step that picks what to order". Both checkers
must read the same column, and a missing file must read as "no data" rather than
silently passing.
"""
import pytest

from conftest import load_module

pytestmark = pytest.mark.unit

mux = load_module("scripts/m11_multiplex_checker.py", "m11_multiplex_checker")
tri = load_module("scripts/m11_triplex_checker.py", "m11_triplex_checker")

CHECKERS = pytest.mark.parametrize(
    "mod", [mux, tri], ids=["multiplex", "triplex"]
)

# The columns M10 actually writes for a populated result (scripts/m10_insilico_pcr.py
# emits raa_detected_pct alongside pcr_/raa_ amplified and detected percentages).
_M10_HEADER = "rank\tregion\tgene\traa_amplified_pct\traa_detected_pct\n"


def _write_summary(path, rows):
    path.write_text(_M10_HEADER + "".join(rows), encoding="utf-8")
    return str(path)


@CHECKERS
def test_coverage_is_read_by_rank(mod, tmp_path):
    summary = _write_summary(tmp_path / "coverage_summary.tsv", [
        "1\tMP\tM1\t99.0\t95.5\n",
        "2\tMP\tM1\t88.0\t70.0\n",
    ])
    coverage = mod.load_coverage(summary)
    assert coverage == {1: 95.5, 2: 70.0}


@CHECKERS
def test_missing_file_reads_as_no_data(mod, tmp_path):
    assert mod.load_coverage(str(tmp_path / "does_not_exist.tsv")) == {}
    assert mod.load_coverage(None) == {}


@CHECKERS
def test_wrong_coverage_column_yields_no_data_not_a_crash(mod, tmp_path):
    # If M10's column were renamed, load_coverage must degrade to "no data"
    # rather than raise -- the gate then reports the absence instead of passing.
    bad = tmp_path / "coverage_summary.tsv"
    bad.write_text("rank\tdetected_pct\n1\t95.5\n", encoding="utf-8")
    assert mod.load_coverage(str(bad)) == {}


def test_attach_coverage_maps_onto_candidates():
    # Only the multiplex checker exposes attach_coverage; the triplex checker
    # threads coverage differently.
    cands = [{"rank": 1}, {"rank": 2}, {"rank": 3}]
    out = mux.attach_coverage(cands, {1: 95.5, 2: 70.0})
    assert [c["coverage"] for c in out] == [95.5, 70.0, None]


def test_both_checkers_read_the_same_column(tmp_path):
    # A behavioural cross-check: feed one file to both, and they must agree.
    # This is what stops one checker's column name from drifting from the other's.
    summary = _write_summary(tmp_path / "coverage_summary.tsv", ["1\tMP\tM1\t99.0\t42.0\n"])
    assert mux.load_coverage(summary) == tri.load_coverage(summary) == {1: 42.0}
