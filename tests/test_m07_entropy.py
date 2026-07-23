"""Regression tests for M07 Shannon-entropy scoring (scripts/m07_bridge_builder.py).

Guards the gap-handling fix: the documented alphabet is {A,T,G,C,-}, so deletions
must raise entropy (lower conservation). The pre-fix code excluded gaps and would
mark a deletion hotspot as perfectly conserved.
"""
import math

import pytest

from conftest import load_module  # noqa: E402

m07 = load_module("scripts/m07_bridge_builder.py", "m07_bridge_builder")


class _Rec:
    """Minimal stand-in for a Bio.SeqRecord (only .seq indexing is used)."""
    def __init__(self, seq):
        self.seq = seq


def _entropy(column_bases):
    prof, _, _ = m07.build_entropy_profile([_Rec(b) for b in column_bases])
    return prof[0]["entropy"]


def _conservation(column_bases):
    return 1.0 - _entropy(column_bases) / m07.MAX_ENTROPY_BITS


def test_max_entropy_constant_is_log2_5():
    assert m07.MAX_ENTROPY_BITS == pytest.approx(math.log2(5))


def test_fully_conserved_column_has_zero_entropy():
    assert _entropy(["A"] * 10) == pytest.approx(0.0)
    assert _conservation(["A"] * 10) == pytest.approx(1.0)


def test_doc_worked_example_8A_1G_1gap():
    # docs/concepts/shannon-entropy.md worked example: H=0.9219, C=0.603
    col = ["A"] * 8 + ["G"] + ["-"]
    assert _entropy(col) == pytest.approx(0.9219, abs=1e-3)
    assert _conservation(col) == pytest.approx(0.603, abs=1e-3)


def test_deletion_hotspot_is_not_perfectly_conserved():
    # 20% deletion, 80% conserved-A. Pre-fix (gap-excluded) gave C=1.0 (a trap).
    col = ["A"] * 8 + ["-"] * 2
    assert _conservation(col) < 0.8


def test_all_gap_column_is_maximally_penalised():
    # An undesignable all-gap column must score C=0, never "perfectly conserved".
    assert _entropy(["-"] * 5) == pytest.approx(m07.MAX_ENTROPY_BITS)
    assert _conservation(["-"] * 5) == pytest.approx(0.0)


def test_more_deletion_lowers_conservation_monotonically():
    c10 = _conservation(["A"] * 9 + ["-"])
    c30 = _conservation(["A"] * 7 + ["-"] * 3)
    c50 = _conservation(["A"] * 5 + ["-"] * 5)
    assert c10 > c30 > c50


def test_ambiguous_bases_folded_into_variability():
    # N is not a reliable base; it should reduce conservation like a gap.
    assert _conservation(["A"] * 8 + ["N"] * 2) < 1.0
