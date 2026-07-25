"""The ordering-sheet recommendation label must follow the selection verdict.

Bug class 5 in the changelog: "BEST RECOMMENDED (RANK 1)" was applied by list
position. A combination that sorted first but failed the selection gates was
therefore presented to an oligo house as the best assay to order. The label is
now derived by recommendation_category(); these tests pin the mapping so it
cannot regress to being position-driven.
"""
import pytest

from conftest import load_module

pytestmark = pytest.mark.unit

es = load_module("scripts/export_supplementary.py", "export_supplementary")


def test_rank_one_that_passed_is_the_best():
    assert es.recommendation_category(1, {"selection_passed": True}) == \
        "BEST RECOMMENDED (RANK 1)"


def test_lower_ranks_are_alternatives():
    assert es.recommendation_category(2, {"selection_passed": True}) == "ALTERNATIVE RANK 2"
    assert es.recommendation_category(3, {"selection_passed": True}) == "ALTERNATIVE RANK 3"


def test_a_failed_combo_is_never_recommended_even_at_rank_one():
    # The exact inversion the bug produced: rank 1, but it failed the gates.
    cat = es.recommendation_category(1, {"selection_passed": False})
    assert cat == "NOT RECOMMENDED — failed selection gates"
    assert "BEST" not in cat


def test_a_failed_combo_at_a_lower_rank_is_also_not_recommended():
    assert es.recommendation_category(4, {"selection_passed": False}) == \
        "NOT RECOMMENDED — failed selection gates"


def test_missing_verdict_defaults_to_passed():
    # Older combo JSON without the field is treated as passed, matching the
    # get(..., True) default the exporter has always used.
    assert es.recommendation_category(1, {}) == "BEST RECOMMENDED (RANK 1)"


def test_relaxed_floors_are_disclosed_on_the_best_row():
    cat = es.recommendation_category(
        1, {"selection_passed": True, "selection_relaxed": True, "selection_level": "L3"}
    )
    assert cat.startswith("BEST RECOMMENDED (RANK 1)")
    assert "relaxed to L3" in cat


def test_relaxation_is_not_disclosed_on_alternatives():
    # The relaxed suffix only applies to the rank-1 row.
    cat = es.recommendation_category(
        2, {"selection_passed": True, "selection_relaxed": True, "selection_level": "L3"}
    )
    assert cat == "ALTERNATIVE RANK 2"
