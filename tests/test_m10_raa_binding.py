"""Regression tests for M10 RAA binding assessment (scripts/m10_insilico_pcr.py).

Guards the reverse-primer 3'-orientation fix. Reverse primers are supplied in
GENOMIC orientation (RC of the 5'->3' oligo), so the oligo's critical 3' terminus
sits at index 0. The pre-fix code always treated the array END as 3', mislabelling
a real 3'-terminal mismatch as a tolerant "internal" mismatch and passing assays
that would fail at the bench.
"""
import pytest

from conftest import load_module  # noqa: E402

m10 = load_module("scripts/m10_insilico_pcr.py", "m10_insilico_pcr")

rc = m10.reverse_complement


def _assessors():
    funcs = [m10.assess_raa_binding]
    if getattr(m10, "_HAS_NUMPY", False):
        funcs.append(m10._assess_raa_numpy)
    return funcs


def test_perfect_forward_primer_binds():
    for f in _assessors():
        ok, t3, ti, tot = f("ACGTACGTACGTACGT", "ACGTACGTACGTACGT" + "TTTT",
                            critical_3prime_bp=8, max_3prime_mm=1, max_internal_mm=4)
        assert ok and tot == 0


def test_forward_3prime_mismatch_counted_at_right_end():
    # Forward primer given 5'->3'; a mismatch at its last base is 3'-critical.
    primer = "AAAAAAAAAAAAAAAG"          # 3' terminus = final G
    amp = "AAAAAAAAAAAAAAAA" + "TTTTTTTT"  # last base A mismatches the primer's G
    for f in _assessors():
        ok, t3, ti, tot = f(primer, amp, critical_3prime_bp=8,
                            max_3prime_mm=0, max_internal_mm=4, three_prime_at="right")
        assert t3 == 1 and ti == 0
        assert not ok  # max_3prime_mm=0 -> a 3' mismatch must reject


def test_reverse_primer_3prime_terminus_at_left_index0():
    # Reverse oligo (5'->3') with 3' end ...C G; supplied to M10 as its genomic RC.
    rev_oligo = "AAAAAAAAAAAAAAAACG"
    genomic = rc(rev_oligo)               # what M10 actually passes
    # Introduce ONE mismatch at the oligo's true 3' terminus = index 0 of genomic.
    amp = "A" + genomic[1:] + "TTTTTTTTTTTTTTTTTTTT"
    for f in _assessors():
        ok, t3, ti, tot = f(genomic, amp, critical_3prime_bp=8,
                            max_3prime_mm=1, max_internal_mm=4, three_prime_at="left")
        assert t3 == 1 and ti == 0  # correctly attributed to the 3' region


def test_reverse_primer_two_3prime_mismatches_reject():
    # The decisive case: 2 mismatches at the true 3' terminus must REJECT.
    rev_oligo = "AAAAAAAAAAAAAAAACG"
    genomic = rc(rev_oligo)
    amp = "AA" + genomic[2:] + "TTTTTTTTTTTTTTTTTTTT"
    for f in _assessors():
        ok, t3, ti, tot = f(genomic, amp, critical_3prime_bp=8,
                            max_3prime_mm=1, max_internal_mm=4, three_prime_at="left")
        assert t3 == 2 and not ok


def test_old_orientation_would_have_masked_the_defect():
    # Documents the bug: with the WRONG side ('right'), the same 3' mismatches are
    # miscounted as internal and the bad assay is (incorrectly) accepted.
    rev_oligo = "AAAAAAAAAAAAAAAACG"
    genomic = rc(rev_oligo)
    amp = "AA" + genomic[2:] + "TTTTTTTTTTTTTTTTTTTT"
    ok_wrong, t3_wrong, ti_wrong, _ = m10.assess_raa_binding(
        genomic, amp, critical_3prime_bp=8, max_3prime_mm=1, max_internal_mm=4,
        three_prime_at="right")
    assert ti_wrong == 2 and t3_wrong == 0 and ok_wrong  # the masked (buggy) behaviour
