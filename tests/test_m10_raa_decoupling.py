"""Guard the M10 RAA-vs-PCR decoupling.

M10 used to gate RAA amplification by the PCR seqkit `-m max_mismatch` result, so
`raa_detected_pct` could never exceed PCR coverage — backwards for RAA, which
tolerates more mismatches (Higgins 2022; TwistAmp 2.7). The fix runs seqkit at the
RAA budget and derives PCR (post-filter to max_mismatch) and RAA (3'/internal gate)
from the same amplicons. These tests pin the decoupling at the gate level: a
sequence RAA would amplify but PCR would not must now be RAA-amplified, and a
PCR-amplifiable sequence must stay PCR-amplified.

The seqkit call itself is not exercised (it needs the binary + a long run); this
validates the per-sequence gate logic that seqkit feeds, which is where the bug was.
"""
import pytest

from conftest import load_module

pytestmark = pytest.mark.unit

m10 = load_module("scripts/m10_insilico_pcr.py", "m10_insilico_pcr")

PRIMER = "ACGTACGTACGTACGTACGT"   # 20 bp; 3' critical window = last 4 (positions 16-19)
MAX_MISMATCH = 3                   # PCR budget
CRIT, MAX_3P, MAX_INT = 4, 1, 4    # RAA budget


def _amplicon(mismatch_positions):
    """PRIMER with the given positions substituted to a non-matching base, padded."""
    swap = {"A": "C", "C": "A", "G": "T", "T": "G"}
    bases = list(PRIMER)
    for p in mismatch_positions:
        bases[p] = swap[bases[p]]
    return "GGGGG" + "".join(bases) + "GGGGG"


def _gates(amp):
    """Reproduce M10's per-sequence PCR and RAA gate decisions for one primer."""
    pcr_mm = m10.min_mismatches_sliding_window(PRIMER, amp)
    pcr_amplified = pcr_mm <= MAX_MISMATCH
    raa_ok, _t3, _ti, _tot = m10.assess_raa_binding(
        PRIMER, amp, critical_3prime_bp=CRIT, max_3prime_mm=MAX_3P, max_internal_mm=MAX_INT)
    return pcr_mm, pcr_amplified, raa_ok


def test_raa_budget_is_wider_than_pcr():
    # raa_budget = max(max_mismatch, raa_max_3prime_mm + raa_max_internal_mm)
    assert max(MAX_MISMATCH, MAX_3P + MAX_INT) == 5 > MAX_MISMATCH


def test_four_internal_mismatches_raa_amplifies_pcr_does_not():
    # 4 internal mismatches, none in the 3' anchor: PCR (>3) drops it, RAA keeps it.
    pcr_mm, pcr_amp, raa_ok = _gates(_amplicon([2, 5, 8, 11]))
    assert pcr_mm == 4
    assert pcr_amp is False       # PCR would NOT amplify (was the old cap)
    assert raa_ok is True         # RAA DOES — decoupled from the PCR gate


def test_low_mismatch_stays_pcr_amplified():
    # 2 internal mismatches: within both budgets. PCR must remain amplified
    # (the fix must not lower PCR coverage).
    pcr_mm, pcr_amp, raa_ok = _gates(_amplicon([5, 11]))
    assert pcr_mm == 2
    assert pcr_amp is True
    assert raa_ok is True


def test_two_3prime_mismatches_raa_rejects_pcr_accepts():
    # 2 mismatches inside the 3' anchor: RAA (max_3prime_mm=1) vetoes, PCR (<=3) accepts.
    # RAA stays STRICTER than PCR where it should — decoupling is not blanket leniency.
    pcr_mm, pcr_amp, raa_ok = _gates(_amplicon([18, 19]))
    assert pcr_mm == 2
    assert pcr_amp is True
    assert raa_ok is False


def test_perfect_match_amplifies_both():
    pcr_mm, pcr_amp, raa_ok = _gates(_amplicon([]))
    assert pcr_mm == 0 and pcr_amp is True and raa_ok is True
