"""Guard the M10 sliding-window refactor.

best_binding() was introduced to expose WHERE mismatches fall (for coverage-
failure attribution) while min_mismatches_sliding_window() became a thin wrapper.
The refactor must not change the mismatch COUNT — coverage numbers depend on it,
and re-running M10 to check takes an hour. This test pins the count against an
independent reference copy of the original algorithm over many random inputs, and
checks that the exposed positions/aligned bases are correct.
"""
import random

import pytest

from conftest import load_module

pytestmark = pytest.mark.unit

m10 = load_module("scripts/m10_insilico_pcr.py", "m10_insilico_pcr")


def _reference_min_mismatches(probe_seq, amplicon_seq):
    """Verbatim copy of the ORIGINAL min_mismatches_sliding_window algorithm."""
    probe_len = len(probe_seq)
    amp_len = len(amplicon_seq)
    if amp_len < probe_len:
        return probe_len
    if probe_seq in amplicon_seq or m10.reverse_complement(probe_seq) in amplicon_seq:
        return 0
    best = probe_len
    for p_seq in [probe_seq, m10.reverse_complement(probe_seq)]:
        for i in range(amp_len - probe_len + 1):
            window = amplicon_seq[i:i + probe_len]
            mm = 0
            for j in range(probe_len):
                if not m10.match_base(p_seq[j], window[j]):
                    mm += 1
                    if mm >= best:
                        break
            best = min(best, mm)
            if best == 0:
                return 0
    return best


def _rand_seq(rng, n, alphabet="ACGT"):
    return "".join(rng.choice(alphabet) for _ in range(n))


def test_count_matches_the_original_over_random_inputs():
    rng = random.Random(1234)
    for _ in range(400):
        oligo = _rand_seq(rng, rng.randint(4, 22))
        target = _rand_seq(rng, rng.randint(0, 60))
        assert m10.min_mismatches_sliding_window(oligo, target) == \
            _reference_min_mismatches(oligo, target)


def test_count_matches_with_iupac_degenerate_oligos():
    rng = random.Random(99)
    for _ in range(200):
        oligo = _rand_seq(rng, rng.randint(6, 18), alphabet="ACGTRYSWKMN")
        target = _rand_seq(rng, rng.randint(6, 40))
        assert m10.min_mismatches_sliding_window(oligo, target) == \
            _reference_min_mismatches(oligo, target)


def test_exact_match_reports_zero_and_offset():
    target = "GGGG" + "ACGTACGT" + "TTTT"
    count, offset, strand, positions, aligned = m10.best_binding("ACGTACGT", target)
    assert count == 0 and positions == []
    assert strand == "+" and offset == 4
    assert "".join(aligned) == "ACGTACGT"


def test_positions_are_in_oligo_frame_and_correct():
    # oligo perfectly matches target except one base; the reported position and
    # aligned base must point at that exact mismatch, in the oligo's 5'->3' frame.
    oligo = "ACGTACGT"
    target = "TT" + "ACGTAAGT" + "CC"   # position 5 (0-based) differs: C -> A
    count, offset, strand, positions, aligned = m10.best_binding(oligo, target)
    assert count == 1
    assert positions == [5]
    assert aligned[5] == "A"
    assert offset == 2 and strand == "+"


def test_reverse_strand_match():
    oligo = "AAACCCGTAG"  # not a reverse-complement palindrome
    assert m10.reverse_complement(oligo) != oligo
    target = m10.reverse_complement("GG" + oligo + "TT")
    count, offset, strand, positions, aligned = m10.best_binding(oligo, target)
    assert count == 0
    assert strand == "-"
    assert "".join(aligned) == oligo


def test_too_short_target_returns_len():
    assert m10.min_mismatches_sliding_window("ACGTACGT", "ACG") == 8
    count, offset, *_ = m10.best_binding("ACGTACGT", "ACG")
    assert count == 8 and offset is None
