"""Tests for scripts/m12_coverage_attribution.py.

Pins the statistics (G-test, BH-FDR) against hand values, the driver detection,
and the end-to-end two-lens attribution: a subtype with a dominant driver is
'mechanistic' and names the position; a subtype absent from the design set
short-circuits as 'no_focal_data'.
"""
import json
import math

import pytest

from conftest import load_module

pytestmark = pytest.mark.unit

m12 = load_module("scripts/m12_coverage_attribution.py", "m12_coverage_attribution")


# --- statistics ---------------------------------------------------------------

def test_g_test_perfect_separation_is_significant():
    g, p = m12.g_test_2x2(40, 0, 0, 40)   # group1 all-match, group2 all-mismatch
    assert g > 50 and p < 1e-10


def test_g_test_no_association_is_not_significant():
    g, p = m12.g_test_2x2(20, 20, 20, 20)
    assert g == pytest.approx(0.0, abs=1e-9)
    assert p == pytest.approx(1.0, abs=1e-9)


def test_g_test_matches_hand_computed():
    # 2x2 [[8,2],[2,8]]; G = 2*sum o*ln(o/e), e=5 each -> 2*4*(8ln1.6+2ln0.4)
    g, _p = m12.g_test_2x2(8, 2, 2, 8)
    expected = 2 * (8 * math.log(8 / 5) + 2 * math.log(2 / 5)) * 2
    assert g == pytest.approx(expected, rel=1e-9)


def test_bh_fdr_is_monotone_and_bounded():
    adj = m12.bh_fdr([0.01, 0.02, 0.5, 0.9])
    assert all(0 <= a <= 1 for a in adj)
    assert adj[0] <= adj[1]  # smaller raw p -> smaller/equal adjusted


# --- driver detection ---------------------------------------------------------

def test_attribute_oligo_flags_the_driver_position():
    oligo = "ACGTACGTAC"
    # detected group: all match the oligo; failing group: all differ at position 9.
    detected = [list(oligo) for _ in range(30)]
    failing = [list(oligo[:9] + "A") for _ in range(30)]  # pos 9 C->A
    res = m12.attribute_oligo("reverse", oligo, detected, failing, crit_3prime=4)
    drivers = {d["pos"] for d in res["driver_positions"]}
    assert drivers == {9}
    assert res["has_3prime_driver"] is True   # pos 9 is within the last 4 bases


# --- end-to-end ---------------------------------------------------------------

FWD, REV, PROBE = "ACGTACGTAC", "TTTTGGGGCC", "GATTACAGAT"


def _seq(rev_region):
    return "GGGGG" + FWD + "AAAA" + rev_region + "AAAA" + PROBE + "GGGGG"


def _write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _make_run(root):
    build = root / "influenza-x-mp"
    # coverage summary: HX (has design data) and HZ (absent from design) both fail.
    cov = [{
        "rank": 1, "fwd_seq": FWD, "rev_seq": REV, "probe_seq": PROBE,
        "fwd_degenerate_seq": FWD, "rev_degenerate_seq": REV,
        "raa_params": {"critical_3prime_bp": 4, "max_3prime_mm": 1, "max_internal_mm": 4},
        "breakdowns": {"subtype": [
            {"subtype": "A / HX", "raa_detected_pct": 45.0, "total": 80},
            {"subtype": "A / HZ", "raa_detected_pct": 10.0, "total": 10},
        ]},
    }]
    _write(build / "m10" / "degenerate" / "coverage_summary.json", json.dumps(cov))

    # metadata + FASTA. HX detected (clade SEEN, Indonesia) vs failing (clade FAIL,
    # Japan; two 3'-terminal mismatches in the reverse site). HZ: clade not in design.
    meta = ["strain\tsubtype\tclade\tcountry"]
    raw = []
    for i in range(40):
        s = f"HXdet{i}"
        meta.append(f"{s}\tA / HX\tSEEN\tIndonesia")
        raw.append(f">{s}\n{_seq(REV)}")
    for i in range(40):
        s = f"HXfail{i}"
        meta.append(f"{s}\tA / HX\tFAIL\tJapan")
        raw.append(f">{s}\n{_seq('TTTTGGGGAA')}")   # positions 8,9 C->A (3')
    for i in range(10):
        s = f"HZ{i}"
        meta.append(f"{s}\tA / HZ\tOTHER\tChina")
        raw.append(f">{s}\n{_seq(REV)}")
    _write(build / "m01" / "metadata_enriched.tsv", "\n".join(meta) + "\n")
    _write(build / "m01" / "raw_sequences.fasta", "\n".join(raw) + "\n")
    # design set: only HX clade SEEN -> HX has focal data, HZ does not.
    _write(build / "m02" / "filtered_sequences.fasta",
           "\n".join(f">HXdet{i}\n{_seq(REV)}" for i in range(10)) + "\n")
    return "influenza-x-mp"


def test_end_to_end_two_lens_and_driver(tmp_path):
    build = _make_run(tmp_path)
    result = m12.attribute_build(str(tmp_path), build, rank=1, max_mismatch=3)
    by_subtype = {r["subtype"]: r for r in result["subtypes"]}

    hx = by_subtype["A / HX"]
    assert hx["failure_mode"] == "mechanistic"
    # RAA-viable (detected) well above the reported/strict-PCR number.
    assert hx["raa_viable_pct"] == pytest.approx(50.0, abs=5.0)  # ~40/80 detected
    assert hx["failing_n"] >= 30
    assert ("Japan", 40) in hx["failing_top_countries"]
    drivers = [d for o in hx["oligos"] for d in o["driver_positions"]]
    assert any(d["three_prime"] for d in drivers)     # the 3'-terminal C->A driver
    assert "reverse" in hx["explanation"]

    hz = by_subtype["A / HZ"]
    assert hz["failure_mode"] == "no_focal_data"

    assert "max_mismatch=3" in result["finding_note"]


def test_finding_note_flags_the_pcr_cap(tmp_path):
    build = _make_run(tmp_path)
    result = m12.attribute_build(str(tmp_path), build, rank=1)
    assert "seqkit" in result["finding_note"] and "RAA" in result["finding_note"]
