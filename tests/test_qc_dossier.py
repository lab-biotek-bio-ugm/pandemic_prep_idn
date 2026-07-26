"""Tests for scripts/qc_dossier.py.

The dossier's core job is to classify each human off-target hit exactly as the
selection gate does (candidate_selection.FATAL_HIT) and to spell out why a
tolerated hit is safe. These tests pin that classification and drive the whole
generator end-to-end over a synthetic run, so a change that silently mislabels a
hit -- or drops the coverage-hole surfacing -- fails here.
"""
import json

import pytest

import qc_dossier

pytestmark = pytest.mark.unit


# --- per-hit classification (mirrors candidate_selection.FATAL_HIT) -----------

def test_a_fully_anchored_near_perfect_hit_is_fatal():
    hit = {"three_prime_overhang": 0, "mismatches": 1, "coverage_percent": 95.0,
           "is_significant": True, "evalue": 1e-9}
    verdict, _ = qc_dossier.classify_hit(hit)
    assert verdict == "FATAL"


def test_a_hit_with_a_3prime_overhang_is_tolerated():
    hit = {"three_prime_overhang": 5, "mismatches": 1, "coverage_percent": 95.0,
           "is_significant": False, "evalue": 1.2}
    verdict, reason = qc_dossier.classify_hit(hit)
    assert verdict == "tolerated"
    assert "not anchored" in reason


def test_a_low_coverage_hit_is_tolerated():
    hit = {"three_prime_overhang": 0, "mismatches": 1, "coverage_percent": 68.8,
           "is_significant": False, "evalue": 4.6}
    verdict, reason = qc_dossier.classify_hit(hit)
    assert verdict == "tolerated"
    assert "of the primer binds" in reason


def test_collect_human_hits_spans_all_three_oligos():
    blast = {
        "forward_result": {"offtargets": [{"sequence_id": "f", "three_prime_overhang": 5}]},
        "reverse_result": {"offtargets": [{"sequence_id": "r", "three_prime_overhang": 7}]},
        "probe_result": {"offtargets": []},
    }
    hits = qc_dossier.collect_human_hits(blast)
    assert {h["role"] for h in hits} == {"forward", "reverse"}
    assert qc_dossier.human_worst_hit(hits)["sequence_id"] == "f"  # smallest overhang


# --- end-to-end over a synthetic run ------------------------------------------

def _blast(offtargets, products=0, grade="A", score=92.6):
    return {
        "combined_score": score, "is_specific": True, "grade": grade,
        "potential_products": products,
        "forward_result": {"primer_seq": "ACGT", "specificity_score": 100.0,
                           "offtarget_count": 0, "offtargets": []},
        "reverse_result": {"primer_seq": "TGCA", "specificity_score": 66.0,
                           "offtarget_count": len(offtargets), "offtargets": offtargets},
        "probe_result": {"primer_seq": "GATTACA", "specificity_score": 100.0,
                         "offtarget_count": 0, "offtargets": []},
    }


def _coverage_record(rank, subtypes):
    return [{
        "rank": rank, "total_seqs": 50000, "raa_detected_pct": 91.17,
        "design_db_total": 700, "design_db_raa_det_pct": 95.13,
        "breakdowns": {"subtype": subtypes},
        "design_db_breakdowns": {"subtype": []},
    }]


def _write(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj), encoding="utf-8")


def _synthetic_run(root):
    tolerated = {"sequence_id": "ref|NC_060938.1|",
                 "sequence_title": "Homo sapiens chromosome 14, T2T-CHM13",
                 "identity": 95.5, "mismatches": 1, "coverage_percent": 68.8,
                 "three_prime_involved": True, "three_prime_overhang": 5,
                 "is_significant": False, "evalue": 1.2, "risk_level": "low"}
    for build, rank in (("influenza-a-mp-indo-only", 1), ("influenza-b-mp-indo-only", 2)):
        stem = root / build / "m09" / "degenerate" / "segment-MP_region_1" / f"candidate_{rank}"
        _write(stem / "db_human" / "blast_result.json", _blast([tolerated], products=0))
        _write(stem / "db_virus" / "blast_result.json", _blast([], products=0))
        subs = [{"subtype": "A / H3N2", "total": 19000, "raa_detected_pct": 99.0},
                {"subtype": "A / H5N5", "total": 35, "raa_detected_pct": 0.0}]  # a hole
        _write(root / build / "m10" / "degenerate" / "coverage_summary.json",
               _coverage_record(rank, subs))
    combo = {
        "pair_id": "FluA_Rank1_x_FluB_Rank2", "rank_a": 1, "rank_b": 2,
        "score": 64.5, "status_raa_39c": "PASS", "dye_ok": True,
        "min_dg_raa_39c": -7.19, "delta_tm": 2.23, "max_3prime_anchor_bp": 5,
        "selection_passed": True, "selection_level": "L1",
        "selection_level_label": "coverage vs best achievable",
        "selection_relaxed": True, "selection_rank": 1,
        "oligos_synthesis": {"flu_a": {}, "flu_b": {}},
        "cross_dimers_39c": {},
    }
    mdir = root / "multiplex-influenza-a-mp-indo-only_x_influenza-b-mp-indo-only"
    _write(mdir / "multiplex_candidates.json", [combo])


def test_end_to_end_dossier(tmp_path):
    root = tmp_path / "results"
    _synthetic_run(root)
    out = tmp_path / "qc"
    n = qc_dossier.build_dossiers(str(root), str(out), top_n=5)
    assert n == 1

    html = (out / "FluA_Rank1_x_FluB_Rank2.html").read_text(encoding="utf-8")
    md = (out / "FluA_Rank1_x_FluB_Rank2.md").read_text(encoding="utf-8")

    # The tolerated human hit is spelled out with its accession and a reason.
    assert "NC_060938.1" in html
    assert "tolerated" in md.lower()
    # No amplicon -> the conclusion states there is no false-positive mechanism.
    assert "No false-positive mechanism" in html
    # Both coverage figures are surfaced, and the gate note names the global one.
    assert "95.13" in html and "91.17" in html
    assert "gate judges" in html
    # The subtype hole is flagged, not hidden in the average.
    assert "H5N5" in html


def test_missing_build_does_not_crash(tmp_path):
    # A multiplex file whose member builds are absent must degrade, not raise.
    root = tmp_path / "results"
    mdir = root / "multiplex-influenza-a-mp-indo-only_x_influenza-b-mp-indo-only"
    _write(mdir / "multiplex_candidates.json", [{
        "pair_id": "FluA_Rank1_x_FluB_Rank2", "oligos_synthesis": {"flu_a": {}, "flu_b": {}},
        "selection_rank": 1, "selection_passed": True, "selection_level": "L1",
    }])
    n = qc_dossier.build_dossiers(str(root), str(tmp_path / "qc"), top_n=5)
    assert n == 1
    assert (tmp_path / "qc" / "FluA_Rank1_x_FluB_Rank2.html").exists()
