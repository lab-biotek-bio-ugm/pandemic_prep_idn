"""Tests for the gate-then-rank selection that decides which oligos get ordered.

The behaviour being pinned here was found in the first production run: M11 ranked
pairs on multiplex thermodynamics alone and recommended a pair whose members were
both human-specificity grade C, while grade A candidates with equal coverage went
unused — and M10's coverage, the number that decides whether the assay detects
circulating strains, was read by nobody at all.
"""
import glob

import pytest
import yaml

from conftest import load_module  # noqa: E402

sel = load_module("scripts/tools/candidate_selection.py", "candidate_selection")


def member(**kw):
    """A candidate that passes every gate unless a field is overridden."""
    base = dict(status="PASS", human_grade="A", human_score=95.0, human_products=0,
                virus_offtargets=0, coverage=99.0, final_score=5.0,
                human_worst_hit={"three_prime_overhang": 3, "mismatches": 2,
                                 "coverage_percent": 60.0})
    base.update(kw)
    return base


def combo(score=70.0, status="PASS", a=None, b=None, tag="x"):
    return {"pair_id": tag, "score": score, "status_raa_39c": status,
            "_a": a or member(), "_b": b or member()}


def run(combos, levels=None):
    return sel.select(combos, lambda c: [("A", c["_a"]), ("B", c["_b"])],
                      status_key="status_raa_39c", levels=levels)


# --- Class 1: never relaxed --------------------------------------------------

@pytest.mark.parametrize("field,value,fragment", [
    ("human_products", 2, "off-target amplicon"),
    ("virus_offtargets", 1, "another virus"),
    ("status", "FAIL", "failed BLAST"),
])
def test_class1_gates_survive_every_relaxation_level(field, value, fragment):
    """These make the assay WRONG, not merely worse. No level may admit them —
    otherwise the ladder that prevents an empty result would do so by shipping a
    candidate that produces false positives."""
    c = combo(a=member(**{field: value}))
    result = run([c])
    assert result["selected"] == []
    assert result["level"] is None
    assert any(fragment in r for r in result["rejections"])


def test_a_fully_anchored_human_hit_is_rejected_regardless_of_grade():
    """The specificity score is an accumulating penalty, so one catastrophic
    binding site and twenty trivial ones yield the same grade. The gate must test
    the hit itself, not the aggregate."""
    fatal = member(human_grade="A", human_score=99.0,
                   human_worst_hit={"three_prime_overhang": 0, "mismatches": 0,
                                    "coverage_percent": 98.0})
    result = run([combo(a=fatal)])
    assert result["selected"] == []
    assert any("3' end fully anchored" in r for r in result["rejections"])


def test_dye_collision_veto_is_class1():
    result = run([combo(status="VETO (Dye Channel Collision)")])
    assert result["selected"] == []


# --- The relaxation ladder ---------------------------------------------------

def test_ladder_stops_at_the_first_level_that_yields_candidates():
    """A candidate that only misses the L0 coverage floor must be found at L1 or
    L2, not dropped and not rescued from the bottom of the ladder."""
    c = combo(a=member(coverage=91.0), b=member(coverage=99.0))
    result = run([c])
    assert result["selected"], result["rejections"]
    assert result["level"] in ("L1", "L2")
    assert result["relaxed"] is True


def test_no_relaxation_reported_when_the_target_profile_is_met():
    result = run([combo(a=member(coverage=99.5), b=member(coverage=99.5))])
    assert result["level"] == "L0"
    assert result["relaxed"] is False


def test_empty_result_only_when_every_candidate_trips_class1():
    """The user's main worry: over-strict gates returning nothing. That can only
    happen for Class 1 violations, and then it is a finding worth surfacing."""
    only_class2 = combo(a=member(human_grade="F", coverage=10.0))
    res = run([only_class2])
    assert res["selected"], "Class 2 alone must never empty the result"
    assert res["level"] == "L5", "a candidate failing only quality floors must " \
                                 "surface at the last-resort level, clearly labelled"

    class1 = combo(a=member(human_products=3))
    assert run([class1])["selected"] == []


def test_relative_coverage_keeps_options_for_a_difficult_target():
    """A fixed floor can wipe out every candidate for an intrinsically hard
    target. Judging coverage against the best achievable keeps the best available
    rather than returning nothing."""
    hard = [combo(tag=f"c{i}", a=member(coverage=cov), b=member(coverage=cov))
            for i, cov in enumerate((70.0, 69.0, 68.0))]
    result = run(hard)
    assert result["selected"], "a hard target must still yield the best available"
    assert result["selected"][0]["_a"]["coverage"] == 70.0


# --- Ranking among survivors -------------------------------------------------

def test_pass_outranks_warning_even_at_a_lower_score():
    """The observed defect: a WARNING pair sat at rank 1 on a 2.6-point score
    edge while twelve PASS pairs existed."""
    warn = combo(tag="warn", score=69.1, status="WARNING")
    ok = combo(tag="pass", score=66.5, status="PASS")
    result = run([warn, ok], levels=[dict(sel.DEFAULT_LEVELS[-1], allow_warning=True)])
    assert result["selected"][0]["pair_id"] == "pass"


def test_coverage_outranks_the_compatibility_score():
    low_cov = combo(tag="lowcov", score=99.0,
                    a=member(coverage=91.0), b=member(coverage=91.0))
    high_cov = combo(tag="highcov", score=60.0,
                     a=member(coverage=99.0), b=member(coverage=99.0))
    result = run([low_cov, high_cov])
    assert result["selected"][0]["pair_id"] == "highcov"


def test_coverage_is_bucketed_so_noise_does_not_decide():
    """91.17 vs 91.16 is noise. Left unbucketed it would outrank a real
    difference in specificity further down the key."""
    noisy_better_cov = combo(tag="noisy", score=70.0,
                             a=member(coverage=91.17, human_score=80.0),
                             b=member(coverage=91.17, human_score=80.0))
    better_spec = combo(tag="spec", score=70.0,
                        a=member(coverage=91.16, human_score=95.0),
                        b=member(coverage=91.16, human_score=95.0))
    result = run([noisy_better_cov, better_spec])
    assert result["selected"][0]["pair_id"] == "spec"


def test_missing_coverage_rejects_rather_than_passes():
    """M10 output absent must not read as 'fine'. Sensitivity unknown is a reason
    to exclude, not to assume."""
    # Checked at a level that still has quality floors; the last-resort level
    # deliberately waives them and would admit the candidate with a loud label.
    normal = [lvl for lvl in sel.DEFAULT_LEVELS if not lvl.get("last_resort")]
    result = run([combo(a=member(coverage=None))], levels=normal)
    assert result["selected"] == []
    assert any("no M10 coverage" in r for r in result["rejections"])


# --- Wiring ------------------------------------------------------------------

def test_config_levels_match_the_module_defaults():
    """Config is the reviewable home for clinical thresholds; the module default
    must not silently disagree with it."""
    cfg = yaml.safe_load(open("config/config.yaml"))
    levels = sel.levels_from_config(cfg)
    assert len(levels) == len(sel.DEFAULT_LEVELS)
    for got, want in zip(levels, sel.DEFAULT_LEVELS):
        for key in ("name", "max_grade", "min_coverage", "allow_warning"):
            assert got[key] == want[key], f"{got['name']}.{key}"


def test_malformed_config_falls_back_instead_of_dropping_the_gates():
    broken = {"multiplex": {"selection": {"levels": [{"name": "L0"}]}}}
    assert sel.levels_from_config(broken) is sel.DEFAULT_LEVELS
    assert sel.levels_from_config({}) is sel.DEFAULT_LEVELS


def test_every_config_defines_the_selection_ladder():
    for path in glob.glob("config/config*.yaml"):
        levels = (yaml.safe_load(open(path)).get("multiplex", {})
                  .get("selection", {}).get("levels"))
        assert levels, f"{path} has no multiplex.selection.levels"
        assert levels[-1]["allow_warning"] is True, \
            f"{path}: the last level must accept WARNING, or a valid-but-flagged " \
            f"combination could yield an empty result"


def test_both_checkers_use_the_shared_selection():
    for path in ("scripts/m11_multiplex_checker.py", "scripts/m11_triplex_checker.py"):
        src = open(path).read()
        assert "candidate_selection import" in src and "select" in src, path
        assert 'results.sort(key=lambda x: x["score"], reverse=True)' in src, \
            f"{path} must keep the score sort only as the no-selection fallback"


def test_orchestrator_passes_coverage_to_the_checkers():
    src = open("scripts/run_multiplex_batch.py").read()
    assert "def find_coverage(" in src
    assert '"--coverage-a"' in src and '"--coverage-c"' in src


def test_m09_forwards_the_mechanistic_offtarget_fields():
    src = open("scripts/m09_blast_verification.py").read()
    for field in ("human_products", "human_worst_hit"):
        assert f'"{field}":' in src, f"M09 does not forward {field}"


# --- The verdict must reach the ordering sheet -------------------------------

def test_annotate_marks_pass_fail_and_rank():
    """The export reads the JSON, not the report. Without the verdict in the file
    it labelled whatever sat at index 0 as BEST RECOMMENDED, with no way to know
    whether it had passed the gates or merely sorted first in a fallback."""
    good = combo(tag="good", a=member(coverage=99.5), b=member(coverage=99.5))
    bad = combo(tag="bad", a=member(human_products=2))
    combos = [good, bad]
    result = run(combos)
    sel.annotate(combos, result)

    assert good["selection_passed"] is True
    assert good["selection_rank"] == 1
    assert bad["selection_passed"] is False
    assert bad["selection_rank"] is None
    assert good["selection_level"] == result["level"]


def test_export_labels_come_from_the_verdict_not_the_position():
    src = open("scripts/export_supplementary.py").read()
    assert 'combo.get("selection_passed"' in src, \
        "the export still labels by list position"
    assert "NOT RECOMMENDED" in src
    assert "selection_relaxed" in src, \
        "a recommendation made under relaxed floors must say so on the sheet"


def test_export_does_not_call_the_primer3_penalty_a_score():
    """final_score is a PENALTY: lower is better. Labelling it 'Score' made the
    best candidate (6.0) look like the worst on the overview sheet."""
    src = open("scripts/export_supplementary.py").read()
    assert "Rank 1 Score:" not in src
    assert "lower is better" in src


def test_on_target_virus_hits_are_not_called_pathogen_hits():
    """`Virus_Pathogen_Hits` counts hits to the TARGET pathogen — 200 of them is
    the desired outcome — but reads as hits to other pathogens."""
    src = open("scripts/m09_blast_verification.py").read()
    assert "Virus_Pathogen_Hits" not in src
    assert "Virus_OnTarget_Hits" in src


def test_both_checkers_stamp_the_verdict_before_writing_json():
    for path in ("scripts/m11_multiplex_checker.py", "scripts/m11_triplex_checker.py"):
        src = open(path).read()
        assert "annotate(results, selection)" in src, path
        assert src.index("annotate(results, selection)") < src.index("json.dump(results"), \
            f"{path}: the verdict must be stamped before the JSON is written"


# --- Per-candidate gates run before combination ------------------------------

def test_prefilter_removes_candidates_before_they_are_paired():
    """Specificity, human amplicons and coverage are properties of ONE candidate.
    Applying them after combining wasted thermodynamics on pairs that could never
    be recommended — and was the only reason a cut by Primer3 penalty existed,
    which is what discarded the run's highest-coverage candidate."""
    cands = [member(coverage=99.0, human_grade="A"),
             member(coverage=99.0, human_products=2),
             member(coverage=99.0, human_grade="F")]
    kept, level, reasons = sel.prefilter(cands, "T")
    assert len(kept) == 1
    assert level["name"] == "L0"
    assert any("off-target amplicon" in r for r in reasons)


def test_prefilter_reports_what_each_level_would_admit():
    """The choice between one candidate at 98% coverage and ten at 91% belongs to
    whoever orders the oligos, and cannot be made without both numbers."""
    cands = [member(coverage=98.0)] + [member(coverage=91.0) for _ in range(9)]
    kept, level, _ = sel.prefilter(cands, "T")
    counts = level["level_counts"]
    assert counts["L0"] == 1
    assert counts["L2"] == 10
    assert len(kept) == 1, "the strictest level that yields anything must win"


def test_prefilter_relaxes_per_target_independently():
    easy = [member(coverage=99.0) for _ in range(3)]
    hard = [member(coverage=60.0) for _ in range(3)]
    _, lvl_easy, _ = sel.prefilter(easy, "easy")
    _, lvl_hard, _ = sel.prefilter(hard, "hard")
    assert lvl_easy["name"] == "L0"
    assert lvl_hard["name"] != "L0", \
        "a hard target must relax without dragging an easy one down with it"


def test_budget_trims_by_diagnostic_merit_not_primer3_penalty():
    """The previous cut ranked by Primer3 penalty. When a budget forces a trim it
    must keep the candidates that matter diagnostically."""
    group = [member(coverage=99.0, human_score=95.0, final_score=99.0),   # best assay
             member(coverage=80.0, human_score=95.0, final_score=1.0)]    # best penalty
    (trimmed,), was_trimmed = sel.cap_to_budget([group], budget=1)
    assert was_trimmed is True
    assert trimmed[0]["coverage"] == 99.0


def test_budget_never_starves_a_target_completely():
    groups = [[member() for _ in range(20)] for _ in range(3)]
    capped, was_trimmed = sel.cap_to_budget(groups, budget=8)
    assert was_trimmed is True
    assert all(len(g) >= 1 for g in capped)
    n = 1
    for g in capped:
        n *= len(g)
    assert n <= 8


def test_checkers_prefilter_before_the_combination_loop():
    for path in ("scripts/m11_multiplex_checker.py", "scripts/m11_triplex_checker.py"):
        src = open(path).read()
        assert "prefilter(" in src, path
        assert "cap_to_budget(" in src, path
        # The gates must be applied before the loop that builds combinations.
        assert src.index("prefilter(") < src.index("for cand_a in cands_a"), \
            f"{path}: gates must run before the combination loop, not after"


def test_m10_scopes_by_specificity_not_by_rank():
    """M10's coverage budget was spent on the top N by Primer3 penalty, so the
    highest-coverage candidate barely made the cut while worse ones consumed it."""
    src = open("scripts/m10_insilico_pcr.py").read()
    assert "ELIGIBLE_GRADES" in src
    assert "eligible = [c for c in candidates" in src
    assert "candidates[:args.num_candidates]" in src, \
        "the rank cut must remain only as the no-eligible-candidate fallback"


def test_m10_ships_sequence_metadata_once_per_worker():
    """Passing it per task pickled ~50k rows for every candidate, so the cost grew
    with the number of candidates — penalising the decision to measure more."""
    src = open("scripts/m10_insilico_pcr.py").read()
    assert "initializer=_init_worker" in src
    assert "_WORKER_SEQ_MAPPINGS" in src
