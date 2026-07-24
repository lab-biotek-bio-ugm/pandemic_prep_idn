#!/usr/bin/env python3
"""
Candidate selection: gates first, then ranking.
=============================================================================
Shared by the pairwise checker, the triplex checker and the supplementary
export, so the criterion that decides which oligos get ORDERED exists in exactly
one place. The M11 thresholds previously lived as separate literals in two
checkers and had already drifted apart; this is the same failure mode with a far
more expensive consequence.

Why gates and ranking are separate
----------------------------------
For every criterion, ask whether violating it makes the assay WRONG or merely
WORSE. Wrong belongs in a gate; worse belongs in the ranking. A criterion that is
already gated must NOT also appear in the ranking — that invents a trade-off
which should not exist.

The previous behaviour ranked pairs on multiplex thermodynamics alone and
recommended a pair whose members were both human-specificity grade C, while grade
A candidates with equal coverage sat unused. Coverage — the single number that
decides whether the assay detects circulating strains — was computed by M10 and
read by nobody.

Why the gates are mechanistic, not score-based
----------------------------------------------
PrimerLab's specificity score is an ACCUMULATING penalty: one perfectly
3'-anchored off-target costs 20 points, and twenty trivial ones also cost 20. The
same grade can therefore mean "one catastrophic binding site" or "background
noise", so no grade threshold can rule the first out.

The Class 1 gates below test the mechanism instead. `human_products` counts
off-target amplicons the primer PAIR could actually produce — a single primer
binding somewhere produces nothing on its own; a false positive needs two primers
facing each other within amplicon distance. That cannot be masked by aggregation.

Empirically, over the 100 Influenza A candidates of the first production run, no
grade A or B candidate had a single potential human amplicon, while grade C and
below did. The grade gate happened to be sufficient there — but that is a
correlation in one dataset, not a property of the formula, which is exactly why
the mechanistic gate exists.
"""

from typing import Any, Callable, Dict, List, Optional, Tuple

# Letter grades ordered best to worst; used only for the Class 2 quality floor.
GRADE_ORDER = {"A": 0, "B": 1, "C": 2, "C-": 3, "D": 4, "F": 5}


def grade_rank(grade: Any) -> int:
    """Numeric order for a letter grade. Unknown/missing sorts worst."""
    return GRADE_ORDER.get(str(grade).strip().upper(), 99)


# ---------------------------------------------------------------------------
# Relaxation profiles
# ---------------------------------------------------------------------------
# Class 1 gates are absent from this table on purpose: they are never relaxed.
#
# Each level loosens exactly one Class 2 gate relative to the previous one, so a
# report can state precisely what was given up. Selection walks down the list and
# stops at the first level that yields any candidate, which makes an empty result
# impossible unless every candidate violates a Class 1 gate — and that outcome is
# information worth surfacing, not a failure to paper over.
DEFAULT_LEVELS: List[Dict[str, Any]] = [
    {"name": "L0", "label": "target profile",
     "max_grade": "B", "min_coverage": 95.0, "coverage_within": None, "allow_warning": False},
    {"name": "L1", "label": "coverage judged relative to the best achievable",
     "max_grade": "B", "min_coverage": 95.0, "coverage_within": 2.0, "allow_warning": False},
    {"name": "L2", "label": "absolute coverage floor lowered to 90%",
     "max_grade": "B", "min_coverage": 90.0, "coverage_within": 2.0, "allow_warning": False},
    {"name": "L3", "label": "specificity floor lowered to grade C",
     "max_grade": "C", "min_coverage": 90.0, "coverage_within": 2.0, "allow_warning": False},
    {"name": "L4", "label": "multiplex WARNING accepted",
     "max_grade": "C", "min_coverage": 85.0, "coverage_within": 5.0, "allow_warning": True},
    # Last resort: Class 1 only. Without this the ladder could still return
    # nothing when every candidate merely fails a quality floor, which is the
    # opposite of what the ladder exists for — an operator facing a difficult
    # target needs the least-bad option, loudly labelled, not silence.
    {"name": "L5", "label": "LAST RESORT — every quality floor waived; only the "
                            "non-negotiable safety gates were applied",
     "max_grade": None, "min_coverage": None, "coverage_within": None,
     "allow_warning": True, "last_resort": True},
]

# A hit this tight is a real priming site rather than background: the 3' terminus
# is fully anchored, it is essentially a perfect match, and nearly the whole
# primer is bound. No candidate in the first production run reached it, which is
# the point — the gate costs nothing today and covers the next dataset.
FATAL_HIT = {"max_three_prime_overhang": 0, "max_mismatches": 1, "min_coverage_percent": 90.0}


def levels_from_config(cfg: Optional[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Read the relaxation ladder from config, falling back to DEFAULT_LEVELS.

    Config is the reviewable place for clinical thresholds, but a malformed or
    absent block must not silently drop the ladder — that would mean selecting
    with no quality floors at all. Anything unusable returns the defaults and says so.
    """
    levels = ((cfg or {}).get("multiplex", {}).get("selection", {}) or {}).get("levels")
    if not levels:
        return DEFAULT_LEVELS
    required = {"name", "label", "max_grade", "min_coverage", "allow_warning"}
    for lvl in levels:
        missing = required - set(lvl)
        if missing:
            print(f"[selection] config level {lvl.get('name', '?')} is missing "
                  f"{sorted(missing)}; using built-in defaults instead")
            return DEFAULT_LEVELS
    return levels


class Rejection(Exception):
    """Raised internally to carry the reason a candidate was dropped."""

    def __init__(self, reason: str, gate_class: int):
        super().__init__(reason)
        self.reason = reason
        self.gate_class = gate_class


def _check_class1(member: Dict[str, Any], who: str) -> None:
    """Gates that are never relaxed: violating them makes the assay WRONG."""
    if str(member.get("status", "")).upper() == "FAIL":
        raise Rejection(f"{who}: failed BLAST verification in M09", 1)

    products = member.get("human_products")
    if isinstance(products, (int, float)) and products > 0:
        raise Rejection(
            f"{who}: {int(products)} off-target amplicon(s) possible on human DNA "
            f"— a direct false-positive mechanism", 1)

    virus_off = member.get("virus_offtargets")
    if isinstance(virus_off, (int, float)) and virus_off > 0:
        raise Rejection(
            f"{who}: {int(virus_off)} off-target(s) against another virus "
            f"— risks reporting the wrong pathogen", 1)

    worst = member.get("human_worst_hit") or {}
    if worst:
        if (worst.get("three_prime_overhang", 99) <= FATAL_HIT["max_three_prime_overhang"]
                and worst.get("mismatches", 99) <= FATAL_HIT["max_mismatches"]
                and (worst.get("coverage_percent") or 0) >= FATAL_HIT["min_coverage_percent"]):
            raise Rejection(
                f"{who}: a human off-target binds with the 3' end fully anchored "
                f"({worst.get('mismatches')} mismatch, "
                f"{worst.get('coverage_percent'):.0f}% of the primer bound)", 1)


def _check_class2(member: Dict[str, Any], who: str, level: Dict[str, Any],
                  best_coverage: Optional[float]) -> None:
    """Quality floors. Relaxed step by step when nothing survives."""
    if grade_rank(member.get("human_grade")) > grade_rank(level["max_grade"]):
        raise Rejection(
            f"{who}: human specificity grade {member.get('human_grade')} "
            f"is below the {level['max_grade']} floor", 2)

    cov = member.get("coverage")
    if cov is None:
        raise Rejection(f"{who}: no M10 coverage data — cannot judge sensitivity", 2)

    within = level.get("coverage_within")
    ok_absolute = cov >= level["min_coverage"]
    ok_relative = (within is not None and best_coverage is not None
                   and cov >= best_coverage - within)
    if not (ok_absolute or ok_relative):
        detail = f"{cov:.2f}% < {level['min_coverage']}%"
        if within is not None and best_coverage is not None:
            detail += f" and more than {within} pp below the best achievable ({best_coverage:.2f}%)"
        raise Rejection(f"{who}: coverage {detail}", 2)



def best_coverage_of(candidates: List[Dict[str, Any]]) -> Optional[float]:
    covs = [c.get("coverage") for c in candidates if c.get("coverage") is not None]
    return max(covs) if covs else None


def prefilter(candidates: List[Dict[str, Any]], who: str,
              levels: Optional[List[Dict[str, Any]]] = None,
              min_keep: int = 1) -> Tuple[List[Dict[str, Any]], Optional[Dict[str, Any]], Dict[str, int]]:
    """Apply the per-candidate gates BEFORE combinations are built.

    Specificity grade, human off-target amplicons and coverage are properties of a
    single candidate — they do not depend on what it is paired with. Applying them
    after combining was both wasteful and the reason an arbitrary cut existed at
    all: M11 took the top N by Primer3 penalty purely to keep the pair count down,
    and that cut discarded the highest-coverage candidate in the run.

    Filtering first removes the need for that cut. Pairwise on the survivors is a
    few thousand combinations, which is cheap; the cubic triplex case is handled by
    cap_to_budget() afterwards.

    Each target relaxes independently, because one may be intrinsically harder than
    another and forcing a shared level would penalise the easy one.

    Returns (survivors, level_used, rejection_counts).
    """
    levels = levels or DEFAULT_LEVELS
    best_cov = best_coverage_of(candidates)
    reasons: Dict[str, int] = {}
    # Every level is evaluated, not just up to the first that succeeds, so the
    # report can show what relaxing one more step would have bought. The choice
    # between "one candidate at 98% coverage" and "ten at 91%" belongs to the
    # person ordering the oligos, and they cannot make it without both numbers.
    counts: Dict[str, int] = {}
    chosen = None
    for level in levels:
        keep, reasons = [], {}
        for cand in candidates:
            try:
                _check_class1(cand, who)
                if not level.get("last_resort", False):
                    _check_class2(cand, who, level, best_cov)
            except Rejection as rej:
                reasons[rej.reason] = reasons.get(rej.reason, 0) + 1
            else:
                keep.append(cand)
        counts[level["name"]] = len(keep)
        if chosen is None and len(keep) >= min_keep:
            chosen = (keep, level, reasons)
    if chosen:
        keep, level, reasons = chosen
        level = dict(level, level_counts=counts)
        return keep, level, reasons
    return [], None, reasons


def candidate_sort_key(cand: Dict[str, Any]) -> tuple:
    """Rank one candidate on diagnostic merit: coverage first, then specificity.

    Used only to decide which survivors to keep when a combination budget forces a
    trim. The previous cut ranked by Primer3 penalty — a general oligo-quality
    measure with little bearing on whether the assay detects circulating strains.
    """
    cov = cand.get("coverage")
    spec = cand.get("human_score")
    return (
        -round((cov if cov is not None else 0.0) * 2) / 2,
        -(spec if isinstance(spec, (int, float)) else 0.0),
        cand.get("final_score") if isinstance(cand.get("final_score"), (int, float)) else 1e9,
    )


def cap_to_budget(groups: List[List[Dict[str, Any]]], budget: int) -> Tuple[List[List[Dict[str, Any]]], bool]:
    """Trim each target's survivors until the combination count fits `budget`.

    Only needed for the cubic triplex case: 55 x 45 x 45 is 111,375 combinations,
    each holding a full oligo specification, so the constraint is output size and
    memory rather than raw speed.

    Trimming keeps the best candidates by diagnostic merit, and takes from the
    largest group first so no target is starved. Returns (groups, was_trimmed) so
    the report can say a cut happened rather than hiding it.
    """
    groups = [sorted(g, key=candidate_sort_key) for g in groups]

    def total(gs):
        n = 1
        for g in gs:
            n *= max(len(g), 1)
        return n

    trimmed = False
    while total(groups) > budget and any(len(g) > 1 for g in groups):
        biggest = max(range(len(groups)), key=lambda i: len(groups[i]))
        groups[biggest] = groups[biggest][:-1]
        trimmed = True
    return groups, trimmed


def evaluate_combination(
    members: List[Tuple[str, Dict[str, Any]]],
    combo: Dict[str, Any],
    level: Dict[str, Any],
    best_coverage: Dict[str, float],
    status_key: str,
) -> Tuple[Optional[tuple], Optional[str], Optional[int]]:
    """Apply the gates and, if they pass, return the lexicographic sort key.

    Returns (sort_key, rejection_reason, gate_class). Exactly one of sort_key and
    rejection_reason is None.
    """
    try:
        for who, member in members:
            _check_class1(member, who)

        status = str(combo.get(status_key, ""))
        if status.startswith("VETO"):
            raise Rejection(f"multiplex veto: {status}", 1)

        if not level.get("last_resort", False):
            for who, member in members:
                _check_class2(member, who, level, best_coverage.get(who))

        if status != "PASS" and not level.get("allow_warning", False):
            raise Rejection(f"multiplex status {status or 'unknown'}", 2)

    except Rejection as rej:
        return None, rej.reason, rej.gate_class

    # Lexicographic, not a weighted sum: these quantities are not commensurable,
    # and a weighted sum would let one silently buy off another.
    coverages = [m.get("coverage") or 0.0 for _, m in members]
    spec_scores = [m.get("human_score") for _, m in members
                   if isinstance(m.get("human_score"), (int, float))]
    penalties = [m.get("final_score") for _, m in members
                 if isinstance(m.get("final_score"), (int, float))]

    key = (
        0 if status == "PASS" else 1,
        # Bucketed to 0.5 pp. Without this, 91.17 vs 91.16 — pure noise — would
        # outrank a real difference in specificity further down the key.
        -round(min(coverages) * 2) / 2,
        -(min(spec_scores) if spec_scores else 0.0),
        -float(combo.get("score", 0.0)),
        sum(penalties) if penalties else 0.0,
    )
    return key, None, None


def select(
    combos: List[Dict[str, Any]],
    members_of: Callable[[Dict[str, Any]], List[Tuple[str, Dict[str, Any]]]],
    status_key: str,
    levels: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """Rank combinations, relaxing Class 2 gates only as far as necessary.

    `members_of(combo)` returns [(label, candidate_dict), ...] for that
    combination, where each candidate carries the M09 fields plus `coverage`
    merged in from M10.

    The result always reports which level was used and why everything else was
    dropped, so a reviewer can see what was traded away rather than inferring it.
    """
    levels = levels or DEFAULT_LEVELS

    best_coverage: Dict[str, float] = {}
    for combo in combos:
        for who, member in members_of(combo):
            cov = member.get("coverage")
            if cov is not None and cov > best_coverage.get(who, -1.0):
                best_coverage[who] = cov

    class1_reasons: Dict[str, int] = {}
    last_reasons: Dict[str, int] = {}
    for level in levels:
        ranked, reasons = [], {}
        for combo in combos:
            key, why, gate_class = evaluate_combination(
                members_of(combo), combo, level, best_coverage, status_key)
            if key is None:
                reasons[why] = reasons.get(why, 0) + 1
                if gate_class == 1:
                    class1_reasons[why] = class1_reasons.get(why, 0) + 1
            else:
                ranked.append((key, combo))

        last_reasons = reasons
        if ranked:
            ranked.sort(key=lambda x: x[0])
            return {
                "selected": [c for _, c in ranked],
                "level": level["name"],
                "level_label": level["label"],
                "relaxed": level["name"] != levels[0]["name"],
                "rejections": dict(sorted(reasons.items(), key=lambda x: -x[1])),
                "best_coverage": best_coverage,
                "n_input": len(combos),
                "n_selected": len(ranked),
            }

    # Every level exhausted. This is only reachable when every combination trips a
    # Class 1 gate, which is a finding in itself and must not be hidden.
    # Report the reasons from the most permissive level attempted, not only the
    # Class 1 ones: when nothing passes at all, the operator needs the full
    # picture of what blocked it, and a Class 2 reason may be the whole story
    # (e.g. every candidate missing M10 coverage).
    merged = dict(last_reasons)
    for reason, n in class1_reasons.items():
        merged.setdefault(reason, n)
    return {
        "selected": [],
        "level": None,
        "level_label": "no level produced a candidate",
        "relaxed": True,
        "rejections": dict(sorted(merged.items(), key=lambda x: -x[1])),
        "best_coverage": best_coverage,
        "n_input": len(combos),
        "n_selected": 0,
    }



def annotate(combos: List[Dict[str, Any]], result: Dict[str, Any]) -> None:
    """Stamp the selection verdict onto each combination, in place.

    The verdict has to survive into the JSON, because the supplementary export
    reads that file and previously labelled whatever sat at index 0 as
    "BEST RECOMMENDED" — with no way to know whether it had passed the gates or
    merely come first in a fallback ordering.
    """
    passed = {id(c) for c in result.get("selected", [])}
    for i, combo in enumerate(combos):
        is_sel = id(combo) in passed
        combo["selection_passed"] = is_sel
        combo["selection_level"] = result.get("level")
        combo["selection_level_label"] = result.get("level_label")
        combo["selection_relaxed"] = bool(result.get("relaxed"))
        combo["selection_rank"] = (i + 1) if is_sel else None


def selection_note(result: Dict[str, Any], top_n: int = 8) -> str:
    """Markdown block explaining what the gates did, for the report."""
    lines = ["", "---", "", "## How This Recommendation Was Chosen", ""]

    if not result["selected"]:
        lines += [
            "**No combination passed.** Every candidate violated a gate that is "
            "never relaxed — one that makes the assay wrong rather than merely "
            "worse. Nothing is recommended; the reasons are listed below.", ""]
    else:
        lines += [
            f"Screened {result['n_input']} combinations; {result['n_selected']} "
            f"passed at level **{result['level']}** ({result['level_label']}).", ""]
        if result["relaxed"]:
            lines += [
                "!!! warning \"Quality floors were relaxed\"",
                f"    No combination met the target profile, so the gates were "
                f"loosened stepwise until candidates appeared. The recommendation "
                f"comes from **{result['level']}** — {result['level_label']}. "
                f"Read it as the best available under these data, not as one that "
                f"met the intended standard.", ""]

    lines += [
        "Selection applies gates before ranking. A criterion that makes the assay "
        "*wrong* excludes a candidate outright and never trades against one that "
        "makes it merely *worse*.", "",
        "Gates that are never relaxed: an off-target amplicon possible on human "
        "DNA, an off-target against another virus, a human binding site with the "
        "3' end fully anchored, a dye-channel collision, and a cross-dimer past "
        "the veto limit. These test the mechanism rather than a score, because the "
        "specificity score is an accumulating penalty in which one catastrophic "
        "binding site and twenty trivial ones produce the same grade.", "",
        "Ranking among survivors is lexicographic — multiplex status, then coverage "
        "bucketed to 0.5 pp, then the weakest member's specificity score, then the "
        "compatibility score. Coverage leads not because it outranks specificity "
        "but because specificity is already guaranteed at the gate.", ""]

    if result["rejections"]:
        lines += ["### Why combinations were dropped", "",
                  "| Count | Reason |", "|---|---|"]
        for reason, n in list(result["rejections"].items())[:top_n]:
            lines.append(f"| {n} | {reason} |")
        lines.append("")
    return "\n".join(lines)
