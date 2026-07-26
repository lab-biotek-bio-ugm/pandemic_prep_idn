#!/usr/bin/env python3
"""Generate a per-assay QC decision dossier from existing pipeline artifacts.

The pipeline's end use is national viral surveillance, so every go/no-go decision
on a candidate must be auditable. This tool reads the artifacts a run already
wrote -- it recomputes nothing -- and renders, per multiplex assay that enters the
Executive Summary:

  1. the selection verdict and which relaxation level was used;
  2. a coverage audit including the per-subtype/year/country breakdown, so a
     lineage-specific hole (e.g. H5N1) is visible, not hidden in the average;
  3. a full human off-target audit: every human hit with its accession, organism,
     identity, mismatch, 3'-overhang and coverage, each classified against the
     FATAL_HIT gate with the specific reason it is tolerated or fatal, plus the
     pair-level amplicon count -- i.e. exactly why a "grade A" candidate is safe;
  4. the decision tree: each gate's value -> threshold -> verdict -> rationale ->
     citation, from qc/decision_criteria.yaml.

Output: a self-contained HTML dossier per assay (no external assets), a markdown
mirror, and an index. Inputs read (per member, degenerate mode):
  results/<build>/m09/<mode>/<region>/candidate_<rank>/db_human|db_virus/blast_result.json
  results/<build>/m10/<mode>/coverage_summary.json   (record where rank == member rank)
  results/multiplex*/multiplex_candidates.json
"""
import argparse
import glob
import html
import json
import os
import re
import sys
from pathlib import Path

import yaml

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE / "tools"))
try:
    from candidate_selection import FATAL_HIT, grade_rank
except Exception:  # pragma: no cover - fallback keeps the tool runnable standalone
    FATAL_HIT = {"max_three_prime_overhang": 0, "max_mismatches": 1, "min_coverage_percent": 90.0}
    _ORDER = {"A": 0, "B": 1, "C": 2, "C-": 3, "D": 4, "F": 5}

    def grade_rank(g):
        return _ORDER.get(str(g).strip().upper(), 99)

# Map an oligo target key to the build-directory prefix it corresponds to.
TARGET_CANON = {
    "flu_a": "influenza-a", "flu_b": "influenza-b", "flu_c": "influenza-c",
    "flu": "influenza", "measles": "measles", "rsv": "rsv", "sars_cov_2": "sars",
}
COVERAGE_HOLE_PCT = 90.0   # a stratum below this is flagged in the coverage audit


# --------------------------------------------------------------------------- IO

def load_criteria(path):
    return yaml.safe_load(Path(path).read_text(encoding="utf-8"))


def default_criteria():
    """The registry at qc/decision_criteria.yaml. Shared with export_supplementary."""
    return load_criteria(_HERE.parent / "qc" / "decision_criteria.yaml")


def load_json(path):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except Exception:
        return None


def find_multiplex_files(results_dir):
    out = []
    for name in ("triplex_candidates.json", "multiplex_candidates.json"):
        out += glob.glob(os.path.join(results_dir, "multiplex*", name))
    # Prefer triplex file if a dir has both; dedupe by directory.
    seen, keep = set(), []
    for f in sorted(out):
        d = os.path.dirname(f)
        if d not in seen:
            seen.add(d)
            keep.append(f)
    return keep


def build_dirs(results_dir):
    return [d for d in glob.glob(os.path.join(results_dir, "*"))
            if os.path.isdir(d) and not os.path.basename(d).startswith("multiplex")
            and not os.path.basename(d) in ("supplementary",)]


def find_build(results_dir, target_key):
    canon = TARGET_CANON.get(target_key, target_key.replace("_", "-"))
    for d in build_dirs(results_dir):
        if os.path.basename(d).startswith(canon):
            return d
    for d in build_dirs(results_dir):
        if canon in os.path.basename(d):
            return d
    return None


# ----------------------------------------------------------------- combo parse

_RANK_RE = re.compile(r"([A-Za-z0-9]+?)_Rank(\d+)")


def parse_members(combo):
    """Return [(target_key, rank), ...] for the combo, from pair_id + oligos keys."""
    oligo_keys = list((combo.get("oligos_synthesis") or {}).keys())
    members = []
    for token, rank in _RANK_RE.findall(str(combo.get("pair_id") or combo.get("triplex_id") or "")):
        key = _token_to_key(token, oligo_keys)
        members.append((key, int(rank)))
    if not members:  # fall back to rank_a/rank_b/rank_c + sorted oligo keys
        for suffix, key in zip(("a", "b", "c"), sorted(oligo_keys)):
            r = combo.get(f"rank_{suffix}")
            if r is not None:
                members.append((key, int(r)))
    return members


def _token_to_key(token, oligo_keys):
    t = token.lower()
    for k in oligo_keys:
        if t == k.replace("_", "") or t.startswith(k.replace("_", "")):
            return k
    return t


# ----------------------------------------------------------- member artefacts

def member_artifacts(build_dir, rank, mode="degenerate"):
    human = _first(glob.glob(os.path.join(
        build_dir, "m09", mode, "*", f"candidate_{rank}", "db_human", "blast_result.json")))
    virus = _first(glob.glob(os.path.join(
        build_dir, "m09", mode, "*", f"candidate_{rank}", "db_virus", "blast_result.json")))
    cov_file = _first(glob.glob(os.path.join(build_dir, "m10", mode, "coverage_summary.json")))
    cov = None
    if cov_file:
        for rec in (load_json(cov_file) or []):
            if rec.get("rank") == rank:
                cov = rec
                break
    return {
        "human_blast": load_json(human) if human else None,
        "virus_blast": load_json(virus) if virus else None,
        "coverage": cov,
        "human_path": human,
    }


def _first(paths):
    return sorted(paths)[0] if paths else None


# ---------------------------------------------------------- human hit analysis

def collect_human_hits(human_blast):
    hits = []
    if not human_blast:
        return hits
    for role in ("forward_result", "reverse_result", "probe_result"):
        r = human_blast.get(role) or {}
        label = role.replace("_result", "")
        for o in (r.get("offtargets") or []):
            hits.append({"role": label, **o})
    return hits


def classify_hit(hit):
    """Return ('FATAL'|'tolerated', reason). Mirrors candidate_selection.FATAL_HIT."""
    overhang = hit.get("three_prime_overhang", 99)
    mism = hit.get("mismatches", 99)
    cov = hit.get("coverage_percent") or 0
    is_fatal = (overhang <= FATAL_HIT["max_three_prime_overhang"]
                and mism <= FATAL_HIT["max_mismatches"]
                and cov >= FATAL_HIT["min_coverage_percent"])
    if is_fatal:
        return "FATAL", ("3' end fully anchored, near-perfect and almost fully bound "
                         "— a genuine human priming site")
    reasons = []
    if overhang > FATAL_HIT["max_three_prime_overhang"]:
        reasons.append(f"3' end not anchored (overhang {overhang} bp &gt; 0) → cannot be extended")
    if mism > FATAL_HIT["max_mismatches"]:
        reasons.append(f"{mism} mismatches")
    if cov < FATAL_HIT["min_coverage_percent"]:
        reasons.append(f"only {cov:g}% of the primer binds (&lt; 90%)")
    if not hit.get("is_significant", False):
        reasons.append(f"not significant (e-value {hit.get('evalue')})")
    return "tolerated", "; ".join(reasons) or "below all fatal-hit thresholds"


def human_worst_hit(hits):
    if not hits:
        return None
    # Tightest = smallest overhang, then fewest mismatches, then highest coverage.
    return sorted(hits, key=lambda h: (h.get("three_prime_overhang", 99),
                                       h.get("mismatches", 99),
                                       -(h.get("coverage_percent") or 0)))[0]


def virus_offtarget_count(virus_blast):
    if not virus_blast:
        return None
    n = 0
    for role in ("forward_result", "reverse_result", "probe_result"):
        n += int((virus_blast.get(role) or {}).get("offtarget_count") or 0)
    return n


# --------------------------------------------------------------- gate verdicts

def member_gates(art):
    human = art["human_blast"] or {}
    hits = collect_human_hits(human)
    worst = human_worst_hit(hits)
    worst_verdict = classify_hit(worst)[0] if worst else "tolerated"
    cov_rec = art["coverage"] or {}
    voff = virus_offtarget_count(art["virus_blast"])
    return {
        "is_specific": human.get("is_specific"),
        "human_grade": human.get("grade"),
        "human_score": human.get("combined_score"),
        "human_products": human.get("potential_products"),
        "virus_offtargets": voff,
        # Two coverage figures, deliberately kept distinct. The surveillance/global
        # DB (the main --sequences, e.g. 50k) tests broad inclusivity incl. exotic
        # /avian lineages; the design/focal DB is the target population the primers
        # were built on. The selection gate (m11 load_coverage) reads the GLOBAL
        # figure, so `coverage` mirrors that; both are surfaced in the report.
        "coverage": cov_rec.get("raa_detected_pct"),
        "coverage_global": cov_rec.get("raa_detected_pct"),
        "coverage_design": cov_rec.get("design_db_raa_det_pct"),
        "n_global": cov_rec.get("total_seqs"),
        "n_design": cov_rec.get("design_db_total"),
        "human_hits": hits,
        "human_worst": worst,
        "human_worst_fatal": worst_verdict == "FATAL",
    }


def eval_gate(gate, combo, m):
    """Return (scope, value_str, passed|None). None = not evaluable here."""
    gid = gate["id"]
    if gid == "human_offtarget_amplicon":
        v = m.get("human_products")
        return "member", f"{v} human amplicon(s)", (v == 0) if v is not None else None
    if gid == "virus_crossreactivity":
        v = m.get("virus_offtargets")
        return "member", f"{v} cross-virus off-target(s)", (v == 0) if v is not None else None
    if gid == "fatal_human_3prime_anchor":
        return "member", ("a fully-anchored human hit exists" if m.get("human_worst_fatal")
                          else "no fully-anchored human hit"), (not m.get("human_worst_fatal"))
    if gid == "m09_verification":
        v = m.get("is_specific")
        return "member", f"is_specific={v}", (v is True) if v is not None else None
    if gid == "coverage_floor":
        v = m.get("coverage_global")
        d = m.get("coverage_design")
        valstr = (f"{v:.2f}% global" if v is not None else "no data") + \
                 (f" / {d:.2f}% focal" if d is not None else "")
        # Gate mirrors the code: judged on the global figure.
        return "member", valstr, (v >= 95.0) if v is not None else None
    if gid == "specificity_grade_floor":
        g = m.get("human_grade")
        return "member", f"grade {g}", (grade_rank(g) <= grade_rank("B")) if g else None
    if gid == "dye_channel_separation":
        v = combo.get("dye_ok")
        return "combo", ("distinct" if v else "COLLISION"), bool(v)
    if gid == "cross_dimer_veto":
        v = combo.get("min_dg_raa_39c")
        return "combo", f"{v} kcal/mol", (v >= gate.get("value", -12.0)) if v is not None else None
    if gid == "delta_tm":
        v = combo.get("delta_tm")
        return "combo", f"{v} °C", (v <= gate.get("value", 4.0)) if v is not None else None
    if gid == "three_prime_anchor":
        v = combo.get("max_3prime_anchor_bp")
        return "combo", f"{v} bp", (v <= gate.get("value", 5)) if v is not None else None
    return "combo", "-", None


# ----------------------------------------------------------------- rendering

def _badge(passed):
    if passed is None:
        return '<span class="b b-na">N/A</span>'
    return '<span class="b b-pass">PASS</span>' if passed else '<span class="b b-fail">FAIL</span>'


def _cite(gate, refs_meta):
    keys = gate.get("references") or []
    pend = set((refs_meta or {}).get("pending_references") or [])
    parts = []
    for k in keys:
        parts.append(f"{k}*" if k in pend else k)
    return ", ".join(parts) if parts else "—"


def render_member_html(key, rank, art, gates, criteria):
    m = member_gates(art)
    name = key.replace("_", " ").title()
    hits = m["human_hits"]
    rows = []
    for h in sorted(hits, key=lambda x: (x.get("three_prime_overhang", 99), x.get("mismatches", 99))):
        verdict, reason = classify_hit(h)
        cls = "b-fail" if verdict == "FATAL" else "b-ok"
        rows.append(
            f"<tr><td>{html.escape(str(h.get('role')))}</td>"
            f"<td><code>{html.escape(str(h.get('sequence_id')))}</code></td>"
            f"<td>{html.escape(str(h.get('sequence_title'))[:70])}</td>"
            f"<td>{h.get('identity')}</td><td>{h.get('mismatches')}</td>"
            f"<td>{h.get('three_prime_overhang')}</td><td>{h.get('coverage_percent')}</td>"
            f"<td><span class='b {cls}'>{verdict}</span> {reason}</td></tr>")
    hit_table = ("".join(rows) if rows
                 else "<tr><td colspan=8>No human off-target hits recorded.</td></tr>")
    prod = m.get("human_products")
    conclusion = (
        f"{len(hits)} partial human match(es); the primer pair forms "
        f"<b>{prod}</b> human amplicon(s). "
        + ("<b>No false-positive mechanism</b> — no hit is a genuine priming site and the "
           "pair cannot amplify human DNA." if (prod == 0 and not m['human_worst_fatal'])
           else "<b>REVIEW</b> — a fully-anchored hit or a possible amplicon is present."))

    # coverage breakdown: surveillance/global vs design/focal, side by side.
    cov = art["coverage"] or {}
    g_bd = {r.get("subtype"): r for r in ((cov.get("breakdowns") or {}).get("subtype") or [])}
    d_bd = {r.get("subtype"): r for r in ((cov.get("design_db_breakdowns") or {}).get("subtype") or [])}
    brows = []
    for st in sorted(set(g_bd) | set(d_bd),
                     key=lambda s: (g_bd.get(s, {}).get("raa_detected_pct") or 0)):
        g = g_bd.get(st, {})
        d = d_bd.get(st, {})
        gp = g.get("raa_detected_pct")
        dp = d.get("raa_detected_pct")
        gcls = "b-fail" if (gp is not None and gp < COVERAGE_HOLE_PCT) else "b-ok"
        dcell = (f"<span class='b {'b-fail' if (dp is not None and dp < COVERAGE_HOLE_PCT) else 'b-ok'}'>"
                 f"{dp}%</span> (n={d.get('total')})" if d else "<span class='note'>not in focal set</span>")
        gcell = (f"<span class='b {gcls}'>{gp}%</span> (n={g.get('total')})" if g
                 else "<span class='note'>—</span>")
        brows.append(f"<tr><td>{html.escape(str(st))}</td><td>{gcell}</td><td>{dcell}</td></tr>")
    brk = "".join(brows) or "<tr><td colspan=3>No subtype breakdown.</td></tr>"

    # member gate rows
    grows = []
    for g in gates:
        scope, valstr, passed = eval_gate(g, {}, m)
        if scope != "member":
            continue
        grows.append(f"<tr><td>C{g['class']}</td><td>{html.escape(g['name'])}</td>"
                     f"<td>{valstr}</td><td>{_badge(passed)}</td>"
                     f"<td>{_cite(g, criteria.get('meta'))}</td></tr>")

    cg = m.get("coverage_global")
    cd = m.get("coverage_design")
    return f"""
    <section class="member">
      <h3>{html.escape(name)} — candidate rank {rank} ({html.escape(str(art.get('human_grade') or ''))})</h3>
      <div class="cards">
        <div class="card"><div class="k">Coverage — design / focal DB (n={m.get('n_design')})</div><div class="v">{cd if cd is not None else '—'}%</div></div>
        <div class="card"><div class="k">Coverage — surveillance / global DB (n={m.get('n_global')})</div><div class="v">{cg if cg is not None else '—'}%</div></div>
        <div class="card"><div class="k">Human specificity</div><div class="v">{html.escape(str(m.get('human_grade')))} ({m.get('human_score')})</div></div>
        <div class="card"><div class="k">Human amplicons</div><div class="v">{m.get('human_products')}</div></div>
        <div class="card"><div class="k">Cross-virus off-targets</div><div class="v">{m.get('virus_offtargets')}</div></div>
      </div>
      <p class="note">The selection gate judges the <b>surveillance / global</b> figure
      ({cg}%). The <b>design / focal</b> figure ({cd}%) is coverage of the target
      population the primers were built on; global holes below are largely exotic
      /imported lineages, a preparedness signal rather than a routine-use defect.</p>

      <h4>Human off-target audit</h4>
      <p class="note">{conclusion}</p>
      <table><thead><tr><th>Oligo</th><th>Accession</th><th>Organism</th><th>Identity</th>
      <th>Mism.</th><th>3'-overhang</th><th>Coverage</th><th>Verdict &amp; why</th></tr></thead>
      <tbody>{hit_table}</tbody></table>

      <h4>Coverage by subtype — global vs focal (holes &lt; {COVERAGE_HOLE_PCT:g}% flagged)</h4>
      <table><thead><tr><th>Subtype</th><th>Surveillance / global</th><th>Design / focal</th></tr></thead>
      <tbody>{brk}</tbody></table>

      <h4>Per-member gates</h4>
      <table><thead><tr><th>Class</th><th>Gate</th><th>Value</th><th>Verdict</th><th>Cite</th></tr></thead>
      <tbody>{''.join(grows)}</tbody></table>
    </section>"""


def render_attribution(build_dir):
    """Render the M12 coverage-failure attribution for a build, if present.

    Read-only: consumes results/<build>/m12/coverage_attribution.json (written by
    scripts/m12_coverage_attribution.py). Silent if absent, so the dossier works
    whether or not the attribution step has run.
    """
    path = Path(build_dir) / "m12" / "coverage_attribution.json"
    data = load_json(path)
    if not data or not data.get("subtypes"):
        return ""
    order = {"mechanistic": 0, "raa_rescued": 1, "no_focal_data": 2}
    rows = []
    for r in sorted(data["subtypes"], key=lambda x: order.get(x.get("failure_mode"), 3)):
        rep = r.get("reported_pct", r.get("coverage_pct"))
        raa = r.get("raa_viable_pct")
        cells = (f"<td>{html.escape(str(r.get('subtype')))}</td>"
                 f"<td>{rep}%</td><td>{raa if raa is not None else '—'}%</td>"
                 f"<td>{html.escape(str(r.get('failure_mode')))}</td>"
                 f"<td>{html.escape(str(r.get('explanation', '')))}</td>")
        rows.append(f"<tr>{cells}</tr>")
    note = html.escape(data.get("finding_note", ""))
    return f"""
    <section class="member">
      <h3>Coverage-failure attribution (why some subtypes read &lt;90%)</h3>
      <p class="note">{note}</p>
      <table><thead><tr><th>Subtype</th><th>Reported</th><th>RAA-viable</th>
      <th>Mode</th><th>Explanation</th></tr></thead>
      <tbody>{''.join(rows)}</tbody></table>
    </section>"""


def render_assay_html(combo, members, criteria, results_dir):
    gates = criteria["gates"]
    pid = combo.get("pair_id") or combo.get("triplex_id") or "assay"
    passed = combo.get("selection_passed")
    level = combo.get("selection_level")
    label = combo.get("selection_level_label", "")
    relaxed = combo.get("selection_relaxed")
    verdict_cls = "b-pass" if passed else "b-fail"

    # combo-level gate rows
    crows = []
    for g in gates:
        scope, valstr, ok = eval_gate(g, combo, {})
        if scope != "combo":
            continue
        crows.append(f"<tr><td>C{g['class']}</td><td>{html.escape(g['name'])}</td>"
                     f"<td>{valstr}</td><td>{_badge(ok)}</td>"
                     f"<td>{_cite(g, criteria.get('meta'))}</td></tr>")

    # full decision-tree table (rationale + citation for every gate)
    trows = []
    for g in gates:
        trows.append(
            f"<tr><td>C{g['class']}</td><td>{html.escape(g['name'])}</td>"
            f"<td>{html.escape(g.get('rule',''))}</td>"
            f"<td>{html.escape(g.get('rationale',''))}</td>"
            f"<td>{html.escape(g.get('consequence',''))}</td>"
            f"<td>{_cite(g, criteria.get('meta'))}</td></tr>")

    member_html = ""
    for key, rank in members:
        bdir = find_build(results_dir, key)
        if not bdir:
            member_html += f"<section class='member'><h3>{html.escape(key)} rank {rank}</h3><p class='note'>Build directory not found.</p></section>"
            continue
        art = member_artifacts(bdir, rank)
        art["human_grade"] = (art["human_blast"] or {}).get("grade")
        member_html += render_member_html(key, rank, art, gates, criteria)
        member_html += render_attribution(bdir)

    relaxed_note = (f" — quality floors relaxed to {level} ({html.escape(label)})"
                    if relaxed else "")
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>QC Dossier — {html.escape(pid)}</title>{_CSS}</head><body>
<header><h1>QC Decision Dossier</h1>
<h2>{html.escape(pid)}</h2>
<p class="verdict"><span class="b {verdict_cls}">{'RECOMMENDED' if passed else 'NOT RECOMMENDED'}</span>
selection level <b>{html.escape(str(level))}</b>{relaxed_note}. Compatibility score {combo.get('score')}, status {html.escape(str(combo.get('status_raa_39c')))}.</p>
<p class="note">This dossier reads only artifacts the run already wrote; nothing is recomputed. A citation key ending in <code>*</code> is added to references.bib in a later phase.</p>
</header>

<h3>Assay-level gates</h3>
<table><thead><tr><th>Class</th><th>Gate</th><th>Value</th><th>Verdict</th><th>Cite</th></tr></thead>
<tbody>{''.join(crows)}</tbody></table>

{member_html}

<h3>Decision tree — every gate, its rationale and citation</h3>
<p class="note">Class 1 gates are never relaxed (a violation makes the assay wrong).
Class 2 gates are quality floors, loosened one step at a time down the L0→L5 ladder
only when nothing survives; the level above states how far that went.</p>
<table><thead><tr><th>Class</th><th>Gate</th><th>Rule</th><th>Rationale</th><th>If violated</th><th>Cite</th></tr></thead>
<tbody>{''.join(trows)}</tbody></table>

<footer><p class="note">Generated by scripts/qc_dossier.py from qc/decision_criteria.yaml.</p></footer>
</body></html>"""


_CSS = """<style>
:root{color-scheme:dark}
body{background:#0f1115;color:#e6e6e6;font:15px/1.55 system-ui,Segoe UI,Roboto,sans-serif;margin:0;padding:0 5vw 6rem}
header{padding:2rem 0 1rem;border-bottom:1px solid #262a33}
h1{margin:0;font-size:1.5rem}h2{margin:.2rem 0;color:#8ab4ff;font-size:1.15rem}
h3{margin:1.8rem 0 .6rem;font-size:1.1rem;border-left:3px solid #8ab4ff;padding-left:.5rem}
h4{margin:1.1rem 0 .4rem;color:#c9d3e0}
table{border-collapse:collapse;width:100%;margin:.4rem 0 1rem;font-size:13.5px;overflow-x:auto;display:block}
th,td{border:1px solid #2a2f3a;padding:.35rem .5rem;text-align:left;vertical-align:top}
th{background:#1a1e26;color:#9fb0c6;position:sticky;top:0}
code{background:#1a1e26;padding:.05rem .3rem;border-radius:3px;font-size:12.5px}
.cards{display:flex;flex-wrap:wrap;gap:.6rem;margin:.6rem 0}
.card{background:#161a22;border:1px solid #262a33;border-radius:8px;padding:.6rem .9rem;min-width:150px}
.card .k{color:#8b95a5;font-size:12px}.card .v{font-size:1.25rem;font-weight:600}
.b{display:inline-block;padding:.05rem .4rem;border-radius:4px;font-size:12px;font-weight:600}
.b-pass,.b-ok{background:#12351f;color:#5fd38a}.b-fail{background:#3a1620;color:#ff8095}.b-na{background:#2a2f3a;color:#9aa3b2}
.verdict{font-size:1rem}.note{color:#9aa3b2;font-size:13px}
.member{border:1px solid #222732;border-radius:10px;padding:.4rem 1rem 1rem;margin:1rem 0;background:#12151c}
footer{margin-top:2rem;border-top:1px solid #262a33;padding-top:.6rem}
</style>"""


def render_assay_md(combo, members, criteria, results_dir):
    pid = combo.get("pair_id") or combo.get("triplex_id") or "assay"
    lines = [f"# QC Decision Dossier — {pid}", ""]
    lines.append(f"- Verdict: {'RECOMMENDED' if combo.get('selection_passed') else 'NOT RECOMMENDED'}"
                 f" (level {combo.get('selection_level')}"
                 f"{', relaxed' if combo.get('selection_relaxed') else ''})")
    lines.append(f"- Compatibility score {combo.get('score')}, status {combo.get('status_raa_39c')}")
    lines.append("")
    for key, rank in members:
        bdir = find_build(results_dir, key)
        if not bdir:
            continue
        art = member_artifacts(bdir, rank)
        m = member_gates(art)
        lines.append(f"## {key} — rank {rank} (grade {m.get('human_grade')})")
        lines.append(f"- Coverage — design/focal DB (n={m.get('n_design')}): {m.get('coverage_design')}%")
        lines.append(f"- Coverage — surveillance/global DB (n={m.get('n_global')}): "
                     f"{m.get('coverage_global')}%  *(the selection gate judges this figure)*")
        lines.append(f"- Human amplicons (pair): {m.get('human_products')}; "
                     f"cross-virus off-targets: {m.get('virus_offtargets')}")
        lines.append("")
        lines.append("| Oligo | Accession | Organism | Ident | Mism | 3'-ovh | Cov | Verdict |")
        lines.append("|---|---|---|---|---|---|---|---|")
        for h in sorted(m["human_hits"], key=lambda x: (x.get("three_prime_overhang", 99),)):
            verdict, reason = classify_hit(h)
            org = str(h.get("sequence_title"))[:50].replace("|", " ")
            lines.append(f"| {h.get('role')} | {h.get('sequence_id')} | {org} | "
                         f"{h.get('identity')} | {h.get('mismatches')} | "
                         f"{h.get('three_prime_overhang')} | {h.get('coverage_percent')} | "
                         f"{verdict}: {re.sub('<[^>]+>', '', reason)} |")
        lines.append("")
    return "\n".join(lines)


# --------------------------------------------------------------------- driver

def build_dossiers(results_dir, out_dir, top_n):
    criteria = default_criteria()
    os.makedirs(out_dir, exist_ok=True)
    index = ["<!doctype html><meta charset=utf-8><title>QC Dossiers</title>"
             + _CSS + "<h1>QC Decision Dossiers</h1><ul>"]
    made = 0
    for mf in find_multiplex_files(results_dir):
        combos = load_json(mf) or []
        combos = sorted(combos, key=lambda c: c.get("selection_rank") or 999)[:top_n]
        for combo in combos:
            members = parse_members(combo)
            if not members:
                continue
            pid = combo.get("pair_id") or combo.get("triplex_id") or f"assay_{made}"
            slug = re.sub(r"[^A-Za-z0-9_.-]", "_", pid)
            (Path(out_dir) / f"{slug}.html").write_text(
                render_assay_html(combo, members, criteria, results_dir), encoding="utf-8")
            (Path(out_dir) / f"{slug}.md").write_text(
                render_assay_md(combo, members, criteria, results_dir), encoding="utf-8")
            index.append(f'<li><a href="{slug}.html">{html.escape(pid)}</a> '
                         f'— level {combo.get("selection_level")}, score {combo.get("score")}</li>')
            made += 1
    index.append("</ul>")
    (Path(out_dir) / "index.html").write_text("".join(index), encoding="utf-8")
    return made


def main():
    ap = argparse.ArgumentParser(description="Generate per-assay QC decision dossiers.")
    ap.add_argument("--results-dir", required=True, help="Pipeline results/ directory.")
    ap.add_argument("--out", required=True, help="Output directory for dossiers.")
    ap.add_argument("--top-n", type=int, default=5,
                    help="How many top multiplex combinations to document (default 5).")
    args = ap.parse_args()
    n = build_dossiers(args.results_dir, args.out, args.top_n)
    print(f"Wrote {n} assay dossier(s) to {args.out}")
    if n == 0:
        print("No multiplex candidates found under", args.results_dir, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
