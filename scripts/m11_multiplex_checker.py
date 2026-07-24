#!/usr/bin/env python3
"""
Module 11: Multiplex Assay Compatibility Screening Engine (Influenza A + B MP)
Evaluates pairwise cross-dimer thermodynamics (Primer3), spectral dye separation,
and kinetics alignment for dual-target assays at 39.0°C (RAA) and 55.0°C (qPCR).
Outputs full sequence synthesis specs and cross-dimer details.

IUPAC expansion, Tm calculation, and heterodimerisation logic are provided by
the shared `scripts/tools/iupac_utils.py` module.
"""

import os
import sys
import json
import argparse
from pathlib import Path

# Allow importing shared utilities from scripts/tools/
sys.path.insert(0, str(Path(__file__).resolve().parent))
from tools.iupac_utils import (
    expand_iupac,
    calculate_gc,
    calculate_tm,
    evaluate_pair_heterodimer,
    evaluate_three_prime_dimer,
    three_prime_anchor_length,
    compatibility_score,
    check_dye_separation,
    screening_methodology_note,
    probe_synthesis_spec,
    MAX_DELTA_TM_C,
    MIN_CROSS_DIMER_DG,
    VETO_CROSS_DIMER_DG,
    MAX_THREE_PRIME_ANCHOR_BP,
    THREE_PRIME_WINDOW_NT,
)

from tools.candidate_selection import (
    select, selection_note, annotate, prefilter, cap_to_budget)

# MAX_VARIANTS is set from CLI --max-variants (default 16, matching config.yaml).
# Injected into module-level functions at startup by main().
_MAX_VARIANTS = 16


def load_coverage(path):
    """Map candidate rank -> RAA detection coverage from an M10 summary TSV.

    M10 computes this and, until now, nothing read it: the number that decides
    whether an assay detects circulating strains never reached the step that
    picks what to order. Returns {} when the file is absent, and the coverage
    gate then reports "no data" rather than quietly passing.
    """
    import csv
    if not path or not os.path.exists(path):
        return {}
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                out[int(row["rank"])] = float(row["raa_detected_pct"])
            except (KeyError, TypeError, ValueError):
                continue
    return out


def attach_coverage(cands, coverage):
    for c in cands:
        c["coverage"] = coverage.get(c.get("rank"))
    return cands

def run_multiplex_screening(file_a, file_b, out_dir, num_candidates=10,
                            dye_a="FAM", dye_b="HEX",
                            label_a="Target A", label_b="Target B",
                            coverage_a=None, coverage_b=None,
                            combination_budget=50000):
    os.makedirs(out_dir, exist_ok=True)

    # Real spectral-separation check (not a hardcoded string): two targets read on
    # the same detection channel cannot be distinguished in a multiplex.
    dye_status, dye_ok = check_dye_separation(dye_a, dye_b)
    print(f"[M11] Dye separation {dye_a} vs {dye_b}: {dye_status}")

    print(f"[M11] Loading {label_a} candidates from: {file_a}")
    with open(file_a, "r") as f:
        cands_a = json.load(f)[:num_candidates]

    print(f"[M11] Loading {label_b} candidates from: {file_b}")
    with open(file_b, "r") as f:
        cands_b = json.load(f)[:num_candidates]

    cov_a, cov_b = load_coverage(coverage_a), load_coverage(coverage_b)
    print(f"[M11] Coverage data: {label_a} {len(cov_a)} candidates, "
          f"{label_b} {len(cov_b)} candidates"
          + ("" if (cov_a and cov_b) else "  [WARNING: missing — the coverage gate "
                                          "cannot be applied and will reject]"))
    attach_coverage(cands_a, cov_a)
    attach_coverage(cands_b, cov_b)
    by_rank_a = {c.get("rank"): c for c in cands_a}
    by_rank_b = {c.get("rank"): c for c in cands_b}

    # Per-candidate gates BEFORE the pair loop. Specificity, human off-target
    # amplicons and coverage do not depend on what a candidate is paired with, so
    # applying them afterwards meant computing thermodynamics for pairs that could
    # never be recommended — and it was the sole reason a cut by Primer3 penalty
    # existed, which is what discarded the highest-coverage candidate of the run.
    kept_a, lvl_a, rej_a = prefilter(cands_a, label_a)
    kept_b, lvl_b, rej_b = prefilter(cands_b, label_b)
    for who, kept, before, lvl in ((label_a, kept_a, cands_a, lvl_a),
                                   (label_b, kept_b, cands_b, lvl_b)):
        print(f"[M11] {who}: {len(kept)}/{len(before)} candidates passed the "
              f"per-candidate gates"
              + (f" at {lvl['name']}" if lvl else " (none passed at any level)"))
        if lvl and lvl.get("level_counts"):
            print(f"[M11]   what each level would admit: "
                  + ", ".join(f"{k}={v}" for k, v in lvl["level_counts"].items()))
    (kept_a, kept_b), trimmed = cap_to_budget([kept_a, kept_b], combination_budget)
    if trimmed:
        print(f"[M11] Trimmed to {len(kept_a)} x {len(kept_b)} to stay within the "
              f"{combination_budget:,}-combination budget, keeping the best by "
              f"coverage then specificity")
    cands_a, cands_b = kept_a, kept_b
    prefilter_rejections = {**rej_a, **rej_b}
        
    print(f"[M11] Screening {len(cands_a)} x {len(cands_b)} = {len(cands_a) * len(cands_b)} candidate pairs...")
    
    results = []
    
    for cand_a in cands_a:
        rank_a = cand_a.get("rank", 1)
        fwd_a_deg = cand_a.get("fwd_degenerate_seq") or cand_a.get("fwd_seq")
        rev_a_deg = cand_a.get("rev_degenerate_seq") or cand_a.get("rev_seq")
        prb_a_deg = cand_a.get("probe_degenerate_seq") or cand_a.get("probe_seq")
        
        fwd_a_ref = cand_a.get("fwd_seq")
        rev_a_ref = cand_a.get("rev_seq")
        prb_a_ref = cand_a.get("probe_seq")
        
        tm_fwd_a = calculate_tm(fwd_a_deg)
        tm_rev_a = calculate_tm(rev_a_deg)
        tm_prb_a = calculate_tm(prb_a_deg)
        
        for cand_b in cands_b:
            rank_b = cand_b.get("rank", 1)
            fwd_b_deg = cand_b.get("fwd_degenerate_seq") or cand_b.get("fwd_seq")
            rev_b_deg = cand_b.get("rev_degenerate_seq") or cand_b.get("rev_seq")
            prb_b_deg = cand_b.get("probe_degenerate_seq") or cand_b.get("probe_seq")
            
            fwd_b_ref = cand_b.get("fwd_seq")
            rev_b_ref = cand_b.get("rev_seq")
            prb_b_ref = cand_b.get("probe_seq")
            
            tm_fwd_b = calculate_tm(fwd_b_deg)
            tm_rev_b = calculate_tm(rev_b_deg)
            tm_prb_b = calculate_tm(prb_b_deg)
            
            all_primer_tms = [tm_fwd_a, tm_rev_a, tm_fwd_b, tm_rev_b]
            delta_tm = max(all_primer_tms) - min(all_primer_tms)
            
            oligos_a = [("Fwd_A", fwd_a_deg), ("Rev_A", rev_a_deg), ("Probe_A", prb_a_deg)]
            oligos_b = [("Fwd_B", fwd_b_deg), ("Rev_B", rev_b_deg), ("Probe_B", prb_b_deg)]
            
            # --- Temp 1: RAA Mode (39.0°C) ---
            cross_dimers_raa = {}
            min_dg_raa = 0.0
            for name_a, seq_a in oligos_a:
                for name_b, seq_b in oligos_b:
                    if seq_a and seq_b:
                        dg, _ = evaluate_pair_heterodimer(seq_a, seq_b, temp_c=39.0)
                        cross_dimers_raa[f"{name_a}_vs_{name_b}"] = round(dg, 2)
                        if dg < min_dg_raa:
                            min_dg_raa = dg
                            
            # --- Temp 2: qPCR Mode (55.0°C) ---
            cross_dimers_qpcr = {}
            min_dg_qpcr = 0.0
            for name_a, seq_a in oligos_a:
                for name_b, seq_b in oligos_b:
                    if seq_a and seq_b:
                        dg, _ = evaluate_pair_heterodimer(seq_a, seq_b, temp_c=55.0)
                        cross_dimers_qpcr[f"{name_a}_vs_{name_b}"] = round(dg, 2)
                        if dg < min_dg_qpcr:
                            min_dg_qpcr = dg
                            
            # --- 3'-end complementarity (primer noise precursor) ---
            # A global ΔG can look acceptable while it is the 3' ends that
            # interact; only an annealed 3' terminus gets extended, and those
            # extension products are what amplify as primer noise.
            min_dg_3p = 0.0
            max_anchor_3p = 0
            three_prime_pairs = {}
            for name_a, seq_a in oligos_a:
                for name_b, seq_b in oligos_b:
                    if seq_a and seq_b:
                        dg3 = evaluate_three_prime_dimer(seq_a, seq_b, temp_c=39.0)
                        anchor = three_prime_anchor_length(seq_a, seq_b)
                        three_prime_pairs[f"{name_a}_vs_{name_b}"] = {
                            "anchor_bp": anchor, "dg_39c": round(dg3, 2)}
                        min_dg_3p = min(min_dg_3p, dg3)
                        max_anchor_3p = max(max_anchor_3p, anchor)

            status_raa = ("PASS" if (min_dg_raa >= MIN_CROSS_DIMER_DG
                                     and delta_tm <= MAX_DELTA_TM_C
                                     and max_anchor_3p <= MAX_THREE_PRIME_ANCHOR_BP)
                          else "WARNING")
            if min_dg_raa < VETO_CROSS_DIMER_DG:
                status_raa = "VETO (High Dimer Risk)"
            # A dye-channel collision means the two targets are optically
            # inseparable — no thermodynamic quality can rescue the pair.
            if not dye_ok:
                status_raa = "VETO (Dye Channel Collision)"

            score = compatibility_score(min_dg_raa, delta_tm, dye_ok)

            results.append({
                "pair_id": f"FluA_Rank{rank_a}_x_FluB_Rank{rank_b}",
                "rank_a": rank_a,
                "rank_b": rank_b,
                "score": score,
                "status_raa_39c": status_raa,
                "min_dg_raa_39c": round(min_dg_raa, 2),
                "min_dg_qpcr_55c": round(min_dg_qpcr, 2),
                "max_3prime_anchor_bp": max_anchor_3p,
                "min_dg_3prime_39c": round(min_dg_3p, 2),
                "three_prime_window_nt": THREE_PRIME_WINDOW_NT,
                "delta_tm": round(delta_tm, 2),
                "dye_a": dye_a,
                "dye_b": dye_b,
                "dye_ok": dye_ok,
                "dye_status": dye_status,
                "oligos_synthesis": {
                    "flu_a": {
                        "fwd": {"deg_seq": fwd_a_deg, "ref_seq": fwd_a_ref, "len": len(fwd_a_deg), "tm": round(tm_fwd_a, 1), "gc": calculate_gc(fwd_a_deg), "pos": f"{cand_a.get('fwd_start')}-{cand_a.get('fwd_end')}"},
                        "rev": {"deg_seq": rev_a_deg, "ref_seq": rev_a_ref, "len": len(rev_a_deg), "tm": round(tm_rev_a, 1), "gc": calculate_gc(rev_a_deg), "pos": f"{cand_a.get('rev_start')}-{cand_a.get('rev_end')}"},
                        "probe": {"deg_seq": prb_a_deg, "ref_seq": prb_a_ref, "len": len(prb_a_deg), "tm": round(tm_prb_a, 1), "gc": calculate_gc(prb_a_deg), "pos": f"{cand_a.get('probe_start')}-{cand_a.get('probe_end')}", "spec": probe_synthesis_spec(cand_a, dye_a)}
                    },
                    "flu_b": {
                        "fwd": {"deg_seq": fwd_b_deg, "ref_seq": fwd_b_ref, "len": len(fwd_b_deg), "tm": round(tm_fwd_b, 1), "gc": calculate_gc(fwd_b_deg), "pos": f"{cand_b.get('fwd_start')}-{cand_b.get('fwd_end')}"},
                        "rev": {"deg_seq": rev_b_deg, "ref_seq": rev_b_ref, "len": len(rev_b_deg), "tm": round(tm_rev_b, 1), "gc": calculate_gc(rev_b_deg), "pos": f"{cand_b.get('rev_start')}-{cand_b.get('rev_end')}"},
                        "probe": {"deg_seq": prb_b_deg, "ref_seq": prb_b_ref, "len": len(prb_b_deg), "tm": round(tm_prb_b, 1), "gc": calculate_gc(prb_b_deg), "pos": f"{cand_b.get('probe_start')}-{cand_b.get('probe_end')}", "spec": probe_synthesis_spec(cand_b, dye_b)}
                    }
                },
                "cross_dimers_39c": cross_dimers_raa,
                "cross_dimers_55c": cross_dimers_qpcr,
                "three_prime_dimers_39c": three_prime_pairs
            })
            
    # Gates before ranking. Sorting on the compatibility score alone previously
    # recommended a pair whose members were both human-specificity grade C while
    # grade A candidates with equal coverage went unused.
    selection = select(
        results,
        lambda r: [(label_a, by_rank_a.get(r["rank_a"], {})),
                   (label_b, by_rank_b.get(r["rank_b"], {}))],
        status_key="status_raa_39c",
    )
    # Candidates dropped before pairing never became combinations, so their
    # reasons would otherwise vanish from the report entirely.
    for reason, n in prefilter_rejections.items():
        selection["rejections"].setdefault(f"(pre-pairing) {reason}", n)
    if selection["selected"]:
        chosen = selection["selected"]
        dropped = [r for r in results if r not in chosen]
        dropped.sort(key=lambda x: x["score"], reverse=True)
        results = chosen + dropped
        print(f"[M11] Selection level {selection['level']} "
              f"({selection['level_label']}): {selection['n_selected']} of "
              f"{selection['n_input']} combinations passed the gates")
    else:
        results.sort(key=lambda x: x["score"], reverse=True)
        print("[M11] WARNING: no combination passed the non-negotiable gates. "
              "Nothing is recommended; see the report for reasons.")
    
    # Stamp the selection verdict onto each combination so it survives into the
    # JSON: the supplementary export reads that file and must not label whatever
    # sits at index 0 as "BEST RECOMMENDED" without knowing it passed the gates.
    annotate(results, selection)

    json_path = os.path.join(out_dir, "multiplex_candidates.json")
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"[M11] Saved multiplex results JSON to: {json_path}")
    
    # Write Markdown Report
    md_path = os.path.join(out_dir, "multiplex_report.md")
    with open(md_path, "w") as f:
        f.write("# Dual-Target Multiplex Assay Compatibility & Synthesis Report\n")
        f.write(f"## Target Pair: {label_a} ({dye_a}) + {label_b} ({dye_b})\n\n")
        f.write(f"**Spectral dye separation:** {dye_status}\n\n")
        f.write(f"Screened **{len(results)} candidate combinations** ($10 \\times 10$) across two temperature regimes ($39.0^\\circ\\text{{C}}$ RT-RAA & $55.0^\\circ\\text{{C}}$ qPCR).\n\n")
        f.write("### Top 5 Recommended Multiplex Combinations\n\n")
        f.write("| Rank | Candidate Pair | Compatibility Score | RAA ΔG (39°C) | 3' Anchor | qPCR ΔG (55°C) | ΔTm (°C) | Status (39°C) |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        for idx, r in enumerate(results[:5], 1):
            f.write(f"| {idx} | **{r['pair_id']}** | **{r['score']}** | {r['min_dg_raa_39c']} kcal/mol | {r['max_3prime_anchor_bp']} bp | {r['min_dg_qpcr_55c']} kcal/mol | {r['delta_tm']}°C | `{r['status_raa_39c']}` |\n")

        f.write(f"\n\n### Rank 1 Oligo Synthesis Specifications (`{results[0]['pair_id']}`)\n\n")
        r1 = results[0]
        oligos = r1['oligos_synthesis']
        
        f.write(f"#### 🟩 Target A: {label_a} ({dye_a} Channel)\n\n")
        f.write("| Oligo Name | Type | Degenerate Sequence (5' -> 3') | Length | Tm (°C) | GC (%) | Modification / Label |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        f.write(f"| FluA_Fwd | Fwd Primer | `{oligos['flu_a']['fwd']['deg_seq']}` | {oligos['flu_a']['fwd']['len']} nt | {oligos['flu_a']['fwd']['tm']}°C | {oligos['flu_a']['fwd']['gc']}% | Unlabeled |\n")
        f.write(f"| FluA_Rev | Rev Primer | `{oligos['flu_a']['rev']['deg_seq']}` | {oligos['flu_a']['rev']['len']} nt | {oligos['flu_a']['rev']['tm']}°C | {oligos['flu_a']['rev']['gc']}% | Unlabeled |\n")
        _sp = oligos['flu_a']['probe']['spec']
        f.write(f"| FluA_Probe | {_sp['probe_type']}-Probe | `{_sp['order_sequence']}` | {oligos['flu_a']['probe']['len']} nt | {oligos['flu_a']['probe']['tm']}°C | {oligos['flu_a']['probe']['gc']}% | {_sp['modification']} |\n\n")
        if _sp['thf_index'] is not None:
            f.write(f"Abasic site at position {_sp['thf_index'] + 1} of the FluA probe (1-based).\n\n")
        if _sp['notes']:
            f.write(f"!!! note\n    FluA probe: {_sp['notes']}.\n\n")

        f.write(f"#### 🟨 Target B: {label_b} ({dye_b} Channel)\n\n")
        f.write("| Oligo Name | Type | Degenerate Sequence (5' -> 3') | Length | Tm (°C) | GC (%) | Modification / Label |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        f.write(f"| FluB_Fwd | Fwd Primer | `{oligos['flu_b']['fwd']['deg_seq']}` | {oligos['flu_b']['fwd']['len']} nt | {oligos['flu_b']['fwd']['tm']}°C | {oligos['flu_b']['fwd']['gc']}% | Unlabeled |\n")
        f.write(f"| FluB_Rev | Rev Primer | `{oligos['flu_b']['rev']['deg_seq']}` | {oligos['flu_b']['rev']['len']} nt | {oligos['flu_b']['rev']['tm']}°C | {oligos['flu_b']['rev']['gc']}% | Unlabeled |\n")
        _sp = oligos['flu_b']['probe']['spec']
        f.write(f"| FluB_Probe | {_sp['probe_type']}-Probe | `{_sp['order_sequence']}` | {oligos['flu_b']['probe']['len']} nt | {oligos['flu_b']['probe']['tm']}°C | {oligos['flu_b']['probe']['gc']}% | {_sp['modification']} |\n\n")
        if _sp['thf_index'] is not None:
            f.write(f"Abasic site at position {_sp['thf_index'] + 1} of the FluB probe (1-based).\n\n")
        if _sp['notes']:
            f.write(f"!!! note\n    FluB probe: {_sp['notes']}.\n\n")

        f.write("### Cross-Dimer Inter-Assay Interactions\n\n")
        f.write("| Oligo Pair | Cross-Dimer ΔG (39°C) | 3' Anchor | 3'-End ΔG (39°C) | Cross-Dimer ΔG (55°C) |\n")
        f.write("|---|---|---|---|---|\n")
        for pair_name, dg in r1['cross_dimers_39c'].items():
            dg_55 = r1['cross_dimers_55c'].get(pair_name, "-")
            tp = r1.get('three_prime_dimers_39c', {}).get(pair_name) or {}
            f.write(f"| {pair_name} | {dg} kcal/mol | {tp.get('anchor_bp', '-')} bp "
                    f"| {tp.get('dg_39c', '-')} kcal/mol | {dg_55} kcal/mol |\n")

        f.write(selection_note(selection))
        f.write(screening_methodology_note())

    print(f"[M11] Saved Markdown report to: {md_path}")

def main():
    parser = argparse.ArgumentParser(
        description="Dual-Target Multiplex Compatibility Checker (Influenza A + B MP)"
    )
    parser.add_argument(
        "--file-a",
        default="results/influenza-a-mp-indo-only/m09/degenerate/segment-MP_region_1/final_candidates.json",
        help="Path to Target A final_candidates.json (Influenza A)",
    )
    parser.add_argument(
        "--file-b",
        default="results/influenza-b-mp-indo-only/m09/degenerate/segment-MP_region_1/final_candidates.json",
        help="Path to Target B final_candidates.json (Influenza B)",
    )
    parser.add_argument("--outdir", default="results/multiplex-indo-only", help="Output directory")
    parser.add_argument("--num-candidates", type=int, default=10, help="Number of top candidates to screen per target")
    parser.add_argument("--dye-a", default="FAM", help="Detection fluorophore for Target A (e.g. FAM)")
    parser.add_argument("--dye-b", default="HEX", help="Detection fluorophore for Target B (e.g. HEX)")
    parser.add_argument("--label-a", default="Target A", help="Human-readable name for Target A")
    parser.add_argument("--label-b", default="Target B", help="Human-readable name for Target B")
    parser.add_argument("--combination-budget", type=int, default=50000,
                        help="Maximum combinations to evaluate. Survivors of the "
                             "per-candidate gates are trimmed by coverage then "
                             "specificity if the product exceeds this. Pairwise "
                             "rarely reaches it; the cubic triplex case does.")
    parser.add_argument("--coverage-a", default=None,
                        help="Path to Target A M10 coverage_summary.tsv. Without it the "
                             "coverage gate cannot be evaluated and every candidate is "
                             "rejected for missing data rather than silently passed.")
    parser.add_argument("--coverage-b", default=None,
                        help="Path to Target B M10 coverage_summary.tsv")
    parser.add_argument(
        "--max-variants",
        type=int,
        default=16,
        help="IUPAC expansion cap per oligo (default: 16, matches config.yaml degenerate.max_iupac_variants)",
    )
    args = parser.parse_args()

    # Propagate --max-variants into the shared utility module caches
    import tools.iupac_utils as _iu
    _iu.DEFAULT_MAX_VARIANTS = args.max_variants

    run_multiplex_screening(args.file_a, args.file_b, args.outdir, args.num_candidates,
                            dye_a=args.dye_a, dye_b=args.dye_b,
                            label_a=args.label_a, label_b=args.label_b,
                            coverage_a=args.coverage_a, coverage_b=args.coverage_b,
                            combination_budget=args.combination_budget)

if __name__ == "__main__":
    main()
