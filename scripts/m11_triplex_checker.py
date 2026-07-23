#!/usr/bin/env python3
"""
Module 11: Triplex Multiplex Assay Compatibility Screening Engine
Evaluates three-way cross-dimer thermodynamics, spectral dye separation,
and Tm alignment for triplex assays (Influenza A + Influenza B + Measles).
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

# MAX_VARIANTS is set from CLI --max-variants (default 16, matching config.yaml).
_MAX_VARIANTS = 16


def _safe_tm(seq):
    """Tm of an oligo, or 0.0 if the oligo is absent (e.g. a candidate with no probe)."""
    return calculate_tm(seq) if seq else 0.0


def run_triplex_screening(file_a, file_b, file_c, out_dir, num_candidates=10,
                          dye_a="FAM", dye_b="HEX", dye_c="Cy5",
                          label_a="Target A", label_b="Target B", label_c="Target C"):
    os.makedirs(out_dir, exist_ok=True)

    # Real 3-way spectral-separation check (not hardcoded): all three targets must
    # occupy distinct detection channels or they cannot be told apart.
    dye_status, dye_ok = check_dye_separation(dye_a, dye_b, dye_c)
    print(f"[Triplex] Dye separation {dye_a}/{dye_b}/{dye_c}: {dye_status}")

    print(f"[Triplex] Loading {label_a} candidates from: {file_a}")
    with open(file_a, "r") as f:
        cands_a = json.load(f)[:num_candidates]

    print(f"[Triplex] Loading {label_b} candidates from: {file_b}")
    with open(file_b, "r") as f:
        cands_b = json.load(f)[:num_candidates]

    print(f"[Triplex] Loading {label_c} candidates from: {file_c}")
    with open(file_c, "r") as f:
        cands_c = json.load(f)[:num_candidates]
        
    total_combs = len(cands_a) * len(cands_b) * len(cands_c)
    print(f"[Triplex] Screening {len(cands_a)} x {len(cands_b)} x {len(cands_c)} = {total_combs} triplex combinations...")
    
    results = []
    
    for cand_a in cands_a:
        rank_a = cand_a.get("rank", 1)
        fwd_a = cand_a.get("fwd_degenerate_seq") or cand_a.get("fwd_seq")
        rev_a = cand_a.get("rev_degenerate_seq") or cand_a.get("rev_seq")
        prb_a = cand_a.get("probe_degenerate_seq") or cand_a.get("probe_seq")
        
        tm_fwd_a = _safe_tm(fwd_a)
        tm_rev_a = _safe_tm(rev_a)
        tm_prb_a = _safe_tm(prb_a)
        
        for cand_b in cands_b:
            rank_b = cand_b.get("rank", 1)
            fwd_b = cand_b.get("fwd_degenerate_seq") or cand_b.get("fwd_seq")
            rev_b = cand_b.get("rev_degenerate_seq") or cand_b.get("rev_seq")
            prb_b = cand_b.get("probe_degenerate_seq") or cand_b.get("probe_seq")
            
            tm_fwd_b = _safe_tm(fwd_b)
            tm_rev_b = _safe_tm(rev_b)
            tm_prb_b = _safe_tm(prb_b)
            
            for cand_c in cands_c:
                rank_c = cand_c.get("rank", 1)
                fwd_c = cand_c.get("fwd_degenerate_seq") or cand_c.get("fwd_seq")
                rev_c = cand_c.get("rev_degenerate_seq") or cand_c.get("rev_seq")
                prb_c = cand_c.get("probe_degenerate_seq") or cand_c.get("probe_seq")
                
                tm_fwd_c = _safe_tm(fwd_c)
                tm_rev_c = _safe_tm(rev_c)
                tm_prb_c = _safe_tm(prb_c)
                
                all_primer_tms = [tm_fwd_a, tm_rev_a, tm_fwd_b, tm_rev_b, tm_fwd_c, tm_rev_c]
                delta_tm = max(all_primer_tms) - min(all_primer_tms)
                
                oligos_a = [("Fwd_A", fwd_a), ("Rev_A", rev_a), ("Probe_A", prb_a)]
                oligos_b = [("Fwd_B", fwd_b), ("Rev_B", rev_b), ("Probe_B", prb_b)]
                oligos_c = [("Fwd_C", fwd_c), ("Rev_C", rev_c), ("Probe_C", prb_c)]
                
                # Check cross-dimers for all three pairs (A-B, A-C, B-C)
                cross_dimers_39c = {}
                cross_dimers_55c = {}
                min_dg_39c = 0.0
                min_dg_55c = 0.0
                
                # Pair A-B
                for name_a, seq_a in oligos_a:
                    for name_b, seq_b in oligos_b:
                        if not (seq_a and seq_b):
                            continue  # skip absent oligos (e.g. probe-less candidate)
                        dg_39, _ = evaluate_pair_heterodimer(seq_a, seq_b, temp_c=39.0)
                        dg_55, _ = evaluate_pair_heterodimer(seq_a, seq_b, temp_c=55.0)
                        cross_dimers_39c[f"{name_a}_vs_{name_b}"] = round(dg_39, 2)
                        cross_dimers_55c[f"{name_a}_vs_{name_b}"] = round(dg_55, 2)
                        if dg_39 < min_dg_39c:
                            min_dg_39c = dg_39
                        if dg_55 < min_dg_55c:
                            min_dg_55c = dg_55
                
                # Pair A-C
                for name_a, seq_a in oligos_a:
                    for name_c, seq_c in oligos_c:
                        if not (seq_a and seq_c):
                            continue
                        dg_39, _ = evaluate_pair_heterodimer(seq_a, seq_c, temp_c=39.0)
                        dg_55, _ = evaluate_pair_heterodimer(seq_a, seq_c, temp_c=55.0)
                        cross_dimers_39c[f"{name_a}_vs_{name_c}"] = round(dg_39, 2)
                        cross_dimers_55c[f"{name_a}_vs_{name_c}"] = round(dg_55, 2)
                        if dg_39 < min_dg_39c:
                            min_dg_39c = dg_39
                        if dg_55 < min_dg_55c:
                            min_dg_55c = dg_55
                            
                # Pair B-C
                for name_b, seq_b in oligos_b:
                    for name_c, seq_c in oligos_c:
                        if not (seq_b and seq_c):
                            continue
                        dg_39, _ = evaluate_pair_heterodimer(seq_b, seq_c, temp_c=39.0)
                        dg_55, _ = evaluate_pair_heterodimer(seq_b, seq_c, temp_c=55.0)
                        cross_dimers_39c[f"{name_b}_vs_{name_c}"] = round(dg_39, 2)
                        cross_dimers_55c[f"{name_b}_vs_{name_c}"] = round(dg_55, 2)
                        if dg_39 < min_dg_39c:
                            min_dg_39c = dg_39
                        if dg_55 < min_dg_55c:
                            min_dg_55c = dg_55
                            
                # 3'-end complementarity across all three target pairs. Only an
                # annealed 3' terminus gets extended, and those products are what
                # amplify as primer noise (TwistAmp manual 2.6).
                min_dg_3p = 0.0
                max_anchor_3p = 0
                three_prime_dimers = {}
                for (grp_x, grp_y) in ((oligos_a, oligos_b), (oligos_a, oligos_c), (oligos_b, oligos_c)):
                    for name_x, seq_x in grp_x:
                        for name_y, seq_y in grp_y:
                            if not (seq_x and seq_y):
                                continue
                            dg3 = evaluate_three_prime_dimer(seq_x, seq_y, temp_c=39.0)
                            anchor = three_prime_anchor_length(seq_x, seq_y)
                            three_prime_dimers[f"{name_x}_vs_{name_y}"] = {
                                "anchor_bp": anchor, "dg_39c": round(dg3, 2)}
                            min_dg_3p = min(min_dg_3p, dg3)
                            max_anchor_3p = max(max_anchor_3p, anchor)

                status_39c = ("PASS" if (min_dg_39c >= MIN_CROSS_DIMER_DG
                                         and delta_tm <= MAX_DELTA_TM_C
                                         and max_anchor_3p <= MAX_THREE_PRIME_ANCHOR_BP)
                              else "WARNING")
                if min_dg_39c < VETO_CROSS_DIMER_DG:
                    status_39c = "VETO"
                # A dye-channel collision makes the three targets optically
                # inseparable regardless of thermodynamic quality.
                if not dye_ok:
                    status_39c = "VETO (Dye Channel Collision)"

                score = compatibility_score(min_dg_39c, delta_tm, dye_ok)

                results.append({
                    "triplex_id": f"{label_a}_R{rank_a}_{label_b}_R{rank_b}_{label_c}_R{rank_c}",
                    "score": score,
                    "status_39c": status_39c,
                    "dye_ok": dye_ok,
                    "dye_status": dye_status,
                    "min_dg_39c": round(min_dg_39c, 2),
                    "min_dg_55c": round(min_dg_55c, 2),
                    "max_3prime_anchor_bp": max_anchor_3p,
                    "min_dg_3prime_39c": round(min_dg_3p, 2),
                    "three_prime_window_nt": THREE_PRIME_WINDOW_NT,
                    "delta_tm": round(delta_tm, 2),
                    "oligos_synthesis": {
                        "flu_a": {
                            "fwd": fwd_a, "rev": rev_a, "probe": prb_a,
                            "fwd_tm": round(tm_fwd_a, 1), "rev_tm": round(tm_rev_a, 1), "probe_tm": round(tm_prb_a, 1),
                            "probe_spec": probe_synthesis_spec(cand_a, dye_a, quencher="BHQ1")
                        },
                        "flu_b": {
                            "fwd": fwd_b, "rev": rev_b, "probe": prb_b,
                            "fwd_tm": round(tm_fwd_b, 1), "rev_tm": round(tm_rev_b, 1), "probe_tm": round(tm_prb_b, 1),
                            "probe_spec": probe_synthesis_spec(cand_b, dye_b, quencher="BHQ1")
                        },
                        "measles": {
                            "fwd": fwd_c, "rev": rev_c, "probe": prb_c,
                            "fwd_tm": round(tm_fwd_c, 1), "rev_tm": round(tm_rev_c, 1), "probe_tm": round(tm_prb_c, 1),
                            "probe_spec": probe_synthesis_spec(cand_c, dye_c, quencher="BHQ2")
                        }
                    },
                    "cross_dimers_39c": cross_dimers_39c,
                    "cross_dimers_55c": cross_dimers_55c,
                    "three_prime_dimers_39c": three_prime_dimers
                })
                
    results.sort(key=lambda x: x["score"], reverse=True)
    
    # Save JSON
    json_path = os.path.join(out_dir, "triplex_candidates.json")
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"[Triplex] Saved JSON candidates to: {json_path}")
    
    # Write Markdown Report
    md_path = os.path.join(out_dir, "triplex_report.md")
    with open(md_path, "w") as f:
        f.write("# Triplex Multiplex Assay Compatibility & Synthesis Report\n")
        f.write(f"## Targets: {label_a} ({dye_a}) + {label_b} ({dye_b}) + {label_c} ({dye_c})\n\n")
        f.write(f"**Spectral dye separation:** {dye_status}\n\n")
        f.write(f"Screened **{len(results)} triplex combinations** ({len(cands_a)} x {len(cands_b)} x {len(cands_c)}) at RAA operating temperature ($39.0^\\circ\\text{{C}}$) and qPCR benchmark ($55.0^\\circ\\text{{C}}$).\n\n")
        f.write("### Top 5 Recommended Triplex Combinations\n\n")
        f.write("| Rank | Triplex Candidate Set | Compatibility Score | Worst RAA ΔG (39°C) | Longest 3' Anchor | Worst qPCR ΔG (55°C) | Max ΔTm (°C) | Status (39°C) |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        for idx, r in enumerate(results[:5], 1):
            f.write(f"| {idx} | **{r['triplex_id']}** | **{r['score']}** | {r['min_dg_39c']} kcal/mol | {r['max_3prime_anchor_bp']} bp | {r['min_dg_55c']} kcal/mol | {r['delta_tm']}°C | `{r['status_39c']}` |\n")
            
        # Detailed specs of Rank 1
        r1 = results[0]
        f.write(f"\n\n### Rank 1 Triplex Specifications (`{r1['triplex_id']}`)\n\n")
        
        # Target A (Flu A)
        f.write(f"#### 🟩 Target A: {label_a} ({dye_a} Channel)\n")
        f.write("| Oligo Name | Type | Sequence (5' -> 3') | Tm (°C) | Modification / Label |\n")
        f.write("|---|---|---|---|---|\n")
        f.write(f"| FluA_Fwd | Fwd Primer | `{r1['oligos_synthesis']['flu_a']['fwd']}` | {r1['oligos_synthesis']['flu_a']['fwd_tm']}°C | Unlabeled |\n")
        f.write(f"| FluA_Rev | Rev Primer | `{r1['oligos_synthesis']['flu_a']['rev']}` | {r1['oligos_synthesis']['flu_a']['rev_tm']}°C | Unlabeled |\n")
        _sp = r1['oligos_synthesis']['flu_a']['probe_spec']
        f.write(f"| FluA_Probe | {_sp['probe_type']}-Probe | `{_sp['order_sequence']}` | {r1['oligos_synthesis']['flu_a']['probe_tm']}°C | {_sp['modification']} |\n\n")
        if _sp['thf_index'] is not None:
            f.write(f"Abasic site at position {_sp['thf_index'] + 1} of the FluA probe (1-based).\n\n")
        if _sp['notes']:
            f.write(f"!!! note\n    FluA probe: {_sp['notes']}.\n\n")
        
        # Target B (Flu B)
        f.write(f"#### 🟨 Target B: {label_b} ({dye_b} Channel)\n")
        f.write("| Oligo Name | Type | Sequence (5' -> 3') | Tm (°C) | Modification / Label |\n")
        f.write("|---|---|---|---|---|\n")
        f.write(f"| FluB_Fwd | Fwd Primer | `{r1['oligos_synthesis']['flu_b']['fwd']}` | {r1['oligos_synthesis']['flu_b']['fwd_tm']}°C | Unlabeled |\n")
        f.write(f"| FluB_Rev | Rev Primer | `{r1['oligos_synthesis']['flu_b']['rev']}` | {r1['oligos_synthesis']['flu_b']['rev_tm']}°C | Unlabeled |\n")
        _sp = r1['oligos_synthesis']['flu_b']['probe_spec']
        f.write(f"| FluB_Probe | {_sp['probe_type']}-Probe | `{_sp['order_sequence']}` | {r1['oligos_synthesis']['flu_b']['probe_tm']}°C | {_sp['modification']} |\n\n")
        if _sp['thf_index'] is not None:
            f.write(f"Abasic site at position {_sp['thf_index'] + 1} of the FluB probe (1-based).\n\n")
        if _sp['notes']:
            f.write(f"!!! note\n    FluB probe: {_sp['notes']}.\n\n")
        
        # Target C (Measles)
        f.write(f"#### 🟦 Target C: {label_c} ({dye_c} Channel)\n")
        f.write("| Oligo Name | Type | Sequence (5' -> 3') | Tm (°C) | Modification / Label |\n")
        f.write("|---|---|---|---|---|\n")
        f.write(f"| Measles_Fwd | Fwd Primer | `{r1['oligos_synthesis']['measles']['fwd']}` | {r1['oligos_synthesis']['measles']['fwd_tm']}°C | Unlabeled |\n")
        f.write(f"| Measles_Rev | Rev Primer | `{r1['oligos_synthesis']['measles']['rev']}` | {r1['oligos_synthesis']['measles']['rev_tm']}°C | Unlabeled |\n")
        _sp = r1['oligos_synthesis']['measles']['probe_spec']
        f.write(f"| Measles_Probe | {_sp['probe_type']}-Probe | `{_sp['order_sequence']}` | {r1['oligos_synthesis']['measles']['probe_tm']}°C | {_sp['modification']} |\n\n")
        if _sp['thf_index'] is not None:
            f.write(f"Abasic site at position {_sp['thf_index'] + 1} of the Measles probe (1-based).\n\n")
        if _sp['notes']:
            f.write(f"!!! note\n    Measles probe: {_sp['notes']}.\n\n")
        
        # Cross dimer detail. The listing threshold is deliberately more permissive
        # than the PASS threshold so a reader can see what sits just short of it;
        # at 39°C a clean random pair already sits near -4 kcal/mol (see the
        # methodology note below), so -5 is "worth a look", not "a problem".
        f.write("### Notable Oligo Interactions (ΔG < -5.0 kcal/mol at 39°C, or a 3' anchor over the limit)\n\n")
        f.write("| Oligo Pair | Cross-Dimer ΔG (39°C) | 3' Anchor | 3'-End ΔG (39°C) | Cross-Dimer ΔG (55°C) |\n")
        f.write("|---|---|---|---|---|\n")
        listed = 0
        for pair_name, dg in r1['cross_dimers_39c'].items():
            tp = r1.get('three_prime_dimers_39c', {}).get(pair_name) or {}
            if dg < -5.0 or tp.get('anchor_bp', 0) > MAX_THREE_PRIME_ANCHOR_BP:
                dg_55 = r1['cross_dimers_55c'].get(pair_name, "-")
                f.write(f"| {pair_name} | {dg} kcal/mol | {tp.get('anchor_bp', '-')} bp "
                        f"| {tp.get('dg_39c', '-')} kcal/mol | {dg_55} kcal/mol |\n")
                listed += 1
        if not listed:
            f.write("| _(none)_ | | | | |\n")

        f.write(screening_methodology_note())

    print(f"[Triplex] Saved Markdown report to: {md_path}")

def main():
    parser = argparse.ArgumentParser(
        description="Triplex Multiplex Compatibility Checker (Influenza A + Influenza B + Measles)"
    )
    parser.add_argument("--file-a", required=True, help="Path to Target A final_candidates.json (Influenza A)")
    parser.add_argument("--file-b", required=True, help="Path to Target B final_candidates.json (Influenza B)")
    parser.add_argument("--file-c", required=True, help="Path to Target C final_candidates.json (Measles)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--num-candidates", type=int, default=10, help="Number of top candidates to screen per target")
    parser.add_argument("--dye-a", default="FAM", help="Detection fluorophore for Target A")
    parser.add_argument("--dye-b", default="HEX", help="Detection fluorophore for Target B")
    parser.add_argument("--dye-c", default="Cy5", help="Detection fluorophore for Target C")
    parser.add_argument("--label-a", default="Target A", help="Human-readable name for Target A")
    parser.add_argument("--label-b", default="Target B", help="Human-readable name for Target B")
    parser.add_argument("--label-c", default="Target C", help="Human-readable name for Target C")
    parser.add_argument(
        "--max-variants",
        type=int,
        default=16,
        help="IUPAC expansion cap per oligo (default: 16, matches config.yaml degenerate.max_iupac_variants)",
    )
    args = parser.parse_args()

    # Propagate --max-variants into the shared utility module
    import tools.iupac_utils as _iu
    _iu.DEFAULT_MAX_VARIANTS = args.max_variants

    run_triplex_screening(args.file_a, args.file_b, args.file_c, args.outdir, args.num_candidates,
                          dye_a=args.dye_a, dye_b=args.dye_b, dye_c=args.dye_c,
                          label_a=args.label_a, label_b=args.label_b, label_c=args.label_c)

if __name__ == "__main__":
    main()
