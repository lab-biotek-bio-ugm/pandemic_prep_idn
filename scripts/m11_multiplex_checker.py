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
    check_dye_separation,
    MAX_DELTA_TM_C,
    MIN_CROSS_DIMER_DG,
    VETO_CROSS_DIMER_DG,
)

# MAX_VARIANTS is set from CLI --max-variants (default 16, matching config.yaml).
# Injected into module-level functions at startup by main().
_MAX_VARIANTS = 16

def run_multiplex_screening(file_a, file_b, out_dir, num_candidates=10,
                            dye_a="FAM", dye_b="HEX",
                            label_a="Target A", label_b="Target B"):
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
                            
            status_raa = ("PASS" if (min_dg_raa >= MIN_CROSS_DIMER_DG and delta_tm <= MAX_DELTA_TM_C)
                          else "WARNING")
            if min_dg_raa < VETO_CROSS_DIMER_DG:
                status_raa = "VETO (High Dimer Risk)"
            # A dye-channel collision means the two targets are optically
            # inseparable — no thermodynamic quality can rescue the pair.
            if not dye_ok:
                status_raa = "VETO (Dye Channel Collision)"

            score = 100.0
            if min_dg_raa < 0:
                score += min_dg_raa * 3.5
            score -= delta_tm * 4.0
            if not dye_ok:
                score = 0.0
            score = max(0.0, round(score, 1))
            
            results.append({
                "pair_id": f"FluA_Rank{rank_a}_x_FluB_Rank{rank_b}",
                "rank_a": rank_a,
                "rank_b": rank_b,
                "score": score,
                "status_raa_39c": status_raa,
                "min_dg_raa_39c": round(min_dg_raa, 2),
                "min_dg_qpcr_55c": round(min_dg_qpcr, 2),
                "delta_tm": round(delta_tm, 2),
                "dye_a": dye_a,
                "dye_b": dye_b,
                "dye_ok": dye_ok,
                "dye_status": dye_status,
                "oligos_synthesis": {
                    "flu_a": {
                        "fwd": {"deg_seq": fwd_a_deg, "ref_seq": fwd_a_ref, "len": len(fwd_a_deg), "tm": round(tm_fwd_a, 1), "gc": calculate_gc(fwd_a_deg), "pos": f"{cand_a.get('fwd_start')}-{cand_a.get('fwd_end')}"},
                        "rev": {"deg_seq": rev_a_deg, "ref_seq": rev_a_ref, "len": len(rev_a_deg), "tm": round(tm_rev_a, 1), "gc": calculate_gc(rev_a_deg), "pos": f"{cand_a.get('rev_start')}-{cand_a.get('rev_end')}"},
                        "probe": {"deg_seq": prb_a_deg, "ref_seq": prb_a_ref, "len": len(prb_a_deg), "tm": round(tm_prb_a, 1), "gc": calculate_gc(prb_a_deg), "pos": f"{cand_a.get('probe_start')}-{cand_a.get('probe_end')}", "label": f"5'-{dye_a} / 3'-BHQ1 / Internal-THF"}
                    },
                    "flu_b": {
                        "fwd": {"deg_seq": fwd_b_deg, "ref_seq": fwd_b_ref, "len": len(fwd_b_deg), "tm": round(tm_fwd_b, 1), "gc": calculate_gc(fwd_b_deg), "pos": f"{cand_b.get('fwd_start')}-{cand_b.get('fwd_end')}"},
                        "rev": {"deg_seq": rev_b_deg, "ref_seq": rev_b_ref, "len": len(rev_b_deg), "tm": round(tm_rev_b, 1), "gc": calculate_gc(rev_b_deg), "pos": f"{cand_b.get('rev_start')}-{cand_b.get('rev_end')}"},
                        "probe": {"deg_seq": prb_b_deg, "ref_seq": prb_b_ref, "len": len(prb_b_deg), "tm": round(tm_prb_b, 1), "gc": calculate_gc(prb_b_deg), "pos": f"{cand_b.get('probe_start')}-{cand_b.get('probe_end')}", "label": f"5'-{dye_b} / 3'-BHQ1 / Internal-THF"}
                    }
                },
                "cross_dimers_39c": cross_dimers_raa,
                "cross_dimers_55c": cross_dimers_qpcr
            })
            
    results.sort(key=lambda x: x["score"], reverse=True)
    
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
        f.write("| Rank | Candidate Pair | Compatibility Score | RAA ΔG (39°C) | qPCR ΔG (55°C) | ΔTm (°C) | Status (39°C) |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for idx, r in enumerate(results[:5], 1):
            f.write(f"| {idx} | **{r['pair_id']}** | **{r['score']}%** | {r['min_dg_raa_39c']} kcal/mol | {r['min_dg_qpcr_55c']} kcal/mol | {r['delta_tm']}°C | `{r['status_raa_39c']}` |\n")
            
        f.write("\n\n### Rank 1 Oligo Synthesis Specifications (`FluA_Rank7_x_FluB_Rank6`)\n\n")
        r1 = results[0]
        oligos = r1['oligos_synthesis']
        
        f.write(f"#### 🟩 Target A: {label_a} ({dye_a} Channel)\n\n")
        f.write("| Oligo Name | Type | Degenerate Sequence (5' -> 3') | Length | Tm (°C) | GC (%) | Modification / Label |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        f.write(f"| FluA_Fwd | Fwd Primer | `{oligos['flu_a']['fwd']['deg_seq']}` | {oligos['flu_a']['fwd']['len']} nt | {oligos['flu_a']['fwd']['tm']}°C | {oligos['flu_a']['fwd']['gc']}% | Unlabeled |\n")
        f.write(f"| FluA_Rev | Rev Primer | `{oligos['flu_a']['rev']['deg_seq']}` | {oligos['flu_a']['rev']['len']} nt | {oligos['flu_a']['rev']['tm']}°C | {oligos['flu_a']['rev']['gc']}% | Unlabeled |\n")
        f.write(f"| FluA_Probe | Exo-Probe | `{oligos['flu_a']['probe']['deg_seq']}` | {oligos['flu_a']['probe']['len']} nt | {oligos['flu_a']['probe']['tm']}°C | {oligos['flu_a']['probe']['gc']}% | {oligos['flu_a']['probe']['label']} |\n\n")

        f.write(f"#### 🟨 Target B: {label_b} ({dye_b} Channel)\n\n")
        f.write("| Oligo Name | Type | Degenerate Sequence (5' -> 3') | Length | Tm (°C) | GC (%) | Modification / Label |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        f.write(f"| FluB_Fwd | Fwd Primer | `{oligos['flu_b']['fwd']['deg_seq']}` | {oligos['flu_b']['fwd']['len']} nt | {oligos['flu_b']['fwd']['tm']}°C | {oligos['flu_b']['fwd']['gc']}% | Unlabeled |\n")
        f.write(f"| FluB_Rev | Rev Primer | `{oligos['flu_b']['rev']['deg_seq']}` | {oligos['flu_b']['rev']['len']} nt | {oligos['flu_b']['rev']['tm']}°C | {oligos['flu_b']['rev']['gc']}% | Unlabeled |\n")
        f.write(f"| FluB_Probe | Exo-Probe | `{oligos['flu_b']['probe']['deg_seq']}` | {oligos['flu_b']['probe']['len']} nt | {oligos['flu_b']['probe']['tm']}°C | {oligos['flu_b']['probe']['gc']}% | {oligos['flu_b']['probe']['label']} |\n\n")

        f.write("### Cross-Dimer Inter-Assay Interactions\n\n")
        f.write("| Oligo Pair | Cross-Dimer ΔG (39°C) | Cross-Dimer ΔG (55°C) |\n")
        f.write("|---|---|---|\n")
        for pair_name, dg in r1['cross_dimers_39c'].items():
            dg_55 = r1['cross_dimers_55c'].get(pair_name, "-")
            f.write(f"| {pair_name} | {dg} kcal/mol | {dg_55} kcal/mol |\n")
            
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
                            label_a=args.label_a, label_b=args.label_b)

if __name__ == "__main__":
    main()
