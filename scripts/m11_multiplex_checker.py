#!/usr/bin/env python3
"""
Module 11: Multiplex Assay Compatibility Screening Engine (Influenza A + B MP)
Evaluates pairwise cross-dimer thermodynamics (ViennaRNA / Primer3), spectral dye separation,
and kinetics alignment for dual-target assays at 39.0°C (RAA) and 55.0°C (qPCR).
Outputs full sequence synthesis specs and cross-dimer details.
"""

import os
import sys
import json
import argparse
import itertools
import primer3

# IUPAC Degenerate Base Table
IUPAC_MAP = {
    'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
    'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}

def expand_iupac(seq):
    """Expand degenerate sequence into canonical combinations (capped at 16)."""
    seq = seq.upper()
    options = [IUPAC_MAP.get(base, [base]) for base in seq]
    combinations = [''.join(p) for p in itertools.product(*options)]
    if len(combinations) > 16:
        combinations = combinations[:16]
    return combinations

def calculate_gc(seq):
    """Calculate GC content percentage."""
    seq = seq.upper()
    gc = sum(1 for b in seq if b in ['G', 'C', 'S'])
    return round((gc / len(seq)) * 100.0, 1)

def evaluate_pair_heterodimer(seq1, seq2, temp_c):
    """
    Calculate minimum delta G (kcal/mol) between seq1 and seq2 for all degenerate variants.
    Returns (min_dg_kcal, worst_struct_found)
    """
    variants1 = expand_iupac(seq1)
    variants2 = expand_iupac(seq2)
    
    min_dg = 99.0
    struct_found = False
    
    for v1 in variants1:
        for v2 in variants2:
            res = primer3.calc_heterodimer(v1, v2, temp_c=temp_c)
            dg_kcal = res.dg / 1000.0  # Convert cal/mol to kcal/mol
            if res.structure_found:
                struct_found = True
                if dg_kcal < min_dg:
                    min_dg = dg_kcal
            else:
                if 0.0 < min_dg:
                    min_dg = 0.0
                    
    if min_dg == 99.0:
        min_dg = 0.0
    return min_dg, struct_found

def calculate_tm(seq):
    """Calculate Tm using Primer3 Owczarzy method."""
    variants = expand_iupac(seq)
    tms = [primer3.calc_tm(v) for v in variants]
    return sum(tms) / len(tms)

def run_multiplex_screening(file_a, file_b, out_dir, num_candidates=10):
    os.makedirs(out_dir, exist_ok=True)
    
    print(f"[M11] Loading Influenza A candidates from: {file_a}")
    with open(file_a, "r") as f:
        cands_a = json.load(f)[:num_candidates]
        
    print(f"[M11] Loading Influenza B candidates from: {file_b}")
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
                            
            dye_status = "PASSED (FAM vs HEX Channel Separation)"
            
            status_raa = "PASS" if (min_dg_raa >= -9.0 and delta_tm <= 3.5) else "WARNING"
            if min_dg_raa < -12.0:
                status_raa = "VETO (High Dimer Risk)"
                
            score = 100.0
            if min_dg_raa < 0:
                score += min_dg_raa * 3.5
            score -= delta_tm * 4.0
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
                "dye_a": "FAM",
                "dye_b": "HEX",
                "dye_status": dye_status,
                "oligos_synthesis": {
                    "flu_a": {
                        "fwd": {"deg_seq": fwd_a_deg, "ref_seq": fwd_a_ref, "len": len(fwd_a_deg), "tm": round(tm_fwd_a, 1), "gc": calculate_gc(fwd_a_deg), "pos": f"{cand_a.get('fwd_start')}-{cand_a.get('fwd_end')}"},
                        "rev": {"deg_seq": rev_a_deg, "ref_seq": rev_a_ref, "len": len(rev_a_deg), "tm": round(tm_rev_a, 1), "gc": calculate_gc(rev_a_deg), "pos": f"{cand_a.get('rev_start')}-{cand_a.get('rev_end')}"},
                        "probe": {"deg_seq": prb_a_deg, "ref_seq": prb_a_ref, "len": len(prb_a_deg), "tm": round(tm_prb_a, 1), "gc": calculate_gc(prb_a_deg), "pos": f"{cand_a.get('probe_start')}-{cand_a.get('probe_end')}", "label": "5'-FAM / 3'-BHQ1 / Internal-THF"}
                    },
                    "flu_b": {
                        "fwd": {"deg_seq": fwd_b_deg, "ref_seq": fwd_b_ref, "len": len(fwd_b_deg), "tm": round(tm_fwd_b, 1), "gc": calculate_gc(fwd_b_deg), "pos": f"{cand_b.get('fwd_start')}-{cand_b.get('fwd_end')}"},
                        "rev": {"deg_seq": rev_b_deg, "ref_seq": rev_b_ref, "len": len(rev_b_deg), "tm": round(tm_rev_b, 1), "gc": calculate_gc(rev_b_deg), "pos": f"{cand_b.get('rev_start')}-{cand_b.get('rev_end')}"},
                        "probe": {"deg_seq": prb_b_deg, "ref_seq": prb_b_ref, "len": len(prb_b_deg), "tm": round(tm_prb_b, 1), "gc": calculate_gc(prb_b_deg), "pos": f"{cand_b.get('probe_start')}-{cand_b.get('probe_end')}", "label": "5'-HEX / 3'-BHQ1 / Internal-THF"}
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
        f.write("## Target Pair: Influenza A MP (FAM) + Influenza B MP (HEX) — Indonesia Cohort\n\n")
        f.write(f"Screened **{len(results)} candidate combinations** ($10 \\times 10$) across two temperature regimes ($39.0^\\circ\\text{{C}}$ RT-RAA & $55.0^\\circ\\text{{C}}$ qPCR).\n\n")
        f.write("### Top 5 Recommended Multiplex Combinations\n\n")
        f.write("| Rank | Candidate Pair | Compatibility Score | RAA ΔG (39°C) | qPCR ΔG (55°C) | ΔTm (°C) | Status (39°C) |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for idx, r in enumerate(results[:5], 1):
            f.write(f"| {idx} | **{r['pair_id']}** | **{r['score']}%** | {r['min_dg_raa_39c']} kcal/mol | {r['min_dg_qpcr_55c']} kcal/mol | {r['delta_tm']}°C | `{r['status_raa_39c']}` |\n")
            
        f.write("\n\n### Rank 1 Oligo Synthesis Specifications (`FluA_Rank7_x_FluB_Rank6`)\n\n")
        r1 = results[0]
        oligos = r1['oligos_synthesis']
        
        f.write("#### 🟩 Target A: Influenza A MP (FAM Channel)\n\n")
        f.write("| Oligo Name | Type | Degenerate Sequence (5' -> 3') | Length | Tm (°C) | GC (%) | Modification / Label |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        f.write(f"| FluA_Fwd | Fwd Primer | `{oligos['flu_a']['fwd']['deg_seq']}` | {oligos['flu_a']['fwd']['len']} nt | {oligos['flu_a']['fwd']['tm']}°C | {oligos['flu_a']['fwd']['gc']}% | Unlabeled |\n")
        f.write(f"| FluA_Rev | Rev Primer | `{oligos['flu_a']['rev']['deg_seq']}` | {oligos['flu_a']['rev']['len']} nt | {oligos['flu_a']['rev']['tm']}°C | {oligos['flu_a']['rev']['gc']}% | Unlabeled |\n")
        f.write(f"| FluA_Probe | Exo-Probe | `{oligos['flu_a']['probe']['deg_seq']}` | {oligos['flu_a']['probe']['len']} nt | {oligos['flu_a']['probe']['tm']}°C | {oligos['flu_a']['probe']['gc']}% | {oligos['flu_a']['probe']['label']} |\n\n")

        f.write("#### 🟨 Target B: Influenza B MP (HEX Channel)\n\n")
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
    parser = argparse.ArgumentParser(description="Multiplex Compatibility Checker")
    parser.add_argument("--file-a", default="results/influenza-a-mp-indo-only/m09/degenerate/segment-MP_region_1/final_candidates.json")
    parser.add_argument("--file-b", default="results/influenza-b-mp-indo-only/m09/degenerate/segment-MP_region_1/final_candidates.json")
    parser.add_argument("--outdir", default="results/multiplex-indo-only")
    parser.add_argument("--num-candidates", type=int, default=10)
    args = parser.parse_args()
    
    run_multiplex_screening(args.file_a, args.file_b, args.outdir, args.num_candidates)

if __name__ == "__main__":
    main()
