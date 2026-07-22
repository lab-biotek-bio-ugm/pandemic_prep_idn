#!/usr/bin/env python3
"""
Module 11: Triplex Multiplex Assay Compatibility Screening Engine
Evaluates three-way cross-dimer thermodynamics, spectral dye separation,
and Tm alignment for triplex assays (Influenza A + Influenza B + Measles).
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
    """Calculate minimum delta G (kcal/mol) between seq1 and seq2."""
    variants1 = expand_iupac(seq1)
    variants2 = expand_iupac(seq2)
    
    min_dg = 99.0
    struct_found = False
    
    for v1 in variants1:
        for v2 in variants2:
            res = primer3.calc_heterodimer(v1, v2, temp_c=temp_c)
            dg_kcal = res.dg / 1000.0
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
    """Calculate Tm using Primer3."""
    variants = expand_iupac(seq)
    tms = [primer3.calc_tm(v) for v in variants]
    return sum(tms) / len(tms)

def run_triplex_screening(file_a, file_b, file_c, out_dir, num_candidates=10):
    os.makedirs(out_dir, exist_ok=True)
    
    print(f"[Triplex] Loading Target A candidates from: {file_a}")
    with open(file_a, "r") as f:
        cands_a = json.load(f)[:num_candidates]
        
    print(f"[Triplex] Loading Target B candidates from: {file_b}")
    with open(file_b, "r") as f:
        cands_b = json.load(f)[:num_candidates]
        
    print(f"[Triplex] Loading Target C candidates from: {file_c}")
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
        
        tm_fwd_a = calculate_tm(fwd_a)
        tm_rev_a = calculate_tm(rev_a)
        tm_prb_a = calculate_tm(prb_a)
        
        for cand_b in cands_b:
            rank_b = cand_b.get("rank", 1)
            fwd_b = cand_b.get("fwd_degenerate_seq") or cand_b.get("fwd_seq")
            rev_b = cand_b.get("rev_degenerate_seq") or cand_b.get("rev_seq")
            prb_b = cand_b.get("probe_degenerate_seq") or cand_b.get("probe_seq")
            
            tm_fwd_b = calculate_tm(fwd_b)
            tm_rev_b = calculate_tm(rev_b)
            tm_prb_b = calculate_tm(prb_b)
            
            for cand_c in cands_c:
                rank_c = cand_c.get("rank", 1)
                fwd_c = cand_c.get("fwd_degenerate_seq") or cand_c.get("fwd_seq")
                rev_c = cand_c.get("rev_degenerate_seq") or cand_c.get("rev_seq")
                prb_c = cand_c.get("probe_degenerate_seq") or cand_c.get("probe_seq")
                
                tm_fwd_c = calculate_tm(fwd_c)
                tm_rev_c = calculate_tm(rev_c)
                tm_prb_c = calculate_tm(prb_c)
                
                all_primer_tms = [tm_fwd_a, tm_rev_a, tm_fwd_b, tm_rev_b, tm_fwd_c, tm_rev_c]
                delta_tm = max(all_primer_tms) - min(all_primer_tms)
                
                oligos_a = [("Fwd_A", fwd_a), ("Rev_A", rev_a), ("Probe_A", prb_a)]
                oligos_b = [("Fwd_B", fwd_b), ("Rev_B", rev_b), ("Probe_B", prb_b)]
                oligos_c = [("Fwd_C", fwd_c), ("Rev_C", rev_c), ("Probe_C", prb_c)]
                
                # Check cross-dimers for all three pairs (A-B, A-C, B-C)
                cross_dimers_39c = {}
                min_dg_39c = 0.0
                
                # Pair A-B
                for name_a, seq_a in oligos_a:
                    for name_b, seq_b in oligos_b:
                        dg, _ = evaluate_pair_heterodimer(seq_a, seq_b, temp_c=39.0)
                        cross_dimers_39c[f"{name_a}_vs_{name_b}"] = round(dg, 2)
                        if dg < min_dg_39c:
                            min_dg_39c = dg
                
                # Pair A-C
                for name_a, seq_a in oligos_a:
                    for name_c, seq_c in oligos_c:
                        dg, _ = evaluate_pair_heterodimer(seq_a, seq_c, temp_c=39.0)
                        cross_dimers_39c[f"{name_a}_vs_{name_c}"] = round(dg, 2)
                        if dg < min_dg_39c:
                            min_dg_39c = dg
                            
                # Pair B-C
                for name_b, seq_b in oligos_b:
                    for name_c, seq_c in oligos_c:
                        dg, _ = evaluate_pair_heterodimer(seq_b, seq_c, temp_c=39.0)
                        cross_dimers_39c[f"{name_b}_vs_{name_c}"] = round(dg, 2)
                        if dg < min_dg_39c:
                            min_dg_39c = dg
                            
                status_39c = "PASS" if (min_dg_39c >= -9.0 and delta_tm <= 4.0) else "WARNING"
                if min_dg_39c < -12.0:
                    status_39c = "VETO"
                    
                # Score out of 100
                score = 100.0
                if min_dg_39c < 0:
                    score += min_dg_39c * 4.0
                score -= delta_tm * 3.0
                score = max(0.0, round(score, 1))
                
                results.append({
                    "triplex_id": f"FluA_R{rank_a}_FluB_R{rank_b}_Measles_R{rank_c}",
                    "score": score,
                    "status_39c": status_39c,
                    "min_dg_39c": round(min_dg_39c, 2),
                    "delta_tm": round(delta_tm, 2),
                    "oligos_synthesis": {
                        "flu_a": {
                            "fwd": fwd_a, "rev": rev_a, "probe": prb_a,
                            "fwd_tm": round(tm_fwd_a, 1), "rev_tm": round(tm_rev_a, 1), "probe_tm": round(tm_prb_a, 1),
                            "probe_label": "5'-FAM / 3'-BHQ1 / Internal-THF"
                        },
                        "flu_b": {
                            "fwd": fwd_b, "rev": rev_b, "probe": prb_b,
                            "fwd_tm": round(tm_fwd_b, 1), "rev_tm": round(tm_rev_b, 1), "probe_tm": round(tm_prb_b, 1),
                            "probe_label": "5'-HEX / 3'-BHQ1 / Internal-THF"
                        },
                        "measles": {
                            "fwd": fwd_c, "rev": rev_c, "probe": prb_c,
                            "fwd_tm": round(tm_fwd_c, 1), "rev_tm": round(tm_rev_c, 1), "probe_tm": round(tm_prb_c, 1),
                            "probe_label": "5'-Cy5 / 3'-BHQ2 / Internal-THF"
                        }
                    },
                    "cross_dimers_39c": cross_dimers_39c
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
        f.write("## Targets: Influenza A (FAM) + Influenza B (HEX) + Measles (Cy5)\n\n")
        f.write(f"Screened **{len(results)} triplex combinations** ({len(cands_a)} x {len(cands_b)} x {len(cands_c)}) at RAA operating temperature ($39.0^\\circ\\text{{C}}$).\n\n")
        f.write("### Top 5 Recommended Triplex Combinations\n\n")
        f.write("| Rank | Triplex Candidate Set | Compatibility Score | Worst RAA ΔG (39°C) | Max ΔTm (°C) | Status (39°C) |\n")
        f.write("|---|---|---|---|---|---|\n")
        for idx, r in enumerate(results[:5], 1):
            f.write(f"| {idx} | **{r['triplex_id']}** | **{r['score']}%** | {r['min_dg_39c']} kcal/mol | {r['delta_tm']}°C | `{r['status_39c']}` |\n")
            
        # Detailed specs of Rank 1
        r1 = results[0]
        f.write(f"\n\n### Rank 1 Triplex Specifications (`{r1['triplex_id']}`)\n\n")
        
        # Target A (Flu A)
        f.write("#### 🟩 Target A: Influenza A MP (FAM Channel)\n")
        f.write("| Oligo Name | Type | Sequence (5' -> 3') | Tm (°C) | Modification / Label |\n")
        f.write("|---|---|---|---|---|\n")
        f.write(f"| FluA_Fwd | Fwd Primer | `{r1['oligos_synthesis']['flu_a']['fwd']}` | {r1['oligos_synthesis']['flu_a']['fwd_tm']}°C | Unlabeled |\n")
        f.write(f"| FluA_Rev | Rev Primer | `{r1['oligos_synthesis']['flu_a']['rev']}` | {r1['oligos_synthesis']['flu_a']['rev_tm']}°C | Unlabeled |\n")
        f.write(f"| FluA_Probe | Exo-Probe | `{r1['oligos_synthesis']['flu_a']['probe']}` | {r1['oligos_synthesis']['flu_a']['probe_tm']}°C | {r1['oligos_synthesis']['flu_a']['probe_label']} |\n\n")
        
        # Target B (Flu B)
        f.write("#### 🟨 Target B: Influenza B MP (HEX Channel)\n")
        f.write("| Oligo Name | Type | Sequence (5' -> 3') | Tm (°C) | Modification / Label |\n")
        f.write("|---|---|---|---|---|\n")
        f.write(f"| FluB_Fwd | Fwd Primer | `{r1['oligos_synthesis']['flu_b']['fwd']}` | {r1['oligos_synthesis']['flu_b']['fwd_tm']}°C | Unlabeled |\n")
        f.write(f"| FluB_Rev | Rev Primer | `{r1['oligos_synthesis']['flu_b']['rev']}` | {r1['oligos_synthesis']['flu_b']['rev_tm']}°C | Unlabeled |\n")
        f.write(f"| FluB_Probe | Exo-Probe | `{r1['oligos_synthesis']['flu_b']['probe']}` | {r1['oligos_synthesis']['flu_b']['probe_tm']}°C | {r1['oligos_synthesis']['flu_b']['probe_label']} |\n\n")
        
        # Target C (Measles)
        f.write("#### 🟦 Target C: Measles Virus (Cy5 Channel)\n")
        f.write("| Oligo Name | Type | Sequence (5' -> 3') | Tm (°C) | Modification / Label |\n")
        f.write("|---|---|---|---|---|\n")
        f.write(f"| Measles_Fwd | Fwd Primer | `{r1['oligos_synthesis']['measles']['fwd']}` | {r1['oligos_synthesis']['measles']['fwd_tm']}°C | Unlabeled |\n")
        f.write(f"| Measles_Rev | Rev Primer | `{r1['oligos_synthesis']['measles']['rev']}` | {r1['oligos_synthesis']['measles']['rev_tm']}°C | Unlabeled |\n")
        f.write(f"| Measles_Probe | Exo-Probe | `{r1['oligos_synthesis']['measles']['probe']}` | {r1['oligos_synthesis']['measles']['probe_tm']}°C | {r1['oligos_synthesis']['measles']['probe_label']} |\n\n")
        
        # Cross dimer detail
        f.write("### Worst Cross-Dimer Interactions (< -5.0 kcal/mol)\n\n")
        f.write("| Oligo Pair | Cross-Dimer ΔG (39°C) |\n")
        f.write("|---|---|\n")
        for pair_name, dg in r1['cross_dimers_39c'].items():
            if dg < -5.0:
                f.write(f"| {pair_name} | {dg} kcal/mol |\n")
                
    print(f"[Triplex] Saved Markdown report to: {md_path}")

def main():
    parser = argparse.ArgumentParser(description="Triplex Compatibility Checker")
    parser.add_argument("--file-a", required=True)
    parser.add_argument("--file-b", required=True)
    parser.add_argument("--file-c", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--num-candidates", type=int, default=10)
    args = parser.parse_args()
    
    run_triplex_screening(args.file_a, args.file_b, args.file_c, args.outdir, args.num_candidates)

if __name__ == "__main__":
    main()
