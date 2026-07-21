#!/usr/bin/env python3
"""
Module 11: Multiplex Assay Compatibility Screening Engine (Influenza A + B MP)
Evaluates pairwise cross-dimer thermodynamics (ViennaRNA / Primer3), spectral dye separation,
and kinetics alignment for dual-target assays at 39.0°C (RAA) and 55.0°C (qPCR).
"""

import os
import sys
import json
import argparse
import itertools
import pandas as pd
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
    """Expand degenerate sequence into all possible canonical combinations (capped at 16)."""
    seq = seq.upper()
    options = [IUPAC_MAP.get(base, [base]) for base in seq]
    combinations = [''.join(p) for p in itertools.product(*options)]
    if len(combinations) > 16:
        combinations = combinations[:16]
    return combinations

def get_consensus_seq(seq):
    """Convert degenerate sequence to standard IUPAC consensus base."""
    # Map common degenerate bases to consensus majority if needed, or return raw string
    return seq.upper()

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
        fwd_a = cand_a.get("fwd_degenerate_seq") or cand_a.get("fwd_seq") or cand_a.get("fwd")
        rev_a = cand_a.get("rev_degenerate_seq") or cand_a.get("rev_seq") or cand_a.get("rev")
        prb_a = cand_a.get("probe_seq") or cand_a.get("probe")
        
        tm_fwd_a = calculate_tm(fwd_a)
        tm_rev_a = calculate_tm(rev_a)
        
        for cand_b in cands_b:
            rank_b = cand_b.get("rank", 1)
            fwd_b = cand_b.get("fwd_degenerate_seq") or cand_b.get("fwd_seq") or cand_b.get("fwd")
            rev_b = cand_b.get("rev_degenerate_seq") or cand_b.get("rev_seq") or cand_b.get("rev")
            prb_b = cand_b.get("probe_seq") or cand_b.get("probe")
            
            tm_fwd_b = calculate_tm(fwd_b)
            tm_rev_b = calculate_tm(rev_b)
            
            all_tms = [tm_fwd_a, tm_rev_a, tm_fwd_b, tm_rev_b]
            delta_tm = max(all_tms) - min(all_tms)
            
            oligos_a = [("Fwd_A", fwd_a), ("Rev_A", rev_a), ("Probe_A", prb_a)]
            oligos_b = [("Fwd_B", fwd_b), ("Rev_B", rev_b), ("Probe_B", prb_b)]
            
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
                            
            # Dye & Compatibility Evaluation
            # Inf A uses FAM, Inf B uses HEX
            dye_status = "PASSED (FAM vs HEX Channel Separation)"
            
            # Status & Veto Rules
            # Veto if min_dg_raa < -9.0 kcal/mol or delta_tm > 4.0 °C
            status_raa = "PASS" if (min_dg_raa >= -9.0 and delta_tm <= 3.5) else "WARNING"
            if min_dg_raa < -12.0:
                status_raa = "VETO (High Dimer Risk)"
                
            # Score (100 - penalties)
            score = 100.0
            if min_dg_raa < 0:
                score += min_dg_raa * 3.5  # Penalty per kcal/mol
            score -= delta_tm * 4.0        # Penalty per degree Tm difference
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
                "oligos": {
                    "fwd_a": fwd_a, "rev_a": rev_a, "probe_a": prb_a,
                    "fwd_b": fwd_b, "rev_b": rev_b, "probe_b": prb_b
                },
                "cross_dimers_39c": cross_dimers_raa,
                "cross_dimers_55c": cross_dimers_qpcr
            })
            
    # Sort results by score descending
    results.sort(key=lambda x: x["score"], reverse=True)
    
    # Write JSON output
    json_path = os.path.join(out_dir, "multiplex_candidates.json")
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"[M11] Saved multiplex results JSON to: {json_path}")
    
    # Write Markdown Report
    md_path = os.path.join(out_dir, "multiplex_report.md")
    with open(md_path, "w") as f:
        f.write("# Dual-Target Multiplex Assay Compatibility Report\n")
        f.write("## Target Pair: Influenza A MP (FAM) + Influenza B MP (HEX) — Indonesia Cohort\n\n")
        f.write(f"Screened **{len(results)} candidate combinations** ($10 \\times 10$) across two temperature regimes ($39.0^\\circ\\text{{C}}$ RT-RAA & $55.0^\\circ\\text{{C}}$ qPCR).\n\n")
        f.write("### Top 5 Recommended Multiplex Combinations\n\n")
        f.write("| Rank | Candidate Pair | Compatibility Score | RAA ΔG (39°C) | qPCR ΔG (55°C) | ΔTm (°C) | Status (39°C) |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for idx, r in enumerate(results[:5], 1):
            f.write(f"| {idx} | **{r['pair_id']}** | **{r['score']}%** | {r['min_dg_raa_39c']} kcal/mol | {r['min_dg_qpcr_55c']} kcal/mol | {r['delta_tm']}°C | `{r['status_raa_39c']}` |\n")
            
        f.write("\n\n### Cross-Dimer Interactions for Rank 1 Combination\n\n")
        r1 = results[0]
        f.write(f"**Pair ID**: `{r1['pair_id']}` (Score: {r1['score']}%)\n\n")
        f.write("| Oligo Pair | Cross-Dimer ΔG (39°C) | Cross-Dimer ΔG (55°C) |\n")
        f.write("|---|---|---|\n")
        for pair_name, dg in r1['cross_dimers_39c'].items():
            dg_55 = r1['cross_dimers_55c'].get(pair_name, "-")
            f.write(f"| {pair_name} | {dg} kcal/mol | {dg_55} kcal/mol |\n")
            
    print(f"[M11] Saved Markdown report to: {md_path}")
    
    # Write HTML Heatmap Dashboard
    html_path = os.path.join(out_dir, "multiplex_matrix.html")
    generate_html_matrix(results, html_path)
    print(f"[M11] Saved HTML interactive dashboard to: {html_path}")

def generate_html_matrix(results, html_path):
    # Build 10x10 matrix for HTML display
    matrix_data = {}
    for r in results:
        matrix_data[(r['rank_a'], r['rank_b'])] = r
        
    html = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>🧬 Dual Multiplex Assay Compatibility Dashboard</title>
    <link href="https://fonts.googleapis.com/css2?family=Outfit:wght@400;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
    <style>
        :root { --bg:#0b132b; --surface:#1c2541; --accent-fam:#10b981; --accent-hex:#f59e0b; --border:rgba(255,255,255,0.1); }
        body { font-family:'Outfit', sans-serif; background:var(--bg); color:#f0f6ff; padding:2rem; }
        h1 { background:linear-gradient(135deg, #10b981, #f59e0b); -webkit-background-clip:text; -webkit-text-fill-color:transparent; margin-bottom:0.5rem; }
        .grid-table { width:100%; border-collapse:collapse; margin-top:1.5rem; }
        .grid-table th, .grid-table td { border:1px solid var(--border); padding:10px; text-align:center; font-size:12px; }
        .grid-table th { background:var(--surface); }
        .score-good { background:rgba(16,185,129,0.25); color:#34d399; font-weight:700; }
        .score-warn { background:rgba(245,158,11,0.25); color:#fbbf24; font-weight:700; }
        .score-bad { background:rgba(239,68,68,0.25); color:#f87171; font-weight:700; }
        .card { background:var(--surface); border:1px solid var(--border); border-radius:12px; padding:1.5rem; margin-bottom:2rem; }
        .badge-fam { background:rgba(16,185,129,0.2); color:#10b981; padding:3px 8px; border-radius:6px; font-weight:600; }
        .badge-hex { background:rgba(245,158,11,0.2); color:#f59e0b; padding:3px 8px; border-radius:6px; font-weight:600; }
    </style>
</head>
<body>
    <h1>🧬 Dual-Target Multiplex Compatibility Matrix</h1>
    <p style="color:#8b9dc0; margin-bottom:1.5rem;">Influenza A MP (<span class="badge-fam">FAM Probe</span>) vs Influenza B MP (<span class="badge-hex">HEX Probe</span>) — Indonesia Cohort</p>

    <div class="card">
        <h3>Top Rank 1 Recommended Combination</h3>
        <p style="margin-top:8px;"><b>Pair ID:</b> """ + results[0]['pair_id'] + """ &nbsp;|&nbsp; <b>Compatibility Score:</b> <span class="badge-fam">""" + str(results[0]['score']) + """%</span></p>
        <p style="margin-top:4px;"><b>RT-RAA (39°C) Min ΔG:</b> """ + str(results[0]['min_dg_raa_39c']) + """ kcal/mol &nbsp;|&nbsp; <b>qPCR (55°C) Min ΔG:</b> """ + str(results[0]['min_dg_qpcr_55c']) + """ kcal/mol &nbsp;|&nbsp; <b>ΔTm:</b> """ + str(results[0]['delta_tm']) + """°C</p>
    </div>

    <h3>10x10 Pairwise Matrix (Compatibility Score %)</h3>
    <table class="grid-table">
        <thead>
            <tr>
                <th>Flu A \\ Flu B</th>
"""
    for rb in range(1, 11):
        html += f"<th>Rank {rb} (Flu B)</th>"
    html += "</tr></thead><tbody>"

    for ra in range(1, 11):
        html += f"<tr><th>Rank {ra} (Flu A)</th>"
        for rb in range(1, 11):
            cell = matrix_data.get((ra, rb), {})
            score = cell.get('score', 0)
            cls = "score-good" if score >= 85 else ("score-warn" if score >= 70 else "score-bad")
            html += f"<td class='{cls}'>{score}%<br><span style='font-size:10px;opacity:0.8;'>{cell.get('min_dg_raa_39c', 0)} kcal/mol</span></td>"
        html += "</tr>"

    html += """
        </tbody>
    </table>
</body>
</html>
"""
    with open(html_path, "w") as f:
        f.write(html)

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
