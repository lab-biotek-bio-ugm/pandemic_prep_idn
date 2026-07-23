#!/usr/bin/env python3
import math
import argparse
import sys
import os
import json
from pathlib import Path
import pandas as pd
import yaml
from Bio import SeqIO

# Add the scripts directory to path to import m07 logic if needed, 
# but for safety and speed, we'll redefine the core scoring math here.
# This prevents potential side-effects from the main script's GFF/JSON parsing loops.

def calculate_composite_score(cons_score, fubar_score, meme_score, recomb_penalty, weights):
    """Calculate composite score based on weights."""
    score = (weights['conservation'] * cons_score +
             weights['purifying'] * fubar_score +
             weights['episodic_penalty'] * meme_score -
             weights['recombination'] * recomb_penalty)
    return max(0.0, min(1.0, score))

def main():
    parser = argparse.ArgumentParser(description="M07 Weight Sensitivity Sweep")
    parser.add_argument("--report", required=True, help="Path to a baseline bridge_report.json")
    parser.add_argument("--outdir", required=True, help="Output directory for sensitivity results")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. Load baseline data
    if not os.path.exists(args.report):
        sys.exit(f"ERROR: Report file {args.report} not found.")
    
    with open(args.report, "r") as f:
        master_data = json.load(f)
    
    regions = master_data.get("regions", [])
    if not regions:
        sys.exit("ERROR: No regions found in report.")

    print(f"[Sensitivity] Loaded {len(regions)} regions from baseline report.")

    # 2. Define Sweep Range (Conservation weight from 0.3 to 0.7)
    # Original Ratios (P:E:R = 0.3 : 0.1 : 0.2 -> 3:1:2)
    cons_weights = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70]
    
    all_results = []

    for cw in cons_weights:
        remainder = 1.0 - cw
        # Distribute remainder in 3:1:2 ratio
        pw = remainder * (0.3 / 0.6)
        ew = remainder * (0.1 / 0.6)
        rw = remainder * (0.2 / 0.6)
        
        current_weights = {
            'conservation': cw,
            'purifying': pw,
            'episodic_penalty': ew,
            'recombination': rw
        }
        
        print(f"  > Simulating W_cons={cw:.2f} (P={pw:.2f}, E={ew:.2f}, R={rw:.2f})")

        iteration_top_regions = []
        for reg in regions:
            # Prepare scores
            # Normalise by the 5-symbol alphabet maximum (log2(5)), matching
            # MAX_ENTROPY_BITS in m07_bridge_builder.py. Without this, regions whose
            # mean entropy exceeds 1.0 (reachable now that gaps count) produce a
            # NEGATIVE conservation score and an uninterpretable composite.
            cons_score = max(0.0, 1.0 - reg['mean_entropy'] / math.log2(5))
            fubar_score = reg['avg_fubar_purifying_prob']
            meme_score = 1.0 - reg['meme_danger_pct']
            recomb_penalty = 1.0 if reg['flags'].get('near_recombination', False) else 0.0
            
            score = calculate_composite_score(cons_score, fubar_score, meme_score, recomb_penalty, current_weights)
            
            # Grade
            if score >= 0.90: grade = "A"
            elif score >= 0.75: grade = "B"
            elif score >= 0.60: grade = "C"
            else: grade = "D"
            
            iteration_top_regions.append({
                "id": f"{reg['gene']}_{reg['start']}_{reg['end']}",
                "score": score,
                "grade": grade,
                "cw": cw
            })
        
        # Rank them
        iteration_top_regions.sort(key=lambda x: x['score'], reverse=True)
        for rank, item in enumerate(iteration_top_regions):
            item['rank'] = rank + 1
            all_results.append(item)

    # 3. Analyze Stability
    df = pd.DataFrame(all_results)
    
    # Calculate Frequency in Top 5
    top5_counts = df[df['rank'] <= 5].groupby('id').size().reset_index(name='top5_frequency')
    top5_counts['stability_score'] = (top5_counts['top5_frequency'] / len(cons_weights)) * 100
    
    # Calculate Average Rank and Score
    avg_stats = df.groupby('id').agg({
        'score': 'mean',
        'rank': 'mean'
    }).reset_index()
    avg_stats.rename(columns={'score': 'avg_score', 'rank': 'avg_rank'}, inplace=True)

    # Merge
    stability_report = pd.merge(top5_counts, avg_stats, on='id', how='outer').fillna(0)
    stability_report = stability_report.sort_values(by=['stability_score', 'avg_score'], ascending=False)

    # 4. Save results
    report_path = outdir / "sensitivity_stability_report.tsv"
    stability_report.to_csv(report_path, sep="\t", index=False)
    
    detailed_path = outdir / "sensitivity_detailed_iterations.csv"
    df.to_csv(detailed_path, index=False)

    print(f"\n[Sensitivity] SUCCESS!")
    print(f"[Sensitivity] Stability Report: {report_path}")
    print(f"[Sensitivity] Detailed Logs: {detailed_path}")
    
    # Quick Print of Top Stable Regions
    print("\nTOP STABLE REGIONS (Persistence Score):")
    print(stability_report.head(10).to_string(index=False))

if __name__ == "__main__":
    main()
