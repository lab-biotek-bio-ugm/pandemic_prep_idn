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

def durability_multiplier(rel_rate, floor, rate_cap=3.0, no_evidence=0.95):
    """Mirror of compute_durability_multiplier() in m07_bridge_builder.py."""
    if rel_rate is None:
        return no_evidence
    if rel_rate <= 1.0:
        return 1.0
    if rel_rate >= rate_cap:
        return floor
    return 1.0 - ((rel_rate - 1.0) / (rate_cap - 1.0)) * (1.0 - floor)


def calculate_composite_score(cons_score, rel_rate, meme_danger_sites, near_recomb,
                              floor, meme_cfg=None, recomb_penalty=0.5):
    """Reproduce the scorer's actual model: Core x Multipliers.

    This tool previously evaluated a DIFFERENT function from the one M07 uses -- a
    four-term purely additive score (conservation + purifying + episodic -
    recombination) with its own grade cut-offs. Any sensitivity result it produced
    described a model the pipeline does not run, which made it useless for exactly
    the job it exists for: justifying the parameters.
    """
    m_dur = durability_multiplier(rel_rate, floor)
    m_meme = 1.0
    for t in (meme_cfg or [{"max_sites": 0, "multiplier": 1.0},
                           {"max_sites": 1, "multiplier": 0.9},
                           {"max_sites": 3, "multiplier": 0.7}]):
        if meme_danger_sites <= t["max_sites"]:
            m_meme = t["multiplier"]
            break
    else:
        m_meme = 0.0  # veto_above
    m_rec = recomb_penalty if near_recomb else 1.0
    return max(0.0, min(1.0, cons_score * m_dur * m_meme * m_rec))

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

    # 2. Sweep the durability floor -- the model's only free parameter now that the
    #    conservation/purifying weights are gone. floor = 1.0 means "ignore
    #    substitution rate entirely", the supervisor's conservation-only baseline;
    #    lower values let nucleotide turnover discount a region more heavily.
    floors = [1.00, 0.99, 0.97, 0.95, 0.92, 0.90, 0.85, 0.80]

    all_results = []

    for floor in floors:
        print(f"  > Simulating durability floor={floor:.2f} "
              f"(may overturn a conservation gap of up to {(1-floor)*100:.0f} pp)")

        iteration_top_regions = []
        for reg in regions:
            # Prepare scores
            # Normalise by the 5-symbol alphabet maximum (log2(5)), matching
            # MAX_ENTROPY_BITS in m07_bridge_builder.py. Without this, regions whose
            # mean entropy exceeds 1.0 (reachable now that gaps count) produce a
            # NEGATIVE conservation score and an uninterpretable composite.
            cons_score = max(0.0, 1.0 - reg['mean_entropy'] / math.log2(5))
            rel_rate = reg.get('rel_subst_rate')
            if rel_rate in ("", None):
                rel_rate = None
            near_recomb = reg['flags'].get('near_recombination', False)
            meme_sites = reg.get('meme_episodic_sites', 0)

            score = calculate_composite_score(cons_score, rel_rate, meme_sites,
                                              near_recomb, floor)

            # Same cut-offs as the scorer, so grades here mean what they mean there.
            if score >= 0.95: grade = "S"
            elif score >= 0.90: grade = "A+"
            elif score >= 0.85: grade = "A"
            elif score >= 0.80: grade = "A-"
            elif score >= 0.75: grade = "B+"
            elif score >= 0.70: grade = "B"
            elif score >= 0.60: grade = "C"
            else: grade = "D"

            iteration_top_regions.append({
                "id": f"{reg['gene']}_{reg['start']}_{reg['end']}",
                "score": score,
                "grade": grade,
                "floor": floor
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
    top5_counts['stability_score'] = (top5_counts['top5_frequency'] / len(floors)) * 100
    
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
