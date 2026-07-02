import os
import json
import yaml

def get_top_3_candidates(run_dir):
    degenerate_dir = os.path.join(run_dir, 'm09', 'degenerate')
    if not os.path.exists(degenerate_dir):
        return None
        
    regions = os.listdir(degenerate_dir)
    results = {}
    for r in regions:
        candidates_path = os.path.join(degenerate_dir, r, 'final_candidates.json')
        if os.path.exists(candidates_path):
            with open(candidates_path, 'r') as f:
                cands = json.load(f)
            # Let's get the top 3
            results[r] = []
            for cand in cands[:3]:
                results[r].append({
                    "rank": cand.get("rank"),
                    "fwd": cand.get("fwd_degenerate_seq"),
                    "rev": cand.get("rev_degenerate_seq"),
                    "probe": cand.get("probe_degenerate_seq"),
                    "fwd_tm": cand.get("fwd_tm"),
                    "rev_tm": cand.get("rev_tm"),
                    "probe_tm": cand.get("probe_tm"),
                    "amplicon_length": cand.get("amplicon_length"),
                    "human_score": cand.get("human_score"),
                    "virus_score": cand.get("virus_score"),
                    "hairpin_dg": cand.get("hairpin_dg"),
                    "cross_dimer_dg": cand.get("cross_dimer_dg")
                })
    return results

def main():
    base_dir = 'ssh-output/results'
    runs = ['influenza-a-mp-indo-only', 'influenza-a-mp-sea-only', 'influenza-a-mp-all-partial']
    all_data = {}
    for r in runs:
        run_path = os.path.join(base_dir, r)
        all_data[r] = get_top_3_candidates(run_path)
        
    with open('scratch/top_3_candidates.json', 'w') as f:
        json.dump(all_data, f, indent=2)
    print("Top 3 candidates extracted successfully!")

if __name__ == '__main__':
    main()
