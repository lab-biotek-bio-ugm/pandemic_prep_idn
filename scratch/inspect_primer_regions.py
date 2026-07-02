import os
import json

def get_primers_for_run(run_dir):
    degenerate_dir = os.path.join(run_dir, 'm09', 'degenerate')
    if not os.path.exists(degenerate_dir):
        return "Directory not found"
        
    regions = os.listdir(degenerate_dir)
    results = {}
    for r in regions:
        candidates_path = os.path.join(degenerate_dir, r, 'final_candidates.json')
        if os.path.exists(candidates_path):
            with open(candidates_path, 'r') as f:
                cands = json.load(f)
            if cands:
                top_cand = cands[0]
                results[r] = {
                    "fwd_start": top_cand.get("fwd_start"),
                    "fwd_end": top_cand.get("fwd_end"),
                    "rev_start": top_cand.get("rev_start"),
                    "rev_end": top_cand.get("rev_end"),
                    "amplicon_start": top_cand.get("amplicon_start"),
                    "amplicon_end": top_cand.get("amplicon_end"),
                    "amplicon_length": top_cand.get("amplicon_length")
                }
    return results

def main():
    base_dir = 'ssh-output/results'
    runs = ['influenza-a-mp-indo-only', 'influenza-a-mp-sea-only', 'influenza-a-mp-all-partial']
    for r in runs:
        print(f"\n=== Primers designed for {r} ===")
        run_path = os.path.join(base_dir, r)
        primers = get_primers_for_run(run_path)
        print(json.dumps(primers, indent=2))

if __name__ == '__main__':
    main()
