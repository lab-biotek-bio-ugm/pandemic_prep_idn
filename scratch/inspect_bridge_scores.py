import os
import csv
import yaml

def get_bridge_details(run_dir):
    summary_path = os.path.join(run_dir, 'm07', 'candidate_targets_summary.csv')
    manifest_path = os.path.join(run_dir, 'm07', 'primer_targets_manifest.yaml')
    
    details = {}
    if os.path.exists(summary_path):
        with open(summary_path, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            if rows:
                # Let's get the first row details
                r0 = rows[0]
                details['summary'] = {
                    "region_id": r0.get("region_id"),
                    "composite_score": r0.get("composite_score"),
                    "grade": r0.get("grade"),
                    "mean_entropy": r0.get("mean_entropy"),
                    "conservation_pct": r0.get("conservation_pct")
                }
    if os.path.exists(manifest_path):
        with open(manifest_path, 'r') as f:
            manifest = yaml.safe_load(f)
        candidate_regions = manifest.get("candidate_regions", [])
        if candidate_regions:
            r0 = candidate_regions[0]
            details['manifest'] = {
                "region_id": r0.get("region_id"),
                "composite_score": r0.get("composite_score"),
                "grade": r0.get("grade"),
                "start": r0.get("conserved_zone", {}).get("start"),
                "end": r0.get("conserved_zone", {}).get("end"),
                "mean_entropy": r0.get("biomechanics", {}).get("mean_entropy"),
                "conservation_pct": r0.get("biomechanics", {}).get("conservation_pct")
            }
    return details

def main():
    base_dir = 'ssh-output/results'
    runs = ['influenza-a-mp-indo-only', 'influenza-a-mp-sea-only', 'influenza-a-mp-all-partial']
    for r in runs:
        print(f"\n=== Bridge details for {r} ===")
        run_path = os.path.join(base_dir, r)
        details = get_bridge_details(run_path)
        print(yaml.dump(details))

if __name__ == '__main__':
    main()
