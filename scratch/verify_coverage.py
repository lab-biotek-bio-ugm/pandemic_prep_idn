import json

def check_run(run):
    dg = json.load(open(f'ssh-output/results/{run}/m10/degenerate/coverage_summary.json'))
    nd = json.load(open(f'ssh-output/results/{run}/m10/non_degenerate/coverage_summary.json'))
    cand1_dg = next((x for x in dg if x['rank'] == 1), dg[0])
    cand1_nd = next((x for x in nd if x['rank'] == 1), nd[0])
    
    print(f"=== {run} ===")
    print(f"  NON-DEGENERATE: PCR={cand1_nd['pcr_amplified_pct']}%  RAA={cand1_nd['raa_amplified_pct']}%")
    print(f"  DEGENERATE:     PCR={cand1_dg['pcr_amplified_pct']}%  RAA={cand1_dg['raa_amplified_pct']}%")
    
    uplift_pcr = round(cand1_dg['pcr_amplified_pct'] - cand1_nd['pcr_amplified_pct'], 2)
    uplift_raa = round(cand1_dg['raa_amplified_pct'] - cand1_nd['raa_amplified_pct'], 2)
    print(f"  UPLIFT:         PCR={uplift_pcr}%  RAA={uplift_raa}%")
    
    print("  Subtype breakdown (DEGENERATE):")
    for st in cand1_dg.get('design_db_breakdowns', {}).get('subtype', []):
        print(f"    {st['subtype']}: PCR={st['pcr_amplified_pct']}%, RAA_amp={st['raa_amplified_pct']}%, RAA_det={st['raa_detected_pct']}%")
    print()

for run in ['influenza-a-mp-indo-only', 'influenza-a-mp-sea-only', 'influenza-a-mp-all-partial']:
    check_run(run)
