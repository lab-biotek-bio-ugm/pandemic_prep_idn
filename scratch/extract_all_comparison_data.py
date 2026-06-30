import os
import sys
import re
import json

if sys.platform.startswith('win'):
    sys.stdout.reconfigure(encoding='utf-8')

def parse_m10_html(path):
    if not os.path.exists(path):
        return f"File {path} not found"
        
    with open(path, 'r', encoding='utf-8') as f:
        html = f.read()
        
    # Extract nonDegData
    non_deg_match = re.search(r'const nonDegData\s*=\s*(.*?);', html, re.DOTALL)
    deg_match = re.search(r'const degData\s*=\s*(.*?);', html, re.DOTALL)
    
    if not non_deg_match or not deg_match:
        # Try checking different spelling
        non_deg_match = re.search(r'nonDegData\s*=\s*(.*?);', html, re.DOTALL)
        deg_match = re.search(r'degData\s*=\s*(.*?);', html, re.DOTALL)
        
    if not non_deg_match or not deg_match:
        return "Could not find nonDegData or degData variables in HTML script"
        
    try:
        non_deg_data = json.loads(non_deg_match.group(1).strip())
        deg_data = json.loads(deg_match.group(1).strip())
    except Exception as e:
        return f"JSON Parse error: {e}"
        
    # Get details for Rank 1 (usually index 0 or cand_id ending in _Cand_1)
    cand_1_non_deg = next((c for c in non_deg_data if c.get('rank') == 1 or 'Cand_1' in c.get('cand_id', '')), None)
    cand_1_deg = next((c for c in deg_data if c.get('rank') == 1 or 'Cand_1' in c.get('cand_id', '')), None)
    
    if not cand_1_non_deg or not cand_1_deg:
        return "Could not find Candidate 1 in non-degenerate or degenerate lists"
        
    res = {
        "non_degenerate": {
            "pcr_pct": cand_1_non_deg.get("pcr_amplified_pct"),
            "raa_pct": cand_1_non_deg.get("raa_amplified_pct"),
            "design_db_pcr_pct": cand_1_non_deg.get("design_db_amplified_pct"),
            "design_db_raa_pct": cand_1_non_deg.get("design_db_raa_amp_pct"),
            "subtypes": cand_1_non_deg.get("design_db_breakdowns", {}).get("subtype", [])
        },
        "degenerate": {
            "pcr_pct": cand_1_deg.get("pcr_amplified_pct"),
            "raa_pct": cand_1_deg.get("raa_amplified_pct"),
            "design_db_pcr_pct": cand_1_deg.get("design_db_amplified_pct"),
            "design_db_raa_pct": cand_1_deg.get("design_db_raa_amp_pct"),
            "subtypes": cand_1_deg.get("design_db_breakdowns", {}).get("subtype", [])
        }
    }
    return res

def main():
    runs = ['influenza-a-mp-indo-only', 'influenza-a-mp-sea-only', 'influenza-a-mp-all-partial']
    for r in runs:
        path = f"ssh-output/results/{r}/m10/comparison_report.html"
        print(f"\n==========================================")
        print(f"JALAN RUN: {r}")
        print(f"==========================================")
        res = parse_m10_html(path)
        if isinstance(res, str):
            print(res)
        else:
            print("OVERALL COVERAGE ON GLOBAL DATASET (N=50,750):")
            print(f"  - Non-Degenerate: PCR {res['non_degenerate']['pcr_pct']}% | RAA {res['non_degenerate']['raa_pct']}%")
            print(f"  - Degenerate:     PCR {res['degenerate']['pcr_pct']}% | RAA {res['degenerate']['raa_pct']}%")
            print(f"  - Uplift (Gain):  PCR {round(res['degenerate']['pcr_pct'] - res['non_degenerate']['pcr_pct'], 2)}%")
            
            print("LOCAL DESIGN DATABASE COVERAGE:")
            print(f"  - Non-Degenerate: PCR {res['non_degenerate']['design_db_pcr_pct']}% | RAA {res['non_degenerate']['design_db_raa_pct']}%")
            print(f"  - Degenerate:     PCR {res['degenerate']['design_db_pcr_pct']}% | RAA {res['degenerate']['design_db_raa_pct']}%")
            
            print("SUBTYPE COVERAGE (LOCAL DESIGN DATABASE):")
            for sub in res['degenerate']['subtypes'][:4]: # Print top 4 subtypes
                sub_name = sub.get('subtype')
                sub_total = sub.get('total')
                sub_deg_pct = sub.get('pcr_amplified_pct')
                
                # Find matching non-deg subtype
                non_deg_sub = next((s for s in res['non_degenerate']['subtypes'] if s.get('subtype') == sub_name), {})
                sub_non_deg_pct = non_deg_sub.get('pcr_amplified_pct', 0.0)
                
                print(f"  - {sub_name} (N={sub_total}): Non-Degen {sub_non_deg_pct}% -> Degen {sub_deg_pct}%")

if __name__ == '__main__':
    main()
