import pandas as pd
import json
import os
import yaml

def main():
    # 1. Paths
    profile_path = "ssh-output/results/influenza-a-mp-indo-only/m07/entropy_profile.tsv"
    primers_path = "ssh-output/results/influenza-a-mp-indo-only/m08/degenerate/segment-MP_region_1/best_primers.json"
    manifest_path = "ssh-output/results/influenza-a-mp-indo-only/m07/primer_targets_manifest.yaml"
    
    if not os.path.exists(profile_path) or not os.path.exists(primers_path):
        print("Error: Missing input files.")
        return
        
    # 2. Read entropy profile
    profile = pd.read_csv(profile_path, sep="\t")
    
    # 3. Read best primers
    with open(primers_path, "r") as f:
        primers = json.load(f)
        
    # Let's inspect the forward primer
    fwd = primers.get("forward")
    if not fwd:
        print("Error: Forward primer not found.")
        return
        
    start = fwd["start"]
    end = fwd["end"]
    seq = fwd["sequence"]
    maj_seq = fwd["majority_sequence"]
    
    print("=" * 80)
    print(f"DEGENERATE PRIMER AUDIT REPORT")
    print(f"Primer ID: {fwd['id']} (Forward)")
    print(f"Coordinates: {start} - {end} bp (0-based)")
    print(f"Rescued Sequence:   {seq}")
    print(f"Consensus Sequence: {maj_seq}")
    print("=" * 80)
    
    # IUPAC lookup table
    iupac_map = {
        'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
        'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 
        'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T']
    }
    
    # Extract matching slice of the profile
    # The start/end in the primer dictionary represents the 0-based coordinate on the consensus sequence
    # Let's read the positions
    slice_df = profile[(profile["position"] >= start) & (profile["position"] < end)].copy()
    
    # Let's print out positions where IUPAC code differs from the majority consensus base
    audit_rows = []
    
    for idx, row in slice_df.iterrows():
        pos = int(row["position"])
        cons = row["consensus_base"]
        iupac = row["iupac_base"]
        
        freqs = {
            "A": row["freq_A"],
            "T": row["freq_T"],
            "G": row["freq_G"],
            "C": row["freq_C"]
        }
        
        # Check if degenerate
        is_degenerate = iupac != cons
        
        # Position relative to primer (0-based)
        primer_pos = pos - start
        
        # We want to show the audit
        audit_rows.append({
            "pos": pos,
            "primer_pos": primer_pos,
            "cons": cons,
            "iupac": iupac,
            "A%": freqs["A"] * 100,
            "T%": freqs["T"] * 100,
            "G%": freqs["G"] * 100,
            "C%": freqs["C"] * 100,
            "is_degenerate": is_degenerate
        })
        
    print(f"{'Pos':<6} | {'PrimerPos':<9} | {'Majority':<8} | {'IUPAC':<6} | {'A %':<7} | {'T %':<7} | {'G %':<7} | {'C %':<7} | {'Status/Rationale'}")
    print("-" * 100)
    for r in audit_rows:
        rationale = ""
        if r["is_degenerate"]:
            # Find which bases cleared the 15% threshold
            cleared = []
            for base in ["A", "T", "G", "C"]:
                pct = r[f"{base}%"]
                if pct >= 15.0:
                    cleared.append(f"{base} ({pct:.1f}%)")
            rationale = f"Degenerate! Cleared >=15% threshold: " + ", ".join(cleared)
        else:
            # Check if there is any variation but below 15%
            vars_present = []
            for base in ["A", "T", "G", "C"]:
                pct = r[f"{base}%"]
                if 0.0 < pct < 15.0:
                    vars_present.append(f"{base} ({pct:.1f}%)")
            if vars_present:
                rationale = f"Conserved (minor variation <15% ignored: " + ", ".join(vars_present) + ")"
            else:
                rationale = "100% Conserved"
                
        status_str = f"*{r['iupac']}*" if r["is_degenerate"] else r["iupac"]
        print(f"{r['pos']:<6} | {r['primer_pos']:<9} | {r['cons']:<8} | {status_str:<6} | {r['A%']:<7.1f} | {r['T%']:<7.1f} | {r['G%']:<7.1f} | {r['C%']:<7.1f} | {rationale}")
        
if __name__ == "__main__":
    main()
