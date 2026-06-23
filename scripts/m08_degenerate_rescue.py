#!/usr/bin/env python3
"""
M08 Degenerate Rescue Script
Back-maps non-degenerate primers designed on majority consensus sequence to the IUPAC consensus sequence.
Calculates worst-case thermodynamic parameters (Tm range, hairpin, dimers) using primer3-py.
"""

import argparse
import json
import os
import sys
import yaml
import itertools
from pathlib import Path
import primer3

# IUPAC degeneracy nucleotide mapping
IUPAC_DICT = {
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
    'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 
    'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}

def get_degeneracy_multiplier(sequence: str) -> int:
    multiplier = 1
    for char in sequence.upper():
        if char in IUPAC_DICT:
            multiplier *= len(IUPAC_DICT[char])
    return multiplier

def enumerate_sequences(sequence: str) -> list:
    bases_lists = [IUPAC_DICT.get(c, [c]) for c in sequence.upper()]
    return [''.join(p) for p in itertools.product(*bases_lists)]

def resolve_excess_degeneracy(degenerate_seq: str, original_seq: str, max_multiplier: int) -> str:
    """
    If degeneracy multiplier is greater than max_multiplier, resolve the positions with the highest
    multiplicity back to the original non-degenerate (majority consensus) base, until multiplier is safe.
    """
    if get_degeneracy_multiplier(degenerate_seq) <= max_multiplier:
        return degenerate_seq

    # Find degenerate positions
    degenerate_positions = []
    for idx, (deg_char, orig_char) in enumerate(zip(degenerate_seq, original_seq)):
        if deg_char in IUPAC_DICT:
            mult = len(IUPAC_DICT[deg_char])
            degenerate_positions.append((idx, deg_char, orig_char, mult))

    # Sort degenerate positions by multiplicity descending, then index
    # We resolve highest multiplicity first (e.g. N -> 4, B/D/H/V -> 3, then R/Y/S/W/K/M -> 2)
    degenerate_positions.sort(key=lambda x: (-x[3], x[0]))

    seq_list = list(degenerate_seq)
    for idx, deg_char, orig_char, mult in degenerate_positions:
        # Resolve to original base
        seq_list[idx] = orig_char
        new_seq = "".join(seq_list)
        if get_degeneracy_multiplier(new_seq) <= max_multiplier:
            return new_seq

    return "".join(seq_list)

def calculate_gc(sequence: str) -> float:
    seq_upper = sequence.upper()
    g_c_count = seq_upper.count('G') + seq_upper.count('C')
    total = len(seq_upper)
    return round((g_c_count / total) * 100, 2) if total > 0 else 0.0

def evaluate_thermodynamics(possible_seqs: list, thermo_mode: str = "worst_case") -> dict:
    """
    Calculates Tm range, GC range, and worst-case hairpin / self-dimer dG for a list of possible sequences.
    """
    tms = []
    gcs = []
    hairpin_dgs = []
    homodimer_dgs = []

    for seq in possible_seqs:
        # Tm calculation using primer3
        try:
            tm = primer3.calc_tm(seq)
            tms.append(tm)
        except Exception:
            # Fallback Tm
            w_a = seq.count('A') + seq.count('T')
            g_c = seq.count('G') + seq.count('C')
            tms.append(2 * w_a + 4 * g_c)

        gcs.append(calculate_gc(seq))

        # Hairpin using primer3 (cal to kcal)
        try:
            hp = primer3.calc_hairpin(seq).dg / 1000.0
            hairpin_dgs.append(hp)
        except Exception:
            hairpin_dgs.append(0.0)

        # Homodimer using primer3
        try:
            hd = primer3.calc_homodimer(seq).dg / 1000.0
            homodimer_dgs.append(hd)
        except Exception:
            homodimer_dgs.append(0.0)

    # Heterodimers between different sequences in the same pool (if any)
    cross_self_dgs = []
    if len(possible_seqs) > 1:
        for seq1, seq2 in itertools.combinations(possible_seqs, 2):
            try:
                ht = primer3.calc_heterodimer(seq1, seq2).dg / 1000.0
                cross_self_dgs.append(ht)
            except Exception:
                cross_self_dgs.append(0.0)

    # Determine aggregated values
    min_tm = min(tms)
    max_tm = max(tms)
    avg_tm = sum(tms) / len(tms)

    min_gc = min(gcs)
    max_gc = max(gcs)
    avg_gc = sum(gcs) / len(gcs)

    # For dG, the "worst" is the most negative (most stable secondary structure)
    worst_hairpin = min(hairpin_dgs) if hairpin_dgs else 0.0
    worst_homodimer = min(homodimer_dgs) if homodimer_dgs else 0.0
    if cross_self_dgs:
        worst_homodimer = min(worst_homodimer, min(cross_self_dgs))

    return {
        "tm_min": min_tm,
        "tm_max": max_tm,
        "tm_avg": avg_tm,
        "gc_min": min_gc,
        "gc_max": max_gc,
        "gc_avg": avg_gc,
        "worst_hairpin_dg": worst_hairpin,
        "worst_homodimer_dg": worst_homodimer,
        "tm_str": f"{min_tm:.2f} - {max_tm:.2f} °C" if min_tm != max_tm else f"{min_tm:.2f} °C",
        "gc_str": f"{min_gc:.2f} - {max_gc:.2f} %" if min_gc != max_gc else f"{min_gc:.2f} %"
    }

def calculate_cross_heterodimer_worst(seqs1: list, seqs2: list) -> float:
    """Calculates the worst-case (most negative) heterodimer dG between two pools."""
    worst_dg = 0.0
    for seq1 in seqs1:
        for seq2 in seqs2:
            try:
                hd = primer3.calc_heterodimer(seq1, seq2).dg / 1000.0
                if hd < worst_dg:
                    worst_dg = hd
            except Exception:
                pass
    return worst_dg

def rescue_primer(primer_dict: dict, iupac_seq: str, consensus_seq: str, max_multiplier: int, thermo_mode: str) -> dict:
    if primer_dict is None:
        return None

    start = primer_dict.get("start")
    end = primer_dict.get("end")
    if start is None or end is None:
        return primer_dict

    original_seq = primer_dict["sequence"]
    
    # Get slice from IUPAC consensus sequence
    iupac_slice = iupac_seq[start:end]
    consensus_slice = consensus_seq[start:end]

    # Detect orientation of original primer relative to the reference genome slice
    # If original_seq is closer to the reverse complement of consensus_slice, it's reverse-oriented.
    rc_consensus = reverse_complement(consensus_slice)
    min_len = min(len(original_seq), len(consensus_slice))
    
    direct_matches = sum(1 for c1, c2 in zip(original_seq[:min_len], consensus_slice[:min_len]) if c1 == c2)
    rc_matches = sum(1 for c1, c2 in zip(original_seq[:min_len], rc_consensus[:min_len]) if c1 == c2)
    
    is_reverse = rc_matches > direct_matches

    # Resolve excess degeneracy if any
    rescued_seq = resolve_excess_degeneracy(iupac_slice, consensus_slice, max_multiplier)
    
    # If reverse oriented, we must reverse-complement the rescued degenerate sequence
    # to keep it in the 5'-3' primer direction.
    if is_reverse:
        rescued_seq = reverse_complement(rescued_seq)
    
    multiplier = get_degeneracy_multiplier(rescued_seq)
    possible_seqs = enumerate_sequences(rescued_seq) if multiplier > 1 else [rescued_seq]

    # Calculate thermodynamics
    thermo = evaluate_thermodynamics(possible_seqs, thermo_mode)

    # Update primer dict
    updated = dict(primer_dict)
    updated["majority_sequence"] = original_seq
    updated["sequence"] = rescued_seq
    updated["degeneracy_multiplier"] = multiplier
    updated["possible_sequences"] = possible_seqs if multiplier > 1 else []
    updated["tm"] = thermo["tm_str"]
    updated["gc"] = thermo["gc_str"]
    updated["hairpin_dg"] = f"{thermo['worst_hairpin_dg']:.2f} kcal/mol"
    updated["homodimer_dg"] = f"{thermo['worst_homodimer_dg']:.2f} kcal/mol"
    
    # Keep numerical values for calculations
    updated["_tm_min"] = thermo["tm_min"]
    updated["_tm_max"] = thermo["tm_max"]
    updated["_tm_avg"] = thermo["tm_avg"]
    updated["_gc_min"] = thermo["gc_min"]
    updated["_gc_max"] = thermo["gc_max"]
    updated["_gc_avg"] = thermo["gc_avg"]
    updated["_hairpin_dg"] = thermo["worst_hairpin_dg"]
    updated["_homodimer_dg"] = thermo["worst_homodimer_dg"]

    return updated

def main():
    parser = argparse.ArgumentParser(description="M08 Degenerate Rescue Tool")
    parser.add_argument("--manifest", required=True, help="Path to M07 primer_targets_manifest.yaml")
    parser.add_argument("--alternatives", required=True, help="Path to M08 alternatives.json")
    parser.add_argument("--summary", required=True, help="Path to M08 summary.json")
    parser.add_argument("--best-primers", required=True, help="Path to M08 best_primers.json")
    parser.add_argument("--outdir", required=True, help="Output directory to write rescued JSON files")
    parser.add_argument("--maf-threshold", type=float, default=0.10, help="MAF threshold for degenerate base selection")
    parser.add_argument("--max-multiplier", type=int, default=16, help="Maximum allowed degeneracy combinations")
    parser.add_argument("--thermo-mode", default="worst_case", choices=["worst_case", "average", "best_case"], help="Thermodynamic evaluation mode")
    parser.add_argument("--vienna-temp", type=float, default=37.0, help="Temperature for thermodynamic calculations")

    args = parser.parse_args()

    # 1. Read Manifest
    print(f"[Rescue] Reading manifest from {args.manifest}...")
    with open(args.manifest, 'r') as f:
        manifest = yaml.safe_load(f)

    # Build sequence mapping by region_id
    regions_map = {}
    for region in manifest.get("candidate_regions", []):
        regions_map[region["region_id"]] = {
            "iupac_seq": region["sequences"]["iupac_consensus_sequence"],
            "consensus_seq": region["sequences"]["consensus_sequence"]
        }

    # Discover region_id from paths or manifest
    # For subdirectories, the region_id is part of the output folder or we read it from summary
    # Let's read summary first to identify the region_id
    with open(args.summary, 'r') as f:
        summary_data = json.load(f)

    region_id = summary_data.get("metadata", {}).get("region_id") or summary_data.get("m07_metadata", {}).get("region_id")
    if not region_id:
        # Fallback: parse region_id from alternatives path
        region_id = Path(args.alternatives).parent.name
        
    print(f"[Rescue] Identified region: {region_id}")
    if region_id not in regions_map:
        print(f"[Rescue] ERROR: Region {region_id} not found in manifest candidate regions.")
        sys.exit(1)

    iupac_seq = regions_map[region_id]["iupac_seq"]
    consensus_seq = regions_map[region_id]["consensus_seq"]

    # 2. Read Alternatives
    print(f"[Rescue] Reading alternatives from {args.alternatives}...")
    with open(args.alternatives, 'r') as f:
        alternatives = json.load(f)

    # 3. Rescue Alternatives
    rescued_alternatives = []
    for idx, candidate in enumerate(alternatives):
        primers = candidate.get("primers", {})
        fwd = primers.get("forward")
        rev = primers.get("reverse")
        prb = primers.get("probe")

        # Rescue each primer type
        fwd_res = rescue_primer(fwd, iupac_seq, consensus_seq, args.max_multiplier, args.thermo_mode)
        rev_res = rescue_primer(rev, iupac_seq, consensus_seq, args.max_multiplier, args.thermo_mode)
        prb_res = rescue_primer(prb, iupac_seq, consensus_seq, args.max_multiplier, args.thermo_mode)

        # Recalculate cross heterodimerization
        fwd_seqs = fwd_res["possible_sequences"] if fwd_res["degeneracy_multiplier"] > 1 else [fwd_res["sequence"]]
        rev_seqs = rev_res["possible_sequences"] if rev_res["degeneracy_multiplier"] > 1 else [rev_res["sequence"]]
        
        worst_fwd_rev_dg = calculate_cross_heterodimer_worst(fwd_seqs, rev_seqs)

        worst_prb_fwd_dg = 0.0
        worst_prb_rev_dg = 0.0
        if prb_res:
            prb_seqs = prb_res["possible_sequences"] if prb_res["degeneracy_multiplier"] > 1 else [prb_res["sequence"]]
            worst_prb_fwd_dg = calculate_cross_heterodimer_worst(prb_seqs, fwd_seqs)
            worst_prb_rev_dg = calculate_cross_heterodimer_worst(prb_seqs, rev_seqs)

        # Update heterodimer fields
        fwd_res["heterodimer_dg"] = f"{worst_fwd_rev_dg:.2f} kcal/mol"
        fwd_res["_heterodimer_dg"] = worst_fwd_rev_dg
        rev_res["heterodimer_dg"] = f"{worst_fwd_rev_dg:.2f} kcal/mol"
        rev_res["_heterodimer_dg"] = worst_fwd_rev_dg

        if prb_res:
            prb_res["heterodimer_dg"] = f"{min(worst_prb_fwd_dg, worst_prb_rev_dg):.2f} kcal/mol"
            prb_res["_heterodimer_dg"] = min(worst_prb_fwd_dg, worst_prb_rev_dg)

        # Update candidate QC block
        qc = dict(candidate.get("qc", {}))
        qc["hairpin_dg"] = min(fwd_res["_hairpin_dg"], rev_res["_hairpin_dg"], (prb_res["_hairpin_dg"] if prb_res else 0.0))
        qc["homodimer_dg"] = min(fwd_res["_homodimer_dg"], rev_res["_homodimer_dg"], (prb_res["_homodimer_dg"] if prb_res else 0.0))
        qc["heterodimer_dg"] = worst_fwd_rev_dg
        
        cross_dg = worst_fwd_rev_dg
        if prb_res:
            cross_dg = min(worst_fwd_rev_dg, worst_prb_fwd_dg, worst_prb_rev_dg)
        qc["cross_dimer_dg"] = cross_dg

        # Tm Difference
        tm_diff = abs(fwd_res["_tm_avg"] - rev_res["_tm_avg"])
        qc["tm_diff"] = round(tm_diff, 2)

        # Check QC criteria (thermo mode check)
        qc["hairpin_ok"] = qc["hairpin_dg"] >= -2.0
        qc["homodimer_ok"] = qc["homodimer_dg"] >= -6.0
        qc["heterodimer_ok"] = qc["heterodimer_dg"] >= -7.0
        qc["tm_balance_ok"] = tm_diff <= 3.0

        # Construct new candidate object
        new_candidate = dict(candidate)
        new_candidate["primers"] = {
            "forward": fwd_res,
            "reverse": rev_res,
            "probe": prb_res
        }
        new_candidate["qc"] = qc

        # Update amplicon sequence to reflect degenerate bases if desired
        # Standard amplicon sequence can be left as-is, but Tm calculations on amplicon are updated
        new_candidate["amplicon"]["tm_forward"] = round(fwd_res["_tm_avg"], 2)
        new_candidate["amplicon"]["tm_reverse"] = round(rev_res["_tm_avg"], 2)

        rescued_alternatives.append(new_candidate)

    # Sort rescued alternatives:
    # First priority: passes QC checks
    # Second priority: lower/better final score (we can preserve primerlab's final score or re-sort based on worst-case penalties)
    # We will preserve the order, but filter out bad ones if strict, or just keep them sorted by final_score.
    # Let's keep them sorted by final_score.
    rescued_alternatives.sort(key=lambda x: x.get("final_score", 999.0))

    # Re-rank candidates
    for idx, cand in enumerate(rescued_alternatives):
        cand["final_rank"] = idx + 1

    # 4. Write Rescued Files
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[Rescue] Writing rescued alternatives to {out_dir / 'alternatives.json'}...")
    with open(out_dir / "alternatives.json", "w") as f:
        json.dump(rescued_alternatives, f, indent=2)

    # Write best_primers.json
    best_primers = {}
    if rescued_alternatives:
        best_primers = rescued_alternatives[0]["primers"]
    
    print(f"[Rescue] Writing rescued best primers to {out_dir / 'best_primers.json'}...")
    with open(out_dir / "best_primers.json", "w") as f:
        json.dump(best_primers, f, indent=2)

    # Write summary.json
    summary_data["results"] = {
        "best_pair": best_primers,
        "num_alternatives": len(rescued_alternatives)
    }
    
    print(f"[Rescue] Writing rescued summary to {out_dir / 'summary.json'}...")
    with open(out_dir / "summary.json", "w") as f:
        json.dump(summary_data, f, indent=2)

    # Write ranking details CSV
    import csv
    print(f"[Rescue] Writing rescued ranking details CSV to {out_dir / 'ranking_details.csv'}...")
    with open(out_dir / "ranking_details.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "rank", "fwd_seq", "fwd_tm", "fwd_degeneracy",
            "rev_seq", "rev_tm", "rev_degeneracy",
            "probe_seq", "probe_tm", "probe_degeneracy",
            "hairpin_dg", "homodimer_dg", "cross_dimer_dg", "tm_diff"
        ])
        for cand in rescued_alternatives:
            p = cand["primers"]
            q = cand["qc"]
            writer.writerow([
                cand["final_rank"],
                p["forward"]["sequence"], p["forward"]["tm"], p["forward"]["degeneracy_multiplier"],
                p["reverse"]["sequence"], p["reverse"]["tm"], p["reverse"]["degeneracy_multiplier"],
                p["probe"]["sequence"] if p.get("probe") else "N/A", 
                p["probe"]["tm"] if p.get("probe") else "N/A",
                p["probe"]["degeneracy_multiplier"] if p.get("probe") else 0,
                q["hairpin_dg"], q["homodimer_dg"], q["cross_dimer_dg"], q["tm_diff"]
            ])

    print("[Rescue] Degenerate rescue completed successfully.")

if __name__ == "__main__":
    main()
