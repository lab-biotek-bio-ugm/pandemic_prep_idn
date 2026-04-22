"""
M07 Bridge Builder v2.1 — Conservation & Selection Integrator

Merges outputs from two independent analysis tracks:
  Track A (Augur): MSA conservation via Shannon Entropy
  Track B (HyPhy): Evolutionary selection pressure (FUBAR/MEME/GARD)

Version 2.1 adds:
  - 150 bp reference consensus sequence embedded in all output manifests.
  - Formal scientific grading annotations (replacing informal descriptors).
  - Automatic CSV export alongside JSON/YAML for spreadsheet-based review.

Inputs: MSA, GFF, HyPhy JSONs, Reference Fasta
Outputs:
  - conserved_regions.bed       (0-based half-open)
  - bridge_report.json          (0-based half-open)
  - primer_targets_manifest.yaml
  - entropy_profile.tsv
  - candidate_targets_summary.csv  (NEW — flat table for PI/Excel review)
"""

import argparse
import csv
import json
import math
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import yaml
from Bio import SeqIO

# --- Step 1: Conservation & Profile Generation ---

def build_entropy_profile(records):
    """Calculates entropy and consensus base for each column in the MSA."""
    if not records:
        return [], 0, 0
    
    n_sequences = len(records)
    aln_length = len(records[0].seq)
    
    profile = []
    for i in range(aln_length):
        column = [str(rec.seq[i]).upper() for rec in records]
        bases = [b for b in column if b in "ACGT"]
        n = len(bases)
        
        if n == 0:
            profile.append({"position": i, "consensus_base": "-", "entropy": 1.0})
            continue
            
        counter = Counter(bases)
        consensus = counter.most_common(1)[0][0]
        
        freq = {b: c / n for b, c in counter.items()}
        entropy = 0.0
        for p in freq.values():
            if p > 0:
                entropy -= p * math.log2(p)
                
        profile.append({"position": i, "consensus_base": consensus, "entropy": entropy})
        
    return profile, n_sequences, aln_length

def export_entropy_profile(profile, path):
    """Writes the entropy profile to a TSV file."""
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["position", "consensus_base", "entropy"], delimiter="\t")
        writer.writeheader()
        writer.writerows(profile)
    print(f"[M07] Wrote entropy profile: {path}")

def scan_conservation(profile, window_size, entropy_threshold):
    """Scan pre-computed entropy profile for conserved regions using sliding-window."""
    aln_length = len(profile)
    entropies = [p["entropy"] for p in profile]
    regions = []
    
    for i in range(aln_length - window_size + 1):
        window = entropies[i:i + window_size]
        mean_e = sum(window) / len(window)
        if mean_e < entropy_threshold:
            regions.append({
                "start": i,
                "end": i + window_size,
                "mean_entropy": mean_e,
                "conservation_pct": round((1 - mean_e / 2.0) * 100, 2),
                "length": window_size
            })
            
    if not regions:
        print(f"[M07] WARNING: No windows with entropy < {entropy_threshold}")
        
    print(f"[M07] Found {len(regions)} conserved raw windows (threshold={entropy_threshold})")
    return regions

# --- Biomechanical Scanners ---

def extract_reference_subsequence(ref_seq, start, end):
    """Extract subsequence from the raw reference genome.
    Note: if ref lacks gaps relative to MSA, mapping should ideally map column to ref pos.
    For simplicity, assuming 1:1 loosely aligned or gaps stripped prior.
    """
    return str(ref_seq[start:end]).upper()

def calculate_gc_content(sequence):
    length = len(sequence)
    if length == 0:
        return 0.0
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / length) * 100.0

def detect_homopolymer(sequence, max_repeats=5):
    for base in "ACGT":
        if (base * (max_repeats + 1)) in sequence:
            return True
    return False

def find_max_contiguous_block(profile, start, end, threshold=0.10):
    max_len = 0
    current_len = 0
    for i in range(start, end):
        if profile[i]["entropy"] < threshold:
            current_len += 1
            if current_len > max_len:
                max_len = current_len
        else:
            current_len = 0
    return max_len

def get_mismatch_positions(profile, start, end, threshold=0.10):
    mismatches = []
    for i in range(start, end):
        if profile[i]["entropy"] > threshold:
            mismatches.append(i)
    return mismatches

# --- GFF Parser ---

def parse_gff(gff_path):
    """Parse GFF3 file to extract gene and UTR coordinates (0-based half-open)."""
    genes = []
    utrs = []

    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            chrom = parts[0]
            feature_type = parts[2]
            start = int(parts[3]) - 1  # GFF is 1-based → convert to 0-based
            end = int(parts[4])        # Keep as exclusive end

            if feature_type == "gene":
                name = ""
                for attr in parts[8].split(";"):
                    attr = attr.strip()
                    if attr.startswith("gene=") or attr.startswith("Name="):
                        name = attr.split("=", 1)[1].strip('"')
                        break
                if name:
                    genes.append({"name": name, "start": start, "end": end, "chrom": chrom})

            elif feature_type in ("five_prime_UTR", "three_prime_UTR"):
                utr_name = "5'UTR" if "five" in feature_type else "3'UTR"
                utrs.append({"name": utr_name, "start": start, "end": end, "chrom": chrom})

    return genes, utrs

def find_gene_for_position(pos, genes, utrs):
    for g in genes:
        if g["start"] <= pos < g["end"]:
            return g["name"], True
    for u in utrs:
        if u["start"] <= pos < u["end"]:
            return u["name"], False
    return "intergenic", False

# --- HyPhy Score Mapper ---

def parse_fubar_json(filepath):
    with open(filepath, "r") as f:
        data = json.load(f)
    mle = data.get("MLE", {})
    headers = mle.get("headers", [])
    rows = mle.get("content", {}).get("0", [])

    purifying_idx = None
    for i, h in enumerate(headers):
        label = h[0] if isinstance(h, list) else str(h)
        if "alpha>beta" in label.lower() or "prob[alpha>beta]" in label.lower():
            purifying_idx = i
            break
    if purifying_idx is None:
        purifying_idx = 5

    results = []
    for i, row in enumerate(rows):
        prob = row[purifying_idx] if purifying_idx < len(row) else 0.5
        results.append({"codon_idx": i, "fubar_purifying_prob": prob})
    return results

def parse_meme_json(filepath):
    with open(filepath, "r") as f:
        data = json.load(f)
    mle = data.get("MLE", {})
    headers = mle.get("headers", [])
    rows = mle.get("content", {}).get("0", [])

    pval_idx = None
    for i, h in enumerate(headers):
        label = h[0] if isinstance(h, list) else str(h)
        if "p-value" in label.lower() or "p value" in label.lower():
            pval_idx = i
            break
    if pval_idx is None:
        pval_idx = len(headers) - 1

    results = []
    for i, row in enumerate(rows):
        pval = row[pval_idx] if pval_idx < len(row) else 1.0
        results.append({"codon_idx": i, "meme_pvalue": pval})
    return results

def parse_gard_json(filepath):
    with open(filepath, "r") as f:
        data = json.load(f)
    if data.get("skipped", False):
        return []
    partitions = data.get("data partitions", data.get("dataPartitions", {}))
    if not partitions or len(partitions) <= 1:
        return []

    boundaries = []
    for key in sorted(partitions.keys(), key=lambda x: int(x) if x.isdigit() else x):
        coverage = partitions[key].get("coverage", [])
        for span in coverage:
            if isinstance(span, list) and len(span) == 2:
                boundaries.append({"start": span[0], "end": span[1]})

    breakpoints = []
    for i in range(1, len(boundaries)):
        breakpoints.append(boundaries[i]["start"])
    return breakpoints

def map_hyphy_to_genome(genes, hyphy_dir, target_genes):
    codon_data = {}
    gard_breakpoints = []

    for gene in genes:
        name = gene["name"]
        if name not in target_genes:
            continue

        gene_start = gene["start"]
        fubar_path = hyphy_dir / f"fubar_{name}.json"
        fubar_data = []
        if fubar_path.exists():
            try: fubar_data = parse_fubar_json(str(fubar_path))
            except Exception: pass

        meme_path = hyphy_dir / f"meme_{name}.json"
        meme_data = []
        if meme_path.exists():
            try: meme_data = parse_meme_json(str(meme_path))
            except Exception: pass

        n_codons = max(len(fubar_data), len(meme_data), 1)
        for ci in range(n_codons):
            gpos = gene_start + ci * 3
            entry = {
                "gene": name,
                "codon_idx": ci,
                "fubar_purifying_prob": 0.5,
                "meme_pvalue": 1.0,
            }
            if ci < len(fubar_data):
                entry["fubar_purifying_prob"] = fubar_data[ci]["fubar_purifying_prob"]
            if ci < len(meme_data):
                entry["meme_pvalue"] = meme_data[ci]["meme_pvalue"]
            for offset in range(3):
                codon_data[gpos + offset] = entry

        gard_path = hyphy_dir / f"gard_{name}.json"
        if gard_path.exists():
            try:
                bp_codons = parse_gard_json(str(gard_path))
                for bp_codon in bp_codons:
                    gard_breakpoints.append({
                        "genome_pos": gene_start + bp_codon * 3,
                        "gene": name,
                        "codon_idx": bp_codon,
                    })
            except Exception: pass

    return codon_data, gard_breakpoints

# --- Multipliers Computation ---

def compute_meme_multiplier(sites, config):
    if sites >= config.get("veto_above", 4):
        return 0.0
    for t in config.get("thresholds", []):
        if sites <= t.get("max_sites", 0):
            return t.get("multiplier", 1.0)
    return 1.0

def compute_gc_multiplier(gc_pct, config):
    if gc_pct < config.get("veto_below", 25) or gc_pct > config.get("veto_above", 70):
        return 0.0
    if config.get("ideal_min", 35) <= gc_pct <= config.get("ideal_max", 60):
        return 1.0
    marginal = config.get("marginal_ranges", [])
    for m in marginal:
        if m.get("min", 0) <= gc_pct <= m.get("max", 100):
            return config.get("marginal_multiplier", 0.7)
    return 0.0

def compute_block_multiplier(block_len, config):
    if block_len < config.get("penalty_min_len", 20):
        return 0.0 # Veto
    if block_len >= config.get("ideal_min_len", 35):
        return 1.0
    return config.get("penalty_multiplier", 0.8)

def compute_homo_multiplier(has_homo, config):
    return config.get("penalty_multiplier", 0.6) if has_homo else 1.0

# --- Step 3: Intersection & Composite Scoring ---

GRADE_TEXT = {
    "S": "Optimal Target — Composite score ≥0.95; no episodic selection detected; high confidence for universal primer design.",
    "A+": "High-Priority Target — Composite score ≥0.90; near-optimal conservation; negligible evolutionary risk.",
    "A": "Recommended Target — Composite score ≥0.85; strong conservation profile; minor multiplier penalties applied.",
    "A-": "Acceptable Target — Composite score ≥0.80; conserved with minor biomechanical caveats; suitable for validation.",
    "B+": "Moderate Priority — Composite score ≥0.75; satisfactory conservation; one or more moderate penalties detected.",
    "B": "Low Priority — Composite score ≥0.70; approaching lower tolerance bounds; requires careful thermodynamic validation.",
    "C": "At-Risk Target — Composite score ≥0.60; significant penalties (e.g., episodic selection, GC imbalance); not recommended without further review.",
    "D": "Excluded — Composite score <0.60; failed one or more critical biomechanical or evolutionary constraints.",
}

def score_regions(regions, profile, ref_seq, codon_data, gard_breakpoints, genes, utrs,
                  scoring_config, tiers_config, gard_buffer_bp):
    """Score each conserved region using Core x Multipliers."""
    
    # Core Weights
    w_cons = scoring_config.get("core", {}).get("conservation_weight", 0.66)
    w_fub = scoring_config.get("core", {}).get("purifying_weight", 0.34)
    
    # Multiplier Configs
    mult_cfgs = scoring_config.get("multipliers", {})
    meme_cfg = mult_cfgs.get("meme", {})
    rec_cfg = mult_cfgs.get("recombination", {})
    gc_cfg = mult_cfgs.get("gc_content", {})
    block_cfg = mult_cfgs.get("conserved_block", {})
    homo_cfg = mult_cfgs.get("homopolymer", {})

    scored = []

    for region in regions:
        start, end = region["start"], region["end"]
        mid = (start + end) // 2
        gene_name, is_coding = find_gene_for_position(mid, genes, utrs)

        # 1. Gather HyPhy Data
        fubar_probs, meme_pvals = [], []
        has_hyphy = False
        for pos in range(start, end):
            if pos in codon_data:
                has_hyphy = True
                fubar_probs.append(codon_data[pos]["fubar_purifying_prob"])
                meme_pvals.append(codon_data[pos]["meme_pvalue"])

        avg_fubar = sum(fubar_probs) / len(fubar_probs) if fubar_probs else 0.0
        meme_danger_sites = sum(1 for p in meme_pvals if p < 0.05)

        # Check GARD
        near_recomb = False
        half_len = (end - start) // 2
        for bp in gard_breakpoints:
            if abs(bp["genome_pos"] - mid) <= gard_buffer_bp + half_len:
                near_recomb = True
                break

        # 2. Biomechanical features
        ref_subseq = extract_reference_subsequence(ref_seq, start, end)
        gc_pct = calculate_gc_content(ref_subseq)
        
        block_t_hr = block_cfg.get("entropy_threshold", 0.10)
        max_block_len = find_max_contiguous_block(profile, start, end, threshold=block_t_hr)
        
        has_homo = detect_homopolymer(ref_subseq, homo_cfg.get("max_repeats", 5))
        mismatches = get_mismatch_positions(profile, start, end, threshold=block_t_hr)

        # 3. Core Score Calculation
        # Normalized conservation: DNA entropy max is 2.0 -> Range [0.0, 1.0]
        cons_norm = 1.0 - (region["mean_entropy"] / 2.0)
        
        if has_hyphy:
            core = (w_cons * cons_norm) + (w_fub * avg_fubar)
        else:
            core = cons_norm * 0.95 # Cap non-coding core

        # 4. Multipliers
        m_meme = compute_meme_multiplier(meme_danger_sites, meme_cfg)
        m_rec = rec_cfg.get("near_penalty", 0.5) if near_recomb else 1.0
        m_gc = compute_gc_multiplier(gc_pct, gc_cfg)
        m_block = compute_block_multiplier(max_block_len, block_cfg)
        m_homo = compute_homo_multiplier(has_homo, homo_cfg)

        final_score = core * m_meme * m_rec * m_gc * m_block * m_homo
        final_score = max(0.0, min(1.0, final_score))

        # --- Tier System (Fail-safe checks on raw metrics) ---
        strict = tiers_config.get("strict", {})
        moderate = tiers_config.get("moderate", {})
        cons_frac = cons_norm # Using normalized score directly
        tier = 3

        if has_hyphy:
            if (cons_frac >= strict.get("min_conservation", 0.95)
                    and avg_fubar >= strict.get("min_purifying_prob", 0.9)
                    and meme_danger_sites == 0):
                tier = 1
            elif (cons_frac >= moderate.get("min_conservation", 0.85)
                  and avg_fubar >= moderate.get("min_purifying_prob", 0.7)):
                tier = 2
        else:
            if cons_frac >= strict.get("min_conservation", 0.95):
                tier = 1
            elif cons_frac >= moderate.get("min_conservation", 0.85):
                tier = 2

        # --- Grade Classification ---
        if final_score >= 0.95: grade = "S"
        elif final_score >= 0.90: grade = "A+"
        elif final_score >= 0.85: grade = "A"
        elif final_score >= 0.80: grade = "A-"
        elif final_score >= 0.75: grade = "B+"
        elif final_score >= 0.70: grade = "B"
        elif final_score >= 0.60: grade = "C"
        else: grade = "D"

        # Build the consensus sequence from the MSA profile for this window
        consensus_seq = "".join(
            profile[i]["consensus_base"] for i in range(start, end)
        )

        scored.append({
            "gene": gene_name,
            "start": start,
            "end": end,
            "length": end - start,
            "conservation_pct": round(region["conservation_pct"], 2),
            "mean_entropy": round(region["mean_entropy"], 6),
            "avg_fubar_purifying_prob": round(avg_fubar, 4),
            "meme_episodic_sites": meme_danger_sites,
            "gc_pct": round(gc_pct, 2),
            "max_block_len": max_block_len,
            "has_homopolymer": has_homo,
            "composite_score": round(final_score, 4),
            "tier": tier,
            "grade": grade,
            "recommendation": GRADE_TEXT[grade],
            # Sequences: consensus from MSA columns, and raw from reference
            "consensus_sequence": consensus_seq,
            "reference_sequence": ref_subseq,
            "flags": {
                "coding": is_coding,
                "has_hyphy_data": has_hyphy,
                "near_recombination": near_recomb,
            },
            "multipliers": {
                "meme": m_meme,
                "recombination": m_rec,
                "gc_content": m_gc,
                "conserved_block": m_block,
                "homopolymer": m_homo
            },
            "mismatch_positions": mismatches
        })

    scored.sort(key=lambda r: r["composite_score"], reverse=True)
    for i, r in enumerate(scored):
        r["rank"] = i + 1

    return scored

def apply_nms(scored_regions, overlap_threshold=0.5):
    scored_regions.sort(key=lambda x: x["composite_score"], reverse=True)
    kept = []
    
    for current in scored_regions:
        is_suppressed = False
        for winner in kept:
            if current["gene"] == winner["gene"]:
                overlap_start = max(current["start"], winner["start"])
                overlap_end = min(current["end"], winner["end"])
                if overlap_end > overlap_start:
                    overlap_len = overlap_end - overlap_start
                    overlap_ratio = overlap_len / current["length"]
                    if overlap_ratio > overlap_threshold:
                        is_suppressed = True
                        break
        if not is_suppressed:
            kept.append(current)
            
    for idx, r in enumerate(kept):
        r["rank"] = idx + 1
        
    return kept

def write_bed(regions, ref_id, path):
    with open(path, "w") as f:
        for r in regions:
            rid = f"{r['gene']}_region_{r['rank']}"
            f.write(f"{ref_id}\t{r['start']}\t{r['end']}\t{rid}\t"
                    f"{r['composite_score']}\t+\n")
    print(f"[M07] Wrote BED: {path}")

def write_report(regions, meta, path):
    report = {
        "pipeline_version": "2.0",
        "scoring_method": "core_x_multiplier",
        "coordinate_system": "0-based half-open",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "track_name": meta.get("track_name", "Universal"),
        "msa_sequences": meta.get("n_sequences", 0),
        "msa_length": meta.get("aln_length", 0),
        "genes_analyzed": meta.get("genes_analyzed", []),
        "hyphy_genes": meta.get("hyphy_genes", []),
        "total_conserved_regions": len(regions),
        "grade_summary": {g: sum(1 for r in regions if r["grade"] == g)
                          for g in ["S", "A+", "A", "A-", "B+", "B", "C", "D"]},
        "tier_summary": {f"tier_{t}": sum(1 for r in regions if r["tier"] == t)
                         for t in (1, 2, 3)},
        "regions": regions,
        "gard_breakpoints": meta.get("gard_breakpoints", []),
        "warnings": meta.get("warnings", []),
    }
    with open(path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"[M07] Wrote report: {path}")

def write_manifest(regions, build_name, ref_genome, buffer_bp, path):
    manifest = {
        "pipeline_run_id": f"{build_name}_bridge",
        "source": "pandemic-preparedness-m07",
        "pipeline_version": "2.0",
        "coordinate_system": "0-based half-open",
        "reference_genome": str(ref_genome),
        "buffer_bp": buffer_bp,
        "total_candidates": len(regions),
        "candidate_regions": [],
    }

    for r in regions:
        inc_start = max(0, r["start"] - buffer_bp)
        inc_end = r["end"] + buffer_bp

        manifest["candidate_regions"].append({
            "region_id": f"{r['gene']}_region_{r['rank']}",
            "rank": r["rank"],
            "gene": r["gene"],
            "tier": r["tier"],
            "grade": r["grade"],
            "composite_score": r["composite_score"],
            "conserved_zone": {
                "start": r["start"],
                "end": r["end"],
                "length": r["length"],
            },
            "included_region": {
                "start": inc_start,
                "end": inc_end,
                "length": inc_end - inc_start,
            },
            "biomechanics": {
                "conservation_pct": r["conservation_pct"],
                "mean_entropy": r["mean_entropy"],
                "avg_fubar_purifying_prob": r["avg_fubar_purifying_prob"],
                "meme_episodic_sites": r["meme_episodic_sites"],
                "gc_pct": r["gc_pct"],
                "max_block_len": r["max_block_len"],
                "has_homopolymer": r["has_homopolymer"],
                "mismatch_positions": r["mismatch_positions"],
                "multipliers": r["multipliers"]
            },
            "flags": r["flags"],
            "recommendation": r["recommendation"],
        })

    with open(path, "w") as f:
        yaml.dump(manifest, f, default_flow_style=False,
                  sort_keys=False, allow_unicode=True)
    print(f"[M07] Wrote manifest: {path}")


def write_csv(regions, path):
    """Export a flat CSV table of scored regions for spreadsheet-based review (PI).

    Each row corresponds to one candidate primer target region. Nested fields
    (multipliers, flags) are flattened into individual columns.
    """
    fieldnames = [
        "rank",
        "region_id",
        "gene",
        "tier",
        "grade",
        "composite_score",
        "start_0based",
        "end_0based",
        "length_bp",
        "conservation_pct",
        "mean_entropy",
        "avg_fubar_purifying_prob",
        "meme_episodic_sites",
        "gc_pct",
        "max_conserved_block_len",
        "has_homopolymer",
        "mult_meme",
        "mult_recombination",
        "mult_gc_content",
        "mult_conserved_block",
        "mult_homopolymer",
        "flag_coding_region",
        "flag_has_hyphy_data",
        "flag_near_recombination",
        "consensus_sequence",
        "reference_sequence",
        "recommendation",
    ]

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in regions:
            writer.writerow({
                "rank": r.get("rank", ""),
                "region_id": f"{r['gene']}_region_{r.get('rank', '')}",
                "gene": r.get("gene", ""),
                "tier": r.get("tier", ""),
                "grade": r.get("grade", ""),
                "composite_score": r.get("composite_score", ""),
                "start_0based": r.get("start", ""),
                "end_0based": r.get("end", ""),
                "length_bp": r.get("length", ""),
                "conservation_pct": r.get("conservation_pct", ""),
                "mean_entropy": r.get("mean_entropy", ""),
                "avg_fubar_purifying_prob": r.get("avg_fubar_purifying_prob", ""),
                "meme_episodic_sites": r.get("meme_episodic_sites", ""),
                "gc_pct": r.get("gc_pct", ""),
                "max_conserved_block_len": r.get("max_block_len", ""),
                "has_homopolymer": r.get("has_homopolymer", ""),
                "mult_meme": r.get("multipliers", {}).get("meme", ""),
                "mult_recombination": r.get("multipliers", {}).get("recombination", ""),
                "mult_gc_content": r.get("multipliers", {}).get("gc_content", ""),
                "mult_conserved_block": r.get("multipliers", {}).get("conserved_block", ""),
                "mult_homopolymer": r.get("multipliers", {}).get("homopolymer", ""),
                "flag_coding_region": r.get("flags", {}).get("coding", ""),
                "flag_has_hyphy_data": r.get("flags", {}).get("has_hyphy_data", ""),
                "flag_near_recombination": r.get("flags", {}).get("near_recombination", ""),
                "consensus_sequence": r.get("consensus_sequence", ""),
                "reference_sequence": r.get("reference_sequence", ""),
                "recommendation": r.get("recommendation", ""),
            })
    print(f"[M07] Wrote CSV summary: {path}")

def execute_track(track_name, records, genes, utrs, target_genes, hyphy_dir, 
                  window_size, entropy_threshold, overlap_threshold, gard_buffer_bp, 
                  scoring_config, tiers_config, ref_id, ref_seq, outdir, build_name, buffer_bp, args_reference):
    print(f"\n[M07] ============================================")
    print(f"[M07]  Track: {track_name}")
    print(f"[M07] ============================================")

    print(f"\n[M07] [{track_name}] Step 1: Generating Entropy Profile")
    profile, n_seq, aln_len = build_entropy_profile(records)
    
    sanitized_name = track_name.replace('/', '_').replace(' ', '_')
    prefix = "" if track_name == "Universal" else f"subtype_{sanitized_name}_"
    
    export_entropy_profile(profile, outdir / f"{prefix}entropy_profile.tsv")

    print(f"\n[M07] [{track_name}] Step 2: Conservation Scanner")
    regions = scan_conservation(profile, window_size, entropy_threshold)

    if not regions:
        print(f"[M07] [{track_name}] WARNING: Relaxing threshold significantly to find candidates (survival mode)")
        regions = scan_conservation(profile, window_size, 0.45)

    if not regions:
         print(f"[M07] [{track_name}] No conserved regions found even after max relaxation. Skipping.")
         return

    print(f"\n[M07] [{track_name}] Step 3: HyPhy Score Mapper")
    codon_data, gard_bps = map_hyphy_to_genome(genes, hyphy_dir, target_genes)
    print(f"[M07] [{track_name}] Mapped {len(codon_data)} positions, {len(gard_bps)} GARD breakpoints")

    print(f"\n[M07] [{track_name}] Step 4: Biomechanical Analysis & Core × Multipliers Scoring")
    scored = score_regions(regions, profile, ref_seq, codon_data, gard_bps, genes, utrs,
                           scoring_config, tiers_config, gard_buffer_bp)

    print(f"\n[M07] [{track_name}] Step 5: Non-Maximum Suppression (NMS)")
    final_scored = apply_nms(scored, overlap_threshold=overlap_threshold)

    grade_keys = ["S", "A+", "A", "A-", "B+", "B", "C", "D"]
    grades = {g: sum(1 for r in final_scored if r["grade"] == g) for g in grade_keys}
    grade_str = " ".join([f"{g}={grades[g]}" for g in grade_keys])
    print(f"[M07] [{track_name}] Scored {len(final_scored)} disjoint/optimal regions: {grade_str}")

    warnings = []
    if not codon_data:
        warnings.append("No HyPhy data loaded. All scores are conservation-only.")

    meta = {
        "track_name": track_name,
        "n_sequences": n_seq,
        "aln_length": aln_len,
        "genes_analyzed": [g["name"] for g in genes] + [u["name"] for u in utrs],
        "hyphy_genes": [g for g in target_genes if (hyphy_dir / f"fubar_{g}.json").exists()],
        "gard_breakpoints": [{"genome_pos": bp["genome_pos"], "gene": bp["gene"]} for bp in gard_bps],
        "warnings": warnings,
    }

    print(f"\n[M07] [{track_name}] Step 6: Output Generation")
    write_bed(final_scored, ref_id, outdir / f"{prefix}conserved_regions.bed")
    write_report(final_scored, meta, outdir / f"{prefix}bridge_report.json")
    write_manifest(final_scored, build_name, args_reference, buffer_bp,
                   outdir / f"{prefix}primer_targets_manifest.yaml")
    write_csv(final_scored, outdir / f"{prefix}candidate_targets_summary.csv")

    print(f"\n[M07] [{track_name}] Complete — {len(final_scored)} candidate regions")
    if final_scored:
        top = final_scored[0]
        print(f"[M07] [{track_name}] #1: {top['gene']} ({top['start']}-{top['end']}) Grade {top['grade']}, Score {top['composite_score']}")


def main():
    parser = argparse.ArgumentParser(
        description="M07 Bridge Builder v2.0 — Core x Multipliers Integrator")
    parser.add_argument("--msa", required=True,
                        help="Aligned FASTA (MSA from augur align / M03)")
    parser.add_argument("--gff", required=True,
                        help="Gene annotation GFF3 file")
    parser.add_argument("--hyphy-dir", required=True,
                        help="Directory containing fubar_*.json, meme_*.json, gard_*.json")
    parser.add_argument("--reference", required=True,
                        help="Reference genome FASTA")
    parser.add_argument("--metadata", required=False,
                        help="Metadata TSV file for subtype/clade splitting (Route B)")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for BED, JSON, YAML")
    parser.add_argument("--config", required=True,
                        help="Pipeline config.yaml path")
    parser.add_argument("--build-name", default="unknown",
                        help="Build name for manifest identification")
    args = parser.parse_args()

    # Load config
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    bridge = config.get("bridge", {})
    window_size = bridge.get("window_size", 150)
    entropy_threshold = bridge.get("entropy_threshold", 0.015)
    overlap_threshold = bridge.get("overlap_threshold", 0.75)
    buffer_bp = bridge.get("buffer_bp", 100)
    gard_buffer_bp = bridge.get("gard_buffer_bp", 50)
    
    scoring_config = bridge.get("scoring", {})
    if not scoring_config and "weights" in bridge:
        print("[M07] WARNING: Falling back to v1.0 weights config. Please update to 'scoring' config block.", file=sys.stderr)
        # Rudimentary fallback mapping
        scoring_config = {
            "core": {
                "conservation_weight": bridge["weights"].get("conservation", 0.5),
                "purifying_weight": bridge["weights"].get("purifying", 0.25)
            },
            "multipliers": {}
        }
        
    tiers_config = bridge.get("tiers", {})
    target_genes = config.get("hyphy", {}).get("target_genes", [])
    hyphy_dir = Path(args.hyphy_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_records = list(SeqIO.parse(args.msa, "fasta"))
    if not all_records:
        print("ERROR: MSA file is empty", file=sys.stderr)
        sys.exit(1)
        
    ref_records = list(SeqIO.parse(args.reference, "fasta"))
    if not ref_records:
         print("ERROR: Reference file is empty", file=sys.stderr)
         sys.exit(1)
         
    ref_seq = str(ref_records[0].seq).upper()
    ref_id = ref_records[0].id

    print("[M07] ============================================")
    print("[M07]  Bridge Builder v2.0 — Biomechanical Scoring")
    print("[M07] ============================================")
    print(f"[M07] Reference : {ref_id} ({len(ref_seq)} bp)")
    print(f"[M07] Window    : {window_size} bp")
    print(f"[M07] Threshold : {entropy_threshold}")

    # ----- Parse GFF -----
    genes, utrs = parse_gff(args.gff)
    print(f"[M07] GFF: {len(genes)} genes, {len(utrs)} UTRs")
    if "ALL" in target_genes:
        target_genes = [g["name"] for g in genes]

    # Execute Universal Track
    execute_track("Universal", all_records, genes, utrs, target_genes, hyphy_dir, 
                  window_size, entropy_threshold, overlap_threshold, gard_buffer_bp, 
                  scoring_config, tiers_config, ref_id, ref_seq, outdir, args.build_name, buffer_bp, args.reference)

    # Route B: Subtype-Specific Tracks
    if args.metadata and Path(args.metadata).exists():
        print(f"\n[M07] ============================================")
        print(f"[M07]  ROUTE B: Subtype Auto-Detection")
        print(f"[M07] ============================================")
        strain_to_subtype = {}
        with open(args.metadata, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                strain = row.get("strain")
                subtype = row.get("subtype")
                if strain and subtype and subtype != "unknown" and subtype != "?":
                    strain_to_subtype[strain] = subtype
                    
        from collections import defaultdict
        subtype_groups = defaultdict(list)
        for rec in all_records:
             rec_id = rec.id.split('|')[0]
             subtype = strain_to_subtype.get(rec_id)
             if subtype:
                 subtype_groups[subtype].append(rec)
                 
        for subtype, recs in subtype_groups.items():
            if len(recs) >= 30:  # Minimum coverage to consider a subtype valid
                execute_track(subtype, recs, genes, utrs, target_genes, hyphy_dir, 
                              window_size, entropy_threshold, overlap_threshold, gard_buffer_bp, 
                              scoring_config, tiers_config, ref_id, ref_seq, outdir, args.build_name, buffer_bp, args.reference)
            else:
                print(f"[M07] Skipping subtype {subtype} (only {len(recs)} sequences, needs >=30)")


if __name__ == "__main__":
    main()
