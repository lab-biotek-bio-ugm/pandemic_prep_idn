"""
M07 Bridge Builder — Conservation & Selection Pressure Integrator

Merges outputs from two independent analysis tracks:
  Track A (Augur): MSA conservation via Shannon Entropy
  Track B (HyPhy): Evolutionary selection pressure (FUBAR/MEME/GARD)

Produces ranked "Gold Zones" — genomic regions that are both empirically
conserved AND under active purifying selection — ideal for primer design.

Outputs:
  1. conserved_regions.bed          — BED format for genome browser (IGV)
  2. bridge_report.json             — Comprehensive analysis report
  3. primer_targets_manifest.yaml   — Tool-agnostic primer target coordinates

Scoring: composite = w1*(1-entropy) + w2*fubar_purifying + w3*(1-meme_danger) - w4*recomb
Grade:   A (>=0.90), B (>=0.75), C (>=0.60), D (<0.60)
Tier:    1 (strict), 2 (moderate), 3 (survival) — fail-safe protocol
Flags:   coding, has_hyphy_data, near_recombination
"""

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml
from Bio import SeqIO


# Step 1: Conservation Scanner

def calculate_entropy(column):
    """Calculate Shannon entropy for a single MSA column.
    Returns 0.0 for perfectly conserved, up to 2.0 for maximum diversity (DNA).
    """
    bases = [b for b in column if b in "ACGT"]
    n = len(bases)
    if n == 0:
        return 1.0  # Full uncertainty if all gaps/N

    freq = {}
    for b in bases:
        freq[b] = freq.get(b, 0) + 1

    entropy = 0.0
    for count in freq.values():
        p = count / n
        if p > 0:
            entropy -= p * math.log2(p)

    return entropy


def scan_conservation(msa_path, window_size, entropy_threshold, min_region_length):
    """Scan MSA for conserved regions using sliding-window entropy."""
    records = list(SeqIO.parse(str(msa_path), "fasta"))
    if not records:
        print("ERROR: MSA file is empty", file=sys.stderr)
        sys.exit(1)

    n_sequences = len(records)
    aln_length = len(records[0].seq)
    print(f"[M07] MSA: {n_sequences} sequences, {aln_length} positions")

    # Per-position entropy
    entropies = []
    for i in range(aln_length):
        column = [str(rec.seq[i]).upper() for rec in records]
        entropies.append(calculate_entropy(column))

    # Sliding window average
    conserved_windows = []
    for i in range(aln_length - window_size + 1):
        window = entropies[i:i + window_size]
        mean_e = sum(window) / len(window)
        if mean_e < entropy_threshold:
            conserved_windows.append({"start": i, "end": i + window_size, "mean_entropy": mean_e})

    if not conserved_windows:
        print(f"[M07] WARNING: No windows with entropy < {entropy_threshold}")
        return [], n_sequences, aln_length

    # Merge overlapping/adjacent windows
    regions = []
    current = conserved_windows[0].copy()

    for w in conserved_windows[1:]:
        if w["start"] <= current["end"]:
            current["end"] = max(current["end"], w["end"])
            current["mean_entropy"] = min(current["mean_entropy"], w["mean_entropy"])
        else:
            if current["end"] - current["start"] >= min_region_length:
                regions.append(current)
            current = w.copy()

    if current["end"] - current["start"] >= min_region_length:
        regions.append(current)

    # Recalculate mean entropy for merged regions
    for region in regions:
        region_e = entropies[region["start"]:region["end"]]
        region["mean_entropy"] = sum(region_e) / len(region_e)
        # Max DNA entropy = log2(4) = 2.0 → conservation_pct relative to that
        region["conservation_pct"] = round((1 - region["mean_entropy"] / 2.0) * 100, 2)
        region["length"] = region["end"] - region["start"]

    print(f"[M07] Found {len(regions)} conserved regions (threshold={entropy_threshold})")
    return regions, n_sequences, aln_length


# GFF Parser

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
    """Determine which gene/UTR a genome position falls into.
    Returns (name, is_coding).
    """
    for g in genes:
        if g["start"] <= pos < g["end"]:
            return g["name"], True
    for u in utrs:
        if u["start"] <= pos < u["end"]:
            return u["name"], False
    return "intergenic", False


# Step 2: HyPhy Score Mapper

def parse_fubar_json(filepath):
    """Parse HyPhy FUBAR JSON → per-codon purifying probability."""
    with open(filepath, "r") as f:
        data = json.load(f)

    mle = data.get("MLE", {})
    headers = mle.get("headers", [])
    rows = mle.get("content", {}).get("0", [])

    # Find Prob[alpha>beta] (purifying selection) column index
    purifying_idx = None
    for i, h in enumerate(headers):
        label = h[0] if isinstance(h, list) else str(h)
        if "alpha>beta" in label.lower() or "prob[alpha>beta]" in label.lower():
            purifying_idx = i
            break
    if purifying_idx is None:
        purifying_idx = 5  # FUBAR default layout fallback

    results = []
    for i, row in enumerate(rows):
        prob = row[purifying_idx] if purifying_idx < len(row) else 0.5
        results.append({"codon_idx": i, "fubar_purifying_prob": prob})

    return results


def parse_meme_json(filepath):
    """Parse HyPhy MEME JSON → per-codon p-value for episodic selection."""
    with open(filepath, "r") as f:
        data = json.load(f)

    mle = data.get("MLE", {})
    headers = mle.get("headers", [])
    rows = mle.get("content", {}).get("0", [])

    # Find p-value column
    pval_idx = None
    for i, h in enumerate(headers):
        label = h[0] if isinstance(h, list) else str(h)
        if "p-value" in label.lower() or "p value" in label.lower():
            pval_idx = i
            break
    if pval_idx is None:
        pval_idx = len(headers) - 1  # Last column fallback

    results = []
    for i, row in enumerate(rows):
        pval = row[pval_idx] if pval_idx < len(row) else 1.0
        results.append({"codon_idx": i, "meme_pvalue": pval})

    return results


def parse_gard_json(filepath):
    """Parse HyPhy GARD JSON → breakpoint codon positions.
    Breakpoints are boundaries between data partitions.
    """
    with open(filepath, "r") as f:
        data = json.load(f)

    # Check if GARD was skipped by the safe wrapper
    if data.get("skipped", False):
        print(f"[M07]   GARD  was skipped: {data.get('analysis', {}).get('skip_reason', 'unknown')}",
              file=sys.stderr)
        return []

    # GARD stores results in "data partitions" (or variants)
    partitions = data.get("data partitions",
                          data.get("dataPartitions", {}))

    if not partitions or len(partitions) <= 1:
        return []  # Single partition = no recombination detected

    # Collect partition boundaries sorted by key
    boundaries = []
    for key in sorted(partitions.keys(), key=lambda x: int(x) if x.isdigit() else x):
        coverage = partitions[key].get("coverage", [])
        for span in coverage:
            if isinstance(span, list) and len(span) == 2:
                boundaries.append({"start": span[0], "end": span[1]})

    # Breakpoints are at transitions between consecutive partitions
    breakpoints = []
    for i in range(1, len(boundaries)):
        breakpoints.append(boundaries[i]["start"])

    return breakpoints


def map_hyphy_to_genome(genes, hyphy_dir, target_genes):
    """Map per-gene HyPhy codon data back to genome coordinates."""
    codon_data = {}      # genome_pos → {fubar_purifying_prob, meme_pvalue}
    gard_breakpoints = []

    for gene in genes:
        name = gene["name"]
        if name not in target_genes:
            continue

        gene_start = gene["start"]  # 0-based

        # --- FUBAR ---
        fubar_path = hyphy_dir / f"fubar_{name}.json"
        fubar_data = []
        if fubar_path.exists():
            try:
                fubar_data = parse_fubar_json(str(fubar_path))
                print(f"[M07]   FUBAR {name}: {len(fubar_data)} codons")
            except Exception as e:
                print(f"[M07]   WARNING: FUBAR parse failed for {name}: {e}",
                      file=sys.stderr)

        # --- MEME ---
        meme_path = hyphy_dir / f"meme_{name}.json"
        meme_data = []
        if meme_path.exists():
            try:
                meme_data = parse_meme_json(str(meme_path))
                print(f"[M07]   MEME  {name}: {len(meme_data)} codons")
            except Exception as e:
                print(f"[M07]   WARNING: MEME parse failed for {name}: {e}",
                      file=sys.stderr)

        # --- Map codons → genome positions ---
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
            # Spread codon data across its 3 base positions
            for offset in range(3):
                codon_data[gpos + offset] = entry

        # --- GARD ---
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
                print(f"[M07]   GARD  {name}: {len(bp_codons)} breakpoints")
            except Exception as e:
                print(f"[M07]   WARNING: GARD parse failed for {name}: {e}",
                      file=sys.stderr)

    return codon_data, gard_breakpoints


# Step 3: Intersection & Composite Scoring

GRADE_TEXT = {
    "A": "EXCELLENT — Ideal for universal primer design",
    "B": "GOOD — Strong candidate, minor concerns",
    "C": "ACCEPTABLE — Usable with caveats",
    "D": "CAUTION — Consider alternatives",
}


def score_regions(regions, codon_data, gard_breakpoints, genes, utrs,
                  weights, tiers_config, gard_buffer_bp):
    """Score each conserved region with the composite formula and assign grade/tier."""
    w1 = weights.get("conservation", 0.35)
    w2 = weights.get("purifying", 0.35)
    w3 = weights.get("episodic_penalty", 0.15)
    w4 = weights.get("recombination", 0.15)

    scored = []

    for region in regions:
        start, end = region["start"], region["end"]
        mid = (start + end) // 2

        gene_name, is_coding = find_gene_for_position(mid, genes, utrs)

        # Gather HyPhy data overlapping this region
        fubar_probs, meme_pvals = [], []
        has_hyphy = False

        for pos in range(start, end):
            if pos in codon_data:
                has_hyphy = True
                fubar_probs.append(codon_data[pos]["fubar_purifying_prob"])
                meme_pvals.append(codon_data[pos]["meme_pvalue"])

        avg_fubar = sum(fubar_probs) / len(fubar_probs) if fubar_probs else 0.0
        meme_danger_sites = sum(1 for p in meme_pvals if p < 0.05)
        meme_danger_pct = meme_danger_sites / len(meme_pvals) if meme_pvals else 0.0

        # Check GARD proximity
        near_recomb = False
        recomb_penalty = 0.0
        half_len = (end - start) // 2
        for bp in gard_breakpoints:
            if abs(bp["genome_pos"] - mid) <= gard_buffer_bp + half_len:
                near_recomb = True
                recomb_penalty = 1.0
                break

        # --- Composite score ---
        cons_score = 1.0 - region["mean_entropy"]  # Higher is better

        if has_hyphy:
            composite = (w1 * cons_score
                         + w2 * avg_fubar
                         + w3 * (1.0 - meme_danger_pct)
                         - w4 * recomb_penalty)
        else:
            # No HyPhy data (non-coding or gene not in target_genes):
            # Use conservation-only score, capped so it can't beat full-data A
            composite = cons_score * 0.95

        composite = max(0.0, min(1.0, composite))

        # --- Tier (fail-safe protocol) ---
        strict = tiers_config.get("strict", {})
        moderate = tiers_config.get("moderate", {})
        cons_frac = region["conservation_pct"] / 100.0
        tier = 3

        if has_hyphy:
            if (cons_frac >= strict.get("min_conservation", 0.99)
                    and avg_fubar >= strict.get("min_purifying_prob", 0.9)
                    and meme_danger_sites == 0):
                tier = 1
            elif (cons_frac >= moderate.get("min_conservation", 0.95)
                  and avg_fubar >= moderate.get("min_purifying_prob", 0.7)):
                tier = 2
        else:
            if cons_frac >= strict.get("min_conservation", 0.99):
                tier = 1
            elif cons_frac >= moderate.get("min_conservation", 0.95):
                tier = 2

        # --- Grade (quality classification) ---
        if composite >= 0.90:
            grade = "A"
        elif composite >= 0.75:
            grade = "B"
        elif composite >= 0.60:
            grade = "C"
        else:
            grade = "D"

        scored.append({
            "gene": gene_name,
            "start": start,
            "end": end,
            "length": end - start,
            "conservation_pct": round(region["conservation_pct"], 2),
            "mean_entropy": round(region["mean_entropy"], 6),
            "avg_fubar_purifying_prob": round(avg_fubar, 4),
            "meme_episodic_sites": meme_danger_sites,
            "meme_danger_pct": round(meme_danger_pct, 4),
            "composite_score": round(composite, 4),
            "tier": tier,
            "grade": grade,
            "recommendation": GRADE_TEXT[grade],
            "flags": {
                "coding": is_coding,
                "has_hyphy_data": has_hyphy,
                "near_recombination": near_recomb,
            },
        })

    scored.sort(key=lambda r: r["composite_score"], reverse=True)
    for i, r in enumerate(scored):
        r["rank"] = i + 1

    return scored


# Step 4: Output Builder

def write_bed(regions, ref_id, path):
    """Write BED file for genome browser visualization."""
    with open(path, "w") as f:
        for r in regions:
            rid = f"{r['gene']}_region_{r['rank']}"
            f.write(f"{ref_id}\t{r['start']}\t{r['end']}\t{rid}\t"
                    f"{r['composite_score']}\t+\n")
    print(f"[M07] Wrote BED: {path}")


def write_report(regions, meta, path):
    """Write bridge_report.json with full analysis details."""
    report = {
        "pipeline_version": "1.0",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "msa_sequences": meta["n_sequences"],
        "msa_length": meta["aln_length"],
        "genes_analyzed": meta["genes_analyzed"],
        "hyphy_genes": meta["hyphy_genes"],
        "total_conserved_regions": len(regions),
        "grade_summary": {g: sum(1 for r in regions if r["grade"] == g)
                          for g in "ABCD"},
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
    """Write primer_targets_manifest.yaml (tool-agnostic)."""
    manifest = {
        "pipeline_run_id": f"{build_name}_bridge",
        "source": "pandemic-preparedness-m07",
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
            "metrics": {
                "conservation_pct": r["conservation_pct"],
                "mean_entropy": r["mean_entropy"],
                "avg_fubar_purifying_prob": r["avg_fubar_purifying_prob"],
                "meme_episodic_sites": r["meme_episodic_sites"],
            },
            "flags": r["flags"],
            "recommendation": r["recommendation"],
        })

    with open(path, "w") as f:
        yaml.dump(manifest, f, default_flow_style=False,
                  sort_keys=False, allow_unicode=True)
    print(f"[M07] Wrote manifest: {path}")


# Main

def main():
    parser = argparse.ArgumentParser(
        description="M07 Bridge Builder — Conservation & Selection Integrator")
    parser.add_argument("--msa", required=True,
                        help="Aligned FASTA (MSA from augur align / M03)")
    parser.add_argument("--gff", required=True,
                        help="Gene annotation GFF3 file")
    parser.add_argument("--hyphy-dir", required=True,
                        help="Directory containing fubar_*.json, meme_*.json, gard_*.json")
    parser.add_argument("--reference", required=True,
                        help="Reference genome FASTA")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for BED, JSON, YAML")
    parser.add_argument("--config", required=True,
                        help="Pipeline config.yaml path")
    parser.add_argument("--build-name", default="unknown",
                        help="Build name for manifest identification")
    args = parser.parse_args()

    # Load config
    with open(args.config, "r") as f:
        full_config = yaml.safe_load(f)

    bridge = full_config.get("bridge", {})
    window_size = bridge.get("window_size", 150)
    entropy_threshold = bridge.get("entropy_threshold", 0.01)
    min_region_length = bridge.get("min_region_length", 100)
    buffer_bp = bridge.get("buffer_bp", 100)
    gard_buffer_bp = bridge.get("gard_buffer_bp", 50)
    weights = bridge.get("weights", {})
    tiers_config = bridge.get("tiers", {})

    target_genes = full_config.get("hyphy", {}).get("target_genes", [])
    hyphy_dir = Path(args.hyphy_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Reference ID (strip version suffix for matching)
    ref_records = list(SeqIO.parse(args.reference, "fasta"))
    ref_id = ref_records[0].id if ref_records else "unknown"

    print("[M07] ============================================")
    print("[M07]  Bridge Builder — Conservation × Selection")
    print("[M07] ============================================")
    print(f"[M07] Reference : {ref_id}")
    print(f"[M07] Window    : {window_size} bp")
    print(f"[M07] Threshold : {entropy_threshold}")
    print(f"[M07] Weights   : cons={weights.get('conservation', 0.35)}, "
          f"purify={weights.get('purifying', 0.35)}, "
          f"episod={weights.get('episodic_penalty', 0.15)}, "
          f"recomb={weights.get('recombination', 0.15)}")

    # ----- Step 1: Conservation Scanner -----
    print(f"\n[M07] Step 1: Conservation Scanner")
    regions, n_seq, aln_len = scan_conservation(
        args.msa, window_size, entropy_threshold, min_region_length)

    # Fail-safe: relax threshold if nothing found
    if not regions:
        relaxed = 0.05
        print(f"[M07] WARNING: No regions at {entropy_threshold}. "
              f"Relaxing to {relaxed} (survival mode)")
        regions, _, _ = scan_conservation(
            args.msa, window_size, relaxed, min_region_length)

    # ----- Parse GFF -----
    genes, utrs = parse_gff(args.gff)
    print(f"[M07] GFF: {len(genes)} genes, {len(utrs)} UTRs")

    # Resolve "ALL" target genes
    if "ALL" in target_genes:
        target_genes = [g["name"] for g in genes]

    # ----- Step 2: HyPhy Score Mapper -----
    print(f"\n[M07] Step 2: HyPhy Score Mapper")
    codon_data, gard_bps = map_hyphy_to_genome(genes, hyphy_dir, target_genes)
    print(f"[M07] Mapped {len(codon_data)} positions, "
          f"{len(gard_bps)} GARD breakpoints")

    # ----- Step 3: Intersection & Scoring -----
    print(f"\n[M07] Step 3: Intersection & Composite Scoring")
    scored = score_regions(regions, codon_data, gard_bps, genes, utrs,
                           weights, tiers_config, gard_buffer_bp)

    grades = {g: sum(1 for r in scored if r["grade"] == g) for g in "ABCD"}
    print(f"[M07] Scored {len(scored)} regions: "
          f"A={grades['A']} B={grades['B']} C={grades['C']} D={grades['D']}")

    # Warnings
    warnings = []
    if not scored:
        warnings.append("No conserved regions found. Check MSA quality or entropy threshold.")
    if not codon_data:
        warnings.append("No HyPhy data loaded. All scores are conservation-only.")
    if gard_bps:
        warnings.append(f"{len(gard_bps)} recombination breakpoints detected. "
                        f"Regions near breakpoints are penalized.")

    hyphy_genes_found = [g for g in target_genes
                         if (hyphy_dir / f"fubar_{g}.json").exists()]

    meta = {
        "n_sequences": n_seq,
        "aln_length": aln_len,
        "genes_analyzed": [g["name"] for g in genes] + [u["name"] for u in utrs],
        "hyphy_genes": hyphy_genes_found,
        "gard_breakpoints": [{"genome_pos": bp["genome_pos"], "gene": bp["gene"]}
                             for bp in gard_bps],
        "warnings": warnings,
    }

    # ----- Step 4: Output -----
    print(f"\n[M07] Step 4: Output Generation")
    write_bed(scored, ref_id, outdir / "conserved_regions.bed")
    write_report(scored, meta, outdir / "bridge_report.json")
    write_manifest(scored, args.build_name, args.reference, buffer_bp,
                   outdir / "primer_targets_manifest.yaml")

    print(f"\n[M07] ============================================")
    print(f"[M07]  Complete — {len(scored)} candidate regions")
    if scored:
        top = scored[0]
        print(f"[M07]  #1: {top['gene']} ({top['start']}-{top['end']}) "
              f"Grade {top['grade']}, Score {top['composite_score']}")
    print(f"[M07] ============================================")


if __name__ == "__main__":
    main()
