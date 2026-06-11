#!/usr/bin/env python3
"""
M10 — In-silico PCR Coverage Analysis Script
Validates candidate primers/probes against a sequence database using seqkit amplicon,
computes coverage breakdowns (by subtype, collection year, country), and generates reports.
"""

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path
import pandas as pd

# IUPAC degeneracy nucleotide mapping
IUPAC_GROUPS = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'}, 'U': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'}, 'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'}, '-': {'-'}
}

def match_base(p_base: str, a_base: str) -> bool:
    """Check if probe base matches amplicon base (IUPAC-aware)."""
    p_base = p_base.upper()
    a_base = a_base.upper()
    p_set = IUPAC_GROUPS.get(p_base, {p_base})
    a_set = IUPAC_GROUPS.get(a_base, {a_base})
    return not p_set.isdisjoint(a_set)

def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
        'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N', '-': '-'
    }
    return "".join(complement.get(base.upper(), base) for base in reversed(seq))

def min_mismatches_sliding_window(probe_seq: str, amplicon_seq: str) -> int:
    """Find the minimum mismatches of the probe in the amplicon sequence (IUPAC-aware)."""
    probe_len = len(probe_seq)
    amp_len = len(amplicon_seq)
    if amp_len < probe_len:
        return probe_len
    
    # 1. Exact string matching shortcut (works for most cases, saving compute)
    if probe_seq in amplicon_seq or reverse_complement(probe_seq) in amplicon_seq:
        return 0
        
    best_mismatch = probe_len
    
    # Check both original probe and its reverse complement
    for p_seq in [probe_seq, reverse_complement(probe_seq)]:
        for i in range(amp_len - probe_len + 1):
            window = amplicon_seq[i:i + probe_len]
            mismatches = 0
            for j in range(probe_len):
                if not match_base(p_seq[j], window[j]):
                    mismatches += 1
                    if mismatches >= best_mismatch:
                        break  # Early exit for this window
            if mismatches < best_mismatch:
                best_mismatch = mismatches
            if best_mismatch == 0:
                return 0
                
    return best_mismatch

def run_seqkit_amplicon(fwd: str, rev: str, max_mismatch: int, fasta_path: str) -> dict:
    """Run seqkit amplicon and return a dict of seq_id -> amplicon_seq."""
    cmd = [
        "seqkit", "amplicon",
        "-F", fwd,
        "-R", rev,
        "-m", str(max_mismatch),
        "--bed",
        str(fasta_path)
    ]
    
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running seqkit: {e.stderr}", file=sys.stderr)
        raise e
        
    matches = {}
    for line in result.stdout.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) >= 7:
            seq_id = parts[0]
            amp_seq = parts[6]
            matches[seq_id] = amp_seq
    return matches

def get_all_fasta_headers(fasta_path: str) -> list:
    """Fast extraction of all sequence headers from a FASTA file."""
    headers = []
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                # Strip the '>' and whitespace
                headers.append(line[1:].strip())
    return headers

def find_metadata_for_header(header: str, strain_to_meta: dict) -> dict:
    """Match sequence header to a row in metadata dict."""
    if header in strain_to_meta:
        return strain_to_meta[header]
    
    # Try GISAID style pipe splitting
    parts = header.split('|')
    for p in parts:
        p_strip = p.strip()
        if p_strip in strain_to_meta:
            return strain_to_meta[p_strip]
            
    # Try first word matching (NCBI style)
    first_word = header.split()[0]
    if first_word in strain_to_meta:
        return strain_to_meta[first_word]
        
    # Try case-insensitive matching
    for p in parts:
        p_lower = p.strip().lower()
        for k in strain_to_meta:
            if k.lower() == p_lower:
                return strain_to_meta[k]
                
    return None

def parse_year(date_str: str) -> str:
    """Extract 4-digit year from date string."""
    if not date_str or date_str == "Unknown":
        return "Unknown"
    match = re.search(r'\b(19|20)\d{2}\b', str(date_str))
    if match:
        return match.group(0)
    return "Unknown"

def main():
    parser = argparse.ArgumentParser(description="M10 — In-silico PCR Coverage Analysis")
    parser.add_argument("--m09-dirs", nargs="+", required=True,
                        help="List of M09 region output directories (e.g. results/BUILD/m09/M1_region_4). "
                             "Script will discover final_candidates.json inside each directory at runtime.")
    parser.add_argument("--sequences", required=True, help="Path to raw sequence FASTA database")
    parser.add_argument("--metadata", required=True, help="Path to pipeline metadata TSV")
    parser.add_argument("--gff", required=True, help="Path to GFF3 annotation file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--max-mismatch", type=int, default=3, help="Max mismatches allowed (default: 3)")
    parser.add_argument("--num-candidates", type=int, default=10, help="Number of candidates to evaluate per region")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. Load Metadata
    print(f"Loading metadata from {args.metadata}...")
    df_meta = pd.read_csv(args.metadata, sep='\t', dtype=str)
    
    # Fill empty/na values with default values
    df_meta['subtype'] = df_meta['subtype'].fillna('Unknown')
    df_meta['country'] = df_meta['country'].fillna('Unknown')
    df_meta['date'] = df_meta['date'].fillna('Unknown')
    
    strain_to_meta = {row['strain']: row.to_dict() for _, row in df_meta.iterrows()}
    print(f"Loaded {len(strain_to_meta)} metadata records.")

    # 2. Extract FASTA Headers
    print(f"Parsing FASTA headers from {args.sequences}...")
    fasta_headers = get_all_fasta_headers(args.sequences)
    total_seqs = len(fasta_headers)
    print(f"Found {total_seqs} sequences in database.")

    # Map sequences to metadata
    seq_mappings = []
    for header in fasta_headers:
        meta_row = find_metadata_for_header(header, strain_to_meta)
        if meta_row:
            seq_mappings.append({
                "header": header,
                "seq_id": header,
                "subtype": meta_row.get("subtype", "Unknown"),
                "year": parse_year(meta_row.get("date", "Unknown")),
                "country": meta_row.get("country", "Unknown")
            })
        else:
            # Fallback for sequences without metadata
            seq_mappings.append({
                "header": header,
                "seq_id": header,
                "subtype": "Unknown",
                "year": "Unknown",
                "country": "Unknown"
            })
    
    df_seq_map = pd.DataFrame(seq_mappings)

    # 3. Discover and Load Candidates from M09 directories
    # We discover final_candidates.json at runtime to handle cases where some regions
    # had no BLAST-passing candidates (and thus no final_candidates.json was written).
    all_candidates = []
    for m09_dir in args.m09_dirs:
        region_dir = Path(m09_dir)
        cand_path = region_dir / "final_candidates.json"

        if not cand_path.exists():
            print(f"[INFO] No final_candidates.json in {region_dir.name} - no BLAST-passing candidates for this region. Skipping.")
            continue

        # Region name is the directory basename (e.g. M1_region_4)
        region_name = region_dir.name
        print(f"Loading candidates for region '{region_name}' from {cand_path.name}...")

        try:
            with open(cand_path, "r", encoding="utf-8") as f:
                candidates = json.load(f)
        except Exception as e:
            print(f"Warning: Failed to parse {cand_path}: {e}. Skipping.")
            continue

        # Limit to top N per region
        for idx, c in enumerate(candidates[:args.num_candidates]):
            c["region"] = region_name
            c["gene"] = region_name.split("_region_")[0]
            # Create a unique ID for this candidate
            c["cand_id"] = f"{region_name}_Cand_{c['rank']}"
            all_candidates.append(c)

    if not all_candidates:
        print("[WARNING] No candidates found across all M09 regions. Nothing to evaluate.")
        print("    This usually means all regions had BLAST specificity scores below the threshold.")
        # Write empty outputs so Snakemake output targets are satisfied
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "coverage_summary.tsv").write_text("cand_id\tregion\tgene\trank\tfwd_seq\trev_seq\tprobe_seq\ttotal\tamplified\tdetected\tamplified_pct\tdetected_pct\n")
        (outdir / "coverage_summary.json").write_text("[]")
        (outdir / "coverage_report.md").write_text("# M10 In-silico PCR Coverage Report\n\nNo candidates passed BLAST specificity thresholds. No PCR simulation performed.\n")
        (outdir / "coverage_heatmap.html").write_text("<html><body><h1>No candidates available</h1></body></html>")
        (outdir / "m10_candidates.gff3").write_text("##gff-version 3\n")
        sys.exit(0)
        
    print(f"Evaluating {len(all_candidates)} candidates across all regions.")

    # 4. PCR Simulation & Probe Binding Check
    results_list = []
    summary_data = []

    for idx, cand in enumerate(all_candidates):
        cand_id = cand["cand_id"]
        fwd = cand["fwd_seq"]
        rev = cand["rev_seq"]
        probe = cand.get("probe_seq", "")
        region = cand["region"]
        gene = cand["gene"]
        rank = cand["rank"]
        
        print(f"[{idx+1}/{len(all_candidates)}] Simulating PCR for {cand_id}...")
        
        # Run seqkit amplicon for mismatch 0, 1, 2, 3
        # We determine the minimum primer mismatch for each sequence.
        # To optimize, we run mismatch 3 to get all matches, then check matching ones at lower mismatch levels.
        # This reduces seqkit invocations!
        matches_m3 = run_seqkit_amplicon(fwd, rev, 3, args.sequences)
        matches_m2 = run_seqkit_amplicon(fwd, rev, 2, args.sequences)
        matches_m1 = run_seqkit_amplicon(fwd, rev, 1, args.sequences)
        matches_m0 = run_seqkit_amplicon(fwd, rev, 0, args.sequences)
        
        matched_seqs = set(matches_m3.keys())
        
        # Evaluate each sequence's mismatch levels
        cand_results = []
        for seq in seq_mappings:
            seq_id = seq["seq_id"]
            
            # Primer mismatch level
            if seq_id in matches_m0:
                p_mismatch = 0
            elif seq_id in matches_m1:
                p_mismatch = 1
            elif seq_id in matches_m2:
                p_mismatch = 2
            elif seq_id in matches_m3:
                p_mismatch = 3
            else:
                p_mismatch = 99  # No match
                
            # Probe mismatch level
            prb_mismatch = 99
            if seq_id in matched_seqs:
                amp_seq = matches_m3[seq_id]
                if probe:
                    prb_mismatch = min_mismatches_sliding_window(probe, amp_seq)
                else:
                    prb_mismatch = 0  # No probe required → 0 mismatch
            
            is_amplified = (p_mismatch <= args.max_mismatch)
            is_detected = is_amplified and (prb_mismatch <= args.max_mismatch)
            
            cand_results.append({
                "seq_id": seq_id,
                "subtype": seq["subtype"],
                "year": seq["year"],
                "country": seq["country"],
                "primer_mismatch": p_mismatch,
                "probe_mismatch": prb_mismatch,
                "is_amplified": is_amplified,
                "is_detected": is_detected
            })
            
        df_res = pd.DataFrame(cand_results)
        
        # Calculate global metrics
        total = len(df_res)
        amplified = df_res["is_amplified"].sum()
        detected = df_res["is_detected"].sum()
        amp_pct = round((amplified / total * 100), 2) if total > 0 else 0.0
        det_pct = round((detected / total * 100), 2) if total > 0 else 0.0
        
        print(f"  -> Amplified: {amplified}/{total} ({amp_pct}%), Detected: {detected}/{total} ({det_pct}%)")
        
        # Store global results
        cand_summary = {
            "cand_id": cand_id,
            "region": region,
            "gene": gene,
            "rank": rank,
            "fwd_seq": fwd,
            "rev_seq": rev,
            "probe_seq": probe,
            "total_seqs": total,
            "amplified_seqs": int(amplified),
            "detected_seqs": int(detected),
            "amplified_pct": amp_pct,
            "detected_pct": det_pct,
            "breakdowns": {}
        }
        
        # Calculate breakdowns by Subtype, Year, Country
        for col in ["subtype", "year", "country"]:
            grouped = df_res.groupby(col).agg(
                total=('seq_id', 'count'),
                amplified=('is_amplified', 'sum'),
                detected=('is_detected', 'sum')
            ).reset_index()
            grouped['amplified_pct'] = (grouped['amplified'] / grouped['total'] * 100).round(2)
            grouped['detected_pct'] = (grouped['detected'] / grouped['total'] * 100).round(2)
            
            # Sort breakdowns logically
            if col == "year":
                grouped = grouped.sort_values(by="year", ascending=False)
            else:
                grouped = grouped.sort_values(by="total", ascending=False)
                
            cand_summary["breakdowns"][col] = grouped.to_dict(orient="records")
            
        summary_data.append(cand_summary)
        
        # Flattened row for TSV output
        results_list.append({
            "cand_id": cand_id,
            "region": region,
            "gene": gene,
            "rank": rank,
            "fwd_seq": fwd,
            "rev_seq": rev,
            "probe_seq": probe,
            "total": total,
            "amplified": int(amplified),
            "detected": int(detected),
            "amplified_pct": amp_pct,
            "detected_pct": det_pct
        })

    # 5. Output TSV and JSON summaries
    df_summary = pd.DataFrame(results_list)
    df_summary.to_csv(outdir / "coverage_summary.tsv", sep="\t", index=False)
    
    with open(outdir / "coverage_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary_data, f, indent=2)
        
    print(f"Wrote summaries to {outdir / 'coverage_summary.tsv'} and {outdir / 'coverage_summary.json'}.")

    # 6. Generate GFF3 genomic candidate track with coverage attributes
    gff3_path = outdir / "m10_candidates.gff3"
    with open(gff3_path, "w", encoding="utf-8") as f:
        f.write("##gff-version 3\n")
        # Find reference ID from GFF file or metadata (standard NCBI)
        ref_id = "NC_002016.1"
        if os.path.exists(args.gff):
            with open(args.gff, "r") as gff_f:
                for line in gff_f:
                    if line.startswith("#"):
                        if "sequence-region" in line:
                            ref_id = line.strip().split()[1]
                            f.write(line)
                            break
        
        f.write(f"##sequence-region {ref_id} 1 15000\n")
        for cand in all_candidates:
            cand_id = cand["cand_id"]
            gene = cand["gene"]
            rank = cand["rank"]
            
            # Find matching coverage data
            cov_info = next(c for c in summary_data if c["cand_id"] == cand_id)
            det_pct = cov_info["detected_pct"]
            amp_pct = cov_info["amplified_pct"]
            
            amp_start = cand.get("amplicon_start")
            amp_end = cand.get("amplicon_end")
            fwd_start = cand.get("fwd_start")
            fwd_end = cand.get("fwd_end")
            rev_start = cand.get("rev_start")
            rev_end = cand.get("rev_end")
            prb_start = cand.get("probe_start")
            prb_end = cand.get("probe_end")
            
            amp_len = amp_end - amp_start + 1 if (amp_start is not None and amp_end is not None) else 0
            fwd_len = fwd_end - fwd_start + 1 if (fwd_start is not None and fwd_end is not None) else 0
            rev_len = rev_end - rev_start + 1 if (rev_start is not None and rev_end is not None) else 0
            prb_len = prb_end - prb_start + 1 if (prb_start is not None and prb_end is not None) else 0
            
            # GFF coordinates are 1-based, inclusive
            if amp_start is not None and amp_end is not None:
                f.write(f"{ref_id}\tPrimerLab_M10\ttranscript\t{amp_start+1}\t{amp_end+1}\t.\t+\t.\tID={cand_id}_amp;Name={cand_id}_({amp_len}_bp);Coverage={det_pct}%;AmpCoverage={amp_pct}%\n")
            if fwd_start is not None and fwd_end is not None:
                f.write(f"{ref_id}\tPrimerLab_M10\tprimer\t{fwd_start+1}\t{fwd_end+1}\t.\t+\t.\tID={cand_id}_F;Name={cand_id}_Fwd_({fwd_len}_bp);Parent={cand_id}_amp\n")
            if rev_start is not None and rev_end is not None:
                f.write(f"{ref_id}\tPrimerLab_M10\tprimer\t{rev_start+1}\t{rev_end+1}\t.\t-\t.\tID={cand_id}_R;Name={cand_id}_Rev_({rev_len}_bp);Parent={cand_id}_amp\n")
            if prb_start is not None and prb_end is not None:
                f.write(f"{ref_id}\tPrimerLab_M10\tprobe\t{prb_start+1}\t{prb_end+1}\t.\t+\t.\tID={cand_id}_P;Name={cand_id}_Prb_({prb_len}_bp);Parent={cand_id}_amp\n")

    print(f"Wrote GFF3 track to {gff3_path}.")

    # 7. Generate Interactive Heatmap HTML Dashboard
    generate_heatmap_dashboard(summary_data, outdir / "coverage_heatmap.html")

    # 8. Generate coverage_report.md Report
    generate_markdown_report(summary_data, outdir / "coverage_report.md")

def generate_heatmap_dashboard(summary_data: list, output_path: Path):
    """Generate a premium interactive dark-mode HTML coverage dashboard."""
    # Convert summary data to JSON for JS inclusion
    data_json = json.dumps(summary_data, default=str)
    
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>🧬 In-silico PCR Coverage Dashboard</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;600;700&display=swap" rel="stylesheet">
    <style>
        :root {{
            --bg-color: #0f172a;
            --card-bg: rgba(30, 41, 59, 0.7);
            --border-color: rgba(255, 255, 255, 0.1);
            --accent-color: #10b981;
            --accent-glow: rgba(16, 185, 129, 0.15);
            --text-primary: #f8fafc;
            --text-secondary: #94a3b8;
        }}

        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}

        body {{
            font-family: 'Outfit', sans-serif;
            background-color: var(--bg-color);
            color: var(--text-primary);
            padding: 2rem;
            min-height: 100vh;
        }}

        .dashboard-container {{
            max-width: 1400px;
            margin: 0 auto;
        }}

        header {{
            margin-bottom: 2rem;
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 1.5rem;
        }}

        h1 {{
            font-size: 2.2rem;
            font-weight: 700;
            background: linear-gradient(135deg, #34d399 0%, #059669 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            margin-bottom: 0.5rem;
        }}

        .subtitle {{
            color: var(--text-secondary);
            font-size: 1rem;
        }}

        .controls {{
            display: flex;
            align-items: center;
            gap: 1rem;
            margin-bottom: 2rem;
            background: var(--card-bg);
            padding: 1rem;
            border-radius: 12px;
            border: 1px solid var(--border-color);
            backdrop-filter: blur(10px);
        }}

        label {{
            font-weight: 600;
            color: var(--text-secondary);
        }}

        select {{
            background-color: #1e293b;
            color: var(--text-primary);
            border: 1px solid var(--border-color);
            padding: 0.5rem 1rem;
            border-radius: 8px;
            font-family: 'Outfit', sans-serif;
            font-size: 1rem;
            cursor: pointer;
            outline: none;
            transition: all 0.3s ease;
        }}

        select:focus {{
            border-color: var(--accent-color);
            box-shadow: 0 0 10px var(--accent-glow);
        }}

        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 1.5rem;
            margin-bottom: 2rem;
        }}

        .stat-card {{
            background: var(--card-bg);
            border: 1px solid var(--border-color);
            border-radius: 12px;
            padding: 1.5rem;
            backdrop-filter: blur(10px);
            transition: transform 0.3s ease, border-color 0.3s ease;
        }}

        .stat-card:hover {{
            transform: translateY(-5px);
            border-color: var(--accent-color);
        }}

        .stat-label {{
            font-size: 0.85rem;
            color: var(--text-secondary);
            text-transform: uppercase;
            letter-spacing: 0.05em;
            margin-bottom: 0.5rem;
        }}

        .stat-value {{
            font-size: 1.8rem;
            font-weight: 700;
        }}

        .stat-value.highlight {{
            color: var(--accent-color);
            text-shadow: 0 0 10px var(--accent-glow);
        }}

        .seq-info {{
            font-family: monospace;
            font-size: 0.9rem;
            background: #0f172a;
            padding: 0.25rem 0.5rem;
            border-radius: 4px;
            border: 1px solid rgba(255,255,255,0.05);
            word-break: break-all;
        }}

        .main-layout {{
            display: grid;
            grid-template-columns: 1fr;
            gap: 2rem;
        }}

        @media (min-width: 1024px) {{
            .main-layout {{
                grid-template-columns: 1fr 1fr;
            }}
            .heatmap-full {{
                grid-column: span 2;
            }}
        }}

        .card {{
            background: var(--card-bg);
            border: 1px solid var(--border-color);
            border-radius: 16px;
            padding: 1.5rem;
            backdrop-filter: blur(10px);
        }}

        .card-title {{
            font-size: 1.2rem;
            font-weight: 600;
            margin-bottom: 1.5rem;
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 0.75rem;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 0.5rem;
        }}

        th, td {{
            padding: 0.75rem 1rem;
            text-align: left;
            border-bottom: 1px solid rgba(255,255,255,0.05);
        }}

        th {{
            color: var(--text-secondary);
            font-weight: 600;
            font-size: 0.9rem;
            text-transform: uppercase;
        }}

        .heatmap-cell {{
            font-weight: 700;
            text-align: center;
            border-radius: 6px;
            color: #0f172a;
            transition: all 0.2s ease;
        }}

        .heatmap-cell:hover {{
            transform: scale(1.05);
        }}

        .badge {{
            padding: 0.25rem 0.5rem;
            border-radius: 6px;
            font-size: 0.8rem;
            font-weight: 600;
            background: var(--accent-glow);
            color: var(--accent-color);
            border: 1px solid var(--accent-color);
        }}
    </style>
</head>
<body>
    <div class="dashboard-container">
        <header>
            <h1>🦠 PCR In-silico Validation Heatmap</h1>
            <div class="subtitle">Validasi diagnostik komprehensif terhadap database sekuens mentah.</div>
        </header>

        <div class="controls">
            <label for="candidateSelect">Select Candidate:</label>
            <select id="candidateSelect" onchange="updateDashboard(this.value)">
                <!-- Dropdown items dynamically loaded -->
            </select>
        </div>

        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-label">Candidate Region / Rank</div>
                <div id="statRegionRank" class="stat-value">Region - Rank</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Total Sequences Evaluated</div>
                <div id="statTotalSeqs" class="stat-value">0</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Overall Amplification (Primers)</div>
                <div id="statAmpCoverage" class="stat-value highlight">0%</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Overall Detection (Primers + Probe)</div>
                <div id="statDetCoverage" class="stat-value highlight">0%</div>
            </div>
        </div>

        <div class="stats-grid" style="grid-template-columns: 1fr 1fr 1fr; margin-bottom: 2rem;">
            <div class="stat-card">
                <div class="stat-label">Forward Primer Sequence (5' → 3')</div>
                <div id="seqFwd" class="seq-info">-</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Reverse Primer Sequence (5' → 3')</div>
                <div id="seqRev" class="seq-info">-</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Probe Sequence (5' → 3')</div>
                <div id="seqProbe" class="seq-info">-</div>
            </div>
        </div>

        <div class="main-layout">
            <div class="card">
                <div class="card-title">
                    <span>Subtype Coverage Distribution</span>
                    <span class="badge">EpiFlu Subtypes</span>
                </div>
                <table id="subtypeTable">
                    <thead>
                        <tr>
                            <th>Subtype</th>
                            <th>Total Seqs</th>
                            <th>Amplified (%)</th>
                            <th>Detected (%)</th>
                        </tr>
                    </thead>
                    <tbody></tbody>
                </table>
            </div>

            <div class="card">
                <div class="card-title">
                    <span>Collection Year Breakdown</span>
                    <span class="badge">Temporal Coverage</span>
                </div>
                <table id="yearTable">
                    <thead>
                        <tr>
                            <th>Year</th>
                            <th>Total Seqs</th>
                            <th>Amplified (%)</th>
                            <th>Detected (%)</th>
                        </tr>
                    </thead>
                    <tbody></tbody>
                </table>
            </div>

            <div class="card heatmap-full">
                <div class="card-title">
                    <span>Country Coverage Breakdown</span>
                    <span class="badge">Geographical Coverage</span>
                </div>
                <table id="countryTable">
                    <thead>
                        <tr>
                            <th>Country</th>
                            <th>Total Seqs</th>
                            <th>Amplified (%)</th>
                            <th>Detected (%)</th>
                        </tr>
                    </thead>
                    <tbody></tbody>
                </table>
            </div>
        </div>
    </div>

    <script>
        // Load data from python
        const candidatesData = {data_json};

        // Initialize select dropdown
        const select = document.getElementById("candidateSelect");
        candidatesData.forEach((cand, idx) => {{
            const option = document.createElement("option");
            option.value = idx;
            option.textContent = `${{cand.cand_id}} (${{cand.detected_pct}}% det)`;
            select.appendChild(option);
        }});

        function getCellColor(pct) {{
            if (pct >= 95.0) {{
                return `background-color: rgba(16, 185, 129, 0.95); color: #0f172a;`;
            }} else if (pct >= 85.0) {{
                return `background-color: rgba(52, 211, 153, 0.8); color: #0f172a;`;
            }} else if (pct >= 50.0) {{
                return `background-color: rgba(245, 158, 11, 0.9); color: #0f172a;`;
            }} else {{
                return `background-color: rgba(239, 68, 68, 0.9); color: #f8fafc;`;
            }}
        }}

        function populateTable(tableId, data, categoryKey) {{
            const tbody = document.querySelector(`#${{tableId}} tbody`);
            tbody.innerHTML = "";
            
            // Limit geographical list to top 25 countries to prevent screen bloat
            let displayData = data;
            if (tableId === "countryTable" && data.length > 25) {{
                displayData = data.slice(0, 25);
            }}

            displayData.forEach(row => {{
                const tr = document.createElement("tr");
                
                const tdCat = document.createElement("td");
                tdCat.textContent = row[categoryKey];
                tr.appendChild(tdCat);
                
                const tdTotal = document.createElement("td");
                tdTotal.textContent = row.total.toLocaleString();
                tr.appendChild(tdTotal);
                
                const tdAmp = document.createElement("td");
                tdAmp.className = "heatmap-cell";
                tdAmp.style = getCellColor(row.amplified_pct);
                tdAmp.textContent = `${{row.amplified_pct}}% (${{row.amplified.toLocaleString()}})`;
                tr.appendChild(tdAmp);
                
                const tdDet = document.createElement("td");
                tdDet.className = "heatmap-cell";
                tdDet.style = getCellColor(row.detected_pct);
                tdDet.textContent = `${{row.detected_pct}}% (${{row.detected.toLocaleString()}})`;
                tr.appendChild(tdDet);
                
                tbody.appendChild(tr);
            }});
            
            if (tableId === "countryTable" && data.length > 25) {{
                const tr = document.createElement("tr");
                const td = document.createElement("td");
                td.colSpan = 4;
                td.style = "text-align: center; color: var(--text-secondary); font-style: italic;";
                td.textContent = `Showing top 25 countries out of ${{data.length}} total.`;
                tr.appendChild(td);
                tbody.appendChild(tr);
            }}
        }}

        function updateDashboard(index) {{
            const cand = candidatesData[index];
            if (!cand) return;

            // Update stats cards
            document.getElementById("statRegionRank").textContent = `${{cand.region}} - Rank ${{cand.rank}}`;
            document.getElementById("statTotalSeqs").textContent = cand.total_seqs.toLocaleString();
            document.getElementById("statAmpCoverage").textContent = `${{cand.amplified_pct}}%`;
            document.getElementById("statDetCoverage").textContent = `${{cand.detected_pct}}%`;

            // Update detailed sequences
            document.getElementById("seqFwd").textContent = cand.fwd_seq;
            document.getElementById("seqRev").textContent = cand.rev_seq;
            document.getElementById("seqProbe").textContent = cand.probe_seq || "— (No Probe)";

            // Populate breakdown tables
            populateTable("subtypeTable", cand.breakdowns.subtype, "subtype");
            populateTable("yearTable", cand.breakdowns.year, "year");
            populateTable("countryTable", cand.breakdowns.country, "country");
        }}

        // Initialize dashboard with the first candidate
        if (candidatesData.length > 0) {{
            updateDashboard(0);
        }}
    </script>
</body>
</html>
"""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    print(f"Generated Heatmap Dashboard: {output_path}")

def generate_markdown_report(summary_data: list, output_path: Path):
    """Generate an executive summary report of the coverage analysis."""
    lines = [
        "# 📋 In-silico PCR Coverage Analysis Report",
        "",
        "Laporan ini menyajikan hasil validasi *in-silico* PCR dari kandidat primer-probe yang berhasil lolos BLAST specificity checks terhadap seluruh sekuens input database.",
        "",
        "## Overall Summary",
        "",
        "| Candidate ID | Target Region | Forward Sequence (5'-3') | Reverse Sequence (5'-3') | Probe Sequence (5'-3') | Amplification Coverage (%) | Detection Coverage (%) |",
        "|--------------|---------------|--------------------------|--------------------------|------------------------|--------------------------|------------------------|",
    ]
    
    for cand in summary_data:
        probe_str = f"`{cand['probe_seq']}`" if cand['probe_seq'] else "—"
        lines.append(
            f"| **{cand['cand_id']}** | {cand['region']} | `{cand['fwd_seq']}` | `{cand['rev_seq']}` | {probe_str} | **{cand['amplified_pct']}%** | **{cand['detected_pct']}%** |"
        )
        
    lines.append("")
    lines.append("## Detailed Subtype Coverage Breakdown")
    lines.append("")
    
    for cand in summary_data:
        lines.append(f"### {cand['cand_id']} (Region: {cand['region']})")
        lines.append("")
        lines.append("| Subtype | Total Sequences | Amplified Coverage (%) | Detection Coverage (%) |")
        lines.append("|---------|-----------------|------------------------|------------------------|")
        for sub in cand["breakdowns"]["subtype"]:
            lines.append(f"| {sub['subtype']} | {sub['total']:,} | {sub['amplified_pct']}% | {sub['detected_pct']}% |")
        lines.append("")
        
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Generated Markdown Report: {output_path}")

if __name__ == "__main__":
    main()
