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


def assess_raa_binding(
    primer_seq: str,
    amp_seq: str,
    critical_3prime_bp: int = 8,
    max_3prime_mm: int = 1,
    max_internal_mm: int = 4,
) -> tuple:
    """Assess if a primer can bind in RAA/RPA conditions.

    RAA biology:
    - Recombinase (RecA-like) loads on the primer and drives strand invasion.
    - INTERNAL mismatches (5' half) are tolerated because recombinase can bypass them.
    - 3' END mismatches are CRITICAL: polymerase needs accurate 3' annealing for extension.

    Returns:
        (is_raa_bound, 3prime_mm, internal_mm, total_mm)
    """
    primer_len = len(primer_seq)
    amp_len = len(amp_seq)
    if amp_len < primer_len:
        return False, primer_len, primer_len, primer_len

    best_3prime_mm = primer_len
    best_internal_mm = primer_len
    best_total_mm = primer_len

    for p_seq in [primer_seq, reverse_complement(primer_seq)]:
        for i in range(amp_len - primer_len + 1):
            window = amp_seq[i : i + primer_len]
            total_mm = 0
            three_prime_mm = 0
            internal_mm = 0
            internal_len = primer_len - critical_3prime_bp

            for j in range(primer_len):
                mismatch = not match_base(p_seq[j], window[j])
                if mismatch:
                    total_mm += 1
                    # Positions [internal_len .. end] are the 3' end
                    if j >= internal_len:
                        three_prime_mm += 1
                    else:
                        internal_mm += 1

            if total_mm < best_total_mm:
                best_total_mm = total_mm
                best_3prime_mm = three_prime_mm
                best_internal_mm = internal_mm
            if best_total_mm == 0:
                break
        if best_total_mm == 0:
            break

    is_raa_bound = (best_3prime_mm <= max_3prime_mm) and (best_internal_mm <= max_internal_mm)
    return is_raa_bound, best_3prime_mm, best_internal_mm, best_total_mm


import tempfile
import itertools

def is_degenerate(seq: str) -> bool:
    return any(c.upper() not in {'A', 'T', 'G', 'C', 'U', '-'} for c in seq)

def expand_degenerate_sequence(seq: str) -> list:
    """Expand a degenerate DNA sequence containing IUPAC codes into all possible combinations."""
    iupac_dict = {
        'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'], 'U': ['T'],
        'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
        'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 
        'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T']
    }
    bases_lists = [iupac_dict.get(c.upper(), [c]) for c in seq]
    return [''.join(p) for p in itertools.product(*bases_lists)]

def run_seqkit_amplicon(fwd: str, rev: str, max_mismatch: int, fasta_path: str, threads: int = 4) -> dict:
    """Run seqkit amplicon and return a dict of seq_id -> amplicon_seq.
    Supports degenerate primers by expanding them and using the -p/--primer-file flag.
    """
    temp_tsv_path = None
    
    if is_degenerate(fwd) or is_degenerate(rev):
        # Expand degenerate primers
        fwd_list = expand_degenerate_sequence(fwd)
        rev_list = expand_degenerate_sequence(rev)
        
        # Limit combinations to a safe maximum to avoid exploding complexity
        combinations = list(itertools.product(fwd_list, rev_list))
        
        # Create temp TSV file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as temp_f:
            for idx, (f_comb, r_comb) in enumerate(combinations):
                temp_f.write(f"comb_{idx}\t{f_comb}\t{r_comb}\n")
            temp_tsv_path = temp_f.name
            
        cmd = [
            "seqkit", "amplicon",
            "-p", temp_tsv_path,
            "-m", str(max_mismatch),
            "-j", str(threads),
            "--bed",
            str(fasta_path)
        ]
    else:
        cmd = [
            "seqkit", "amplicon",
            "-F", fwd,
            "-R", rev,
            "-m", str(max_mismatch),
            "-j", str(threads),
            "--bed",
            str(fasta_path)
        ]
    
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running seqkit: {e.stderr}", file=sys.stderr)
        if temp_tsv_path and os.path.exists(temp_tsv_path):
            os.remove(temp_tsv_path)
        raise e
    finally:
        if temp_tsv_path and os.path.exists(temp_tsv_path):
            os.remove(temp_tsv_path)
        
    matches = {}
    for line in result.stdout.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) >= 7:
            seq_id = parts[0]
            amp_seq = parts[6]
            # If multiple combinations match, keep the first one (or shortest)
            if seq_id not in matches or len(amp_seq) < len(matches[seq_id]):
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
    parser.add_argument("--max-mismatch", type=int, default=3, help="Max mismatches allowed for PCR mode (default: 3)")
    parser.add_argument("--num-candidates", type=int, default=10, help="Number of candidates to evaluate per region")
    parser.add_argument("--threads", type=int, default=None,
                        help="Total CPU threads to use. Splits between seqkit (-j) and candidate parallelism. Default: auto-detect (os.cpu_count).")
    # RAA-specific binding parameters
    parser.add_argument("--raa-3prime-bp", type=int, default=8,
                        help="Number of 3' end bases considered critical for RAA binding (default: 8)")
    parser.add_argument("--raa-max-3prime-mm", type=int, default=1,
                        help="Max mismatches allowed in the RAA critical 3' region (default: 1)")
    parser.add_argument("--raa-max-internal-mm", type=int, default=4,
                        help="Max mismatches allowed in the RAA internal/5' region (default: 4)")
    args = parser.parse_args()
    
    # Resolve thread allocation
    total_cpus = args.threads or os.cpu_count() or 8
    # Give seqkit at least 4 threads per call, rest for candidate parallelism
    seqkit_threads = max(4, total_cpus // 4)
    worker_threads = max(1, total_cpus // seqkit_threads)
    print(f"[M10] Thread config: {total_cpus} total CPUs → {worker_threads} parallel candidates × {seqkit_threads} seqkit threads each")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. Load Metadata
    print(f"Loading metadata from {args.metadata}...")
    
    # Auto-detect delimiter: sniff first line to decide tab vs comma
    meta_path = Path(args.metadata)
    with open(meta_path, 'r', encoding='utf-8', errors='ignore') as _f:
        _first_line = _f.readline()
    _sep = '\t' if '\t' in _first_line else ','
    print(f"[M10] Detected metadata delimiter: {'TAB' if _sep == chr(9) else 'COMMA'}")
    
    df_meta = pd.read_csv(args.metadata, sep=_sep, dtype=str, on_bad_lines='warn')

    
    # Normalize column names to lowercase to handle case differences (Subtype vs subtype)
    df_meta.columns = [c.strip().lower() for c in df_meta.columns]
    
    # Fill empty/na values with default values (only if column exists)
    if 'subtype' in df_meta.columns:
        df_meta['subtype'] = df_meta['subtype'].fillna('Unknown')
    else:
        df_meta['subtype'] = 'Unknown'
        
    if 'country' in df_meta.columns:
        df_meta['country'] = df_meta['country'].fillna('Unknown')
    else:
        df_meta['country'] = 'Unknown'
        
    if 'date' in df_meta.columns:
        df_meta['date'] = df_meta['date'].fillna('Unknown')
    else:
        df_meta['date'] = 'Unknown'
    
    # Use 'strain' column as key; fall back to first column if not present
    if 'strain' not in df_meta.columns:
        print(f"[WARNING] 'strain' column not found in metadata. Using first column '{df_meta.columns[0]}' as key.")
        df_meta = df_meta.rename(columns={df_meta.columns[0]: 'strain'})
    
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
            # Create a unique ID for this candidate (use rank or index as fallback)
            rank_val = c.get('rank', idx + 1)
            c["cand_id"] = f"{region_name}_Cand_{rank_val}"
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

    # 4. PCR Simulation & Probe Binding Check — fully parallelized
    # Strategy: run seqkit ONCE per candidate (with max_mismatch), then use sliding-window
    # mismatch detection to determine exact mismatch level per-sequence. This cuts seqkit
    # calls from 4× down to 1× per candidate while still reporting per-mismatch breakdown.

    import concurrent.futures

    def _evaluate_candidate(idx_cand):
        """Worker function: evaluate one candidate. Returns (cand_summary, result_row) or None."""
        idx, cand = idx_cand
        cand_id = cand.get("cand_id", f"cand_{idx}")
        fwd = cand.get("fwd_seq") or cand.get("forward_primer_seq") or cand.get("fwd")
        rev = cand.get("rev_seq") or cand.get("reverse_primer_seq") or cand.get("rev")
        probe = cand.get("probe_seq") or cand.get("probe") or ""
        region = cand.get("region", "unknown")
        gene = cand.get("gene", region.split("_region_")[0] if "_region_" in region else region)
        rank = cand.get("rank", idx + 1)

        if not fwd or not rev:
            print(f"[SKIP] {cand_id}: missing fwd_seq or rev_seq.")
            return None

        print(f"[{idx+1}/{len(all_candidates)}] Simulating PCR for {cand_id} (seqkit -j {seqkit_threads})...")

        # ── Single seqkit call with max_mismatch ──────────────────────────────
        # seqkit returns amplicons for sequences where BOTH primers bind within max_mismatch.
        # We then compute per-sequence primer mismatch via Python sliding window for the
        # matched set to determine per-mismatch-level breakdown.
        try:
            matches_max = run_seqkit_amplicon(fwd, rev, args.max_mismatch, args.sequences, seqkit_threads)
        except Exception as e:
            print(f"[ERROR] seqkit failed for {cand_id}: {e}")
            matches_max = {}

        matched_seqs = set(matches_max.keys())

        # ── Per-sequence evaluation ───────────────────────────────────────────
        # CRITICAL: seqkit amplicon is the authoritative source for amplification.
        # Sequences in matched_seqs are DEFINITIVELY amplifiable within max_mismatch.
        # The sliding window classifies mismatch LEVEL for breakdown reporting ONLY —
        # it does NOT re-gate amplification. Previously this caused double-filtering
        # and catastrophically underreported coverage (e.g. H3N2 showed 0.43%).
        cand_results = []
        for seq in seq_mappings:
            seq_id = seq["seq_id"]

            if seq_id not in matched_seqs:
                # seqkit found no amplicon → genuinely not amplifiable
                is_amplified = False
                p_mismatch   = 99
                prb_mismatch = 99
                is_raa_amplified  = False
                raa_3p_mm         = 99
                raa_int_mm        = 99
                is_raa_detected   = False
                raa_probe_mm      = 99
            else:
                # seqkit confirmed amplification → trust it unconditionally for PCR
                is_amplified = True
                amp_seq = matches_max[seq_id]

                # ─ PCR mismatch level (for reporting, not gating) ─
                fwd_mm = min_mismatches_sliding_window(fwd, amp_seq)
                rev_mm = min_mismatches_sliding_window(reverse_complement(rev), amp_seq)
                p_mismatch = min(max(fwd_mm, rev_mm), args.max_mismatch)

                # ─ RAA binding assessment for both primers ─
                fwd_raa_ok, fwd_3p, fwd_int, _ = assess_raa_binding(
                    fwd, amp_seq,
                    critical_3prime_bp=args.raa_3prime_bp,
                    max_3prime_mm=args.raa_max_3prime_mm,
                    max_internal_mm=args.raa_max_internal_mm,
                )
                rev_raa_ok, rev_3p, rev_int, _ = assess_raa_binding(
                    reverse_complement(rev), amp_seq,
                    critical_3prime_bp=args.raa_3prime_bp,
                    max_3prime_mm=args.raa_max_3prime_mm,
                    max_internal_mm=args.raa_max_internal_mm,
                )
                is_raa_amplified = fwd_raa_ok and rev_raa_ok
                raa_3p_mm  = max(fwd_3p, rev_3p)
                raa_int_mm = max(fwd_int, rev_int)

                # ─ Probe mismatch (shared for both PCR & RAA probe detection) ─
                if probe:
                    prb_mismatch = min_mismatches_sliding_window(probe, amp_seq)
                    # For RAA: probe (exo-probe) has similar 3' critical region
                    probe_raa_ok, _, _, _ = assess_raa_binding(
                        probe, amp_seq,
                        critical_3prime_bp=args.raa_3prime_bp,
                        max_3prime_mm=args.raa_max_3prime_mm,
                        max_internal_mm=args.raa_max_internal_mm,
                    )
                    raa_probe_mm = prb_mismatch
                else:
                    prb_mismatch = 0
                    probe_raa_ok = True
                    raa_probe_mm = 0

            is_detected     = is_amplified     and (prb_mismatch <= args.max_mismatch)
            is_raa_detected = is_raa_amplified and probe_raa_ok


            cand_results.append({
                "seq_id":            seq_id,
                "subtype":           seq["subtype"],
                "year":              seq["year"],
                "country":           seq["country"],
                "primer_mismatch":   p_mismatch,
                "probe_mismatch":    prb_mismatch,
                # PCR
                "is_amplified":      is_amplified,
                "is_detected":       is_detected,
                # RAA
                "is_raa_amplified":  is_raa_amplified,
                "raa_3p_mm":         raa_3p_mm,
                "raa_int_mm":        raa_int_mm,
                "is_raa_detected":   is_raa_detected,
            })

        df_res = pd.DataFrame(cand_results)
        total     = len(df_res)
        # PCR metrics
        amplified     = df_res["is_amplified"].sum()
        detected      = df_res["is_detected"].sum()
        amp_pct       = round(amplified / total * 100, 2) if total > 0 else 0.0
        det_pct       = round(detected  / total * 100, 2) if total > 0 else 0.0
        # RAA metrics
        raa_amplified = df_res["is_raa_amplified"].sum()
        raa_detected  = df_res["is_raa_detected"].sum()
        raa_amp_pct   = round(raa_amplified / total * 100, 2) if total > 0 else 0.0
        raa_det_pct   = round(raa_detected  / total * 100, 2) if total > 0 else 0.0

        print(f"  [{cand_id}] PCR → Amplified: {amplified}/{total} ({amp_pct}%), Detected: {detected}/{total} ({det_pct}%)")
        print(f"  [{cand_id}] RAA → Amplified: {raa_amplified}/{total} ({raa_amp_pct}%), Detected: {raa_detected}/{total} ({raa_det_pct}%)")

        # ── Breakdown by subtype / year / country ─────────────────────────────
        breakdowns = {}
        for col in ["subtype", "year", "country"]:
            grp = df_res.groupby(col).agg(
                total=("seq_id", "count"),
                pcr_amplified=("is_amplified", "sum"),
                pcr_detected=("is_detected",  "sum"),
                raa_amplified=("is_raa_amplified", "sum"),
                raa_detected=("is_raa_detected",  "sum"),
            ).reset_index()
            grp["pcr_amplified_pct"] = (grp["pcr_amplified"] / grp["total"] * 100).round(2)
            grp["pcr_detected_pct"]  = (grp["pcr_detected"]  / grp["total"] * 100).round(2)
            grp["raa_amplified_pct"] = (grp["raa_amplified"] / grp["total"] * 100).round(2)
            grp["raa_detected_pct"]  = (grp["raa_detected"]  / grp["total"] * 100).round(2)
            grp = grp.sort_values("year" if col == "year" else "total",
                                   ascending=(col == "year"))
            breakdowns[col] = grp.to_dict(orient="records")

        cand_summary = {
            "cand_id": cand_id, "region": region, "gene": gene, "rank": rank,
            "fwd_seq": fwd, "rev_seq": rev, "probe_seq": probe,
            # PCR
            "total_seqs": total,
            "pcr_amplified_seqs": int(amplified), "pcr_amplified_pct": amp_pct,
            "pcr_detected_seqs":  int(detected),  "pcr_detected_pct":  det_pct,
            # RAA
            "raa_amplified_seqs": int(raa_amplified), "raa_amplified_pct": raa_amp_pct,
            "raa_detected_seqs":  int(raa_detected),  "raa_detected_pct":  raa_det_pct,
            "raa_params": {
                "critical_3prime_bp":  args.raa_3prime_bp,
                "max_3prime_mm":       args.raa_max_3prime_mm,
                "max_internal_mm":     args.raa_max_internal_mm,
            },
            "breakdowns": breakdowns,
        }
        result_row = {
            "cand_id": cand_id, "region": region, "gene": gene, "rank": rank,
            "fwd_seq": fwd, "rev_seq": rev, "probe_seq": probe, "total": total,
            # PCR
            "pcr_amplified": int(amplified), "pcr_detected": int(detected),
            "pcr_amplified_pct": amp_pct,    "pcr_detected_pct": det_pct,
            # RAA
            "raa_amplified": int(raa_amplified), "raa_detected": int(raa_detected),
            "raa_amplified_pct": raa_amp_pct,    "raa_detected_pct": raa_det_pct,
        }
        return (cand_summary, result_row)

    # ── Run candidates in parallel ────────────────────────────────────────────
    results_list = []
    summary_data = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=worker_threads) as pool:
        futures = {pool.submit(_evaluate_candidate, (i, c)): i
                   for i, c in enumerate(all_candidates)}
        for fut in concurrent.futures.as_completed(futures):
            out = fut.result()
            if out is not None:
                cand_summary, result_row = out
                summary_data.append(cand_summary)
                results_list.append(result_row)

    # Restore deterministic order (sorted by region then rank)
    summary_data.sort(key=lambda x: (x["region"], x["rank"]))
    results_list.sort(key=lambda x: (x["region"],  x["rank"]))

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
            cov_info = next((c for c in summary_data if c["cand_id"] == cand_id), {})
            det_pct = cov_info.get("pcr_detected_pct", cov_info.get("detected_pct", 0))
            amp_pct = cov_info.get("pcr_amplified_pct", cov_info.get("amplified_pct", 0))
            raa_det_pct = cov_info.get("raa_detected_pct", 0)
            raa_amp_pct = cov_info.get("raa_amplified_pct", 0)
            
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
            
            if amp_start is not None and amp_end is not None:
                f.write(f"{ref_id}\tPrimerLab_M10\ttranscript\t{amp_start+1}\t{amp_end+1}\t.\t+\t.\t"
                        f"ID={cand_id}_amp;Name={cand_id}_({amp_len}_bp);"
                        f"PCR_AmpCov={amp_pct}%;PCR_DetCov={det_pct}%;"
                        f"RAA_AmpCov={raa_amp_pct}%;RAA_DetCov={raa_det_pct}%;gene={gene}\n")
            if fwd_start is not None and fwd_end is not None:
                f.write(f"{ref_id}\tPrimerLab_M10\tprimer\t{fwd_start+1}\t{fwd_end+1}\t.\t+\t.\tID={cand_id}_F;Name={cand_id}_Fwd_({fwd_len}_bp);Parent={cand_id}_amp;gene={gene}\n")
            if rev_start is not None and rev_end is not None:
                f.write(f"{ref_id}\tPrimerLab_M10\tprimer\t{rev_start+1}\t{rev_end+1}\t.\t-\t.\tID={cand_id}_R;Name={cand_id}_Rev_({rev_len}_bp);Parent={cand_id}_amp;gene={gene}\n")
            if prb_start is not None and prb_end is not None:
                f.write(f"{ref_id}\tPrimerLab_M10\tprobe\t{prb_start+1}\t{prb_end+1}\t.\t+\t.\tID={cand_id}_P;Name={cand_id}_Prb_({prb_len}_bp);Parent={cand_id}_amp;gene={gene}\n")

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
                            <th>PCR Amp (%)</th>
                            <th>PCR Det (%)</th>
                            <th>RAA Amp (%)</th>
                            <th>RAA Det (%)</th>
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
                            <th>PCR Amp (%)</th>
                            <th>PCR Det (%)</th>
                            <th>RAA Amp (%)</th>
                            <th>RAA Det (%)</th>
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
                            <th>PCR Amp (%)</th>
                            <th>PCR Det (%)</th>
                            <th>RAA Amp (%)</th>
                            <th>RAA Det (%)</th>
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
            const detPct = cand.pcr_detected_pct ?? cand.detected_pct ?? 0;
            const raaDetPct = cand.raa_detected_pct ?? "—";
            option.textContent = `${{cand.cand_id}} (PCR Det: ${{detPct}}%, RAA Det: ${{raaDetPct}}%)`;
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
                tdAmp.style = getCellColor(row.pcr_amplified_pct ?? row.amplified_pct ?? 0);
                tdAmp.textContent = `${{(row.pcr_amplified_pct ?? row.amplified_pct ?? 0)}}%`;
                tr.appendChild(tdAmp);

                const tdDet = document.createElement("td");
                tdDet.className = "heatmap-cell";
                tdDet.style = getCellColor(row.pcr_detected_pct ?? row.detected_pct ?? 0);
                tdDet.textContent = `${{(row.pcr_detected_pct ?? row.detected_pct ?? 0)}}%`;
                tr.appendChild(tdDet);

                const tdRaaAmp = document.createElement("td");
                tdRaaAmp.className = "heatmap-cell";
                tdRaaAmp.style = getCellColor(row.raa_amplified_pct ?? 0);
                tdRaaAmp.textContent = `${{(row.raa_amplified_pct ?? "-")}}%`;
                tr.appendChild(tdRaaAmp);

                const tdRaaDet = document.createElement("td");
                tdRaaDet.className = "heatmap-cell";
                tdRaaDet.style = getCellColor(row.raa_detected_pct ?? 0);
                tdRaaDet.textContent = `${{(row.raa_detected_pct ?? "-")}}%`;
                tr.appendChild(tdRaaDet);

                tbody.appendChild(tr);
            }});
            
            if (tableId === "countryTable" && data.length > 25) {{
                const tr = document.createElement("tr");
                const td = document.createElement("td");
                td.colSpan = 6;
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
            const pcrAmp = cand.pcr_amplified_pct ?? cand.amplified_pct ?? 0;
            const raaAmp = cand.raa_amplified_pct ?? "—";
            document.getElementById("statAmpCoverage").textContent = `PCR: ${{pcrAmp}}% | RAA: ${{raaAmp}}%`;
            const pcrDet = cand.pcr_detected_pct ?? cand.detected_pct ?? 0;
            const raaDet = cand.raa_detected_pct ?? "—";
            document.getElementById("statDetCoverage").textContent = `PCR: ${{pcrDet}}% | RAA: ${{raaDet}}%`;

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
    """Generate an executive summary report comparing PCR vs RAA coverage."""
    lines = [
        "# 📋 In-silico Coverage Analysis Report: PCR vs RAA",
        "",
        "Laporan ini membandingkan cakupan amplifikasi antara dua mode:",
        "- **PCR:** Simulasi via `seqkit amplicon`, toleransi mismatch total ≤3",
        "- **RAA:** Model biologis — 3' end (8 bp) ≤1 mismatch, internal ≤4 mismatch",
        "",
        "> RAA lebih toleran mismatch internal (recombinase strand invasion), lebih ketat di 3' end.",
        "",
        "## 📊 Overall Summary — PCR vs RAA",
        "",
        "| Candidate ID | Region | PCR Amp (%) | PCR Det (%) | RAA Amp (%) | RAA Det (%) |",
        "|--------------|--------|-------------|-------------|-------------|-------------|" ,
    ]

    for cand in summary_data:
        pcr_amp = cand.get("pcr_amplified_pct", cand.get("amplified_pct", 0))
        pcr_det = cand.get("pcr_detected_pct",  cand.get("detected_pct", 0))
        raa_amp = cand.get("raa_amplified_pct", "n/a")
        raa_det = cand.get("raa_detected_pct",  "n/a")
        lines.append(
            f"| **{cand['cand_id']}** | {cand['region']} "
            f"| **{pcr_amp}%** | {pcr_det}% "
            f"| **{raa_amp}%** | {raa_det}% |"
        )

    lines.append("")
    lines.append("## 🔬 Primer Sequences")
    lines.append("")
    lines.append("| Candidate ID | Forward (5'→3') | Reverse (5'→3') | Probe (5'→3') |")
    lines.append("|--------------|-----------------|-----------------|---------------|")
    for cand in summary_data:
        probe_str = f"`{cand['probe_seq']}`" if cand.get("probe_seq") else "—"
        lines.append(
            f"| **{cand['cand_id']}** | `{cand['fwd_seq']}` | `{cand['rev_seq']}` | {probe_str} |"
        )

    lines.append("")
    lines.append("## 📈 Detailed Subtype Coverage: PCR vs RAA")
    lines.append("")

    for cand in summary_data:
        pcr_amp = cand.get("pcr_amplified_pct", cand.get("amplified_pct", 0))
        raa_amp = cand.get("raa_amplified_pct", "n/a")
        lines.append(f"### {cand['cand_id']} (Region: {cand['region']})")
        lines.append("")
        lines.append(f"**Overall:** PCR Amp = {pcr_amp}% | RAA Amp = {raa_amp}%")
        lines.append("")
        lines.append("| Subtype | N Sequences | PCR Amp (%) | PCR Det (%) | RAA Amp (%) | RAA Det (%) |")
        lines.append("|---------|-------------|-------------|-------------|-------------|-------------|")
        for sub in cand["breakdowns"]["subtype"]:
            pcr_a = sub.get("pcr_amplified_pct", sub.get("amplified_pct", 0))
            pcr_d = sub.get("pcr_detected_pct",  sub.get("detected_pct", 0))
            raa_a = sub.get("raa_amplified_pct", "n/a")
            raa_d = sub.get("raa_detected_pct",  "n/a")
            lines.append(
                f"| {sub.get('subtype', '?')} | {sub['total']:,} "
                f"| {pcr_a}% | {pcr_d}% "
                f"| {raa_a}% | {raa_d}% |"
            )
        lines.append("")

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Generated Markdown Report: {output_path}")

if __name__ == "__main__":
    main()
