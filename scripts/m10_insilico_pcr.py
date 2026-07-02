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


try:
    import numpy as _np
    _HAS_NUMPY = True
except ImportError:
    _np = None
    _HAS_NUMPY = False

# IUPAC bitmask for vectorized comparison: A=1,C=2,G=4,T=8
_IUPAC_MASK = {
    'A':1,'C':2,'G':4,'T':8,'U':8,'R':5,'Y':10,'S':6,'W':9,'K':12,
    'M':3,'B':14,'D':13,'H':11,'V':7,'N':15,'-':0,
}

def _assess_raa_numpy(primer_seq, amp_seq, critical_3prime_bp=8, max_3prime_mm=1, max_internal_mm=4):
    """Vectorized RAA binding assessment using numpy bitmask (~8x faster than pure Python)."""
    primer_len = len(primer_seq)
    amp_len    = len(amp_seq)
    if amp_len < primer_len:
        return False, primer_len, primer_len, primer_len
    internal_len = primer_len - critical_3prime_bp
    best_3p, best_int, best_tot = primer_len, primer_len, primer_len
    rc = reverse_complement(primer_seq)
    seqs = [primer_seq] if primer_seq == rc else [primer_seq, rc]
    for p_seq in seqs:
        p_arr   = _np.array([_IUPAC_MASK.get(b.upper(), 15) for b in p_seq], dtype=_np.uint8)
        amp_arr = _np.array([_IUPAC_MASK.get(b.upper(), 15) for b in amp_seq], dtype=_np.uint8)
        n_win = amp_len - primer_len + 1
        if n_win <= 0:
            continue
        s = amp_arr.strides[0]
        wins = _np.lib.stride_tricks.as_strided(amp_arr, shape=(n_win, primer_len), strides=(s, s))
        mm = ((p_arr[None, :] & wins) == 0).view(_np.uint8)
        tot_arr = mm.sum(axis=1)
        int_arr = mm[:, :internal_len].sum(axis=1)
        t3_arr  = mm[:, internal_len:].sum(axis=1)
        bi = int(tot_arr.argmin())
        tot = int(tot_arr[bi]); t3 = int(t3_arr[bi]); it = int(int_arr[bi])
        if tot < best_tot:
            best_tot = tot; best_3p = t3; best_int = it
        if best_tot == 0:
            break
    return (best_3p <= max_3prime_mm) and (best_int <= max_internal_mm), best_3p, best_int, best_tot

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
            if seq_id not in matches or len(amp_seq) < len(matches[seq_id]):
                matches[seq_id] = amp_seq
    return matches


def run_seqkit_amplicon_pool(
    fwd_list: list, rev_list: list,
    max_mismatch: int, fasta_path: str, threads: int = 4,
    cand_id: str = "unknown",
    rev_list_is_genomic: bool = False,
) -> dict:
    """Run seqkit amplicon with a pool of fwd/rev combinations via -p TSV.
    Returns a dict of seq_id -> amplicon_seq (UNION across all combinations).

    IMPORTANT: seqkit -p TSV with --bed outputs:
      col0=seqID, col1=start, col2=end, col3=primer_name, col4=score, col5=strand, col6=amplicon_seq
    The primer_name (col3) encodes which combination matched, not the seq_id.
    seq_id is ALWAYS col0.

    rev_list_is_genomic: If True, rev_list sequences are already in genomic/RC orientation
      (i.e., they are the RC of the 5'->3' primer). seqkit amplicon expects the REVERSE PRIMER
      in 5'->3' format and performs RC internally. If rev_list_is_genomic=True, we must RC
      each rev sequence back to 5'->3' before writing to the TSV.
    """
    # Build cross-product TSV: all fwd × rev combinations
    combinations = list(itertools.product(fwd_list, rev_list))
    if not combinations:
        return {}

    tmp_tsv = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False, encoding="utf-8") as tmp:
            for c_idx, (f_comb, r_comb) in enumerate(combinations):
                # FIX: rev_possible_seqs from m08_degenerate_rescue are stored in GENOMIC orientation
                # (RC of the 5'->3' primer) because back-mapping detected a reverse primer.
                # seqkit amplicon -p TSV expects the reverse primer in 5'->3' format and will RC
                # internally. So we must RC the genomic-orientation sequences back to 5'->3'.
                if rev_list_is_genomic:
                    r_comb = reverse_complement(r_comb)
                tmp.write(f"comb_{c_idx}\t{f_comb}\t{r_comb}\n")
            tmp_tsv = tmp.name

        cmd = [
            "seqkit", "amplicon",
            "-p", tmp_tsv,
            "-m", str(max_mismatch),
            "-j", str(threads),
            "--bed",
            str(fasta_path)
        ]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            print(f"[WARN] seqkit pool returned code {res.returncode} for {cand_id}: {res.stderr[:500]}",
                  file=sys.stderr)
            return {}

        if res.stderr.strip():
            # Print first 3 lines of stderr for debugging (seqkit logs progress there)
            for ln in res.stderr.strip().splitlines()[:3]:
                print(f"  [seqkit] {ln}", file=sys.stderr)

        # UNION: collect ALL seq_ids across ALL combination matches
        # BED7 format: seqID start end primerName score strand ampSeq
        matches: dict = {}
        for line in res.stdout.splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 7:
                sid = parts[0]       # Column 1: sequence ID (always correct)
                amp  = parts[6]      # Column 7: amplicon sequence
                # Keep shortest amplicon per sequence (best match across combinations)
                if sid not in matches or len(amp) < len(matches[sid]):
                    matches[sid] = amp
            elif len(parts) == 6:
                # Some seqkit versions omit amp_seq in --bed; fall back without amp_seq
                sid = parts[0]
                if sid not in matches:
                    matches[sid] = ""  # Mark as amplified without stored seq
        return matches

    except Exception as e:
        print(f"[ERROR] seqkit pool exception for {cand_id}: {e}", file=sys.stderr)
        return {}
    finally:
        if tmp_tsv and os.path.exists(tmp_tsv):
            os.remove(tmp_tsv)


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



def _evaluate_candidate(idx_cand, args, seq_mappings, actual_candidates, seqkit_threads):
    """Worker function: evaluate one candidate. Returns (cand_summary, result_row) or None."""
    idx, cand = idx_cand
    cand_id = cand.get("cand_id", f"cand_{idx}")

    # ── Resolve primer sequences ─────────────────────────────────────────
    # fwd_seq / rev_seq = majority (ATGC) sequence used for BLAST / display
    # fwd_degenerate_seq = IUPAC sequence for RAA binding assessment
    # fwd_possible_seqs  = all expanded variants for seqkit union coverage
    fwd   = cand.get("fwd_seq") or cand.get("forward_primer_seq") or cand.get("fwd") or ""
    rev   = cand.get("rev_seq") or cand.get("reverse_primer_seq") or cand.get("rev") or ""
    probe = cand.get("probe_seq") or cand.get("probe") or ""

    # Degenerate IUPAC sequences (fall back to majority if not available)
    fwd_deg   = cand.get("fwd_degenerate_seq") or fwd
    rev_deg   = cand.get("rev_degenerate_seq") or rev
    probe_deg = cand.get("probe_degenerate_seq") or probe

    # Expanded variant lists from m08 degenerate rescue (may be empty for non-degenerate)
    fwd_possible   = cand.get("fwd_possible_seqs") or []
    rev_possible   = cand.get("rev_possible_seqs") or []
    probe_possible = cand.get("probe_possible_seqs") or []

    region = cand.get("region", "unknown")
    gene   = cand.get("gene", region.split("_region_")[0] if "_region_" in region else region)
    rank   = cand.get("rank", idx + 1)

    if not fwd or not rev:
        print(f"[SKIP] {cand_id}: missing fwd_seq or rev_seq.")
        return None

    # Determine if this candidate has degenerate possible_sequences
    has_deg_pool = bool(fwd_possible) and bool(rev_possible)

    # ── ORIENTATION NOTE for degenerate track ────────────────────────────
    # rev_possible_seqs from m08_degenerate_rescue are stored in GENOMIC ORIENTATION
    # (i.e., they are the RC of the 5'->3' reverse primer). This is because m08
    # detects that the reverse primer is reverse-oriented and applies RC during rescue.
    #
    # For seqkit -p TSV: rev primer must be 5'->3' → run_seqkit_amplicon_pool will RC them
    #   via rev_list_is_genomic=True flag.
    # For RAA assessment: rev_possible (genomic orientation) can be used directly to
    #   match against the + strand amplicon without any additional RC.
    #
    # CRITICAL FIX: always include majority sequence in fwd pool.
    # For rev pool: add RC(rev_seq) as the majority variant in GENOMIC orientation
    # to be consistent with rev_possible_seqs orientation.
    if has_deg_pool:
        if fwd not in fwd_possible:
            fwd_possible = [fwd] + fwd_possible
            print(f"  [POOL] Added majority fwd_seq ({len(fwd)}bp) to fwd pool (MSA length drift)")
        # rev_possible is in GENOMIC orientation; add RC(rev) as the majority genomic variant
        rev_majority_genomic = reverse_complement(rev)
        if rev_majority_genomic not in rev_possible:
            rev_possible = [rev_majority_genomic] + rev_possible
            print(f"  [POOL] Added RC(majority rev_seq) ({len(rev_majority_genomic)}bp) to rev pool (MSA length drift)")

    n_combinations = len(fwd_possible) * len(rev_possible) if has_deg_pool else 1
    deg_label = f"({n_combinations} combos incl. majority)" if has_deg_pool else "(non-degenerate)"
    print(f"[{idx+1}/{actual_candidates}] Simulating PCR for {cand_id} {deg_label} (seqkit -j {seqkit_threads})...")


    # ── seqkit amplicon — union-based coverage (surveillance DB) ──────────
    # Strategy: run seqkit against the main --sequences (surveillance/50k DB).
    # For degenerate: pool includes majority + all IUPAC variants.
    if has_deg_pool:
        print(f"  -> [SURVEILLANCE] seqkit pool {n_combinations} combos (rev in genomic orientation, will RC for seqkit)...")
        matches_max = run_seqkit_amplicon_pool(
            fwd_possible, rev_possible,
            args.max_mismatch, args.sequences,
            seqkit_threads, cand_id,
            rev_list_is_genomic=True,  # rev_possible stored as RC of 5'->3' primer
        )
        print(f"  -> [SURVEILLANCE] Union matched: {len(matches_max)} sequences")
    else:
        try:
            matches_max = run_seqkit_amplicon(fwd, rev, args.max_mismatch, args.sequences, seqkit_threads)
        except Exception as e:
            print(f"[ERROR] seqkit failed for {cand_id}: {e}")
            matches_max = {}

    # ── [OPTIONAL] Design-DB coverage ───────────────────────────────────────
    # If --design-sequences is provided, run a second pass against the smaller
    # design database (M02 filtered sequences). This gives a 'fair' coverage
    # metric: how well the primer covers the dataset it was actually designed for.
    design_matches = {}
    if args.design_sequences:
        print(f"  -> [DESIGN-DB] Running against design sequences...")
        if has_deg_pool:
            design_matches = run_seqkit_amplicon_pool(
                fwd_possible, rev_possible,
                args.max_mismatch, args.design_sequences,
                seqkit_threads, cand_id,
                rev_list_is_genomic=True,  # rev_possible stored as RC of 5'->3' primer
            )
        else:
            try:
                design_matches = run_seqkit_amplicon(fwd, rev, args.max_mismatch, args.design_sequences, seqkit_threads)
            except Exception as e:
                print(f"[ERROR] seqkit design-db failed for {cand_id}: {e}")
        print(f"  -> [DESIGN-DB]  Matched: {len(design_matches)} sequences")

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
            # For non-degenerate: use fwd_deg/rev_deg (same as majority) with IUPAC-aware match.
            # For degenerate: rev_deg is in GENOMIC orientation (RC of 5'->3'); use it directly
            # for matching against + strand amplicon (no additional RC needed).
            if has_deg_pool:
                # rev_deg is genomic (RC of 5'->3'), match directly against + strand
                fwd_mm = min_mismatches_sliding_window(fwd_deg, amp_seq)
                rev_mm = min_mismatches_sliding_window(rev_deg, amp_seq)
            else:
                # rev_deg is 5'->3' format, must RC before matching + strand
                fwd_mm = min_mismatches_sliding_window(fwd_deg, amp_seq)
                rev_mm = min_mismatches_sliding_window(reverse_complement(rev_deg), amp_seq)
            p_mismatch = min(max(fwd_mm, rev_mm), args.max_mismatch)

            # ─ RAA binding assessment — best-of-pool strategy ─
            # For DEGENERATE track: evaluate all variants from fwd_possible and rev_possible.
            #   fwd_possible: 5'->3' format → use directly for RAA binding check.
            #   rev_possible: GENOMIC orientation → use directly against + strand amplicon.
            # For NON-DEGENERATE: use majority sequences with correct RC for rev.
            if has_deg_pool:
                # fwd: best-of-pool (fwd_possible is 5'->3')
                # Cap pool to 8 variants for RAA (don't need all 289 combos)
                _fwd_raa_pool = fwd_possible[:8]
                _rev_raa_pool = rev_possible[:8]
                fwd_raa_ok = False
                fwd_3p, fwd_int = args.raa_3prime_bp, args.raa_max_internal_mm + 1
                for fwd_var in _fwd_raa_pool:
                    ok, t3, ti, _ = (_assess_raa_numpy if _HAS_NUMPY else assess_raa_binding)(
                        fwd_var, amp_seq,
                        critical_3prime_bp=args.raa_3prime_bp,
                        max_3prime_mm=args.raa_max_3prime_mm,
                        max_internal_mm=args.raa_max_internal_mm,
                    )
                    if ok:
                        fwd_raa_ok = True
                        fwd_3p, fwd_int = t3, ti
                        break
                    if t3 < fwd_3p:
                        fwd_3p, fwd_int = t3, ti

                # rev: best-of-pool (rev_possible is GENOMIC = matches + strand directly)
                rev_raa_ok = False
                rev_3p, rev_int = args.raa_3prime_bp, args.raa_max_internal_mm + 1
                for rev_var in _rev_raa_pool:
                    ok, t3, ti, _ = (_assess_raa_numpy if _HAS_NUMPY else assess_raa_binding)(
                        rev_var, amp_seq,
                        critical_3prime_bp=args.raa_3prime_bp,
                        max_3prime_mm=args.raa_max_3prime_mm,
                        max_internal_mm=args.raa_max_internal_mm,
                    )
                    if ok:
                        rev_raa_ok = True
                        rev_3p, rev_int = t3, ti
                        break
                    if t3 < rev_3p:
                        rev_3p, rev_int = t3, ti
            else:
                # Non-degenerate: rev_deg is 5'->3', must RC before RAA assessment vs + strand
                fwd_raa_ok, fwd_3p, fwd_int, _ = assess_raa_binding(
                    fwd_deg, amp_seq,
                    critical_3prime_bp=args.raa_3prime_bp,
                    max_3prime_mm=args.raa_max_3prime_mm,
                    max_internal_mm=args.raa_max_internal_mm,
                )
                rev_raa_ok, rev_3p, rev_int, _ = assess_raa_binding(
                    reverse_complement(rev_deg), amp_seq,
                    critical_3prime_bp=args.raa_3prime_bp,
                    max_3prime_mm=args.raa_max_3prime_mm,
                    max_internal_mm=args.raa_max_internal_mm,
                )

            is_raa_amplified = fwd_raa_ok and rev_raa_ok
            raa_3p_mm  = max(fwd_3p, rev_3p)
            raa_int_mm = max(fwd_int, rev_int)

            # ─ Probe mismatch — use IUPAC degenerate probe sequence ─
            # probe_deg and probe_possible are always in 5'->3' format.
            if probe_deg:
                if has_deg_pool and probe_possible:
                    # Best-of-pool for probe (5'->3' format)
                    prb_mismatch = min(
                        min_mismatches_sliding_window(pv, amp_seq) for pv in probe_possible
                    )
                    probe_raa_ok = False
                    for pv in probe_possible:
                        ok, _, _, _ = assess_raa_binding(
                            pv, amp_seq,
                            critical_3prime_bp=args.raa_3prime_bp,
                            max_3prime_mm=args.raa_max_3prime_mm,
                            max_internal_mm=args.raa_max_internal_mm,
                        )
                        if ok:
                            probe_raa_ok = True
                            break
                else:
                    prb_mismatch = min_mismatches_sliding_window(probe_deg, amp_seq)
                    probe_raa_ok, _, _, _ = assess_raa_binding(
                        probe_deg, amp_seq,
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

    print(f"  [{cand_id}] PCR -> Amplified: {amplified}/{total} ({amp_pct}%), Detected: {detected}/{total} ({det_pct}%)")
    print(f"  [{cand_id}] RAA -> Amplified: {raa_amplified}/{total} ({raa_amp_pct}%), Detected: {raa_detected}/{total} ({raa_det_pct}%)")

    # ── Design-DB coverage WITH subtype/year/country breakdown ─────────────
    # seq_mappings already contains metadata for all sequences (including design DB seqs).
    # We filter by seq_ids in design_matches to build breakdown tables.
    design_total = 0
    design_amplified = 0
    design_amp_pct = 0.0
    design_raa_amplified = 0
    design_raa_detected  = 0
    design_raa_amp_pct   = 0.0
    design_raa_det_pct   = 0.0
    design_db_breakdowns = {"subtype": [], "year": [], "country": []}
    if design_matches:
        # Count total design FASTA headers for denominator
        try:
            design_total = sum(1 for l in open(args.design_sequences, encoding='utf-8', errors='ignore') if l.startswith('>'))
        except Exception:
            design_total = 0
        design_amplified = len(design_matches)
        design_amp_pct = round(design_amplified / design_total * 100, 2) if design_total else 0.0
        print(f"  [{cand_id}] DESIGN-DB PCR -> Amplified: {design_amplified}/{design_total} ({design_amp_pct}%)")

        # Build per-row records for design DB using seq_mappings metadata.
        # design_matches keys = amplified seq_ids; join with metadata for breakdowns.
        design_matched_ids = set(design_matches.keys())
        try:
            fasta_design_ids = set()
            with open(args.design_sequences, encoding='utf-8', errors='ignore') as _fh:
                for _ln in _fh:
                    if _ln.startswith('>'):
                        fasta_design_ids.add(_ln[1:].split()[0].strip())
        except Exception:
            fasta_design_ids = design_matched_ids.copy()

        meta_lookup = {s["seq_id"]: s for s in seq_mappings}
        design_rows = []
        for sid in fasta_design_ids:
            meta = meta_lookup.get(sid)
            is_amp = sid in design_matched_ids

            # ── RAA assessment for design DB sequences ──────────────────────
            # Run same sliding-window mismatch check as surveillance DB loop.
            # Only run RAA on amplified sequences (non-amplified → False by definition).
            is_raa_amp_des = False
            is_raa_det_des = False
            if is_amp:
                amp_seq_des = design_matches[sid]  # amplified sequence string from seqkit
                if has_deg_pool:
                    _fwd_raa = fwd_possible[:8]
                    _rev_raa = rev_possible[:8]
                    _fwd_ok = any(
                        (_assess_raa_numpy if _HAS_NUMPY else assess_raa_binding)(
                            fv, amp_seq_des,
                            critical_3prime_bp=args.raa_3prime_bp,
                            max_3prime_mm=args.raa_max_3prime_mm,
                            max_internal_mm=args.raa_max_internal_mm,
                        )[0]
                        for fv in _fwd_raa
                    )
                    _rev_ok = any(
                        (_assess_raa_numpy if _HAS_NUMPY else assess_raa_binding)(
                            rv, amp_seq_des,
                            critical_3prime_bp=args.raa_3prime_bp,
                            max_3prime_mm=args.raa_max_3prime_mm,
                            max_internal_mm=args.raa_max_internal_mm,
                        )[0]
                        for rv in _rev_raa
                    )
                    is_raa_amp_des = _fwd_ok and _rev_ok
                else:
                    _fwd_ok, _, _, _ = assess_raa_binding(
                        fwd_deg, amp_seq_des,
                        critical_3prime_bp=args.raa_3prime_bp,
                        max_3prime_mm=args.raa_max_3prime_mm,
                        max_internal_mm=args.raa_max_internal_mm,
                    )
                    _rev_ok, _, _, _ = assess_raa_binding(
                        reverse_complement(rev_deg), amp_seq_des,
                        critical_3prime_bp=args.raa_3prime_bp,
                        max_3prime_mm=args.raa_max_3prime_mm,
                        max_internal_mm=args.raa_max_internal_mm,
                    )
                    is_raa_amp_des = _fwd_ok and _rev_ok
                # Probe RAA check (only if primer binding passes)
                if is_raa_amp_des and probe:
                    if has_deg_pool and probe_possible:
                        is_raa_det_des = any(
                            assess_raa_binding(
                                pv, amp_seq_des,
                                critical_3prime_bp=args.raa_3prime_bp,
                                max_3prime_mm=args.raa_max_3prime_mm,
                                max_internal_mm=args.raa_max_internal_mm,
                            )[0]
                            for pv in probe_possible
                        )
                    else:
                        is_raa_det_des, _, _, _ = assess_raa_binding(
                            probe_deg, amp_seq_des,
                            critical_3prime_bp=args.raa_3prime_bp,
                            max_3prime_mm=args.raa_max_3prime_mm,
                            max_internal_mm=args.raa_max_internal_mm,
                        )
                elif is_raa_amp_des and not probe:
                    # No probe → detection == amplification
                    is_raa_det_des = True

            design_rows.append({
                "seq_id":           sid,
                "subtype":          meta["subtype"]   if meta else "Unknown",
                "year":             str(meta["year"]) if meta else "Unknown",
                "country":          meta["country"]   if meta else "Unknown",
                "is_amplified":     is_amp,
                "is_raa_amplified": is_raa_amp_des,
                "is_raa_detected":  is_raa_det_des,
            })

        # Compute design DB RAA aggregate metrics
        design_raa_amplified = 0
        design_raa_detected  = 0
        design_raa_amp_pct   = 0.0
        design_raa_det_pct   = 0.0
        if design_rows:
            import pandas as _pd_des
            df_des = _pd_des.DataFrame(design_rows)
            design_raa_amplified = int(df_des["is_raa_amplified"].sum())
            design_raa_detected  = int(df_des["is_raa_detected"].sum())
            design_raa_amp_pct   = round(design_raa_amplified / design_total * 100, 2) if design_total else 0.0
            design_raa_det_pct   = round(design_raa_detected  / design_total * 100, 2) if design_total else 0.0
            print(f"  [{cand_id}] DESIGN-DB RAA -> Amplified: {design_raa_amplified}/{design_total} ({design_raa_amp_pct}%), Detected: {design_raa_detected}/{design_total} ({design_raa_det_pct}%)")
            for _col in ["subtype", "year", "country"]:
                _grp = df_des.groupby(_col).agg(
                    total=("seq_id", "count"),
                    pcr_amplified=("is_amplified", "sum"),
                    raa_amplified=("is_raa_amplified", "sum"),
                    raa_detected=("is_raa_detected", "sum"),
                ).reset_index()
                _grp["pcr_amplified"]     = _grp["pcr_amplified"].astype(int)
                _grp["raa_amplified"]     = _grp["raa_amplified"].astype(int)
                _grp["raa_detected"]      = _grp["raa_detected"].astype(int)
                _grp["pcr_amplified_pct"] = (_grp["pcr_amplified"] / _grp["total"] * 100).round(2)
                _grp["raa_amplified_pct"] = (_grp["raa_amplified"] / _grp["total"] * 100).round(2)
                _grp["raa_detected_pct"]  = (_grp["raa_detected"]  / _grp["total"] * 100).round(2)
                _grp = _grp.sort_values("year" if _col == "year" else "total",
                                        ascending=(_col == "year"))
                design_db_breakdowns[_col] = _grp.to_dict(orient="records")

    # ── Breakdown by subtype / year / country ─────────────────────────────
    # For subtype breakdown: exclude subtypes with < min_subtype_seqs total sequences
    # (they cannot be statistically generalised). Aggregate them as 'Other (rare subtypes)'.
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

        if col == "subtype":
            # Split into major (>= threshold) and rare (< threshold) subtypes
            min_seqs = args.min_subtype_seqs
            major_mask = grp["total"] >= min_seqs
            rare_rows  = grp[~major_mask]
            grp = grp[major_mask].copy()
            if not rare_rows.empty:
                other_row = {
                    "subtype": f"Other (rare, n<{min_seqs})",
                    "total":          int(rare_rows["total"].sum()),
                    "pcr_amplified":  int(rare_rows["pcr_amplified"].sum()),
                    "pcr_detected":   int(rare_rows["pcr_detected"].sum()),
                    "raa_amplified":  int(rare_rows["raa_amplified"].sum()),
                    "raa_detected":   int(rare_rows["raa_detected"].sum()),
                }
                t = other_row["total"]
                other_row["pcr_amplified_pct"] = round(other_row["pcr_amplified"] / t * 100, 2) if t else 0.0
                other_row["pcr_detected_pct"]  = round(other_row["pcr_detected"]  / t * 100, 2) if t else 0.0
                other_row["raa_amplified_pct"] = round(other_row["raa_amplified"] / t * 100, 2) if t else 0.0
                other_row["raa_detected_pct"]  = round(other_row["raa_detected"]  / t * 100, 2) if t else 0.0
                grp = pd.concat([grp, pd.DataFrame([other_row])], ignore_index=True)
            grp = grp.sort_values("total", ascending=False)
        elif col == "year":
            grp = grp.sort_values("year", ascending=True)
        else:
            grp = grp.sort_values("total", ascending=False)

        breakdowns[col] = grp.to_dict(orient="records")

    cand_summary = {
        "cand_id": cand_id, "region": region, "gene": gene, "rank": rank,
        "fwd_seq": fwd, "rev_seq": rev, "probe_seq": probe,
        # Degenerate IUPAC sequences (for display/RAA context)
        "fwd_degenerate_seq": fwd_deg if fwd_deg != fwd else None,
        "rev_degenerate_seq": rev_deg if rev_deg != rev else None,
        "n_combinations": n_combinations,
        # Surveillance DB (full 50k EpiFlu)
        "total_seqs": total,
        "pcr_amplified_seqs": int(amplified), "pcr_amplified_pct": amp_pct,
        "pcr_detected_seqs":  int(detected),  "pcr_detected_pct":  det_pct,
        # RAA (surveillance)
        "raa_amplified_seqs": int(raa_amplified), "raa_amplified_pct": raa_amp_pct,
        "raa_detected_seqs":  int(raa_detected),  "raa_detected_pct":  raa_det_pct,
        # Design DB (M02 filtered — primer was designed from these)
        "design_db_total": design_total,
        "design_db_amplified": design_amplified,
        "design_db_amplified_pct": design_amp_pct,
        # Design DB RAA (same sliding-window binding model, applied to design DB amplicons)
        "design_db_raa_amplified": design_raa_amplified if design_matches else 0,
        "design_db_raa_detected":  design_raa_detected  if design_matches else 0,
        "design_db_raa_amp_pct":   design_raa_amp_pct   if design_matches else 0.0,
        "design_db_raa_det_pct":   design_raa_det_pct   if design_matches else 0.0,
        "design_db_breakdowns": design_db_breakdowns,
        "raa_params": {
            "critical_3prime_bp":  args.raa_3prime_bp,
            "max_3prime_mm":       args.raa_max_3prime_mm,
            "max_internal_mm":     args.raa_max_internal_mm,
        },
        "breakdowns": breakdowns,
    }
    result_row = {
        "cand_id": cand_id, "region": region, "gene": gene, "rank": rank,
        "fwd_seq": fwd, "rev_seq": rev, "probe_seq": probe,
        # Surveillance DB
        "total": total,
        "pcr_amplified": int(amplified), "pcr_detected": int(detected),
        "pcr_amplified_pct": amp_pct,    "pcr_detected_pct": det_pct,
        # RAA (surveillance)
        "raa_amplified": int(raa_amplified), "raa_detected": int(raa_detected),
        "raa_amplified_pct": raa_amp_pct,    "raa_detected_pct": raa_det_pct,
        # Design DB
        "design_db_total": design_total,
        "design_db_amplified": design_amplified,
        "design_db_amplified_pct": design_amp_pct,
        # Design DB RAA
        "design_db_raa_amplified": design_raa_amplified if design_matches else 0,
        "design_db_raa_detected":  design_raa_detected  if design_matches else 0,
        "design_db_raa_amp_pct":   design_raa_amp_pct   if design_matches else 0.0,
        "design_db_raa_det_pct":   design_raa_det_pct   if design_matches else 0.0,
    }
    return (cand_summary, result_row)


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
    parser.add_argument("--min-subtype-seqs", type=int, default=30,
                        help="Minimum total sequences a subtype must have to appear in breakdown report (default: 30). "
                             "Subtypes with fewer sequences are aggregated into 'Other (rare subtypes)' group.")
    parser.add_argument("--design-sequences", default=None,
                        help="[OPTIONAL] Path to the FASTA used for primer design (e.g. M02 filtered sequences). "
                             "If provided, coverage will be computed against BOTH this database (design coverage) "
                             "AND the main --sequences database (surveillance coverage), giving a dual-layer "
                             "view of primer performance.")
    args = parser.parse_args()
    
    # ── Thread allocation strategy ──────────────────────────────────────────
    # Goal: run ALL candidates in parallel, giving each seqkit as many threads as possible.
    # For degenerate candidates with hundreds of primer combinations, seqkit benefits
    # greatly from more -j goroutines, so we balance between candidate count and threads.
    total_cpus = args.threads or os.cpu_count() or 8
    # Pre-count candidates to tune allocation (we haven't loaded them yet, so estimate)
    # Conservative estimate: max_candidates = num_candidates × (number of m09 dirs)
    est_candidates = max(1, args.num_candidates)  # refine below after load
    # Give seqkit at least 2 threads, and run as many candidates in parallel as possible.
    # We'll recompute below once actual candidate count is known.
    seqkit_threads = max(2, min(16, total_cpus // max(1, est_candidates)))
    worker_threads = max(1, total_cpus // seqkit_threads)
    print(f"[M10] Thread pre-config: {total_cpus} total CPUs -> est {worker_threads} parallel candidates x {seqkit_threads} seqkit threads each (will recompute after loading)")

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
    
    # Map GISAID/NCBI/EpiFlu fields to standardized pipeline columns
    if 'collection_date' in df_meta.columns and 'date' not in df_meta.columns:
        df_meta = df_meta.rename(columns={'collection_date': 'date'})
        print("  Mapped 'collection_date' -> 'date' for GISAID/EpiFlu compatibility.")
        
    if 'source_country' in df_meta.columns and 'country' not in df_meta.columns:
        df_meta = df_meta.rename(columns={'source_country': 'country'})
        print("  Mapped 'source_country' -> 'country' for GISAID/EpiFlu compatibility.")
    elif 'location' in df_meta.columns and 'country' not in df_meta.columns:
        # Extract country from GISAID location (Continent / Country / Division / City)
        def parse_gisaid_country(loc_val):
            if pd.isna(loc_val) or not isinstance(loc_val, str) or not loc_val.strip():
                return 'Unknown'
            parts = [p.strip() for p in loc_val.split('/')]
            return parts[1] if len(parts) > 1 else 'Unknown'
        df_meta['country'] = df_meta['location'].apply(parse_gisaid_country)
        print("  Extracted 'country' from 'location' for GISAID compatibility.")

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
        for possible_strain in ['isolate_id', 'isolate_name', 'accession']:
            if possible_strain in df_meta.columns:
                df_meta = df_meta.rename(columns={possible_strain: 'strain'})
                print(f"  Mapped '{possible_strain}' -> 'strain' for compatibility.")
                break
        else:
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
        seq_id = header.split()[0].strip() if header.strip() else ""
        if meta_row:
            seq_mappings.append({
                "header": header,
                "seq_id": seq_id,
                "subtype": meta_row.get("subtype", "Unknown"),
                "year": parse_year(meta_row.get("date", "Unknown")),
                "country": meta_row.get("country", "Unknown")
            })
        else:
            # Fallback for sequences without metadata
            seq_mappings.append({
                "header": header,
                "seq_id": seq_id,
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

    # ── Recompute thread allocation now that we know the actual candidate count ───
    actual_candidates = len(all_candidates)
    # Each candidate runs one seqkit amplicon call (potentially with hundreds of combinations
    # via -p TSV for degenerate). We want to run all candidates in parallel.
    # seqkit itself parallelises internally with -j goroutines per call.
    # Optimal: worker_threads = actual_candidates, seqkit_threads = total / worker_threads
    # ProcessPoolExecutor: each worker = 1 separate Python process (bypasses GIL).
    # Moving _evaluate_candidate to module level allows pickling/serialization.
    # Cap at 16 to avoid memory oversubscription; seqkit handles internal threading.
    worker_threads = min(actual_candidates, total_cpus, 16)
    seqkit_threads = max(2, total_cpus // worker_threads)
    print(f"[M10] Thread final config: {total_cpus} total CPUs -> {worker_threads} parallel candidates x {seqkit_threads} seqkit threads/candidate")
    # Strategy: run seqkit ONCE per candidate (with max_mismatch), then use sliding-window
    # mismatch detection to determine exact mismatch level per-sequence. This cuts seqkit
    # calls from 4× down to 1× per candidate while still reporting per-mismatch breakdown.


    # 4. PCR Simulation & Probe Binding Check — fully parallelized
    # Strategy for DEGENERATE primers:
    #   - candidate carries fwd_possible_seqs / rev_possible_seqs (all ATGC variants from m08)
    #   - build cross-product TSV and run seqkit with -p (single call, all combinations)
    #   - UNION all matched seq_ids across ALL combinations → true degenerate coverage
    #   - RAA assessment uses the IUPAC degenerate sequence (IUPAC-aware match_base)
    # Strategy for NON-DEGENERATE primers:
    #   - single seqkit call with -F / -R as before

    import concurrent.futures


    # ── Run candidates in parallel ────────────────────────────────────────────
    results_list = []
    summary_data = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=worker_threads) as pool:
        futures = {pool.submit(_evaluate_candidate, (i, c), args, seq_mappings, actual_candidates, seqkit_threads): i
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

def generate_heatmap_dashboard(summary_data: list, output_path):
    """Generate a premium interactive dark-mode HTML coverage dashboard with full PCR and RAA data."""
    import json as _json
    data_json = _json.dumps(summary_data, default=str)

    html_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>\U0001f9ec In-silico PCR/RAA Coverage Dashboard</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
    <style>
        :root {
            --bg:#080f1e;--surface:#0f1e35;--surface2:#162036;
            --card-bg:rgba(22,33,62,0.85);--border:rgba(255,255,255,0.07);
            --accent-pcr:#10b981;--accent-pcr-dim:rgba(16,185,129,0.15);
            --accent-raa:#6366f1;--accent-raa-dim:rgba(99,102,241,0.15);
            --accent-design:#a78bfa;--text-1:#f0f6ff;--text-2:#8b9dc0;
        }
        *{box-sizing:border-box;margin:0;padding:0;}
        body{font-family:'Outfit',sans-serif;background:var(--bg);color:var(--text-1);padding:2rem;}
        .dashboard{max-width:1500px;margin:0 auto;}
        header{margin-bottom:2rem;padding-bottom:1.5rem;border-bottom:1px solid var(--border);}
        .header-flex{display:flex;align-items:flex-start;justify-content:space-between;flex-wrap:wrap;gap:1rem;}
        h1{font-size:2rem;font-weight:700;background:linear-gradient(135deg,#34d399 0%,#6366f1 100%);-webkit-background-clip:text;-webkit-text-fill-color:transparent;margin-bottom:0.4rem;}
        .subtitle{color:var(--text-2);font-size:0.95rem;line-height:1.6;}
        .toggle-container{display:inline-flex;background:var(--card-bg);border:1px solid var(--border);border-radius:14px;padding:4px;box-shadow:0 8px 20px rgba(0,0,0,0.3);}
        .toggle-btn{background:none;border:none;color:var(--text-2);padding:10px 22px;border-radius:10px;font-weight:600;cursor:pointer;transition:all 0.25s ease;font-family:inherit;font-size:0.9rem;display:flex;align-items:center;gap:8px;white-space:nowrap;}
        .toggle-btn.active-global{background:var(--accent-pcr);color:#0a1628;box-shadow:0 4px 16px rgba(16,185,129,0.3);}
        .toggle-btn.active-design{background:var(--accent-design);color:#0a1628;box-shadow:0 4px 16px rgba(167,139,250,0.3);}
        .toggle-btn:hover:not(.active-global):not(.active-design){color:var(--text-1);background:rgba(255,255,255,0.05);}
        .controls{display:flex;align-items:center;gap:1rem;margin-bottom:2rem;background:var(--card-bg);padding:1rem 1.5rem;border-radius:14px;border:1px solid var(--border);flex-wrap:wrap;}
        .controls label{font-weight:600;color:var(--text-2);font-size:0.9rem;}
        select{background:var(--surface2);color:var(--text-1);border:1px solid var(--border);padding:0.5rem 1rem;border-radius:8px;font-family:inherit;font-size:0.95rem;cursor:pointer;outline:none;flex:1;min-width:280px;}
        select:focus{border-color:var(--accent-pcr);}
        .dual-metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:1.5rem;margin-bottom:2rem;}
        .metric-card{background:var(--card-bg);border:1px solid var(--border);border-radius:16px;padding:1.5rem;backdrop-filter:blur(10px);transition:border-color 0.3s,transform 0.3s;}
        .metric-card:hover{border-color:rgba(16,185,129,0.35);transform:translateY(-3px);}
        .metric-label{font-size:0.8rem;color:var(--text-2);text-transform:uppercase;letter-spacing:0.08em;margin-bottom:0.8rem;}
        .metric-dual{display:flex;align-items:baseline;gap:1.5rem;}
        .metric-block{flex:1;}
        .metric-mode-label{font-size:0.72rem;font-weight:600;text-transform:uppercase;letter-spacing:0.06em;margin-bottom:0.2rem;}
        .metric-value{font-size:1.7rem;font-weight:700;}
        .pcr-color{color:var(--accent-pcr);}
        .raa-color{color:var(--accent-raa);}
        .divider-v{width:1px;background:var(--border);height:50px;align-self:center;}
        .prog-wrap{height:4px;background:rgba(255,255,255,0.05);border-radius:9px;margin-top:0.6rem;overflow:hidden;}
        .prog-bar{height:100%;border-radius:9px;transition:width 0.6s ease;}
        .seq-grid{display:grid;grid-template-columns:1fr 1fr 1fr;gap:1rem;margin-bottom:2rem;}
        @media(max-width:1100px){.seq-grid{grid-template-columns:1fr 1fr;}}
        @media(max-width:700px){.seq-grid{grid-template-columns:1fr;}}
        .seq-card{background:var(--card-bg);border:1px solid var(--border);border-radius:14px;padding:1.2rem;}
        .seq-title{font-size:0.78rem;font-weight:600;text-transform:uppercase;letter-spacing:0.07em;color:var(--text-2);margin-bottom:0.6rem;}
        .seq-mono{font-family:'JetBrains Mono',monospace;font-size:0.82rem;background:#060d1a;padding:0.6rem 0.8rem;border-radius:8px;border:1px solid rgba(255,255,255,0.04);word-break:break-all;line-height:1.6;color:#a5f3c3;}
        .seq-deg-note{font-size:0.72rem;color:var(--text-2);margin-top:0.4rem;}
        .design-banner{background:linear-gradient(135deg,rgba(139,92,246,0.12),rgba(109,40,217,0.08));border:1px solid rgba(139,92,246,0.3);border-radius:16px;padding:1.2rem 1.5rem;margin-bottom:2rem;}
        .design-banner-header{font-size:0.75rem;font-weight:700;text-transform:uppercase;letter-spacing:0.08em;color:#a78bfa;margin-bottom:0.8rem;}
        .design-content{display:flex;gap:2.5rem;flex-wrap:wrap;align-items:center;}
        .design-stat-label{font-size:0.72rem;color:rgba(167,139,250,0.65);margin-bottom:0.2rem;}
        .design-stat-value{font-size:1.5rem;font-weight:700;color:#c4b5fd;}
        .design-bar-wrap{flex:1;min-width:120px;}
        .design-bar-track{background:rgba(255,255,255,0.05);border-radius:9px;height:8px;overflow:hidden;}
        .design-bar-fill{height:100%;border-radius:9px;transition:width 0.6s ease;}
        .main-layout{display:grid;grid-template-columns:1fr 1fr;gap:1.5rem;}
        @media(max-width:900px){.main-layout{grid-template-columns:1fr;}.full-width{grid-column:span 1!important;}}
        .full-width{grid-column:span 2;}
        .card{background:var(--card-bg);border:1px solid var(--border);border-radius:16px;padding:1.5rem;backdrop-filter:blur(10px);}
        .card-title{font-size:1.05rem;font-weight:600;margin-bottom:1.2rem;padding-bottom:0.75rem;border-bottom:1px solid var(--border);display:flex;justify-content:space-between;align-items:center;}
        .badge{padding:0.2rem 0.6rem;border-radius:6px;font-size:0.75rem;font-weight:600;}
        .badge-pcr{background:var(--accent-pcr-dim);color:var(--accent-pcr);border:1px solid rgba(16,185,129,0.3);}
        .badge-raa{background:var(--accent-raa-dim);color:var(--accent-raa);border:1px solid rgba(99,102,241,0.3);}
        table{width:100%;border-collapse:collapse;font-size:0.88rem;}
        th,td{padding:0.6rem 0.75rem;text-align:left;border-bottom:1px solid rgba(255,255,255,0.04);}
        th{color:var(--text-2);font-size:0.75rem;font-weight:700;text-transform:uppercase;letter-spacing:0.05em;}
        tr:hover td{background:rgba(255,255,255,0.02);}
        .cov-cell{font-weight:700;border-radius:6px;padding:0.35rem 0.55rem;display:inline-block;transition:transform 0.15s;font-size:0.85rem;}
        .cov-cell:hover{transform:scale(1.08);}
        .mini-bar{width:100%;background:rgba(255,255,255,0.06);border-radius:4px;height:4px;margin-top:3px;overflow:hidden;}
        .mini-bar-fill{height:100%;border-radius:4px;}
    </style>
</head>
<body>
<div class="dashboard">
    <header>
        <div class="header-flex">
            <div>
                <h1>\U0001f9ec PCR &amp; RAA Coverage Dashboard</h1>
                <div class="subtitle">
                    Dual-mode in-silico validation: <span style="color:var(--accent-pcr);font-weight:600;">PCR</span> (seqkit amplicon)
                    &amp; <span style="color:var(--accent-raa);font-weight:600;">RAA/RPA</span> (recombinase binding model).
                </div>
            </div>
            <div class="toggle-container">
                <button class="toggle-btn active-global" id="btn-global" onclick="setDatasetMode('global')">\U0001f30d Global Surveillance DB</button>
                <button class="toggle-btn" id="btn-design" onclick="setDatasetMode('design')">\U0001f3af Design DB Only</button>
            </div>
        </div>
    </header>

    <div class="controls">
        <label for="candidateSelect">Select Candidate:</label>
        <select id="candidateSelect" onchange="updateDashboard(this.value)"></select>
    </div>

    <div class="dual-metrics">
        <div class="metric-card">
            <div class="metric-label">Candidate &middot; Region &middot; Rank</div>
            <div id="statRegion" style="font-size:1.2rem;font-weight:700;">&mdash;</div>
        </div>
        <div class="metric-card">
            <div class="metric-label">Total Sequences (DB)</div>
            <div id="statTotal" style="font-size:1.7rem;font-weight:700;">&mdash;</div>
        </div>
        <div class="metric-card">
            <div class="metric-label">\U0001f52c Amplification Coverage</div>
            <div class="metric-dual">
                <div class="metric-block">
                    <div class="metric-mode-label pcr-color">PCR</div>
                    <div class="metric-value pcr-color" id="pcrAmpPct">&mdash;%</div>
                    <div class="prog-wrap"><div class="prog-bar" id="pcrAmpBar" style="background:var(--accent-pcr);width:0%;"></div></div>
                </div>
                <div class="divider-v"></div>
                <div class="metric-block">
                    <div class="metric-mode-label raa-color" id="raaAmpLabel">RAA</div>
                    <div class="metric-value raa-color" id="raaAmpPct">&mdash;</div>
                    <div id="raaAmpNote" style="font-size:0.62rem;color:var(--text-2);margin-top:1px;display:none;">(Global Surveillance DB)</div>
                    <div class="prog-wrap"><div class="prog-bar" id="raaAmpBar" style="background:var(--accent-raa);width:0%;"></div></div>
                </div>
            </div>
        </div>
        <div class="metric-card">
            <div class="metric-label">\U0001f3af Detection Coverage (+Probe)</div>
            <div class="metric-dual">
                <div class="metric-block">
                    <div class="metric-mode-label pcr-color">PCR</div>
                    <div class="metric-value pcr-color" id="pcrDetPct">&mdash;%</div>
                    <div class="prog-wrap"><div class="prog-bar" id="pcrDetBar" style="background:var(--accent-pcr);width:0%;"></div></div>
                </div>
                <div class="divider-v"></div>
                <div class="metric-block">
                    <div class="metric-mode-label raa-color" id="raaDetLabel">RAA</div>
                    <div class="metric-value raa-color" id="raaDetPct">&mdash;</div>
                    <div id="raaDetNote" style="font-size:0.62rem;color:var(--text-2);margin-top:1px;display:none;">(Global Surveillance DB)</div>
                    <div class="prog-wrap"><div class="prog-bar" id="raaDetBar" style="background:var(--accent-raa);width:0%;"></div></div>
                </div>
            </div>
        </div>
    </div>

    <div class="design-banner" id="designBanner" style="display:none;">
        <div class="design-banner-header">\U0001f3af Design DB Coverage (sekuens yang dipakai desain primer &mdash; M02 filtered)</div>
        <div class="design-content">
            <div><div class="design-stat-label">Total Sequences</div><div class="design-stat-value" id="designTotal">&mdash;</div></div>
            <div><div class="design-stat-label">PCR Amplified</div><div class="design-stat-value" id="designAmplified">&mdash;</div></div>
            <div><div class="design-stat-label">PCR Coverage</div><div class="design-stat-value" id="designPct">&mdash;%</div></div>
            <div><div class="design-stat-label">RAA Amplified</div><div class="design-stat-value" id="designRaaAmplified">&mdash;</div></div>
            <div><div class="design-stat-label">RAA Coverage</div><div class="design-stat-value" id="designRaaPct">&mdash;%</div></div>
            <div class="design-bar-wrap">
                <div class="design-stat-label">PCR Coverage Bar</div>
                <div class="design-bar-track"><div class="design-bar-fill" id="designBar" style="width:0%;background:linear-gradient(90deg,#7c3aed,#a78bfa);"></div></div>
            </div>
        </div>
    </div>

    <div class="seq-grid">
        <div class="seq-card">
            <div class="seq-title">&#x2B06; Forward Primer (5&#x27;&#x2192;3&#x27;) &mdash; Majority</div>
            <div class="seq-mono" id="seqFwd">&mdash;</div>
            <div class="seq-deg-note" id="seqFwdDeg" style="display:none;"></div>
        </div>
        <div class="seq-card">
            <div class="seq-title">&#x2B07; Reverse Primer (5&#x27;&#x2192;3&#x27;) &mdash; Majority</div>
            <div class="seq-mono" id="seqRev">&mdash;</div>
            <div class="seq-deg-note" id="seqRevDeg" style="display:none;"></div>
        </div>
        <div class="seq-card">
            <div class="seq-title">&#x1F534; Probe Sequence (5&#x27;&#x2192;3&#x27;)</div>
            <div class="seq-mono" id="seqProbe">&mdash;</div>
            <div class="seq-deg-note" id="seqProbeDeg" style="display:none;"></div>
        </div>
    </div>

    <div id="globalTables">
        <div class="main-layout">
            <div class="card">
                <div class="card-title">
                    Subtype Coverage Distribution
                    <div style="display:flex;gap:0.5rem;">
                        <span class="badge badge-pcr">PCR</span>
                        <span class="badge badge-raa">RAA</span>
                    </div>
                </div>
                <table id="subtypeTable">
                    <thead><tr>
                        <th>Subtype</th><th>N</th>
                        <th style="color:var(--accent-pcr);">PCR Amp</th>
                        <th style="color:var(--accent-pcr);">PCR Det</th>
                        <th style="color:var(--accent-raa);">RAA Amp</th>
                        <th style="color:var(--accent-raa);">RAA Det</th>
                    </tr></thead>
                    <tbody></tbody>
                </table>
            </div>
            <div class="card">
                <div class="card-title">
                    Collection Year Breakdown
                    <span class="badge badge-pcr">Temporal</span>
                </div>
                <table id="yearTable">
                    <thead><tr>
                        <th>Year</th><th>N</th>
                        <th style="color:var(--accent-pcr);">PCR Amp</th>
                        <th style="color:var(--accent-pcr);">PCR Det</th>
                        <th style="color:var(--accent-raa);">RAA Amp</th>
                        <th style="color:var(--accent-raa);">RAA Det</th>
                    </tr></thead>
                    <tbody></tbody>
                </table>
            </div>
            <div class="card full-width">
                <div class="card-title">
                    Country Coverage Breakdown
                    <span class="badge badge-raa">Geographical</span>
                </div>
                <table id="countryTable">
                    <thead><tr>
                        <th>Country</th><th>N</th>
                        <th style="color:var(--accent-pcr);">PCR Amp</th>
                        <th style="color:var(--accent-pcr);">PCR Det</th>
                        <th style="color:var(--accent-raa);">RAA Amp</th>
                        <th style="color:var(--accent-raa);">RAA Det</th>
                    </tr></thead>
                    <tbody></tbody>
                </table>
            </div>
        </div>
    </div>
    <div id="designTables" style="display:none;">
        <div class="main-layout">
            <div class="card">
                <div class="card-title">
                    Subtype Distribution (Design DB)
                    <div style="display:flex;gap:0.5rem;">
                        <span class="badge badge-pcr">PCR</span>
                        <span class="badge badge-raa">RAA</span>
                    </div>
                </div>
                <table id="designSubtypeTable">
                    <thead><tr>
                        <th>Subtype</th><th>N</th>
                        <th style="color:var(--accent-pcr);">PCR Amp</th>
                        <th style="color:var(--accent-pcr);">PCR Det</th>
                        <th style="color:var(--accent-raa);">RAA Amp</th>
                        <th style="color:var(--accent-raa);">RAA Det</th>
                    </tr></thead>
                    <tbody></tbody>
                </table>
            </div>
            <div class="card">
                <div class="card-title">
                    Collection Year (Design DB)
                    <span class="badge badge-pcr">Temporal</span>
                </div>
                <table id="designYearTable">
                    <thead><tr>
                        <th>Year</th><th>N</th>
                        <th style="color:var(--accent-pcr);">PCR Amp</th>
                        <th style="color:var(--accent-pcr);">PCR Det</th>
                        <th style="color:var(--accent-raa);">RAA Amp</th>
                        <th style="color:var(--accent-raa);">RAA Det</th>
                    </tr></thead>
                    <tbody></tbody>
                </table>
            </div>
            <div class="card full-width">
                <div class="card-title">
                    Country Distribution (Design DB)
                    <span class="badge badge-raa">Geographical</span>
                </div>
                <table id="designCountryTable">
                    <thead><tr>
                        <th>Country</th><th>N</th>
                        <th style="color:var(--accent-pcr);">PCR Amp</th>
                        <th style="color:var(--accent-pcr);">PCR Det</th>
                        <th style="color:var(--accent-raa);">RAA Amp</th>
                        <th style="color:var(--accent-raa);">RAA Det</th>
                    </tr></thead>
                    <tbody></tbody>
                </table>
            </div>
        </div>
    </div>

</div>

<script>
""" + "    const candidatesData = " + data_json + """;
    let currentMode = 'global';
    const sel = document.getElementById('candidateSelect');
    candidatesData.forEach((c, i) => {
        const opt = document.createElement('option');
        opt.value = i;
        const pa = c.pcr_amplified_pct ?? 0;
        const ra = c.raa_amplified_pct ?? '\u2014';
        const nc = c.n_combinations > 1 ? ` [${c.n_combinations} combos]` : '';
        opt.textContent = `${c.cand_id}${nc}  |  PCR Amp: ${pa}%  \u00b7  RAA Amp: ${ra}%`;
        sel.appendChild(opt);
    });
    function pctStyle(pct) {
        const v = parseFloat(pct || 0);
        if (v >= 80) return 'background:rgba(52,211,153,0.85);color:#052e1c;';
        if (v >= 60) return 'background:rgba(16,185,129,0.7);color:#052e1c;';
        if (v >= 35) return 'background:rgba(251,191,36,0.85);color:#1a0e00;';
        if (v >= 15) return 'background:rgba(251,146,60,0.85);color:#1a0500;';
        return 'background:rgba(248,113,113,0.85);color:#1a0000;';
    }
    function barColor(pct, type) {
        const v = parseFloat(pct || 0);
        if (type === 'raa') { return v >= 70 ? '#6366f1' : v >= 40 ? '#818cf8' : '#f87171'; }
        return v >= 70 ? '#10b981' : v >= 40 ? '#fbbf24' : '#f87171';
    }
    function makeCell(pct, type) {
        const td = document.createElement('td');
        const span = document.createElement('span');
        span.className = 'cov-cell';
        span.style.cssText = pctStyle(pct);
        span.textContent = pct !== null && pct !== undefined ? `${pct}%` : '\u2014';
        const bar = document.createElement('div'); bar.className = 'mini-bar';
        const fill = document.createElement('div'); fill.className = 'mini-bar-fill';
        fill.style.width = `${Math.min(parseFloat(pct||0), 100)}%`;
        fill.style.background = barColor(pct, type);
        bar.appendChild(fill);
        const w = document.createElement('div');
        w.appendChild(span); w.appendChild(bar);
        td.appendChild(w);
        return td;
    }
    function populateTable(tid, data, catKey, maxRows) {
        const tbody = document.querySelector(`#${tid} tbody`);
        tbody.innerHTML = '';
        if (!data) return;
        const rows = maxRows ? data.slice(0, maxRows) : data;
        rows.forEach(row => {
            const tr = document.createElement('tr');
            const tdCat = document.createElement('td'); tdCat.textContent = row[catKey] || '\u2014'; tr.appendChild(tdCat);
            const tdN = document.createElement('td'); tdN.textContent = (row.total||0).toLocaleString(); tdN.style.color='var(--text-2)'; tr.appendChild(tdN);
            tr.appendChild(makeCell(row.pcr_amplified_pct ?? row.amplified_pct, 'pcr'));
            tr.appendChild(makeCell(row.pcr_detected_pct ?? row.detected_pct, 'pcr'));
            tr.appendChild(makeCell(row.raa_amplified_pct, 'raa'));
            tr.appendChild(makeCell(row.raa_detected_pct, 'raa'));
            tbody.appendChild(tr);
        });
        if (maxRows && data.length > maxRows) {
            const tr = document.createElement('tr');
            const td = document.createElement('td'); td.colSpan = 6;
            td.style.cssText = 'text-align:center;color:var(--text-2);font-style:italic;padding:0.5rem;';
            td.textContent = `Showing ${maxRows} of ${data.length} entries`;
            tr.appendChild(td); tbody.appendChild(tr);
        }
    }

    function populateDesignTable(tid, data, catKey, maxRows) {
        const tbody = document.querySelector(`#${tid} tbody`);
        tbody.innerHTML = '';
        if (!data || data.length === 0) {
            const tr = document.createElement('tr');
            const td = document.createElement('td'); td.colSpan = 6;
            td.style.cssText = 'text-align:center;color:var(--text-2);font-style:italic;padding:1rem;';
            td.textContent = 'No metadata available for design DB sequences';
            tr.appendChild(td); tbody.appendChild(tr);
            return;
        }
        const rows = maxRows ? data.slice(0, maxRows) : data;
        rows.forEach(row => {
            const tr = document.createElement('tr');
            const tdCat = document.createElement('td'); tdCat.textContent = row[catKey] || '\u2014'; tr.appendChild(tdCat);
            const tdN = document.createElement('td'); tdN.textContent = (row.total||0).toLocaleString(); tdN.style.color='var(--text-2)'; tr.appendChild(tdN);
            // PCR Amp (same column structure as global: use pcr_amplified_pct)
            tr.appendChild(makeCell(row.pcr_amplified_pct ?? 0, 'pcr'));
            // PCR Det: design DB PCR has no separate detection step (same as amp for design DB)
            // Fall back to pcr_amplified_pct if pcr_detected_pct is absent
            tr.appendChild(makeCell(row.pcr_amplified_pct ?? 0, 'pcr'));
            // RAA Amp + RAA Det — computed from design DB amplicons by m10_insilico_pcr.py
            tr.appendChild(makeCell(row.raa_amplified_pct ?? null, 'raa'));
            tr.appendChild(makeCell(row.raa_detected_pct  ?? null, 'raa'));
            tbody.appendChild(tr);
        });
        if (maxRows && data.length > maxRows) {
            const tr = document.createElement('tr');
            const td = document.createElement('td'); td.colSpan = 6;
            td.style.cssText = 'text-align:center;color:var(--text-2);font-style:italic;padding:0.5rem;';
            td.textContent = `Showing ${maxRows} of ${data.length} entries`;
            tr.appendChild(td); tbody.appendChild(tr);
        }
    }
    function setDatasetMode(mode) {
        currentMode = mode;
        document.querySelectorAll('.toggle-btn').forEach(b => b.classList.remove('active-global','active-design'));
        if (mode === 'global') {
            document.getElementById('btn-global').classList.add('active-global');
            document.getElementById('globalTables').style.display = 'block';
            document.getElementById('designTables').style.display = 'none';
        } else {
            document.getElementById('btn-design').classList.add('active-design');
            document.getElementById('globalTables').style.display = 'none';
            document.getElementById('designTables').style.display = 'block';
        }
        updateDashboard(sel.value);
    }
    function updateDashboard(idx) {
        const c = candidatesData[idx]; if (!c) return;
        document.getElementById('statRegion').textContent = `${c.region} \u00b7 Rank ${c.rank}`;
        const total = currentMode === 'global' ? (c.total_seqs||0) : (c.design_db_total||0);
        document.getElementById('statTotal').textContent = total.toLocaleString();
        let pcrAmp=0,pcrDet=0,raaAmp=null,raaDet=null;
        if (currentMode === 'global') {
            pcrAmp=c.pcr_amplified_pct??0; pcrDet=c.pcr_detected_pct??0;
            raaAmp=c.raa_amplified_pct; raaDet=c.raa_detected_pct;
        } else {
            // Design DB mode: all metrics (PCR + RAA) now use design DB.
            // RAA binding was computed per-sequence on design DB amplicons (same sliding-window
            // model used for surveillance DB), so design_db_raa_* holds accurate values.
            pcrAmp=c.design_db_amplified_pct??0; pcrDet=pcrAmp;
            raaAmp=c.design_db_raa_amp_pct??null; raaDet=c.design_db_raa_det_pct??null;
        }
        // Hide "(Global Surveillance DB)" note — RAA now from design DB when in design mode
        ['raaAmpNote','raaDetNote'].forEach(id => {
            document.getElementById(id).style.display = 'none';
        });
        document.getElementById('pcrAmpPct').textContent=`${pcrAmp}%`;
        document.getElementById('pcrDetPct').textContent=`${pcrDet}%`;
        document.getElementById('raaAmpPct').textContent=raaAmp!=null?`${raaAmp}%`:'\u2014';
        document.getElementById('raaDetPct').textContent=raaDet!=null?`${raaDet}%`:'\u2014';
        [['pcrAmpBar',pcrAmp,'pcr'],['pcrDetBar',pcrDet,'pcr'],['raaAmpBar',raaAmp,'raa'],['raaDetBar',raaDet,'raa']].forEach(([id,v,t])=>{
            const el=document.getElementById(id);
            el.style.width=`${Math.min(parseFloat(v||0),100)}%`;
            el.style.background=barColor(v,t);
        });
        const banner=document.getElementById('designBanner');
        if (c.design_db_total>0 && currentMode==='global') {
            banner.style.display='block';
            document.getElementById('designTotal').textContent=c.design_db_total.toLocaleString();
            document.getElementById('designAmplified').textContent=(c.design_db_amplified??0).toLocaleString();
            const dp=c.design_db_amplified_pct||0;
            document.getElementById('designPct').textContent=`${dp}%`;
            // RAA stats for design DB
            const dr=c.design_db_raa_amp_pct??null;
            document.getElementById('designRaaAmplified').textContent=dr!==null?(c.design_db_raa_amplified??0).toLocaleString():'\u2014';
            document.getElementById('designRaaPct').textContent=dr!==null?`${dr}%`:'\u2014';
            const db=document.getElementById('designBar');
            db.style.width=`${Math.min(dp,100)}%`;
            db.style.background=dp>=80?'linear-gradient(90deg,#059669,#34d399)':dp>=50?'linear-gradient(90deg,#7c3aed,#a78bfa)':dp>=30?'linear-gradient(90deg,#d97706,#fbbf24)':'linear-gradient(90deg,#dc2626,#f87171)';
        } else { banner.style.display='none'; }
        document.getElementById('seqFwd').textContent=c.fwd_seq||'\u2014';
        document.getElementById('seqRev').textContent=c.rev_seq||'\u2014';
        document.getElementById('seqProbe').textContent=c.probe_seq||'\u2014 (No Probe)';
        const fwdDeg=document.getElementById('seqFwdDeg');
        const revDeg=document.getElementById('seqRevDeg');
        const prbDeg=document.getElementById('seqProbeDeg');
        if (c.fwd_degenerate_seq && c.fwd_degenerate_seq!==c.fwd_seq) {
            fwdDeg.innerHTML=`\U0001f9ec IUPAC: <span style="font-family:'JetBrains Mono',monospace;color:#c4b5fd;font-size:0.8rem;">${c.fwd_degenerate_seq}</span>`+(c.n_combinations>1?` \u00b7 ${c.n_combinations} combos`:'');
            fwdDeg.style.display='block';
        } else { fwdDeg.style.display='none'; }
        if (c.rev_degenerate_seq && c.rev_degenerate_seq!==c.rev_seq) {
            revDeg.innerHTML=`\U0001f9ec IUPAC: <span style="font-family:'JetBrains Mono',monospace;color:#c4b5fd;font-size:0.8rem;">${c.rev_degenerate_seq}</span>`;
            revDeg.style.display='block';
        } else { revDeg.style.display='none'; }
        if (c.probe_seq && /[RYMKSWBDHVN]/i.test(c.probe_seq)) {
            prbDeg.textContent='\U0001f9ec Contains IUPAC degenerate bases';
            prbDeg.style.display='block';
        } else { prbDeg.style.display='none'; }
        if (currentMode==='global' && c.breakdowns) {
            populateTable('subtypeTable',c.breakdowns.subtype,'subtype');
            populateTable('yearTable',c.breakdowns.year,'year');
            populateTable('countryTable',c.breakdowns.country,'country',30);
        }
        if (currentMode==='design' && c.design_db_breakdowns) {
            const ddb = c.design_db_breakdowns;
            populateDesignTable('designSubtypeTable', ddb.subtype, 'subtype');
            populateDesignTable('designYearTable',    ddb.year,    'year');
            populateDesignTable('designCountryTable', ddb.country, 'country', 30);
        }
    }
    if (candidatesData.length > 0) updateDashboard(0);
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
