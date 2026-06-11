#!/usr/bin/env python3
"""
Safe GARD wrapper — runs HyPhy GARD with pre-flight checks and auto-subsampling.

Strategy:
1. Pre-flight: Check if alignment meets GARD's mathematical requirements.
2. If too many sequences for the alignment length, AUTO-SUBSAMPLE to a safe
   number (keeping representative diversity via random seed).
3. Run GARD on the (possibly subsampled) alignment.
4. If GARD crashes at runtime, catch the error and write a valid empty JSON
   so the pipeline continues gracefully (soft-fail).

Usage:
    python run_hyphy_gard_safe.py <alignment.fasta> <output.json> [--rv GDD] [--seed 42]
"""

import json
import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# Safety margin: we target max_seqs = alignment_length / SAFETY_DIVISOR
# HyPhy internally needs ~2*n_seqs*2 sites. Using divisor 5 gives headroom.
SAFETY_DIVISOR = 5
MIN_SEQS_FOR_GARD = 20  # GARD needs reasonable diversity


def assess_alignment(alignment_path):
    """Read alignment and return (records, n_seqs, aln_len, n_variable)."""
    records = list(SeqIO.parse(alignment_path, "fasta"))
    n_seqs = len(records)
    if n_seqs == 0:
        return records, 0, 0, 0

    aln_len = len(records[0].seq)

    # Count variable sites (quick scan, sample columns if too large)
    variable = 0
    step = max(1, aln_len // 2000)  # sample up to 2000 columns for speed
    for i in range(0, aln_len, step):
        col = set(str(r.seq[i]).upper() for r in records) - {"-", "N"}
        if len(col) > 1:
            variable += 1
    if step > 1:
        variable = int(variable * (aln_len / (aln_len // step)))  # extrapolate

    return records, n_seqs, aln_len, variable


def subsample_alignment(records, max_seqs, seed=42):
    """Randomly subsample records to max_seqs, preserving diversity via seed."""
    random.seed(seed)
    indices = random.sample(range(len(records)), max_seqs)
    return [records[i] for i in sorted(indices)]


def write_empty_gard_json(output_path, reason):
    """Write a valid but empty GARD result JSON that downstream tools can parse."""
    result = {
        "analysis": {
            "info": "GARD analysis was skipped (safe mode)",
            "skip_reason": reason,
            "version": "safe-wrapper-2.0"
        },
        "breakpointData": {},
        "improvements": {},
        "singleTreeAICc": 0,
        "totalModelCount": 0,
        "potentialBreakpoints": 0,
        "skipped": True
    }
    with open(output_path, "w") as f:
        json.dump(result, f, indent=2)


def main():
    if len(sys.argv) < 3:
        print("Usage: python run_hyphy_gard_safe.py <alignment> <output> [--rv GDD] [--seed 42]",
              file=sys.stderr)
        sys.exit(1)

    alignment = sys.argv[1]
    output = sys.argv[2]
    rv = "GDD"
    seed = 42

    if "--rv" in sys.argv:
        rv_idx = sys.argv.index("--rv")
        if rv_idx + 1 < len(sys.argv):
            rv = sys.argv[rv_idx + 1]
    if "--seed" in sys.argv:
        seed_idx = sys.argv.index("--seed")
        if seed_idx + 1 < len(sys.argv):
            seed = int(sys.argv[seed_idx + 1])

    print(f"[GARD-SAFE] Assessing alignment: {alignment}", file=sys.stderr)

    # ---- Step 1: Assess alignment ----
    records, n_seqs, aln_len, n_variable = assess_alignment(alignment)

    if n_seqs == 0:
        print("[GARD-SAFE] [SKIPPED]: Alignment is empty.", file=sys.stderr)
        write_empty_gard_json(output, "Alignment is empty (0 sequences)")
        sys.exit(0)

    if n_variable < 10:
        reason = f"Too few variable sites (~{n_variable}). GARD needs sufficient variation."
        print(f"[GARD-SAFE] [SKIPPED]: {reason}", file=sys.stderr)
        write_empty_gard_json(output, reason)
        sys.exit(0)

    max_safe_seqs = aln_len // SAFETY_DIVISOR
    print(f"[GARD-SAFE] Alignment: {n_seqs} seqs × {aln_len} sites, "
          f"~{n_variable} variable sites", file=sys.stderr)
    print(f"[GARD-SAFE] Max safe sequences for GARD: {max_safe_seqs} "
          f"(alignment_length / {SAFETY_DIVISOR})", file=sys.stderr)

    if max_safe_seqs < MIN_SEQS_FOR_GARD:
        reason = (f"Alignment too short for meaningful GARD analysis. "
                  f"{aln_len} sites can support max ~{max_safe_seqs} sequences, "
                  f"but need at least {MIN_SEQS_FOR_GARD}.")
        print(f"[GARD-SAFE] [SKIPPED]: {reason}", file=sys.stderr)
        write_empty_gard_json(output, reason)
        sys.exit(0)

    # ---- Step 2: Subsample if needed ----
    actual_alignment = alignment
    tmp_file = None

    if n_seqs > max_safe_seqs:
        print(f"[GARD-SAFE] [SUBSAMPLING] Auto-subsampling: {n_seqs} -> {max_safe_seqs} sequences "
              f"(seed={seed})", file=sys.stderr)
        subsampled = subsample_alignment(records, max_safe_seqs, seed=seed)

        # Write subsampled alignment to temp file
        tmp_dir = os.path.dirname(output) or "."
        tmp_fd, tmp_path = tempfile.mkstemp(suffix="_gard_subsample.fasta", dir=tmp_dir)
        os.close(tmp_fd)
        SeqIO.write(subsampled, tmp_path, "fasta")
        actual_alignment = tmp_path
        tmp_file = tmp_path
        print(f"[GARD-SAFE] Subsampled alignment written to {tmp_path}", file=sys.stderr)
    else:
        print(f"[GARD-SAFE] [OK] No subsampling needed ({n_seqs} <= {max_safe_seqs})",
              file=sys.stderr)

    threads = 1
    if "--threads" in sys.argv:
        t_idx = sys.argv.index("--threads")
        if t_idx + 1 < len(sys.argv):
            threads = int(sys.argv[t_idx + 1])

    # ---- Step 3: Run GARD ----
    print(f"[GARD-SAFE] Running HyPhy GARD with {threads} threads...", file=sys.stderr)
    cmd = [
        "hyphy", f"CPU={threads}", "gard",
        "--alignment", actual_alignment,
        "--output", output,
        "--rv", rv
    ]

    try:
        result = subprocess.run(cmd, capture_output=False, text=True)
        if result.returncode != 0:
            print(f"[GARD-SAFE] [WARNING] GARD exited with code {result.returncode}",
                  file=sys.stderr)
            if not Path(output).exists() or Path(output).stat().st_size == 0:
                reason = f"GARD runtime failure (exit code {result.returncode})"
                print(f"[GARD-SAFE] Writing empty GARD result", file=sys.stderr)
                write_empty_gard_json(output, reason)
            sys.exit(0)
    except FileNotFoundError:
        reason = "HyPhy binary not found in PATH"
        print(f"[GARD-SAFE] [SKIPPED]: {reason}", file=sys.stderr)
        write_empty_gard_json(output, reason)
        sys.exit(0)
    except Exception as e:
        reason = f"Unexpected error: {str(e)}"
        print(f"[GARD-SAFE] [SKIPPED]: {reason}", file=sys.stderr)
        write_empty_gard_json(output, reason)
        sys.exit(0)
    finally:
        # Clean up temp file
        if tmp_file and os.path.exists(tmp_file):
            os.remove(tmp_file)
            print(f"[GARD-SAFE] Cleaned up temp subsampled file", file=sys.stderr)

    # ---- Step 4: Annotate output if subsampled ----
    if n_seqs > max_safe_seqs and Path(output).exists():
        try:
            with open(output, "r") as f:
                gard_result = json.load(f)
            gard_result["_safe_wrapper"] = {
                "original_sequences": n_seqs,
                "subsampled_to": max_safe_seqs,
                "seed": seed,
                "note": ("GARD was run on a representative subsample. "
                         "Breakpoint positions are in alignment coordinates "
                         "and are valid for the full dataset.")
            }
            with open(output, "w") as f:
                json.dump(gard_result, f, indent=2)
        except Exception:
            pass  # Non-critical annotation, don't fail

    print(f"[GARD-SAFE] [SUCCESS] GARD completed successfully.", file=sys.stderr)


if __name__ == "__main__":
    main()
