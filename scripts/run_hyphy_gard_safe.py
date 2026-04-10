#!/usr/bin/env python3
"""
Safe GARD wrapper — runs HyPhy GARD with pre-flight checks.
If the alignment fails feasibility checks (too short, too few variable sites),
or if GARD itself crashes, this script writes a valid empty JSON output
and exits with code 0 so the pipeline continues gracefully.

Usage:
    python run_hyphy_gard_safe.py <alignment.fasta> <output.json> [--rv GDD]
"""

import json
import subprocess
import sys
from pathlib import Path
from Bio import SeqIO


def preflight_check(alignment_path):
    """Check if alignment meets GARD's minimum requirements.
    GARD needs alignment_length > 2 * num_sequences for c-AIC.
    Returns (ok: bool, reason: str).
    """
    records = list(SeqIO.parse(alignment_path, "fasta"))
    n_seqs = len(records)
    if n_seqs == 0:
        return False, "Alignment is empty (0 sequences)"

    aln_len = len(records[0].seq)
    # GARD minimum: at least 2 * n_seqs sites for a two-partition model
    min_sites = 2 * n_seqs * 2  # conservative estimate matching HyPhy's internal check
    if aln_len < min_sites:
        return False, (
            f"Alignment too short for GARD c-AIC model comparison. "
            f"Have {aln_len} sites with {n_seqs} sequences, "
            f"need at least ~{min_sites} sites. "
            f"This is a mathematical limitation, not a data quality issue."
        )

    # Check variable sites
    variable = 0
    for i in range(aln_len):
        col = set(str(r.seq[i]).upper() for r in records) - {"-", "N"}
        if len(col) > 1:
            variable += 1
    if variable < 10:
        return False, f"Too few variable sites ({variable}). GARD needs sufficient variation."

    return True, "Preflight OK"


def write_empty_gard_json(output_path, reason):
    """Write a valid but empty GARD result JSON that downstream tools can parse."""
    result = {
        "analysis": {
            "info": "GARD analysis was skipped (safe mode)",
            "skip_reason": reason,
            "version": "safe-wrapper-1.0"
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
        print("Usage: python run_hyphy_gard_safe.py <alignment> <output> [--rv GDD]",
              file=sys.stderr)
        sys.exit(1)

    alignment = sys.argv[1]
    output = sys.argv[2]
    rv = "GDD"
    if "--rv" in sys.argv:
        rv_idx = sys.argv.index("--rv")
        if rv_idx + 1 < len(sys.argv):
            rv = sys.argv[rv_idx + 1]

    print(f"[GARD-SAFE] Pre-flight check on {alignment}...", file=sys.stderr)

    # --- Pre-flight check ---
    ok, reason = preflight_check(alignment)
    if not ok:
        print(f"[GARD-SAFE] ⚠ SKIPPED: {reason}", file=sys.stderr)
        print(f"[GARD-SAFE] Writing empty GARD result to {output}", file=sys.stderr)
        write_empty_gard_json(output, reason)
        sys.exit(0)  # Exit 0 so Snakemake continues

    print(f"[GARD-SAFE] Pre-flight passed. Running HyPhy GARD...", file=sys.stderr)

    # --- Run GARD ---
    cmd = [
        "hyphy", "gard",
        "--alignment", alignment,
        "--output", output,
        "--rv", rv
    ]

    try:
        result = subprocess.run(cmd, capture_output=False, text=True)
        if result.returncode != 0:
            print(f"[GARD-SAFE] ⚠ GARD exited with code {result.returncode}",
                  file=sys.stderr)
            # Check if output was partially written
            if not Path(output).exists() or Path(output).stat().st_size == 0:
                reason = f"GARD exited with code {result.returncode}"
                print(f"[GARD-SAFE] Writing empty GARD result due to runtime failure",
                      file=sys.stderr)
                write_empty_gard_json(output, reason)
            # Always exit 0 for soft-fail
            sys.exit(0)
    except FileNotFoundError:
        reason = "HyPhy binary not found in PATH"
        print(f"[GARD-SAFE] ⚠ SKIPPED: {reason}", file=sys.stderr)
        write_empty_gard_json(output, reason)
        sys.exit(0)
    except Exception as e:
        reason = f"Unexpected error: {str(e)}"
        print(f"[GARD-SAFE] ⚠ SKIPPED: {reason}", file=sys.stderr)
        write_empty_gard_json(output, reason)
        sys.exit(0)

    print(f"[GARD-SAFE] ✓ GARD completed successfully.", file=sys.stderr)


if __name__ == "__main__":
    main()
