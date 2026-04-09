"""
QC Prefilter — Sequence quality control before augur filter.

Performs:
  1. N-content filtering (removes sequences with too many ambiguous bases)
  2. Deduplication (removes identical sequences by MD5 hash)
  3. Generates a QC report JSON with statistics

Usage:
  python qc_prefilter.py \\
      --input-fasta raw_sequences.fasta \\
      --output-fasta prefiltered_sequences.fasta \\
      --report qc_report.json \\
      --max-n-pct 5.0
"""

import argparse
import hashlib
import json
import sys
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path

from Bio import SeqIO


def calculate_n_pct(record):
    """Calculate the percentage of ambiguous (N) bases in a sequence."""
    seq_str = str(record.seq).upper()
    if len(seq_str) == 0:
        return 100.0
    n_count = seq_str.count("N")
    return (n_count / len(seq_str)) * 100.0


def seq_hash(record):
    """Generate MD5 hash of the uppercase sequence for deduplication."""
    return hashlib.md5(str(record.seq).upper().encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description="QC Prefilter for genome sequences")
    parser.add_argument("--input-fasta", required=True, help="Input FASTA file")
    parser.add_argument("--output-fasta", required=True, help="Output filtered FASTA file")
    parser.add_argument("--report", required=True, help="Output QC report JSON file")
    parser.add_argument("--max-n-pct", type=float, default=5.0,
                        help="Maximum allowed N-content percentage (default: 5.0)")
    args = parser.parse_args()

    input_path = Path(args.input_fasta)
    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    # --- Phase 1: Read all sequences ---
    print(f"[QC] Reading sequences from: {input_path}")
    records = list(SeqIO.parse(str(input_path), "fasta"))
    total_input = len(records)
    print(f"[QC] Total input sequences: {total_input}")

    if total_input == 0:
        print("ERROR: Input FASTA contains no sequences.", file=sys.stderr)
        sys.exit(1)

    # Collect length stats from the raw input
    lengths = [len(r.seq) for r in records]
    length_stats = {
        "min": min(lengths),
        "max": max(lengths),
        "mean": round(sum(lengths) / len(lengths), 1),
        "median": sorted(lengths)[len(lengths) // 2],
    }

    # --- Phase 2: N-content filter ---
    n_content_removed = 0
    n_content_removed_ids = []
    passed_n_filter = []

    for record in records:
        pct = calculate_n_pct(record)
        if pct > args.max_n_pct:
            n_content_removed += 1
            n_content_removed_ids.append(record.id)
            print(f"[QC]   Removed (N={pct:.1f}%): {record.id}")
        else:
            passed_n_filter.append(record)

    print(f"[QC] After N-content filter ({args.max_n_pct}%): "
          f"{len(passed_n_filter)} kept, {n_content_removed} removed")

    # --- Phase 3: Deduplication by MD5 hash ---
    seen_hashes = OrderedDict()
    duplicates_removed = 0
    duplicate_ids = []

    for record in passed_n_filter:
        h = seq_hash(record)
        if h in seen_hashes:
            duplicates_removed += 1
            duplicate_ids.append(record.id)
        else:
            seen_hashes[h] = record

    deduplicated = list(seen_hashes.values())
    print(f"[QC] After deduplication: "
          f"{len(deduplicated)} unique, {duplicates_removed} duplicates removed")

    # --- Phase 4: Write filtered output ---
    output_path = Path(args.output_fasta)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(deduplicated, str(output_path), "fasta")
    print(f"[QC] Wrote {len(deduplicated)} sequences to: {output_path}")

    # --- Phase 5: Write QC report ---
    report = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "input_file": str(input_path),
        "input_count": total_input,
        "length_stats_raw": length_stats,
        "filters": {
            "n_content": {
                "max_n_pct": args.max_n_pct,
                "removed_count": n_content_removed,
                "removed_ids": n_content_removed_ids[:50],  # Cap at 50 to avoid bloat
            },
            "deduplication": {
                "duplicates_removed": duplicates_removed,
                "duplicate_ids": duplicate_ids[:50],
            },
        },
        "output_count": len(deduplicated),
        "retention_pct": round((len(deduplicated) / total_input) * 100, 1) if total_input > 0 else 0,
    }

    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"[QC] Report written to: {report_path}")
    print(f"[QC] Summary: {total_input} → {len(deduplicated)} sequences "
          f"({report['retention_pct']}% retained)")


if __name__ == "__main__":
    main()
