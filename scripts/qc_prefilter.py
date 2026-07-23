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


def _date_key(record_id: str, id_to_date: dict) -> str:
    """Sort key for choosing which duplicate to keep: earliest real date first.

    Missing/unknown dates map to a value that sorts after any real ISO date, so a
    record with a known collection date always beats an undated duplicate.
    """
    lookup = record_id.split("|")[0].strip()
    value = (id_to_date.get(lookup) or "").strip()
    if not value or value.upper().startswith("X"):
        return "9999-99-99"
    return value


def seq_hash(record):
    """Generate MD5 hash of the uppercase sequence for deduplication."""
    return hashlib.md5(str(record.seq).upper().encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description="QC Prefilter for genome sequences")
    parser.add_argument("--input-fasta", required=True, help="Input FASTA file")
    parser.add_argument("--output-fasta", required=True, help="Output filtered FASTA file")
    parser.add_argument("--report", required=True, help="Output QC report JSON file")
    parser.add_argument("--metadata", help="Optional metadata file (CSV/TSV) for subtype filtering")
    parser.add_argument("--min-subtype-count", type=int, default=30,
                        help="Minimum number of sequences required per subtype (default: 30)")
    parser.add_argument("--max-n-pct", type=float, default=5.0,
                        help="Maximum allowed N-content percentage (default: 5.0)")
    args = parser.parse_args()

    input_path = Path(args.input_fasta)
    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    # --- Phase 1: Read Metadata (Optional but recommended for subtype filtering) ---
    blacklist_ids = set()
    subtype_stats = {}
    id_to_date = {}
    
    if args.metadata:
        print(f"[QC] Loading metadata for subtype filtering: {args.metadata}")
        import csv
        
        # 1. Count subtypes
        counts = {}
        id_to_subtype = {}
        id_to_date = {}
        
        # Determine delimiter
        delim = "\t" if str(args.metadata).endswith(".tsv") else ","
        
        try:
            with open(args.metadata, "r", encoding="utf-8", errors="ignore") as f:
                reader = csv.DictReader(f, delimiter=delim)
                fieldnames = reader.fieldnames or []

                # Case-insensitive column lookup. M02 runs on metadata AFTER
                # sanitize_metadata.py, which renames the ID column to 'strain'
                # and the subtype column to lowercase 'subtype'. We still accept
                # the legacy raw GISAID/NCBI names so pre-sanitize files also work.
                def find_col(candidates):
                    lower_map = {c.lower(): c for c in fieldnames}
                    for cand in candidates:
                        if cand.lower() in lower_map:
                            return lower_map[cand.lower()]
                    # substring fallback (e.g. "MP Segment_Id")
                    for cand in candidates:
                        for c in fieldnames:
                            if cand.lower() in c.lower():
                                return c
                    return None

                id_col = find_col(["strain", "Segment_Id", "Isolate_Id", "Accession", "gisaid_epi_isl"])
                subtype_col = find_col(["subtype", "Subtype", "Subtype_Full", "serotype"])
                date_col = find_col(["date", "Collection_Date", "collection_date", "Collection Date"])

                if not id_col:
                    print("WARNING: Could not find ID column (tried strain/Segment_Id/Isolate_Id/Accession). "
                          "Subtype filtering might be limited.")
                if not subtype_col:
                    print("WARNING: Could not find subtype column (tried subtype/Subtype/serotype). "
                          "Subtype scarcity filter will be a no-op.")
                else:
                    print(f"[QC] Using ID column '{id_col}' and subtype column '{subtype_col}'.")

                for row in reader:
                    st = (row.get(subtype_col, "Unknown") if subtype_col else "Unknown").strip() or "Unknown"
                    counts[st] = counts.get(st, 0) + 1

                    if id_col:
                        # GISAID IDs often look like "EPI123456|NAME"
                        raw_id = row.get(id_col, "")
                        sid = raw_id.split("|")[0].strip() if "|" in raw_id else raw_id.strip()
                        if sid:
                            id_to_subtype[sid] = st
                            if date_col:
                                id_to_date[sid] = (row.get(date_col) or "").strip()
            
            # 2. Identify blacklisted subtypes
            removed_subtypes = []
            for st, count in counts.items():
                if count < args.min_subtype_count:
                    removed_subtypes.append(st)
                subtype_stats[st] = count
            
            print(f"[QC] Subtype Analysis: Found {len(counts)} subtypes.")
            if removed_subtypes:
                print(f"[QC]   Blacklisted (n < {args.min_subtype_count}): {removed_subtypes}")
                
            # 3. Mark IDs for removal
            # We don't know the exact FASTA header format yet, so we store subtype info
            # for the filtering phase below.
        except Exception as e:
            print(f"ERROR: Failed to process metadata: {e}")

    # --- Phase 2: Read and Filter sequences ---
    print(f"[QC] Reading sequences from: {input_path}")
    
    total_input = 0
    n_content_removed = 0
    invalid_char_removed = 0
    subtype_removed = 0
    
    passed_filters = []
    seen_hashes = OrderedDict()
    duplicates_removed = 0
    
    VALID_IUPAC = set("ACGTURSYWKMBDHVN-.")
    
    for record in SeqIO.parse(str(input_path), "fasta"):
        total_input += 1
        rid = record.id
        
        # A. Subtype Data Availability Filter
        if args.metadata:
            # Match FASTA ID to Metadata ID
            # Fasta IDs can be "EPI123456" or "EPI_ISL_123456" or "EPI123|Name"
            lookup_id = rid.split("|")[0].strip()
            st = id_to_subtype.get(lookup_id)
            
            if st and counts.get(st, 0) < args.min_subtype_count:
                subtype_removed += 1
                if subtype_removed <= 5:
                    print(f"[QC]   Removed (subtype {st} n={counts[st]}): {rid}")
                continue

        # B. N-content filter
        pct = calculate_n_pct(record)
        if pct > args.max_n_pct:
            n_content_removed += 1
            if n_content_removed <= 5:
                print(f"[QC]   Removed (N={pct:.1f}%): {rid}")
            continue

        # C. Character validation
        seq_upper = str(record.seq).upper()
        bad_chars = set(seq_upper) - VALID_IUPAC
        if bad_chars:
            invalid_char_removed += 1
            if invalid_char_removed <= 5:
                print(f"[QC]   Removed (invalid chars {bad_chars}): {rid}")
            continue

        # D. Deduplication — keep the EARLIEST collection date.
        # Previously whichever record SeqIO happened to parse first won, so the
        # retained duplicate depended on raw FASTA ordering, contradicting the
        # documented behaviour ("the record with the earliest collection date is
        # retained"). Unknown/'' dates sort last so a dated record always wins.
        h = seq_hash(record)
        if h in seen_hashes:
            duplicates_removed += 1
            if _date_key(record.id, id_to_date) < _date_key(seen_hashes[h].id, id_to_date):
                seen_hashes[h] = record
        else:
            seen_hashes[h] = record

    deduplicated = list(seen_hashes.values())
    
    print(f"[QC] Summary of Filtering:")
    print(f"[QC]   Total Input: {total_input}")
    if args.metadata:
        print(f"[QC]   Removed by Subtype scarcity (<{args.min_subtype_count}): {subtype_removed}")
    print(f"[QC]   Removed by N-content (>{args.max_n_pct}%): {n_content_removed}")
    print(f"[QC]   Removed by Invalid Chars: {invalid_char_removed}")
    print(f"[QC]   Removed by Deduplication: {duplicates_removed}")
    print(f"[QC]   Final Output: {len(deduplicated)}")

    # --- Phase 3: Write filtered output ---
    output_path = Path(args.output_fasta)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(deduplicated, str(output_path), "fasta")
    print(f"[QC] Wrote {len(deduplicated)} sequences to: {output_path}")

    # --- Phase 4: Write QC report ---
    report = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "input_file": str(input_path),
        "input_count": total_input,
        "subtype_filtering": {
            "metadata_file": args.metadata,
            "min_count": args.min_subtype_count,
            "removed_count": subtype_removed,
            "all_subtypes_found": subtype_stats
        },
        "filters": {
            "n_content": {
                "max_n_pct": args.max_n_pct,
                "removed_count": n_content_removed,
            },
            "invalid_characters": {
                "removed_count": invalid_char_removed,
            },
            "deduplication": {
                "duplicates_removed": duplicates_removed,
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
