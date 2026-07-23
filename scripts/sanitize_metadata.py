"""
Metadata Sanitizer — Normalize and validate metadata from any source.

Supports auto-detection of source format (GISAID, NCBI, generic TSV).
Applies a 3-tier column validation system:
  - Tier 1 (Primary):     strain, date, country — essential for pipeline
  - Tier 2 (Secondary):   region, host, division — filled with defaults if empty
  - Tier 3 (Pass-through): all other columns preserved as-is

Usage:
  python sanitize_metadata.py \\
      --input metadata.tsv \\
      --output metadata_sanitized.tsv \\
      --report sanitize_report.json \\
      --format auto
"""

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

# Column mapping per source format

# GISAID Nextstrain-preprocessed format (already has standard columns)
GISAID_COLUMN_MAP = {
    "strain": "strain",
    "virus": "virus",
    "gisaid_epi_isl": "gisaid_epi_isl",
    "date": "date",
    "region": "region",
    "country": "country",
    "division": "division",
    "location": "location",
    "host": "host",
    "age": "age",
    "sex": "sex",
    "pangolin_lineage": "pangolin_lineage",
    "GISAID_clade": "GISAID_clade",
    "originating_lab": "originating_lab",
    "submitting_lab": "submitting_lab",
    "authors": "authors",
    "date_submitted": "date_submitted",
}

# GISAID raw download format (different column names)
GISAID_RAW_COLUMN_MAP = {
    "Virus name": "strain",
    "Accession ID": "gisaid_epi_isl",
    "Collection date": "date",
    "Location": "_location_raw",  # needs parsing
    "Host": "host",
    "Additional location information": "location",
    "Lineage": "pangolin_lineage",
    "Clade": "GISAID_clade",
    "Submission date": "date_submitted",
}

# GISAID EpiFlu download format (XLS format)
# NOTE: EpiFlu FASTA headers use Isolate_Id (EPI_ISL_xxx) as the primary ID,
# so we map that to 'strain' for augur compatibility.
GISAID_EPIFLU_MAP = {
    "Isolate_Id": "strain",
    "Isolate_Name": "virus_name",
    "Collection_Date": "date",
    "Location": "_location_raw",  # uses continent / country / ...
    "Host": "host",
    "Subtype": "subtype",
    "Lineage": "lineage",
    "Clade": "clade",
    "Submission_Date": "date_submitted",
    "Pathogenicity": "pathogenicity",
}

# NCBI (already handled by convert_ncbi_csv, but support as fallback)
NCBI_COLUMN_MAP = {
    "Accession": "strain",
    "Organism_Name": "virus",
    "Collection_Date": "date",
    "Geo_Location": "_geo_location_raw",  # needs parsing
    "Host": "host",
    "Pangolin": "pangolin_lineage",
    "Length": "length",
    "Release_Date": "date_submitted",
}

# Tier 1: Primary columns — pipeline may crash or produce bad results without these
TIER1_DEFAULTS = {
    "strain": None,       # None = DROP row if missing (fatal)
    "date": "XXXX-XX-XX", # augur convention for unknown date
    "country": "Unknown",
}

# Tier 2: Secondary columns — filled silently with defaults
TIER2_DEFAULTS = {
    "region": "Unknown",
    # NOT "Human": for an influenza pipeline, silently labelling a record with a
    # missing host as human misclassifies avian/swine-origin strains in the Auspice
    # host colour-by. "Unknown" is honest and matches the other Tier 2 defaults.
    "host": "Unknown",
    "division": "Unknown",
}


def detect_format(headers: list[str]) -> str:
    """Auto-detect metadata format from column headers."""
    header_set = set(headers)

    # Check for GISAID raw format
    if "Virus name" in header_set and "Accession ID" in header_set:
        return "gisaid_raw"

    # Check for GISAID EpiFlu format
    if "Isolate_Name" in header_set and "Isolate_Id" in header_set:
        return "gisaid_epiflu"

    # Check for GISAID Nextstrain-preprocessed format
    if "strain" in header_set and "gisaid_epi_isl" in header_set:
        return "gisaid_tsv"

    # Check for NCBI CSV format (shouldn't reach here normally,
    # but support as fallback)
    if "Accession" in header_set and "Organism_Name" in header_set:
        return "ncbi_csv"

    # Check for our internal pipeline format
    if "strain" in header_set and "date" in header_set:
        return "tsv"

    # Unknown — try to proceed with generic TSV using first column as strain
    return "generic"


def parse_gisaid_location(location_str: str) -> dict:
    """Parse GISAID Location field: 'Continent / Country / Division / City'."""
    parts = [p.strip() for p in location_str.split("/")]
    return {
        "region": parts[0] if len(parts) > 0 else "Unknown",
        "country": parts[1] if len(parts) > 1 else "Unknown",
        "division": parts[2] if len(parts) > 2 else "Unknown",
        "location": parts[3] if len(parts) > 3 else "",
    }


def parse_ncbi_geo(geo_str: str) -> dict:
    """Parse NCBI Geo_Location field: 'Country: Region'."""
    parts = geo_str.split(":", 1)
    return {
        "country": parts[0].strip() if parts else "Unknown",
        "region": parts[1].strip() if len(parts) > 1 else "Unknown",
    }


def normalize_row(row: dict, fmt: str, column_map: dict) -> dict:
    """Map source columns to internal pipeline columns."""
    out = {}

    # Map known columns
    for src_col, dst_col in column_map.items():
        if src_col in row:
            out[dst_col] = row[src_col]

    # Handle special location parsing
    if fmt in ["gisaid_raw", "gisaid_epiflu"] and "_location_raw" in out:
        loc = parse_gisaid_location(out.pop("_location_raw"))
        out.update(loc)

    if fmt == "ncbi_csv" and "_geo_location_raw" in out:
        geo = parse_ncbi_geo(out.pop("_geo_location_raw"))
        out.update(geo)

    # Pass-through: add any columns not in the map (Tier 3)
    mapped_src_cols = set(column_map.keys())
    for col, val in row.items():
        if col not in mapped_src_cols:
            # Use original column name, sanitize to lowercase
            clean_col = col.strip().lower().replace(" ", "_")
            if clean_col not in out:
                out[clean_col] = val

    return out


def sanitize_value(value: str) -> str:
    """Clean a single metadata value."""
    if value is None:
        return ""
    v = str(value).strip()
    # GISAID uses '?' for unknown values
    if v in ("?", "unknown", "Unknown", "N/A", "n/a", "NA", "na", ""):
        return ""
    return v


def main():
    parser = argparse.ArgumentParser(description="Metadata Sanitizer")
    parser.add_argument("--input", required=True, help="Input metadata file (TSV/CSV)")
    parser.add_argument("--output", required=True, help="Output sanitized metadata (TSV)")
    parser.add_argument("--report", required=True, help="Output sanitization report (JSON)")
    parser.add_argument("--format", default="auto",
                        choices=["auto", "gisaid_tsv", "gisaid_raw", "gisaid_epiflu",
                                 "ncbi_csv", "tsv", "generic"],
                        help="Metadata format (default: auto-detect)")
    parser.add_argument("--orig-file", default="", help="Original file path (for extension detection)")
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    # --- Read input ---
    orig_ext = Path(args.orig_file).suffix.lower() if args.orig_file else input_path.suffix.lower()
    
    if orig_ext in [".xls", ".xlsx"]:
        import pandas as pd
        engine = "openpyxl" if orig_ext == ".xlsx" else "xlrd"
        print(f"[SANITIZE] Detected Excel file format ({orig_ext}). Loading via pandas (engine={engine})...")
        df = pd.read_excel(input_path, engine=engine)
        # GISAID Excel often has NaNs for empty cells, convert them to empty strings
        df = df.fillna("")
        raw_headers = list(df.columns)
        raw_rows = df.to_dict('records')
    else:
        # Detect delimiter (tab for TSV, comma for CSV)
        with open(input_path, "r", encoding="utf-8") as f:
            first_line = f.readline()
            delimiter = "\t" if "\t" in first_line else ","

        with open(input_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            raw_headers = reader.fieldnames or []
            raw_rows = list(reader)

    print(f"[SANITIZE] Read {len(raw_rows)} rows from {input_path}")
    print(f"[SANITIZE] Columns detected: {raw_headers}")

    # --- Detect format ---
    fmt = args.format
    if fmt == "auto":
        fmt = detect_format(raw_headers)
    print(f"[SANITIZE] Format detected/selected: {fmt}")

    # Select column map
    column_maps = {
        "gisaid_tsv": GISAID_COLUMN_MAP,
        "gisaid_raw": GISAID_RAW_COLUMN_MAP,
        "gisaid_epiflu": GISAID_EPIFLU_MAP,
        "ncbi_csv": NCBI_COLUMN_MAP,
        "tsv": {},      # Already in internal format
        "generic": {},   # Try as-is
    }
    column_map = column_maps.get(fmt, {})

    # --- Process rows ---
    sanitized_rows = []
    stats = {
        "total_input": len(raw_rows),
        "dropped_no_strain": 0,
        "dropped_strain_ids": [],
        "filled_date": 0,
        "filled_country": 0,
        "filled_region": 0,
        "filled_host": 0,
    }

    for i, raw_row in enumerate(raw_rows):
        # Step 1: Normalize column names
        if column_map:
            row = normalize_row(raw_row, fmt, column_map)
        else:
            # For TSV/generic: clean column names, keep values
            row = {}
            for k, v in raw_row.items():
                clean_key = k.strip().lower().replace(" ", "_") if k else f"col_{i}"
                row[clean_key] = v

        # Step 2: Sanitize all values
        row = {k: sanitize_value(v) for k, v in row.items()}

        # Step 3: Tier 1 — Check primary columns
        # strain: DROP if empty
        strain = row.get("strain", "")
        if not strain:
            stats["dropped_no_strain"] += 1
            if stats["dropped_no_strain"] <= 50:
                stats["dropped_strain_ids"].append(f"row_{i+2}")  # +2 for 1-indexed + header
            continue

        # date: FILL default if empty
        if not row.get("date", ""):
            row["date"] = TIER1_DEFAULTS["date"]
            stats["filled_date"] += 1

        # country: FILL default if empty
        if not row.get("country", ""):
            row["country"] = TIER1_DEFAULTS["country"]
            stats["filled_country"] += 1

        # Step 4: Tier 2 — Fill secondary defaults silently
        for col, default in TIER2_DEFAULTS.items():
            if not row.get(col, ""):
                row[col] = default
                stats[f"filled_{col}"] = stats.get(f"filled_{col}", 0) + 1

        sanitized_rows.append(row)

    # --- Determine output columns ---
    # Collect all unique columns across all rows, maintaining order:
    # Tier 1 first, then Tier 2, then pass-through alphabetically
    tier1_cols = ["strain", "date", "country"]
    tier2_cols = ["region", "host", "division"]
    seen_cols = set()
    passthrough_cols = []

    for row in sanitized_rows:
        for col in row.keys():
            if col not in seen_cols and col not in tier1_cols and col not in tier2_cols:
                passthrough_cols.append(col)
                seen_cols.add(col)

    output_columns = tier1_cols + tier2_cols + sorted(passthrough_cols)

    # --- Write output ---
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=output_columns, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for row in sanitized_rows:
            writer.writerow(row)

    print(f"[SANITIZE] Wrote {len(sanitized_rows)} sanitized rows to {output_path}")

    # --- Print warnings for Tier 1 fills ---
    if stats["dropped_no_strain"] > 0:
        print(f"[SANITIZE] WARNING: Dropped {stats['dropped_no_strain']} rows "
              f"with empty strain (ID)", file=sys.stderr)
    if stats["filled_date"] > 0:
        print(f"[SANITIZE] WARNING: {stats['filled_date']} rows had "
              f"empty date → filled 'XXXX-XX-XX'", file=sys.stderr)
    if stats["filled_country"] > 0:
        print(f"[SANITIZE] WARNING: {stats['filled_country']} rows had "
              f"empty country → filled 'Unknown'", file=sys.stderr)

    # --- Write report ---
    report = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "input_file": str(input_path),
        "detected_format": fmt,
        "original_columns": raw_headers,
        "output_columns": output_columns,
        "stats": stats,
        "output_count": len(sanitized_rows),
        "retention_pct": round(
            (len(sanitized_rows) / len(raw_rows)) * 100, 1
        ) if raw_rows else 0,
    }

    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"[SANITIZE] Report written to {report_path}")


if __name__ == "__main__":
    main()
