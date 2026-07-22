#!/usr/bin/env python3
"""
Supplementary Data Exporter
============================
Collects all pipeline output files scattered across the results/ directory
and consolidates them into a single supplementary/ folder with:
  - XLSX files (multi-sheet, colour-coded, freeze-panes) for supervisor review
  - CSV backup for every sheet (in supplementary/csv/)
  - README.md summarising all files and column descriptions

Usage
-----
python scripts/export_supplementary.py \\
    --results-dir results/ \\
    --outdir supplementary/

# Specify an external/Windows results path (via WSL):
python scripts/export_supplementary.py \\
    --results-dir /mnt/c/Users/engki/Downloads/results \\
    --outdir supplementary/

# Only export degenerate primers; take top-N candidates per build:
python scripts/export_supplementary.py \\
    --results-dir /mnt/c/Users/engki/Downloads/results \\
    --outdir supplementary/ \\
    --only-degenerate --top-n 20
"""

import argparse
import csv
import json
import os
import re
import sys
from datetime import datetime
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    print("❌ pandas is required:  pip install pandas openpyxl")
    sys.exit(1)

try:
    import openpyxl
    from openpyxl.styles import (
        Alignment, Font, PatternFill, Border, Side, numbers
    )
    from openpyxl.utils import get_column_letter
    from openpyxl.utils.dataframe import dataframe_to_rows
    from openpyxl.worksheet.worksheet import Worksheet
except ImportError:
    print("❌ openpyxl is required:  pip install openpyxl")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Style constants
# ---------------------------------------------------------------------------
FILL_HEADER  = PatternFill("solid", fgColor="1E3A5F")   # dark navy
FILL_PASS    = PatternFill("solid", fgColor="C8E6C9")   # soft green
FILL_WARN    = PatternFill("solid", fgColor="FFF9C4")   # soft yellow
FILL_FAIL    = PatternFill("solid", fgColor="FFCDD2")   # soft red
FILL_SUBHEAD = PatternFill("solid", fgColor="E3F2FD")   # light blue
FONT_HEADER  = Font(bold=True, color="FFFFFF", name="Calibri", size=10)
FONT_SUBHEAD = Font(bold=True, color="1E3A5F", name="Calibri", size=9)
FONT_BODY    = Font(name="Calibri", size=9)
FONT_TITLE   = Font(bold=True, name="Calibri", size=11)
BORDER_THIN  = Border(
    left=Side(border_style="thin", color="BDBDBD"),
    right=Side(border_style="thin", color="BDBDBD"),
    top=Side(border_style="thin", color="BDBDBD"),
    bottom=Side(border_style="thin", color="BDBDBD"),
)
EXPORT_DATE  = datetime.now().strftime("%Y-%m-%d %H:%M")

# ---------------------------------------------------------------------------
# Column mappings — technical key → friendly label (supervisor-readable)
# ---------------------------------------------------------------------------
CANDIDATE_COLUMNS = {
    "rank":                "Rank",
    "status":              "Status",
    "final_score":         "Final Score",
    "tier":                "Tier",
    "amplicon_start":      "Amplicon Start (bp)",
    "amplicon_end":        "Amplicon End (bp)",
    "amplicon_length":     "Amplicon Length (bp)",
    "fwd_degenerate_seq":  "Forward Primer (IUPAC)",
    "rev_degenerate_seq":  "Reverse Primer (IUPAC)",
    "probe_degenerate_seq": "Probe (IUPAC)",
    "fwd_tm":              "Fwd Tm (°C)",
    "rev_tm":              "Rev Tm (°C)",
    "probe_tm":            "Probe Tm (°C)",
    "fwd_gc":              "Fwd GC (%)",
    "rev_gc":              "Rev GC (%)",
    "probe_gc":            "Probe GC (%)",
    "hairpin_dg":          "Hairpin ΔG (kcal/mol)",
    "cross_dimer_dg":      "Cross-Dimer ΔG (kcal/mol)",
    "dimer_dg":            "Self-Dimer ΔG (kcal/mol)",
    "tm_delta":            "ΔTm Fwd-Rev (°C)",
    "human_score":         "Human Specificity Score",
    "human_grade":         "Human Specificity Grade",
    "human_status":        "Human Specificity Status",
    "virus_score":         "Virus Specificity Score",
    "virus_grade":         "Virus Specificity Grade",
    "virus_status":        "Virus Specificity Status",
    "remote_score":        "Remote Specificity Score",
    "remote_grade":        "Remote Specificity Grade",
    "gard_near":           "Near GARD Breakpoint?",
}

MULTIPLEX_COLUMNS = {
    "pair_id":         "Candidate Pair",
    "rank_a":          "Rank Target A",
    "rank_b":          "Rank Target B",
    "score":           "Compatibility Score (%)",
    "status_raa_39c":  "Status (RAA 39°C)",
    "min_dg_raa_39c":  "Worst RAA ΔG 39°C (kcal/mol)",
    "min_dg_qpcr_55c": "Worst qPCR ΔG 55°C (kcal/mol)",
    "delta_tm":        "ΔTm Fwd Primers (°C)",
    "dye_a":           "Fluorophore Target A",
    "dye_b":           "Fluorophore Target B",
    "dye_status":      "Dye Separation",
}

BLAST_COLUMNS = {
    "rank":              "Rank",
    "candidate_id":      "Candidate ID",
    "target":            "Target",
    "specificity_grade": "Specificity Grade",
    "human_hits":        "Human Genome Hits",
    "virus_hits":        "Pathogen Hits",
    "off_targets":       "Off-Targets",
    "score":             "Specificity Score",
}

# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _safe_sheet_name(name: str, maxlen: int = 31) -> str:
    """Sanitize string into a valid Excel sheet name (<= 31 chars, no /\\?*:[])."""
    name = re.sub(r"[\/\\?\*:\[\]]", "_", name)
    return name[:maxlen]


def _build_label(build_dir: str) -> str:
    """Convert a build directory name to a clean label.
    e.g. 'influenza-a-mp-indo-only' → 'FluA_MP_Indonesia'
    """
    s = build_dir.lower()
    s = s.replace("influenza-a", "FluA").replace("influenza-b", "FluB") \
         .replace("measles", "Measles")
    s = s.replace("-mp-", "_MP_").replace("-np-", "_NP_")
    s = s.replace("indo-only", "Indonesia").replace("indonesia-only", "Indonesia")
    s = s.replace("sea-only", "SEAsia").replace("all-partial", "AllPartial")
    s = s.replace("-", "_")
    return s


def _style_header_row(ws: Worksheet, row_idx: int, ncols: int):
    """Apply header styling to a given row."""
    for col in range(1, ncols + 1):
        cell = ws.cell(row=row_idx, column=col)
        cell.fill = FILL_HEADER
        cell.font = FONT_HEADER
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = BORDER_THIN


def _style_body_row(ws: Worksheet, row_idx: int, ncols: int, status_col_idx: int = None):
    """Apply body row styling, optionally colour-coding by status."""
    # Determine fill based on status value (if applicable)
    fill = None
    if status_col_idx:
        cell_val = str(ws.cell(row=row_idx, column=status_col_idx).value or "").upper()
        if "PASS" in cell_val:
            fill = FILL_PASS
        elif "WARN" in cell_val or "BORDER" in cell_val:
            fill = FILL_WARN
        elif "FAIL" in cell_val or "VETO" in cell_val:
            fill = FILL_FAIL

    for col in range(1, ncols + 1):
        cell = ws.cell(row=row_idx, column=col)
        cell.font = FONT_BODY
        cell.border = BORDER_THIN
        cell.alignment = Alignment(wrap_text=False, vertical="center")
        if fill:
            cell.fill = fill


def _auto_width(ws: Worksheet, min_w: int = 8, max_w: int = 45):
    """Set column widths based on content."""
    for col in ws.columns:
        max_len = 0
        col_letter = get_column_letter(col[0].column)
        for cell in col:
            try:
                max_len = max(max_len, len(str(cell.value or "")))
            except Exception:
                pass
        ws.column_dimensions[col_letter].width = max(min_w, min(max_w, max_len + 2))


def _add_metadata_header(ws: Worksheet, title: str, build: str, count: int, ncols: int):
    """Insert two metadata rows at top of sheet before data."""
    ws.insert_rows(1, amount=3)
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=max(ncols, 4))
    ws.merge_cells(start_row=2, start_column=1, end_row=2, end_column=max(ncols, 4))
    ws.merge_cells(start_row=3, start_column=1, end_row=3, end_column=max(ncols, 4))

    ws.cell(1, 1).value = f"🔬 {title}"
    ws.cell(1, 1).font = FONT_TITLE
    ws.cell(1, 1).fill = PatternFill("solid", fgColor="E8F4FD")

    ws.cell(2, 1).value = f"Build: {build}    |    Total records: {count}    |    Exported: {EXPORT_DATE}"
    ws.cell(2, 1).font = Font(italic=True, name="Calibri", size=9, color="555555")

    ws.cell(3, 1).value = ""
    ws.row_dimensions[1].height = 18
    ws.row_dimensions[2].height = 14


def write_df_to_sheet(
    wb: openpyxl.Workbook,
    df: pd.DataFrame,
    sheet_name: str,
    title: str,
    build: str,
    status_col: str = None,
) -> None:
    """Write a DataFrame to a new sheet with full styling."""
    ws = wb.create_sheet(title=_safe_sheet_name(sheet_name))

    # Write data starting at row 1
    for r_idx, row in enumerate(dataframe_to_rows(df, index=False, header=True), 1):
        for c_idx, value in enumerate(row, 1):
            ws.cell(r_idx, c_idx, value)

    ncols = len(df.columns)
    status_col_idx = None
    if status_col and status_col in df.columns:
        status_col_idx = list(df.columns).index(status_col) + 1

    # Style header row (row 1)
    _style_header_row(ws, 1, ncols)
    # Style body rows
    for r in range(2, len(df) + 2):
        _style_body_row(ws, r, ncols, status_col_idx)

    # Insert metadata above data
    _add_metadata_header(ws, title, build, len(df), ncols)

    # Freeze pane below metadata + header
    ws.freeze_panes = ws.cell(row=6, column=1)

    _auto_width(ws)


# ---------------------------------------------------------------------------
# S1 — Primer Candidates
# ---------------------------------------------------------------------------

def export_primer_candidates(
    results_root: Path,
    outdir: Path,
    only_degenerate: bool,
    top_n: int,
) -> openpyxl.Workbook:
    """Build S1_Primer_Candidates.xlsx from all final_candidates.json files."""
    wb = openpyxl.Workbook()
    wb.remove(wb.active)  # remove default sheet
    csv_dir = outdir / "csv"
    csv_dir.mkdir(exist_ok=True)

    modes = ["degenerate"] if only_degenerate else ["degenerate", "non_degenerate"]

    builds = sorted([d for d in results_root.iterdir() if d.is_dir() and not d.name.startswith("multiplex")])
    sheet_count = 0

    for build_dir in builds:
        build_label = _build_label(build_dir.name)
        m09_root = build_dir / "m09"
        if not m09_root.exists():
            continue

        for mode in modes:
            mode_dir = m09_root / mode
            if not mode_dir.exists():
                continue

            for region_dir in sorted(mode_dir.iterdir()):
                json_path = region_dir / "final_candidates.json"
                if not json_path.exists():
                    continue

                with open(json_path, "r") as f:
                    data = json.load(f)

                if top_n:
                    data = data[:top_n]

                if not data:
                    continue

                df = pd.DataFrame(data)
                # Keep only columns that exist in the file
                cols_to_keep = [c for c in CANDIDATE_COLUMNS if c in df.columns]
                df = df[cols_to_keep].rename(columns=CANDIDATE_COLUMNS)

                region_name = region_dir.name
                mode_suffix = "Deg" if mode == "degenerate" else "NonDeg"
                sheet_name = f"{build_label}_{region_name}_{mode_suffix}"

                write_df_to_sheet(
                    wb, df,
                    sheet_name=sheet_name,
                    title=f"Primer Candidates — {build_dir.name} / {region_name} [{mode}]",
                    build=build_dir.name,
                    status_col="Status",
                )

                # CSV backup
                csv_path = csv_dir / f"S1_{sheet_name}.csv"
                df.to_csv(str(csv_path), index=False)
                sheet_count += 1
                print(f"  ✓ S1 sheet: {sheet_name} ({len(df)} candidates)")

    print(f"  → S1: {sheet_count} sheets created")
    return wb


# ---------------------------------------------------------------------------
# S2 — BLAST Specificity
# ---------------------------------------------------------------------------

def export_blast_specificity(
    results_root: Path,
    outdir: Path,
    only_degenerate: bool,
) -> openpyxl.Workbook:
    """Build S2_BLAST_Specificity.xlsx from all blast_specificity_summary.csv files."""
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    csv_dir = outdir / "csv"
    csv_dir.mkdir(exist_ok=True)

    modes = ["degenerate"] if only_degenerate else ["degenerate", "non_degenerate"]
    builds = sorted([d for d in results_root.iterdir() if d.is_dir() and not d.name.startswith("multiplex")])
    sheet_count = 0

    for build_dir in builds:
        build_label = _build_label(build_dir.name)
        m09_root = build_dir / "m09"
        if not m09_root.exists():
            continue

        for mode in modes:
            mode_dir = m09_root / mode
            if not mode_dir.exists():
                continue

            for region_dir in sorted(mode_dir.iterdir()):
                csv_path = region_dir / "blast_specificity_summary.csv"
                if not csv_path.exists():
                    continue

                df = pd.read_csv(str(csv_path))
                region_name = region_dir.name
                mode_suffix = "Deg" if mode == "degenerate" else "NonDeg"
                sheet_name = f"{build_label}_{region_name}_{mode_suffix}"

                write_df_to_sheet(
                    wb, df,
                    sheet_name=sheet_name,
                    title=f"BLAST Specificity — {build_dir.name} / {region_name} [{mode}]",
                    build=build_dir.name,
                    status_col=None,
                )

                # CSV backup
                out_csv = csv_dir / f"S2_{sheet_name}.csv"
                df.to_csv(str(out_csv), index=False)
                sheet_count += 1
                print(f"  ✓ S2 sheet: {sheet_name} ({len(df)} rows)")

    print(f"  → S2: {sheet_count} sheets created")
    return wb


# ---------------------------------------------------------------------------
# S3 — Multiplex Compatibility
# ---------------------------------------------------------------------------

def _flatten_cross_dimers(record: dict, temp_key: str) -> dict:
    """Flatten cross_dimers_39c / cross_dimers_55c sub-dict into flat columns."""
    flat = {}
    dimers = record.get(temp_key, {})
    suffix = "39C" if "39" in temp_key else "55C"
    for pair_name, dg in dimers.items():
        flat[f"ΔG {pair_name} ({suffix})"] = dg
    return flat


def export_multiplex(results_root: Path, outdir: Path) -> openpyxl.Workbook:
    """Build S3_Multiplex_Compatibility.xlsx from multiplex_candidates.json files."""
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    csv_dir = outdir / "csv"
    csv_dir.mkdir(exist_ok=True)

    multiplex_dirs = sorted([
        d for d in results_root.iterdir()
        if d.is_dir() and d.name.startswith("multiplex")
    ])
    sheet_count = 0

    for mux_dir in multiplex_dirs:
        json_path = mux_dir / "multiplex_candidates.json"
        triplex_path = mux_dir / "triplex_candidates.json"

        target_path = triplex_path if triplex_path.exists() else json_path
        if not target_path.exists():
            continue

        with open(target_path, "r") as f:
            data = json.load(f)

        if not data:
            continue

        rows = []
        for rec in data:
            row = {}
            # Main fields
            for key, label in MULTIPLEX_COLUMNS.items():
                if key in rec:
                    row[label] = rec[key]
            # Flatten cross-dimer tables
            row.update(_flatten_cross_dimers(rec, "cross_dimers_39c"))
            row.update(_flatten_cross_dimers(rec, "cross_dimers_55c"))
            rows.append(row)

        df = pd.DataFrame(rows)
        is_triplex = triplex_path.exists()
        prefix = "Triplex" if is_triplex else "Pairwise"
        # Shorten name by removing "multiplex-" prefix
        short_name = mux_dir.name.replace("multiplex-", "").replace("multiplex_", "")
        sheet_name = f"{prefix}_{short_name}"

        write_df_to_sheet(
            wb, df,
            sheet_name=sheet_name,
            title=f"Multiplex Compatibility — {mux_dir.name}",
            build=mux_dir.name,
            status_col="Status (RAA 39°C)",
        )

        out_csv = csv_dir / f"S3_{sheet_name}.csv"
        df.to_csv(str(out_csv), index=False)
        sheet_count += 1
        print(f"  ✓ S3 sheet: {sheet_name} ({len(df)} pairs)")

    print(f"  → S3: {sheet_count} sheets created")
    return wb


# ---------------------------------------------------------------------------
# S4 — QC Summary
# ---------------------------------------------------------------------------

def export_qc_summary(results_root: Path, outdir: Path) -> openpyxl.Workbook:
    """Build S4_QC_Summary.xlsx from all qc_report.json files."""
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    csv_dir = outdir / "csv"
    csv_dir.mkdir(exist_ok=True)

    rows = []
    builds = sorted([d for d in results_root.iterdir() if d.is_dir() and not d.name.startswith("multiplex")])

    for build_dir in builds:
        qc_path = build_dir / "m02" / "qc_report.json"
        if not qc_path.exists():
            continue

        with open(qc_path, "r") as f:
            qc = json.load(f)

        row = {
            "Build": build_dir.name,
            "Export Timestamp": qc.get("timestamp", ""),
            "Input Sequences": qc.get("input_count", ""),
            "Output Sequences (after QC)": qc.get("output_count", ""),
            "Retention Rate (%)": qc.get("retention_pct", ""),
        }
        
        # Safely convert dict/lists to string representation for Excel compatibility
        for k in ["subtype_filtering", "filters"]:
            val = qc.get(k)
            if isinstance(val, dict):
                # Flatten the keys and values nicely
                row[k.replace("_", " ").title()] = "; ".join(f"{sub_k}: {sub_v}" for sub_k, sub_v in val.items())
            elif val is not None:
                row[k.replace("_", " ").title()] = str(val)

        rows.append(row)

    if not rows:
        print("  ⚠ S4: No qc_report.json files found")
        wb.create_sheet("QC_No_Data")
        return wb

    df = pd.DataFrame(rows)
    write_df_to_sheet(
        wb, df,
        sheet_name="QC_All_Builds",
        title="QC & Sequence Filtering Summary — All Builds",
        build="All builds",
        status_col=None,
    )

    out_csv = csv_dir / "S4_QC_All_Builds.csv"
    df.to_csv(str(out_csv), index=False)
    print(f"  ✓ S4: QC summary ({len(df)} builds)")
    return wb


# ---------------------------------------------------------------------------
# S5 — Conservation / Bridge Scoring
# ---------------------------------------------------------------------------

def export_conservation(results_root: Path, outdir: Path) -> openpyxl.Workbook:
    """Build S5_Conservation_Scoring.xlsx from candidate_targets_summary.csv files."""
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    csv_dir = outdir / "csv"
    csv_dir.mkdir(exist_ok=True)

    builds = sorted([d for d in results_root.iterdir() if d.is_dir() and not d.name.startswith("multiplex")])
    sheet_count = 0

    for build_dir in builds:
        csv_path = build_dir / "m07" / "candidate_targets_summary.csv"
        if not csv_path.exists():
            continue

        df = pd.read_csv(str(csv_path))
        build_label = _build_label(build_dir.name)
        sheet_name = f"Conservation_{build_label}"

        write_df_to_sheet(
            wb, df,
            sheet_name=sheet_name,
            title=f"Conservation & Bridge Scoring — {build_dir.name}",
            build=build_dir.name,
            status_col=None,
        )

        out_csv = csv_dir / f"S5_{sheet_name}.csv"
        df.to_csv(str(out_csv), index=False)
        sheet_count += 1
        print(f"  ✓ S5 sheet: {sheet_name} ({len(df)} regions)")

    print(f"  → S5: {sheet_count} sheets created")
    return wb


# ---------------------------------------------------------------------------
# README generator
# ---------------------------------------------------------------------------

README_TEMPLATE = """\
# Supplementary Data — Pandemic Preparedness Primer Design Pipeline
**Exported:** {date}  
**Pipeline source:** https://github.com/engkinandatama/pandemic-preparedness  
**Data source:** {results_dir}

---

## Files in This Folder

| File | Contents |
|---|---|
| `S1_Primer_Candidates.xlsx` | All primer candidate sets (degenerate IUPAC sequences, Tm, GC%, ΔG scores, specificity grades) |
| `S2_BLAST_Specificity.xlsx` | BLAST specificity verification results against human genome, known pathogens, and remote organisms |
| `S3_Multiplex_Compatibility.xlsx` | Cross-dimer thermodynamic compatibility screening for pairwise and triplex assay combinations |
| `S4_QC_Summary.xlsx` | Sequence quality control and filtering summary per analysis build |
| `S5_Conservation_Scoring.xlsx` | Shannon entropy–based conservation scoring and bridge-report region grades |
| `csv/` | CSV backup of every sheet above (for import into R, Python, or other tools) |

---

## Column Glossary

### S1 — Primer Candidates
| Column | Meaning |
|---|---|
| Rank | Priority rank within the build (lower = better) |
| Status | PASS / BORDERLINE / FAIL based on all scoring criteria |
| Final Score | Composite PrimerLab score (thermodynamics + specificity) |
| Forward / Reverse Primer (IUPAC) | Degenerate primer sequence using IUPAC ambiguity codes |
| Probe (IUPAC) | Exonuclease probe sequence |
| Fwd/Rev/Probe Tm | Melting temperature (°C) |
| Fwd/Rev/Probe GC% | GC content percentage |
| Hairpin ΔG | Self-folding free energy (kcal/mol); less negative = safer |
| Cross-Dimer ΔG | Worst-case heterodimerisation energy; ≥ −9 kcal/mol required |
| Human Specificity Grade | A (excellent) → D (poor off-target risk against human genome) |
| Virus Specificity Grade | A (excellent) → D (poor cross-reactivity against non-target viruses) |
| Near GARD Breakpoint? | True if amplicon overlaps an intra-segment recombination breakpoint |

### S3 — Multiplex Compatibility
| Column | Meaning |
|---|---|
| Compatibility Score (%) | 100% = no dimer risk; lower = worse heterodimerisation |
| Worst RAA ΔG 39°C | Most negative cross-dimer ΔG at isothermal 39°C operating temperature |
| Worst qPCR ΔG 55°C | Most negative cross-dimer ΔG at qPCR 55°C operating temperature |
| ΔTm Fwd Primers | Difference in melting temperatures between forward primers (°C); ≤ 3.5°C recommended |
| ΔG [OligoA]_vs_[OligoB] | Individual oligo-pair heterodimerisation ΔG at stated temperature |

---

## Status Colour Coding in Excel Files
- 🟩 **Green** = PASS
- 🟨 **Yellow** = BORDERLINE / WARNING
- 🟥 **Red** = FAIL / VETO
"""


def write_readme(outdir: Path, results_dir: str):
    readme_path = outdir / "README.md"
    with open(readme_path, "w") as f:
        f.write(README_TEMPLATE.format(
            date=datetime.now().strftime("%Y-%m-%d %H:%M"),
            results_dir=results_dir,
        ))
    print(f"  ✓ README.md written to {readme_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Export supplementary data from results/ to a supervisor-friendly folder"
    )
    parser.add_argument(
        "--results-dir",
        required=True,
        help="Path to the results/ directory (can be a WSL /mnt/c/... path)",
    )
    parser.add_argument(
        "--outdir",
        default="supplementary",
        help="Output folder for supplementary files (default: supplementary/)",
    )
    parser.add_argument(
        "--only-degenerate",
        action="store_true",
        help="Only export degenerate primer candidates (skip non-degenerate)",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=0,
        help="Export only the top-N candidates per build (0 = all)",
    )
    args = parser.parse_args()

    results_root = Path(args.results_dir).resolve()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "csv").mkdir(exist_ok=True)

    if not results_root.exists():
        print(f"❌ Results directory not found: {results_root}")
        sys.exit(1)

    print(f"\n📁 Source:  {results_root}")
    print(f"📁 Output:  {outdir.resolve()}")
    print(f"🔬 Mode:    {'Degenerate only' if args.only_degenerate else 'Degenerate + Non-degenerate'}")
    print(f"🔢 Top-N:   {'All' if not args.top_n else args.top_n}")
    print()

    # --- S1: Primer Candidates ---
    print("📊 Exporting S1: Primer Candidates...")
    wb1 = export_primer_candidates(results_root, outdir, args.only_degenerate, args.top_n)
    wb1.save(str(outdir / "S1_Primer_Candidates.xlsx"))
    print()

    # --- S2: BLAST Specificity ---
    print("📊 Exporting S2: BLAST Specificity...")
    wb2 = export_blast_specificity(results_root, outdir, args.only_degenerate)
    wb2.save(str(outdir / "S2_BLAST_Specificity.xlsx"))
    print()

    # --- S3: Multiplex Compatibility ---
    print("📊 Exporting S3: Multiplex Compatibility...")
    wb3 = export_multiplex(results_root, outdir)
    wb3.save(str(outdir / "S3_Multiplex_Compatibility.xlsx"))
    print()

    # --- S4: QC Summary ---
    print("📊 Exporting S4: QC Summary...")
    wb4 = export_qc_summary(results_root, outdir)
    wb4.save(str(outdir / "S4_QC_Summary.xlsx"))
    print()

    # --- S5: Conservation Scoring ---
    print("📊 Exporting S5: Conservation Scoring...")
    wb5 = export_conservation(results_root, outdir)
    wb5.save(str(outdir / "S5_Conservation_Scoring.xlsx"))
    print()

    # --- README ---
    write_readme(outdir, str(results_root))

    print("=" * 60)
    print(f"✅ Export complete! All files saved to: {outdir.resolve()}")
    print(f"   Open XLSX files with Excel or LibreOffice Calc.")
    print("=" * 60)


if __name__ == "__main__":
    main()
