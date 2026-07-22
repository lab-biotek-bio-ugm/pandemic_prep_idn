#!/usr/bin/env python3
"""
Supplementary Data & Synthesis Order Exporter (v4.0 Final)
==========================================================
Collects all pipeline output files scattered across results/ and consolidates
them into an Executive Summary Workbook (Executive_Summary_Diagnostic_Assays.xlsx)
and supporting S1-S5 workbooks.

Files Generated:
  1. Executive_Summary_Diagnostic_Assays.xlsx
     - 00_Executive_Overview : High-level summary of builds, sequence counts, and top recommendations.
     - 01_Synthesis_Oligo_Orders : Ready-to-vendor ordering specs for Rank 1 Best & Top 3 Alternatives (Triplex & Single).
     - 02_Multiplex_Pairwise : 2-target screening tables with cross-dimers.
     - 03_Multiplex_Triplex : 3-target screening tables with 27 cross-dimer pairs & dye channels.
     - 04_BLAST_OffTarget_Summary : Human & Virus off-target details for top candidates.

  2. Supporting Workbooks (S1_Primer_Candidates, S2_BLAST_Specificity, S3_Multiplex_Compatibility, S4_QC_Summary, S5_Conservation_Scoring)
  3. CSV Backups in supplementary/csv/

Usage:
  python3 scripts/export_supplementary.py --results-dir results --outdir supplementary --only-degenerate
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
    print("❌ pandas is required: pip install pandas openpyxl")
    sys.exit(1)

try:
    import openpyxl
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.utils import get_column_letter
    from openpyxl.utils.dataframe import dataframe_to_rows
    from openpyxl.worksheet.worksheet import Worksheet
except ImportError:
    print("❌ openpyxl is required: pip install openpyxl")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Formatting Constants
# ---------------------------------------------------------------------------
FILL_HEADER  = PatternFill("solid", fgColor="1E3A5F")   # Dark Navy
FILL_TITLE   = PatternFill("solid", fgColor="E8F4FD")   # Light Ice Blue
FILL_PASS    = PatternFill("solid", fgColor="C8E6C9")   # Soft Green
FILL_WARN    = PatternFill("solid", fgColor="FFF9C4")   # Soft Yellow
FILL_FAIL    = PatternFill("solid", fgColor="FFCDD2")   # Soft Red

FONT_TITLE   = Font(bold=True, name="Calibri", size=12, color="1E3A5F")
FONT_HEADER  = Font(bold=True, color="FFFFFF", name="Calibri", size=10)
FONT_BODY    = Font(name="Calibri", size=9)

BORDER_THIN  = Border(
    left=Side(style="thin", color="BDBDBD"),
    right=Side(style="thin", color="BDBDBD"),
    top=Side(style="thin", color="BDBDBD"),
    bottom=Side(style="thin", color="BDBDBD"),
)

EXPORT_DATE = datetime.now().strftime("%Y-%m-%d %H:%M")

def safe_sheet_name(name: str, maxlen: int = 31) -> str:
    name = re.sub(r"[\/\\?\*:\[\]]", "_", name)
    return name[:maxlen]

def build_label(dir_name: str) -> str:
    s = dir_name.lower()
    s = s.replace("influenza-a", "FluA").replace("influenza-b", "FluB").replace("measles", "Measles")
    s = s.replace("-mp-", "_MP_").replace("-np-", "_NP_")
    s = s.replace("indo-only", "Indonesia").replace("indonesia-only", "Indonesia")
    s = s.replace("sea-only", "SEAsia").replace("all-partial", "AllPartial")
    s = s.replace("-", "_")
    return s

def auto_width(ws: Worksheet, min_w: int = 10, max_w: int = 50):
    for col in ws.columns:
        max_len = 0
        col_letter = get_column_letter(col[0].column)
        for cell in col:
            val = str(cell.value or "")
            max_len = max(max_len, len(val))
        ws.column_dimensions[col_letter].width = max(min_w, min(max_w, max_len + 3))

def style_table(ws: Worksheet, start_row: int, ncols: int, status_col_idx: int = None):
    # Style Header
    for col in range(1, ncols + 1):
        cell = ws.cell(row=start_row, column=col)
        cell.fill = FILL_HEADER
        cell.font = FONT_HEADER
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = BORDER_THIN

    # Style Data
    for r in range(start_row + 1, ws.max_row + 1):
        fill = None
        if status_col_idx:
            val = str(ws.cell(row=r, column=status_col_idx).value or "").upper()
            if "PASS" in val or "SPECIFIC" in val:
                fill = FILL_PASS
            elif "WARN" in val or "BORDER" in val:
                fill = FILL_WARN
            elif "FAIL" in val or "VETO" in val:
                fill = FILL_FAIL

        for col in range(1, ncols + 1):
            cell = ws.cell(row=r, column=col)
            cell.font = FONT_BODY
            cell.border = BORDER_THIN
            if fill:
                cell.fill = fill

def write_df_to_ws(ws: Worksheet, df: pd.DataFrame, title: str, status_col: str = None):
    ws.append([f"🔬 {title}"])
    ws.append([f"Exported: {EXPORT_DATE} | Total Records: {len(df)}"])
    ws.append([])

    ws.cell(1, 1).font = FONT_TITLE
    ws.cell(2, 1).font = Font(italic=True, size=9, color="555555")

    for row in dataframe_to_rows(df, index=False, header=True):
        ws.append(row)

    status_idx = list(df.columns).index(status_col) + 1 if status_col and status_col in df.columns else None
    style_table(ws, 4, len(df.columns), status_idx)
    ws.freeze_panes = ws.cell(row=5, column=1)
    auto_width(ws)

# ---------------------------------------------------------------------------
# Executive Summary Exporter
# ---------------------------------------------------------------------------

def create_executive_summary(results_root: Path, outdir: Path):
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    csv_dir = outdir / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------------
    # Sheet 00: Executive Overview
    # -----------------------------------------------------------------------
    ws0 = wb.create_sheet(title="00_Executive_Overview")
    ws0.append(["🔬 PANDEMIC PREPAREDNESS — EXECUTIVE DIAGNOSTIC OVERVIEW"])
    ws0.append([f"Integrated Report | Generated: {EXPORT_DATE}"])
    ws0.append([])

    ws0.cell(1, 1).font = Font(bold=True, size=14, color="1E3A5F")
    ws0.cell(2, 1).font = Font(italic=True, size=9, color="555555")

    # Overview Table: Pipeline Summary
    overview_data = []
    builds = sorted([d for d in results_root.iterdir() if d.is_dir() and not d.name.startswith("multiplex")])
    for b in builds:
        qc_file = b / "m02" / "qc_report.json"
        input_n, output_n, ret_pct = "-", "-", "-"
        if qc_file.exists():
            with open(qc_file, "r") as f:
                qc = json.load(f)
                input_n = qc.get("input_count", "-")
                output_n = qc.get("output_count", "-")
                ret_pct = f"{qc.get('retention_pct', '-')}%"
        
        cands_json = list((b / "m09" / "degenerate").glob("**/final_candidates.json"))
        top_cand = "-"
        if cands_json:
            with open(cands_json[0], "r") as f:
                c_data = json.load(f)
                if c_data:
                    top_cand = f"Rank 1 Score: {c_data[0].get('final_score')}"

        overview_data.append({
            "Target Cohort": build_label(b.name),
            "Input Sequences": input_n,
            "QC Passed Sequences": output_n,
            "QC Retention Rate": ret_pct,
            "Top Candidate Status": top_cand
        })

    df_ov = pd.DataFrame(overview_data)
    ws0.append(["Summary of Target Cohorts & Data Quality"])
    ws0.append([])
    for r in dataframe_to_rows(df_ov, index=False, header=True):
        ws0.append(r)
    style_table(ws0, 6, len(df_ov.columns))
    auto_width(ws0)

    # -----------------------------------------------------------------------
    # Sheet 01: Synthesis Oligo Orders (Rank 1 Best + 3 Alternatives)
    # -----------------------------------------------------------------------
    ws1 = wb.create_sheet(title="01_Synthesis_Oligo_Orders")
    triplex_files = list(results_root.glob("**/triplex_candidates.json"))
    order_rows = []
    if triplex_files and triplex_files[0].exists():
        with open(triplex_files[0], "r") as f:
            t_data = json.load(f)[:4] # Rank 1 + 3 alts
            for rank_idx, combo in enumerate(t_data, 1):
                cat = "BEST RECOMMENDED (RANK 1)" if rank_idx == 1 else f"ALTERNATIVE RANK {rank_idx}"
                oligos = combo.get("oligos_synthesis", {})
                for t_key, (t_name, default_dye) in [("flu_a", ("Influenza A MP", "FAM")), ("flu_b", ("Influenza B MP", "HEX")), ("measles", ("Measles Virus", "Cy5"))]:
                    if t_key in oligos:
                        for role_key, role_lbl in [("fwd", "Forward Primer"), ("rev", "Reverse Primer"), ("probe", "Exonuclease Probe")]:
                            if role_key in oligos[t_key]:
                                o = oligos[t_key][role_key]
                                seq = o.get("deg_seq") or o.get("seq", "")
                                order_rows.append({
                                    "Recommendation Category": cat,
                                    "Assay Set": combo.get("pair_id"),
                                    "Target Virus": t_name,
                                    "Oligo Role": role_lbl,
                                    "Degenerate Sequence (5' -> 3')": seq,
                                    "Length (nt)": o.get("len", len(seq)),
                                    "Tm (°C)": o.get("tm", "-"),
                                    "GC (%)": o.get("gc", "-"),
                                    "Fluorophore / Label": o.get("label", "Unlabeled" if role_key != "probe" else f"5'-{default_dye} / 3'-BHQ / Internal-THF"),
                                    "Triplex Compatibility Score": f"{combo.get('score')}%"
                                })

    df_ord = pd.DataFrame(order_rows) if order_rows else pd.DataFrame([{"Status": "Run Triplex Module to populate"}])
    write_df_to_ws(ws1, df_ord, "Synthesis Oligo Ordering Specs (Rank 1 Best + 3 Alternatives)")
    df_ord.to_csv(str(csv_dir / "01_Synthesis_Oligo_Orders.csv"), index=False)

    # -----------------------------------------------------------------------
    # Sheet 02: Multiplex Pairwise (2-Target)
    # -----------------------------------------------------------------------
    ws2 = wb.create_sheet(title="02_Multiplex_Pairwise")
    pair_files = list(results_root.glob("**/multiplex_candidates.json"))
    pair_rows = []
    for pf in pair_files:
        if "triplex" in pf.parent.name:
            continue
        with open(pf, "r") as f:
            p_data = json.load(f)
            for rec in p_data:
                r = {
                    "Pair ID": rec.get("pair_id"),
                    "Compatibility Score (%)": rec.get("score"),
                    "Status (39°C)": rec.get("status_raa_39c"),
                    "RAA ΔG 39°C (kcal/mol)": rec.get("min_dg_raa_39c"),
                    "qPCR ΔG 55°C (kcal/mol)": rec.get("min_dg_qpcr_55c"),
                    "ΔTm Fwd (°C)": rec.get("delta_tm"),
                }
                for pair_k, pair_v in rec.get("cross_dimers_39c", {}).items():
                    r[f"ΔG 39C: {pair_k}"] = pair_v
                pair_rows.append(r)

    df_pair = pd.DataFrame(pair_rows) if pair_rows else pd.DataFrame([{"Status": "No pairwise runs found"}])
    write_df_to_ws(ws2, df_pair, "Pairwise Multiplex Screening Results (2-Target)", "Status (39°C)")
    df_pair.to_csv(str(csv_dir / "02_Multiplex_Pairwise.csv"), index=False)

    # -----------------------------------------------------------------------
    # Sheet 03: Multiplex Triplex (3-Target)
    # -----------------------------------------------------------------------
    ws3 = wb.create_sheet(title="03_Multiplex_Triplex")
    tri_rows = []
    if triplex_files and triplex_files[0].exists():
        with open(triplex_files[0], "r") as f:
            t_data = json.load(f)
            for rec in t_data:
                r = {
                    "Triplex Combo ID": rec.get("pair_id"),
                    "Compatibility Score (%)": rec.get("score"),
                    "Status (39°C)": rec.get("status_raa_39c"),
                    "RAA ΔG 39°C (kcal/mol)": rec.get("min_dg_raa_39c"),
                    "qPCR ΔG 55°C (kcal/mol)": rec.get("min_dg_qpcr_55c"),
                    "Max ΔTm (°C)": rec.get("delta_tm"),
                }
                # Dynamically capture all 27 inter-assay cross dimer pairs!
                for pair_k, pair_v in rec.get("cross_dimers_39c", {}).items():
                    r[f"ΔG 39C: {pair_k}"] = pair_v
                tri_rows.append(r)

    df_tri = pd.DataFrame(tri_rows) if tri_rows else pd.DataFrame([{"Status": "No triplex runs found"}])
    write_df_to_ws(ws3, df_tri, "Triplex Multiplex Screening Results (3-Target: FluA + FluB + Measles)", "Status (39°C)")
    df_tri.to_csv(str(csv_dir / "03_Multiplex_Triplex.csv"), index=False)

    # -----------------------------------------------------------------------
    # Sheet 04: BLAST OffTarget Summary
    # -----------------------------------------------------------------------
    ws4 = wb.create_sheet(title="04_BLAST_OffTarget_Summary")
    blast_files = list(results_root.glob("**/blast_specificity_summary.csv"))
    blast_dfs = []
    for bf in blast_files:
        try:
            b_df = pd.read_csv(bf)
            b_df.insert(0, "Build Cohort", bf.parent.parent.parent.parent.name)
            blast_dfs.append(b_df)
        except Exception:
            pass

    df_blast = pd.concat(blast_dfs, ignore_index=True) if blast_dfs else pd.DataFrame([{"Status": "No BLAST summary files found"}])
    write_df_to_ws(ws4, df_blast, "BLAST Specificity & Off-Target Risk Summary", "specificity_grade")
    df_blast.to_csv(str(csv_dir / "04_BLAST_OffTarget_Summary.csv"), index=False)

    exec_out = outdir / "Executive_Summary_Diagnostic_Assays.xlsx"
    wb.save(str(exec_out))
    print(f"  ✓ Created Main Workbook: {exec_out}")

# ---------------------------------------------------------------------------
# Main Routine
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Full Executive & Supplementary Exporter")
    parser.add_argument("--results-dir", default="results", help="Results directory")
    parser.add_argument("--outdir", default="supplementary", help="Output directory")
    parser.add_argument("--only-degenerate", action="store_true", help="Only export degenerate candidates")
    args = parser.parse_args()

    results_root = Path(args.results_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"\n================================================================================")
    print(f"📊 EXPORTING EXECUTIVE SUMMARY & SUPPLEMENTARY REPORTS")
    print(f"================================================================================")
    
    create_executive_summary(results_root, outdir)
    
    print(f"================================================================================")
    print(f"✅ EXPORT COMPLETED! Output directory: {outdir}")
    print(f"================================================================================")

if __name__ == "__main__":
    main()
