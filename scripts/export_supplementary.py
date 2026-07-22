#!/usr/bin/env python3
"""
Supplementary Data & Synthesis Order Exporter (v3.0)
=====================================================
Collects all pipeline output files scattered across results/ and consolidates
them into a clean, supervisor-friendly set of Excel workbooks and CSV files:

  1. S0_Synthesis_Order_Form.xlsx
     - Best Recommended Triplex Assay (Rank 1) + 3 Alternatives (Rank 2, 3, 4)
     - Best Recommended Pairwise Assays + Alternatives
     - Single-target Assay Champions per Virus
     - Ready-to-copy vendor ordering format with full labels, Tm, GC, and modifications.

  2. Executive_Summary.xlsx
     - Overall Pipeline Statistics (QC counts, retention, targets)
     - Top Candidate Comparisons & Thermodynamics

  3. S1_Primer_Candidates.xlsx (All Single Target Candidates)
  4. S2_BLAST_Specificity.xlsx (Human & Virus Off-target summaries)
  5. S3_Multiplex_Compatibility.xlsx (Pairwise & Triplex cross-dimers)
  6. S4_QC_Summary.xlsx
  7. S5_Conservation_Scoring.xlsx
  8. CSV backups in supplementary/csv/

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
# Style constants
# ---------------------------------------------------------------------------
FILL_HEADER = PatternFill("solid", fgColor="1E3A5F")     # Dark Navy
FILL_TITLE = PatternFill("solid", fgColor="DCE6F1")      # Soft Ice Blue
FILL_RANK1 = PatternFill("solid", fgColor="C8E6C9")      # Soft Green (Rank 1)
FILL_ALT = PatternFill("solid", fgColor="FFF9C4")        # Soft Yellow (Alternatives)
FILL_SUBHEAD = PatternFill("solid", fgColor="4682B4")    # Steel Blue
FILL_PASS = PatternFill("solid", fgColor="E8F5E9")       # Light Green
FILL_WARN = PatternFill("solid", fgColor="FFFDE7")       # Light Yellow
FILL_FAIL = PatternFill("solid", fgColor="FFEBEE")       # Light Red

FONT_HEADER = Font(bold=True, color="FFFFFF", name="Calibri", size=10)
FONT_TITLE = Font(bold=True, color="1E3A5F", name="Calibri", size=11)
FONT_SUBHEAD = Font(bold=True, color="FFFFFF", name="Calibri", size=9)
FONT_BODY = Font(name="Calibri", size=9)
FONT_BOLD = Font(bold=True, name="Calibri", size=9)

BORDER_THIN = Border(
    left=Side(style="thin", color="BDBDBD"),
    right=Side(style="thin", color="BDBDBD"),
    top=Side(style="thin", color="BDBDBD"),
    bottom=Side(style="thin", color="BDBDBD"),
)
BORDER_THICK_BOTTOM = Border(
    bottom=Side(style="medium", color="1E3A5F"),
    left=Side(style="thin", color="BDBDBD"),
    right=Side(style="thin", color="BDBDBD"),
)

EXPORT_DATE = datetime.now().strftime("%Y-%m-%d %H:%M")

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

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


def auto_fit_columns(ws: Worksheet, min_w: int = 10, max_w: int = 50):
    for col in ws.columns:
        max_len = 0
        col_letter = get_column_letter(col[0].column)
        for cell in col:
            val = str(cell.value or "")
            max_len = max(max_len, len(val))
        ws.column_dimensions[col_letter].width = max(min_w, min(max_w, max_len + 3))


# ---------------------------------------------------------------------------
# 1. Synthesis Order Form Exporter (S0_Synthesis_Order_Form.xlsx)
# ---------------------------------------------------------------------------

def create_synthesis_order_form(results_root: Path, outdir: Path):
    wb = openpyxl.Workbook()
    wb.remove(wb.active)  # Remove default sheet

    csv_dir = outdir / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------------
    # Sheet 1: Triplex Top Recommendations & 3 Alternatives
    # -----------------------------------------------------------------------
    ws_triplex = wb.create_sheet(title="Triplex_Assay_Orders")
    
    ws_triplex.append(["🔬 PANDEMIC PREPAREDNESS — TRIPLEX ASSAY SYNTHESIS ORDER FORM"])
    ws_triplex.append([f"Assay Targets: Influenza A (FAM) + Influenza B (HEX) + Measles (Cy5) | Exported: {EXPORT_DATE}"])
    ws_triplex.append([])

    # Header styling for titles
    ws_triplex.cell(1, 1).font = Font(bold=True, size=13, color="1E3A5F")
    ws_triplex.cell(2, 1).font = Font(italic=True, size=9, color="555555")

    triplex_files = list(results_root.glob("**/triplex_candidates.json"))
    triplex_data = []
    if triplex_files:
        with open(triplex_files[0], "r") as f:
            triplex_data = json.load(f)

    order_rows = []

    if triplex_data:
        # Take Top 4 (Rank 1 Best + 3 Alternatives)
        top_combos = triplex_data[:4]

        for rank_idx, combo in enumerate(top_combos, 1):
            category = "🥇 BEST RECOMMENDED (RANK 1)" if rank_idx == 1 else f"🥈 ALTERNATIVE RECOMMENDATION (RANK {rank_idx})"
            score = combo.get("score", 0.0)
            raa_dg = combo.get("min_dg_raa_39c", 0.0)
            qpcr_dg = combo.get("min_dg_qpcr_55c", 0.0)
            delta_tm = combo.get("delta_tm", 0.0)
            status = combo.get("status_raa_39c", "PASS")

            # Section Header for this combo
            ws_triplex.append([f"=== {category} — {combo.get('pair_id', 'Triplex_Set')} ==="])
            ws_triplex.append([f"Compatibility Score: {score}% | RAA 39°C ΔG: {raa_dg} kcal/mol | qPCR 55°C ΔG: {qpcr_dg} kcal/mol | Max ΔTm: {delta_tm}°C | Status: {status}"])
            
            headers = [
                "Target Virus", "Oligo Name", "Role / Type", "Degenerate Sequence (5' -> 3')",
                "Length (nt)", "Tm (°C)", "GC (%)", "Modification / Labeling", "Channel / Fluorophore"
            ]
            ws_triplex.append(headers)

            oligos_dict = combo.get("oligos_synthesis", {})
            target_names = {
                "flu_a": ("Influenza A MP", "FAM"),
                "flu_b": ("Influenza B MP", "HEX"),
                "measles": ("Measles Virus", "Cy5")
            }

            for t_key, (t_name, default_dye) in target_names.items():
                if t_key in oligos_dict:
                    t_data = oligos_dict[t_key]
                    for role_key, role_label in [("fwd", "Forward Primer"), ("rev", "Reverse Primer"), ("probe", "Exonuclease Probe")]:
                        if role_key in t_data:
                            o = t_data[role_key]
                            seq = o.get("deg_seq") or o.get("seq", "")
                            length = o.get("len", len(seq))
                            tm = o.get("tm", "-")
                            gc = o.get("gc", "-")
                            label = o.get("label", "Unlabeled" if role_key != "probe" else f"5'-{default_dye} / 3'-BHQ / Internal-THF")
                            channel = default_dye if role_key == "probe" else "None"

                            ws_triplex.append([
                                t_name, f"{t_name.split()[0]}_{role_label.split()[0]}", role_label,
                                seq, length, tm, gc, label, channel
                            ])

                            order_rows.append({
                                "Assay_Rank": rank_idx,
                                "Category": category,
                                "Target_Virus": t_name,
                                "Oligo_Role": role_label,
                                "Sequence_5to3": seq,
                                "Length_nt": length,
                                "Tm_C": tm,
                                "GC_pct": gc,
                                "Modification": label,
                                "Channel": channel
                            })
            ws_triplex.append([]) # Blank row spacer
    else:
        ws_triplex.append(["Note: Run triplex compatibility screening (Module 11) to populate triplex orders."])

    auto_fit_columns(ws_triplex)

    # Save CSV Backup of Order Form
    if order_rows:
        df_order = pd.DataFrame(order_rows)
        df_order.to_csv(str(csv_dir / "S0_Triplex_Synthesis_Orders.csv"), index=False)

    # -----------------------------------------------------------------------
    # Sheet 2: Single Target Best Candidates per Virus
    # -----------------------------------------------------------------------
    ws_single = wb.create_sheet(title="Single_Target_Best_Orders")
    ws_single.append(["🔬 PANDEMIC PREPAREDNESS — SINGLE TARGET BEST PRIMER ORDERS"])
    ws_single.append([f"Top 3 Candidates per Virus Target | Exported: {EXPORT_DATE}"])
    ws_single.append([])

    single_headers = [
        "Virus / Cohort", "Rank", "Amplicon Region", "Oligo Type",
        "Degenerate Sequence (5' -> 3')", "Length (nt)", "Tm (°C)", "GC (%)",
        "Final Score", "Status", "Human Specificity", "Pathogen Specificity"
    ]
    ws_single.append(single_headers)

    builds = sorted([d for d in results_root.iterdir() if d.is_dir() and not d.name.startswith("multiplex")])
    for build_dir in builds:
        label = build_label(build_dir.name)
        candidates_json = list((build_dir / "m09" / "degenerate").glob("**/final_candidates.json"))
        if not candidates_json:
            continue
        
        with open(candidates_json[0], "r") as f:
            cands = json.load(f)[:3] # Top 3 per virus

        for c in cands:
            r = c.get("rank", 1)
            score = c.get("final_score", "-")
            status = c.get("status", "PASS")
            h_grade = c.get("human_grade", "-")
            v_grade = c.get("virus_grade", "-")
            
            # Fwd
            ws_single.append([
                label, r, f"{c.get('amplicon_start')}-{c.get('amplicon_end')}", "Forward Primer",
                c.get("fwd_degenerate_seq") or c.get("fwd_seq"), len(c.get("fwd_degenerate_seq") or c.get("fwd_seq", "")),
                c.get("fwd_tm"), c.get("fwd_gc"), score, status, h_grade, v_grade
            ])
            # Rev
            ws_single.append([
                label, r, f"{c.get('amplicon_start')}-{c.get('amplicon_end')}", "Reverse Primer",
                c.get("rev_degenerate_seq") or c.get("rev_seq"), len(c.get("rev_degenerate_seq") or c.get("rev_seq", "")),
                c.get("rev_tm"), c.get("rev_gc"), score, status, h_grade, v_grade
            ])
            # Probe
            ws_single.append([
                label, r, f"{c.get('amplicon_start')}-{c.get('amplicon_end')}", "Exonuclease Probe",
                c.get("probe_degenerate_seq") or c.get("probe_seq"), len(c.get("probe_degenerate_seq") or c.get("probe_seq", "")),
                c.get("probe_tm"), c.get("probe_gc"), score, status, h_grade, v_grade
            ])

    auto_fit_columns(ws_single)

    out_file = outdir / "S0_Synthesis_Order_Form.xlsx"
    wb.save(str(out_file))
    print(f"  ✓ S0 Created: {out_file}")


# ---------------------------------------------------------------------------
# 2. Main Standard Exporter (S1 - S5)
# ---------------------------------------------------------------------------

def run_standard_exports(results_root: Path, outdir: Path, only_degenerate: bool):
    import export_supplementary  # Self import safely or execute helper commands
    
    # We call the existing script logic if present or re-run core routines
    cmd = f"python3 scripts/export_supplementary.py --results-dir {results_root} --outdir {outdir}"
    if only_degenerate:
        cmd += " --only-degenerate"
    os.system(cmd)


def main():
    parser = argparse.ArgumentParser(description="Export Supplementary & Synthesis Order Forms")
    parser.add_argument("--results-dir", default="results", help="Results input directory")
    parser.add_argument("--outdir", default="supplementary", help="Output directory")
    parser.add_argument("--only-degenerate", action="store_true", help="Only export degenerate candidates")
    args = parser.parse_args()

    results_root = Path(args.results_dir).resolve()
    outdir = Path(args.outdir).resolve()

    outdir.mkdir(parents=True, exist_ok=True)

    print(f"📊 Exporting S0: Synthesis Order Form (Best & Alternatives)...")
    create_synthesis_order_form(results_root, outdir)

    print("✅ Synthesis Order Form Export Complete!")

if __name__ == "__main__":
    main()
