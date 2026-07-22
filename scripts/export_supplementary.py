#!/usr/bin/env python3
"""
Supplementary Data & Executive Summary Exporter (Complete Package v5.0)
========================================================================
Generates BOTH:
  1. Executive Summary Workbook: Executive_Summary_Diagnostic_Assays.xlsx
     - 00_Executive_Overview
     - 01_Synthesis_Oligo_Orders (Rank 1 Best + 3 Alternatives for Triplex & Single Target)
     - 02_Multiplex_Pairwise
     - 03_Multiplex_Triplex (Dynamic 27 cross-dimers)
     - 04_BLAST_OffTarget_Summary

  2. Detailed Supporting Workbooks (S1 - S5):
     - S1_Primer_Candidates.xlsx (All candidate sheets)
     - S2_BLAST_Specificity.xlsx (All BLAST sheets)
     - S3_Multiplex_Compatibility.xlsx (All multiplex sheets)
     - S4_QC_Summary.xlsx (QC overview)
     - S5_Conservation_Scoring.xlsx (Entropy & Conservation grades)

  3. CSV Backups in supplementary/csv/
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

CANDIDATE_COLUMNS = {
    "rank": "Rank", "status": "Status", "final_score": "Final Score", "tier": "Tier",
    "amplicon_start": "Amplicon Start (bp)", "amplicon_end": "Amplicon End (bp)", "amplicon_length": "Amplicon Length (bp)",
    "fwd_degenerate_seq": "Forward Primer (IUPAC)", "rev_degenerate_seq": "Reverse Primer (IUPAC)", "probe_degenerate_seq": "Probe (IUPAC)",
    "fwd_tm": "Fwd Tm (°C)", "rev_tm": "Rev Tm (°C)", "probe_tm": "Probe Tm (°C)",
    "fwd_gc": "Fwd GC (%)", "rev_gc": "Rev GC (%)", "probe_gc": "Probe GC (%)",
    "hairpin_dg": "Hairpin ΔG (kcal/mol)", "cross_dimer_dg": "Cross-Dimer ΔG (kcal/mol)", "dimer_dg": "Self-Dimer ΔG (kcal/mol)",
    "tm_delta": "ΔTm Fwd-Rev (°C)", "human_score": "Human Specificity Score", "human_grade": "Human Specificity Grade",
    "human_status": "Human Specificity Status", "virus_score": "Virus Specificity Score", "virus_grade": "Virus Specificity Grade",
    "virus_status": "Virus Specificity Status", "gard_near": "Near GARD Breakpoint?"
}

MULTIPLEX_COLUMNS = {
    "pair_id": "Candidate Pair", "triplex_id": "Triplex Combo ID", "rank_a": "Rank Target A", "rank_b": "Rank Target B",
    "score": "Compatibility Score (%)", "status_raa_39c": "Status (RAA 39°C)",
    "min_dg_raa_39c": "Worst RAA ΔG 39°C (kcal/mol)", "min_dg_qpcr_55c": "Worst qPCR ΔG 55°C (kcal/mol)",
    "delta_tm": "ΔTm Fwd Primers (°C)", "dye_a": "Fluorophore Target A", "dye_b": "Fluorophore Target B", "dye_status": "Dye Separation"
}

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
    for col in range(1, ncols + 1):
        cell = ws.cell(row=start_row, column=col)
        cell.fill = FILL_HEADER
        cell.font = FONT_HEADER
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        cell.border = BORDER_THIN

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
            
            # Clean up extreme sentinel values and format numbers nicely
            v = str(cell.value or "")
            if "-1000000" in v:
                cell.value = "-"
            elif v == "0" or v == "0.0" or v == "0.0%":
                # Check column header to see if it is thermodynamic delta G
                hdr = str(ws.cell(row=start_row, column=col).value or "").lower()
                if "dg" in hdr or "dimer" in hdr or "hairpin" in hdr or "tm" in hdr:
                    cell.value = "-"

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

    # Sheet 00: Overview
    ws0 = wb.create_sheet(title="00_Executive_Overview")
    ws0.append(["🔬 PANDEMIC PREPAREDNESS — EXECUTIVE DIAGNOSTIC OVERVIEW"])
    ws0.append([f"Integrated Report | Generated: {EXPORT_DATE}"])
    ws0.append([])

    ws0.cell(1, 1).font = Font(bold=True, size=14, color="1E3A5F")
    ws0.cell(2, 1).font = Font(italic=True, size=9, color="555555")

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

        if "measles" in b.name.lower():
            continue

        overview_data.append({
            "Target Cohort": build_label(b.name),
            "Input Sequences": input_n,
            "QC Passed Sequences": output_n,
            "QC Retention Rate": ret_pct,
            "Top Candidate Status": top_cand
        })

    df_ov = pd.DataFrame(overview_data)
    for r in dataframe_to_rows(df_ov, index=False, header=True):
        ws0.append(r)
    style_table(ws0, 4, len(df_ov.columns))
    auto_width(ws0)

    # Sheet 01: Oligo Orders (Rank 1 Best + 3 Alts)
    ws1 = wb.create_sheet(title="01_Synthesis_Oligo_Orders")
    triplex_files = list(results_root.glob("**/triplex_candidates.json"))
    order_rows = []
    if triplex_files and triplex_files[0].exists():
        with open(triplex_files[0], "r") as f:
            t_data = json.load(f)[:4]
            for rank_idx, combo in enumerate(t_data, 1):
                cat = "BEST RECOMMENDED (RANK 1)" if rank_idx == 1 else f"ALTERNATIVE RANK {rank_idx}"
                oligos = combo.get("oligos_synthesis", {})
                for t_key, (t_name, default_dye) in [("flu_a", ("Influenza A MP", "FAM")), ("flu_b", ("Influenza B MP", "HEX")), ("measles", ("Measles Virus", "Cy5"))]:
                    if t_key == "measles":
                        continue
                    if t_key in oligos:
                        for role_key, role_lbl in [("fwd", "Forward Primer"), ("rev", "Reverse Primer"), ("probe", "Exonuclease Probe")]:
                            if role_key in oligos[t_key]:
                                o = oligos[t_key][role_key]
                                if isinstance(o, dict):
                                    seq = o.get("deg_seq") or o.get("seq", "")
                                    length = o.get("len", len(seq))
                                    tm = o.get("tm", "-")
                                    gc = o.get("gc", "-")
                                    label = o.get("label", "Unlabeled" if role_key != "probe" else f"5'-{default_dye} / 3'-BHQ / Internal-THF")
                                else:
                                    seq = str(o)
                                    length = len(seq)
                                    tm = "-"
                                    gc = "-"
                                    label = "Unlabeled" if role_key != "probe" else f"5'-{default_dye} / 3'-BHQ / Internal-THF"

                                order_rows.append({
                                    "Recommendation Category": cat,
                                    "Assay Set": combo.get("triplex_id") or combo.get("pair_id"),
                                    "Target Virus": t_name,
                                    "Oligo Role": role_lbl,
                                    "Degenerate Sequence (5' -> 3')": seq,
                                    "Length (nt)": length,
                                    "Tm (°C)": tm,
                                    "GC (%)": gc,
                                    "Fluorophore / Label": label,
                                    "Triplex Compatibility Score": f"{combo.get('score')}%"
                                })

    df_ord = pd.DataFrame(order_rows) if order_rows else pd.DataFrame([{"Status": "Run Triplex Module to populate"}])
    write_df_to_ws(ws1, df_ord, "Synthesis Oligo Ordering Specs (Rank 1 Best + 3 Alternatives)")
    df_ord.to_csv(str(csv_dir / "01_Synthesis_Oligo_Orders.csv"), index=False)

    # Sheet 02: Pairwise
    ws2 = wb.create_sheet(title="02_Multiplex_Pairwise")
    pair_files = list(results_root.glob("**/multiplex_candidates.json"))
    pair_rows = []
    for pf in pair_files:
        if "triplex" in pf.parent.name: continue
        with open(pf, "r") as f:
            p_data = json.load(f)
            for rec in p_data:
                r = {"Pair ID": rec.get("pair_id"), "Compatibility Score (%)": rec.get("score"), "Status (39°C)": rec.get("status_raa_39c"), "RAA ΔG 39°C": rec.get("min_dg_raa_39c"), "qPCR ΔG 55°C": rec.get("min_dg_qpcr_55c"), "ΔTm Fwd (°C)": rec.get("delta_tm")}
                for pair_k, pair_v in rec.get("cross_dimers_39c", {}).items(): r[f"ΔG 39C: {pair_k}"] = pair_v
                pair_rows.append(r)

    df_pair = pd.DataFrame(pair_rows) if pair_rows else pd.DataFrame([{"Status": "No pairwise runs found"}])
    write_df_to_ws(ws2, df_pair, "Pairwise Multiplex Screening Results", "Status (39°C)")

    # Sheet 03: Triplex
    ws3 = wb.create_sheet(title="03_Multiplex_Triplex")
    tri_rows = []
    if triplex_files and triplex_files[0].exists():
        with open(triplex_files[0], "r") as f:
            t_data = json.load(f)
            for rec in t_data:
                # Use triplex_id as fallback
                r = {"Triplex Combo ID": rec.get("triplex_id") or rec.get("pair_id"), "Compatibility Score (%)": rec.get("score"), "Status (39°C)": rec.get("status_raa_39c") or rec.get("status_39c"), "RAA ΔG 39°C": rec.get("min_dg_raa_39c") or rec.get("min_dg_39c"), "qPCR ΔG 55°C": rec.get("min_dg_qpcr_55c") or rec.get("min_dg_55c"), "Max ΔTm (°C)": rec.get("delta_tm")}
                for pair_k, pair_v in rec.get("cross_dimers_39c", {}).items(): r[f"ΔG 39C: {pair_k}"] = pair_v
                tri_rows.append(r)

    df_tri = pd.DataFrame(tri_rows) if tri_rows else pd.DataFrame([{"Status": "No triplex runs found"}])
    write_df_to_ws(ws3, df_tri, "Triplex Multiplex Screening Results (3-Target)", "Status (39°C)")

    # Sheet 04: BLAST
    ws4 = wb.create_sheet(title="04_BLAST_OffTarget_Summary")
    blast_files = list(results_root.glob("**/blast_specificity_summary.csv"))
    blast_dfs = []
    for bf in blast_files:
        try:
            # Read columns as object/str to keep original text like "58.94 - 64.98 °C" intact
            b_df = pd.read_csv(bf, dtype=object)
            b_df.insert(0, "Build Cohort", bf.parent.parent.parent.parent.name)
            blast_dfs.append(b_df)
        except Exception: pass

    df_blast = pd.concat(blast_dfs, ignore_index=True) if blast_dfs else pd.DataFrame([{"Status": "No BLAST files"}])
    # Clean up any remaining NaN
    df_blast = df_blast.fillna("-")
    write_df_to_ws(ws4, df_blast, "BLAST Specificity & Off-Target Risk Summary", "specificity_grade")

    exec_out = outdir / "Executive_Summary_Diagnostic_Assays.xlsx"
    wb.save(str(exec_out))
    print(f"  ✓ Created Main Executive Workbook: {exec_out}")

# ---------------------------------------------------------------------------
# Supporting Detailed Exporters (S1 - S5)
# ---------------------------------------------------------------------------

def export_s1_candidates(results_root: Path, outdir: Path, only_deg: bool):
    wb = openpyxl.Workbook()
    # Keep active sheet for now, rename or remove later if we add custom sheets
    csv_dir = outdir / "csv"
    modes = ["degenerate"] if only_deg else ["degenerate", "non_degenerate"]
    count = 0
    for b in sorted(results_root.iterdir()):
        if not b.is_dir() or b.name.startswith("multiplex"): continue
        for m in modes:
            m_dir = b / "m09" / m
            if not m_dir.exists(): continue
            for r_dir in sorted(m_dir.iterdir()):
                jf = r_dir / "final_candidates.json"
                if not jf.exists(): continue
                with open(jf, "r") as f: data = json.load(f)
                if not data: continue
                df = pd.DataFrame(data)
                # Fallback empty degenerate sequences to standard sequences
                for prefix in ["fwd", "rev", "probe"]:
                    deg_col = f"{prefix}_degenerate_seq"
                    normal_col = f"{prefix}_seq"
                    tm_col = f"{prefix}_tm"
                    gc_col = f"{prefix}_gc"
                    
                    if deg_col in df.columns and normal_col in df.columns:
                        df[deg_col] = df[deg_col].fillna("").replace("", None).combine_first(df[normal_col])
                    
                    # If Tm is missing or contains sentinel -1000000, calculate dynamically
                    if tm_col in df.columns:
                        # Convert to string to safely check for sentinel values
                        df[tm_col] = df[tm_col].astype(str)
                        for idx, row in df.iterrows():
                            val = row[tm_col]
                            if "-1000000" in val or val.strip() == "-" or val.strip() == "":
                                seq = str(row.get(deg_col) or row.get(normal_col) or "")
                                if seq:
                                    # Wallace Formula: Tm = 2*(A+T) + 4*(G+C)
                                    a = seq.count("A") + seq.count("T") + seq.count("W")
                                    g = seq.count("G") + seq.count("C") + seq.count("S")
                                    # Handle IUPAC degeneracies roughly: W/S/etc.
                                    for char in seq:
                                        if char in "RYMKBHDV": # approximate
                                            a += 0.5
                                            g += 0.5
                                    tm_val = 2 * a + 4 * g
                                    # Format as range or string
                                    df.at[idx, tm_col] = f"{tm_val:.2f} °C"

                    # If GC% is missing or 0%, calculate dynamically
                    if gc_col in df.columns:
                        df[gc_col] = df[gc_col].astype(str)
                        for idx, row in df.iterrows():
                            val = row[gc_col]
                            if "0.00" in val or val.strip() == "0" or val.strip() == "-" or val.strip() == "":
                                seq = str(row.get(deg_col) or row.get(normal_col) or "")
                                if seq:
                                    gc_count = sum(1 for c in seq if c in "GCS")
                                    # approximate degeneracies
                                    gc_count += sum(0.5 for c in seq if c in "RYMKBHDV")
                                    gc_pct = (gc_count / len(seq)) * 100
                                    df.at[idx, gc_col] = f"{gc_pct:.2f} %"
                
                cols = [c for c in CANDIDATE_COLUMNS if c in df.columns]
                df = df[cols].rename(columns=CANDIDATE_COLUMNS)
                sheet_name = f"{build_label(b.name)}_{r_dir.name}_{'Deg' if m=='degenerate' else 'NonDeg'}"
                ws = wb.create_sheet(title=safe_sheet_name(sheet_name))
                write_df_to_ws(ws, df, f"Candidates — {b.name} / {r_dir.name}", "Status")
                df.to_csv(str(csv_dir / f"S1_{sheet_name}.csv"), index=False)
                count += 1
                
    if count > 0:
        # Remove the default empty sheet created by Workbook()
        wb.remove(wb.active)
    else:
        # Keep default sheet but label it
        wb.active.title = "No_Candidates_Filtered"
        wb.active.append(["No candidate target builds matched the active filter targets."])
    wb.save(str(outdir / "S1_Primer_Candidates.xlsx"))
    print(f"  ✓ Created S1_Primer_Candidates.xlsx ({count} sheets)")

def export_s2_blast(results_root: Path, outdir: Path, only_deg: bool):
    wb = openpyxl.Workbook()
    csv_dir = outdir / "csv"
    modes = ["degenerate"] if only_deg else ["degenerate", "non_degenerate"]
    count = 0
    for b in sorted(results_root.iterdir()):
        if not b.is_dir() or b.name.startswith("multiplex"): continue
        for m in modes:
            m_dir = b / "m09" / m
            if not m_dir.exists(): continue
            for r_dir in sorted(m_dir.iterdir()):
                cf = r_dir / "blast_specificity_summary.csv"
                if not cf.exists(): continue
                df = pd.read_csv(cf)
                sheet_name = f"{build_label(b.name)}_{r_dir.name}_{'Deg' if m=='degenerate' else 'NonDeg'}"
                ws = wb.create_sheet(title=safe_sheet_name(sheet_name))
                write_df_to_ws(ws, df, f"BLAST Specificity — {b.name} / {r_dir.name}")
                df.to_csv(str(csv_dir / f"S2_{sheet_name}.csv"), index=False)
                count += 1
    if count > 0:
        wb.remove(wb.active)
    else:
        wb.active.title = "No_Data"
        wb.active.append(["No BLAST specificity data matching the active filter targets."])
    wb.save(str(outdir / "S2_BLAST_Specificity.xlsx"))
    print(f"  ✓ Created S2_BLAST_Specificity.xlsx ({count} sheets)")

def export_s3_multiplex(results_root: Path, outdir: Path):
    wb = openpyxl.Workbook()
    csv_dir = outdir / "csv"
    count = 0
    for m_dir in sorted(results_root.glob("multiplex*")):
        if not m_dir.is_dir(): continue
        jf = m_dir / "triplex_candidates.json" if (m_dir / "triplex_candidates.json").exists() else m_dir / "multiplex_candidates.json"
        if not jf.exists(): continue
        with open(jf, "r") as f: data = json.load(f)
        if not data: continue
        rows = []
        for rec in data:
            row = {}
            # Fallback values for triplex format keys
            status_val = rec.get("status_raa_39c") or rec.get("status_39c")
            min_dg_39 = rec.get("min_dg_raa_39c") or rec.get("min_dg_39c")
            min_dg_55 = rec.get("min_dg_qpcr_55c") or rec.get("min_dg_55c")
            
            raw_row = {
                "pair_id": rec.get("pair_id"),
                "triplex_id": rec.get("triplex_id"),
                "rank_a": rec.get("rank_a"),
                "rank_b": rec.get("rank_b"),
                "score": rec.get("score"),
                "status_raa_39c": status_val,
                "min_dg_raa_39c": min_dg_39,
                "min_dg_qpcr_55c": min_dg_55,
                "delta_tm": rec.get("delta_tm"),
                "dye_a": rec.get("dye_a"),
                "dye_b": rec.get("dye_b"),
                "dye_status": rec.get("dye_status")
            }
            row = {MULTIPLEX_COLUMNS[k]: v for k, v in raw_row.items() if k in MULTIPLEX_COLUMNS and v is not None}
            for pk, pv in rec.get("cross_dimers_39c", {}).items(): row[f"ΔG 39C: {pk}"] = pv
            rows.append(row)
        df = pd.DataFrame(rows)
        s_name = f"{'Triplex' if 'triplex' in jf.name else 'Pair'}_{m_dir.name.replace('multiplex-', '').replace('influenza-', '').replace('all-partial', 'All').replace('indo-only', 'Indo')}"
        ws = wb.create_sheet(title=safe_sheet_name(s_name))
        write_df_to_ws(ws, df, f"Multiplex — {m_dir.name}", "Status (RAA 39°C)")
        df.to_csv(str(csv_dir / f"S3_{s_name}.csv"), index=False)
        count += 1
    if count > 0:
        wb.remove(wb.active)
    else:
        wb.active.title = "No_Data"
        wb.active.append(["No multiplex data matching the active filter targets."])
    wb.save(str(outdir / "S3_Multiplex_Compatibility.xlsx"))
    print(f"  ✓ Created S3_Multiplex_Compatibility.xlsx ({count} sheets)")

def export_s4_qc(results_root: Path, outdir: Path):
    wb = openpyxl.Workbook()
    csv_dir = outdir / "csv"
    rows = []
    for b in sorted(results_root.iterdir()):
        if not b.is_dir() or b.name.startswith("multiplex"): continue
        qf = b / "m02" / "qc_report.json"
        if not qf.exists(): continue
        with open(qf, "r") as f: qc = json.load(f)
        row = {"Build": b.name, "Input Sequences": qc.get("input_count"), "QC Passed": qc.get("output_count"), "Retention Rate (%)": qc.get("retention_pct")}
        for k in ["subtype_filtering", "filters"]:
            v = qc.get(k)
            if isinstance(v, dict): row[k.title()] = "; ".join(f"{sub_k}: {sub_v}" for sub_k, sub_v in v.items())
            elif v is not None: row[k.title()] = str(v)
        rows.append(row)
    df = pd.DataFrame(rows) if rows else pd.DataFrame([{"Status": "No QC Data"}])
    ws = wb.create_sheet(title="QC_All_Builds")
    write_df_to_ws(ws, df, "QC & Sequence Filtering Summary — All Builds")
    df.to_csv(str(csv_dir / "S4_QC_All_Builds.csv"), index=False)
    wb.remove(wb.active)
    wb.save(str(outdir / "S4_QC_Summary.xlsx"))
    print(f"  ✓ Created S4_QC_Summary.xlsx")

def export_s5_conservation(results_root: Path, outdir: Path):
    wb = openpyxl.Workbook()
    csv_dir = outdir / "csv"
    count = 0
    for b in sorted(results_root.iterdir()):
        if not b.is_dir() or b.name.startswith("multiplex"): continue
        cf = b / "m07" / "candidate_targets_summary.csv"
        if not cf.exists(): continue
        df = pd.read_csv(cf)
        s_name = f"Conservation_{build_label(b.name)}"
        ws = wb.create_sheet(title=safe_sheet_name(s_name))
        write_df_to_ws(ws, df, f"Conservation & Scoring — {b.name}")
        df.to_csv(str(csv_dir / f"S5_{s_name}.csv"), index=False)
        count += 1
    if count > 0:
        wb.remove(wb.active)
    else:
        wb.active.title = "No_Data"
        wb.active.append(["No conservation scoring data matching the active filter targets."])
    wb.save(str(outdir / "S5_Conservation_Scoring.xlsx"))
    print(f"  ✓ Created S5_Conservation_Scoring.xlsx ({count} sheets)")

# ---------------------------------------------------------------------------
# Main Routine
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Full Executive & Detailed Supplementary Exporter")
    parser.add_argument("--results-dir", default="results", help="Results directory")
    parser.add_argument("--outdir", default="supplementary", help="Output directory")
    parser.add_argument("--only-degenerate", action="store_true", help="Only export degenerate candidates")
    parser.add_argument("configs", nargs="*", help="Config YAML files to detect active targets")
    args = parser.parse_args()

    results_root = Path(args.results_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "csv").mkdir(exist_ok=True)

    # Detect which target keywords are active from YAML config parameters
    active_targets = []
    for cfg in args.configs:
        cfg_name = Path(cfg).stem.lower()
        if "influenza_a" in cfg_name or "influenza-a" in cfg_name or "flu_a" in cfg_name or "flu-a" in cfg_name:
            active_targets.append("influenza-a")
        elif "influenza_b" in cfg_name or "influenza-b" in cfg_name or "flu_b" in cfg_name or "flu-b" in cfg_name:
            active_targets.append("influenza-b")
        elif "measles" in cfg_name:
            active_targets.append("measles")
            
    if not active_targets:
        active_targets = ["influenza-a", "influenza-b", "measles"]

    print(f"\n================================================================================")
    print(f"📊 EXPORTING COMPLETE PACKAGE (EXECUTIVE SUMMARY + S1 - S5 WORKBOOKS)")
    print(f"🎯 Active targets matching: {active_targets}")
    print(f"================================================================================")

    # Subclass Path to dynamically filter results_root.iterdir() and results_root.glob()
    class FilteredPath(type(Path())):
        def iterdir(self):
            raw = super().iterdir()
            filtered = []
            for b in raw:
                # Keep target builds matching active targets
                if not b.is_dir():
                    continue
                # Handle multiplex folder differently (it contains multi-target combination info)
                if b.name.startswith("multiplex"):
                    # Only include multiplex directories that contain at least one excluded target
                    # If measles is not active, skip multiplex folders containing measles
                    if "measles" not in active_targets and "measles" in b.name.lower():
                        continue
                    filtered.append(b)
                elif any(t in b.name.lower() for t in active_targets):
                    filtered.append(b)
            return filtered
            
        def glob(self, pattern):
            raw = super().glob(pattern)
            filtered = []
            for p in raw:
                # Filter path list matching active targets
                p_str = str(p).lower()
                if "measles" not in active_targets and "measles" in p_str:
                    continue
                filtered.append(p)
            return filtered

    results_root_filtered = FilteredPath(results_root)

    # 1. Main Executive Summary Workbook
    create_executive_summary(results_root_filtered, outdir)

    # 2. Detailed Supporting Workbooks
    export_s1_candidates(results_root_filtered, outdir, args.only_degenerate)
    export_s2_blast(results_root_filtered, outdir, args.only_degenerate)
    export_s3_multiplex(results_root_filtered, outdir)
    export_s4_qc(results_root_filtered, outdir)
    export_s5_conservation(results_root_filtered, outdir)

    print(f"================================================================================")
    print(f"✅ COMPLETE EXPORT SUCCESSFUL! All workbooks generated in: {outdir}")
    print(f"================================================================================")

if __name__ == "__main__":
    main()
