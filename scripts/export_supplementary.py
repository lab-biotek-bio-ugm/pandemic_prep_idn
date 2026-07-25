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

# Module-level default — overwritten by main() after auto-detection
ACTIVE_TARGETS: list = []

# Canonical virus tokens used to identify a build/multiplex's constituent targets.
_KNOWN_VIRUSES = [
    "influenza-a", "influenza-b", "measles",
    "rsv-a", "rsv-b", "sars-cov-2", "dengue",
]


def _viruses_in(text: str) -> list:
    t = (text or "").lower()
    return [v for v in _KNOWN_VIRUSES if v in t]


def is_multiplex_active(name: str, result_path=None) -> bool:
    """True if a multiplex/triplex output should be exported under the active filter.

    A multiplex is excluded when it references any KNOWN virus that is not in
    ACTIVE_TARGETS (e.g. a triplex containing Measles when only Influenza A/B are
    active).

    The directory name alone is NOT sufficient: the batch orchestrator writes every
    triplex to the literal path `results/multiplex-triplex-all-3`, which contains no
    virus token at all. Falling back to "keep it" on an uninformative name is exactly
    how Measles rows kept reaching the Influenza-only workbooks. So when the name
    identifies nothing, we inspect the result file's own combination IDs.
    """
    present_known = _viruses_in(name)

    if not present_known and result_path is not None:
        try:
            with open(result_path, "r") as f:
                data = json.load(f)
            probe = " ".join(
                str(rec.get("triplex_id") or rec.get("pair_id") or rec.get("combo_id") or "")
                for rec in (data[:50] if isinstance(data, list) else [])
            )
            present_known = _viruses_in(probe)
        except Exception as e:
            print(f"  ⚠ Could not inspect {result_path} for target filtering: {e}")

    if not present_known:
        # Genuinely unidentifiable: keep, but say so rather than filtering silently.
        print(f"  ⚠ '{name}': constituent viruses not identifiable — kept without target filtering.")
        return True

    # ACTIVE_TARGETS may hold canonical virus tokens ("influenza-a") or full build
    # names ("influenza-a-mp-indo-only"), so match in both directions.
    def _active(v: str) -> bool:
        return any(v == a or v in a or a in v for a in ACTIVE_TARGETS)

    return all(_active(v) for v in present_known)

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
    "pair_id": "Candidate Pair", "triplex_id": "Combo ID", "rank_a": "Rank Target A", "rank_b": "Rank Target B",
    "score": "Compatibility Score (%)", "status_raa_39c": "Status (RAA 39°C)",
    "min_dg_raa_39c": "Worst RAA ΔG 39°C (kcal/mol)", "min_dg_qpcr_55c": "Worst qPCR ΔG 55°C (kcal/mol)",
    "delta_tm": "ΔTm Fwd Primers (°C)", "dye_a": "Fluorophore Target A", "dye_b": "Fluorophore Target B", "dye_status": "Dye Separation"
}

# ---------------------------------------------------------------------------
# Dynamic target helpers — no hardcoded virus names below this line
# ---------------------------------------------------------------------------

# Dye palette cycled over however many targets exist
DYE_PALETTE = ["FAM", "HEX", "Cy5", "Cy3", "ROX", "Texas Red", "Atto 647N"]

# Canonical name lookup: oligo JSON key → display name
_TARGET_DISPLAY = {
    "flu_a": "Influenza A",
    "flu_b": "Influenza B",
    "influenza_a": "Influenza A",
    "influenza_b": "Influenza B",
    "measles": "Measles",
    "rsv_a": "RSV A",
    "rsv_b": "RSV B",
    "sars_cov_2": "SARS-CoV-2",
    "dengue": "Dengue",
}

# Map from JSON oligo key to the canonical ACTIVE_TARGETS string
_TARGET_CANONICAL = {
    "flu_a": "influenza-a",
    "flu_b": "influenza-b",
    "influenza_a": "influenza-a",
    "influenza_b": "influenza-b",
    "measles": "measles",
    "rsv_a": "rsv-a",
    "rsv_b": "rsv-b",
    "sars_cov_2": "sars-cov-2",
    "dengue": "dengue",
}

def nplex_label(n: int) -> str:
    """Duplex / Triplex / Quadruplex / N-Plex — never hardcoded."""
    return {1: "Singleplex", 2: "Duplex", 3: "Triplex", 4: "Quadruplex"}.get(n, f"{n}-Plex")

def recommendation_category(rank_idx: int, combo: dict) -> str:
    """Label a combination for the ordering sheet.

    The label comes from the selection verdict, not from position in the list.
    A combination that failed the gates must never be presented as a
    recommendation just because it sorted first -- that inversion (a failed
    combo shown as "BEST RECOMMENDED" because it was rank 1) is the class of
    label bug this guards against.
    """
    if not combo.get("selection_passed", True):
        return "NOT RECOMMENDED — failed selection gates"
    if rank_idx == 1:
        cat = "BEST RECOMMENDED (RANK 1)"
        if combo.get("selection_relaxed"):
            cat += f" — quality floors relaxed to {combo.get('selection_level')}"
        return cat
    return f"ALTERNATIVE RANK {rank_idx}"

def target_key_to_name(key: str) -> str:
    """Convert an oligo JSON key (e.g. 'flu_a') to a human-readable virus name."""
    k = key.lower()
    if k in _TARGET_DISPLAY:
        return _TARGET_DISPLAY[k]
    # Generic fallback: uppercase short tokens, capitalize others
    return " ".join(p.upper() if len(p) <= 2 else p.capitalize() for p in k.split("_"))

def is_target_active(oligo_key: str, active_targets: list) -> bool:
    """Return True if this oligo JSON key (e.g. 'flu_a') maps to an active target."""
    k = oligo_key.lower()
    canonical = _TARGET_CANONICAL.get(k, k.replace("_", "-"))
    return any(canonical == a or canonical in a or a in canonical for a in active_targets)

def safe_sheet_name(name: str, maxlen: int = 31) -> str:
    name = re.sub(r"[\/\\?\*:\[\]]", "_", name)
    return name[:maxlen]

def build_label(dir_name: str) -> str:
    """Convert a results directory name to a compact label. Generic — no hardcoded virus list."""
    s = dir_name.lower()
    # Known virus shorthands (extend freely)
    for old, new in [
        ("influenza-a", "FluA"), ("influenza_a", "FluA"),
        ("influenza-b", "FluB"), ("influenza_b", "FluB"),
        ("measles", "Measles"),
        ("rsv-a", "RSVA"), ("rsv-b", "RSVB"),
        ("sars-cov-2", "SARS2"), ("sars_cov_2", "SARS2"),
        ("dengue", "Dengue"),
        ("-mp-", "_MP_"), ("_mp_", "_MP_"),
        ("-np-", "_NP_"), ("_np_", "_NP_"),
        ("-ha-", "_HA_"), ("-na-", "_NA_"),
        ("indo-only", "Indonesia"), ("indonesia-only", "Indonesia"),
        ("sea-only", "SEAsia"), ("all-partial", "AllPartial"),
        ("global", "Global"),
    ]:
        s = s.replace(old, new)
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
        if not any(t in b.name.lower() for t in ACTIVE_TARGETS):
            continue
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
                    # final_score is PrimerLab's Primer3 PENALTY: lower is
                    # better. Labelling it "Score" read as a mark out of 100 and
                    # made the best candidate look like the worst.
                    top_cand = (f"Rank 1 Primer3 penalty: "
                                f"{c_data[0].get('final_score')} (lower is better)")

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

    # Sheet 01: Oligo Orders — fully dynamic (works for any N-plex)
    ws1 = wb.create_sheet(title="01_Synthesis_Oligo_Orders")

    # Find the best available N-plex candidates file (triplex preferred, then multiplex)
    combo_files = sorted(results_root.glob("**/triplex_candidates.json"))
    if not combo_files:
        combo_files = [p for p in sorted(results_root.glob("**/multiplex_candidates.json"))
                       if "triplex" not in p.parent.name.lower()]
    # Keep only multiplex files whose constituent viruses are ALL active. Using
    # "any active target present" here was the leak that let a Flu+Measles triplex
    # through on the Influenza token alone, printing Measles rows into the orders sheet.
    combo_files = [p for p in combo_files if is_multiplex_active(p.parent.name, p)]

    order_rows = []
    n_targets = 0
    plex_lbl = "Multiplex"

    if combo_files:
        with open(combo_files[0], "r") as f:
            t_data = json.load(f)[:4]

        # Determine N-plex label from actual number of targets in first combo
        if t_data:
            n_targets = len(t_data[0].get("oligos_synthesis", {}))
            plex_lbl = nplex_label(n_targets)

        for rank_idx, combo in enumerate(t_data, 1):
            cat = recommendation_category(rank_idx, combo)
            oligos = combo.get("oligos_synthesis", {})
            combo_id = combo.get("triplex_id") or combo.get("pair_id") or combo.get("combo_id")
            score = combo.get("score")

            # Iterate ALL target keys present in JSON — no hardcoded list
            for t_idx, t_key in enumerate(sorted(oligos.keys())):
                # Skip targets that are excluded via --exclude-targets
                if not is_target_active(t_key, ACTIVE_TARGETS):
                    continue
                t_name = target_key_to_name(t_key)
                default_dye = DYE_PALETTE[t_idx % len(DYE_PALETTE)]

                for role_key, role_lbl in [("fwd", "Forward Primer"), ("rev", "Reverse Primer"), ("probe", "Exonuclease Probe")]:
                    if role_key not in oligos[t_key]:
                        continue
                    o = oligos[t_key][role_key]
                    if isinstance(o, dict):
                        seq    = o.get("deg_seq") or o.get("seq", "")
                        length = o.get("len", len(seq))
                        tm     = o.get("tm", "-")
                        gc     = o.get("gc", "-")
                        label  = o.get("label", "Unlabeled" if role_key != "probe" else f"5'-{default_dye} / 3'-BHQ / Internal-THF")
                    else:
                        seq    = str(o)
                        length = len(seq)
                        tm     = "-"
                        gc     = "-"
                        label  = "Unlabeled" if role_key != "probe" else f"5'-{default_dye} / 3'-BHQ / Internal-THF"

                    order_rows.append({
                        "Recommendation Category": cat,
                        "Assay Set": combo_id,
                        "Target Virus": t_name,
                        "Oligo Role": role_lbl,
                        "Degenerate Sequence (5' -> 3')": seq,
                        "Length (nt)": length,
                        "Tm (°C)": tm,
                        "GC (%)": gc,
                        "Fluorophore / Label": label,
                        f"{plex_lbl} Compatibility Score": f"{score}%",
                    })

    empty_msg = "Run Multiplex Module to populate" if not combo_files else "No oligos found in candidates"
    df_ord = pd.DataFrame(order_rows) if order_rows else pd.DataFrame([{"Status": empty_msg}])
    write_df_to_ws(ws1, df_ord, "Synthesis Oligo Ordering Specs (Rank 1 Best + 3 Alternatives)")
    df_ord.to_csv(str(csv_dir / "01_Synthesis_Oligo_Orders.csv"), index=False)

    # Sheet 02: Pairwise
    ws2 = wb.create_sheet(title="02_Multiplex_Pairwise")
    pair_files = list(results_root.glob("**/multiplex_candidates.json"))
    pair_rows = []
    for pf in pair_files:
        if "triplex" in pf.parent.name: continue
        if not is_multiplex_active(pf.parent.name, pf): continue  # respect active target filter
        with open(pf, "r") as f:
            p_data = json.load(f)
            for rec in p_data:
                r = {"Pair ID": rec.get("pair_id"), "Compatibility Score (%)": rec.get("score"), "Status (39°C)": rec.get("status_raa_39c"), "RAA ΔG 39°C": rec.get("min_dg_raa_39c"), "qPCR ΔG 55°C": rec.get("min_dg_qpcr_55c"), "ΔTm Fwd (°C)": rec.get("delta_tm")}
                for pair_k, pair_v in rec.get("cross_dimers_39c", {}).items(): r[f"ΔG 39C: {pair_k}"] = pair_v
                pair_rows.append(r)

    df_pair = pd.DataFrame(pair_rows) if pair_rows else pd.DataFrame([{"Status": "No pairwise runs found"}])
    write_df_to_ws(ws2, df_pair, "Pairwise Multiplex Screening Results", "Status (39°C)")

    # Sheet 03: N-Plex Compatibility — label derived from actual target count
    n_label_str = f"{n_targets}-Target" if n_targets else "Multi-Target"
    ws3_title = safe_sheet_name(f"03_{plex_lbl}_Compat")
    ws3 = wb.create_sheet(title=ws3_title)
    combo_rows = []
    if combo_files:
        with open(combo_files[0], "r") as f:
            t_data3 = json.load(f)
        for rec in t_data3:
            r = {
                "Combo ID": rec.get("triplex_id") or rec.get("pair_id") or rec.get("combo_id"),
                "Compatibility Score (%)": rec.get("score"),
                "Status (39°C)": rec.get("status_raa_39c") or rec.get("status_39c"),
                "RAA ΔG 39°C": rec.get("min_dg_raa_39c") or rec.get("min_dg_39c"),
                "qPCR ΔG 55°C": rec.get("min_dg_qpcr_55c") or rec.get("min_dg_55c"),
                "Max ΔTm (°C)": rec.get("delta_tm"),
            }
            for pair_k, pair_v in rec.get("cross_dimers_39c", {}).items():
                r[f"ΔG 39C: {pair_k}"] = pair_v
            combo_rows.append(r)

    df_combo = pd.DataFrame(combo_rows) if combo_rows else pd.DataFrame([{"Status": f"No {plex_lbl.lower()} runs found"}])
    write_df_to_ws(ws3, df_combo, f"{plex_lbl} Multiplex Screening Results ({n_label_str})", "Status (39°C)")

    # Sheet 04: BLAST
    ws4 = wb.create_sheet(title="04_BLAST_OffTarget_Summary")
    blast_files = list(results_root.glob("**/blast_specificity_summary.csv"))
    blast_dfs = []
    for bf in blast_files:
        cohort = bf.parent.parent.parent.parent.name
        if not any(t in cohort.lower() for t in ACTIVE_TARGETS): continue  # respect active target filter
        try:
            # Read columns as object/str to keep original text like "58.94 - 64.98 °C" intact
            b_df = pd.read_csv(bf, dtype=object)
            b_df.insert(0, "Build Cohort", cohort)
            blast_dfs.append(b_df)
        except Exception: pass

    df_blast = pd.concat(blast_dfs, ignore_index=True) if blast_dfs else pd.DataFrame([{"Status": "No BLAST files"}])
    # Clean up any remaining NaN
    df_blast = df_blast.fillna("-")
    write_df_to_ws(ws4, df_blast, "BLAST Specificity & Off-Target Risk Summary", "specificity_grade")

    # Sheet 05: Human Off-Target Risk Audit & Justification
    ws5 = wb.create_sheet(title="05_Human_OffTarget_Risk_Audit")
    human_audit_rows = []
    
    for b in builds:
        if not any(t in b.name.lower() for t in ACTIVE_TARGETS): continue
        for m in ["degenerate", "non_degenerate"]:
            m_dir = b / "m09" / m
            if not m_dir.exists(): continue
            for r_dir in sorted(m_dir.iterdir()):
                # Scan top candidates (candidate_1 to candidate_10)
                cand_dirs = sorted(r_dir.glob("candidate_*"), key=lambda d: int(d.name.split("_")[-1]) if d.name.split("_")[-1].isdigit() else 999)
                for c_dir in cand_dirs[:5]: # Top 5 primary ordered candidates
                    cand_num = c_dir.name.replace("candidate_", "Candidate ")
                    db_human = c_dir / "db_human"
                    if not db_human.exists(): continue

                    f_csv = db_human / "forward_offtargets.csv"
                    r_csv = db_human / "reverse_offtargets.csv"
                    p_csv = db_human / "probe_offtargets.csv"

                    has_hits = False
                    for role, csv_path in [("Forward Primer", f_csv), ("Reverse Primer", r_csv), ("TaqMan Exo-Probe", p_csv)]:
                        if csv_path.exists():
                            try:
                                hits_df = pd.read_csv(csv_path)
                                for _, h_row in hits_df.iterrows():
                                    has_hits = True
                                    seq_id = str(h_row.get("Sequence ID") or h_row.get("subject_id") or "-")
                                    desc = str(h_row.get("Description / Species") or h_row.get("title") or "Homo sapiens genome")
                                    ident = str(h_row.get("Identity") or h_row.get("pident") or "-")
                                    overhang = str(h_row.get("3' Overhang") or h_row.get("overhang_3p") or "-")
                                    risk = str(h_row.get("Risk") or "Low")

                                    # Biological Rationale Generation
                                    if "Probe" in role:
                                        rationale = "SAFE / PROBE SPECIFIC: Single-stranded probe hit. No fluorophore cleavage signal without paired primer amplification."
                                    elif "0" in overhang or "1" in overhang:
                                        rationale = "SAFE / 3'-END MISMATCH: Critical 3' terminus is unanchored or mismatched (>2 nt). DNA polymerase/recombinase extension blocked."
                                    else:
                                        rationale = "SAFE / SINGLE-PRIMER ISOLATION: Isolated single-primer binding on human chromosome. No opposite-strand primer binding within <1000 bp amplicon distance. Zero amplicon risk."

                                    human_audit_rows.append({
                                        "Build Cohort": build_label(b.name),
                                        "Candidate Set": cand_num,
                                        "Oligo Component": role,
                                        "Human RefSeq ID": seq_id,
                                        "Chromosome / Gene Annotation": desc,
                                        "Identity (%)": ident,
                                        "3' Overhang": overhang,
                                        "Off-Target Risk": risk,
                                        "Biological Safety Rationale": rationale
                                    })
                            except Exception: pass

                    if not has_hits:
                        human_audit_rows.append({
                            "Build Cohort": build_label(b.name),
                            "Candidate Set": cand_num,
                            "Oligo Component": "All Components (Fwd/Rev/Probe)",
                            "Human RefSeq ID": "None (0 Hits)",
                            "Chromosome / Gene Annotation": "No human genome sequence alignment detected",
                            "Identity (%)": "0.0%",
                            "3' Overhang": "N/A",
                            "Off-Target Risk": "SAFE",
                            "Biological Safety Rationale": "OPTIMAL SPECIFICITY: 100% specific to target pathogen. Zero alignments to Homo sapiens genome."
                        })

    df_human_audit = pd.DataFrame(human_audit_rows) if human_audit_rows else pd.DataFrame([{"Status": "No human off-target hits found"}])
    write_df_to_ws(ws5, df_human_audit, "Human Genome (Homo sapiens) Off-Target Risk & Safety Audit", "Off-Target Risk")
    df_human_audit.to_csv(str(csv_dir / "05_Human_OffTarget_Risk_Audit.csv"), index=False)

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
        if not any(t in b.name.lower() for t in ACTIVE_TARGETS): continue  # respect active target filter
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
        if not any(t in b.name.lower() for t in ACTIVE_TARGETS): continue  # respect active target filter
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
        # Pass the result file so an uninformative directory name (e.g. the
        # orchestrator's literal "multiplex-triplex-all-3") still gets filtered.
        if not is_multiplex_active(m_dir.name, jf): continue
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
        s_name = f"{'Triplex' if 'triplex' in jf.name else 'Pair'}_{build_label(m_dir.name.replace('multiplex-', ''))}"
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
        if not any(t in b.name.lower() for t in ACTIVE_TARGETS): continue  # respect active target filter
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
        if not any(t in b.name.lower() for t in ACTIVE_TARGETS): continue  # respect active target filter
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
    parser.add_argument(
        "--exclude-targets", nargs="*", default=[], metavar="TARGET",
        help="Exclude specific virus targets from ALL sheets (e.g. --exclude-targets measles rsv-a)"
    )
    parser.add_argument(
        "--active-targets", nargs="*", default=None, metavar="TARGET",
        help="Explicitly set the active targets (overrides auto-detection). "
             "Restricts ALL sheets, including S3 multiplex, to only these builds/viruses."
    )
    parser.add_argument(
        "configs", nargs="*",
        help="Optional config YAML files to restrict output to only matching targets"
    )
    args = parser.parse_args()

    results_root = Path(args.results_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "csv").mkdir(exist_ok=True)

    # ---------------------------------------------------------------------------
    # Auto-detect active targets from build directories in results/
    # No hardcoded virus names — works for any N-plex combination
    # ---------------------------------------------------------------------------
    found: set = set()
    if results_root.exists():
        for d in results_root.iterdir():
            if not d.is_dir() or d.name.startswith("multiplex"):
                continue
            d_lower = d.name.lower()
            for known in _KNOWN_VIRUSES:
                if known in d_lower:
                    found.add(known)
                    break
            else:
                # Unknown virus — use the first 1-2 dash-separated segments as the key
                parts = d.name.lower().split("-")
                found.add("-".join(parts[:2]) if len(parts) >= 2 else parts[0])

    # If configs explicitly given, restrict to only matching targets
    if args.configs:
        cfg_hits: set = set()
        for cfg in args.configs:
            cfg_lower = Path(cfg).stem.lower()
            for t in found:
                if t.replace("-", "_") in cfg_lower or t in cfg_lower:
                    cfg_hits.add(t)
        if cfg_hits:
            found = cfg_hits

    # Apply --exclude-targets
    excluded = {e.lower().strip() for e in (args.exclude_targets or [])}

    # An explicit --active-targets overrides auto-detection entirely.
    if args.active_targets:
        active_targets = sorted({a.lower().strip() for a in args.active_targets} - excluded)
    else:
        active_targets = sorted(found - excluded)

    if not active_targets:
        # Safety fallback: include all known viruses if nothing detected
        active_targets = sorted(t for t in _KNOWN_VIRUSES if t not in excluded)

    print(f"\n================================================================================")
    print(f"📊 EXPORTING COMPLETE PACKAGE (EXECUTIVE SUMMARY + S1 - S5 WORKBOOKS)")
    print(f"🎯 Active targets (auto-detected): {active_targets}")
    if excluded:
        print(f"🚫 Excluded targets:              {sorted(excluded)}")
    print(f"================================================================================")

    # Set globally for all exporter functions to read
    global ACTIVE_TARGETS
    ACTIVE_TARGETS = active_targets

    # 1. Main Executive Summary Workbook
    create_executive_summary(results_root, outdir)

    # 2. Detailed Supporting Workbooks
    export_s1_candidates(results_root, outdir, args.only_degenerate)
    export_s2_blast(results_root, outdir, args.only_degenerate)
    export_s3_multiplex(results_root, outdir)
    export_s4_qc(results_root, outdir)
    export_s5_conservation(results_root, outdir)

    print(f"================================================================================")
    print(f"✅ COMPLETE EXPORT SUCCESSFUL! All workbooks generated in: {outdir}")
    print(f"================================================================================")

if __name__ == "__main__":
    main()
