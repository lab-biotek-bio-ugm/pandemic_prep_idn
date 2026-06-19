#!/usr/bin/env python3
"""
M10 Run Comparison Tool
Compares non-degenerate vs degenerate primer design runs and generates a beautiful HTML report.
"""

import argparse
import json
import os
import sys
from pathlib import Path

def generate_html_report(non_deg_data: list, deg_data: list, output_path: Path, build_name: str, has_degenerate: bool):
    # Parse candidates
    non_deg_map = {c["cand_id"]: c for c in non_deg_data}
    deg_map = {c["cand_id"]: c for c in deg_data} if deg_data else {}

    # Compile table rows
    all_cand_ids = sorted(list(set(list(non_deg_map.keys()) + list(deg_map.keys()))))
    
    rows_html = []
    
    # We want to display side-by-side comparison for each candidate slot (by rank or ID)
    for idx, cand_id in enumerate(all_cand_ids):
        nd_cand = non_deg_map.get(cand_id)
        d_cand = deg_map.get(cand_id)
        
        # Check if we have both or just one
        nd_seqs = f"<b>F:</b> {nd_cand['fwd_seq']}<br><b>R:</b> {nd_cand['rev_seq']}"
        if nd_cand.get("probe_seq"):
            nd_seqs += f"<br><b>P:</b> {nd_cand['probe_seq']}"
            
        nd_cov = f"{nd_cand['raa_detected_pct']}% ({nd_cand['raa_detected_seqs']}/{nd_cand['total_seqs']})"
        nd_amp = f"{nd_cand['raa_amplified_pct']}%"
        
        if has_degenerate and d_cand:
            d_seqs = f"<b>F:</b> {d_cand['fwd_seq']}<br><b>R:</b> {d_cand['rev_seq']}"
            if d_cand.get("probe_seq"):
                d_seqs += f"<br><b>P:</b> {d_cand['probe_seq']}"
                
            d_cov = f"{d_cand['raa_detected_pct']}% ({d_cand['raa_detected_seqs']}/{d_cand['total_seqs']})"
            d_amp = f"{d_cand['raa_amplified_pct']}%"
            
            # Highlight coverage increase
            diff = d_cand['raa_detected_pct'] - nd_cand['raa_detected_pct']
            if diff > 0:
                badge = f'<span class="badge positive">+{diff:.2f}%</span>'
            elif diff < 0:
                badge = f'<span class="badge negative">{diff:.2f}%</span>'
            else:
                badge = '<span class="badge neutral">0.00%</span>'
        else:
            d_seqs = '<span class="disabled">Degenerate Mode Disabled</span>'
            d_cov = '-'
            d_amp = '-'
            badge = '-'

        rows_html.append(f"""
        <tr>
            <td class="cand-id">Rank {idx + 1}<br><small>{cand_id}</small></td>
            <td class="seq-cell">{nd_seqs}</td>
            <td class="center-cell highlight-nd">{nd_amp}</td>
            <td class="center-cell highlight-nd"><b>{nd_cov}</b></td>
            <td class="seq-cell">{d_seqs}</td>
            <td class="center-cell highlight-d">{d_amp}</td>
            <td class="center-cell highlight-d"><b>{d_cov}</b></td>
            <td class="center-cell">{badge}</td>
        </tr>
        """)

    # Subtype breakdown comparison for Rank 1
    subtype_rows_html = []
    if all_cand_ids:
        r1_id = all_cand_ids[0]
        nd_r1 = non_deg_map.get(r1_id)
        d_r1 = deg_map.get(r1_id) if has_degenerate else None
        
        nd_subtypes = {s["subtype"]: s for s in nd_r1["breakdowns"]["subtype"]} if nd_r1 else {}
        d_subtypes = {s["subtype"]: s for s in d_r1["breakdowns"]["subtype"]} if d_r1 else {}
        
        all_subtypes = sorted(list(set(list(nd_subtypes.keys()) + list(d_subtypes.keys()))))
        for st in all_subtypes:
            nd_st = nd_subtypes.get(st)
            d_st = d_subtypes.get(st)
            
            nd_pct = f"{nd_st['raa_detected_pct']}%" if nd_st else "0.0%"
            nd_abs = f"({nd_st['raa_detected']}/{nd_st['total']})" if nd_st else "-"
            
            if has_degenerate and d_st:
                d_pct = f"{d_st['raa_detected_pct']}%"
                d_abs = f"({d_st['raa_detected']}/{d_st['total']})"
                diff = d_st['raa_detected_pct'] - nd_st['raa_detected_pct'] if nd_st else d_st['raa_detected_pct']
                if diff > 0:
                    badge = f'<span class="badge positive">+{diff:.2f}%</span>'
                elif diff < 0:
                    badge = f'<span class="badge negative">{diff:.2f}%</span>'
                else:
                    badge = '<span class="badge neutral">0.00%</span>'
            else:
                d_pct = "-"
                d_abs = "-"
                badge = "-"

            subtype_rows_html.append(f"""
            <tr>
                <td><b>{st}</b></td>
                <td class="center-cell highlight-nd">{nd_pct} <small class="text-muted">{nd_abs}</small></td>
                <td class="center-cell highlight-d">{d_pct} <small class="text-muted">{d_abs}</small></td>
                <td class="center-cell">{badge}</td>
            </tr>
            """)

    # HTML Page Template
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Universal Influenza A Diagnostics — Comparison Report</title>
    <style>
        :root {{
            --bg-color: #0f172a;
            --panel-bg: #1e293b;
            --border-color: #334155;
            --text-color: #f8fafc;
            --text-muted: #94a3b8;
            --primary: #3b82f6;
            --success: #10b981;
            --warning: #f59e0b;
            --danger: #ef4444;
            --nd-color: #6366f1;
            --d-color: #10b981;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            background-color: var(--bg-color);
            color: var(--text-color);
            margin: 0;
            padding: 30px;
            line-height: 1.5;
        }}
        .container {{
            max-width: 1300px;
            margin: 0 auto;
        }}
        header {{
            margin-bottom: 40px;
            border-bottom: 2px solid var(--border-color);
            padding-bottom: 20px;
        }}
        h1 {{
            margin: 0 0 10px 0;
            font-size: 2.2rem;
            font-weight: 800;
            background: linear-gradient(to right, #60a5fa, #34d399);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }}
        .meta {{
            color: var(--text-muted);
            font-size: 0.95rem;
        }}
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .card {{
            background-color: var(--panel-bg);
            border: 1px solid var(--border-color);
            border-radius: 12px;
            padding: 20px;
        }}
        .card-title {{
            font-size: 1.1rem;
            font-weight: 700;
            margin-top: 0;
            margin-bottom: 15px;
            color: var(--text-muted);
            text-transform: uppercase;
            letter-spacing: 0.05em;
        }}
        .metric-val {{
            font-size: 2.5rem;
            font-weight: 800;
            line-height: 1;
            margin-bottom: 10px;
        }}
        .metric-sub {{
            font-size: 0.85rem;
            color: var(--text-muted);
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 30px;
            background-color: var(--panel-bg);
            border: 1px solid var(--border-color);
            border-radius: 12px;
            overflow: hidden;
        }}
        th, td {{
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid var(--border-color);
        }}
        th {{
            background-color: #151f32;
            font-weight: 700;
            color: var(--text-muted);
            font-size: 0.85rem;
            text-transform: uppercase;
        }}
        tr:last-child td {{
            border-bottom: none;
        }}
        .seq-cell {{
            font-family: "Courier New", Courier, monospace;
            font-size: 0.9rem;
            white-space: nowrap;
        }}
        .center-cell {{
            text-align: center;
        }}
        .cand-id {{
            font-weight: 700;
        }}
        .disabled {{
            color: var(--text-muted);
            font-style: italic;
        }}
        .badge {{
            display: inline-block;
            padding: 4px 8px;
            border-radius: 6px;
            font-size: 0.8rem;
            font-weight: 700;
        }}
        .badge.positive {{
            background-color: rgba(16, 185, 129, 0.15);
            color: var(--success);
        }}
        .badge.negative {{
            background-color: rgba(239, 68, 68, 0.15);
            color: var(--danger);
        }}
        .badge.neutral {{
            background-color: rgba(148, 163, 184, 0.15);
            color: var(--text-muted);
        }}
        .highlight-nd {{
            background-color: rgba(99, 102, 241, 0.05);
            border-left: 1px solid rgba(99, 102, 241, 0.2);
            border-right: 1px solid rgba(99, 102, 241, 0.2);
        }}
        .highlight-d {{
            background-color: rgba(16, 185, 129, 0.05);
            border-left: 1px solid rgba(16, 185, 129, 0.2);
            border-right: 1px solid rgba(16, 185, 129, 0.2);
        }}
        .footer {{
            margin-top: 50px;
            text-align: center;
            color: var(--text-muted);
            font-size: 0.85rem;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>Universal Influenza A Diagnostics Comparison Report</h1>
            <div class="meta">
                Build: <b>{build_name}</b> | Status: <b>Dual-Track Completed</b> | Target segment: <b>Matrix Protein (MP)</b>
            </div>
        </header>

        <div class="grid">
            <div class="card" style="border-top: 4px solid var(--nd-color);">
                <div class="card-title">Non-Degenerate Best Coverage</div>
                <div class="metric-val" style="color: #818cf8;">{non_deg_data[0]['raa_detected_pct'] if non_deg_data else 0.0}%</div>
                <div class="metric-sub">Total Covered: {non_deg_data[0]['raa_detected_seqs'] if non_deg_data else 0} / {non_deg_data[0]['total_seqs'] if non_deg_data else 0} sequences</div>
            </div>
            
            <div class="card" style="border-top: 4px solid var(--d-color);">
                <div class="card-title">Degenerate Best Coverage</div>
                <div class="metric-val" style="color: #34d399;">
                    {deg_data[0]['raa_detected_pct'] if (has_degenerate and deg_data) else "N/A"}%
                </div>
                <div class="metric-sub">
                    {f"Total Covered: {deg_data[0]['raa_detected_seqs']} / {deg_data[0]['total_seqs']} sequences" if (has_degenerate and deg_data) else "Degenerate simulation disabled"}
                </div>
            </div>

            <div class="card" style="border-top: 4px solid var(--warning);">
                <div class="card-title">Coverage Uplift</div>
                <div class="metric-val" style="color: #fbbf24;">
                    {f"+{round(deg_data[0]['raa_detected_pct'] - non_deg_data[0]['raa_detected_pct'], 2)}%" if (has_degenerate and deg_data and non_deg_data) else "0.00%"}
                </div>
                <div class="metric-sub">Net increase in diagnostic universality</div>
            </div>
        </div>

        <h2>Candidate Comparison Table</h2>
        <table>
            <thead>
                <tr>
                    <th rowspan="2">Candidate</th>
                    <th colspan="3" class="center-cell" style="background-color: rgba(99, 102, 241, 0.2);">Non-Degenerate Track</th>
                    <th colspan="3" class="center-cell" style="background-color: rgba(16, 185, 129, 0.2);">Degenerate Track</th>
                    <th rowspan="2" class="center-cell">Uplift</th>
                </tr>
                <tr>
                    <th style="background-color: rgba(99, 102, 241, 0.1);">Sequences</th>
                    <th class="center-cell" style="background-color: rgba(99, 102, 241, 0.1);">Amplified %</th>
                    <th class="center-cell" style="background-color: rgba(99, 102, 241, 0.1);">Detected %</th>
                    <th style="background-color: rgba(16, 185, 129, 0.1);">Sequences</th>
                    <th class="center-cell" style="background-color: rgba(16, 185, 129, 0.1);">Amplified %</th>
                    <th class="center-cell" style="background-color: rgba(16, 185, 129, 0.1);">Detected %</th>
                </tr>
            </thead>
            <tbody>
                {"".join(rows_html)}
            </tbody>
        </table>

        <h2>Subtype Coverage Breakdown (Rank 1 Candidate)</h2>
        <table>
            <thead>
                <tr>
                    <th>Subtype</th>
                    <th class="center-cell highlight-nd">Non-Degenerate Coverage</th>
                    <th class="center-cell highlight-d">Degenerate Coverage</th>
                    <th class="center-cell">Uplift</th>
                </tr>
            </thead>
            <tbody>
                {"".join(subtype_rows_html)}
            </tbody>
        </table>

        <div class="footer">
            Pandemic Preparedness Pipeline &bull; Automated Diagnostic Design Verification
        </div>
    </div>
</body>
</html>
"""

    with open(output_path, "w") as f:
        f.write(html_content)

def main():
    parser = argparse.ArgumentParser(description="M10 Run Comparison Tool")
    parser.add_argument("--non-degenerate-summary", required=True, help="Path to non-degenerate coverage_summary.json")
    parser.add_argument("--degenerate-summary", default=None, help="Path to degenerate coverage_summary.json (optional)")
    parser.add_argument("--output-report", required=True, help="Path to write the comparison HTML report")
    parser.add_argument("--output-marker", required=True, help="Path to write verification complete marker file")
    parser.add_argument("--build-name", required=True, help="Run build name metadata")

    args = parser.parse_args()

    # Read Non-Degenerate Summary
    print(f"[Compare] Reading non-degenerate summary from {args.non_degenerate_summary}...")
    with open(args.non_degenerate_summary, 'r') as f:
        non_deg_data = json.load(f)

    # Read Degenerate Summary if available
    deg_data = None
    has_degenerate = False
    if args.degenerate_summary and os.path.exists(args.degenerate_summary):
        print(f"[Compare] Reading degenerate summary from {args.degenerate_summary}...")
        with open(args.degenerate_summary, 'r') as f:
            deg_data = json.load(f)
            has_degenerate = True

    # Generate Laporan HTML
    output_path = Path(args.output_report)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"[Compare] Generating comparison HTML report at {output_path}...")
    generate_html_report(non_deg_data, deg_data, output_path, args.build_name, has_degenerate)

    # Tulis Marker File
    print(f"[Compare] Writing complete marker to {args.output_marker}...")
    with open(args.output_marker, "w") as f:
        f.write(f"Comparison completed successfully.\n")
        f.write(f"Non-degenerate summary: {args.non_degenerate_summary}\n")
        if has_degenerate:
            f.write(f"Degenerate summary: {args.degenerate_summary}\n")
        else:
            f.write("Degenerate summary: disabled\n")

    print("[Compare] Run comparison completed successfully.")

if __name__ == "__main__":
    main()
