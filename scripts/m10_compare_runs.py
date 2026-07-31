#!/usr/bin/env python3
"""
M10 Run Comparison Tool
Compares non-degenerate vs degenerate primer design runs and generates a beautiful HTML report
with dynamic toggle between Global Surveillance DB and Design DB.
"""

import argparse
import json
import os
import sys
from pathlib import Path

def generate_html_report(non_deg_data: list, deg_data: list, output_path: Path, build_name: str, has_degenerate: bool):
    # Serialise the lists to JSON for embedding in the HTML
    non_deg_json = json.dumps(non_deg_data)
    deg_json = json.dumps(deg_data) if has_degenerate else "null"
    has_deg_js = "true" if has_degenerate else "false"

    # Compile the subtype breakdown table statically for Rank 1 (Global mode only, but initialized statically)
    subtype_rows_html = []
    nd_r1 = non_deg_data[0] if (non_deg_data and len(non_deg_data) > 0) else None
    d_r1 = deg_data[0] if (has_degenerate and deg_data and len(deg_data) > 0) else None
    
    nd_subtypes = {s["subtype"]: s for s in nd_r1["breakdowns"]["subtype"]} if (nd_r1 and "breakdowns" in nd_r1 and "subtype" in nd_r1["breakdowns"]) else {}
    d_subtypes = {s["subtype"]: s for s in d_r1["breakdowns"]["subtype"]} if (d_r1 and "breakdowns" in d_r1 and "subtype" in d_r1["breakdowns"]) else {}
    
    all_subtypes = sorted(list(set(list(nd_subtypes.keys()) + list(d_subtypes.keys()))))
    for st in all_subtypes:
        nd_st = nd_subtypes.get(st)
        d_st = d_subtypes.get(st)
        
        nd_pct = f"{nd_st['raa_detected_pct']}%" if nd_st else "0.0%"
        nd_abs = f"({nd_st['raa_detected']}/{nd_st['total']})" if nd_st else "-"
        
        if has_degenerate and d_st:
            d_pct = f"{d_st['raa_detected_pct']}%"
            d_abs = f"({d_st['raa_detected']}/{d_st['total']})"
            
            # Hitung selisih
            nd_pct_val = nd_st['raa_detected_pct'] if nd_st else 0.0
            diff = d_st['raa_detected_pct'] - nd_pct_val
            
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
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Universal Influenza A Diagnostics — Comparison Report</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;600;700;800&family=JetBrains+Mono:wght@400;600&display=swap" rel="stylesheet">
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
            font-family: 'Outfit', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            background-color: var(--bg-color);
            color: var(--text-color);
            margin: 0;
            padding: 30px;
            line-height: 1.5;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
        }}
        header {{
            margin-bottom: 30px;
            border-bottom: 2px solid var(--border-color);
            padding-bottom: 20px;
        }}
        h1 {{
            margin: 0 0 10px 0;
            font-size: 2.4rem;
            font-weight: 800;
            background: linear-gradient(to right, #60a5fa, #34d399);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }}
        .meta {{
            color: var(--text-muted);
            font-size: 0.95rem;
        }}
        
        /* Toggle Switch Styling */
        .toggle-container {{
            display: inline-flex;
            background-color: #1e293b;
            border: 1px solid var(--border-color);
            border-radius: 12px;
            padding: 4px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
        }}
        .toggle-btn {{
            background: none;
            border: none;
            color: var(--text-muted);
            padding: 10px 20px;
            border-radius: 8px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.25s cubic-bezier(0.4, 0, 0.2, 1);
            font-family: inherit;
            font-size: 0.95rem;
            display: flex;
            align-items: center;
            gap: 8px;
        }}
        .toggle-btn.active {{
            background-color: var(--primary);
            color: #ffffff;
            box-shadow: 0 4px 12px -2px rgba(59, 130, 246, 0.3);
        }}
        .toggle-btn:hover:not(.active) {{
            color: var(--text-color);
            background-color: rgba(255, 255, 255, 0.05);
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
            border-radius: 16px;
            padding: 24px;
            transition: transform 0.3s ease, border-color 0.3s ease;
            box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1);
        }}
        .card:hover {{
            transform: translateY(-2px);
            border-color: rgba(255,255,255,0.15);
        }}
        .card-title {{
            font-size: 1rem;
            font-weight: 700;
            margin-top: 0;
            margin-bottom: 12px;
            color: var(--text-muted);
            text-transform: uppercase;
            letter-spacing: 0.08em;
        }}
        .metric-val {{
            font-size: 2.8rem;
            font-weight: 800;
            line-height: 1;
            margin-bottom: 8px;
            transition: color 0.3s ease;
        }}
        .metric-sub {{
            font-size: 0.85rem;
            color: var(--text-muted);
        }}
        h2 {{
            font-size: 1.5rem;
            margin-top: 40px;
            margin-bottom: 20px;
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 30px;
            background-color: var(--panel-bg);
            border: 1px solid var(--border-color);
            border-radius: 16px;
            overflow: hidden;
            box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1);
        }}
        th, td {{
            padding: 14px 18px;
            text-align: left;
            border-bottom: 1px solid var(--border-color);
        }}
        th {{
            background-color: #151f32;
            font-weight: 700;
            color: var(--text-muted);
            font-size: 0.85rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
        }}
        tr:last-child td {{
            border-bottom: none;
        }}
        .seq-cell {{
            font-family: 'JetBrains Mono', monospace;
            font-size: 0.85rem;
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
            background-color: rgba(99, 102, 241, 0.03);
            border-left: 1px solid rgba(99, 102, 241, 0.15);
            border-right: 1px solid rgba(99, 102, 241, 0.15);
        }}
        .highlight-d {{
            background-color: rgba(16, 185, 129, 0.03);
            border-left: 1px solid rgba(16, 185, 129, 0.15);
            border-right: 1px solid rgba(16, 185, 129, 0.15);
        }}
        .footer {{
            margin-top: 60px;
            text-align: center;
            color: var(--text-muted);
            font-size: 0.85rem;
            border-top: 1px solid var(--border-color);
            padding-top: 20px;
        }}
        
        .tab-alert {{
            background-color: rgba(59, 130, 246, 0.1);
            border: 1px solid rgba(59, 130, 246, 0.2);
            padding: 12px 18px;
            border-radius: 12px;
            font-size: 0.9rem;
            color: #93c5fd;
            margin-bottom: 20px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>Universal {f"Influenza B" if "influenza-b" in build_name.lower() else f"Measles" if "measles" in build_name.lower() else "Influenza A"} Diagnostics Comparison Report</h1>
            <div class="meta">
                Build: <b>{build_name}</b> | Status: <b>Dual-Track Completed</b> | Target segment: <b>Matrix Protein (MP) / Target Gene</b>
            </div>
        </header>

        <!-- Toggle Switch -->
        <div class="toggle-container">
            <button class="toggle-btn active" id="btn-global" onclick="setDatasetMode('global')">
                🌍 Global Surveillance DB
            </button>
            <button class="toggle-btn" id="btn-design" onclick="setDatasetMode('design')">
                🎯 Design DB
            </button>
        </div>

        <div class="grid">
            <div class="card" style="border-top: 4px solid var(--nd-color);">
                <div class="card-title">Non-Degenerate Best Coverage</div>
                <div class="metric-val" id="card-nd-val" style="color: #818cf8;">—%</div>
                <div class="metric-sub" id="card-nd-sub">—</div>
            </div>
            
            <div class="card" style="border-top: 4px solid var(--d-color);">
                <div class="card-title">Degenerate Best Coverage</div>
                <div class="metric-val" id="card-d-val" style="color: #34d399;">—%</div>
                <div class="metric-sub" id="card-d-sub">—</div>
            </div>

            <div class="card" style="border-top: 4px solid var(--warning);">
                <div class="card-title">Coverage Uplift</div>
                <div class="metric-val" id="card-uplift-val" style="color: #fbbf24;">—%</div>
                <div class="metric-sub">Net increase in diagnostic universality</div>
            </div>
        </div>

        <h2>Candidate Comparison Table</h2>
        <table>
            <thead id="table-header">
                <!-- Loaded dynamically by JS -->
            </thead>
            <tbody id="table-body">
                <!-- Loaded dynamically by JS -->
            </tbody>
        </table>

        <div id="subtype-section">
            <h2>Subtype Coverage Breakdown (Rank 1 Candidate)</h2>
            <div class="tab-alert">
                💡 Subtype breakdowns are generated using GISAID EpiFlu metadata, which is only available for the full global Surveillance DB.
            </div>
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
        </div>

        <div class="footer">
            Pandemic Preparedness Pipeline &bull; Automated Diagnostic Design Verification &bull; Dual-Path Analysis
        </div>
    </div>

    <script>
        const nonDegData = {non_deg_json};
        const degData = {deg_json};
        const hasDegenerate = {has_deg_js};

        // Label the dataset-toggle buttons from the actual data — the design DB
        // size varies per build (e.g. indo ~759, all-partial ~6,306) and the
        // global set is not always 50k (Flu B ~7k), so these must not be hardcoded.
        (function() {{
            const c = (nonDegData && nonDegData[0]) || {{}};
            const g = c.total_seqs, d = c.design_db_total;
            if (g != null) document.getElementById('btn-global').innerHTML =
                '🌍 Global Surveillance DB (' + g.toLocaleString() + ')';
            if (d != null) document.getElementById('btn-design').innerHTML =
                '🎯 Design DB (' + d.toLocaleString() + ' seqs)';
        }})();

        let currentMode = 'global';

        function getBadge(diff) {{
            if (diff > 0) {{
                return `<span class="badge positive">+${{diff.toFixed(2)}}%</span>`;
            }} else if (diff < 0) {{
                return `<span class="badge negative">${{diff.toFixed(2)}}%</span>`;
            }} else {{
                return `<span class="badge neutral">0.00%</span>`;
            }}
        }}

        function setDatasetMode(mode) {{
            currentMode = mode;
            
            // Toggle active state on buttons
            document.querySelectorAll('.toggle-btn').forEach(btn => btn.classList.remove('active'));
            if (mode === 'global') {{
                document.getElementById('btn-global').classList.add('active');
                document.getElementById('subtype-section').style.display = 'block';
            }} else {{
                document.getElementById('btn-design').classList.add('active');
                document.getElementById('subtype-section').style.display = 'none';
            }}
            
            renderCards();
            renderTable();
        }}

        function renderCards() {{
            const nd0 = nonDegData[0];
            const d0 = hasDegenerate && degData ? degData[0] : null;
            
            const cardNdVal = document.getElementById('card-nd-val');
            const cardNdSub = document.getElementById('card-nd-sub');
            const cardDVal = document.getElementById('card-d-val');
            const cardDSub = document.getElementById('card-d-sub');
            const cardUpliftVal = document.getElementById('card-uplift-val');
            
            if (currentMode === 'global') {{
                const ndPct = nd0 ? nd0.raa_detected_pct : 0.0;
                const ndSeqs = nd0 ? nd0.raa_detected_seqs : 0;
                const total = nd0 ? nd0.total_seqs : 0;
                cardNdVal.textContent = ndPct + '%';
                cardNdSub.textContent = `Total Covered: ${{ndSeqs.toLocaleString()}} / ${{total.toLocaleString()}} sequences (RAA Det)`;
                
                if (d0) {{
                    const dPct = d0.raa_detected_pct;
                    const dSeqs = d0.raa_detected_seqs;
                    cardDVal.textContent = dPct + '%';
                    cardDSub.textContent = `Total Covered: ${{dSeqs.toLocaleString()}} / ${{total.toLocaleString()}} sequences (RAA Det)`;
                    
                    const diff = dPct - ndPct;
                    cardUpliftVal.textContent = (diff > 0 ? '+' : '') + diff.toFixed(2) + '%';
                }} else {{
                    cardDVal.textContent = 'N/A';
                    cardDSub.textContent = 'Degenerate simulation disabled';
                    cardUpliftVal.textContent = '0.00%';
                }}
            }} else {{
                // Design DB — use design_db_raa_det_pct (computed from design DB sequences)
                const ndPct = nd0 ? (nd0.design_db_raa_det_pct ?? nd0.design_db_amplified_pct) : 0.0;
                const ndSeqs = nd0 ? (nd0.design_db_raa_detected ?? nd0.design_db_amplified) : 0;
                const total = nd0 ? nd0.design_db_total : 0;
                cardNdVal.textContent = ndPct + '%';
                cardNdSub.textContent = `Total Covered: ${{ndSeqs.toLocaleString()}} / ${{total.toLocaleString()}} sequences (RAA Det)`;
                
                if (d0) {{
                    const dPct = d0.design_db_raa_det_pct ?? d0.design_db_amplified_pct;
                    const dSeqs = d0.design_db_raa_detected ?? d0.design_db_amplified;
                    cardDVal.textContent = dPct + '%';
                    cardDSub.textContent = `Total Covered: ${{dSeqs.toLocaleString()}} / ${{total.toLocaleString()}} sequences (RAA Det)`;
                    
                    const diff = dPct - ndPct;
                    cardUpliftVal.textContent = (diff > 0 ? '+' : '') + diff.toFixed(2) + '%';
                }} else {{
                    cardDVal.textContent = 'N/A';
                    cardDSub.textContent = 'Degenerate simulation disabled';
                    cardUpliftVal.textContent = '0.00%';
                }}
            }}
        }}

        function renderTable() {{
            const tableHeader = document.getElementById('table-header');
            const tableBody = document.getElementById('table-body');
            tableBody.innerHTML = '';
            
            // Set headers — always 4 metric columns (PCR Amp, PCR Det, RAA Amp, RAA Det) in both modes
            const modeLabel = currentMode === 'global' ? 'Global' : 'Design DB';
            tableHeader.innerHTML = `
                <tr>
                    <th rowspan="2">Candidate</th>
                    <th colspan="5" class="center-cell" style="background-color: rgba(99, 102, 241, 0.2);">Non-Degenerate Track (${{modeLabel}})</th>
                    <th colspan="5" class="center-cell" style="background-color: rgba(16, 185, 129, 0.2);">Degenerate Track (${{modeLabel}})</th>
                    <th rowspan="2" class="center-cell">Uplift</th>
                </tr>
                <tr>
                    <th style="background-color: rgba(99, 102, 241, 0.1); width: 260px;">Sequences</th>
                    <th class="center-cell" style="background-color: rgba(99, 102, 241, 0.1);">PCR Amp %</th>
                    <th class="center-cell" style="background-color: rgba(99, 102, 241, 0.1);">PCR Det %</th>
                    <th class="center-cell" style="background-color: rgba(99, 102, 241, 0.1);">RAA Amp %</th>
                    <th class="center-cell" style="background-color: rgba(99, 102, 241, 0.1);">RAA Det %</th>
                    <th style="background-color: rgba(16, 185, 129, 0.1); width: 260px;">Sequences</th>
                    <th class="center-cell" style="background-color: rgba(16, 185, 129, 0.1);">PCR Amp %</th>
                    <th class="center-cell" style="background-color: rgba(16, 185, 129, 0.1);">PCR Det %</th>
                    <th class="center-cell" style="background-color: rgba(16, 185, 129, 0.1);">RAA Amp %</th>
                    <th class="center-cell" style="background-color: rgba(16, 185, 129, 0.1);">RAA Det %</th>
                </tr>
            `;
            
            const maxRanks = Math.max(nonDegData.length, hasDegenerate && degData ? degData.length : 0);
            for (let idx = 0; idx < maxRanks; idx++) {{
                const nd = nonDegData[idx];
                const d = hasDegenerate && degData ? degData[idx] : null;
                
                const ndId = nd ? nd.cand_id : '-';
                const dId = d ? d.cand_id : '-';
                
                let ndSeqs = '-';
                let ndPcrAmp = '-', ndPcrDet = '-', ndRaaAmp = '-', ndRaaDet = '-', ndCountStr = '';
                if (nd) {{
                    ndSeqs = `<b>F:</b> ${{nd.fwd_seq}}<br><b>R:</b> ${{nd.rev_seq}}`;
                    if (nd.probe_seq) ndSeqs += `<br><b>P:</b> ${{nd.probe_seq}}`;
                    if (currentMode === 'global') {{
                        ndPcrAmp = (nd.pcr_amplified_pct ?? '-') + '%';
                        ndPcrDet = (nd.pcr_detected_pct ?? '-') + '%';
                        ndRaaAmp = (nd.raa_amplified_pct ?? '-') + '%';
                        ndRaaDet = (nd.raa_detected_pct ?? '-') + '%';
                        ndCountStr = `(${{(nd.raa_detected_seqs ?? 0).toLocaleString()}}/${{(nd.total_seqs ?? 0).toLocaleString()}})`;
                    }} else {{
                        ndPcrAmp = (nd.design_db_amplified_pct ?? '-') + '%';
                        ndPcrDet = ndPcrAmp;  // design DB PCR has no separate det (probe is same run)
                        ndRaaAmp = (nd.design_db_raa_amp_pct ?? '-') + '%';
                        ndRaaDet = (nd.design_db_raa_det_pct ?? '-') + '%';
                        ndCountStr = `(${{(nd.design_db_raa_detected ?? nd.design_db_amplified ?? 0).toLocaleString()}}/${{(nd.design_db_total ?? 0).toLocaleString()}})`;
                    }}
                }}
                
                let dSeqs = '-';
                let dPcrAmp = '-', dPcrDet = '-', dRaaAmp = '-', dRaaDet = '-', dCountStr = '', badge = '-';
                if (hasDegenerate && d) {{
                    dSeqs = `<b>F:</b> ${{d.fwd_seq}}<br><b>R:</b> ${{d.rev_seq}}`;
                    if (d.probe_seq) dSeqs += `<br><b>P:</b> ${{d.probe_seq}}`;
                    if (currentMode === 'global') {{
                        dPcrAmp = (d.pcr_amplified_pct ?? '-') + '%';
                        dPcrDet = (d.pcr_detected_pct ?? '-') + '%';
                        dRaaAmp = (d.raa_amplified_pct ?? '-') + '%';
                        dRaaDet = (d.raa_detected_pct ?? '-') + '%';
                        dCountStr = `(${{(d.raa_detected_seqs ?? 0).toLocaleString()}}/${{(d.total_seqs ?? 0).toLocaleString()}})`;
                        if (nd) {{
                            const diff = (d.raa_detected_pct ?? 0) - (nd.raa_detected_pct ?? 0);
                            badge = getBadge(diff);
                        }}
                    }} else {{
                        dPcrAmp = (d.design_db_amplified_pct ?? '-') + '%';
                        dPcrDet = dPcrAmp;
                        dRaaAmp = (d.design_db_raa_amp_pct ?? '-') + '%';
                        dRaaDet = (d.design_db_raa_det_pct ?? '-') + '%';
                        dCountStr = `(${{(d.design_db_raa_detected ?? d.design_db_amplified ?? 0).toLocaleString()}}/${{(d.design_db_total ?? 0).toLocaleString()}})`;
                        if (nd) {{
                            const diff = (d.design_db_raa_det_pct ?? d.design_db_amplified_pct ?? 0)
                                       - (nd.design_db_raa_det_pct ?? nd.design_db_amplified_pct ?? 0);
                            badge = getBadge(diff);
                        }}
                    }}
                }} else if (!hasDegenerate) {{
                    dSeqs = '<span class="disabled">Degenerate Mode Disabled</span>';
                }}
                
                const tr = document.createElement('tr');
                tr.innerHTML = `
                    <td class="cand-id">Rank ${{idx + 1}}<br><small class="text-muted">ND: ${{ndId}}<br>D: ${{dId}}</small></td>
                    <td class="seq-cell">${{ndSeqs}}</td>
                    <td class="center-cell highlight-nd">${{ndPcrAmp}}</td>
                    <td class="center-cell highlight-nd">${{ndPcrDet}}</td>
                    <td class="center-cell highlight-nd">${{ndRaaAmp}}</td>
                    <td class="center-cell highlight-nd"><b>${{ndRaaDet}}</b><br><small class="text-muted">${{ndCountStr}}</small></td>
                    <td class="seq-cell">${{dSeqs}}</td>
                    <td class="center-cell highlight-d">${{dPcrAmp}}</td>
                    <td class="center-cell highlight-d">${{dPcrDet}}</td>
                    <td class="center-cell highlight-d">${{dRaaAmp}}</td>
                    <td class="center-cell highlight-d"><b>${{dRaaDet}}</b><br><small class="text-muted">${{dCountStr}}</small></td>
                    <td class="center-cell">${{badge}}</td>
                `;
                tableBody.appendChild(tr);
            }}
        }}

        // Initialize display
        setDatasetMode('global');
    </script>
</body>
</html>
"""

    with open(output_path, "w", encoding="utf-8") as f:
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
