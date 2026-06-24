#!/usr/bin/env python3
"""
Helper: Replace generate_heatmap_dashboard in m10_insilico_pcr.py with the premium version.
Run once from the project root, then delete this file.
"""
import sys
from pathlib import Path

target = Path("scripts/m10_insilico_pcr.py")
lines = target.read_text(encoding="utf-8").splitlines(keepends=True)

# Find start and end of the old function (1-indexed → 0-indexed)
start_line = None  # def generate_heatmap_dashboard(...)
end_line   = None  # last line before def generate_markdown_report

for i, l in enumerate(lines):
    if l.startswith("def generate_heatmap_dashboard("):
        start_line = i
    if start_line is not None and l.startswith("def generate_markdown_report("):
        end_line = i
        break

if start_line is None or end_line is None:
    print(f"ERROR: could not find function boundaries (start={start_line}, end={end_line})")
    sys.exit(1)

print(f"Replacing lines {start_line+1}..{end_line} (0-indexed) with new function")

NEW_FUNC = r'''def generate_heatmap_dashboard(summary_data: list, output_path):
    """Generate a premium interactive dark-mode HTML coverage dashboard with full PCR and RAA data."""
    import json as _json
    data_json = _json.dumps(summary_data, default=str)

    html_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>\U0001f9ec In-silico PCR/RAA Coverage Dashboard</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
    <style>
        :root {
            --bg:#080f1e;--surface:#0f1e35;--surface2:#162036;
            --card-bg:rgba(22,33,62,0.85);--border:rgba(255,255,255,0.07);
            --accent-pcr:#10b981;--accent-pcr-dim:rgba(16,185,129,0.15);
            --accent-raa:#6366f1;--accent-raa-dim:rgba(99,102,241,0.15);
            --accent-design:#a78bfa;--text-1:#f0f6ff;--text-2:#8b9dc0;
        }
        *{box-sizing:border-box;margin:0;padding:0;}
        body{font-family:'Outfit',sans-serif;background:var(--bg);color:var(--text-1);padding:2rem;}
        .dashboard{max-width:1500px;margin:0 auto;}
        header{margin-bottom:2rem;padding-bottom:1.5rem;border-bottom:1px solid var(--border);}
        .header-flex{display:flex;align-items:flex-start;justify-content:space-between;flex-wrap:wrap;gap:1rem;}
        h1{font-size:2rem;font-weight:700;background:linear-gradient(135deg,#34d399 0%,#6366f1 100%);-webkit-background-clip:text;-webkit-text-fill-color:transparent;margin-bottom:0.4rem;}
        .subtitle{color:var(--text-2);font-size:0.95rem;line-height:1.6;}
        .toggle-container{display:inline-flex;background:var(--card-bg);border:1px solid var(--border);border-radius:14px;padding:4px;box-shadow:0 8px 20px rgba(0,0,0,0.3);}
        .toggle-btn{background:none;border:none;color:var(--text-2);padding:10px 22px;border-radius:10px;font-weight:600;cursor:pointer;transition:all 0.25s ease;font-family:inherit;font-size:0.9rem;display:flex;align-items:center;gap:8px;white-space:nowrap;}
        .toggle-btn.active-global{background:var(--accent-pcr);color:#0a1628;box-shadow:0 4px 16px rgba(16,185,129,0.3);}
        .toggle-btn.active-design{background:var(--accent-design);color:#0a1628;box-shadow:0 4px 16px rgba(167,139,250,0.3);}
        .toggle-btn:hover:not(.active-global):not(.active-design){color:var(--text-1);background:rgba(255,255,255,0.05);}
        .controls{display:flex;align-items:center;gap:1rem;margin-bottom:2rem;background:var(--card-bg);padding:1rem 1.5rem;border-radius:14px;border:1px solid var(--border);flex-wrap:wrap;}
        .controls label{font-weight:600;color:var(--text-2);font-size:0.9rem;}
        select{background:var(--surface2);color:var(--text-1);border:1px solid var(--border);padding:0.5rem 1rem;border-radius:8px;font-family:inherit;font-size:0.95rem;cursor:pointer;outline:none;flex:1;min-width:280px;}
        select:focus{border-color:var(--accent-pcr);}
        .dual-metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:1.5rem;margin-bottom:2rem;}
        .metric-card{background:var(--card-bg);border:1px solid var(--border);border-radius:16px;padding:1.5rem;backdrop-filter:blur(10px);transition:border-color 0.3s,transform 0.3s;}
        .metric-card:hover{border-color:rgba(16,185,129,0.35);transform:translateY(-3px);}
        .metric-label{font-size:0.8rem;color:var(--text-2);text-transform:uppercase;letter-spacing:0.08em;margin-bottom:0.8rem;}
        .metric-dual{display:flex;align-items:baseline;gap:1.5rem;}
        .metric-block{flex:1;}
        .metric-mode-label{font-size:0.72rem;font-weight:600;text-transform:uppercase;letter-spacing:0.06em;margin-bottom:0.2rem;}
        .metric-value{font-size:1.7rem;font-weight:700;}
        .pcr-color{color:var(--accent-pcr);}
        .raa-color{color:var(--accent-raa);}
        .divider-v{width:1px;background:var(--border);height:50px;align-self:center;}
        .prog-wrap{height:4px;background:rgba(255,255,255,0.05);border-radius:9px;margin-top:0.6rem;overflow:hidden;}
        .prog-bar{height:100%;border-radius:9px;transition:width 0.6s ease;}
        .seq-grid{display:grid;grid-template-columns:1fr 1fr 1fr;gap:1rem;margin-bottom:2rem;}
        @media(max-width:1100px){.seq-grid{grid-template-columns:1fr 1fr;}}
        @media(max-width:700px){.seq-grid{grid-template-columns:1fr;}}
        .seq-card{background:var(--card-bg);border:1px solid var(--border);border-radius:14px;padding:1.2rem;}
        .seq-title{font-size:0.78rem;font-weight:600;text-transform:uppercase;letter-spacing:0.07em;color:var(--text-2);margin-bottom:0.6rem;}
        .seq-mono{font-family:'JetBrains Mono',monospace;font-size:0.82rem;background:#060d1a;padding:0.6rem 0.8rem;border-radius:8px;border:1px solid rgba(255,255,255,0.04);word-break:break-all;line-height:1.6;color:#a5f3c3;}
        .seq-deg-note{font-size:0.72rem;color:var(--text-2);margin-top:0.4rem;}
        .design-banner{background:linear-gradient(135deg,rgba(139,92,246,0.12),rgba(109,40,217,0.08));border:1px solid rgba(139,92,246,0.3);border-radius:16px;padding:1.2rem 1.5rem;margin-bottom:2rem;}
        .design-banner-header{font-size:0.75rem;font-weight:700;text-transform:uppercase;letter-spacing:0.08em;color:#a78bfa;margin-bottom:0.8rem;}
        .design-content{display:flex;gap:2.5rem;flex-wrap:wrap;align-items:center;}
        .design-stat-label{font-size:0.72rem;color:rgba(167,139,250,0.65);margin-bottom:0.2rem;}
        .design-stat-value{font-size:1.5rem;font-weight:700;color:#c4b5fd;}
        .design-bar-wrap{flex:1;min-width:120px;}
        .design-bar-track{background:rgba(255,255,255,0.05);border-radius:9px;height:8px;overflow:hidden;}
        .design-bar-fill{height:100%;border-radius:9px;transition:width 0.6s ease;}
        .main-layout{display:grid;grid-template-columns:1fr 1fr;gap:1.5rem;}
        @media(max-width:900px){.main-layout{grid-template-columns:1fr;}.full-width{grid-column:span 1!important;}}
        .full-width{grid-column:span 2;}
        .card{background:var(--card-bg);border:1px solid var(--border);border-radius:16px;padding:1.5rem;backdrop-filter:blur(10px);}
        .card-title{font-size:1.05rem;font-weight:600;margin-bottom:1.2rem;padding-bottom:0.75rem;border-bottom:1px solid var(--border);display:flex;justify-content:space-between;align-items:center;}
        .badge{padding:0.2rem 0.6rem;border-radius:6px;font-size:0.75rem;font-weight:600;}
        .badge-pcr{background:var(--accent-pcr-dim);color:var(--accent-pcr);border:1px solid rgba(16,185,129,0.3);}
        .badge-raa{background:var(--accent-raa-dim);color:var(--accent-raa);border:1px solid rgba(99,102,241,0.3);}
        table{width:100%;border-collapse:collapse;font-size:0.88rem;}
        th,td{padding:0.6rem 0.75rem;text-align:left;border-bottom:1px solid rgba(255,255,255,0.04);}
        th{color:var(--text-2);font-size:0.75rem;font-weight:700;text-transform:uppercase;letter-spacing:0.05em;}
        tr:hover td{background:rgba(255,255,255,0.02);}
        .cov-cell{font-weight:700;border-radius:6px;padding:0.35rem 0.55rem;display:inline-block;transition:transform 0.15s;font-size:0.85rem;}
        .cov-cell:hover{transform:scale(1.08);}
        .mini-bar{width:100%;background:rgba(255,255,255,0.06);border-radius:4px;height:4px;margin-top:3px;overflow:hidden;}
        .mini-bar-fill{height:100%;border-radius:4px;}
    </style>
</head>
<body>
<div class="dashboard">
    <header>
        <div class="header-flex">
            <div>
                <h1>\U0001f9ec PCR &amp; RAA Coverage Dashboard</h1>
                <div class="subtitle">
                    Dual-mode in-silico validation: <span style="color:var(--accent-pcr);font-weight:600;">PCR</span> (seqkit amplicon)
                    &amp; <span style="color:var(--accent-raa);font-weight:600;">RAA/RPA</span> (recombinase binding model).
                </div>
            </div>
            <div class="toggle-container">
                <button class="toggle-btn active-global" id="btn-global" onclick="setDatasetMode('global')">\U0001f30d Global Surveillance DB</button>
                <button class="toggle-btn" id="btn-design" onclick="setDatasetMode('design')">\U0001f3af Design DB Only</button>
            </div>
        </div>
    </header>

    <div class="controls">
        <label for="candidateSelect">Select Candidate:</label>
        <select id="candidateSelect" onchange="updateDashboard(this.value)"></select>
    </div>

    <div class="dual-metrics">
        <div class="metric-card">
            <div class="metric-label">Candidate &middot; Region &middot; Rank</div>
            <div id="statRegion" style="font-size:1.2rem;font-weight:700;">&mdash;</div>
        </div>
        <div class="metric-card">
            <div class="metric-label">Total Sequences (DB)</div>
            <div id="statTotal" style="font-size:1.7rem;font-weight:700;">&mdash;</div>
        </div>
        <div class="metric-card">
            <div class="metric-label">\U0001f52c Amplification Coverage</div>
            <div class="metric-dual">
                <div class="metric-block">
                    <div class="metric-mode-label pcr-color">PCR</div>
                    <div class="metric-value pcr-color" id="pcrAmpPct">&mdash;%</div>
                    <div class="prog-wrap"><div class="prog-bar" id="pcrAmpBar" style="background:var(--accent-pcr);width:0%;"></div></div>
                </div>
                <div class="divider-v"></div>
                <div class="metric-block">
                    <div class="metric-mode-label raa-color">RAA</div>
                    <div class="metric-value raa-color" id="raaAmpPct">&mdash;</div>
                    <div class="prog-wrap"><div class="prog-bar" id="raaAmpBar" style="background:var(--accent-raa);width:0%;"></div></div>
                </div>
            </div>
        </div>
        <div class="metric-card">
            <div class="metric-label">\U0001f3af Detection Coverage (+Probe)</div>
            <div class="metric-dual">
                <div class="metric-block">
                    <div class="metric-mode-label pcr-color">PCR</div>
                    <div class="metric-value pcr-color" id="pcrDetPct">&mdash;%</div>
                    <div class="prog-wrap"><div class="prog-bar" id="pcrDetBar" style="background:var(--accent-pcr);width:0%;"></div></div>
                </div>
                <div class="divider-v"></div>
                <div class="metric-block">
                    <div class="metric-mode-label raa-color">RAA</div>
                    <div class="metric-value raa-color" id="raaDetPct">&mdash;</div>
                    <div class="prog-wrap"><div class="prog-bar" id="raaDetBar" style="background:var(--accent-raa);width:0%;"></div></div>
                </div>
            </div>
        </div>
    </div>

    <div class="design-banner" id="designBanner" style="display:none;">
        <div class="design-banner-header">\U0001f3af Design DB Coverage (sekuens yang dipakai desain primer &mdash; M02 filtered)</div>
        <div class="design-content">
            <div><div class="design-stat-label">Total Sequences</div><div class="design-stat-value" id="designTotal">&mdash;</div></div>
            <div><div class="design-stat-label">PCR Amplified</div><div class="design-stat-value" id="designAmplified">&mdash;</div></div>
            <div><div class="design-stat-label">Coverage</div><div class="design-stat-value" id="designPct">&mdash;%</div></div>
            <div class="design-bar-wrap">
                <div class="design-stat-label">Coverage Bar</div>
                <div class="design-bar-track"><div class="design-bar-fill" id="designBar" style="width:0%;background:linear-gradient(90deg,#7c3aed,#a78bfa);"></div></div>
            </div>
        </div>
    </div>

    <div class="seq-grid">
        <div class="seq-card">
            <div class="seq-title">&#x2B06; Forward Primer (5&#x27;&#x2192;3&#x27;) &mdash; Majority</div>
            <div class="seq-mono" id="seqFwd">&mdash;</div>
            <div class="seq-deg-note" id="seqFwdDeg" style="display:none;"></div>
        </div>
        <div class="seq-card">
            <div class="seq-title">&#x2B07; Reverse Primer (5&#x27;&#x2192;3&#x27;) &mdash; Majority</div>
            <div class="seq-mono" id="seqRev">&mdash;</div>
            <div class="seq-deg-note" id="seqRevDeg" style="display:none;"></div>
        </div>
        <div class="seq-card">
            <div class="seq-title">&#x1F534; Probe Sequence (5&#x27;&#x2192;3&#x27;)</div>
            <div class="seq-mono" id="seqProbe">&mdash;</div>
            <div class="seq-deg-note" id="seqProbeDeg" style="display:none;"></div>
        </div>
    </div>

    <div id="globalTables">
        <div class="main-layout">
            <div class="card">
                <div class="card-title">
                    Subtype Coverage Distribution
                    <div style="display:flex;gap:0.5rem;">
                        <span class="badge badge-pcr">PCR</span>
                        <span class="badge badge-raa">RAA</span>
                    </div>
                </div>
                <table id="subtypeTable">
                    <thead><tr>
                        <th>Subtype</th><th>N</th>
                        <th style="color:var(--accent-pcr);">PCR Amp</th>
                        <th style="color:var(--accent-pcr);">PCR Det</th>
                        <th style="color:var(--accent-raa);">RAA Amp</th>
                        <th style="color:var(--accent-raa);">RAA Det</th>
                    </tr></thead>
                    <tbody></tbody>
                </table>
            </div>
            <div class="card">
                <div class="card-title">
                    Collection Year Breakdown
                    <span class="badge badge-pcr">Temporal</span>
                </div>
                <table id="yearTable">
                    <thead><tr>
                        <th>Year</th><th>N</th>
                        <th style="color:var(--accent-pcr);">PCR Amp</th>
                        <th style="color:var(--accent-pcr);">PCR Det</th>
                        <th style="color:var(--accent-raa);">RAA Amp</th>
                        <th style="color:var(--accent-raa);">RAA Det</th>
                    </tr></thead>
                    <tbody></tbody>
                </table>
            </div>
            <div class="card full-width">
                <div class="card-title">
                    Country Coverage Breakdown
                    <span class="badge badge-raa">Geographical</span>
                </div>
                <table id="countryTable">
                    <thead><tr>
                        <th>Country</th><th>N</th>
                        <th style="color:var(--accent-pcr);">PCR Amp</th>
                        <th style="color:var(--accent-pcr);">PCR Det</th>
                        <th style="color:var(--accent-raa);">RAA Amp</th>
                        <th style="color:var(--accent-raa);">RAA Det</th>
                    </tr></thead>
                    <tbody></tbody>
                </table>
            </div>
        </div>
    </div>
</div>

<script>
""" + "    const candidatesData = " + data_json + """;
    let currentMode = 'global';
    const sel = document.getElementById('candidateSelect');
    candidatesData.forEach((c, i) => {
        const opt = document.createElement('option');
        opt.value = i;
        const pa = c.pcr_amplified_pct ?? 0;
        const ra = c.raa_amplified_pct ?? '\u2014';
        const nc = c.n_combinations > 1 ? ` [${c.n_combinations} combos]` : '';
        opt.textContent = `${c.cand_id}${nc}  |  PCR Amp: ${pa}%  \u00b7  RAA Amp: ${ra}%`;
        sel.appendChild(opt);
    });
    function pctStyle(pct) {
        const v = parseFloat(pct || 0);
        if (v >= 80) return 'background:rgba(52,211,153,0.85);color:#052e1c;';
        if (v >= 60) return 'background:rgba(16,185,129,0.7);color:#052e1c;';
        if (v >= 35) return 'background:rgba(251,191,36,0.85);color:#1a0e00;';
        if (v >= 15) return 'background:rgba(251,146,60,0.85);color:#1a0500;';
        return 'background:rgba(248,113,113,0.85);color:#1a0000;';
    }
    function barColor(pct, type) {
        const v = parseFloat(pct || 0);
        if (type === 'raa') { return v >= 70 ? '#6366f1' : v >= 40 ? '#818cf8' : '#f87171'; }
        return v >= 70 ? '#10b981' : v >= 40 ? '#fbbf24' : '#f87171';
    }
    function makeCell(pct, type) {
        const td = document.createElement('td');
        const span = document.createElement('span');
        span.className = 'cov-cell';
        span.style.cssText = pctStyle(pct);
        span.textContent = pct !== null && pct !== undefined ? `${pct}%` : '\u2014';
        const bar = document.createElement('div'); bar.className = 'mini-bar';
        const fill = document.createElement('div'); fill.className = 'mini-bar-fill';
        fill.style.width = `${Math.min(parseFloat(pct||0), 100)}%`;
        fill.style.background = barColor(pct, type);
        bar.appendChild(fill);
        const w = document.createElement('div');
        w.appendChild(span); w.appendChild(bar);
        td.appendChild(w);
        return td;
    }
    function populateTable(tid, data, catKey, maxRows) {
        const tbody = document.querySelector(`#${tid} tbody`);
        tbody.innerHTML = '';
        if (!data) return;
        const rows = maxRows ? data.slice(0, maxRows) : data;
        rows.forEach(row => {
            const tr = document.createElement('tr');
            const tdCat = document.createElement('td'); tdCat.textContent = row[catKey] || '\u2014'; tr.appendChild(tdCat);
            const tdN = document.createElement('td'); tdN.textContent = (row.total||0).toLocaleString(); tdN.style.color='var(--text-2)'; tr.appendChild(tdN);
            tr.appendChild(makeCell(row.pcr_amplified_pct ?? row.amplified_pct, 'pcr'));
            tr.appendChild(makeCell(row.pcr_detected_pct ?? row.detected_pct, 'pcr'));
            tr.appendChild(makeCell(row.raa_amplified_pct, 'raa'));
            tr.appendChild(makeCell(row.raa_detected_pct, 'raa'));
            tbody.appendChild(tr);
        });
        if (maxRows && data.length > maxRows) {
            const tr = document.createElement('tr');
            const td = document.createElement('td'); td.colSpan = 6;
            td.style.cssText = 'text-align:center;color:var(--text-2);font-style:italic;padding:0.5rem;';
            td.textContent = `Showing ${maxRows} of ${data.length} entries`;
            tr.appendChild(td); tbody.appendChild(tr);
        }
    }
    function setDatasetMode(mode) {
        currentMode = mode;
        document.querySelectorAll('.toggle-btn').forEach(b => b.classList.remove('active-global','active-design'));
        if (mode === 'global') {
            document.getElementById('btn-global').classList.add('active-global');
            document.getElementById('globalTables').style.display = 'block';
        } else {
            document.getElementById('btn-design').classList.add('active-design');
            document.getElementById('globalTables').style.display = 'none';
        }
        updateDashboard(sel.value);
    }
    function updateDashboard(idx) {
        const c = candidatesData[idx]; if (!c) return;
        document.getElementById('statRegion').textContent = `${c.region} \u00b7 Rank ${c.rank}`;
        const total = currentMode === 'global' ? (c.total_seqs||0) : (c.design_db_total||0);
        document.getElementById('statTotal').textContent = total.toLocaleString();
        let pcrAmp=0,pcrDet=0,raaAmp=null,raaDet=null;
        if (currentMode === 'global') {
            pcrAmp=c.pcr_amplified_pct??0; pcrDet=c.pcr_detected_pct??0;
            raaAmp=c.raa_amplified_pct; raaDet=c.raa_detected_pct;
        } else { pcrAmp=c.design_db_amplified_pct??0; pcrDet=pcrAmp; }
        document.getElementById('pcrAmpPct').textContent=`${pcrAmp}%`;
        document.getElementById('pcrDetPct').textContent=`${pcrDet}%`;
        document.getElementById('raaAmpPct').textContent=raaAmp!=null?`${raaAmp}%`:'\u2014';
        document.getElementById('raaDetPct').textContent=raaDet!=null?`${raaDet}%`:'\u2014';
        [['pcrAmpBar',pcrAmp,'pcr'],['pcrDetBar',pcrDet,'pcr'],['raaAmpBar',raaAmp,'raa'],['raaDetBar',raaDet,'raa']].forEach(([id,v,t])=>{
            const el=document.getElementById(id);
            el.style.width=`${Math.min(parseFloat(v||0),100)}%`;
            el.style.background=barColor(v,t);
        });
        const banner=document.getElementById('designBanner');
        if (c.design_db_total>0 && currentMode==='global') {
            banner.style.display='block';
            document.getElementById('designTotal').textContent=c.design_db_total.toLocaleString();
            document.getElementById('designAmplified').textContent=c.design_db_amplified.toLocaleString();
            const dp=c.design_db_amplified_pct||0;
            document.getElementById('designPct').textContent=`${dp}%`;
            const db=document.getElementById('designBar');
            db.style.width=`${Math.min(dp,100)}%`;
            db.style.background=dp>=80?'linear-gradient(90deg,#059669,#34d399)':dp>=50?'linear-gradient(90deg,#7c3aed,#a78bfa)':dp>=30?'linear-gradient(90deg,#d97706,#fbbf24)':'linear-gradient(90deg,#dc2626,#f87171)';
        } else { banner.style.display='none'; }
        document.getElementById('seqFwd').textContent=c.fwd_seq||'\u2014';
        document.getElementById('seqRev').textContent=c.rev_seq||'\u2014';
        document.getElementById('seqProbe').textContent=c.probe_seq||'\u2014 (No Probe)';
        const fwdDeg=document.getElementById('seqFwdDeg');
        const revDeg=document.getElementById('seqRevDeg');
        const prbDeg=document.getElementById('seqProbeDeg');
        if (c.fwd_degenerate_seq && c.fwd_degenerate_seq!==c.fwd_seq) {
            fwdDeg.innerHTML=`\U0001f9ec IUPAC: <span style="font-family:'JetBrains Mono',monospace;color:#c4b5fd;font-size:0.8rem;">${c.fwd_degenerate_seq}</span>`+(c.n_combinations>1?` \u00b7 ${c.n_combinations} combos`:'');
            fwdDeg.style.display='block';
        } else { fwdDeg.style.display='none'; }
        if (c.rev_degenerate_seq && c.rev_degenerate_seq!==c.rev_seq) {
            revDeg.innerHTML=`\U0001f9ec IUPAC: <span style="font-family:'JetBrains Mono',monospace;color:#c4b5fd;font-size:0.8rem;">${c.rev_degenerate_seq}</span>`;
            revDeg.style.display='block';
        } else { revDeg.style.display='none'; }
        if (c.probe_seq && /[RYMKSWBDHVN]/i.test(c.probe_seq)) {
            prbDeg.textContent='\U0001f9ec Contains IUPAC degenerate bases';
            prbDeg.style.display='block';
        } else { prbDeg.style.display='none'; }
        if (currentMode==='global' && c.breakdowns) {
            populateTable('subtypeTable',c.breakdowns.subtype,'subtype');
            populateTable('yearTable',c.breakdowns.year,'year');
            populateTable('countryTable',c.breakdowns.country,'country',30);
        }
    }
    if (candidatesData.length > 0) updateDashboard(0);
</script>
</body>
</html>
"""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    print(f"Generated Heatmap Dashboard: {output_path}")


'''

# Write to file: keep lines before start_line and from end_line onward
new_lines = lines[:start_line] + [NEW_FUNC + "\n"] + lines[end_line:]
target.write_text("".join(new_lines), encoding="utf-8")
print(f"Done. New file has {len(new_lines)} lines.")
print("Verifying syntax...")
import ast
ast.parse(target.read_text(encoding="utf-8"))
print("Syntax OK.")
