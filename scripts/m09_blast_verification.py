import os
import sys
import json
import argparse
import subprocess
import tempfile
import csv
from pathlib import Path
from datetime import datetime
from typing import Optional

def main():
    parser = argparse.ArgumentParser(description="M09 Specificity Verification for Top N Candidates")
    parser.add_argument("--alternatives", required=True, help="Path to alternatives.json")
    parser.add_argument("--num-candidates", type=int, default=3, help="Number of candidates to verify")
    parser.add_argument("--database", required=True, help="BLAST database path (or nt/nr)")
    parser.add_argument("--mode", default="remote", choices=["remote", "local"], help="BLAST mode")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--target", default=None, help="Expected target sequence ID (to distinguish from off-targets)")
    parser.add_argument("--pathogen", default=None, help="Pathogen name keyword (e.g. 'Influenza') to group family hits")
    parser.add_argument(
        "--db-label", default=None,
        help="Label for this database run (e.g. 'human' or 'virus'). "
             "Used to name output subdirs when querying multiple databases in local mode. "
             "If not set, files use default names (blast_result.json, specificity_report.md)."
    )
    
    args = parser.parse_args()
    
    if not os.path.exists(args.alternatives):
        print(f"❌ Alternatives file not found: {args.alternatives}")
        sys.exit(1)
        
    with open(args.alternatives, "r") as f:
        candidates = json.load(f)
        
    if not candidates:
        print("⚠️ No candidate primers found in alternatives.json.")
        if not args.db_label:
            marker_path = Path(args.outdir) / "verification_complete.txt"
            marker_path.parent.mkdir(parents=True, exist_ok=True)
            with open(marker_path, "w") as f:
                f.write("No candidates found.\n")
        sys.exit(0)
        
    # Limit number of candidates to process
    to_process = candidates[:args.num_candidates]
    db_label_str = f" [{args.db_label} db]" if args.db_label else ""
    print(f"🚀 Found {len(candidates)} candidates. Processing top {len(to_process)}...{db_label_str}")
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    for idx, cand in enumerate(to_process):
        rank = cand.get("final_rank") or (idx + 1)
        print(f"\n--- Processing Candidate #{rank}{db_label_str} ---")
        
        primers = cand.get("primers", {})
        fwd = primers.get("forward", {})
        rev = primers.get("reverse", {})
        probe = primers.get("probe", {})
        
        if not fwd or not rev:
            print(f"⚠️ Candidate #{rank} has incomplete primer data. Skipping.")
            continue
            
        # Create a temp file containing the primers dict for this candidate
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as temp_f:
            json.dump(primers, temp_f)
            temp_json_path = temp_f.name
            
        try:
            cand_out_dir = outdir / f"candidate_{rank}"
            cand_out_dir.mkdir(parents=True, exist_ok=True)
            
            # When db_label is set, use a subdirectory to avoid overwriting previous db results
            if args.db_label:
                label_out_dir = cand_out_dir / f"db_{args.db_label}"
                label_out_dir.mkdir(parents=True, exist_ok=True)
                output_arg = str(label_out_dir)
            else:
                output_arg = str(cand_out_dir)
            
            # Build command line
            cmd = [
                "primerlab", "blast",
                "--primers", temp_json_path,
                "--output", output_arg
            ]
            
            if args.target:
                cmd.extend(["--target", args.target])
            if args.pathogen:
                cmd.extend(["--pathogen", args.pathogen])
            
            if args.mode == "remote":
                cmd.append("--online")
                cmd.extend(["--database", args.database])
            else:
                cmd.extend(["--database", args.database])
                
            print(f"Executing: {' '.join(cmd)}")
            
            # Run command and stream output in real-time
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1
            )
            
            if process.stdout:
                for line in process.stdout:
                    print(f"   {line.strip()}")
                    sys.stdout.flush()
                
            process.wait()
            
            if process.returncode == 0:
                print(f"✅ Candidate #{rank}{db_label_str} verification complete.")
            else:
                print(f"❌ Candidate #{rank}{db_label_str} BLAST verification failed with exit code {process.returncode}.")
                
        except Exception as e:
            print(f"❌ Error verifying candidate #{rank}{db_label_str}: {e}")
        finally:
            if os.path.exists(temp_json_path):
                os.remove(temp_json_path)
                
    # Write marker file only if this is a standalone run (no db_label / remote mode).
    # In local dual-db mode, the Snakemake shell writes the marker after both db runs.
    if not args.db_label:
        marker_path = outdir / "verification_complete.txt"
        with open(marker_path, "w") as f:
            f.write(f"Completed verification for {len(to_process)} candidates at {Path(args.alternatives).name}\n")
        print(f"\n🎉 Verification complete. Marker file written: {marker_path}")
    else:
        print(f"\n✅ DB run [{args.db_label}] complete. Marker will be written by Snakemake after all DB runs.")

    # Compile overall summary report and CSV for downstream decisions
    try:
        compile_summary(str(outdir), args.alternatives, args.target)
    except Exception as e:
        import traceback
        print(f"⚠️ Warning: Failed to compile specificity summary: {e}")
        traceback.print_exc()

def get_grade(score):
    if score >= 90: return "A"
    elif score >= 80: return "B"
    elif score >= 70: return "C"
    elif score >= 60: return "D"
    else: return "F"

def compile_summary(outdir_path: str, alternatives_path: Optional[str] = None, target_id: Optional[str] = None):
    """
    Scans output directory, loads JSON results for all candidates,
    and generates master CSV and Markdown reports for comparison.
    Also produces dedicated BED, GFF3, JSON and MD reports specifically for PASSing candidates.
    """
    outdir = Path(outdir_path)
    if not outdir.exists():
        return

    # Load alternatives sequences as fallback
    alt_sequences = {}
    alts = []
    if alternatives_path and os.path.exists(alternatives_path):
        try:
            with open(alternatives_path, "r") as f:
                alts = json.load(f)
            for idx, alt in enumerate(alts):
                rank = alt.get("final_rank") or (idx + 1)
                primers = alt.get("primers", {})
                fwd_seq = primers.get("forward", {}).get("sequence", "")
                rev_seq = primers.get("reverse", {}).get("sequence", "")
                probe_seq = primers.get("probe", {}).get("sequence", "")
                alt_sequences[rank] = {
                    "fwd": fwd_seq,
                    "rev": rev_seq,
                    "probe": probe_seq,
                    "full_alt": alt
                }
        except Exception as e:
            print(f"Warning: Failed to parse alternatives fallback: {e}")

    summary_data = []

    # Find all candidate subdirectories
    for cand_dir in sorted(outdir.glob("candidate_*"), key=lambda p: int(p.name.split("_")[1])):
        rank = int(cand_dir.name.split("_")[1])
        
        fwd_seq, rev_seq, probe_seq = "", "", ""
        
        human_score, human_grade, human_status = "-", "-", "-"
        human_offtargets, human_dangerous = "-", "-"
        
        virus_score, virus_grade, virus_status = "-", "-", "-"
        virus_pathogen, virus_offtargets = "-", "-"
        
        remote_score, remote_grade, remote_status = "-", "-", "-"
        remote_offtargets, remote_dangerous, remote_pathogen = "-", "-", "-"

        # Load local human result
        human_json = cand_dir / "db_human" / "blast_result.json"
        if human_json.exists():
            try:
                with open(human_json, "r") as f:
                    data = json.load(f)
                fwd_res = data.get("forward_result") or data.get("forward") or {}
                rev_res = data.get("reverse_result") or data.get("reverse") or {}
                probe_res = data.get("probe_result") or data.get("probe") or {}
                fwd_seq = fwd_res.get("primer_seq", "")
                rev_seq = rev_res.get("primer_seq", "")
                probe_seq = probe_res.get("primer_seq", "") if probe_res else ""
                
                human_score = data.get("combined_score", 0.0)
                human_grade = get_grade(human_score)
                human_status = "Specific" if data.get("is_specific", True) else "Low Specificity"
                
                if "forward_result" in data:
                    fwd_ots = data.get("forward_result", {}).get("offtargets", [])
                    rev_ots = data.get("reverse_result", {}).get("offtargets", [])
                else:
                    fwd_ots = data.get("forward", {}).get("offtargets", 0)
                    rev_ots = data.get("reverse", {}).get("offtargets", 0)

                fwd_cnt = len(fwd_ots) if isinstance(fwd_ots, list) else int(fwd_ots)
                rev_cnt = len(rev_ots) if isinstance(rev_ots, list) else int(rev_ots)
                human_offtargets = fwd_cnt + rev_cnt
                
                fwd_list = fwd_ots if isinstance(fwd_ots, list) else []
                rev_list = rev_ots if isinstance(rev_ots, list) else []
                human_dangerous = sum(1 for ot in fwd_list + rev_list if ot.get("three_prime_involved", False))
            except Exception as e:
                print(f"Warning: Failed to parse human JSON for candidate {rank}: {e}")

        # Load local virus result
        virus_json = cand_dir / "db_virus" / "blast_result.json"
        if virus_json.exists():
            try:
                with open(virus_json, "r") as f:
                    data = json.load(f)
                if not fwd_seq:
                    fwd_res = data.get("forward_result") or data.get("forward") or {}
                    rev_res = data.get("reverse_result") or data.get("reverse") or {}
                    probe_res = data.get("probe_result") or data.get("probe") or {}
                    fwd_seq = fwd_res.get("primer_seq", "")
                    rev_seq = rev_res.get("primer_seq", "")
                    probe_seq = probe_res.get("primer_seq", "") if probe_res else ""
                
                virus_score = data.get("combined_score", 0.0)
                virus_grade = get_grade(virus_score)
                virus_status = "Specific" if data.get("is_specific", True) else "Low Specificity"
                
                if "forward_result" in data:
                    fwd_ots = data.get("forward_result", {}).get("offtargets", [])
                    rev_ots = data.get("reverse_result", {}).get("offtargets", [])
                    fwd_path = data.get("forward_result", {}).get("pathogen_hits", [])
                    rev_path = data.get("reverse_result", {}).get("pathogen_hits", [])
                else:
                    fwd_ots = data.get("forward", {}).get("offtargets", 0)
                    rev_ots = data.get("reverse", {}).get("offtargets", 0)
                    fwd_path = data.get("forward", {}).get("pathogen_family_hits", 0)
                    rev_path = data.get("reverse", {}).get("pathogen_family_hits", 0)

                fwd_cnt = len(fwd_ots) if isinstance(fwd_ots, list) else int(fwd_ots)
                rev_cnt = len(rev_ots) if isinstance(rev_ots, list) else int(rev_ots)
                virus_offtargets = fwd_cnt + rev_cnt
                
                fwd_path_cnt = len(fwd_path) if isinstance(fwd_path, list) else int(fwd_path)
                rev_path_cnt = len(rev_path) if isinstance(rev_path, list) else int(rev_path)
                virus_pathogen = fwd_path_cnt + rev_path_cnt
            except Exception as e:
                print(f"Warning: Failed to parse virus JSON for candidate {rank}: {e}")

        # Load remote result
        remote_json = cand_dir / "blast_result.json"
        if remote_json.exists():
            try:
                with open(remote_json, "r") as f:
                    data = json.load(f)
                fwd_res = data.get("forward_result") or data.get("forward") or {}
                rev_res = data.get("reverse_result") or data.get("reverse") or {}
                probe_res = data.get("probe_result") or data.get("probe") or {}
                fwd_seq = fwd_res.get("primer_seq", "")
                rev_seq = rev_res.get("primer_seq", "")
                probe_seq = probe_res.get("primer_seq", "") if probe_res else ""
                
                remote_score = data.get("combined_score", 0.0)
                remote_grade = get_grade(remote_score)
                remote_status = "Specific" if data.get("is_specific", True) else "Low Specificity"
                
                if "forward_result" in data:
                    fwd_ots = data.get("forward_result", {}).get("offtargets", [])
                    rev_ots = data.get("reverse_result", {}).get("offtargets", [])
                    fwd_path = data.get("forward_result", {}).get("pathogen_hits", [])
                    rev_path = data.get("reverse_result", {}).get("pathogen_hits", [])
                else:
                    fwd_ots = data.get("forward", {}).get("offtargets", 0)
                    rev_ots = data.get("reverse", {}).get("offtargets", 0)
                    fwd_path = data.get("forward", {}).get("pathogen_family_hits", 0)
                    rev_path = data.get("reverse", {}).get("pathogen_family_hits", 0)

                fwd_cnt = len(fwd_ots) if isinstance(fwd_ots, list) else int(fwd_ots)
                rev_cnt = len(rev_ots) if isinstance(rev_ots, list) else int(rev_ots)
                remote_offtargets = fwd_cnt + rev_cnt
                
                fwd_list = fwd_ots if isinstance(fwd_ots, list) else []
                rev_list = rev_ots if isinstance(rev_ots, list) else []
                remote_dangerous = sum(1 for ot in fwd_list + rev_list if ot.get("three_prime_involved", False))
                
                fwd_path_cnt = len(fwd_path) if isinstance(fwd_path, list) else int(fwd_path)
                rev_path_cnt = len(rev_path) if isinstance(rev_path, list) else int(rev_path)
                remote_pathogen = fwd_path_cnt + rev_path_cnt
            except Exception as e:
                print(f"Warning: Failed to parse remote JSON for candidate {rank}: {e}")

        # Fallback to alternatives.json if sequences are still empty
        fallback = alt_sequences.get(rank, {})
        full_alt = fallback.get("full_alt", {})
        
        if not fwd_seq or not rev_seq:
            if not fwd_seq:
                fwd_seq = fallback.get("fwd", "")
            if not rev_seq:
                rev_seq = fallback.get("rev", "")
            if not probe_seq:
                probe_seq = fallback.get("probe", "")

        if human_score == "-" and virus_score == "-" and remote_score == "-":
            continue

        is_pass = True
        if human_score != "-" and isinstance(human_score, (int, float)) and human_score < 70.0:
            is_pass = False
        if virus_score != "-" and isinstance(virus_score, (int, float)) and virus_score < 70.0:
            is_pass = False
        if remote_score != "-" and isinstance(remote_score, (int, float)) and remote_score < 70.0:
            is_pass = False
        
        status_text = "PASS" if is_pass else "FAIL"

        # Load metadata from alternatives.json
        p3_penalty = full_alt.get("p3_penalty", "-")
        vienna_penalty = full_alt.get("vienna_penalty", "-")
        final_score = full_alt.get("final_score", "-")
        tier = full_alt.get("tier", "-")
        
        primers_dict = full_alt.get("primers", {})
        fwd_info = primers_dict.get("forward", {})
        rev_info = primers_dict.get("reverse", {})
        probe_info = primers_dict.get("probe", {})
        
        fwd_tm = fwd_info.get("tm", "-")
        rev_tm = rev_info.get("tm", "-")
        probe_tm = probe_info.get("tm", "-") if probe_info else "-"
        
        fwd_gc = fwd_info.get("gc", "-")
        rev_gc = rev_info.get("gc", "-")
        probe_gc = probe_info.get("gc", "-") if probe_info else "-"
        
        qc_dict = full_alt.get("qc", {})
        cross_dimer_dg = qc_dict.get("cross_dimer_dg", "-")
        dimer_dg = qc_dict.get("dimer_dg", "-")
        hairpin_dg = qc_dict.get("hairpin_dg", "-")
        
        amplicon_dict = full_alt.get("amplicon", {})
        amplicon_length = amplicon_dict.get("length", "-")

        summary_data.append({
            "rank": rank,
            "fwd_seq": fwd_seq,
            "rev_seq": rev_seq,
            "probe_seq": probe_seq,
            "fwd_start": fwd_info.get("start"),
            "fwd_end": fwd_info.get("end"),
            "rev_start": rev_info.get("start"),
            "rev_end": rev_info.get("end"),
            "probe_start": probe_info.get("start") if probe_info else None,
            "probe_end": probe_info.get("end") if probe_info else None,
            "amplicon_start": amplicon_dict.get("start"),
            "amplicon_end": amplicon_dict.get("end"),
            "p3_penalty": p3_penalty,
            "vienna_penalty": vienna_penalty,
            "final_score": final_score,
            "tier": tier,
            "fwd_tm": fwd_tm,
            "rev_tm": rev_tm,
            "probe_tm": probe_tm,
            "fwd_gc": fwd_gc,
            "rev_gc": rev_gc,
            "probe_gc": probe_gc,
            "cross_dimer_dg": cross_dimer_dg,
            "dimer_dg": dimer_dg,
            "hairpin_dg": hairpin_dg,
            "amplicon_length": amplicon_length,
            "human_score": human_score,
            "human_grade": human_grade,
            "human_status": human_status,
            "human_offtargets": human_offtargets,
            "human_dangerous": human_dangerous,
            "virus_score": virus_score,
            "virus_grade": virus_grade,
            "virus_status": virus_status,
            "virus_pathogen": virus_pathogen,
            "virus_offtargets": virus_offtargets,
            "remote_score": remote_score,
            "remote_grade": remote_grade,
            "remote_status": remote_status,
            "remote_offtargets": remote_offtargets,
            "remote_dangerous": remote_dangerous,
            "remote_pathogen": remote_pathogen,
            "status": status_text
        })

    if not summary_data:
        return

    # Generate CSV
    csv_path = outdir / "blast_specificity_summary.csv"
    csv_headers = [
        "Candidate_Rank", "Forward_Sequence", "Reverse_Sequence", "Probe_Sequence",
        "Human_Score", "Human_Grade", "Human_Status", "Human_Total_Offtargets", "Human_Dangerous_3Prime_Offtargets",
        "Virus_Score", "Virus_Grade", "Virus_Status", "Virus_Pathogen_Hits", "Virus_Total_Offtargets",
        "Remote_Score", "Remote_Grade", "Remote_Status", "Remote_Total_Offtargets", "Remote_Dangerous_3Prime_Offtargets", "Remote_Pathogen_Hits",
        "Overall_Decision"
    ]
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(csv_headers)
        for d in summary_data:
            writer.writerow([
                f"Candidate_{d['rank']}", d["fwd_seq"], d["rev_seq"], d["probe_seq"],
                d["human_score"], d["human_grade"], d["human_status"], d["human_offtargets"], d["human_dangerous"],
                d["virus_score"], d["virus_grade"], d["virus_status"], d["virus_pathogen"], d["virus_offtargets"],
                d["remote_score"], d["remote_grade"], d["remote_status"], d["remote_offtargets"], d["remote_dangerous"], d["remote_pathogen"],
                d["status"]
            ])
    print(f"📊 Overall Summary CSV generated: {csv_path}")

    # Generate Markdown Report
    md_path = outdir / "blast_specificity_summary.md"
    md_lines = [
        "# 📋 M09 BLAST Specificity Verification Summary",
        "",
        f"*Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*",
        "",
        "Laporan ringkasan ini membandingkan spesifisitas seluruh kandidat primer hasil pencarian BLAST untuk mempermudah keputusan pemilihan primer selanjutnya.",
        "",
        "## Summary Table",
        ""
    ]

    is_local_run = any(d["human_score"] != "-" or d["virus_score"] != "-" for d in summary_data)
    
    if is_local_run:
        md_lines.append("| Candidate | Overall | Human Score (Grade) | Human 3' Off-targets | Virus Score (Grade) | Virus Pathogen Hits |")
        md_lines.append("|-----------|---------|---------------------|----------------------|---------------------|---------------------|")
        for d in summary_data:
            h_score_str = f"{d['human_score']} ({d['human_grade']})" if d['human_score'] != "-" else "-"
            v_score_str = f"{d['virus_score']} ({d['virus_grade']})" if d['virus_score'] != "-" else "-"
            status_icon = "✅ PASS" if d["status"] == "PASS" else "❌ FAIL"
            md_lines.append(
                f"| **Candidate {d['rank']}** | {status_icon} | {h_score_str} | {d['human_dangerous']} | {v_score_str} | {d['virus_pathogen']} |"
            )
    else:
        md_lines.append("| Candidate | Overall | Specificity Score (Grade) | Total Off-targets | Pathogen Hits |")
        md_lines.append("|-----------|---------|---------------------------|-------------------|---------------|")
        for d in summary_data:
            r_score_str = f"{d['remote_score']} ({d['remote_grade']})" if d['remote_score'] != "-" else "-"
            status_icon = "✅ PASS" if d["status"] == "PASS" else "❌ FAIL"
            md_lines.append(
                f"| **Candidate {d['rank']}** | {status_icon} | {r_score_str} | {d['remote_offtargets']} | {d['remote_pathogen']} |"
            )

    md_lines.append("")
    md_lines.append("## Detailed Sequences")
    md_lines.append("")
    md_lines.append("| Candidate | Forward Primer (5' → 3') | Reverse Primer (5' → 3') | Probe (5' → 3') |")
    md_lines.append("|-----------|--------------------------|--------------------------|-----------------|")
    for d in summary_data:
        probe_seq_str = f"`{d['probe_seq']}`" if d['probe_seq'] else "-"
        md_lines.append(f"| **Candidate {d['rank']}** | `{d['fwd_seq']}` | `{d['rev_seq']}` | {probe_seq_str} |")

    md_lines.append("")
    md_lines.append("---")
    md_lines.append("*Report generated by M09 Spec Summary Compiler*")

    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(md_lines))
    print(f"📝 Overall Summary Markdown report generated: {md_path}")

    # --- GENERATE DOWNSTREAM CHOSEN CANDIDATES (PASS) FILE AND GENOMIC TRACKS (GFF3 / BED) ---
    pass_candidates = [d for d in summary_data if d["status"] == "PASS"]
    if not pass_candidates:
        print("⚠️ No candidates passed the BLAST specificity threshold (Score >= 70). Skipping final reports.")
        return

    # 1. Save final_candidates.json
    json_path = outdir / "final_candidates.json"
    try:
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(pass_candidates, f, indent=2)
        print(f"💾 final_candidates.json generated: {json_path}")
    except Exception as e:
        print(f"Warning: Failed to generate final_candidates.json: {e}")

    # 2. Save final_candidates.gff3 (Genomic track for IGV/UCSC)
    target = target_id if target_id else "NC_002016.1"
    gff_path = outdir / "final_candidates.gff3"
    try:
        with open(gff_path, "w", encoding="utf-8") as f:
            f.write("##gff-version 3\n")
            f.write(f"##sequence-region {target} 1 15000\n")
            for d in pass_candidates:
                rank = d["rank"]
                amp_start = d["amplicon_start"]
                amp_end = d["amplicon_end"]
                # Write Amplicon
                if amp_start is not None and amp_end is not None:
                    # GFF3 is 1-based, inclusive coordinates
                    f.write(f"{target}\tPrimerLab\tamplicon\t{amp_start+1}\t{amp_end}\t.\t+\t.\tID=cand_{rank}_amp;Name=Candidate_{rank}_Amplicon;Score={d['final_score']};Tier={d['tier']}\n")
                
                # Write Forward
                fwd_start = d["fwd_start"]
                fwd_end = d["fwd_end"]
                if fwd_start is not None and fwd_end is not None:
                    f.write(f"{target}\tPrimerLab\tprimer\t{fwd_start+1}\t{fwd_end}\t.\t+\t.\tID=cand_{rank}_F;Name=Candidate_{rank}_Forward;Parent=cand_{rank}_amp;Tm={d['fwd_tm']};GC={d['fwd_gc']}\n")
                
                # Write Reverse
                rev_start = d["rev_start"]
                rev_end = d["rev_end"]
                if rev_start is not None and rev_end is not None:
                    f.write(f"{target}\tPrimerLab\tprimer\t{rev_start+1}\t{rev_end}\t.\t-\t.\tID=cand_{rank}_R;Name=Candidate_{rank}_Reverse;Parent=cand_{rank}_amp;Tm={d['rev_tm']};GC={d['rev_gc']}\n")
                
                # Write Probe
                prb_start = d["probe_start"]
                prb_end = d["probe_end"]
                if prb_start is not None and prb_end is not None:
                    f.write(f"{target}\tPrimerLab\tprobe\t{prb_start+1}\t{prb_end}\t.\t+\t.\tID=cand_{rank}_P;Name=Candidate_{rank}_Probe;Parent=cand_{rank}_amp;Tm={d['probe_tm']};GC={d['probe_gc']}\n")
        print(f"🧬 GFF3 Genomic track generated: {gff_path}")
    except Exception as e:
        print(f"Warning: Failed to generate GFF3: {e}")

    # 3. Save final_candidates.bed (Track for IGV/UCSC)
    bed_path = outdir / "final_candidates.bed"
    try:
        with open(bed_path, "w", encoding="utf-8") as f:
            f.write('track name="PrimerLab RAA Candidates" description="Verified RAA Primers and Probes" useScore=1\n')
            for d in pass_candidates:
                rank = d["rank"]
                # BED is 0-based start, 1-based exclusive end
                amp_start = d["amplicon_start"]
                amp_end = d["amplicon_end"]
                if amp_start is not None and amp_end is not None:
                    f.write(f"{target}\t{amp_start}\t{amp_end}\tCandidate_{rank}_Amplicon\t1000\t+\n")
                
                fwd_start = d["fwd_start"]
                fwd_end = d["fwd_end"]
                if fwd_start is not None and fwd_end is not None:
                    f.write(f"{target}\t{fwd_start}\t{fwd_end}\tCandidate_{rank}_Forward\t800\t+\n")
                
                rev_start = d["rev_start"]
                rev_end = d["rev_end"]
                if rev_start is not None and rev_end is not None:
                    f.write(f"{target}\t{rev_start}\t{rev_end}\tCandidate_{rank}_Reverse\t800\t-\n")
                
                prb_start = d["probe_start"]
                prb_end = d["probe_end"]
                if prb_start is not None and prb_end is not None:
                    f.write(f"{target}\t{prb_start}\t{prb_end}\tCandidate_{rank}_Probe\t900\t+\n")
        print(f"🧬 BED track generated: {bed_path}")
    except Exception as e:
        print(f"Warning: Failed to generate BED: {e}")

    # 4. Save final_candidates_report.md (Premium markdown report for human viewing)
    report_path = outdir / "final_candidates_report.md"
    try:
        rep_lines = [
            "# 🏆 Final Verified RAA Candidates Report",
            "",
            f"*Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*",
            "",
            "Dokumen ini memuat daftar kandidat primer dan probe yang **LOLOS verifikasi spesifisitas BLAST** (tidak memicu off-target berbahaya) dan memenuhi semua kriteria QC termodinamika.",
            "",
            "## Summary of Passed Candidates",
            "",
            "| Rank | Tier | Penalty Score | Vienna Dimer dG | Human Score | Virus Pathogen Hits | Amplicon Size | Decision |",
            "|------|------|---------------|-----------------|-------------|---------------------|---------------|----------|"
        ]
        for d in pass_candidates:
            v_dg = d["cross_dimer_dg"]
            v_dg_str = f"{v_dg} kcal/mol" if v_dg != "-" else "-"
            rep_lines.append(
                f"| Rank {d['rank']} | {d['tier']} | {d['final_score']} | {v_dg_str} | {d['human_score']} | {d['virus_pathogen']} | {d['amplicon_length']} bp | ✅ PASS |"
            )
        
        rep_lines.append("")
        rep_lines.append("## Detailed Coordinates & Sequences")
        rep_lines.append("")
        for d in pass_candidates:
            rep_lines.extend([
                f"### 📍 Candidate Rank {d['rank']} (Penalty: {d['final_score']})",
                "",
                f"* **Tier**: {d['tier']}",
                f"* **Amplicon Size**: {d['amplicon_length']} bp (Position: {d['amplicon_start']}-{d['amplicon_end']})",
                f"* **Thermodynamics**: Dimer dG: {d['dimer_dg']} | Hairpin dG: {d['hairpin_dg']} | Cross-dimer dG: {d['cross_dimer_dg']}",
                f"* **BLAST Specificity**: Human Score: {d['human_score']} (Off-targets: {d['human_offtargets']}) | Virus Pathogen Hits: {d['virus_pathogen']}",
                "",
                "**Sequences:**",
                f"* **Forward Primer** ({d['fwd_start']}-{d['fwd_end']}): `{d['fwd_seq']}` (Tm: {d['fwd_tm']}°C, GC: {d['fwd_gc']}%)",
                f"* **Reverse Primer** ({d['rev_start']}-{d['rev_end']}): `{d['rev_seq']}` (Tm: {d['rev_tm']}°C, GC: {d['rev_gc']}%)",
            ])
            if d['probe_seq'] and d['probe_seq'] != "-":
                rep_lines.append(f"* **Probe** ({d['probe_start']}-{d['probe_end']}): `{d['probe_seq']}` (Tm: {d['probe_tm']}°C, GC: {d['probe_gc']}%)")
            rep_lines.append("")
            
        rep_lines.extend([
            "---",
            "### 🧬 IGV Integration & Visualization",
            "Kandidat di atas juga diekspor ke dalam format standar anotasi genomik yang bisa langsung dimasukkan ke genome browser seperti **IGV** atau **UCSC**:",
            f"* **GFF3 Track**: [final_candidates.gff3](file:///{str(gff_path.absolute()).replace(chr(92), '/')})",
            f"* **BED Track**: [final_candidates.bed](file:///{str(bed_path.absolute()).replace(chr(92), '/')})",
            "",
            "Untuk melakukan visualisasi publikasi jurnal atau di IGV:",
            "1. Buka **IGV**, load genome/consensus target Anda.",
            "2. Klik `File -> Load from File...` lalu pilih file `final_candidates.gff3` atau `final_candidates.bed`.",
            "3. Primer dan Probe akan dirender secara visual sesuai posisi binding situsnya di atas genom target."
        ])
        
        with open(report_path, "w", encoding="utf-8") as f:
            f.write("\n".join(rep_lines))
        print(f"🏆 Premium Final Candidates Report generated: {report_path}")
    except Exception as e:
        print(f"Warning: Failed to generate premium report: {e}")

if __name__ == "__main__":
    main()
