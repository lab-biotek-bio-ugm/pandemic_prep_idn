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
    parser.add_argument(
        "--compile-only", action="store_true",
        help="Skip BLAST execution; only compile summary reports from existing results."
    )
    parser.add_argument(
        "--bridge-report", default=None,
        help="Path to M07 bridge_report.json for GARD breakpoint cross-check."
    )
    parser.add_argument(
        "--ref-fasta", default=None,
        help="Path to reference genome FASTA file."
    )
    
    args = parser.parse_args()

    # --- Compile-only mode: skip BLAST, just compile reports ---
    if args.compile_only:
        print("📊 Compile-only mode: generating summary reports from existing results...")
        try:
            compile_summary(
                str(args.outdir), args.alternatives, args.target,
                bridge_report_path=args.bridge_report,
                ref_fasta_path=args.ref_fasta,
            )
        except Exception as e:
            import traceback
            print(f"⚠️ Warning: Failed to compile specificity summary: {e}")
            traceback.print_exc()
        return
    
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
    
    skipped_count = 0
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

        # --- GAP-1: Resume/skip for already-BLASTed candidates ---
        cand_out_dir = outdir / f"candidate_{rank}"
        if args.db_label:
            result_check_dir = cand_out_dir / f"db_{args.db_label}"
        else:
            result_check_dir = cand_out_dir
        
        result_check_file = result_check_dir / "blast_result.json"
        if result_check_file.exists():
            print(f"⏭️ Candidate #{rank}{db_label_str} already verified. Skipping (result exists at {result_check_file}).")
            skipped_count += 1
            continue
            
        # Create a temp file containing the primers dict for this candidate
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as temp_f:
            json.dump(primers, temp_f)
            temp_json_path = temp_f.name
            
        try:
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

    if skipped_count > 0:
        print(f"\n📌 {skipped_count} candidate(s) skipped (already verified from a previous run).")
                
    # Write marker file only if this is a standalone run (no db_label / remote mode).
    # In local dual-db mode, the Snakemake shell writes the marker after both db runs.
    if not args.db_label:
        marker_path = outdir / "verification_complete.txt"
        with open(marker_path, "w") as f:
            f.write(f"Completed verification for {len(to_process)} candidates at {Path(args.alternatives).name}\n")
        print(f"\n🎉 Verification complete. Marker file written: {marker_path}")
    else:
        print(f"\n✅ DB run [{args.db_label}] complete. Marker will be written by Snakemake after all DB runs.")

    # --- BUG-1 FIX: Only compile summary if this is a standalone run (no db_label).
    # In local dual-db mode, compile_summary is called ONCE by Snakemake via --compile-only
    # AFTER both human and virus runs finish. This prevents partial/duplicate compilation.
    if not args.db_label:
        try:
            compile_summary(
                str(outdir), args.alternatives, args.target,
                bridge_report_path=args.bridge_report,
                ref_fasta_path=args.ref_fasta,
            )
        except Exception as e:
            import traceback
            print(f"⚠️ Warning: Failed to compile specificity summary: {e}")
            traceback.print_exc()


def get_grade(score):
    """Assign letter grade based on BLAST specificity score.
    Higher score = better specificity (0-100 scale).
    """
    if score >= 90: return "A"
    elif score >= 80: return "B"
    elif score >= 70: return "C"
    elif score >= 65: return "C-"   # POLISH-1: BORDERLINE zone
    elif score >= 60: return "D"
    else: return "F"


def _safe_coord(val):
    """BUG-3: Validate and return coordinate, or None if invalid."""
    if val is None:
        return None
    try:
        v = int(val)
        return v if v >= 0 else None
    except (ValueError, TypeError):
        return None


def _clean_gff_val(val):
    """Clean value for GFF3 attribute writing by removing unit suffixes
    and percent-encoding any reserved characters.
    """
    if val is None or val == "-":
        return "."
    s = str(val).strip()
    # Remove common unit suffixes that cause URLDecoder issues in IGV
    s = s.replace("°C", "").replace("\u00b0C", "").replace("%", "").strip()
    import urllib.parse
    return urllib.parse.quote(s)


def _load_gard_breakpoints(bridge_report_path: Optional[str]) -> list:
    """POLISH-4: Load GARD breakpoints from bridge_report.json."""
    if not bridge_report_path or not os.path.exists(bridge_report_path):
        return []
    try:
        with open(bridge_report_path, "r") as f:
            data = json.load(f)
        return data.get("gard_breakpoints", [])
    except Exception:
        return []


def _check_gard_overlap(start, end, gard_breakpoints, buffer=50):
    """POLISH-4: Check if a coordinate range overlaps with any GARD breakpoint."""
    if not gard_breakpoints or start is None or end is None:
        return False, []
    overlapping = []
    for bp in gard_breakpoints:
        gp = bp.get("genome_pos", -1)
        if start - buffer <= gp <= end + buffer:
            overlapping.append(bp)
    return len(overlapping) > 0, overlapping


def compile_summary(outdir_path: str, alternatives_path: Optional[str] = None,
                    target_id: Optional[str] = None,
                    bridge_report_path: Optional[str] = None,
                    ref_fasta_path: Optional[str] = None):
    """
    Scans output directory, loads JSON results for all candidates,
    and generates master CSV and Markdown reports for comparison.
    Also produces dedicated BED, GFF3, JSON, IGV session and MD reports specifically for PASSing candidates.
    """
    outdir = Path(outdir_path)
    if not outdir.exists():
        return

    # POLISH-4: Load GARD breakpoints for cross-check
    gard_breakpoints = _load_gard_breakpoints(bridge_report_path)
    if gard_breakpoints:
        print(f"📍 Loaded {len(gard_breakpoints)} GARD breakpoints for cross-check.")

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

        # --- POLISH-1: BORDERLINE status for scores 65-69.9 ---
        is_pass = True
        is_borderline = False
        for score_val in [human_score, virus_score, remote_score]:
            if score_val != "-" and isinstance(score_val, (int, float)):
                if score_val < 65.0:
                    is_pass = False
                elif score_val < 70.0:
                    is_borderline = True

        if not is_pass:
            status_text = "FAIL"
        elif is_borderline:
            status_text = "BORDERLINE"
        else:
            status_text = "PASS"

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

        # --- POLISH-5: Tm_Delta computation ---
        tm_delta = "-"
        if isinstance(fwd_tm, (int, float)) and isinstance(rev_tm, (int, float)):
            tm_delta = round(abs(fwd_tm - rev_tm), 2)

        # --- BUG-3: Validate coordinates ---
        fwd_start_raw = fwd_info.get("start")
        fwd_end_raw = fwd_info.get("end")
        rev_start_raw = rev_info.get("start")
        rev_end_raw = rev_info.get("end")
        probe_start_raw = probe_info.get("start") if probe_info else None
        probe_end_raw = probe_info.get("end") if probe_info else None
        amp_start_raw = amplicon_dict.get("start")
        amp_end_raw = amplicon_dict.get("end")

        fwd_start_v = _safe_coord(fwd_start_raw)
        fwd_end_v = _safe_coord(fwd_end_raw)
        rev_start_v = _safe_coord(rev_start_raw)
        rev_end_v = _safe_coord(rev_end_raw)
        probe_start_v = _safe_coord(probe_start_raw)
        probe_end_v = _safe_coord(probe_end_raw)
        amp_start_v = _safe_coord(amp_start_raw)
        amp_end_v = _safe_coord(amp_end_raw)

        # --- POLISH-4: GARD breakpoint cross-check ---
        gard_near = False
        gard_details = []
        if amp_start_v is not None and amp_end_v is not None:
            gard_near, gard_details = _check_gard_overlap(amp_start_v, amp_end_v, gard_breakpoints)

        summary_data.append({
            "rank": rank,
            "fwd_seq": fwd_seq,
            "rev_seq": rev_seq,
            "probe_seq": probe_seq,
            "fwd_start": fwd_start_v,
            "fwd_end": fwd_end_v,
            "rev_start": rev_start_v,
            "rev_end": rev_end_v,
            "probe_start": probe_start_v,
            "probe_end": probe_end_v,
            "amplicon_start": amp_start_v,
            "amplicon_end": amp_end_v,
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
            "tm_delta": tm_delta,
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
            "status": status_text,
            "gard_near": gard_near,
            "gard_details": gard_details,
        })

    if not summary_data:
        return

    # Generate CSV — POLISH-5: includes Tm_Delta column
    csv_path = outdir / "blast_specificity_summary.csv"
    csv_headers = [
        "Candidate_Rank", "Forward_Sequence", "Reverse_Sequence", "Probe_Sequence",
        "Fwd_Tm", "Rev_Tm", "Tm_Delta", "Fwd_GC", "Rev_GC", "Probe_Tm", "Probe_GC",
        "PrimerLab_Penalty_Score", "Vienna_Penalty", "Tier",
        "Amplicon_Length", "Cross_Dimer_dG", "Dimer_dG", "Hairpin_dG",
        "Human_Specificity_Score", "Human_Grade", "Human_Status", "Human_Total_Offtargets", "Human_Dangerous_3Prime_Offtargets",
        "Virus_Specificity_Score", "Virus_Grade", "Virus_Status", "Virus_Pathogen_Hits", "Virus_Total_Offtargets",
        "Remote_Specificity_Score", "Remote_Grade", "Remote_Status", "Remote_Total_Offtargets", "Remote_Dangerous_3Prime_Offtargets", "Remote_Pathogen_Hits",
        "GARD_Breakpoint_Nearby",
        "Overall_Decision"
    ]
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(csv_headers)
        for d in summary_data:
            writer.writerow([
                f"Candidate_{d['rank']}", d["fwd_seq"], d["rev_seq"], d["probe_seq"],
                d["fwd_tm"], d["rev_tm"], d["tm_delta"], d["fwd_gc"], d["rev_gc"], d["probe_tm"], d["probe_gc"],
                d["final_score"], d["vienna_penalty"], d["tier"],
                d["amplicon_length"], d["cross_dimer_dg"], d["dimer_dg"], d["hairpin_dg"],
                d["human_score"], d["human_grade"], d["human_status"], d["human_offtargets"], d["human_dangerous"],
                d["virus_score"], d["virus_grade"], d["virus_status"], d["virus_pathogen"], d["virus_offtargets"],
                d["remote_score"], d["remote_grade"], d["remote_status"], d["remote_offtargets"], d["remote_dangerous"], d["remote_pathogen"],
                "YES" if d["gard_near"] else "NO",
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
        # --- BUG-2: Score unit legend ---
        "> **Scoring Legend:**",
        "> - **BLAST Specificity Score** (0-100): Makin **tinggi** = makin **spesifik** (baik). ✅",
        "> - **PrimerLab Penalty Score** (0+): Makin **rendah** = makin **baik** (sedikit penalti). ✅",
        "> - **PASS** = semua skor BLAST ≥ 70 | **BORDERLINE** = ada skor 65-69.9 | **FAIL** = ada skor < 65",
        "",
        "## Summary Table",
        ""
    ]

    is_local_run = any(d["human_score"] != "-" or d["virus_score"] != "-" for d in summary_data)
    
    if is_local_run:
        md_lines.append("| Candidate | Overall | PrimerLab Penalty | Tm Δ | Human Score (Grade) | Human 3' Off-targets | Virus Score (Grade) | Virus Pathogen Hits | GARD |")
        md_lines.append("|-----------|---------|-------------------|------|---------------------|----------------------|---------------------|---------------------|------|")
        for d in summary_data:
            h_score_str = f"{d['human_score']} ({d['human_grade']})" if d['human_score'] != "-" else "-"
            v_score_str = f"{d['virus_score']} ({d['virus_grade']})" if d['virus_score'] != "-" else "-"
            if d["status"] == "PASS":
                status_icon = "✅ PASS"
            elif d["status"] == "BORDERLINE":
                status_icon = "⚠️ BORDERLINE"
            else:
                status_icon = "❌ FAIL"
            gard_flag = "⚠️ YES" if d["gard_near"] else "—"
            md_lines.append(
                f"| **Candidate {d['rank']}** | {status_icon} | {d['final_score']} | {d['tm_delta']} | {h_score_str} | {d['human_dangerous']} | {v_score_str} | {d['virus_pathogen']} | {gard_flag} |"
            )
    else:
        md_lines.append("| Candidate | Overall | PrimerLab Penalty | Tm Δ | Specificity Score (Grade) | Total Off-targets | Pathogen Hits | GARD |")
        md_lines.append("|-----------|---------|-------------------|------|---------------------------|-------------------|---------------|------|")
        for d in summary_data:
            r_score_str = f"{d['remote_score']} ({d['remote_grade']})" if d['remote_score'] != "-" else "-"
            if d["status"] == "PASS":
                status_icon = "✅ PASS"
            elif d["status"] == "BORDERLINE":
                status_icon = "⚠️ BORDERLINE"
            else:
                status_icon = "❌ FAIL"
            gard_flag = "⚠️ YES" if d["gard_near"] else "—"
            md_lines.append(
                f"| **Candidate {d['rank']}** | {status_icon} | {d['final_score']} | {d['tm_delta']} | {r_score_str} | {d['remote_offtargets']} | {d['remote_pathogen']} | {gard_flag} |"
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

    # --- GENERATE DOWNSTREAM CHOSEN CANDIDATES (PASS + BORDERLINE) FILE AND GENOMIC TRACKS ---
    pass_candidates = [d for d in summary_data if d["status"] in ("PASS", "BORDERLINE")]
    if not pass_candidates:
        print("⚠️ No candidates passed the BLAST specificity threshold (Score >= 65). Skipping final reports.")
        return

    # 1. Save final_candidates.json
    json_path = outdir / "final_candidates.json"
    try:
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(pass_candidates, f, indent=2, default=str)
        print(f"💾 final_candidates.json generated: {json_path}")
    except Exception as e:
        print(f"Warning: Failed to generate final_candidates.json: {e}")

    # 2. Save final_candidates.gff3 — BUG-3: coordinate validation applied
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
                # Compute lengths of primer elements
                fwd_len = d["fwd_end"] - d["fwd_start"] + 1 if (d["fwd_start"] is not None and d["fwd_end"] is not None) else 0
                rev_len = d["rev_end"] - d["rev_start"] + 1 if (d["rev_start"] is not None and d["rev_end"] is not None) else 0
                prb_len = d["probe_end"] - d["probe_start"] + 1 if (d["probe_start"] is not None and d["probe_end"] is not None) else 0
                amp_len = amp_end - amp_start + 1 if (amp_start is not None and amp_end is not None) else 0

                # BUG-3: Only write if coordinates are valid (>= 0 and not None)
                if amp_start is not None and amp_end is not None and amp_start < amp_end:
                    # GFF3 is 1-based, inclusive coordinates
                    clean_score = _clean_gff_val(d['final_score'])
                    clean_tier = _clean_gff_val(d['tier'])
                    f.write(f"{target}\tPrimerLab\tmRNA\t{amp_start+1}\t{amp_end+1}\t{clean_score}\t+\t.\tID=cand_{rank}_amp;Name=Candidate_{rank}_({amp_len}_bp);Score={clean_score};Tier={clean_tier}\n")
                
                fwd_start = d["fwd_start"]
                fwd_end = d["fwd_end"]
                if fwd_start is not None and fwd_end is not None and fwd_start < fwd_end:
                    clean_tm = _clean_gff_val(d['fwd_tm'])
                    clean_gc = _clean_gff_val(d['fwd_gc'])
                    f.write(f"{target}\tPrimerLab\texon\t{fwd_start+1}\t{fwd_end+1}\t.\t+\t.\tID=cand_{rank}_F;Name=Fwd_({fwd_len}_bp);Parent=cand_{rank}_amp;Tm={clean_tm};GC={clean_gc}\n")
                
                rev_start = d["rev_start"]
                rev_end = d["rev_end"]
                if rev_start is not None and rev_end is not None and rev_start < rev_end:
                    clean_tm = _clean_gff_val(d['rev_tm'])
                    clean_gc = _clean_gff_val(d['rev_gc'])
                    f.write(f"{target}\tPrimerLab\texon\t{rev_start+1}\t{rev_end+1}\t.\t-\t.\tID=cand_{rank}_R;Name=Rev_({rev_len}_bp);Parent=cand_{rank}_amp;Tm={clean_tm};GC={clean_gc}\n")
                
                prb_start = d["probe_start"]
                prb_end = d["probe_end"]
                if prb_start is not None and prb_end is not None and prb_start < prb_end:
                    clean_tm = _clean_gff_val(d['probe_tm'])
                    clean_gc = _clean_gff_val(d['probe_gc'])
                    f.write(f"{target}\tPrimerLab\texon\t{prb_start+1}\t{prb_end+1}\t.\t+\t.\tID=cand_{rank}_P;Name=Prb_({prb_len}_bp);Parent=cand_{rank}_amp;Tm={clean_tm};GC={clean_gc}\n")
        print(f"🧬 GFF3 Genomic track generated: {gff_path}")
    except Exception as e:
        print(f"Warning: Failed to generate GFF3: {e}")

    # 3. Save final_candidates.bed — BUG-3: coordinate validation applied
    bed_path = outdir / "final_candidates.bed"
    try:
        with open(bed_path, "w", encoding="utf-8") as f:
            f.write('track name="PrimerLab RAA Candidates" description="Verified RAA Primers and Probes" useScore=1\n')
            for d in pass_candidates:
                rank = d["rank"]
                amp_start = d["amplicon_start"]
                amp_end = d["amplicon_end"]
                if amp_start is not None and amp_end is not None and amp_start < amp_end:
                    f.write(f"{target}\t{amp_start}\t{amp_end+1}\tCandidate_{rank}_Amplicon\t1000\t+\n")
                
                fwd_start = d["fwd_start"]
                fwd_end = d["fwd_end"]
                if fwd_start is not None and fwd_end is not None and fwd_start < fwd_end:
                    f.write(f"{target}\t{fwd_start}\t{fwd_end+1}\tCandidate_{rank}_Forward\t800\t+\n")
                
                rev_start = d["rev_start"]
                rev_end = d["rev_end"]
                if rev_start is not None and rev_end is not None and rev_start < rev_end:
                    f.write(f"{target}\t{rev_start}\t{rev_end+1}\tCandidate_{rank}_Reverse\t800\t-\n")
                
                prb_start = d["probe_start"]
                prb_end = d["probe_end"]
                if prb_start is not None and prb_end is not None and prb_start < prb_end:
                    f.write(f"{target}\t{prb_start}\t{prb_end+1}\tCandidate_{rank}_Probe\t900\t+\n")
        print(f"🧬 BED track generated: {bed_path}")
    except Exception as e:
        print(f"Warning: Failed to generate BED: {e}")

    # --- GAP-4: Generate IGV session XML ---
    igv_session_path = outdir / "igv_session.xml"
    try:
        _generate_igv_session(igv_session_path, gff_path, bed_path, target, outdir, ref_fasta_path=ref_fasta_path)
        print(f"🖥️ IGV session file generated: {igv_session_path}")
    except Exception as e:
        print(f"Warning: Failed to generate IGV session: {e}")

    # 4. Save final_candidates_report.md — POLISH-2: Table format per candidate
    report_path = outdir / "final_candidates_report.md"
    try:
        rep_lines = [
            "# 🏆 Final Verified RAA Candidates Report",
            "",
            f"*Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*",
            "",
            "Dokumen ini memuat daftar kandidat primer dan probe yang **LOLOS verifikasi spesifisitas BLAST** dan memenuhi semua kriteria QC termodinamika.",
            "",
            # BUG-2: Score legend
            "> **Scoring Guide:**",
            "> - **PrimerLab Penalty Score** — Makin **rendah** = makin **baik** (0 = sempurna, sedikit penalti).",
            "> - **BLAST Specificity Score** — Makin **tinggi** = makin **spesifik** (100 = tidak ada off-target).",
            "> - **Tm Delta** — Selisih Tm antara forward dan reverse (ideal < 2°C).",
            "",
            "## Summary of Verified Candidates",
            "",
            "| Rank | Status | Tier | Penalty Score ↓ | Tm Δ (°C) | Human Score ↑ | Virus Pathogen | Amplicon | GARD |",
            "|------|--------|------|-----------------|-----------|---------------|----------------|----------|------|"
        ]
        for d in pass_candidates:
            if d["status"] == "PASS":
                status_icon = "✅ PASS"
            else:
                status_icon = "⚠️ BORDER"
            gard_flag = "⚠️" if d["gard_near"] else "—"
            rep_lines.append(
                f"| {d['rank']} | {status_icon} | {d['tier']} | {d['final_score']} | {d['tm_delta']} | {d['human_score']} | {d['virus_pathogen']} | {d['amplicon_length']} bp | {gard_flag} |"
            )
        
        rep_lines.append("")
        rep_lines.append("## Detailed Candidate Profiles")
        rep_lines.append("")
        for d in pass_candidates:
            # --- POLISH-2: Compact table per candidate ---
            rep_lines.extend([
                f"### 📍 Candidate Rank {d['rank']}",
                "",
                "| Property | Value |",
                "|----------|-------|",
                f"| **Status** | {'✅ PASS' if d['status'] == 'PASS' else '⚠️ BORDERLINE'} |",
                f"| **Tier** | {d['tier']} |",
                f"| **PrimerLab Penalty Score** ↓ | {d['final_score']} |",
                f"| **Vienna Penalty** | {d['vienna_penalty']} |",
                f"| **Amplicon Size** | {d['amplicon_length']} bp ({d['amplicon_start']}–{d['amplicon_end']}) |",
                f"| **Tm Delta (Fwd-Rev)** | {d['tm_delta']}°C |",
                f"| **Dimer dG** | {d['dimer_dg']} kcal/mol |",
                f"| **Hairpin dG** | {d['hairpin_dg']} kcal/mol |",
                f"| **Cross-dimer dG** | {d['cross_dimer_dg']} kcal/mol |",
                f"| **Human BLAST Score** ↑ | {d['human_score']} ({d['human_grade']}) |",
                f"| **Human Off-targets (3' involved)** | {d['human_dangerous']} |",
                f"| **Virus BLAST Score** ↑ | {d['virus_score']} ({d['virus_grade']}) |",
                f"| **Virus Pathogen Hits** | {d['virus_pathogen']} |",
            ])
            # POLISH-4: GARD warning
            if d["gard_near"]:
                rep_lines.append(f"| **⚠️ GARD Breakpoint** | Recombination breakpoint detected nearby: {[bp.get('genome_pos') for bp in d['gard_details']]} |")
            else:
                rep_lines.append(f"| **GARD Breakpoint** | — (none nearby) |")
            rep_lines.extend([
                "",
                "**Primer Sequences:**",
                "",
                "| Component | Sequence (5'→3') | Position | Tm (°C) | GC (%) |",
                "|-----------|------------------|----------|---------|--------|",
                f"| Forward | `{d['fwd_seq']}` | {d['fwd_start']}–{d['fwd_end']} | {d['fwd_tm']} | {d['fwd_gc']} |",
                f"| Reverse | `{d['rev_seq']}` | {d['rev_start']}–{d['rev_end']} | {d['rev_tm']} | {d['rev_gc']} |",
            ])
            if d['probe_seq'] and d['probe_seq'] != "-":
                rep_lines.append(f"| Probe | `{d['probe_seq']}` | {d['probe_start']}–{d['probe_end']} | {d['probe_tm']} | {d['probe_gc']} |")
            rep_lines.append("")
            
        # --- POLISH-6: Use relative paths ---
        rep_lines.extend([
            "---",
            "### 🧬 IGV Integration & Visualization",
            "Kandidat di atas juga diekspor ke dalam format standar anotasi genomik yang bisa langsung dimasukkan ke genome browser seperti **IGV** atau **UCSC**:",
            "",
            "| Format | File | Deskripsi |",
            "|--------|------|-----------|",
            "| GFF3 | `final_candidates.gff3` | Anotasi primer/probe/amplicon (1-based) |",
            "| BED | `final_candidates.bed` | Track BED untuk IGV/UCSC (0-based) |",
            "| IGV Session | `igv_session.xml` | Session file siap buka di IGV |",
            "| JSON | `final_candidates.json` | Data mentah untuk analisis lanjutan |",
            "",
            "**Cara memuat di IGV:**",
            "1. Buka **IGV**, load genome/consensus target Anda.",
            "2. Klik `File → Load from File...` lalu pilih `final_candidates.gff3` atau `final_candidates.bed`.",
            "3. **Atau** gunakan `File → Open Session...` untuk membuka `igv_session.xml` yang sudah menyertakan semua track.",
            "4. Primer dan Probe akan dirender secara visual sesuai posisi binding situsnya di atas genom target."
        ])
        
        with open(report_path, "w", encoding="utf-8") as f:
            f.write("\n".join(rep_lines))
        print(f"🏆 Premium Final Candidates Report generated: {report_path}")
    except Exception as e:
        print(f"Warning: Failed to generate premium report: {e}")


def _generate_igv_session(session_path: Path, gff_path: Path, bed_path: Path,
                          reference_id: str, outdir: Path, ref_fasta_path: Optional[str] = None):
    """GAP-4: Generate an IGV session XML file that pre-loads all genomic tracks."""
    # Use relative paths from the session file location
    gff_rel = gff_path.name
    
    # Check if M07 BED is available (one directory up, in m07 folder)
    m07_bed = outdir.parent / "m07" / "conserved_regions.bed"
    m07_bed_str = ""
    if m07_bed.exists():
        m07_bed_str = f'''
        <Resource path="../m07/conserved_regions.bed"/>'''

    # Determine genome value
    genome_val = reference_id if reference_id else "NC_002016.1"
    if ref_fasta_path and os.path.exists(ref_fasta_path):
        try:
            # Compute relative path from the session file directory (outdir)
            rel_ref = os.path.relpath(ref_fasta_path, start=outdir).replace("\\", "/")
            genome_val = rel_ref
        except Exception:
            pass

    xml_content = f'''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="{genome_val}" hasGeneTrack="false" version="8">
    <Resources>
        <Resource path="{gff_rel}"/>{m07_bed_str}
    </Resources>
    <Panel name="DataPanel" height="300">
        <Track id="{gff_rel}" name="Primer Candidates" renderer="BASIC_FEATURE" visible="true" displayMode="EXPANDED">
            <RenderOptions colorByAttribute="Name"/>
        </Track>
    </Panel>
</Session>
'''
    with open(session_path, "w", encoding="utf-8") as f:
        f.write(xml_content.strip())


if __name__ == "__main__":
    main()
