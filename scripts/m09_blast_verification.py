import os
import sys
import json
import argparse
import subprocess
import tempfile
from pathlib import Path

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

if __name__ == "__main__":
    main()
