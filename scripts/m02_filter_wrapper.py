import subprocess
import os
import sys

subsampling = snakemake.params.subsampling
input_seq = snakemake.input.sequences
input_meta = snakemake.input.metadata
output_seq = snakemake.output.sequences
log_file = snakemake.log[0]

def build_cmd(query=None, max_seq=None, seq_per_group=None, group_by=None, out_fasta=output_seq):
    cmd = [
        "augur", "filter",
        "--sequences", input_seq,
        "--metadata", input_meta,
        "--min-length", str(subsampling.get("min_length", 800)),
        "--subsample-seed", str(subsampling.get("seed", 42)),
        "--exclude-ambiguous-dates-by", "any",
        "--output-sequences", out_fasta
    ]
    
    # Optional parameters based on configuration
    if group_by:
        # Augur group-by accepts space-separated columns, but from config it might be a single string
        # e.g., "subtype country year month". In subprocess we pass it as separate args.
        cmd.extend(["--group-by"])
        cmd.extend(group_by.split())
    
    if query:
        cmd.extend(["--query", query])
        
    if max_seq:
        cmd.extend(["--subsample-max-sequences", str(max_seq)])
    elif seq_per_group:
        cmd.extend(["--sequences-per-group", str(seq_per_group)])
        
    return cmd

with open(log_file, "w") as f_log:
    focal_config = subsampling.get("focal", {})
    if focal_config.get("enabled", False):
        print("Running Focal + Contextual Sampling...", file=f_log)
        contextual_config = subsampling.get("contextual", {})
        
        # 1. Focal
        focal_out = output_seq + ".focal.fasta"
        cmd_focal = build_cmd(
            query=focal_config.get("query"),
            max_seq=focal_config.get("max_sequences"),
            seq_per_group=focal_config.get("sequences_per_group"),
            group_by=focal_config.get("group_by", subsampling.get("group_by")),
            out_fasta=focal_out
        )
        print("FOCAL CMD:", " ".join(cmd_focal), file=f_log)
        sys.stdout.flush()
        try:
            subprocess.run(cmd_focal, stdout=f_log, stderr=subprocess.STDOUT, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Focal filter failed with error code {e.returncode}", file=f_log)
            sys.exit(1)
            
        # 2. Contextual
        ctx_out = output_seq + ".ctx.fasta"
        cmd_ctx = build_cmd(
            query=contextual_config.get("query"),
            max_seq=contextual_config.get("max_sequences"),
            seq_per_group=contextual_config.get("sequences_per_group"),
            group_by=contextual_config.get("group_by", subsampling.get("group_by")),
            out_fasta=ctx_out
        )
        print("\nCONTEXTUAL CMD:", " ".join(cmd_ctx), file=f_log)
        sys.stdout.flush()
        try:
            subprocess.run(cmd_ctx, stdout=f_log, stderr=subprocess.STDOUT, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Contextual filter failed with error code {e.returncode}", file=f_log)
            sys.exit(1)
        
        # 3. Combine deduplicated
        print("\nCombining fasta files...", file=f_log)
        seen = set()
        with open(output_seq, "w") as out_f:
            for fname in [focal_out, ctx_out]:
                if not os.path.exists(fname): 
                    print(f"Warning: {fname} does not exist. Skipping.", file=f_log)
                    continue
                with open(fname, "r") as in_f:
                    write_mode = False
                    for line in in_f:
                        if line.startswith(">"):
                            header = line.strip()
                            if header not in seen:
                                seen.add(header)
                                write_mode = True
                                out_f.write(line)
                            else:
                                write_mode = False
                        elif write_mode:
                            out_f.write(line)
                            
        # Cleanup
        if os.path.exists(focal_out): os.remove(focal_out)
        if os.path.exists(ctx_out): os.remove(ctx_out)
        print(f"Done. Unique sequences retained: {len(seen)}", file=f_log)
        
    else:
        print("Running Standard Global Subsampling...", file=f_log)
        cmd = build_cmd(
            query=subsampling.get("query"),
            max_seq=subsampling.get("max_sequences"),
            seq_per_group=subsampling.get("sequences_per_group"),
            group_by=subsampling.get("group_by")
        )
        print("CMD:", " ".join(cmd), file=f_log)
        sys.stdout.flush()
        try:
            subprocess.run(cmd, stdout=f_log, stderr=subprocess.STDOUT, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Global filter failed with error code {e.returncode}", file=f_log)
            sys.exit(1)
