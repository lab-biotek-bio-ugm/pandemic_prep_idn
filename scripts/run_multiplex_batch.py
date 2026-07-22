#!/usr/bin/env python3
"""
Master Multi-Target Multiplex Orchestrator
Batch-executes Snakemake pipelines across multiple virus targets/configs,
and triggers Module 11 Multiplex Compatibility Screening across all target pairs.
"""

import os
import sys
import yaml
import argparse
import subprocess

def run_cmd(cmd):
    print(f"\n[Batch Orchestrator] Running: {cmd}")
    res = subprocess.run(cmd, shell=True)
    if res.returncode != 0:
        print(f"[Error] Command failed with code {res.returncode}: {cmd}")
        sys.exit(res.returncode)

def parse_build_name(config_path):
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)
    return cfg.get("build_name", os.path.basename(config_path).replace(".yaml", ""))

def main():
    parser = argparse.ArgumentParser(description="Multi-Target Batching & Multiplex Orchestrator")
    parser.add_argument("--configs", nargs="+", required=True, help="List of target config files (e.g., config_influenza_a_mp.yaml config_influenza_b_mp.yaml)")
    parser.add_argument("--cores", type=int, default=64, help="Number of CPU cores for Snakemake execution")
    parser.add_argument("--skip-snakemake", action="store_true", help="Skip Snakemake runs if result folders already exist")
    parser.add_argument("--snakemake-args", type=str, default="--use-conda --rerun-incomplete", help="Additional flags to pass directly to Snakemake (e.g., '--use-conda --unlock --rerun-incomplete')")
    parser.add_argument("--num-candidates", type=int, default=10, help="Number of top candidates per target to screen for multiplexing")
    args = parser.parse_args()

    build_names = []
    
    print("=" * 80)
    print("🚀 PANDEMIC PREPAREDNESS MULTI-TARGET BATCH ORCHESTRATOR")
    print("=" * 80)

    # Step 1: Run Snakemake for each target config
    for idx, cfg in enumerate(args.configs, 1):
        build_name = parse_build_name(cfg)
        build_names.append((cfg, build_name))
        print(f"\n[Target {idx}/{len(args.configs)}] Config: {cfg} ---> Build: {build_name}")
        
        if args.skip_snakemake:
            print(f"[Info] --skip-snakemake flag set. Skipping Snakemake execution for {build_name}.")
        else:
            extra_args = args.snakemake_args if args.snakemake_args else ""
            snakemake_cmd = f"snakemake -c {args.cores} --configfile {cfg} {extra_args}".strip()
            run_cmd(snakemake_cmd)

    # Step 2: Execute Module 11 Multiplex Compatibility for all unique target pairs
    print("\n" + "=" * 80)
    print("🧬 STARTING MODULE 11 MULTIPLEX COMPATIBILITY SCREENING ACROSS TARGET PAIRS")
    print("=" * 80)

    for i in range(len(build_names)):
        for j in range(i + 1, len(build_names)):
            cfg_a, build_a = build_names[i]
            cfg_b, build_b = build_names[j]
            
            cand_a = f"results/{build_a}/m09/degenerate/segment-MP_region_1/final_candidates.json"
            cand_b = f"results/{build_b}/m09/degenerate/segment-MP_region_1/final_candidates.json"
            
            if not os.path.exists(cand_a):
                for root, dirs, files in os.walk(f"results/{build_a}"):
                    if "final_candidates.json" in files:
                        cand_a = os.path.join(root, "final_candidates.json")
                        break

            if not os.path.exists(cand_b):
                for root, dirs, files in os.walk(f"results/{build_b}"):
                    if "final_candidates.json" in files:
                        cand_b = os.path.join(root, "final_candidates.json")
                        break

            out_dir = f"results/multiplex-{build_a}_x_{build_b}"
            
            print(f"\n[Multiplex Pair] {build_a}  <--->  {build_b}")
            print(f"  Input A: {cand_a}")
            print(f"  Input B: {cand_b}")
            print(f"  Output: {out_dir}")

            m11_cmd = f"python scripts/m11_multiplex_checker.py --file-a {cand_a} --file-b {cand_b} --outdir {out_dir} --num-candidates {args.num_candidates}"
            run_cmd(m11_cmd)

    # Step 3: Execute Module 11 Triplex Compatibility for 3 targets if present
    if len(build_names) == 3:
        print("\n" + "=" * 80)
        print("🧬 DETECTED 3 TARGETS — RUNNING TRIPLEX COMPATIBILITY SCREENING (ALL 3 VIRUSES)")
        print("=" * 80)
        
        cfg_a, build_a = build_names[0]
        cfg_b, build_b = build_names[1]
        cfg_c, build_c = build_names[2]
        
        def find_candidates(build):
            cand = f"results/{build}/m09/degenerate/segment-MP_region_1/final_candidates.json"
            if not os.path.exists(cand):
                for root, dirs, files in os.walk(f"results/{build}"):
                    if "final_candidates.json" in files:
                        return os.path.join(root, "final_candidates.json")
            return cand

        cand_a = find_candidates(build_a)
        cand_b = find_candidates(build_b)
        cand_c = find_candidates(build_c)
        
        out_dir = "results/multiplex-triplex-all-3"
        print(f"\n[Triplex Combo] {build_a} + {build_b} + {build_c}")
        print(f"  Input A: {cand_a}")
        print(f"  Input B: {cand_b}")
        print(f"  Input C: {cand_c}")
        print(f"  Output: {out_dir}")
        
        triplex_cmd = f"python scripts/m11_triplex_checker.py --file-a {cand_a} --file-b {cand_b} --file-c {cand_c} --outdir {out_dir} --num-candidates {args.num_candidates}"
        run_cmd(triplex_cmd)

    print("\n" + "=" * 80)
    print("✅ MULTI-TARGET BATCHING & MULTIPLEX SCREENING COMPLETED SUCCESSFULLY!")
    print("=" * 80)

if __name__ == "__main__":
    main()
