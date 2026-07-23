#!/usr/bin/env python3
"""
Master Multi-Target Multiplex Orchestrator
Batch-executes Snakemake pipelines across multiple virus targets/configs,
and triggers Module 11 Multiplex Compatibility Screening across all target pairs.
"""

import glob
import os
import re
import shlex
import sys
import yaml
import argparse
import subprocess

def run_cmd(cmd):
    """Run a shell command string. Only for the Snakemake invocation, whose
    --snakemake-args are a free-form user-supplied flag string."""
    print(f"\n[Batch Orchestrator] Running: {cmd}")
    res = subprocess.run(cmd, shell=True)
    if res.returncode != 0:
        print(f"[Error] Command failed with code {res.returncode}: {cmd}")
        sys.exit(res.returncode)


def run_argv(argv, fatal=True):
    """Run a command given as an argument LIST — no shell, so values containing
    spaces or shell metacharacters (target_label: "Influenza A (H3N2)", a dye like
    "Texas Red") are passed through intact instead of being word-split or
    interpreted. Returns True on success.

    Building these as f-string shell commands is the same defect class as the old
    bash orchestrator's `awk '{print $2}'` build_name parsing.
    """
    print(f"\n[Batch Orchestrator] Running: {' '.join(shlex.quote(a) for a in argv)}")
    res = subprocess.run(argv)
    if res.returncode != 0:
        msg = f"[Error] Command failed with code {res.returncode}"
        if fatal:
            print(msg)
            sys.exit(res.returncode)
        print(f"{msg} — continuing with the remaining targets.")
        return False
    return True

# Sentinel used when a target's detection fluorophore cannot be determined from
# config. Passing this (rather than letting the checker fall back to FAM/HEX)
# makes check_dye_separation report UNKNOWN instead of silently asserting that
# two targets are spectrally separated.
UNRESOLVED_DYE = "UNRESOLVED"


def sanitize_build_name(name: str) -> str:
    """Make a config-supplied build_name safe to interpolate into result paths.

    build_name comes from user YAML and is used to build `results/<name>/...` and
    `results/multiplex-<a>_x_<b>/`. A value like `../../etc` would write outside
    results/. Keep only path-safe characters and reject empties.
    """
    cleaned = re.sub(r"[^A-Za-z0-9._-]", "-", str(name or "").strip())
    cleaned = cleaned.strip(".-") or "unnamed-build"
    if cleaned != str(name):
        print(f"[Batch Orchestrator] build_name {name!r} sanitized to {cleaned!r}")
    return cleaned


def find_candidates(build: str) -> str:
    """Resolve a build's degenerate M09 final_candidates.json.

    Mirrors _find_m09_candidates() in rules/m11_multiplex.smk: prefer the
    lexicographically first `segment-*_region_1` match so the Snakemake and
    orchestrator entrypoints screen the SAME segment. Previously this preferred
    segment-MP while the rule took sorted(glob)[0], so a multi-segment influenza
    build was screened on different segments depending on how M11 was launched.
    """
    base = f"results/{build}/m09/degenerate"
    hits = sorted(glob.glob(f"{base}/segment-*_region_1/final_candidates.json"))
    if hits:
        return hits[0]
    # Fall back to any final_candidates.json under the build.
    for root, _dirs, files in os.walk(f"results/{build}"):
        if "final_candidates.json" in files:
            return os.path.join(root, "final_candidates.json")
    return f"{base}/segment-UNKNOWN_region_1/final_candidates.json"


def parse_build_name(config_path):
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)
    return sanitize_build_name(
        cfg.get("build_name", os.path.basename(config_path).replace(".yaml", ""))
    )


def resolve_dye(cfg, config_path):
    """Resolve a target's probe fluorophore.

    Order: explicit `dye` / `primerlab.dye` on the target config, else the
    `parameters.probe.labels.fluorophore` of the referenced primerlab template.
    Returns UNRESOLVED_DYE if it cannot be determined, so the multiplex dye-channel
    check flags it rather than assuming a default.
    """
    primerlab = cfg.get("primerlab") or {}
    explicit = cfg.get("dye") or primerlab.get("dye")
    if explicit:
        return str(explicit)

    template = primerlab.get("template")
    if not template:
        print(f"[WARN] {config_path}: no primerlab.template — cannot resolve dye.")
        return UNRESOLVED_DYE
    try:
        with open(template, "r") as f:
            tpl = yaml.safe_load(f) or {}
        dye = (((tpl.get("parameters") or {}).get("probe") or {})
               .get("labels") or {}).get("fluorophore")
        if dye:
            return str(dye)
        print(f"[WARN] {template}: parameters.probe.labels.fluorophore missing.")
    except Exception as e:
        print(f"[WARN] Could not read dye from template '{template}': {e}")
    return UNRESOLVED_DYE


def parse_target_meta(config_path):
    """Return (build_name, dye, label) for one target config."""
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f) or {}
    build = sanitize_build_name(
        cfg.get("build_name", os.path.basename(config_path).replace(".yaml", ""))
    )
    label = cfg.get("target_label") or build
    return build, resolve_dye(cfg, config_path), label

def main():
    parser = argparse.ArgumentParser(description="Multi-Target Batching & Multiplex Orchestrator")
    parser.add_argument("--configs", nargs="+", required=True, help="List of target config files (e.g., config_influenza_a_mp.yaml config_influenza_b_mp.yaml)")
    parser.add_argument("--cores", type=int, default=64, help="Number of CPU cores for Snakemake execution")
    parser.add_argument("--skip-snakemake", action="store_true", help="Skip Snakemake runs if result folders already exist")
    parser.add_argument("--snakemake-args", type=str, default="--use-conda --rerun-incomplete", help="Additional flags to pass directly to Snakemake (e.g., '--use-conda --unlock --rerun-incomplete')")
    parser.add_argument("--num-candidates", type=int, default=10, help="Number of top candidates per target to screen for multiplexing")
    parser.add_argument("--outdir", default="supplementary", help="Output directory for the supplementary workbooks")
    parser.add_argument("--skip-export", action="store_true", help="Skip the supplementary export step")
    parser.add_argument("--mirror-dir", default=os.environ.get("SUPPLEMENTARY_MIRROR_DIR", ""),
                        help="Optional directory to mirror the supplementary folder into "
                             "(default: $SUPPLEMENTARY_MIRROR_DIR, unset = no mirroring)")
    args = parser.parse_args()

    build_names = []
    targets = {}  # build_name -> {"dye": ..., "label": ...}

    print("=" * 80)
    print("🚀 PANDEMIC PREPAREDNESS MULTI-TARGET BATCH ORCHESTRATOR")
    print("=" * 80)

    # Step 1: Run Snakemake for each target config
    for idx, cfg in enumerate(args.configs, 1):
        build_name, dye, label = parse_target_meta(cfg)
        # Two configs sharing a build_name would overwrite each other's dye/label,
        # so a pair would be screened and reported with the WRONG fluorophore.
        if build_name in targets:
            print(f"[Error] Duplicate build_name '{build_name}' (config: {cfg}). "
                  f"Each target config must declare a unique build_name — otherwise "
                  f"their results directories and dye assignments collide.")
            sys.exit(1)
        build_names.append((cfg, build_name))
        targets[build_name] = {"dye": dye, "label": label}
        print(f"\n[Target {idx}/{len(args.configs)}] Config: {cfg} ---> Build: {build_name} "
              f"(dye: {dye})")
        
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
            
            cand_a = find_candidates(build_a)
            cand_b = find_candidates(build_b)

            out_dir = f"results/multiplex-{build_a}_x_{build_b}"
            meta_a = targets.get(build_a, {})
            meta_b = targets.get(build_b, {})
            dye_a = meta_a.get("dye", UNRESOLVED_DYE)
            dye_b = meta_b.get("dye", UNRESOLVED_DYE)

            print(f"\n[Multiplex Pair] {build_a} ({dye_a})  <--->  {build_b} ({dye_b})")
            print(f"  Input A: {cand_a}")
            print(f"  Input B: {cand_b}")
            print(f"  Output: {out_dir}")

            # Fail soft: a build without M09 output must not abort the whole batch
            # (the remaining pairs, the triplex and the export step would all be lost).
            missing = [p for p in (cand_a, cand_b) if not os.path.exists(p)]
            if missing:
                print(f"  [SKIP] No final_candidates.json for: {', '.join(missing)}")
                continue

            run_argv([
                sys.executable, "scripts/m11_multiplex_checker.py",
                "--file-a", cand_a, "--file-b", cand_b, "--outdir", out_dir,
                "--num-candidates", str(args.num_candidates),
                "--dye-a", dye_a, "--dye-b", dye_b,
                "--label-a", str(meta_a.get("label", build_a)),
                "--label-b", str(meta_b.get("label", build_b)),
            ], fatal=False)

    # Step 3: Execute Module 11 Triplex Compatibility for 3 targets if present
    if len(build_names) == 3:
        print("\n" + "=" * 80)
        print("🧬 DETECTED 3 TARGETS — RUNNING TRIPLEX COMPATIBILITY SCREENING (ALL 3 VIRUSES)")
        print("=" * 80)
        
        cfg_a, build_a = build_names[0]
        cfg_b, build_b = build_names[1]
        cfg_c, build_c = build_names[2]
        
        cand_a = find_candidates(build_a)
        cand_b = find_candidates(build_b)
        cand_c = find_candidates(build_c)
        
        out_dir = "results/multiplex-triplex-all-3"
        meta_a = targets.get(build_a, {})
        meta_b = targets.get(build_b, {})
        meta_c = targets.get(build_c, {})
        dye_a = meta_a.get("dye", UNRESOLVED_DYE)
        dye_b = meta_b.get("dye", UNRESOLVED_DYE)
        dye_c = meta_c.get("dye", UNRESOLVED_DYE)

        print(f"\n[Triplex Combo] {build_a} ({dye_a}) + {build_b} ({dye_b}) + {build_c} ({dye_c})")
        print(f"  Input A: {cand_a}")
        print(f"  Input B: {cand_b}")
        print(f"  Input C: {cand_c}")
        print(f"  Output: {out_dir}")

        missing = [p for p in (cand_a, cand_b, cand_c) if not os.path.exists(p)]
        if missing:
            print(f"  [SKIP] Triplex needs all three targets; missing: {', '.join(missing)}")
        else:
            run_argv([
                sys.executable, "scripts/m11_triplex_checker.py",
                "--file-a", cand_a, "--file-b", cand_b, "--file-c", cand_c,
                "--outdir", out_dir,
                "--num-candidates", str(args.num_candidates),
                "--dye-a", dye_a, "--dye-b", dye_b, "--dye-c", dye_c,
                "--label-a", str(meta_a.get("label", build_a)),
                "--label-b", str(meta_b.get("label", build_b)),
                "--label-c", str(meta_c.get("label", build_c)),
            ], fatal=False)

    # Step 4: Consolidate supplementary workbooks, restricted to THIS run's targets
    # so builds from previous runs cannot leak into the exported reports.
    if not args.skip_export:
        print("\n" + "=" * 80)
        print("📊 CONSOLIDATING SUPPLEMENTARY EXCEL REPORTS")
        print("=" * 80)

        run_argv([
            sys.executable, "scripts/export_supplementary.py",
            "--results-dir", "results", "--outdir", args.outdir,
            "--only-degenerate",
            "--active-targets", *[b for _, b in build_names],
        ])

        if args.mirror_dir:
            print(f"[Orchestrator] Mirroring supplementary folder to {args.mirror_dir} ...")
            os.makedirs(args.mirror_dir, exist_ok=True)
            import shutil
            shutil.copytree(args.outdir, args.mirror_dir, dirs_exist_ok=True)

    print("\n" + "=" * 80)
    print("✅ MULTI-TARGET BATCHING & MULTIPLEX SCREENING COMPLETED SUCCESSFULLY!")
    print("=" * 80)

if __name__ == "__main__":
    main()
