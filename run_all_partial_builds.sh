#!/usr/bin/env bash
#
# Run the all-partial builds for Influenza A MP and B MP.
#
# "all-partial" = focal + contextual subsampling: dense focal Indonesian strains
# PLUS a sparse global background (every other country, few sequences per group).
# The subsampling.focal / subsampling.contextual blocks are enabled in these
# configs; the plain top-level subsampling.query is ignored in this mode. Every
# other threshold and design parameter is identical to the indo-only and sea-only
# builds, so all three are directly comparable.
#
# Usage:
#   ./run_all_partial_builds.sh [CORES] [extra snakemake args...]
#
# Examples:
#   ./run_all_partial_builds.sh                 # 16 cores, real run
#   ./run_all_partial_builds.sh 32              # 32 cores
#   ./run_all_partial_builds.sh 16 --dry-run    # preview the DAG only
#   DRYRUN=1 ./run_all_partial_builds.sh        # same as passing --dry-run
#
# Env vars:
#   SNAKEMAKE_ENV   conda env to activate (default: snakemake); set to "" to skip
#
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"

CORES="${1:-16}"
if [[ $# -gt 0 ]]; then shift; fi
EXTRA_ARGS=("$@")
if [[ "${DRYRUN:-0}" == "1" ]]; then EXTRA_ARGS+=("--dry-run"); fi

CONFIGS=(
  "config/config_influenza_a_mp_all_partial.yaml"
  "config/config_influenza_b_mp_all_partial.yaml"
)

# --- Activate the snakemake conda env (skip if SNAKEMAKE_ENV="") ---------------
SNAKEMAKE_ENV="${SNAKEMAKE_ENV-snakemake}"
if [[ -n "$SNAKEMAKE_ENV" ]]; then
  if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$SNAKEMAKE_ENV" || {
      echo "WARN: could not activate conda env '$SNAKEMAKE_ENV'; assuming snakemake is already on PATH." >&2
    }
  else
    echo "WARN: conda not found; assuming snakemake is already on PATH." >&2
  fi
fi

command -v snakemake >/dev/null 2>&1 || { echo "ERROR: snakemake not on PATH." >&2; exit 1; }
echo "snakemake: $(snakemake --version)  |  cores: $CORES"
echo

# --- Pre-flight: how many focal vs contextual sequences will feed the build? ---
# Prints the focal (Indonesia) cohort size and the global contextual pool size
# from the metadata, so the focal/contextual mix is visible before a long run.
# Never fatal -- it only warns.
preflight() {
  local cfg="$1"
  python3 - "$cfg" <<'PY' || echo "  (pre-flight skipped: metadata not readable yet)"
import sys, os, yaml
cfg = yaml.safe_load(open(sys.argv[1]))
build = cfg.get("build_name", "")
raw = cfg.get("local_input", {}).get("metadata")
sub = cfg["subsampling"]
fq = sub.get("focal", {}).get("query")
cq = sub.get("contextual", {}).get("query")
import pandas as pd
# The real filter (augur_filter, M02) runs on the M01-enriched metadata, which
# has the normalized 'country' column. The raw GISAID export only has 'Location',
# so prefer the enriched file when a previous run already produced it.
enriched = f"results/{build}/m01/metadata_enriched.tsv"
meta, sep = (enriched, "\t") if os.path.exists(enriched) else (raw, ",")
try:
    df = pd.read_csv(meta, sep=sep, dtype=str, low_memory=False)
except Exception as e:
    print(f"  metadata not loaded ({meta}): {e}"); sys.exit(0)
if "country" not in df.columns:
    print(f"  {meta} has no normalized 'country' column yet (raw GISAID export "
          f"uses 'Location'). M01 derives 'country' at run time, so the focal/"
          f"contextual country filter applies during the run. Skipping pre-count.")
    sys.exit(0)
print(f"  metadata: {meta}  (total {len(df)} rows)")
for label, q in [("focal", fq), ("contextual", cq)]:
    if not q:
        continue
    try:
        n = len(df.query(q))
    except Exception as e:
        print(f"  {label} query failed ({q}): {e}")
        continue
    print(f"  {label:11s} '{q}' -> {n} sequences")
top = df["country"].value_counts().head(8)
print("  top countries: " + ", ".join(f"{c}={int(n)}" for c, n in top.items()))
PY
}

echo "=== Pre-flight focal/contextual check ==="
for cfg in "${CONFIGS[@]}"; do
  echo "[$cfg]"
  preflight "$cfg"
  echo
done

# --- Run the builds -----------------------------------------------------------
for cfg in "${CONFIGS[@]}"; do
  echo "=================================================================="
  echo "Running: $cfg"
  echo "=================================================================="
  snakemake \
    --use-conda \
    --cores "$CORES" \
    --rerun-incomplete \
    --printshellcmds \
    --configfile "$cfg" \
    "${EXTRA_ARGS[@]}"
  echo
done

# --- Multiplex (A x B) + supplementary package --------------------------------
# `rule all` above only builds the per-target track (m01-m10). The cross-build
# multiplex screen and the supplementary workbooks come from the orchestrator,
# which pairs the two builds and exports to a scenario-specific folder. The
# builds are already done above, so --skip-snakemake goes straight to M12 +
# multiplex + supplementary export. Skipped on a dry run.
if [[ "${DRYRUN:-0}" != "1" && " ${EXTRA_ARGS[*]} " != *" --dry-run "* ]]; then
  echo "=================================================================="
  echo "Multiplex + supplementary (Flu A x Flu B, all-partial) -> supplementary_all_partial/"
  echo "=================================================================="
  python scripts/run_multiplex_batch.py \
    --configs "${CONFIGS[@]}" \
    --cores "$CORES" \
    --skip-snakemake \
    --outdir supplementary_all_partial
  echo
fi

echo "Done. Outputs:"
echo "  results/influenza-a-mp-all-partial/   results/influenza-b-mp-all-partial/"
echo "  results/multiplex-influenza-a-mp-all-partial_x_influenza-b-mp-all-partial/"
echo "  supplementary_all_partial/   (Excel workbooks + qc_dossier)"
echo
echo "Compare against the indo-only and sea-only builds (m10 coverage_summary.json,"
echo "m07 bridge_report.json, m08 best_primers.json) to see how a global focal+"
echo "contextual database shifts the primers and their coverage."
