#!/usr/bin/env bash
#
# Run the Southeast-Asia (sea-only) builds for Influenza A MP and B MP.
#
# These configs differ from the indo-only builds ONLY in the geographic
# subsampling filter (subsampling.query) and the build_name (output directory);
# every threshold and design parameter is identical, so the sea-only vs indo-only
# comparison is apples-to-apples.
#
# Usage:
#   ./run_sea_builds.sh [CORES] [extra snakemake args...]
#
# Examples:
#   ./run_sea_builds.sh                 # 16 cores, real run
#   ./run_sea_builds.sh 32              # 32 cores
#   ./run_sea_builds.sh 16 --dry-run    # preview the DAG only
#   DRYRUN=1 ./run_sea_builds.sh        # same as passing --dry-run
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
  "config/config_influenza_a_mp_sea.yaml"
  "config/config_influenza_b_mp_sea.yaml"
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

# --- Pre-flight: which SEA countries actually exist in each metadata file? -----
# The country strings in the query are best-effort GISAID spellings. This prints
# what the filter will actually match so you can adjust the config before a long
# run. Never fatal -- it only warns.
preflight() {
  local cfg="$1"
  python3 - "$cfg" <<'PY' || echo "  (pre-flight skipped: metadata not readable yet)"
import sys, os, yaml
cfg = yaml.safe_load(open(sys.argv[1]))
build = cfg.get("build_name", "")
raw = cfg.get("local_input", {}).get("metadata")
q = cfg["subsampling"]["query"]
# extract the quoted country names from the "in [...]" list
import re
wanted = re.findall(r"'([^']+)'", q)
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
          f"uses 'Location'). M01 derives 'country' at run time, so the country "
          f"filter ({', '.join(wanted)}) applies during the run. Skipping pre-count.")
    sys.exit(0)
present = df["country"].value_counts()
hit = [(c, int(present.get(c, 0))) for c in wanted if present.get(c, 0) > 0]
miss = [c for c in wanted if present.get(c, 0) == 0]
total = sum(n for _, n in hit)
print(f"  metadata: {meta}")
print(f"  matched SEA countries ({total} seqs): " + ", ".join(f"{c}={n}" for c, n in hit))
if miss:
    print(f"  NOT present (spelling differs or absent): {', '.join(miss)}")
    # surface nearby country names that look South-East Asian but were not matched
    import difflib
    others = [c for c in present.index if c not in wanted]
    near = [c for w in miss for c in difflib.get_close_matches(w, others, n=1, cutoff=0.7)]
    if near:
        print(f"  possible spelling matches to add to the query: {sorted(set(near))}")
PY
}

echo "=== Pre-flight country check ==="
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
  echo "Multiplex + supplementary (Flu A x Flu B, sea-only) -> supplementary_sea/"
  echo "=================================================================="
  python scripts/run_multiplex_batch.py \
    --configs "${CONFIGS[@]}" \
    --cores "$CORES" \
    --skip-snakemake \
    --outdir supplementary_sea
  echo
fi

echo "Done. Outputs:"
echo "  results/influenza-a-mp-sea-only/   results/influenza-b-mp-sea-only/"
echo "  results/multiplex-influenza-a-mp-sea-only_x_influenza-b-mp-sea-only/"
echo "  supplementary_sea/   (Excel workbooks + qc_dossier)"
echo
echo "Compare against the indo-only builds (e.g. m10 coverage_summary.json,"
echo "m07 bridge_report.json, m08 best_primers.json) to see how the design DB"
echo "changes the primers and their coverage."
