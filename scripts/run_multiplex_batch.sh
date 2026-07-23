#!/usr/bin/env bash
# ==============================================================================
# PANDEMIC PREPAREDNESS: MULTI-TARGET BATCH & MULTIPLEX ORCHESTRATOR (wrapper)
# ==============================================================================
# Thin wrapper around scripts/run_multiplex_batch.py.
#
# This script previously re-implemented the whole orchestration in bash. Keeping
# two copies caused them to drift (e.g. the pairwise ΔTm threshold diverged from
# the documented value, and only one copy learned about dye-channel wiring), and
# the bash copy parsed build_name with `awk '{print $2}'`, which truncates any
# quoted value containing spaces. All orchestration now lives in the Python
# implementation; this wrapper only preserves the documented CLI.
#
# Usage:
#   bash scripts/run_multiplex_batch.sh <config1.yaml> <config2.yaml> [snakemake_flags...]
#
# Example:
#   bash scripts/run_multiplex_batch.sh config/config_influenza_a_mp.yaml config/config_influenza_b_mp.yaml --use-conda -c 64
#   bash scripts/run_multiplex_batch.sh config/config_influenza_a_mp.yaml config/config_influenza_b_mp.yaml config/config_measles_all_partial.yaml -c 64
# ==============================================================================

set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage: bash scripts/run_multiplex_batch.sh <config1.yaml> <config2.yaml> [extra snakemake flags...]"
  echo "Example: bash scripts/run_multiplex_batch.sh config/config_influenza_a_mp.yaml config/config_influenza_b_mp.yaml --use-conda -c 64"
  exit 1
fi

CONFIG_FILES=()
SNAKEMAKE_FLAGS=()
CORES=""

# Separate config files (.yaml/.yml) from extra snakemake flags.
# `-c N` / `--cores N` is forwarded to the orchestrator's --cores.
prev=""
for arg in "$@"; do
  if [[ "$arg" == *.yaml ]] || [[ "$arg" == *.yml ]]; then
    CONFIG_FILES+=("$arg")
  elif [[ "$prev" == "-c" || "$prev" == "--cores" ]]; then
    CORES="$arg"
  elif [[ "$arg" == "-c" || "$arg" == "--cores" ]]; then
    :  # consumed via $prev on the next iteration
  else
    SNAKEMAKE_FLAGS+=("$arg")
  fi
  prev="$arg"
done

if [ "${#CONFIG_FILES[@]}" -lt 2 ]; then
  echo "[Error] At least 2 config YAML files are required for multiplex screening."
  exit 1
fi

CMD=(python3 scripts/run_multiplex_batch.py --configs "${CONFIG_FILES[@]}")

if [ -n "$CORES" ]; then
  CMD+=(--cores "$CORES")
fi

if [ "${#SNAKEMAKE_FLAGS[@]}" -gt 0 ]; then
  # Must be the --opt=value form. Passed as two argv elements, argparse sees a
  # value that begins with '-' (e.g. "--use-conda"), treats it as another option
  # and aborts with "expected one argument". The '=' form has no such ambiguity.
  CMD+=("--snakemake-args=${SNAKEMAKE_FLAGS[*]}")
fi

echo "================================================================================"
echo "🚀 MULTI-TARGET BATCH & MULTIPLEX ORCHESTRATOR (delegating to Python)"
echo "================================================================================"
echo "Configs: ${CONFIG_FILES[*]}"
[ -n "$CORES" ] && echo "Cores:   $CORES"
[ "${#SNAKEMAKE_FLAGS[@]}" -gt 0 ] && echo "Flags:   ${SNAKEMAKE_FLAGS[*]}"
echo "================================================================================"

exec "${CMD[@]}"
