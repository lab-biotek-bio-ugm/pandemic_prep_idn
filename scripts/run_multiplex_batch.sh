#!/usr/bin/env bash
# ==============================================================================
# PANDEMIC PREPAREDNESS: MULTI-TARGET BATCH & MULTIPLEX ORCHESTRATOR
# ==============================================================================
# Batch-executes Snakemake across multiple virus configs, then runs Module 11
# Multiplex Compatibility Screening across all generated candidate sets.
#
# Usage:
#   bash scripts/run_multiplex_batch.sh <config1.yaml> <config2.yaml> [snakemake_flags...]
#
# Example:
#   bash scripts/run_multiplex_batch.sh config/config_influenza_a_mp.yaml config/config_influenza_b_mp.yaml --use-conda -c 64
#   bash scripts/run_multiplex_batch.sh config/config_influenza_a_mp.yaml config/config_influenza_b_mp.yaml --unlock
# ==============================================================================

set -e

if [ "$#" -lt 2 ]; then
  echo "Usage: bash scripts/run_multiplex_batch.sh <config1.yaml> <config2.yaml> [extra snakemake flags...]"
  echo "Example: bash scripts/run_multiplex_batch.sh config/config_influenza_a_mp.yaml config/config_influenza_b_mp.yaml --use-conda -c 64"
  exit 1
fi

CONFIG_FILES=()
SNAKEMAKE_FLAGS=()

# Separate config files (.yaml) from extra snakemake flags
for arg in "$@"; do
  if [[ "$arg" == *.yaml ]] || [[ "$arg" == *.yml ]]; then
    CONFIG_FILES+=("$arg")
  else
    SNAKEMAKE_FLAGS+=("$arg")
  fi
done

if [ "${#CONFIG_FILES[@]}" -lt 2 ]; then
  echo "[Error] At least 2 config YAML files are required for multiplex screening."
  exit 1
fi

echo "================================================================================"
echo "🚀 MULTI-TARGET BATCH & MULTIPLEX ORCHESTRATOR"
echo "================================================================================"
echo "Configs to execute: ${CONFIG_FILES[*]}"
echo "Extra Snakemake flags: ${SNAKEMAKE_FLAGS[*]}"
echo "================================================================================"

BUILD_NAMES=()

# Step 1: Run Snakemake for each config
for cfg in "${CONFIG_FILES[@]}"; do
  echo -e "\n[Batch] Processing target config: $cfg"
  
  # Extract build_name from YAML config
  BUILD_NAME=$(grep -E "^build_name:" "$cfg" | awk '{print $2}' | tr -d '"' | tr -d "'")
  if [ -z "$BUILD_NAME" ]; then
    BUILD_NAME=$(basename "$cfg" .yaml)
  fi
  BUILD_NAMES+=("$BUILD_NAME")

  echo "[Batch] Target Build Name: $BUILD_NAME"
  echo "[Batch] Executing Snakemake..."
  
  snakemake --configfile "$cfg" "${SNAKEMAKE_FLAGS[@]}"
done

echo -e "\n================================================================================"
echo "🧬 STARTING MODULE 11 MULTIPLEX COMPATIBILITY SCREENING"
echo "================================================================================"

# Step 2: Run Module 11 across all target build pairs
NUM_BUILDS=${#BUILD_NAMES[@]}
for (( i=0; i<NUM_BUILDS; i++ )); do
  for (( j=i+1; j<NUM_BUILDS; j++ )); do
    BUILD_A="${BUILD_NAMES[$i]}"
    BUILD_B="${BUILD_NAMES[$j]}"
    
    CAND_A=$(find "results/$BUILD_A" -name "final_candidates.json" | head -n 1)
    CAND_B=$(find "results/$BUILD_B" -name "final_candidates.json" | head -n 1)
    
    if [ -z "$CAND_A" ] || [ ! -f "$CAND_A" ]; then
      echo "[Error] Could not find final_candidates.json in results/$BUILD_A"
      continue
    fi
    if [ -z "$CAND_B" ] || [ ! -f "$CAND_B" ]; then
      echo "[Error] Could not find final_candidates.json in results/$BUILD_B"
      continue
    fi
    
    OUT_DIR="results/multiplex-${BUILD_A}_x_${BUILD_B}"
    
    echo -e "\n[Multiplex Pair] $BUILD_A <---> $BUILD_B"
    echo "  Input A: $CAND_A"
    echo "  Input B: $CAND_B"
    echo "  Output:  $OUT_DIR"
    
    python scripts/m11_multiplex_checker.py --file-a "$CAND_A" --file-b "$CAND_B" --outdir "$OUT_DIR" --num-candidates 10
  done
done

echo -e "\n================================================================================"
echo "✅ BATCHING RUN AND MULTIPLEX COMPATIBILITY COMPLETED!"
echo "================================================================================"
