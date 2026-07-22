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

# Step 3: Run Triplex Compatibility if 3 target builds are present
if [ "$NUM_BUILDS" -eq 3 ]; then
  BUILD_A="${BUILD_NAMES[0]}"
  BUILD_B="${BUILD_NAMES[1]}"
  BUILD_C="${BUILD_NAMES[2]}"

  CAND_A=$(find "results/$BUILD_A" -name "final_candidates.json" | head -n 1)
  CAND_B=$(find "results/$BUILD_B" -name "final_candidates.json" | head -n 1)
  CAND_C=$(find "results/$BUILD_C" -name "final_candidates.json" | head -n 1)

  if [ -n "$CAND_A" ] && [ -n "$CAND_B" ] && [ -n "$CAND_C" ]; then
    OUT_DIR="results/multiplex-triplex-all-3"
    echo -e "\n================================================================================"
    echo "🧬 RUNNING TRIPLEX COMPATIBILITY SCREENING (ALL 3 VIRUSES)"
    echo "================================================================================"
    echo "  Input A: $CAND_A"
    echo "  Input B: $CAND_B"
    echo "  Input C: $CAND_C"
    echo "  Output:  $OUT_DIR"

    python scripts/m11_triplex_checker.py --file-a "$CAND_A" --file-b "$CAND_B" --file-c "$CAND_C" --outdir "$OUT_DIR" --num-candidates 10
  fi
fi

# Step 4: Automatically export supplementary data
echo -e "\n================================================================================"
echo "📊 CONSOLIDATING SUPPLEMENTARY EXCEL REPORTS"
echo "================================================================================"
python3 scripts/export_supplementary.py --results-dir results --outdir supplementary --only-degenerate

# Mirror to Windows host Downloads results folder if available (WSL fallback)
if [ -d "/mnt/c/Users/engki/Downloads/results" ]; then
  echo "[Orchestrator] Mirroring supplementary folder to Windows host..."
  mkdir -p "/mnt/c/Users/engki/Downloads/results/supplementary"
  cp -r supplementary/* "/mnt/c/Users/engki/Downloads/results/supplementary/"
fi

echo -e "\n================================================================================"
echo "✅ MULTI-TARGET BATCHING & MULTIPLEX SCREENING COMPLETED SUCCESSFULLY!"
echo "================================================================================"
