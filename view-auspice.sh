#!/bin/bash

# 1. Parse build name from config (defaulting to influenza-a-mp-sea)
CONFIG_FILE="config/config_influenza_a_mp.yaml"
BUILD_NAME="influenza-a-mp-sea" # fallback default

if [ -f "$CONFIG_FILE" ]; then
    # Try to parse build_name using grep and sed
    PARSED_BUILD=$(grep -E "^build_name:" "$CONFIG_FILE" | sed -E 's/build_name:[[:space:]]*"([^"]+)".*/\1/' | sed -E "s/build_name:[[:space:]]*'([^']+)'.*/\1/" | sed -E 's/build_name:[[:space:]]*([^[:space:]]+).*/\1/')
    if [ ! -z "$PARSED_BUILD" ]; then
        BUILD_NAME="$PARSED_BUILD"
    fi
fi

DATASET_DIR="results/${BUILD_NAME}/nextstrain"

echo "========================================================"
echo " Auspice Visualization Launcher"
echo "========================================================"
echo "Detected Build: $BUILD_NAME"
echo "Dataset Dir   : $DATASET_DIR"
echo "========================================================"

# 2. Check if pandemic_nextstrain environment exists
CONDA_ENV_NAME="pandemic_nextstrain"
ENV_EXISTS=$(conda env list | grep -E "^${CONDA_ENV_NAME}[[:space:]]")

if [ -z "$ENV_EXISTS" ]; then
    echo "WARNING: Conda environment '$CONDA_ENV_NAME' not found."
    echo "To create it, run:"
    echo "  conda env create -f envs/nextstrain.yaml"
    echo ""
    echo "Alternatively, if Snakemake has run before using --use-conda,"
    echo "you can find and activate the automatic Snakemake nextstrain conda env."
    exit 1
fi

# 3. Check if dataset directory exists
if [ ! -d "$DATASET_DIR" ]; then
    echo "WARNING: Dataset directory '$DATASET_DIR' does not exist."
    echo "Please run the pipeline first to generate the auspice dataset:"
    echo "  snakemake --use-conda --cores 8"
    echo ""
    echo "We will still attempt to launch Auspice, but it may be empty."
    echo "--------------------------------------------------------"
fi

# 4. Activate conda environment and run auspice
echo "Activating conda environment '$CONDA_ENV_NAME'..."
# Source conda sh to allow 'conda activate' inside script
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV_NAME"

echo "Launching Auspice view server..."
echo "If you are using VS Code, it will automatically forward port 4000."
echo "If not, remember to tunnel using: ssh -L 4000:localhost:4000 user@hpc"
echo "Open your browser at: http://localhost:4000"
echo "--------------------------------------------------------"
auspice view --dataset-dir "$DATASET_DIR"
