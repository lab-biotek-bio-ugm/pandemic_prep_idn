#!/bin/bash
#SBATCH --job-name=pandemic_pipeline
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err
#SBATCH --partition=normal         # Adjust to your HPC partition name (e.g., batch, compute, etc.)
#SBATCH --nodes=1
#SBATCH --ntasks=1                 # Fixed at 1: Snakemake acts as the single orchestrator
#SBATCH --cpus-per-task=16         # Total cores requested, passed to Snakemake
#SBATCH --mem=32GB                 # Optimal for MSA & Tree inference (FastTree/IQ-TREE)
#SBATCH --time=08:00:00            # Wall-clock limit: 8 hours

# SLURM Execution Template for Pandemic Preparedness Pipeline
# Note: Use this script if your HPC environment supports the Slurm Workload Manager.

# Standard I/O optimization
export TMPDIR=/tmp

# Load Miniconda/Anaconda if required by your admin
# module load miniconda3

echo "----------------------------------------------------------------------"
echo "Starting Snakemake on Node: $(hostname)"
echo "Allocated Cores: $SLURM_CPUS_PER_TASK"
echo "Job ID: $SLURM_JOB_ID"
echo "----------------------------------------------------------------------"

# Execute Pipeline
# Snakemake automatically utilizes the number of cores allocated by SLURM
snakemake \
  --use-conda \
  --cores $SLURM_CPUS_PER_TASK \
  --keep-going \
  --printshellcmds

echo "----------------------------------------------------------------------"
echo "Pipeline Execution Completed."
echo "----------------------------------------------------------------------"
