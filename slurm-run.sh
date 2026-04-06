#!/bin/bash
#SBATCH --job-name=pandemic_run
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err
#SBATCH --partition=normal         # Sesuaikan dengan nama partisi HPC-mu (batch/compute/dll)
#SBATCH --nodes=1
#SBATCH --ntasks=1                 # Harus 1, Snakemake bertindak sebagai orchestrator tunggal
#SBATCH --cpus-per-task=16         # Total Cores yang di-request, diteruskan ke Snakemake
#SBATCH --mem=32GB                 # Optimal untuk MSA & Tree inference (FastTree/IQ-TREE)
#SBATCH --time=08:00:00            # Batas waktu wall-clock: 8 Jam

# ==============================================================================
# Pandemic Preparedness Pipeline - SLURM Job Script
# ==============================================================================

# Pastikan environment kondusif untuk IO speed
export TMPDIR=/tmp

# Load modul Miniconda (jika HPC menggunakan Lmod environment modules)
# Hilangkan komentar pada baris di bawah jika direkomendasikan admin
# module load miniconda3

echo "======================================================================"
echo "Starting Snakemake on Node: $(hostname)"
echo "Allocated Cores: $SLURM_CPUS_PER_TASK"
echo "Job ID: $SLURM_JOB_ID"
echo "======================================================================"

# Eksekusi Pipeline (Snakemake otomatis menggunakan jumlah Core dari SLURM)
snakemake \
  --use-conda \
  --cores $SLURM_CPUS_PER_TASK \
  --keep-going \
  --printshellcmds

echo "======================================================================"
echo "Pipeline Execution Completed."
echo "======================================================================"
