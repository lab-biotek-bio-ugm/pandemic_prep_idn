# Quickstart Guide

This guide will help you set up your environment and run the pipeline for the first time.

---

## Prerequisites

*   **Operating System**: Linux (Ubuntu 20.04+ recommended) or macOS. For Windows users, **WSL2** is required.
*   **Package Manager**: [Miniforge3/Mambaforge](https://github.com/conda-forge/miniforge) (recommended) or Miniconda. Mamba is strongly preferred due to its faster dependency resolution.

---

## 1. Setup Conda & Clone Repository

Clone the project repository and navigate into the workspace:

```bash
git clone https://github.com/lab-biotek-bio-ugm/pandemic_prep_idn.git
cd pandemic_prep_idn
```

Configure your conda shell (if not already done) and set up the channels:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

---

## 2. Install Snakemake

We recommend installing Snakemake in a dedicated environment using `mamba`:

```bash
mamba create -c conda-forge -c bioconda -n snakemake snakemake=8.10.8 python=3.10 -y
conda activate snakemake
```

---

## 3. Configure Pathogen Metadata & BLAST

The pipeline integrates with **local BLAST** for primer specificity checks. To set up a local database (e.g., human genome or targeted pathogen genomes):

1.  Download your reference DB (e.g., NCBI RefSeq).
2.  Index it using `makeblastdb`:
    ```bash
    makeblastdb -in reference.fasta -dbtype nucl -out data/blast_db/pathogen_db
    ```
3.  Configure the path in your profile under the `blast_verification` section (see the [Configuration Guide](configuration)).

---

## 4. Run the Pipeline

Always run a **dry-run** first to verify that Snakemake can construct the execution graph without syntax or path errors:

```bash
snakemake --configfile config/config_influenza_a_mp.yaml -n
```

To run the pipeline and allow Snakemake to automatically download and manage separate conda environments for each step, use:

```bash
snakemake --use-conda --cores 16 --configfile config/config_influenza_a_mp.yaml
```

### High-Performance Computing (HPC/Slurm)

If you are running the pipeline on a shared server with a Slurm scheduler:

```bash
sbatch slurm-run.sh
```
This script allocates resources, activates the conda environment, and executes the Snakemake workflow in the background.
