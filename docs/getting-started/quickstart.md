# Quickstart Guide

This guide will help you set up your local environment and execute your first diagnostic target discovery run on Windows WSL2 or Ubuntu.

---

## 1. Setup Windows Subsystem for Linux (WSL2)

If you are using Windows, you must configure a Linux environment to run bioconda packages (such as ViennaRNA, HyPhy, and MAFFT).

1. Open PowerShell as Administrator and enable WSL:
   ```powershell
   wsl --install -d Ubuntu-22.04
   ```
2. Restart your computer if prompted.
3. Launch Ubuntu from your Windows Start Menu, configure your username and password, and update packages:
   ```bash
   sudo apt update && sudo apt upgrade -y
   sudo apt install git unzip wget graphviz -y
   ```

---

## 2. Install Mambaforge (Conda)

Within your WSL2/Linux terminal, download and install Miniforge3 (with `mamba` pre-configured):

```bash
# Download the installer
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

# Run the installer script
bash Miniforge3-Linux-x86_64.sh -b

# Initialise mamba shell
~/miniforge3/bin/conda init bash
~/miniforge3/bin/mamba init bash

# Restart your terminal session
exec bash
```

Set up channels for bioinformatics dependency resolution:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

---

## 3. Clone Repository & Setup Snakemake

Clone the pipeline repository:

```bash
git clone https://github.com/lab-biotek-bio-ugm/pandemic_prep_idn.git
cd pandemic_prep_idn
```

Create an isolated Snakemake environment:

```bash
mamba create -c conda-forge -c bioconda -n snakemake snakemake=8.10.8 python=3.10 -y
conda activate snakemake
```

---

## 4. Run a Validation Dry-Run

Always perform a dry-run first to check that your paths, configurations, and DAG dependencies are properly resolved:

```bash
snakemake --configfile config/config_influenza_a_mp.yaml --dryrun --printshellcmds
```

If the dry-run succeeds, launch the complete workflow. Snakemake will automatically create isolated conda environments for every rule (this takes 10-15 minutes on the first run to build all environments):

```bash
snakemake --use-conda --cores 8 --configfile config/config_influenza_a_mp.yaml
```

Once completed, the final validated assays, specificity reports, and interactive HTML coverage heatmaps will be available in `results/influenza-a-mp-sea/m10/`.
