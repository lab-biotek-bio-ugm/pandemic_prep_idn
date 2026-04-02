"""
Pandemic Preparedness Diagnostic Pipeline
Main Snakemake orchestrator for the Nextstrain-HyPhy-PrimerLab hybrid architecture.
"""

configfile: "config/config.yaml"

# Global namespace for run-specific build directory
BUILD = config["build_name"]

rule all:
    input:
        f"results/{BUILD}/nextstrain/auspice_dataset.json",  # Clade-annotated phylogeny (M6)
        f"results/{BUILD}/nextstrain/aligned.fasta"          # Standardized MSA for HyPhy (M3)

include: "rules/m01_acquisition.smk"
include: "rules/m02_m06_augur_core.smk"
# include: "rules/m05_hyphy.smk"
# include: "rules/m07_bridge.smk"
# include: "rules/m08_primerlab.smk"
