"""
Pandemic Preparedness Diagnostic Pipeline
Main Snakemake orchestrator for the Nextstrain-HyPhy-PrimerLab hybrid architecture.
"""

configfile: "config/config.yaml"

# Global namespace for run-specific build directory
BUILD = config["build_name"]

TARGET_GENES = config["hyphy"].get("target_genes", [])

rule all:
    input:
        # M06 Target
        f"results/{BUILD}/nextstrain/auspice_dataset.json",
        # M05 Targets (HyPhy)
        [f"results/{BUILD}/m05/fubar_{gene}.json" for gene in TARGET_GENES],
        [f"results/{BUILD}/m05/meme_{gene}.json" for gene in TARGET_GENES]

include: "rules/m01_acquisition.smk"
include: "rules/m02_m06_augur_core.smk"
include: "rules/m05_hyphy.smk"
# include: "rules/m07_bridge.smk"
# include: "rules/m08_primerlab.smk"
