"""
Pandemic Preparedness Diagnostic Pipeline
Main Snakemake orchestrator for the Nextstrain-HyPhy-PrimerLab hybrid architecture.
"""

configfile: "config/config.yaml"

# Global namespace for run-specific build directory
BUILD = config["build_name"]

def get_target_genes(config):
    genes = config["hyphy"].get("target_genes", ["ALL"])
    if "ALL" in genes:
        genes = []
        gff_path = config["references"]["annotation"]
        with open(gff_path, "r") as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split("\t")
                if len(parts) >= 9 and parts[2] == "gene":
                    for attr in parts[8].split(";"):
                        attr = attr.strip()
                        if attr.startswith("gene="):
                            name = attr.split("=")[1].strip('"')
                            if name not in genes:
                                genes.append(name)
    return genes

TARGET_GENES = get_target_genes(config)

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
