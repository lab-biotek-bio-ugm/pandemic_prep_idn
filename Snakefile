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
HYPHY_METHODS = config["hyphy"].get("methods", ["FUBAR"])

rule all:
    input:
        # M02 Target (QC Report)
        f"results/{BUILD}/m02/qc_report.json",
        # M06 Target
        f"results/{BUILD}/nextstrain/auspice_dataset.json",
        # M05 Targets (HyPhy — FUBAR and MEME always, GARD conditional)
        [f"results/{BUILD}/m05/fubar_{gene}.json" for gene in TARGET_GENES],
        [f"results/{BUILD}/m05/meme_{gene}.json" for gene in TARGET_GENES],
        ([f"results/{BUILD}/m05/gard_{gene}.json" for gene in TARGET_GENES]
         if "GARD" in HYPHY_METHODS else []),
        # M07 Target (Bridge — primer target manifest)
        f"results/{BUILD}/m07/primer_targets_manifest.yaml",

include: "rules/m01_acquisition.smk"
include: "rules/m02_m06_augur_core.smk"
include: "rules/m05_hyphy.smk"
include: "rules/m07_bridge.smk"
# include: "rules/m08_primerlab.smk"

# DEFERRED DEVELOPMENT ITEMS — Tracked here to avoid loss across audits.
#
# TODO [DEV-1] Route B: Clade-Specific Primer Targets
#   What: Add a second route in M07 that identifies clade-defining SNPs
#         and generates SEQUENCE_TARGET constraints for clade-specific
#         primer design (genotyping: "which variant?", not just "which virus?").
#   Why deferred: Primary milestone is universal detection + multiplex RAA.
#                 Route B requires clade_snp_finder.py + consensus per clade.
#   Prereqs: DEV-4 (consensus per clade) must be implemented first.
#   Priority: Phase 2 — after M07 Route A is validated.
#
# TODO [DEV-3] Protein Phylogeny (Antigen/Serological Track)
#   What: Build a separate protein-level phylogenetic tree from MACSE AA
#         output using IQ-TREE2 with protein substitution models.
#   Why deferred: Current protein/codon data is only used for mutation
#                 pressure analysis (dN/dS via HyPhy). Protein phylogeny
#                 is needed for serological/antigen-based applications,
#                 which is not the focus of Phase 1 (nucleotide diagnostics).
#   Priority: Phase 2 — when antigen/epitope prediction is in scope.
#
# TODO [DEV-4] Consensus Sequence per Clade
#   What: Generate one representative consensus FASTA per clade by
#         majority-vote at each position across all sequences in that clade.
#         Output: clade_X_consensus.fasta for each clade detected by Nextclade.
#   Why deferred: Only needed for DEV-1 (Route B clade-specific primers).
#                 Universal primer design does not require per-clade consensus.
#   Prereqs: Clade labels from Nextclade (already available in metadata).
#   Priority: Phase 2 — together with DEV-1.
#
# TODO [DEV-5] Multi-Target Batch Run
#   What: Support a targets[] array in config.yaml to run the full pipeline
#         for multiple viruses in parallel (e.g., Influenza A, B, SARS-CoV-2)
#         in a single snakemake invocation. Each virus is analyzed independently;
#         results merge only at M08 (multiplex compatibility check).
#   Why deferred: Manual per-virus runs are sufficient for development phase.
#                 Multi-target adds config complexity (templates, per-virus
#                 reference genomes, segment handling) without functional gain.
#   Priority: Phase 3 — after single-virus pipeline is fully validated.
#
