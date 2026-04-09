"""
Module 2-6: Augur Core (Nextstrain Ecosystem)
Executes the Nextstrain spatio-temporal filtering, multiple sequence alignment, 
phylogenetic tree inference, mutational trait analysis, and JSON payload export.
"""

# Import the global BUILD variable defined in the main Snakefile
BUILD = config["build_name"]

rule qc_prefilter:
    """Prefilters raw sequences: removes high N-content and duplicate sequences.
    Produces a QC report JSON with retention statistics.
    """
    input:
        fasta=f"results/{BUILD}/m01/raw_sequences.fasta"
    output:
        fasta=f"results/{BUILD}/m02/prefiltered_sequences.fasta",
        report=f"results/{BUILD}/m02/qc_report.json"
    params:
        max_n_pct=config.get("qc", {}).get("max_n_content", 5.0)
    conda:
        "../envs/nextstrain.yaml"
    log:
        f"logs/{BUILD}/m02_qc_prefilter.log"
    shell:
        """
        mkdir -p results/{BUILD}/m02
        python scripts/qc_prefilter.py \
            --input-fasta {input.fasta} \
            --output-fasta {output.fasta} \
            --report {output.report} \
            --max-n-pct {params.max_n_pct} 2>&1 | tee {log}
        """

rule augur_filter:
    input:
        sequences=f"results/{BUILD}/m02/prefiltered_sequences.fasta",
        metadata=f"results/{BUILD}/m01/metadata_enriched.tsv"
    output:
        sequences=f"results/{BUILD}/m02/filtered_sequences.fasta"
    params:
        group_by=config["subsampling"]["group_by"],
        seq_per_group=config["subsampling"]["sequences_per_group"],
        min_length=config["subsampling"]["min_length"]
    conda:
        "../envs/nextstrain.yaml"
    log:
        f"logs/{BUILD}/m02_augur_filter.log"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --group-by {params.group_by} \
            --sequences-per-group {params.seq_per_group} \
            --min-length {params.min_length} \
            --exclude-ambiguous-dates-by any \
            --output-sequences {output.sequences} 2> {log}
        """

rule augur_align:
    input:
        sequences=f"results/{BUILD}/m02/filtered_sequences.fasta",
        reference=config["references"]["genome"]
    output:
        alignment=f"results/{BUILD}/m03/aligned.fasta"
    threads:
        8
    conda:
        "../envs/nextstrain.yaml"
    log:
        f"logs/{BUILD}/m03_augur_align.log"
    shell:
        """
        mkdir -p results/{BUILD}/m03
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --nthreads {threads} 2> {log}
        """

rule augur_tree:
    input:
        alignment=f"results/{BUILD}/m03/aligned.fasta"
    output:
        tree=f"results/{BUILD}/m04/tree_raw.nwk"
    threads:
        8
    conda:
        "../envs/nextstrain.yaml"
    log:
        f"logs/{BUILD}/m04_augur_tree.log"
    shell:
        """
        mkdir -p results/{BUILD}/m04
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method iqtree \
            --nthreads {threads} 2> {log}
        """

rule augur_refine:
    input:
        tree=f"results/{BUILD}/m04/tree_raw.nwk",
        alignment=f"results/{BUILD}/m03/aligned.fasta",
        metadata=f"results/{BUILD}/m01/metadata_enriched.tsv"
    output:
        tree=f"results/{BUILD}/m04/tree.nwk",
        node_data=f"results/{BUILD}/m04/branch_lengths.json"
    conda:
        "../envs/nextstrain.yaml"
    log:
        f"logs/{BUILD}/m04_augur_refine.log"
    params:
        clock_rate=config["timetree"]["clock_rate"],
        clock_std=config["timetree"]["clock_std_dev"]
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std} \
            --stochastic-resolve 2> {log}
        """

rule augur_ancestral:
    input:
        tree=f"results/{BUILD}/m04/tree.nwk",
        alignment=f"results/{BUILD}/m03/aligned.fasta"
    output:
        node_data=f"results/{BUILD}/m06/nt_muts.json"
    conda:
        "../envs/nextstrain.yaml"
    log:
        f"logs/{BUILD}/m06_augur_ancestral.log"
    shell:
        """
        mkdir -p results/{BUILD}/m06
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference joint 2> {log}
        """

rule augur_translate:
    input:
        tree=f"results/{BUILD}/m04/tree.nwk",
        node_data=f"results/{BUILD}/m06/nt_muts.json",
        reference=config["references"]["annotation"]
    output:
        node_data=f"results/{BUILD}/m06/aa_muts.json"
    conda:
        "../envs/nextstrain.yaml"
    log:
        f"logs/{BUILD}/m06_augur_translate.log"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --validation-mode warn \
            --output-node-data {output.node_data} 2> {log}
        """

rule augur_export:
    input:
        tree=f"results/{BUILD}/m04/tree.nwk",
        metadata=f"results/{BUILD}/m01/metadata_enriched.tsv",
        branch_lengths=f"results/{BUILD}/m04/branch_lengths.json",
        nt_muts=f"results/{BUILD}/m06/nt_muts.json",
        aa_muts=f"results/{BUILD}/m06/aa_muts.json"
    output:
        auspice_json=f"results/{BUILD}/nextstrain/auspice_dataset.json"
    params:
        build_name=BUILD,
        color_by=(
            "clade lineage country region date"
            if config.get("nextclade", {}).get("enabled", True)
            else "subtype country region host date"
        )
    conda:
        "../envs/nextstrain.yaml"
    log:
        f"logs/{BUILD}/m06_augur_export.log"
    shell:
        """
        mkdir -p results/{BUILD}/nextstrain
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} \
            --color-by-metadata {params.color_by} \
            --output {output.auspice_json} \
            --title "Pandemic Preparedness: {params.build_name}" \
            --build-url "https://github.com/engkinandatama/pandemic-preparedness" 2> {log}
        """

# ==========================================================================
# DEFERRED / DROPPED ITEMS for M02-M06
# ==========================================================================
#
# DROPPED: trimAl (M03 — Alignment Trimming)
#   Original plan included trimAl to remove poorly-aligned columns from MSA.
#   Decision: Dropped. augur align (--fill-gaps) performs reference-based
#   alignment which inherently produces cleaner columns than de novo alignment.
#   Additionally, trimming risks removing conserved regions that are critical
#   for primer design. For M05 HyPhy input, clean_alignment_for_hyphy.py
#   already handles codon-level sanitization.
#
# TODO [M04-1] Protein Phylogeny (Antigen/Serological Track)
#   What: Build a separate protein phylogenetic tree from MACSE AA output
#         using IQ-TREE2 with protein substitution models (e.g., WAG, LG).
#   Why deferred: Phase 1 focuses on nucleotide-level diagnostics (primer design).
#                 Protein phylogeny is relevant for serological/antigen-based
#                 applications (epitope prediction, vaccine strain selection)
#                 which are out of scope for the current milestone.
#   Note: Protein/codon data IS used in Phase 1 — but only for HyPhy dN/dS
#         analysis (mutation pressure scoring), not for building protein trees.
#   Priority: Phase 2 — when antigen/serological objectives are in scope.
#
# ==========================================================================
