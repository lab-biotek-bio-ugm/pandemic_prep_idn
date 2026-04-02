"""
Module 2-6: Augur Core (Nextstrain Ecosystem)
Executes the Nextstrain spatio-temporal filtering, multiple sequence alignment, 
phylogenetic tree inference, mutational trait analysis, and JSON payload export.
"""

# Import the global BUILD variable defined in the main Snakefile
BUILD = config["build_name"]

rule augur_filter:
    input:
        sequences=f"results/{BUILD}/m01/raw_sequences.fasta",
        metadata=f"results/{BUILD}/m01/metadata.tsv"
    output:
        sequences=f"results/{BUILD}/m02/filtered_sequences.fasta"
    params:
        group_by=config["subsampling"]["group_by"],
        seq_per_group=config["subsampling"]["sequences_per_group"],
        max_seq=config["subsampling"]["max_sequences"],
        min_length=config["ncbi_query"]["min_length"]
    conda:
        "../envs/nextstrain.yaml"
    log:
        f"logs/{BUILD}/m02_augur_filter.log"
    shell:
        """
        mkdir -p results/{BUILD}/m02
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --group-by {params.group_by} \
            --sequences-per-group {params.seq_per_group} \
            --subsample-max-sequences {params.max_seq} \
            --min-length {params.min_length} \
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
        metadata=f"results/{BUILD}/m01/metadata.tsv"
    output:
        tree=f"results/{BUILD}/m04/tree.nwk",
        node_data=f"results/{BUILD}/m04/branch_lengths.json"
    conda:
        "../envs/nextstrain.yaml"
    log:
        f"logs/{BUILD}/m04_augur_refine.log"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --clock-rate 0.0008 \
            --clock-std-dev 0.0004 2> {log}
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
            --output-node-data {output.node_data} 2> {log}
        """

rule augur_export:
    input:
        tree=f"results/{BUILD}/m04/tree.nwk",
        metadata=f"results/{BUILD}/m01/metadata.tsv",
        branch_lengths=f"results/{BUILD}/m04/branch_lengths.json",
        nt_muts=f"results/{BUILD}/m06/nt_muts.json",
        aa_muts=f"results/{BUILD}/m06/aa_muts.json"
    output:
        auspice_json=f"results/{BUILD}/nextstrain/auspice_dataset.json"
    params:
        build_name=BUILD
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
            --output {output.auspice_json} \
            --title "Pandemic Preparedness: {params.build_name}" \
            --build-url "https://github.com/engkinandatama/pandemic-preparedness" 2> {log}
        """
