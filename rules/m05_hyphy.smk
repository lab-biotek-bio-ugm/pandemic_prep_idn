"""
Module 5: Evolutionary Selection Pressure Engine
Uses MACSE for codon-aware alignment and HyPhy for dN/dS selection analysis.
"""

BUILD = config["build_name"]
HYPHY_METHODS = config["hyphy"].get("methods", ["FUBAR"])

rule m05_extract_genes:
    """Extracts unaligned gene nucleotide sequences based on GFF coordinates inside the MSA."""
    input:
        alignment=f"results/{BUILD}/m03/aligned.fasta",
        gff=config["references"]["annotation"]
    output:
        fastas=[f"results/{BUILD}/m05/unaligned_{gene}.fasta" for gene in TARGET_GENES]
    params:
        outdir=f"results/{BUILD}/m05",
        ref_id=config["references"].get("reference_id", "MN908947")
    conda:
        "../envs/hyphy_macse.yaml"
    log:
        f"logs/{BUILD}/m05_extract_genes.log"
    shell:
        """
        python scripts/extract_genes.py \
            --alignment {input.alignment} \
            --gff {input.gff} \
            --reference-id {params.ref_id} \
            --outdir {params.outdir} 2>&1 | tee {log}
        """

rule m05_macse_align:
    """Uses MACSE to perform strict codon-aware alignment of the extracted gene."""
    input:
        fasta=f"results/{BUILD}/m05/unaligned_{{gene}}.fasta"
    output:
        nt_fasta=f"results/{BUILD}/m05/macse_{{gene}}_NT.fasta",
        aa_fasta=f"results/{BUILD}/m05/macse_{{gene}}_AA.fasta"
    threads:
        4
    conda:
        "../envs/hyphy_macse.yaml"
    log:
        f"logs/{BUILD}/m05_macse_{{gene}}.log"
    shell:
        """
        macse -prog alignSequences \
            -seq {input.fasta} \
            -out_NT {output.nt_fasta} \
            -out_AA {output.aa_fasta} 2>&1 | tee {log}
        """

rule m05_clean_alignment:
    """Sanitizes MACSE outputs by masking stop codons and frameshifts so HyPhy doesn't crash."""
    input:
        alignment=f"results/{BUILD}/m05/macse_{{gene}}_NT.fasta"
    output:
        cleaned=f"results/{BUILD}/m05/cleaned_{{gene}}_NT.fasta"
    conda:
        "../envs/hyphy_macse.yaml"
    log:
        f"logs/{BUILD}/m05_clean_{{gene}}.log"
    shell:
        """
        python scripts/clean_alignment_for_hyphy.py {input.alignment} {output.cleaned} 2>&1 | tee {log}
        """

rule m05_prune_tree:
    """Prunes the full phylogenetic tree to only contain tips present in the gene alignment."""
    input:
        tree=f"results/{BUILD}/m04/tree.nwk",
        alignment=f"results/{BUILD}/m05/macse_{{gene}}_NT.fasta"
    output:
        tree=f"results/{BUILD}/m05/tree_{{gene}}.nwk"
    conda:
        "../envs/hyphy_macse.yaml"
    log:
        f"logs/{BUILD}/m05_prune_tree_{{gene}}.log"
    script:
        "../scripts/prune_tree.py"

rule m05_hyphy_fubar:
    """Runs HyPhy FUBAR to detect pervasive positive and negative selection."""
    input:
        alignment=f"results/{BUILD}/m05/cleaned_{{gene}}_NT.fasta",
        tree=f"results/{BUILD}/m05/tree_{{gene}}.nwk"
    output:
        json=f"results/{BUILD}/m05/fubar_{{gene}}.json"
    threads:
        4
    conda:
        "../envs/hyphy_macse.yaml"
    log:
        f"logs/{BUILD}/m05_hyphy_fubar_{{gene}}.log"
    shell:
        """
        hyphy fubar \
            --alignment {input.alignment} \
            --tree {input.tree} \
            --output {output.json} 2>&1 | tee {log}
        """

rule m05_hyphy_meme:
    """Runs HyPhy MEME to detect episodic positive selection."""
    input:
        alignment=f"results/{BUILD}/m05/cleaned_{{gene}}_NT.fasta",
        tree=f"results/{BUILD}/m05/tree_{{gene}}.nwk"
    output:
        json=f"results/{BUILD}/m05/meme_{{gene}}.json"
    threads:
        4
    conda:
        "../envs/hyphy_macse.yaml"
    log:
        f"logs/{BUILD}/m05_hyphy_meme_{{gene}}.log"
    shell:
        """
        hyphy meme \
            --alignment {input.alignment} \
            --tree {input.tree} \
            --output {output.json} 2>&1 | tee {log}
        """

rule m05_hyphy_gard:
    """Runs HyPhy GARD to detect recombination breakpoints within each gene.
    Critical for Coronaviruses; less so for segmented viruses (Influenza).
    GARD infers its own phylogenies per partition, so no input tree is needed.
    """
    input:
        alignment=f"results/{BUILD}/m05/cleaned_{{gene}}_NT.fasta"
    output:
        json=f"results/{BUILD}/m05/gard_{{gene}}.json"
    threads:
        4
    conda:
        "../envs/hyphy_macse.yaml"
    log:
        f"logs/{BUILD}/m05_hyphy_gard_{{gene}}.log"
    shell:
        """
        hyphy gard \
            --alignment {input.alignment} \
            --output {output.json} \
            --rv GDD 2>&1 | tee {log}
        """
