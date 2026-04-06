"""
Module 5: Evolurioanry Selection Pressure Engine
Uses MACSE for codon-aware alignment and HyPhy for dN/dS selection analysis.
"""

BUILD = config["build_name"]
HYPHY_METHODS = config["hyphy"].get("methods", ["FUBAR"])

rule m05_extract_genes:
    """Extracts unaligned gene nucleotide sequences based on GFF coordinates inside the MSA."""
    input:
        alignment=f"results/{BUILD}/m03/aligned.fasta",
        gff=config["references"]["annotation"],
        # Provide the reference ID used in the MSA (typically MN908947.3 for SARS-CoV-2)
        # In augur align, the reference is usually injected at the start with its original ID.
    output:
        # Dynamically expect files for all target genes defined in config
        fastas=[f"results/{BUILD}/m05/unaligned_{gene}.fasta" for gene in TARGET_GENES]
    params:
        outdir=f"results/{BUILD}/m05",
        # We assume the reference ID contains 'MN908947' or similar based on typical nextstrain
        ref_id="MN908947" 
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
            
        # Rename outputs to match expected wildcard template
        for gene in {TARGET_GENES}; do
            if [ -f "{params.outdir}/$gene.fasta" ]; then
                mv "{params.outdir}/$gene.fasta" "{params.outdir}/unaligned_$gene.fasta"
            fi
        done
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
        # MACSE v2 is called via java wrapper or conda command 'macse'
        macse -prog alignSequences \
            -seq {input.fasta} \
            -out_NT {output.nt_fasta} \
            -out_AA {output.aa_fasta} 2>&1 | tee {log}
        """

rule m05_hyphy_fubar:
    """Runs HyPhy FUBAR to detect pervasive positive and negative selection."""
    input:
        alignment=f"results/{BUILD}/m05/macse_{{gene}}_NT.fasta",
        tree=f"results/{BUILD}/m04/tree.nwk"
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
        alignment=f"results/{BUILD}/m05/macse_{{gene}}_NT.fasta",
        tree=f"results/{BUILD}/m04/tree.nwk"
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
