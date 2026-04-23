"""
Module 5: Evolutionary Selection Pressure Engine
Subsamples the M03 MSA for diversity, slices coding genes via GFF coordinates,
applies harsh masking, then runs HyPhy (FUBAR/MEME/GARD) for dN/dS analysis.
"""

BUILD = config["build_name"]
HYPHY_METHODS = config["hyphy"].get("methods", ["FUBAR"])

_SUB_CFG = config["hyphy"].get("subsample", {})
_SUB_ENABLED = _SUB_CFG.get("enabled", False)
_SUB_GRP = _SUB_CFG.get("group_by", {})

rule m05_hyphy_subsample:
    """Diversity-aware subsampling of the global MSA for Track B (HyPhy).
    Runs ONCE before gene extraction — all genes share the same consistent subset.
    Adaptive grouping: clade x country x year -> clade x country -> clade -> random.
    Track A (M07 Bridge entropy scan) still uses the full aligned.fasta.
    If subsample.enabled is false, copies the full MSA as-is (bypass).
    """
    input:
        alignment=f"results/{BUILD}/m03/aligned.fasta",
        metadata=f"results/{BUILD}/m01/metadata_enriched.tsv"
    output:
        alignment=f"results/{BUILD}/m05/aligned_hyphy.fasta",
        report=f"results/{BUILD}/m05/hyphy_subsample_report.json"
    params:
        enabled=_SUB_ENABLED,
        max_seq=_SUB_CFG.get("max_sequences", 300),
        min_seq=_SUB_CFG.get("min_sequences", 50),
        strategy=_SUB_CFG.get("strategy", "balanced"),
        seed=_SUB_CFG.get("seed", 42),
        subtype_flag="" if _SUB_GRP.get("subtype", True) else "--no-group-by-subtype",
        clade_flag="" if _SUB_GRP.get("clade", True) else "--no-group-by-clade",
        country_flag="" if _SUB_GRP.get("country", True) else "--no-group-by-country",
        year_flag="" if _SUB_GRP.get("year", True) else "--no-group-by-year"
    conda:
        "../envs/hyphy_macse.yaml"
    log:
        f"logs/{BUILD}/m05_hyphy_subsample.log"
    shell:
        """
        if [ "{params.enabled}" = "True" ]; then
            python scripts/hyphy_subsample.py \
                --alignment {input.alignment} \
                --metadata {input.metadata} \
                --output {output.alignment} \
                --report {output.report} \
                --max-sequences {params.max_seq} \
                --min-sequences {params.min_seq} \
                --strategy {params.strategy} \
                --seed {params.seed} \
                {params.subtype_flag} \
                {params.clade_flag} \
                {params.country_flag} \
                {params.year_flag} 2>&1 | tee {log}
        else
            echo '[SUBSAMPLE] DISABLED — copying full MSA' | tee {log}
            cp {input.alignment} {output.alignment}
            echo '{{"bypass": true, "reason": "subsample.enabled=false"}}' > {output.report}
        fi
        """

rule m05_extract_genes:
    """Slices genes from the subsampled MSA and applies Harsh Masking for reading-frame integrity."""
    input:
        alignment=f"results/{BUILD}/m05/aligned_hyphy.fasta",
        gff=config["references"]["annotation"]
    output:
        fastas=[f"results/{BUILD}/m05/cleaned_{gene}_NT.fasta" for gene in TARGET_GENES]
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

rule m05_prune_tree:
    """Prunes the full phylogenetic tree to only contain tips present in the gene alignment."""
    input:
        tree=f"results/{BUILD}/m04/tree.nwk",
        alignment=f"results/{BUILD}/m05/cleaned_{{gene}}_NT.fasta"
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
        config.get("hyphy", {}).get("threads", 32)
    conda:
        "../envs/hyphy_macse.yaml"
    log:
        f"logs/{BUILD}/m05_hyphy_fubar_{{gene}}.log"
    shell:
        """
        hyphy CPU={threads} fubar \
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
        config.get("hyphy", {}).get("threads", 32)
    conda:
        "../envs/hyphy_macse.yaml"
    log:
        f"logs/{BUILD}/m05_hyphy_meme_{{gene}}.log"
    shell:
        """
        hyphy CPU={threads} meme \
            --alignment {input.alignment} \
            --tree {input.tree} \
            --output {output.json} 2>&1 | tee {log}
        """

rule m05_hyphy_gard:
    """Runs HyPhy GARD to detect recombination breakpoints within each gene.
    Uses a safe wrapper that performs pre-flight checks:
    - If alignment is too short for c-AIC (alignment_len < ~4xnum_seqs),
      writes empty JSON and exits 0 (soft-fail, pipeline continues).
    - If GARD crashes at runtime, catches error and soft-fails.
    Critical for long genomes (SARS-CoV-2); auto-skipped for short segments (Flu).
    """
    input:
        alignment=f"results/{BUILD}/m05/cleaned_{{gene}}_NT.fasta"
    output:
        json=f"results/{BUILD}/m05/gard_{{gene}}.json"
    threads:
        config.get("hyphy", {}).get("threads", 32)
    conda:
        "../envs/hyphy_macse.yaml"
    log:
        f"logs/{BUILD}/m05_hyphy_gard_{{gene}}.log"
    shell:
        """
        python scripts/run_hyphy_gard_safe.py \
            {input.alignment} \
            {output.json} \
            --threads {threads} \
            --rv GDD 2>&1 | tee {log}
        """

# DEFERRED / DROPPED ITEMS for M05
#
# TODO [M05-1] Clade-Defining SNP Finder (linked to DEV-1 in Snakefile)
#   What: Script to identify mutations unique to each clade from nt_muts.json
#         combined with clade labels from Nextclade metadata. Output: clade_snps.tsv
#         containing per-position SNPs that distinguish clades from each other.
#   Why deferred: Only needed for DEV-1 (Route B clade-specific primer targets).
#                 Primary milestone is universal primer detection (Route A).
#   Priority: Phase 2 — together with DEV-1 and DEV-4 (consensus per clade).
#
# DROPPED: HyPhy FEL (Fixed Effects Likelihood)
#   Original plan included FEL alongside FUBAR for pervasive selection detection.
#   Decision: Dropped. FUBAR (Bayesian) and FEL (Frequentist) detect the same
#   biological signal (pervasive selection). FUBAR is faster and performs equally
#   well for datasets with >20 sequences (our typical range is 50-1000).
#   FEL is only advantageous for very small datasets (<20 sequences).
#
# TODO [M05-2] SNP Matrix (linked to DEV-1 and M05-1)
#   What: Generate a position x sequence allele matrix from the MSA. Each row
#         is a genome position, each column is a sequence, each cell is the
#         nucleotide at that position. This is the raw data the Clade SNP
#         Finder (M05-1) needs to identify clade-defining mutations.
#   Why deferred: Only needed when Route B (clade-specific primers) is active.
#                 Note: augur ancestral (nt_muts.json) gives mutations per tree
#                 branch, NOT per-position alleles — they are complementary,
#                 not redundant.
#   Priority: Phase 2 — together with M05-1 and DEV-4.
#
