"""
Module 7: The Bridge — Conservation & Selection Pressure Integrator
Merges MSA conservation (Track A) with HyPhy evolutionary data (Track B)
to identify and rank Gold Zones for primer design.
"""

BUILD = config["build_name"]

rule m07_bridge:
    """Runs the Bridge Builder to produce ranked primer target zones.
    Reads MSA + HyPhy outputs → produces BED, JSON report, YAML manifest.
    """
    input:
        msa=f"results/{BUILD}/m03/aligned.fasta",
        gff=config["references"]["annotation"],
        reference=config["references"]["genome"],
        # Ensure all HyPhy outputs are available before running
        fubar=[f"results/{BUILD}/m05/fubar_{gene}.json" for gene in TARGET_GENES],
        meme=[f"results/{BUILD}/m05/meme_{gene}.json" for gene in TARGET_GENES],
        gard=[f"results/{BUILD}/m05/gard_{gene}.json" for gene in TARGET_GENES],
    output:
        bed=f"results/{BUILD}/m07/conserved_regions.bed",
        report=f"results/{BUILD}/m07/bridge_report.json",
        manifest=f"results/{BUILD}/m07/primer_targets_manifest.yaml",
    params:
        hyphy_dir=f"results/{BUILD}/m05",
        build_name=BUILD,
    conda:
        "../envs/hyphy_macse.yaml"
    log:
        f"logs/{BUILD}/m07_bridge.log"
    shell:
        """
        python scripts/m07_bridge_builder.py \
            --msa {input.msa} \
            --gff {input.gff} \
            --hyphy-dir {params.hyphy_dir} \
            --reference {input.reference} \
            --outdir results/{params.build_name}/m07 \
            --config config/config.yaml \
            --build-name {params.build_name} 2>&1 | tee {log}
        """
