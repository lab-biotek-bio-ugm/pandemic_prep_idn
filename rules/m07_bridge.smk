"""
Module 7: The Bridge — Conservation & Selection Pressure Integrator
Merges MSA conservation (Track A) with HyPhy evolutionary data (Track B)
to identify and rank Gold Zones for primer design.
"""

BUILD = config["build_name"]
_HYPHY_METHODS = config["hyphy"].get("methods", ["FUBAR"])

# Build GARD input list conditionally
_GARD_INPUTS = ([f"results/{BUILD}/m05/gard_{gene}.json" for gene in TARGET_GENES]
                if "GARD" in _HYPHY_METHODS else [])

rule m07_export_config:
    """Exports the merged Snakemake configuration to a YAML file for python scripts."""
    output:
        yaml=f"results/{BUILD}/m07/active_config.yaml"
    run:
        import yaml
        import os
        os.makedirs(os.path.dirname(output.yaml), exist_ok=True)
        with open(output.yaml, "w") as f:
            yaml.dump(config, f, default_flow_style=False)

rule m07_bridge:
    """Runs the Bridge Builder to produce ranked primer target zones.
    Reads MSA + HyPhy outputs → produces BED, JSON report, YAML manifest.
    GARD inputs are conditional — if GARD is not in hyphy.methods,
    M07 still runs using only FUBAR + MEME (recombination penalty = 0).
    """
    input:
        msa=f"results/{BUILD}/m03/aligned.fasta",
        gff=config["references"]["annotation"],
        reference=config["references"]["genome"],
        config_yaml=f"results/{BUILD}/m07/active_config.yaml",
        metadata=f"results/{BUILD}/m01/metadata_enriched.tsv",
        fubar=[f"results/{BUILD}/m05/fubar_{gene}.json" for gene in TARGET_GENES],
        meme=[f"results/{BUILD}/m05/meme_{gene}.json" for gene in TARGET_GENES],
        gard=_GARD_INPUTS,
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
            --metadata {input.metadata} \
            --outdir results/{params.build_name}/m07 \
            --config {input.config_yaml} \
            --build-name {params.build_name} 2>&1 | tee {log}
        """
