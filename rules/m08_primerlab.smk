"""
Module 8: PrimerLab Integration
Generates PrimerLab payloads from M07 manifest and executes the design engine.
"""

BUILD = config["build_name"]

rule m08_generate_payload:
    """
    Converts M07 manifest into PrimerLab-compatible YAML configuration.
    Injects consensus sequence and MEME-based excluded regions.
    Also overrides design parameters from global config.
    """
    input:
        manifest=f"results/{BUILD}/m07/primer_targets_manifest.yaml",
        template=config["primerlab"]["template"]
    output:
        config=f"results/{BUILD}/m08/primerlab_config_Whole_Segment_region_1.yaml"
    params:
        outdir=f"results/{BUILD}/m08",
        num_candidates=config["primerlab"]["design"].get("num_candidates"),
        vienna_limit=config["primerlab"]["design"].get("vienna_ranking_limit"),
        window_size=config["primerlab"]["design"].get("window_size"),
        overlap=config["primerlab"]["design"].get("overlap")
    conda:
        "../envs/primerlab.yaml"
    shell:
        """
        # Build command with optional arguments
        CMD="python scripts/m08_payload_generator.py \
            --manifest {input.manifest} \
            --template {input.template} \
            --outdir {params.outdir}"
        
        [[ -n "{params.num_candidates}" ]] && CMD="$CMD --num-candidates {params.num_candidates}"
        [[ -n "{params.vienna_limit}" ]] && CMD="$CMD --vienna-limit {params.vienna_limit}"
        [[ -n "{params.window_size}" ]] && CMD="$CMD --window-size {params.window_size}"
        [[ -n "{params.overlap}" ]] && CMD="$CMD --overlap {params.overlap}"
        
        eval $CMD
        """

rule m08_run_primerlab:
    """
    Executes PrimerLab design engine using the generated config.
    Installs primerlab in editable mode inside the environment if needed.
    """
    input:
        config=f"results/{BUILD}/m08/primerlab_config_Whole_Segment_region_1.yaml"
    output:
        report=f"results/{BUILD}/m08/Whole_Segment_region_1/raa_result.json"
    params:
        repo_path=config["primerlab"]["repo_path"]
    conda:
        "../envs/primerlab.yaml"
    shell:
        """
        # Ensure primerlab is installed in this environment
        pip install -e {params.repo_path}
        
        # Run design
        primerlab raa --config {input.config} --output-dir results/{BUILD}/m08
        """
