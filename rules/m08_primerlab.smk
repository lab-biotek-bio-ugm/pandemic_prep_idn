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
    log:
        f"logs/{BUILD}/m08/generate_payload_Whole_Segment_region_1.log"
    params:
        outdir=f"results/{BUILD}/m08",
        num_candidates=config["primerlab"]["design"].get("num_candidates"),
        vienna_limit=config["primerlab"]["design"].get("vienna_ranking_limit"),
        window_size=config["primerlab"]["design"].get("window_size"),
        overlap=config["primerlab"]["design"].get("overlap"),
        search_strategy=config["primerlab"]["design"].get("search_strategy"),
        hard_qc_filter=config["primerlab"]["design"].get("hard_qc_filter"),
        max_search_iterations=config["primerlab"]["design"].get("max_search_iterations"),
        validate=config["primerlab"]["design"].get("validate"),
        blast_check=config["primerlab"]["design"].get("blast_offtarget_check"),
        blast_limit=config["primerlab"]["design"].get("blast_ranking_limit"),
        blast_db=config["primerlab"]["design"].get("blast_database")
    conda:
        "../envs/primerlab.yaml"
    shell:
        """
        #!/usr/bin/env bash
        set -euo pipefail
        
        # Build command with optional arguments
        CMD="python scripts/m08_payload_generator.py \
            --manifest {input.manifest} \
            --template {input.template} \
            --outdir {params.outdir}"
        
        [[ -n "{params.num_candidates}" ]] && CMD="$CMD --num-candidates {params.num_candidates}"
        [[ -n "{params.vienna_limit}" ]] && CMD="$CMD --vienna-limit {params.vienna_limit}"
        [[ -n "{params.window_size}" ]] && CMD="$CMD --window-size {params.window_size}"
        [[ -n "{params.overlap}" ]] && CMD="$CMD --overlap {params.overlap}"
        [[ -n "{params.search_strategy}" ]] && CMD="$CMD --search-strategy {params.search_strategy}"
        [[ -n "{params.hard_qc_filter}" ]] && CMD="$CMD --hard-qc-filter {params.hard_qc_filter}"
        [[ -n "{params.max_search_iterations}" ]] && CMD="$CMD --max-search-iterations {params.max_search_iterations}"
        [[ -n "{params.validate}" ]] && CMD="$CMD --validate {params.validate}"
        [[ -n "{params.blast_check}" ]] && CMD="$CMD --blast-check {params.blast_check}"
        [[ -n "{params.blast_limit}" ]] && CMD="$CMD --blast-limit {params.blast_limit}"
        [[ -n "{params.blast_db}" ]] && CMD="$CMD --blast-db {params.blast_db}"
        
        eval $CMD > {log} 2>&1
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
    log:
        f"logs/{BUILD}/m08/run_primerlab_Whole_Segment_region_1.log"
    params:
        repo_path=config["primerlab"]["repo_path"],
        validate=config["primerlab"]["design"].get("validate", False)
    conda:
        "../envs/primerlab.yaml"
    shell:
        """
        #!/usr/bin/env bash
        set -euo pipefail
        (
        # Install primerlab if repo_path is provided
        if [ -n "{params.repo_path}" ] && [ "{params.repo_path}" != "conda" ] && [ "{params.repo_path}" != "None" ]; then
            if [[ "{params.repo_path}" == git+* ]]; then
                pip install "{params.repo_path}"
            else
                pip install -e "{params.repo_path}"
            fi
        fi
        
        # Build validation flag
        VALIDATE_FLAG=""
        if [ "{params.validate}" = "True" ] || [ "{params.validate}" = "true" ]; then
            VALIDATE_FLAG="--validate"
        fi
        
        # Run design
        primerlab raa --config {input.config} --output-dir results/{BUILD}/m08 $VALIDATE_FLAG
        ) > {log} 2>&1
        """
