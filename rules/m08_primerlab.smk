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
        blast_db=config["primerlab"]["design"].get("blast_database"),
        blast_mode=config["primerlab"]["design"].get("blast_mode"),
        blast_local_db_path=config["primerlab"]["design"].get("blast_local_db_path"),
        blast_targeted_check=config["primerlab"]["design"].get("blast_targeted_check"),
        blast_targeted_db_path=config["primerlab"]["design"].get("blast_targeted_db_path")
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
        
        [[ -n "{params.num_candidates}" && "{params.num_candidates}" != "None" ]] && CMD="$CMD --num-candidates {params.num_candidates}"
        [[ -n "{params.vienna_limit}" && "{params.vienna_limit}" != "None" ]] && CMD="$CMD --vienna-limit {params.vienna_limit}"
        [[ -n "{params.window_size}" && "{params.window_size}" != "None" ]] && CMD="$CMD --window-size {params.window_size}"
        [[ -n "{params.overlap}" && "{params.overlap}" != "None" ]] && CMD="$CMD --overlap {params.overlap}"
        [[ -n "{params.search_strategy}" && "{params.search_strategy}" != "None" ]] && CMD="$CMD --search-strategy {params.search_strategy}"
        [[ -n "{params.hard_qc_filter}" && "{params.hard_qc_filter}" != "None" ]] && CMD="$CMD --hard-qc-filter {params.hard_qc_filter}"
        [[ -n "{params.max_search_iterations}" && "{params.max_search_iterations}" != "None" ]] && CMD="$CMD --max-search-iterations {params.max_search_iterations}"
        [[ -n "{params.validate}" && "{params.validate}" != "None" ]] && CMD="$CMD --validate {params.validate}"
        [[ -n "{params.blast_check}" && "{params.blast_check}" != "None" ]] && CMD="$CMD --blast-check {params.blast_check}"
        [[ -n "{params.blast_limit}" && "{params.blast_limit}" != "None" ]] && CMD="$CMD --blast-limit {params.blast_limit}"
        [[ -n "{params.blast_db}" && "{params.blast_db}" != "None" ]] && CMD="$CMD --blast-db {params.blast_db}"
        [[ -n "{params.blast_mode}" && "{params.blast_mode}" != "None" ]] && CMD="$CMD --blast-mode {params.blast_mode}"
        [[ -n "{params.blast_local_db_path}" && "{params.blast_local_db_path}" != "None" ]] && CMD="$CMD --blast-local-db-path {params.blast_local_db_path}"
        [[ -n "{params.blast_targeted_check}" && "{params.blast_targeted_check}" != "None" ]] && CMD="$CMD --blast-targeted-check {params.blast_targeted_check}"
        [[ -n "{params.blast_targeted_db_path}" && "{params.blast_targeted_db_path}" != "None" ]] && CMD="$CMD --blast-targeted-db-path {params.blast_targeted_db_path}"
        
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
        report=f"results/{BUILD}/m08/Whole_Segment_region_1/summary.json"
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
        REPO_PATH="{params.repo_path}"
        if [ -n "$REPO_PATH" ] && [ "$REPO_PATH" != "conda" ] && [ "$REPO_PATH" != "None" ]; then
            if [[ "$REPO_PATH" == git+* ]]; then
                pip install --upgrade --force-reinstall --no-cache-dir "$REPO_PATH"
            else
                # If running on Linux and path has Windows drive letter (e.g., e:/Project/...)
                if [[ "$REPO_PATH" =~ ^[a-zA-Z]:/ ]]; then
                    DIR_NAME=$(basename "$REPO_PATH")
                    PARENT_DIR=$(dirname "$(pwd)")
                    LOCAL_LINUX_PATH="$PARENT_DIR/$DIR_NAME"
                    if [ -d "$LOCAL_LINUX_PATH" ]; then
                        REPO_PATH="$LOCAL_LINUX_PATH"
                    else
                        echo "Warning: Windows path $REPO_PATH not found, and fallback $LOCAL_LINUX_PATH does not exist."
                    fi
                fi
                
                if [ -d "$REPO_PATH" ]; then
                    pip install -e "$REPO_PATH"
                else
                    echo "Skipping editable install because directory does not exist: $REPO_PATH"
                fi
            fi
        fi
        
        # Build validation flag
        VALIDATE_FLAG=""
        if [ "{params.validate}" = "True" ] || [ "{params.validate}" = "true" ]; then
            VALIDATE_FLAG="--validate"
        fi
        
        # Run design
        primerlab raa --config {input.config} --out results/{BUILD}/m08 $VALIDATE_FLAG
        ) > {log} 2>&1
        """
