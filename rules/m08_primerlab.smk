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
        config=f"results/{BUILD}/m08/{{mode}}/primerlab_config_{{region}}.yaml"
    log:
        f"logs/{BUILD}/m08/generate_payload_{{mode}}_{{region}}.log"
    params:
        outdir=f"results/{BUILD}/m08/{{mode}}",
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
        blast_targeted_db_path=config["primerlab"]["design"].get("blast_targeted_db_path"),
        entropy_threshold=config["primerlab"]["design"].get("entropy_threshold"),
        min_exclusion_block_length=config["primerlab"]["design"].get("min_exclusion_block_length", 10)
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
            --outdir {params.outdir} \
            --mode {wildcards.mode}"
        
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
        [[ -n "{params.entropy_threshold}" && "{params.entropy_threshold}" != "None" ]] && CMD="$CMD --entropy-threshold {params.entropy_threshold}"
        [[ -n "{params.min_exclusion_block_length}" && "{params.min_exclusion_block_length}" != "None" ]] && CMD="$CMD --min-exclusion-block-length {params.min_exclusion_block_length}"
        
        eval "$CMD" > "{log}" 2>&1
        """

rule m08_install_primerlab:
    """
    Installs PrimerLab once in editable mode or from remote git repo
    to prevent race conditions when running parallel designs.
    """
    output:
        sentinel=f"results/{BUILD}/m08/.primerlab_installed"
    params:
        repo_path=config["primerlab"]["repo_path"]
    conda:
        "../envs/primerlab.yaml"
    log:
        f"logs/{BUILD}/m08/install_primerlab.log"
    shell:
        """
        #!/usr/bin/env bash
        set -euo pipefail
        (
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
        else
            echo "Skipping installation as repo_path is not specified or set to conda/None."
        fi
        touch "{output.sentinel}"
        ) > "{log}" 2>&1
        """

rule m08_run_primerlab:
    """
    Executes PrimerLab design engine using the generated config.
    """
    input:
        config=f"results/{BUILD}/m08/{{mode}}/primerlab_config_{{region}}.yaml",
        installed=f"results/{BUILD}/m08/.primerlab_installed"
    output:
        report=f"results/{BUILD}/m08/{{mode}}/{{region}}/summary.json",
        alternatives=f"results/{BUILD}/m08/{{mode}}/{{region}}/alternatives.json",
        best_primers=f"results/{BUILD}/m08/{{mode}}/{{region}}/best_primers.json"
    log:
        f"logs/{BUILD}/m08/run_primerlab_{{mode}}_{{region}}.log"
    params:
        outdir=f"results/{BUILD}/m08/{{mode}}",
        validate=config["primerlab"]["design"].get("validate", False)
    conda:
        "../envs/primerlab.yaml"
    shell:
        """
        #!/usr/bin/env bash
        set -euo pipefail
        (
        # Build validation flag
        VALIDATE_FLAG=""
        if [ "{params.validate}" = "True" ] || [ "{params.validate}" = "true" ]; then
            VALIDATE_FLAG="--validate"
        fi
        
        # Run design
        primerlab raa --config "{input.config}" --out "{params.outdir}" $VALIDATE_FLAG
        ) > "{log}" 2>&1
        """

rule m08_degenerate_rescue:
    """
    Rescues degenerate positions in design candidates based on alignment data.
    Runs only for degenerate mode.
    """
    input:
        manifest=f"results/{BUILD}/m07/primer_targets_manifest.yaml",
        summary=f"results/{BUILD}/m08/degenerate/{{region}}/summary.json",
        alternatives=f"results/{BUILD}/m08/degenerate/{{region}}/alternatives.json",
        best_primers=f"results/{BUILD}/m08/degenerate/{{region}}/best_primers.json"
    output:
        marker=f"results/{BUILD}/m08/degenerate/{{region}}/.degenerate_rescued"
    log:
        f"logs/{BUILD}/m08/degenerate_rescue_{{region}}.log"
    params:
        outdir=f"results/{BUILD}/m08/degenerate/{{region}}",
        maf_threshold=config["primerlab"].get("design", {}).get("degenerate", {}).get("maf_threshold", 0.10),
        max_multiplier=config["primerlab"].get("design", {}).get("degenerate", {}).get("max_multiplier", 16),
        thermo_mode=config["primerlab"].get("design", {}).get("degenerate", {}).get("thermo_eval_mode", "worst_case"),
        vienna_temp=config["primerlab"].get("design", {}).get("degenerate", {}).get("vienna_temp", 37.0),
        # QC gate thresholds — defaults match the PrimerLab engine templates so the
        # degenerate track can never be more permissive than the non-degenerate one.
        hairpin_dg_min=config["primerlab"].get("design", {}).get("qc", {}).get("hairpin_dg_min", -2.0),
        homodimer_dg_min=config["primerlab"].get("design", {}).get("qc", {}).get("homodimer_dg_min", -5.0),
        heterodimer_dg_min=config["primerlab"].get("design", {}).get("qc", {}).get("heterodimer_dg_min", -6.0),
        tm_diff_max=config["primerlab"].get("design", {}).get("qc", {}).get("tm_diff_max", 3.0),
        strict_qc="--strict-qc" if config["primerlab"].get("design", {}).get("qc", {}).get("strict", False) else ""
    conda:
        "../envs/primerlab.yaml"
    shell:
        """
        python scripts/m08_degenerate_rescue.py \
            --manifest {input.manifest} \
            --summary {input.summary} \
            --alternatives {input.alternatives} \
            --best-primers {input.best_primers} \
            --outdir {params.outdir} \
            --maf-threshold {params.maf_threshold} \
            --max-multiplier {params.max_multiplier} \
            --thermo-mode {params.thermo_mode} \
            --hairpin-dg-min {params.hairpin_dg_min} \
            --homodimer-dg-min {params.homodimer_dg_min} \
            --heterodimer-dg-min {params.heterodimer_dg_min} \
            --tm-diff-max {params.tm_diff_max} \
            {params.strict_qc} \
            --vienna-temp {params.vienna_temp} > "{log}" 2>&1
        
        touch {output.marker}
        """

def get_m08_finalize_inputs(wildcards):
    if wildcards.mode == "degenerate":
        return f"results/{BUILD}/m08/degenerate/{wildcards.region}/.degenerate_rescued"
    else:
        return f"results/{BUILD}/m08/non_degenerate/{wildcards.region}/summary.json"

rule m08_finalize_design:
    """
    Finalizes the design step for the track.
    Combines or routes degenerate rescue outputs appropriately.
    """
    input:
        get_m08_finalize_inputs
    output:
        marker=f"results/{BUILD}/m08/{{mode}}/{{region}}/.design_complete"
    shell:
        "touch {output.marker}"
