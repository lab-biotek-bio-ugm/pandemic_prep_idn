"""
Module 9: BLAST Specificity Verification
Runs detailed NCBI or local BLAST specificity queries for top N candidate primer sets.
"""

BUILD = config["build_name"]

rule m09_blast_verification:
    """
    Runs BLAST off-target queries for top N candidates in alternatives.json
    and exports dedicated specificity_report.md files.
    """
    input:
        alternatives=f"results/{BUILD}/m08/Whole_Segment_region_1/summary.json" # Wait for m08 completion
    output:
        marker=f"results/{BUILD}/m09/verification_complete.txt"
    log:
        f"logs/{BUILD}/m09/blast_verification.log"
    params:
        alternatives_file=f"results/{BUILD}/m08/Whole_Segment_region_1/alternatives.json",
        num_candidates=config["blast_verification"].get("num_candidates", 3),
        database=config["blast_verification"].get("database", "nt"),
        mode=config["blast_verification"].get("mode", "remote"),
        local_db_path=config["blast_verification"].get("local_db_path"),
        outdir=f"results/{BUILD}/m09"
    conda:
        "../envs/primerlab.yaml"
    shell:
        """
        #!/usr/bin/env bash
        set -euo pipefail
        
        DB_TARGET="{params.database}"
        if [ "{params.mode}" = "local" ] && [ -n "{params.local_db_path}" ] && [ "{params.local_db_path}" != "None" ]; then
            DB_TARGET="{params.local_db_path}"
        fi
        
        python scripts/m09_blast_verification.py \
            --alternatives {params.alternatives_file} \
            --num-candidates {params.num_candidates} \
            --database "$DB_TARGET" \
            --mode {params.mode} \
            --outdir {params.outdir} > {log} 2>&1
        """
