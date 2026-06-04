"""
Module 10 (supporting): Custom BLAST Database Builder
Builds human_db and virus_db locally for M09 specificity verification.
Only runs if custom_db_build.enabled = true and mode = local.
"""

import os

BUILD = config["build_name"]

# ─── Helper ──────────────────────────────────────────────────────────────────
def _custom_db_enabled():
    """Return True if local custom DB build is configured and enabled."""
    bv = config.get("blast_verification", {})
    return (
        bv.get("mode", "remote") == "local"
        and bv.get("custom_db_build", {}).get("enabled", False)
    )

# ─── Rule ─────────────────────────────────────────────────────────────────────
if _custom_db_enabled():
    _DB_BUILD = config["blast_verification"]["custom_db_build"]
    _OUTPUT_DIR = _DB_BUILD.get("output_dir", "resources/blast_db")

    rule build_custom_blast_db:
        """
        Build a custom local BLAST database from host genomes (GRCh38, T2T-CHM13)
        and respiratory virus sequences (NCBI Datasets CLI).

        Database is stored in resources/blast_db/ at the project root, persisted
        across all Snakemake runs. Hash-based caching skips rebuild if config
        has not changed.
        """
        output:
            marker=f"{_OUTPUT_DIR}/build_complete.txt"
        params:
            config_file=lambda wildcards, input: os.path.abspath(workflow.configfiles[0])
                        if workflow.configfiles else "config/config.yaml"
        log:
            "logs/build_db/build_custom_blast_db.log"
        conda:
            "../envs/blast_db.yaml"
        shell:
            """
            #!/usr/bin/env bash
            set -euo pipefail

            echo "Building custom BLAST database..."
            echo "Config file: {params.config_file}"

            python -u scripts/build_custom_blast_db.py \
                --config-file "{params.config_file}" \
                > {log} 2>&1

            echo "Build complete marker: {output.marker}"
            """
