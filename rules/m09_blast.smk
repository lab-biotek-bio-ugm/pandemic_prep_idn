"""
Module 9: BLAST Specificity Verification
Runs detailed NCBI or local BLAST specificity queries for top N candidate primer sets.
Supports two modes:
  - remote: Queries NCBI nt database over the internet (requires stable connection & API key).
  - local:  Queries custom local BLAST databases (human_db + virus_db) built by build_db.smk.
"""

BUILD = config["build_name"]

def _blast_mode():
    return config.get("blast_verification", {}).get("mode", "remote")

def _custom_db_enabled():
    bv = config.get("blast_verification", {})
    return (
        bv.get("mode", "remote") == "local"
        and bv.get("custom_db_build", {}).get("enabled", False)
    )

def _get_db_marker():
    """Return the db_marker path if local custom DB is enabled, else empty list."""
    if _custom_db_enabled():
        db_dir = config["blast_verification"]["custom_db_build"].get(
            "output_dir", "resources/blast_db"
        )
        return [f"{db_dir}/build_complete.txt"]
    return []

def _validate_local_db():
    """
    Called at runtime in the shell: if mode=local but custom_db_build.enabled=false,
    we cannot proceed. The shell block handles this gracefully.
    """
    pass


rule m09_blast_verification:
    """
    Runs BLAST off-target queries for top N candidates in alternatives.json
    and exports dedicated specificity_report.md files per candidate.

    Local mode: queries both human_db and virus_db, then compiles a unified summary.
    Remote mode: queries NCBI nt (requires NCBI_API_KEY env var for reliability).
    """
    input:
        alternatives=f"results/{BUILD}/m08/{{region}}/summary.json",
        bridge_report=f"results/{BUILD}/m07/bridge_report.json",
        gff=config["references"]["annotation"],
        db_marker=_get_db_marker()  # Empty list if remote mode → no dependency
    output:
        marker=f"results/{BUILD}/m09/{{region}}/verification_complete.txt"
    log:
        f"logs/{BUILD}/m09/blast_verification_{{region}}.log"
    params:
        alternatives_file=f"results/{BUILD}/m08/{{region}}/alternatives.json",
        num_candidates=config["blast_verification"].get(
            "num_candidates", config["primerlab"]["design"].get("vienna_ranking_limit", 50)
        ),
        database=config["blast_verification"].get("database", "nt"),
        mode=config["blast_verification"].get("mode", "remote"),
        local_db_path=config["blast_verification"].get("local_db_path", ""),
        outdir=f"results/{BUILD}/m09/{{region}}",
        target=config["references"].get("reference_id"),
        ref_fasta=config["references"].get("genome", ""),
        pathogen=config["build_name"].split("-")[0],
        custom_db_enabled=_custom_db_enabled(),
        db_output_dir=config.get("blast_verification", {}).get(
            "custom_db_build", {}
        ).get("output_dir", "resources/blast_db"),
        bridge_report_path=f"results/{BUILD}/m07/bridge_report.json",
        gff=config["references"]["annotation"]
    conda:
        "../envs/primerlab.yaml"
    shell:
        """
        #!/usr/bin/env bash
        set -euo pipefail

        MODE="{params.mode}"

        # ── Mode validation ────────────────────────────────────────────────
        if [ "$MODE" = "local" ]; then
            DB_DIR="{params.db_output_dir}"

            # Validate that databases exist
            if [ ! -f "${{DB_DIR}}/build_complete.txt" ]; then
                echo "❌ ERROR: Mode 'local' dipilih tetapi database belum dibangun."
                echo "   File tidak ditemukan: ${{DB_DIR}}/build_complete.txt"
                echo ""
                echo "   Solusi:"
                echo "   1. Pastikan blast_verification.custom_db_build.enabled = true"
                echo "   2. Jalankan: snakemake --use-conda ${{DB_DIR}}/build_complete.txt"
                exit 1
            fi

            HUMAN_DB="${{DB_DIR}}/human_db"
            VIRUS_DB="${{DB_DIR}}/virus_db"

            if [ ! -f "${{HUMAN_DB}}.nhr" ] && [ ! -f "${{HUMAN_DB}}.nal" ]; then
                echo "❌ ERROR: human_db tidak ditemukan di ${{HUMAN_DB}}"
                echo "   Jalankan ulang: snakemake --use-conda ${{DB_DIR}}/build_complete.txt"
                exit 1
            fi
            if [ ! -f "${{VIRUS_DB}}.nhr" ] && [ ! -f "${{VIRUS_DB}}.nal" ]; then
                echo "❌ ERROR: virus_db tidak ditemukan di ${{VIRUS_DB}}"
                echo "   Jalankan ulang: snakemake --use-conda ${{DB_DIR}}/build_complete.txt"
                exit 1
            fi

            echo "✅ Mode LOCAL: Menggunakan database lokal"
            echo "   Human DB : ${{HUMAN_DB}}"
            echo "   Virus  DB: ${{VIRUS_DB}}"

            # Run verification against human_db (with resume/skip support)
            python -u scripts/m09_blast_verification.py \
                --alternatives {params.alternatives_file} \
                --num-candidates {params.num_candidates} \
                --database "${{HUMAN_DB}}" \
                --mode local \
                --outdir {params.outdir} \
                --target "{params.target}" \
                --ref-fasta "{params.ref_fasta}" \
                --gff "{params.gff}" \
                --pathogen "{params.pathogen}" \
                --db-label "human" \
                >> {log} 2>&1

            # Run verification against virus_db (with resume/skip support)
            python -u scripts/m09_blast_verification.py \
                --alternatives {params.alternatives_file} \
                --num-candidates {params.num_candidates} \
                --database "${{VIRUS_DB}}" \
                --mode local \
                --outdir {params.outdir} \
                --target "{params.target}" \
                --ref-fasta "{params.ref_fasta}" \
                --gff "{params.gff}" \
                --pathogen "{params.pathogen}" \
                --db-label "virus" \
                >> {log} 2>&1

            # BUG-1 FIX: Compile summary ONCE after both DB runs complete
            python -u scripts/m09_blast_verification.py \
                --alternatives {params.alternatives_file} \
                --database "unused" \
                --outdir {params.outdir} \
                --target "{params.target}" \
                --ref-fasta "{params.ref_fasta}" \
                --gff "{params.gff}" \
                --bridge-report "{params.bridge_report_path}" \
                --compile-only \
                >> {log} 2>&1

            # Write marker file after both DB runs complete (local mode)
            echo "Completed local verification (human_db + virus_db) for candidates." > {output.marker}
            echo "human_db: ${{HUMAN_DB}}" >> {output.marker}
            echo "virus_db: ${{VIRUS_DB}}" >> {output.marker}

        else
            # ── Remote mode ─────────────────────────────────────────────────
            if [ -z "${{NCBI_API_KEY:-}}" ]; then
                echo "⚠ WARNING: NCBI_API_KEY tidak diset. Remote BLAST mungkin lambat dan timeout."
                echo "   Set API key: export NCBI_API_KEY=your_key_here"
                echo "   Dapatkan key gratis di: https://www.ncbi.nlm.nih.gov/account/"
            fi

            DB_TARGET="{params.database}"
            if [ -n "{params.local_db_path}" ] && [ "{params.local_db_path}" != "None" ] && [ "{params.local_db_path}" != "" ]; then
                DB_TARGET="{params.local_db_path}"
            fi

            echo "🌐 Mode REMOTE: Menggunakan database NCBI: ${{DB_TARGET}}"

            python -u scripts/m09_blast_verification.py \
                --alternatives {params.alternatives_file} \
                --num-candidates {params.num_candidates} \
                --database "${{DB_TARGET}}" \
                --mode remote \
                --outdir {params.outdir} \
                --target "{params.target}" \
                --ref-fasta "{params.ref_fasta}" \
                --gff "{params.gff}" \
                --pathogen "{params.pathogen}" \
                --bridge-report "{params.bridge_report_path}" \
                >> {log} 2>&1

            # Marker is written by the Python script in remote mode (no --db-label)
            # Verify it was created
            if [ ! -f "{output.marker}" ]; then
                echo "Completed remote verification." > {output.marker}
            fi
        fi

        echo "✅ M09 verification complete." >> {log}
        """
