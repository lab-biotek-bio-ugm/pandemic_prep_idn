"""
Module 1: Data Acquisition & QC
Handles fetching raw viral sequences from NCBI or local files, running Nextclade
for clade/lineage classification + quality control, and formatting
the metadata into the enriched TSV structure required by Nextstrain Augur.

Data source is selected via config["database"]: "ncbi" or "local"
"""

BUILD = config["build_name"]
DATABASE = config["database"]  # "ncbi" or "local"

# Source-specific download rules (only one will be active based on config)

if DATABASE == "ncbi":
    # ------------------------------------------------------------------
    # NCBI Datasets API pathway
    # ------------------------------------------------------------------
    rule download_ncbi_data:
        output:
            zip_file=temp(f"results/{BUILD}/m01/ncbi_dataset.zip")
        params:
            organism=config["ncbi_query"]["organism"],
            location=config["ncbi_query"]["location"],
            released_after=config["ncbi_query"].get("released_after", "")
        conda:
            "../envs/nextstrain.yaml"
        log:
            f"logs/{BUILD}/m01_download.log"
        shell:
            """
            mkdir -p results/{BUILD}/m01

            echo "================================================================" | tee {log}
            echo "  WARNING: The NCBI Datasets API pathway is currently"            | tee -a {log}
            echo "           experimental and might not work as expected."           | tee -a {log}
            echo "           Consider using 'database: local' instead."              | tee -a {log}
            echo "================================================================" | tee -a {log}

            RELEASED_FLAG=""
            if [ -n "{params.released_after}" ]; then
                RELEASED_FLAG="--released-after {params.released_after}"
            fi

            datasets download virus genome taxon "{params.organism}" \
                --geo-location "{params.location}" \
                $RELEASED_FLAG \
                --include genome \
                --filename {output.zip_file} \
                2>&1 | tee {log}
            """

    rule extract_ncbi_data:
        input:
            zip_file=f"results/{BUILD}/m01/ncbi_dataset.zip"
        output:
            fasta=f"results/{BUILD}/m01/raw_sequences.fasta",
            jsonl=f"results/{BUILD}/m01/raw_metadata.jsonl"
        conda:
            "../envs/nextstrain.yaml"
        log:
            f"logs/{BUILD}/m01_extract.log"
        shell:
            """
            TMPDIR=$(mktemp -d)

            unzip -o {input.zip_file} -d $TMPDIR 2>&1 | tee {log}

            # Copy the genomic FASTA
            cp $TMPDIR/ncbi_dataset/data/genomic.fna {output.fasta}

            # Generate TSV metadata from the data report
            dataformat tsv virus-genome \
                --inputfile $TMPDIR/ncbi_dataset/data/data_report.jsonl \
                --fields accession,virus-name,geo-location,geo-region,collection-date,release-date,length,host-name \
                > {output.jsonl}

            rm -rf $TMPDIR
            """

    rule convert_metadata_tsv:
        input:
            raw_meta=f"results/{BUILD}/m01/raw_metadata.jsonl"
        output:
            tsv=f"results/{BUILD}/m01/metadata.tsv"
        run:
            import csv

            FIELD_MAP = {
                "Accession": "strain",
                "Virus Name": "virus",
                "Geographic Location": "location",
                "Geographic Region": "region",
                "Collection Date": "date",
                "Release Date": "date_submitted",
                "Length": "length",
                "Host Name": "host"
            }

            with open(input.raw_meta, "r") as fin, open(output.tsv, "w", newline="") as fout:
                reader = csv.DictReader(fin, delimiter="\t")
                writer = None
                for row in reader:
                    mapped = {}
                    for old_key, new_key in FIELD_MAP.items():
                        mapped[new_key] = row.get(old_key, "")
                    mapped.setdefault("country", mapped.get("location", "").split(":")[0].strip() if mapped.get("location") else "")

                    if writer is None:
                        writer = csv.DictWriter(fout, fieldnames=list(mapped.keys()), delimiter="\t")
                        writer.writeheader()
                    writer.writerow(mapped)
            if writer is None:
                with open(output.tsv, "w") as fout:
                    fout.write("strain\tvirus\tlocation\tregion\tdate\tdate_submitted\tlength\thost\tcountry\n")

elif DATABASE == "local":
    # ------------------------------------------------------------------
    # Local / Manual data pathway
    # Supports files downloaded from NCBI Virus web portal or any
    # manually prepared FASTA + metadata.
    # ------------------------------------------------------------------
    rule copy_local_sequences:
        """Copy manually-placed FASTA sequences to pipeline results directory."""
        input:
            fasta=config["local_input"]["sequences"]
        output:
            fasta=f"results/{BUILD}/m01/raw_sequences.fasta"
        log:
            f"logs/{BUILD}/m01_download.log"
        shell:
            """
            mkdir -p results/{BUILD}/m01
            cp {input.fasta} {output.fasta}
            echo "Copied local sequences: $(grep -c '^>' {output.fasta}) records" | tee {log}
            """

    rule convert_local_metadata:
        """Copy or convert local metadata to pipeline working directory.
        For CSV formats (ncbi_virus_csv): converts to TSV with column mapping.
        For TSV formats (gisaid_tsv, gisaid_raw, tsv, auto): copies as-is.
        Column normalization is handled downstream by sanitize_metadata rule.
        """
        input:
            csv=config["local_input"]["metadata"]
        output:
            tsv=f"results/{BUILD}/m01/metadata.tsv"
        params:
            fmt=config["local_input"].get("format", "auto")
        log:
            f"logs/{BUILD}/m01_convert_metadata.log"
        run:
            import csv as csvmod
            import shutil
            import sys

            fmt = params.fmt
            log_path = str(log[0])

            with open(log_path, "w") as logf:
                if fmt == "ncbi_virus_csv":
                    # --- NCBI Virus web portal CSV conversion ---
                    COLUMN_MAP = {
                        "Accession": "strain",
                        "Organism_Name": "virus",
                        "Pangolin": "pangolin_lineage",
                        "Isolate": "isolate",
                        "Length": "length",
                        "Nuc_Completeness": "completeness",
                        "Geo_Location": "location",
                        "Host": "host",
                        "Collection_Date": "date",
                        "Release_Date": "date_submitted",
                    }

                    with open(str(input.csv), "r", encoding="utf-8") as fin, \
                         open(str(output.tsv), "w", newline="", encoding="utf-8") as fout:
                        reader = csvmod.DictReader(fin)
                        out_fields = list(COLUMN_MAP.values()) + ["country", "region"]
                        writer = csvmod.DictWriter(fout, fieldnames=out_fields, delimiter="\t")
                        writer.writeheader()

                        count = 0
                        for row in reader:
                            out_row = {}
                            for csv_col, tsv_col in COLUMN_MAP.items():
                                out_row[tsv_col] = row.get(csv_col, "")
                            geo = row.get("Geo_Location", "")
                            parts = geo.split(":", 1)
                            out_row["country"] = parts[0].strip() if parts else ""
                            out_row["region"] = parts[1].strip() if len(parts) > 1 else ""
                            writer.writerow(out_row)
                            count += 1

                    logf.write(f"Converted {count} NCBI Virus CSV records → {output.tsv}\n")
                    print(f"Converted {count} records from NCBI Virus CSV", file=sys.stderr)

                elif fmt in ("tsv", "gisaid_tsv", "gisaid_raw", "auto"):
                    # --- TSV-based formats: copy as-is ---
                    # Column normalization is handled by sanitize_metadata.py
                    shutil.copy2(str(input.csv), str(output.tsv))
                    logf.write(f"Copied {fmt} metadata → {output.tsv}\n")
                    print(f"Copied {fmt} metadata (normalization deferred to sanitizer)",
                          file=sys.stderr)

                elif fmt == "gisaid_epiflu":
                    # --- GISAID EpiFlu XLS format: copy as-is ---
                    # sanitize_metadata.py reads Excel directly via --orig-file
                    shutil.copy2(str(input.csv), str(output.tsv))
                    logf.write(f"Copied EpiFlu XLS metadata → {output.tsv}\n")
                    print(f"Copied EpiFlu XLS metadata (Excel parsing deferred to sanitizer)",
                          file=sys.stderr)

                else:
                    logf.write(f"ERROR: Unknown metadata format: {fmt}\n")
                    raise ValueError(
                        f"Unknown local_input format: '{fmt}'. "
                        f"Use 'auto', 'ncbi_virus_csv', 'gisaid_tsv', 'gisaid_raw', "
                        f"'gisaid_epiflu', or 'tsv'."
                    )

else:
    raise ValueError(f"Unknown database provider: '{DATABASE}'. Must be 'ncbi' or 'local'.")


# Shared rules (run regardless of data source)

NEXTCLADE_ENABLED = config.get("nextclade", {}).get("enabled", True)

# Nextclade dataset download + run (only when enabled)
if NEXTCLADE_ENABLED:
    rule nextclade_dataset:
        output:
            zip_file=f"results/{BUILD}/m01/nextclade_dataset.zip"
        params:
            dataset_name=config["nextclade"]["dataset_name"]
        conda:
            "../envs/nextstrain.yaml"
        log:
            f"logs/{BUILD}/m01_nextclade_dataset.log"
        shell:
            """
            nextclade dataset get \
                --name '{params.dataset_name}' \
                --output-zip {output.zip_file} \
                2>&1 | tee {log}
            """

    rule run_nextclade:
        input:
            sequences=f"results/{BUILD}/m01/early_filtered_sequences.fasta",
            dataset=f"results/{BUILD}/m01/nextclade_dataset.zip"
        output:
            tsv=f"results/{BUILD}/m01/nextclade.tsv"
        conda:
            "../envs/nextstrain.yaml"
        log:
            f"logs/{BUILD}/m01_nextclade_run.log"
        shell:
            """
            nextclade run \
                --input-dataset {input.dataset} \
                --output-tsv {output.tsv} \
                {input.sequences} \
                2>&1 | tee {log}
                2>&1 | tee {log}
            """

EARLY_FILTER_ENABLED = config.get("early_filter", {}).get("enabled", False)

rule early_filter_data:
    """Pre-filters raw sequences and sanitized metadata based on user queries
    before passing them to downstream systems like Nextclade.
    """
    input:
        sequences=f"results/{BUILD}/m01/raw_sequences.fasta",
        metadata=f"results/{BUILD}/m01/metadata_sanitized.tsv"
    output:
        sequences=f"results/{BUILD}/m01/early_filtered_sequences.fasta",
        metadata=f"results/{BUILD}/m01/early_filtered_metadata.tsv"
    params:
        query=config.get("early_filter", {}).get("query", ""),
        min_date=config.get("early_filter", {}).get("min_date", ""),
        max_date=config.get("early_filter", {}).get("max_date", "")
    conda:
        "../envs/nextstrain.yaml"
    log:
        f"logs/{BUILD}/m01_early_filter.log"
    shell:
        """
        if [ "{EARLY_FILTER_ENABLED}" = "True" ]; then
            python scripts/early_filter.py \
                --sequences {input.sequences} \
                --metadata {input.metadata} \
                --output-sequences {output.sequences} \
                --output-metadata {output.metadata} \
                --query "{params.query}" \
                --min-date "{params.min_date}" \
                --max-date "{params.max_date}" 2>&1 | tee {log}
        else
            echo "[EARLY FILTER] Disabled in config. Copying raw data as-is." | tee {log}
            cp {input.sequences} {output.sequences}
            cp {input.metadata} {output.metadata}
        fi
        """

# Sanitize metadata: normalize columns, fill defaults, drop empty IDs
rule sanitize_metadata:
    """Normalizes metadata from any source (GISAID, NCBI, generic TSV) into
    a consistent internal format. Applies 3-tier column validation:
    Tier 1 (strain/date/country) → drop or fill with warning.
    Tier 2 (region/host/division) → fill silently.
    Tier 3 (everything else) → pass-through.
    """
    input:
        tsv=f"results/{BUILD}/m01/metadata.tsv"
    output:
        tsv=f"results/{BUILD}/m01/metadata_sanitized.tsv",
        report=f"results/{BUILD}/m01/sanitize_report.json"
    params:
        fmt=config["local_input"].get("format", "auto"),
        orig_file=lambda w: config["local_input"]["metadata"] if DATABASE == "local" else ""
    conda:
        "../envs/nextstrain.yaml"
    log:
        f"logs/{BUILD}/m01_sanitize_metadata.log"
    shell:
        """
        python scripts/sanitize_metadata.py \
            --input {input.tsv} \
            --output {output.tsv} \
            --report {output.report} \
            --format {params.fmt} \
            --orig-file "{params.orig_file}" 2>&1 | tee {log}
        """

# Enrich metadata: with Nextclade (if enabled) or passthrough (if disabled)
if NEXTCLADE_ENABLED:
    rule enrich_metadata:
        input:
            base=f"results/{BUILD}/m01/early_filtered_metadata.tsv",
            nc=f"results/{BUILD}/m01/nextclade.tsv"
        output:
            tsv=f"results/{BUILD}/m01/metadata_enriched.tsv"
        conda:
            "../envs/nextstrain.yaml"
        log:
            f"logs/{BUILD}/m01_enrich_metadata.log"
        shell:
            """
            python scripts/enrich_metadata.py \
                --base-metadata {input.base} \
                --nextclade-tsv {input.nc} \
                --output {output.tsv} \
                2>&1 | tee {log}
            """
else:
    rule enrich_metadata_passthrough:
        """When Nextclade is disabled, copy sanitized metadata directly as
        the enriched output. All downstream rules read metadata_enriched.tsv
        so this maintains full compatibility.
        """
        input:
            base=f"results/{BUILD}/m01/early_filtered_metadata.tsv"
        output:
            tsv=f"results/{BUILD}/m01/metadata_enriched.tsv"
        log:
            f"logs/{BUILD}/m01_enrich_metadata.log"
        run:
            import shutil
            shutil.copy2(str(input.base), str(output.tsv))
            with open(str(log[0]), "w") as f:
                f.write("[ENRICH] Nextclade DISABLED — copied sanitized metadata as enriched.\n")
            print("[ENRICH] Nextclade DISABLED — metadata passthrough (no clade/lineage columns).")

# DEFERRED ITEMS for M01
#
# TODO [M01-1] Virus Template System (config/templates/)
#   What: Create per-virus YAML presets (e.g., influenza_a_h1n1.yaml,
#         sars_cov_2.yaml, rsv_a.yaml) containing pre-filled values for:
#         - references.genome, references.annotation, references.reference_id
#         - nextclade.dataset_name
#         - timetree.clock_rate, timetree.clock_std_dev
#         - hyphy.target_genes
#         - subsampling.min_length
#   Why:  Currently the user must manually change 6+ config parameters when
#         switching target virus. Templates would make it plug-and-play.
#   Prereqs: Pipeline must be validated on at least 2 different viruses first.
#   Priority: Phase 2 — after Influenza A and SARS-CoV-2 runs succeed.
#
# DROPPED: BV-BRC / ViPR Database Support
#   Original plan included BV-BRC (Bacterial and Viral Bioinformatics Resource
#   Center, successor to ViPR) as a third data source.
#   Decision: Dropped. BV-BRC data for respiratory viruses is largely a mirror
#   of NCBI GenBank. Adding it would increase code complexity without providing
#   unique sequences. NCBI + local files are sufficient for all use cases.
#
