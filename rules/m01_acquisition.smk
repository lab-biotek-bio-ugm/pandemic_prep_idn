"""
Module 1: Data Acquisition & QC
Handles fetching raw viral sequences from NCBI or GISAID, running Nextclade
for clade/lineage classification + quality control, and formatting
the metadata into the enriched TSV structure required by Nextstrain Augur.

Data source is selected via config["database"]: "ncbi" or "gisaid"
"""

BUILD = config["build_name"]
DATABASE = config["database"]  # "ncbi" or "gisaid"

# ==========================================================================
# Source-specific download rules (only one will be active based on config)
# ==========================================================================

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

elif DATABASE == "gisaid":
    # ------------------------------------------------------------------
    # GISAID EpiCoV pathway (skeleton — requires credentials)
    # ------------------------------------------------------------------
    rule download_gisaid_data:
        output:
            fasta=f"results/{BUILD}/m01/raw_sequences.fasta",
            tsv=f"results/{BUILD}/m01/metadata.tsv"
        params:
            virus_name=config["gisaid_query"]["virus_name"],
            location=config["gisaid_query"]["location"],
            date_from=config["gisaid_query"].get("collection_date_from", "")
        log:
            f"logs/{BUILD}/m01_download.log"
        shell:
            """
            echo "================================================================" | tee {log}
            echo "  GISAID Data Source Selected" | tee -a {log}
            echo "================================================================" | tee -a {log}
            echo "" | tee -a {log}
            echo "  GISAID requires manual data placement or API credentials." | tee -a {log}
            echo "  Please place your files at:" | tee -a {log}
            echo "    - {output.fasta}" | tee -a {log}
            echo "    - {output.tsv}" | tee -a {log}
            echo "" | tee -a {log}
            echo "  Download from: https://epicov.org/epi3" | tee -a {log}
            echo "================================================================" | tee -a {log}

            if [ ! -f "{output.fasta}" ] || [ ! -f "{output.tsv}" ]; then
                echo "ERROR: GISAID files not found." | tee -a {log}
                exit 1
            fi
            """

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
        """Convert NCBI Virus CSV metadata to pipeline TSV format."""
        input:
            csv=config["local_input"]["metadata"]
        output:
            tsv=f"results/{BUILD}/m01/metadata.tsv"
        params:
            fmt=config["local_input"].get("format", "ncbi_virus_csv")
        log:
            f"logs/{BUILD}/m01_convert_metadata.log"
        run:
            import csv as csvmod
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

                elif fmt == "tsv":
                    # --- Already in TSV format, just copy ---
                    import shutil
                    shutil.copy2(str(input.csv), str(output.tsv))
                    logf.write(f"Copied pre-formatted TSV → {output.tsv}\n")

                else:
                    logf.write(f"ERROR: Unknown metadata format: {fmt}\n")
                    raise ValueError(f"Unknown local_input format: '{fmt}'. Use 'ncbi_virus_csv' or 'tsv'.")

else:
    raise ValueError(f"Unknown database provider: '{DATABASE}'. Must be 'ncbi', 'gisaid', or 'local'.")


# ==========================================================================
# Shared rules (run regardless of data source)
# ==========================================================================

# ---------------------------------------------------------------------------
# Nextclade dataset download (isolated to specific build directory)
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Run Nextclade for QC + Clade/Lineage classification
# ---------------------------------------------------------------------------
rule run_nextclade:
    input:
        sequences=f"results/{BUILD}/m01/raw_sequences.fasta",
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
        """

# ---------------------------------------------------------------------------
# Enrich metadata with Nextclade clade/lineage/QC columns
# ---------------------------------------------------------------------------
rule enrich_metadata:
    input:
        base=f"results/{BUILD}/m01/metadata.tsv",
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
