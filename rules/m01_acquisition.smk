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
    # GISAID does not have a public CLI like NCBI. Data is typically:
    #   1. Downloaded manually from https://epicov.org/epi3
    #   2. Or via the GISAID API (requires institutional access)
    #   3. Or via augur curate (Nextstrain's GISAID ingest pipeline)
    #
    # This skeleton expects the user to place raw files manually at:
    #   results/{BUILD}/m01/raw_sequences.fasta
    #   results/{BUILD}/m01/metadata.tsv
    #
    # The user can also set up automated download via the Nextstrain
    # GISAID ingest workflow: https://github.com/nextstrain/ncov-ingest

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
            echo "  Expected metadata TSV columns:" | tee -a {log}
            echo "    strain, virus, date, region, country, location, host," | tee -a {log}
            echo "    date_submitted, length" | tee -a {log}
            echo "" | tee -a {log}
            echo "  Download from: https://epicov.org/epi3" | tee -a {log}
            echo "  Or use: https://github.com/nextstrain/ncov-ingest" | tee -a {log}
            echo "================================================================" | tee -a {log}

            # Check if files already exist (manually placed by user)
            if [ ! -f "{output.fasta}" ] || [ ! -f "{output.tsv}" ]; then
                echo "ERROR: GISAID files not found. Please download and place them manually." | tee -a {log}
                exit 1
            fi
            """

else:
    raise ValueError(f"Unknown database provider: '{DATABASE}'. Must be 'ncbi' or 'gisaid'.")


# ==========================================================================
# Shared rules (run regardless of data source)
# ==========================================================================

# ---------------------------------------------------------------------------
# Nextclade dataset download (cached locally)
# ---------------------------------------------------------------------------
rule nextclade_dataset:
    output:
        dir=directory("data/nextclade_dataset")
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
            --output-dir {output.dir} \
            2>&1 | tee {log}
        """

# ---------------------------------------------------------------------------
# Run Nextclade for QC + Clade/Lineage classification
# ---------------------------------------------------------------------------
rule run_nextclade:
    input:
        sequences=f"results/{BUILD}/m01/raw_sequences.fasta",
        dataset="data/nextclade_dataset"
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
