"""
Module 1: Data Acquisition & QC
Handles fetching raw viral sequences from NCBI GenBank, running Nextclade
for clade/lineage classification + quality control, and formatting
the metadata into the enriched TSV structure required by Nextstrain Augur.
"""

BUILD = config["build_name"]

# ---------------------------------------------------------------------------
# M1a: Download raw sequences and metadata from NCBI Datasets
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# M1b: Extract FASTA and JSONL from the downloaded zip
# ---------------------------------------------------------------------------
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

        # Generate JSONL metadata from the data report
        dataformat tsv virus-genome \
            --inputfile $TMPDIR/ncbi_dataset/data/data_report.jsonl \
            --fields accession,virus-name,geo-location,geo-region,collection-date,release-date,length,host-name \
            > {output.jsonl}

        rm -rf $TMPDIR
        """

# ---------------------------------------------------------------------------
# M1c: Convert NCBI metadata JSONL/TSV into Augur-compatible metadata TSV
# ---------------------------------------------------------------------------
rule convert_metadata_tsv:
    input:
        raw_meta=f"results/{BUILD}/m01/raw_metadata.jsonl"
    output:
        tsv=f"results/{BUILD}/m01/metadata.tsv"
    run:
        import csv

        # Map NCBI field names -> Augur/Nextstrain standard column names
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
                # Set defaults for required Augur columns that NCBI may not provide
                mapped.setdefault("country", mapped.get("location", "").split(":")[0].strip() if mapped.get("location") else "")

                if writer is None:
                    writer = csv.DictWriter(fout, fieldnames=list(mapped.keys()), delimiter="\t")
                    writer.writeheader()
                writer.writerow(mapped)
        # Handle empty metadata edge case
        if writer is None:
            with open(output.tsv, "w") as fout:
                fout.write("strain\tvirus\tlocation\tregion\tdate\tdate_submitted\tlength\thost\tcountry\n")

# ---------------------------------------------------------------------------
# M1d: Download Nextclade reference dataset (cached locally)
# ---------------------------------------------------------------------------
rule nextclade_dataset:
    output:
        dir=directory("data/nextclade_dataset")
    conda:
        "../envs/nextstrain.yaml"
    log:
        f"logs/{BUILD}/m01_nextclade_dataset.log"
    shell:
        """
        nextclade dataset get \
            --name 'nextstrain/sars-cov-2/wuhan-hu-1/orfs' \
            --output-dir {output.dir} \
            2>&1 | tee {log}
        """

# ---------------------------------------------------------------------------
# M1e: Run Nextclade for QC + Clade/Lineage classification
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
# M1f: Enrich metadata with Nextclade clade/lineage/QC columns
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
