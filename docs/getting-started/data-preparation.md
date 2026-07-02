---
title: "Data Preparation"
description: "Comprehensive guide to preparing, formatting, and sourcing raw genomic data and metadata for the pipeline. Covers GISAID EpiFlu, NCBI Virus, and custom laboratory datasets."
---

# Data Preparation

The pipeline accepts genomic data from three primary sources: GISAID EpiFlu, NCBI Virus, and custom laboratory collections. This page documents the required file formats, metadata schemas, and acquisition workflows for each source.

---

## Overview: Supported Data Sources

| Source | Recommended For | Format | Config `database` |
|---|---|---|---|
| GISAID EpiFlu | Influenza A/B surveillance | XLS (EpiFlu export) or TSV | `local` |
| NCBI Virus (web portal) | Any pathogen with NCBI records | CSV | `local` |
| NCBI Datasets API | Automated large-scale retrieval | Auto (via `ncbi-datasets-cli`) | `ncbi` |
| Custom laboratory | In-house sequencing | FASTA + any TSV | `local` |

---

## Source 1: GISAID EpiFlu

GISAID is the primary recommended source for Influenza A/B sequences due to its comprehensive coverage and standardised metadata.

### Downloading from GISAID

1. Log in at [https://gisaid.org](https://gisaid.org) → Select **EpiFlu™ Database**.
2. Apply filters:
   - **Type**: `A` (Influenza A) or `B` (Influenza B)
   - **Segment**: Select the target segment (e.g., `MP` for M segment)
   - **Host**: `Human`
   - **Location**: Your target region
3. Select all results → Click **Download** → Choose **Sequences (XLS/CSV) + Metadata**.
4. Download format: **EpiFlu™ Package** — this produces a `.xls` metadata file paired with a `.fasta` sequences file.

### GISAID Metadata Schema (EpiFlu XLS)

The pipeline's `sanitize_metadata.py` reads EpiFlu XLS directly via `openpyxl`. Key column mapping:

| GISAID Column | Internal Column | Notes |
|---|---|---|
| `Isolate_Id` | `strain` | Primary join key |
| `Type` | `virus` | e.g., `A / H3N2` |
| `Subtype` | `subtype` | e.g., `H3N2`, `H1N1` |
| `Lineage` | `lineage` | Clade/lineage if assigned |
| `Collection_Date` | `date` | ISO format coerced from GISAID partial dates |
| `Location` | `location` | Full location string (e.g., `Asia / Indonesia / Jakarta`) |
| `Host` | `host` | e.g., `Human` |
| `Passage` | `passage` | e.g., `Original`, `E1`, `MDCK` |

### Configuration for GISAID EpiFlu

```yaml
database: "local"
local_input:
  sequences: "dataset/EpiFlu/MP/EpiFlu_A_MP_Combined_Sequence.fasta"
  metadata:  "dataset/EpiFlu/MP/EpiFlu_A_MP_Combined_Metadata.csv"
  format:    "gisaid_epiflu"
```

<Info>
**GISAID Data Access Agreement**: All data downloaded from GISAID is subject to the GISAID Data Access Agreement. Do not redistribute raw GISAID sequences or metadata. The pipeline stores processed outputs only; raw GISAID files should be retained locally and not committed to public repositories.
</Info>

---

## Source 2: NCBI Virus (Web Portal)

NCBI Virus provides open-access sequence data for all pathogens. Download limits are higher than GISAID, but metadata is less standardised for Influenza.

### Downloading from NCBI Virus

1. Navigate to [https://www.ncbi.nlm.nih.gov/labs/virus/vssi/](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/).
2. Search by taxon (e.g., **Influenza A virus**) and apply filters:
   - Nucleotide completeness: `Complete`
   - Host: `Homo sapiens`
   - Geographic Location: Target country
   - Collection Date: Start and end year range
3. Click **Download** → Select both:
   - **Nucleotide** (FASTA format)
   - **Current table view result** (CSV metadata)
4. Save both files and note their paths.

### NCBI Virus Metadata Schema (CSV)

| NCBI Column | Internal Column | Notes |
|---|---|---|
| `Accession` | `strain` | GenBank accession (e.g., `KF509985`) |
| `Organism_Name` | `virus` | e.g., `Influenza A virus (A/Indonesia/...)` |
| `Geo_Location` | `location` / `country` | Country extracted by splitting on `:` |
| `Host` | `host` | e.g., `Homo sapiens` |
| `Collection_Date` | `date` | ISO 8601 or partial date |
| `Release_Date` | `date_submitted` | NCBI public release date |
| `Length` | `length` | Sequence length in bp |

### Configuration for NCBI Virus CSV

```yaml
database: "local"
local_input:
  sequences: "data/ncbi/ncbi_sequences.fasta"
  metadata:  "data/ncbi/ncbi_metadata.csv"
  format:    "ncbi_virus_csv"
```

---

## Source 3: NCBI Datasets API (Automated)

For large-scale or automated retrieval, use the `ncbi-datasets-cli` pathway.

<Warning>
The NCBI API pathway is currently marked **experimental** in the pipeline. The API is evolving and may return inconsistent results for some organism queries. For production surveillance runs, use `database: local` with manually downloaded files.
</Warning>

```yaml
database: "ncbi"
ncbi_query:
  organism:      "Influenza A virus"
  location:      "Indonesia"
  released_after: "2020-01-01"   # Optional
```

The pipeline will call `ncbi-datasets-cli` to download sequences and auto-convert metadata.

**Prerequisite:** Set your NCBI API key for higher rate limits:

```bash
export NCBI_API_KEY="your_api_key_here"
# Get a key at: https://www.ncbi.nlm.nih.gov/account/settings/
```

---

## Sequence File Requirements (FASTA)

Regardless of source, the input FASTA must meet these requirements:

| Requirement | Value | Notes |
|---|---|---|
| Format | Standard FASTA | Headers starting with `>` |
| Header format | `>AccessionID` or `>StrainName` | Must match `strain` column in metadata |
| Character set | A, T, G, C, N, `-` | IUPAC degenerate bases accepted |
| Minimum length | Configured via `subsampling.min_length` | Default 800 bp for Influenza segments |
| Encoding | UTF-8 | No special characters in headers |

### FASTA Header Alignment with Metadata

The `strain` field in metadata **must exactly match** the FASTA sequence header (after the `>` character). Mismatches cause sequences to be silently dropped during metadata enrichment.

```bash
# Verify header-to-metadata alignment
grep "^>" data/sequences.fasta | sed 's/>//' | sort > /tmp/fasta_ids.txt
awk -F'\t' 'NR>1 {print $1}' data/metadata.tsv | sort > /tmp/meta_ids.txt
comm -23 /tmp/fasta_ids.txt /tmp/meta_ids.txt  # Sequences without metadata
comm -13 /tmp/fasta_ids.txt /tmp/meta_ids.txt  # Metadata without sequences
```

---

## Combining Multiple Downloads

When combining multiple GISAID or NCBI downloads:

```bash
# Combine FASTA files (ensure no duplicate headers)
cat download1.fasta download2.fasta > combined.fasta

# Check for duplicate strain IDs
grep "^>" combined.fasta | sort | uniq -d

# Combine metadata CSVs (ensure identical column structure)
# For CSV: use header from first file only
head -n 1 metadata1.csv > combined_metadata.csv
tail -n +2 metadata1.csv >> combined_metadata.csv
tail -n +2 metadata2.csv >> combined_metadata.csv
```

The pipeline's MD5 deduplication in M02 will handle any exact-duplicate sequences that slip through the combination step.

---

## Reference Files

The pipeline requires two reference files per target segment:

| File | Format | Description | Example |
|---|---|---|---|
| `references.genome` | FASTA or GenBank `.gb` | Reference genome for alignment | `NC_002016.1` (Flu A MP) |
| `references.annotation` | GFF3 | CDS annotation for gene extraction | `flu_a_mp_annotation.gff` |

Reference files are stored in `config/references/`. To add a new pathogen:

1. Download the reference GenBank record: `efetch -db nuccore -id NC_002016.1 -format gb > config/references/flu_a_mp_reference.gb`
2. Convert to GFF3 using `bp_seqfeature_load.pl` or manually create from the GenBank CDS annotations.
3. Update `config.references.genome` and `config.references.annotation` to point to the new files.
4. Update `config.references.reference_id` to the accession number (e.g., `NC_002016.1`).

---

## Dataset Sizing Recommendations

| Dataset Size | Recommended `max_sequences` | Notes |
|---|---|---|
| < 500 sequences | No limit needed | Run full dataset |
| 500 – 3,000 | `3000` | Default, balances accuracy vs. speed |
| 3,000 – 10,000 | `3000` with `group_by: "country year month"` | Stratified subsampling essential |
| > 10,000 | `2000 – 3000` + focal/contextual mode | Required for tractable tree inference |

For HyPhy selection analysis (M05), independently cap the input via `hyphy.subsample.max_sequences: 300` regardless of the total dataset size.
