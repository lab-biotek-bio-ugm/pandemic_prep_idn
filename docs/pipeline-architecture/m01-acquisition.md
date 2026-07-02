---
title: "M01: Data Acquisition & Standardisation"
description: "Fetches raw viral sequences from NCBI or local repositories, runs Nextclade for clade/lineage classification, and normalises metadata into the enriched TSV structure required by downstream Augur modules."
---

# Module M01: Data Acquisition & Standardisation

M01 is the pipeline entry point. It ingests raw viral genomic data from two selectable pathways (NCBI API or local files), executes Nextclade for quality-controlled clade classification, and outputs a standardised metadata TSV and cleaned FASTA that all downstream modules depend upon.

---

## Input Requirements

| Parameter | Config Key | Description |
|---|---|---|
| Data source selector | `database` | `"ncbi"` or `"local"` |
| NCBI organism query | `ncbi_query.organism` | Taxon name (e.g., `"Influenza A virus"`) |
| NCBI location filter | `ncbi_query.location` | Geographic filter (e.g., `"Indonesia"`) |
| Local FASTA path | `local_input.sequences` | Path to pre-downloaded FASTA file |
| Local metadata path | `local_input.metadata` | Path to metadata file |
| Local metadata format | `local_input.format` | See format table below |

### Supported Local Metadata Formats

| Format Value | Source | Parser Behaviour |
|---|---|---|
| `ncbi_virus_csv` | NCBI Virus web portal | Column remapping + country extraction from `Geo_Location` |
| `gisaid_tsv` | GISAID EpiFlu TSV export | Copied as-is; normalisation deferred to `sanitize_metadata.py` |
| `gisaid_epiflu` | GISAID EpiFlu XLS export | Copied as-is; Excel parsing handled by `sanitize_metadata.py` |
| `gisaid_raw` | GISAID raw download | Treated identically to `gisaid_tsv` |
| `tsv` | Any tab-separated file | Copied as-is |
| `auto` | Detection by extension | Falls back to TSV copy behaviour |

<Info>
The NCBI Datasets API pathway (`database: ncbi`) is marked **experimental** in the rule code. For reproducible production runs, use `database: local` with files manually downloaded from NCBI Virus or GISAID.
</Info>

---

## Processing Stages

### Stage 1A — NCBI Pathway (when `database: ncbi`)

**Rule: `download_ncbi_data`**

Calls the `ncbi-datasets-cli` to download a ZIP archive of genomic sequences matching the taxon and location query.

```bash
datasets download virus genome taxon "Influenza A virus" \
    --geo-location "Indonesia" \
    --released-after "2020-01-01" \
    --include genome \
    --filename results/{BUILD}/m01/ncbi_dataset.zip
```

**Rule: `extract_ncbi_data`**

Unpacks the ZIP, extracts `genomic.fna`, and converts the JSONL data report to a raw TSV metadata table using the `dataformat` utility:

```bash
dataformat tsv virus-genome \
    --inputfile ncbi_dataset/data/data_report.jsonl \
    --fields accession,virus-name,geo-location,geo-region,collection-date,release-date,length,host-name \
    > raw_metadata.jsonl
```

**Rule: `convert_metadata_tsv`** *(Python in-process)*

Remaps raw NCBI column names to the Nextstrain standard schema using a hardcoded `FIELD_MAP` dict. Parses `country` by splitting `Geographic Location` on `":"`.

### Stage 1B — Local Pathway (when `database: local`)

**Rule: `copy_local_sequences`** — Simple file copy from `local_input.sequences` path.

**Rule: `convert_local_metadata`** — Dispatches to the appropriate parser based on `local_input.format`. For `ncbi_virus_csv`, applies a full `COLUMN_MAP` conversion including `country`/`region` extraction from `Geo_Location`. For all TSV/GISAID formats, copies as-is.

### Stage 2 — Metadata Sanitisation (both pathways)

**Rule: `sanitize_metadata`** → Script: `scripts/sanitize_metadata.py`

Normalises column names across heterogeneous sources. Handles:
- GISAID EpiFlu XLS (direct Excel reading via `openpyxl`)
- Divergent date formats (ISO 8601, GISAID partial dates `YYYY-XX-XX` → `YYYY`)
- Strain name deduplication and whitespace stripping
- Country/region extraction from compound location strings

### Stage 3 — Nextclade Classification (optional, `nextclade.enabled: true`)

**Rule: `nextclade_dataset`** — Downloads the Nextclade reference dataset:

```bash
nextclade dataset get \
    --name 'flu_h3n2_ha' \
    --output-zip results/{BUILD}/m01/nextclade_dataset.zip
```

**Rule: `run_nextclade`** — Assigns clade and lineage labels to all sequences in the input FASTA:

```bash
nextclade run \
    --input-dataset results/{BUILD}/m01/nextclade_dataset.zip \
    --input-fasta results/{BUILD}/m01/raw_sequences.fasta \
    --output-tsv results/{BUILD}/m01/nextclade_output.tsv \
    --output-all results/{BUILD}/m01/nextclade/
```

**Rule: `enrich_metadata`** → Script: `scripts/enrich_metadata.py`

Left-joins the Nextclade output TSV onto the sanitised metadata on the `strain` column. Appends `clade`, `lineage`, `aaSubstitutions`, `qcOverallStatus` fields to produce the final `metadata_enriched.tsv`.

---

## Outputs

| File | Path | Description |
|---|---|---|
| Raw sequences | `results/{BUILD}/m01/raw_sequences.fasta` | Unfiltered input sequences |
| Sanitised metadata | `results/{BUILD}/m01/metadata.tsv` | Column-normalised TSV |
| **Enriched metadata** | `results/{BUILD}/m01/metadata_enriched.tsv` | **Primary M01 output** — with Nextclade clade/lineage/QC columns |
| Nextclade TSV | `results/{BUILD}/m01/nextclade_output.tsv` | Raw per-sequence Nextclade output |
| Download log | `logs/{BUILD}/m01_download.log` | NCBI or copy operation log |

### Enriched Metadata Schema (Key Columns)

| Column | Source | Notes |
|---|---|---|
| `strain` | All sources | Primary join key; must be unique |
| `virus` | All sources | Organism/taxon string |
| `date` | All sources | ISO 8601 (`YYYY-MM-DD`); partial dates coerced |
| `country` | All sources | Extracted from location string |
| `region` | All sources | Sub-national location if available |
| `clade` | Nextclade | e.g., `3C.2a1b.2a.2`, `XBB.1.5` |
| `lineage` | Nextclade | Unified nomenclature lineage |
| `qcOverallStatus` | Nextclade | `good` / `mediocre` / `bad` |

---

## Conda Environment

```
envs/nextstrain.yaml
```

Key packages: `augur`, `nextclade`, `ncbi-datasets-cli`, `dataformat`, `pandas`, `openpyxl`

---

## Debugging & Manual Invocation

Run Nextclade independently on a FASTA subset to verify clade assignment before launching the full pipeline:

```bash
nextclade run \
    --input-dataset results/{BUILD}/m01/nextclade_dataset.zip \
    --input-fasta data/test_subset.fasta \
    --output-tsv results/debug_m01/nextclade_output.tsv \
    --output-all results/debug_m01/
```

Force-re-run only M01 rules without rebuilding the rest of the DAG:

```bash
snakemake --use-conda \
    results/{BUILD}/m01/metadata_enriched.tsv \
    --forcerun sanitize_metadata enrich_metadata \
    --configfile config/config_influenza_a_mp.yaml
```

---

## Known Limitations & Edge Cases

- **Partial dates**: GISAID sequences with collection dates like `2024-XX-XX` are coerced to `2024`. This intentional loss of precision may introduce small biases in `augur refine` clock fitting. Sequences with `date: unknown` are excluded.
- **NCBI API throttling**: The NCBI API enforces rate limits. Set the `NCBI_API_KEY` environment variable to raise the limit from 3 to 10 requests/second.
- **Nextclade dataset versioning**: Nextclade datasets evolve. Pin `nextclade.dataset_name` to a versioned identifier in long-term surveillance studies to ensure reproducible clade calls across runs.
- **Subtype underrepresentation**: If a subtype has fewer than `global_params.min_subtype_count` sequences in the input, all sequences of that subtype are excluded by `qc_prefilter` in M02 to prevent tree artifacts from low-count groups.
