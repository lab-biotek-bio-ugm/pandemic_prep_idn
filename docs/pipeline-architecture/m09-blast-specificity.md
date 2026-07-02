---
title: "M09: BLAST Specificity Verification"
description: "Verifies primer and probe specificity by running BLAST searches of top candidate sequences against human and pathogen-specific databases. Supports both remote NCBI nt queries and local custom databases for offline or HPC workflows."
---

# Module M09: BLAST Specificity Verification

M09 is the off-target contamination screen for candidate primer sets output by M08. For each top-ranked candidate, it runs BLAST to confirm that the primer/probe sequences match the target viral genome and **do not cross-hybridise** to the human genome or closely related non-target pathogens.

This module operates in two modes controlled by `blast_verification.mode`:
- **`remote`**: Queries the NCBI `nt` database over the internet via BLAST+ API
- **`local`**: Queries custom local databases (`human_db` + `virus_db`) built by `build_db.smk`

---

## Input Requirements

| File | Source | Description |
|---|---|---|
| `results/{BUILD}/m08/{mode}/{region}/summary.json` | M08 | Ranked candidate primer list per region |
| `results/{BUILD}/m08/{mode}/{region}/alternatives.json` | M08 | Full candidate details including sequences |
| `results/{BUILD}/m07/bridge_report.json` | M07 | Bridge report for positional context |
| `config.references.annotation` | Config | GFF3 annotation for coding context |
| `resources/blast_db/build_complete.txt` | build_db.smk | *(Local mode only)* DB build sentinel |

### Configuration Parameters

| Config Key | Type | Default | Description |
|---|---|---|---|
| `blast_verification.mode` | String | `"remote"` | `"remote"` (NCBI) or `"local"` (custom DB) |
| `blast_verification.num_candidates` | Integer | From PrimerLab `vienna_ranking_limit` | Number of candidates to verify per region |
| `blast_verification.database` | String | `"nt"` | NCBI database (remote mode) |
| `blast_verification.local_db_path` | String | `""` | Path prefix to local BLAST DB files |
| `blast_verification.evalue_cutoff` | Float | `0.01` | E-value threshold for significant hits |
| `blast_verification.custom_db_build.enabled` | Boolean | `false` | Auto-build local DB via `build_db.smk` |
| `blast_verification.custom_db_build.output_dir` | String | `"resources/blast_db"` | Output directory for built local DBs |

---

## Processing

### Rule: `m09_blast_verification`

**Script:** `scripts/m09_blast_verification.py`  
**Environment:** `envs/primerlab.yaml`

The rule shell block performs upfront validation before calling the Python script:

**Local mode validation sequence:**

```bash
# 1. Check that local DB was built
if [ ! -f "${DB_DIR}/build_complete.txt" ]; then
    echo "ERROR: Mode 'local' selected but database not built."
    exit 1
fi

# 2. Validate human_db and virus_db exist
if [ ! -f "${HUMAN_DB}.nhr" ] && [ ! -f "${HUMAN_DB}.nal" ]; then
    echo "ERROR: human_db not found at ${HUMAN_DB}"
    exit 1
fi
```

**Local mode — dual database query:**

```bash
# Query 1: Check for human off-target cross-hybridisation
python -u scripts/m09_blast_verification.py \
    --alternatives results/{BUILD}/m08/{mode}/{region}/alternatives.json \
    --num-candidates 5 \
    --database resources/blast_db/human_db \
    --target {reference_id} \
    --ref-fasta config/reference.fasta \
    --output-dir results/{BUILD}/m09/{mode}/{region} \
    --bridge-report results/{BUILD}/m07/bridge_report.json \
    --gff config/reference.gff \
    --mode local

# Query 2: Check for specificity within pathogen-related sequences  
python -u scripts/m09_blast_verification.py \
    --alternatives results/{BUILD}/m08/{mode}/{region}/alternatives.json \
    --num-candidates 5 \
    --database resources/blast_db/virus_db \
    --target {reference_id} \
    --ref-fasta config/reference.fasta \
    --output-dir results/{BUILD}/m09/{mode}/{region} \
    --bridge-report results/{BUILD}/m07/bridge_report.json \
    --gff config/reference.gff \
    --mode local
```

**Remote mode query:**

```bash
python -u scripts/m09_blast_verification.py \
    --alternatives results/{BUILD}/m08/{mode}/{region}/alternatives.json \
    --num-candidates 5 \
    --database nt \
    --target {reference_id} \
    --ref-fasta config/reference.fasta \
    --output-dir results/{BUILD}/m09/{mode}/{region} \
    --bridge-report results/{BUILD}/m07/bridge_report.json \
    --gff config/reference.gff \
    --mode remote
```

<Note>
Remote mode requires a stable internet connection and preferably an NCBI API key for reliable throughput. Set `NCBI_API_KEY` as an environment variable to raise the query rate limit from 3 to 10 requests/second.
</Note>

---

## Building Local BLAST Databases (`build_db.smk`)

The `build_db.smk` rule file handles automated construction of the local databases. To trigger a database build:

```bash
# Build both human_db and virus_db
snakemake --use-conda \
    resources/blast_db/build_complete.txt \
    --configfile config/config_influenza_a_mp.yaml
```

**Databases created:**

| Database | Contents | Use |
|---|---|---|
| `human_db` | Human genome reference sequences | Off-target human cross-hybridisation check |
| `virus_db` | Target pathogen + related virus sequences | Within-genus specificity verification |

Each database is indexed with `makeblastdb -dbtype nucl`. For multi-file databases (e.g., human genome across multiple chromosomes), an alias file (`.nal`) is created to allow querying across all files.

**Manual database creation:**

```bash
# Build a simple pathogen-specific database
makeblastdb \
    -in data/pathogen_reference.fasta \
    -dbtype nucl \
    -out resources/blast_db/virus_db \
    -title "Pathogen reference DB"

# Build a multi-file human database alias
makeblastdb \
    -in "GRCh38_chr1.fna GRCh38_chr2.fna GRCh38_chr3.fna" \
    -dbtype nucl \
    -out resources/blast_db/human_db \
    -title "Human genome GRCh38"
```

---

## Outputs

| File | Path | Description |
|---|---|---|
| Verification sentinel | `results/{BUILD}/m09/{mode}/{region}/verification_complete.txt` | Pipeline completion marker |
| Per-candidate specificity report | `results/{BUILD}/m09/{mode}/{region}/specificity_report_{CANDIDATE}.md` | Markdown BLAST hit summary per candidate |
| BLAST log | `logs/{BUILD}/m09/blast_verification_{mode}_{region}.log` | Full script stderr |

### Reading the Specificity Report

Each `specificity_report_{CANDIDATE}.md` contains:

```markdown
## Candidate: flu_mp_raa_001 | Region: region_0 | Grade: A

### Forward Primer: CTGACTAARGGRATTTTRGGRTTTGTATTCA
- Top BLAST hit: Influenza A segment 7 (NC_002016.1) | E-value: 1e-12 | Identity: 98.2%
- Human genome hits: NONE
- Non-target virus hits: NONE → PASS

### Reverse Primer: GCAGTTAAACTATATARGAARCTGAARAGRG
- Top BLAST hit: Influenza A segment 7 (NC_002016.1) | E-value: 5e-11 | Identity: 96.7%
- Human genome hits: NONE
- Non-target virus hits: NONE → PASS

### Probe: CTTCTTGTTCTTGCTCGAAGACATTTCGTTCTG
- Top BLAST hit: Influenza A segment 7 (NC_002016.1) | E-value: 2e-13 | Identity: 100%
- Human genome hits: NONE → PASS

**Overall Specificity: PASS**
```

A candidate passes verification if:
1. Forward primer, reverse primer, and probe all have top hits to the target sequence
2. No significant hits (E-value `&lt;` `evalue_cutoff`) to the human genome
3. No significant hits to non-target pathogen sequences in `virus_db`

---

## Conda Environment

```
envs/primerlab.yaml
```

Key packages: `blast+`, `python`, `pandas`

---

## Known Limitations & Edge Cases

- **Remote mode rate limits**: NCBI enforces 3 BLAST API requests/second without an API key. For batches of >10 candidate sets, remote mode may take hours. Use `NCBI_API_KEY` or switch to local mode for HPC workflows.
- **Human genome database size**: A full human genome database (GRCh38) is approximately 3 GB compressed. Ensure sufficient disk space in `resources/blast_db/` before enabling `custom_db_build`.
- **Degenerate primer BLAST hits**: IUPAC degenerate bases (e.g., `R`, `Y`, `W`) in primer sequences are interpreted as ambiguities by BLAST. A primer with high degeneracy may produce fewer hits than expected. Review BLAST output manually for highly degenerate candidates.
- **Resume/skip logic**: The BLAST verification script implements a resume mechanism. If a candidate-level verification file already exists, it is skipped on re-runs. This prevents expensive re-BLAST of already-verified candidates after pipeline interruption.
- **E-value threshold sensitivity**: The default `evalue_cutoff: 0.01` is conservative. For short primers (~30 bp), an E-value of `1e-3` to `1e-5` is more appropriate for detecting meaningful off-target hits. Adjust based on primer length and acceptable specificity threshold.
