---
title: "Data Management & GISAID Compliance"
description: "Best practices for managing viral genomic data downloaded from GISAID and NCBI Virus, ensuring database compliance, and tracking reference data changes."
---

# Data Management & GISAID Compliance

This pipeline operates on large genomic datasets sourced from public repositories, principally **GISAID (Global Initiative on Sharing All Influenza Data)** and **NCBI GenBank**. Managing these datasets requires careful adherence to data access agreements, data integrity validation, and efficient local storage practices.

---

## 1. GISAID Data Access & Compliance

All data downloaded from GISAID (via the EpiFlu or EpiCoV portals) is subject to the **GISAID Database Access Agreement (DAA)**. 

### Core DAA Requirements:
- **No Redistribution**: You must **never** upload raw GISAID FASTA sequences or full metadata CSV/TSV files to public code repositories (e.g., GitHub, GitLab) or public data sharing services.
- **Acknowledge Lab Contributions**: Publications using data from GISAID must acknowledge the Originating and Submitting laboratories. The pipeline does not store these lists automatically, so keep a copy of the downloaded acknowledgment table.
- **Fair Use**: GISAID data should be accessed only through authorized credentials. Sharing accounts or API access keys is prohibited.

### Compliance in the Codebase:
- The `.gitignore` file of this repository is configured to exclude `dataset/` and `results/` folders where GISAID files are stored.
- **DO NOT** modify the `.gitignore` to track sequence FASTAs or metadata CSVs.
- Only share configuration profiles (`config/*.yaml`) containing public parameters, paths, and queries — never the sequences themselves.

---

## 2. Managing NCBI GenBank Data

GenBank data (accessed via NCBI Virus or the NCBI Datasets API) is public domain. However, large downloads can consume significant network bandwidth and local disk space.

### NCBI API Key Configuration:
When running with `database: ncbi`, the pipeline invokes `ncbi-datasets-cli` to pull records. NCBI throttles anonymous requests to 3 requests per second. To prevent timeout failures during large downloads, get an API key and set it in your environment:

```bash
# Get an API key from: https://www.ncbi.nlm.nih.gov/account/settings/
export NCBI_API_KEY="your_ncbi_api_key_here"

# Verify the key is active in your terminal session
echo $NCBI_API_KEY
```

---

## 3. Data Integrity & Verification

Large sequence files are prone to truncation during download or extraction. Always verify the checksums of your local datasets before initiating a run:

```bash
# Calculate MD5 checksum of input FASTA
md5sum dataset/EpiFlu/MP/EpiFlu_A_MP_Combined_Sequence.fasta > dataset/EpiFlu/MP/sequences.md5

# Verify in subsequent runs or cluster transfers
md5sum -c dataset/EpiFlu/MP/sequences.md5
```

### Checking for Header Mismatches:
A common cause of M01 enrichment failure is a mismatch between the FASTA sequence headers and the metadata `strain` (or accession) identifiers. Run this quick diagnostic script before launching:

```bash
# Extract FASTA headers (removing the '>' symbol)
grep "^>" dataset/EpiFlu/MP/EpiFlu_A_MP_Combined_Sequence.fasta | sed 's/>//' | sort > fasta_headers.txt

# Extract Metadata Strain names (assuming column 1, ignoring header row)
awk -F',' 'NR>1 {print $1}' dataset/EpiFlu/MP/EpiFlu_A_MP_Combined_Metadata.csv | sort > metadata_strains.txt

# Compare the lists
comm -23 fasta_headers.txt metadata_strains.txt > unmatched_sequences.txt
echo "Unmatched sequences (FASTA has it, Metadata doesn't): $(wc -l < unmatched_sequences.txt)"
```

If the count is high, check if your metadata format matches `local_input.format` in the config profile.

---

## 4. Local Database Directory Conventions

To keep the workspace clean and compliant, follow these directory conventions for datasets and references:

```
pandemic-preparedness/
├── dataset/                     # Large raw sequences & metadata (gitignored)
│   ├── EpiFlu/                  # Influenza data
│   │   ├── MP/
│   │   └── NP/
│   └── GISAID_SARS2/            # SARS-CoV-2 data
├── config/references/           # Small, public reference files (tracked in git)
│   ├── flu_a_mp_reference.fasta # Reference genome
│   └── flu_a_mp_annotation.gff  # GFF3 annotation
└── resources/blast_db/          # Local BLAST databases (gitignored)
    ├── build_complete.txt       # Build sentinel
    ├── human_db.*               # Indexed human chromosomes
    └── virus_db.*               # Indexed respiratory viruses
```

Keep the `dataset/` folder strictly outside the code tracks, and never commit clinical sample FASTAs or raw patient metadata to the repository.
