# Modules M01 - M03: Data Preparation

These three modules prepare raw genomic data for phylogenomic analysis by mapping metadata, applying quality control thresholds, balancing sampling, and performing multiple sequence alignment.

---

## M01: Data Acquisition & Standardisation

**Objective**: Standardize metadata and sequence records across heterogeneous repositories (GISAID, NCBI Virus, and local laboratories).

### Technical Details
*   **Metadata Standardization**: Source headers are converted to Nextstrain's standard columns: `strain`, `virus`, `gisaid_epi_isl`, `date`, `region`, `country`, `division`, `host`, `passage`, `authors`.
*   **Clade & Lineage Profiling**: Sequences are processed through Nextclade to assign unified lineages. This ensures that metadata records carry synchronized nomenclature before trees are built.

### Command Line Debugging
To run Nextclade manually for debugging:
```bash
nextclade run \
  --input-fasta data/raw_sequences.fasta \
  --dataset-name flu_h3n2_ha \
  --output-all results/debug_m01/
```

---

## M02: Quality Control & Subsampling

**Objective**: Eliminate low-quality reads and prevent spatial or temporal sampling bias from skewing tree construction.

### Technical Details
*   **Hard QC Filters**: Discards sequences that exceed the `max_n_content` (e.g., 5%) or do not meet `min_length` requirements.
*   **MD5 Deduplication**: Identical sequences are identified and collapsed to reduce down-stream computational complexity while preserving regional representation in metadata.
*   **Spatio-Temporal Subsampling**: The pipeline uses `augur filter` to balance sequence counts across groups. For example, if country A has 10,000 sequences and country B has 10, `augur filter` draws a balanced sample (e.g., max 10 per country per month) to construct a representative evolutionary history.

### Rule Excerpt
```yaml
# Snakemake subsampling configuration
subsampling:
  group_by: "country year month"
  max_sequences: 3000
```

---

## M03: Multiple Sequence Alignment (MSA)

**Objective**: Align all QC-passed sequences against a fixed reference sequence to establish a unified coordinate system.

### Technical Details
*   **Reference-Based Alignment**: Uses **MAFFT** (integrated via `augur align`) to align sequences against a predefined RefSeq (e.g., NC_002016.1 for Influenza A MP segment).
*   **Gap Preservation**: The `--fill-gaps` parameter is enforced to fill missing bases at coordinates relative to the reference sequence. This ensures that every position corresponds exactly to RefSeq coordinates across the dataset.

### Technical Commands
To check alignment performance on a local subset:
```bash
augur align \
  --sequences results/m02/filtered.fasta \
  --reference-sequence config/reference.gb \
  --output results/m03/aligned.fasta \
  --fill-gaps
```
