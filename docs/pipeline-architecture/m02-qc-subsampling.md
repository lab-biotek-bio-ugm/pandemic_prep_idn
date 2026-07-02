---
title: "M02: Quality Control & Subsampling"
description: "Applies hard QC filters to remove low-quality sequences and performs adaptive spatio-temporal subsampling to produce a balanced, representative dataset for phylogenetic analysis."
---

# Module M02: Quality Control & Subsampling

M02 operates on the raw FASTA and enriched metadata from M01 and produces a filtered, balanced sequence set. It executes two serial stages: (1) hard pre-filtering via a custom script, and (2) adaptive spatio-temporal subsampling via `augur filter`. The output feeds directly into M03 alignment.

---

## Input Requirements

| File | Source | Description |
|---|---|---|
| `results/{BUILD}/m01/raw_sequences.fasta` | M01 | All sequences from acquisition |
| `results/{BUILD}/m01/metadata_enriched.tsv` | M01 | Enriched metadata with clade/QC columns |

### Configuration Parameters

| Config Key | Type | Default | Description |
|---|---|---|---|
| `qc.max_n_content` | Float | `5.0` | Maximum % of ambiguous `N` bases allowed |
| `global_params.min_subtype_count` | Integer | `30` | Minimum sequences required per subtype to keep the subtype |
| `subsampling.group_by` | String | `"country year month"` | Metadata columns used to stratify sampling |
| `subsampling.max_sequences` | Integer | `3000` | Hard upper bound on total sequences passed to alignment |
| `subsampling.sequences_per_group` | Integer | — | Max sequences per stratum (overrides `max_sequences` if set) |

---

## Processing Stages

### Stage 1 — Pre-filter (`qc_prefilter`)

**Script:** `scripts/qc_prefilter.py`  
**Rule:** `qc_prefilter`  
**Environment:** `envs/nextstrain.yaml`

Performs three sequential filters on the raw FASTA:

1. **N-content filter**: Sequences exceeding `max_n_pct` (default 5%) ambiguous bases are discarded. Calculated as `count('N') / len(sequence)`.
2. **MD5 deduplication**: Exact-duplicate sequences (identical nucleotide strings) are collapsed to one representative record. Metadata record with the earliest collection date is retained.
3. **Subtype underrepresentation filter**: Subtypes with fewer than `min_subtype_count` total records after N-filtering are removed entirely. This prevents low-sample-count subtypes from introducing spurious long branches in tree inference.

**Outputs:**
- `results/{BUILD}/m02/prefiltered_sequences.fasta` — Filtered FASTA
- `results/{BUILD}/m02/qc_report.json` — Retention statistics (total input, removed per filter, output count)

```bash
# Manual invocation for debugging
python scripts/qc_prefilter.py \
    --input-fasta results/{BUILD}/m01/raw_sequences.fasta \
    --metadata results/{BUILD}/m01/metadata_enriched.tsv \
    --output-fasta results/{BUILD}/m02/prefiltered_sequences.fasta \
    --report results/{BUILD}/m02/qc_report.json \
    --min-subtype-count 30 \
    --max-n-pct 5.0
```

### Stage 2 — Adaptive Subsampling (`augur_filter`)

**Script:** `scripts/m02_filter_wrapper.py` (Snakemake `script:` directive — runs as Python)  
**Rule:** `augur_filter`  
**Environment:** `envs/nextstrain.yaml`

This rule wraps `augur filter` with adaptive logic that supports two sampling schemes:

**Mode A — Global uniform subsampling:**
Used when `subsampling` block contains only `group_by` and `max_sequences`. Draws a balanced sample globally across all defined strata.

```yaml
subsampling:
  group_by: "country year month"
  max_sequences: 3000
```

**Mode B — Focal + Contextual subsampling:**
Used when `subsampling.focal` and `subsampling.contextual` blocks are defined. Ensures that the target country/region is fully represented, with a smaller contextual sample from other regions for phylogenetic placement.

```yaml
subsampling:
  focal:
    group_by: "country year month"
    max_sequences: 500
    query: "country == 'Indonesia'"
  contextual:
    group_by: "region year"
    max_sequences: 2500
```

The wrapper script detects which mode is active from the config and constructs the appropriate `augur filter` command. If neither focal nor contextual keys are present, it falls back to global mode.

**Output:**
- `results/{BUILD}/m02/filtered_sequences.fasta`

---

## Outputs

| File | Path | Description |
|---|---|---|
| Pre-filtered FASTA | `results/{BUILD}/m02/prefiltered_sequences.fasta` | After N-content, deduplication, subtype filters |
| **Subsampled FASTA** | `results/{BUILD}/m02/filtered_sequences.fasta` | **Primary M02 output** — used by M03 alignment |
| QC report | `results/{BUILD}/m02/qc_report.json` | Filter retention statistics |
| Pre-filter log | `logs/{BUILD}/m02_qc_prefilter.log` | qc_prefilter.py stderr |
| Augur filter log | `logs/{BUILD}/m02_augur_filter.log` | augur filter stderr |

### Reading the QC Report

```json
// results/{BUILD}/m02/qc_report.json
{
  "input_count": 12847,
  "removed_high_n": 231,
  "removed_duplicates": 1563,
  "removed_low_subtype_count": 88,
  "output_count": 10965
}
```

---

## Conda Environment

```
envs/nextstrain.yaml
```

Key packages: `augur`, `pandas`, `biopython`

---

## Known Limitations & Edge Cases

- **Date-stratified subsampling requires complete dates**: `augur filter` with `group_by: "country year month"` silently drops sequences where `date` is missing or only partially known (e.g., `2024`). Check `qc_report.json` and the filter log for the count of sequences excluded due to missing date fields.
- **Subtype filter vs. focal sampling**: When using focal+contextual mode with a narrow focal region, the pre-filter's `min_subtype_count` may remove subtypes that are present only in the focal region. Lower `min_subtype_count` or disable it for single-country analyses.
- **Global `max_sequences` cap**: The `max_sequences` limit is applied globally after stratified sampling. If the number of sequences after stratification is already below `max_sequences`, no additional downsampling occurs.
- **Temporal bias**: Heavy bias toward recent sequences is common in GISAID. To counteract, include `year` in `group_by` and set `sequences_per_group` to cap per-year density.
