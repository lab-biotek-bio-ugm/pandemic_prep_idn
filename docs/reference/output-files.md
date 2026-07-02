---
title: "Output Files Catalog"
description: "Complete catalog of all files produced by each pipeline module, including their schema, location pattern, and downstream consumers."
---

# Output Files Catalog

This page catalogs every file produced by the pipeline, organised by module. All paths use the `{BUILD}` wildcard which resolves to the value of `config.build_name` (e.g., `influenza-a-mp-sea`).

---

## Module M01: Data Acquisition

| File | Path | Format | Description |
|---|---|---|---|
| Raw sequences | `results/{BUILD}/m01/raw_sequences.fasta` | FASTA | Unfiltered sequences from NCBI or local input |
| Raw metadata | `results/{BUILD}/m01/metadata.tsv` | TSV | Column-normalised metadata before enrichment |
| Nextclade output | `results/{BUILD}/m01/nextclade_output.tsv` | TSV | Per-sequence Nextclade clade/lineage/QC output |
| **Enriched metadata** | `results/{BUILD}/m01/metadata_enriched.tsv` | TSV | **Primary M01 output.** Sanitised metadata + Nextclade clade/lineage columns. Used by all downstream modules. |
| Nextclade dataset | `results/{BUILD}/m01/nextclade_dataset.zip` | ZIP | Downloaded Nextclade reference dataset (temp) |

### `metadata_enriched.tsv` Key Columns

| Column | Type | Description |
|---|---|---|
| `strain` | String | Unique sequence identifier (primary join key) |
| `virus` | String | Organism/taxon string |
| `date` | String | Collection date (ISO 8601 or partial) |
| `country` | String | Country of collection |
| `region` | String | Sub-national region |
| `subtype` | String | Influenza subtype (e.g., `H3N2`) |
| `clade` | String | Nextclade clade assignment |
| `lineage` | String | Unified lineage designation |
| `qcOverallStatus` | String | Nextclade QC status: `good` / `mediocre` / `bad` |

---

## Module M02: Quality Control & Subsampling

| File | Path | Format | Description |
|---|---|---|---|
| Pre-filtered FASTA | `results/{BUILD}/m02/prefiltered_sequences.fasta` | FASTA | After N-content, dedup, and subtype filters |
| **Filtered FASTA** | `results/{BUILD}/m02/filtered_sequences.fasta` | FASTA | **Primary M02 output.** Spatio-temporally balanced subset. |
| QC report | `results/{BUILD}/m02/qc_report.json` | JSON | Filter statistics (input count, removed per filter, output count) |

### `qc_report.json` Schema

```json
{
  "input_count": 12847,
  "removed_high_n": 231,
  "removed_duplicates": 1563,
  "removed_low_subtype_count": 88,
  "output_count": 10965
}
```

---

## Module M03: Alignment

| File | Path | Format | Description |
|---|---|---|---|
| **Aligned FASTA** | `results/{BUILD}/m03/aligned.fasta` | FASTA | **Primary M03 output.** Reference-coordinate MSA used by M04, M05, M06, and M07. |

All sequences in `aligned.fasta` have identical length equal to the reference genome length. Gap characters (`-`) represent deletions or assembly gaps relative to the reference; positions filled with `N` (via `--fill-gaps`) represent missing sequence at sequence termini.

---

## Module M04: Phylogenetics

| File | Path | Format | Description |
|---|---|---|---|
| Raw ML tree | `results/{BUILD}/m04/tree_raw.nwk` | Newick | Unrooted maximum-likelihood topology from IQ-TREE2 |
| **Time-calibrated tree** | `results/{BUILD}/m04/tree.nwk` | Newick | **Primary M04 output.** Rooted, time-resolved phylogeny. |
| Branch-length data | `results/{BUILD}/m04/branch_lengths.json` | JSON | Per-node clock data (numdate, mutation lengths) for Auspice |

### `branch_lengths.json` Schema

```json
{
  "nodes": {
    "A/Indonesia/PDDV0001/2024": {
      "numdate": 2024.673,
      "clock_length": 0.00341,
      "mutation_length": 0.00289,
      "raw_date": "2024-09-02"
    },
    "NODE_0001234": {
      "numdate": 2020.15,
      "clock_length": 0.00189,
      "mutation_length": 0.00201
    }
  }
}
```

---

## Module M05: Selection Analysis

Files are generated per gene (e.g., `M1`, `M2` for Flu A MP segment).

| File | Path | Format | Description |
|---|---|---|---|
| HyPhy input MSA | `results/{BUILD}/m05/aligned_hyphy.fasta` | FASTA | Diversity-subsampled alignment for HyPhy |
| Subsample report | `results/{BUILD}/m05/hyphy_subsample_report.json` | JSON | Subsampling statistics |
| Gene FASTA | `results/{BUILD}/m05/cleaned_{GENE}_NT.fasta` | FASTA | Codon-sanitised per-gene alignment |
| Pruned tree | `results/{BUILD}/m05/tree_{GENE}.nwk` | Newick | Tip-pruned tree for each gene |
| **FUBAR results** | `results/{BUILD}/m05/fubar_{GENE}.json` | JSON | Per-codon purifying/positive selection probabilities |
| **MEME results** | `results/{BUILD}/m05/meme_{GENE}.json` | JSON | Per-codon episodic selection p-values |
| GARD results | `results/{BUILD}/m05/gard_{GENE}.json` | JSON | Recombination breakpoint positions (if enabled) |

### `fubar_{GENE}.json` Key Fields

```json
{
  "MLE": {
    "headers": ["alpha", "beta", "beta-alpha", "Prob[alpha>beta]", "Prob[beta>alpha]"],
    "content": {
      "0": [[0.123, 0.001, -0.122, 0.97, 0.03]]
    }
  }
}
```

`Prob[alpha>beta]` (index 3) is the **purifying selection posterior** consumed by M07 as the $P$ score.

### `meme_{GENE}.json` Key Fields

```json
{
  "MLE": {
    "headers": ["alpha", "beta-", "p-", "beta+", "p+", "LRT", "p-value", "# branches"],
    "content": {
      "0": [[0.12, 0.0, 1.0, 0.89, 0.23, 4.67, 0.031, 2]]
    }
  }
}
```

`p-value` (index 6) is compared against `hyphy.p_value_threshold`. Sites below the threshold trigger the $E$ penalty in M07.

---

## Module M06: Ancestral Reconstruction & Export

| File | Path | Format | Description |
|---|---|---|---|
| **NT mutation data** | `results/{BUILD}/m06/nt_muts.json` | JSON | Per-branch nucleotide mutation lists across all tree nodes |
| **AA mutation data** | `results/{BUILD}/m06/aa_muts.json` | JSON | Per-branch amino acid mutation lists per annotated gene |
| **Auspice dataset** | `results/{BUILD}/nextstrain/auspice_dataset.json` | JSON | Complete Nextstrain/Auspice v2 visualisation dataset |

### `nt_muts.json` Schema

```json
{
  "nodes": {
    "A/Indonesia/PDDV0001/2024": {
      "muts": ["A187T", "G452A"],
      "aa_muts": {"M1": ["I297V"]}
    }
  }
}
```

---

## Module M07: Bridge Builder

| File | Path | Format | Description |
|---|---|---|---|
| **Bridge report** | `results/{BUILD}/m07/bridge_report.json` | JSON | **Primary M07 output.** Ranked candidate regions with composite scores and grade classifications. |
| **Primer manifest** | `results/{BUILD}/m07/primer_targets_manifest.yaml` | YAML | Input manifest for M08 PrimerLab — specifies coordinates, consensus sequence, and exclusion zones per design region. |
| Conserved regions BED | `results/{BUILD}/m07/conserved_regions.bed` | BED | Genomic browser track of all candidate regions (loadable in IGV or UCSC). |
| Entropy profile | `results/{BUILD}/m07/entropy_profile.tsv` | TSV | Per-nucleotide Shannon entropy, base frequencies, IUPAC consensus. |
| Active config | `results/{BUILD}/m07/active_config.yaml` | YAML | Snapshot of merged Snakemake config used for this run (for reproducibility). |

### `bridge_report.json` Schema

```json
{
  "build_name": "influenza-a-mp-sea",
  "regions": [
    {
      "region_id": "region_0",
      "start": 187,
      "end": 218,
      "length": 31,
      "gene": "M1",
      "avg_conservation": 0.985,
      "avg_purifying_prob": 0.921,
      "meme_sites": 0,
      "near_recombination": false,
      "composite_score": 0.968,
      "grade": "A",
      "tier": "strict",
      "coding": true
    }
  ]
}
```

### `entropy_profile.tsv` Schema

| Column | Description |
|---|---|
| `position` | 1-based reference coordinate |
| `consensus` | Most frequent nucleotide |
| `iupac` | IUPAC ambiguity code if MAF > `iupac_threshold` |
| `entropy` | Shannon entropy $H$ (0 = fully conserved, log₂(5) ≈ 2.32 = fully random) |
| `freq_A` / `freq_T` / `freq_G` / `freq_C` / `freq_gap` | Per-nucleotide frequency |

---

## Module M08: Primer Design

Files are organised by `{mode}` (e.g., `RAA`, `PCR`) and `{region}` (e.g., `region_0`).

| File | Path | Format | Description |
|---|---|---|---|
| PrimerLab config | `results/{BUILD}/m08/{mode}/primerlab_config_{region}.yaml` | YAML | Auto-generated PrimerLab input config for each region |
| Install sentinel | `results/{BUILD}/m08/.primerlab_installed` | TXT | Marker confirming PrimerLab was installed successfully |
| Design sentinel | `results/{BUILD}/m08/{mode}/{region}/.design_complete` | TXT | Marker confirming design run completed |
| **Alternatives** | `results/{BUILD}/m08/{mode}/{region}/alternatives.json` | JSON | All generated candidate primer sets with sequences and thermodynamic properties |
| **Summary** | `results/{BUILD}/m08/{mode}/{region}/summary.json` | JSON | Ranked top candidates after all PrimerLab filters |

### `summary.json` Schema (PrimerLab output)

```json
{
  "candidates": [
    {
      "candidate_id": "flu_mp_raa_001",
      "region_id": "region_0",
      "forward_primer": "CTGACTAARGGRATTTTRGGRTTTGTATTCA",
      "reverse_primer": "GCAGTTAAACTATATARGAARCTGAARAGRG",
      "probe": "CTTCTTGTTCTTGCTCGAAGACATTTCGTTCTG",
      "multiplier": 8,
      "delta_G_hairpin": -1.2,
      "delta_G_homodimer": -3.1,
      "delta_G_heterodimer": -4.5,
      "amplicon_length": 187,
      "gc_forward": 45.2,
      "gc_reverse": 42.8,
      "rank": 1
    }
  ]
}
```

---

## Module M09: BLAST Specificity

| File | Path | Format | Description |
|---|---|---|---|
| Verification sentinel | `results/{BUILD}/m09/{mode}/{region}/verification_complete.txt` | TXT | Pipeline completion marker |
| **Specificity reports** | `results/{BUILD}/m09/{mode}/{region}/specificity_report_{ID}.md` | Markdown | Per-candidate BLAST hit summary with PASS/FAIL verdict |
| Final candidates | `results/{BUILD}/m09/{mode}/{region}/final_candidates.json` | JSON | BLAST-verified candidates that passed all specificity checks |

---

## Module M10: In-Silico Validation

| File | Path | Format | Description |
|---|---|---|---|
| **Coverage summary TSV** | `results/{BUILD}/m10/{mode}/coverage_summary.tsv` | TSV | Per-candidate, per-year, per-subtype coverage percentages |
| **Coverage summary JSON** | `results/{BUILD}/m10/{mode}/coverage_summary.json` | JSON | Structured version of the coverage summary |
| **Coverage report** | `results/{BUILD}/m10/{mode}/coverage_report.md` | Markdown | Human-readable report with coverage table and interpretation |
| **Coverage heatmap** | `results/{BUILD}/m10/{mode}/coverage_heatmap.html` | HTML | Interactive heatmap of temporal × subtype coverage |
| GFF3 track | `results/{BUILD}/m10/{mode}/m10_candidates.gff3` | GFF3 | Genomic feature file for all validated candidates (loadable in IGV) |

### `coverage_summary.tsv` Schema

| Column | Description |
|---|---|
| `candidate_id` | Candidate identifier from M08/M09 |
| `year` | Collection year |
| `subtype` | Viral subtype (e.g., `H3N2`, `H1N1`) |
| `total_sequences` | Total sequences evaluated for this year/subtype |
| `amplified` | Sequences where both primers matched within mismatch tolerance |
| `coverage_pct` | `amplified / total_sequences × 100` |
| `3prime_failures` | Sequences failing due to 3' mismatch criterion |
| `internal_failures` | Sequences failing due to internal mismatch count |

---

## Nextstrain Visualisation

| File | Path | Format | Description |
|---|---|---|---|
| **Auspice JSON** | `results/{BUILD}/nextstrain/auspice_dataset.json` | JSON | Complete Auspice v2 dataset for interactive phylogenetic visualisation |

Load with: `auspice view --datasetDir results/{BUILD}/nextstrain/`

---

## Log Files

All rule logs are written to `logs/{BUILD}/`:

| Log | Rule |
|---|---|
| `logs/{BUILD}/m01_download.log` | Data acquisition |
| `logs/{BUILD}/m02_qc_prefilter.log` | QC pre-filter |
| `logs/{BUILD}/m02_augur_filter.log` | Subsampling |
| `logs/{BUILD}/m03_augur_align.log` | MAFFT alignment |
| `logs/{BUILD}/m04_augur_tree.log` | IQ-TREE2 |
| `logs/{BUILD}/m04_augur_refine.log` | TreeTime clock fitting |
| `logs/{BUILD}/m05_hyphy_subsample.log` | HyPhy subsampling |
| `logs/{BUILD}/m05_hyphy_fubar_{GENE}.log` | FUBAR per gene |
| `logs/{BUILD}/m05_hyphy_meme_{GENE}.log` | MEME per gene |
| `logs/{BUILD}/m05_hyphy_gard_{GENE}.log` | GARD per gene |
| `logs/{BUILD}/m06_augur_ancestral.log` | Ancestral reconstruction |
| `logs/{BUILD}/m06_augur_translate.log` | AA translation |
| `logs/{BUILD}/m06_augur_export.log` | Auspice export |
| `logs/{BUILD}/m07_bridge.log` | Bridge Builder |
| `logs/{BUILD}/m08/generate_payload_{mode}_{region}.log` | PrimerLab payload generation |
| `logs/{BUILD}/m09/blast_verification_{mode}_{region}.log` | BLAST verification |
| `logs/{BUILD}/m10/insilico_pcr_{mode}.log` | In-silico PCR |
