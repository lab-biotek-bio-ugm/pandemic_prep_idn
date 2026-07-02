---
title: "Configuration Reference"
description: "Complete reference for all configuration keys in the pipeline YAML profiles. Covers all modules from M01 data acquisition through M10 in-silico validation."
---

# Configuration Reference

This page documents every configuration key accepted by the pipeline. All settings are defined in a YAML profile file (e.g., `config/config_influenza_a_mp.yaml`) and passed to Snakemake via `--configfile`.

---

## Global Parameters

```yaml
global_params:
  min_subtype_count: 30
  iupac_threshold: 0.15
```

| Key | Type | Default | Description |
|---|---|---|---|
| `global_params.min_subtype_count` | Integer | `30` | Minimum number of sequences a subtype must have across the entire dataset to be included in any downstream analysis. Subtypes below this threshold are removed in M02 pre-filter and collapsed to `"Other (rare)"` in M10 reports. |
| `global_params.iupac_threshold` | Float | `0.15` | Minor allele frequency (MAF) threshold that triggers IUPAC degenerate base encoding in the consensus sequence. At position i, if any non-consensus allele exceeds this frequency, the position is encoded with the appropriate IUPAC code (e.g., `R` for A/G). |

---

## Pipeline Identity

| Key | Type | Example | Description |
|---|---|---|---|
| `build_name` | String | `"influenza-a-mp-sea"` | Unique identifier for this run. All outputs are written to `results/{build_name}/`. Use distinct names for each pathogen/segment/region combination. |
| `database` | String | `"local"` | Data acquisition mode. `"local"` reads pre-downloaded files; `"ncbi"` fetches via `ncbi-datasets-cli` API. |

---

## References

```yaml
references:
  reference_id: "NC_002016.1"
  genome: "config/references/flu_a_mp_reference.fasta"
  annotation: "config/references/flu_a_mp_annotation.gff"
```

| Key | Type | Description |
|---|---|---|
| `references.reference_id` | String | NCBI accession number of the reference genome. Used by M09 BLAST to identify the correct target hit and by M05 `extract_genes.py` as the reference tip ID in the alignment. |
| `references.genome` | String | Path to reference genome (FASTA or GenBank `.gb`). Used by `augur align` (M03) and `m07_bridge_builder.py`. |
| `references.annotation` | String | Path to GFF3 annotation file. Used by M05 for gene slicing and by M06 `augur translate` for amino acid annotation. |

---

## M01: Data Acquisition

### Local Input Pathway (`database: local`)

```yaml
local_input:
  sequences: "dataset/EpiFlu/MP/EpiFlu_A_MP_Combined_Sequence.fasta"
  metadata:  "dataset/EpiFlu/MP/EpiFlu_A_MP_Combined_Metadata.csv"
  format:    "gisaid_epiflu"
```

| Key | Type | Description |
|---|---|---|
| `local_input.sequences` | String | Path to input FASTA file. |
| `local_input.metadata` | String | Path to metadata file. |
| `local_input.format` | String | Parser mode. Options: `gisaid_epiflu`, `gisaid_tsv`, `gisaid_raw`, `ncbi_virus_csv`, `tsv`, `auto`. |

### NCBI API Pathway (`database: ncbi`)

```yaml
ncbi_query:
  organism: "Influenza A virus"
  location: "Indonesia"
  released_after: "2020-01-01"
```

| Key | Type | Description |
|---|---|---|
| `ncbi_query.organism` | String | NCBI taxon name for the sequence query. |
| `ncbi_query.location` | String | Geographic location filter string. |
| `ncbi_query.released_after` | String | Optional ISO date to filter sequences released after this date. |

### Early Filter (pre-Nextclade)

```yaml
early_filter:
  enabled: true
  query: null           # e.g. "subtype == 'A / H1N1'"
  min_date: null
  max_date: "2026-04-01"
```

| Key | Type | Default | Description |
|---|---|---|---|
| `early_filter.enabled` | Boolean | `false` | Activate pre-Nextclade filtering. |
| `early_filter.query` | String/null | `null` | Pandas query string applied to metadata before Nextclade classification. |
| `early_filter.min_date` | String/null | `null` | Discard sequences collected before this ISO date. |
| `early_filter.max_date` | String/null | `null` | Discard sequences collected after this ISO date. |

### Nextclade

```yaml
nextclade:
  enabled: false
  dataset_name: "nextstrain/flu/h1n1pdm/ha/MW626062"
```

| Key | Type | Default | Description |
|---|---|---|---|
| `nextclade.enabled` | Boolean | `true` | Run Nextclade for clade/lineage assignment. Disable if the reference dataset is not available or for non-standard pathogens. |
| `nextclade.dataset_name` | String | — | Nextclade dataset name (e.g., `flu_h3n2_ha`, `sars-cov-2`). Run `nextclade dataset list` to see available datasets. |

---

## M02: Quality Control & Subsampling

```yaml
qc:
  max_n_content: 5.0
  dedup: true

subsampling:
  group_by: "subtype country year month"
  sequences_per_group: 10
  max_sequences: null
  min_length: 800
  seed: 42
  query: "country in ['Indonesia', 'Singapore', ...]"

  focal:
    enabled: false
    query: "country == 'Indonesia'"
    max_sequences: null
    sequences_per_group: 50
    group_by: "subtype year month"

  contextual:
    query: "country != 'Indonesia'"
    max_sequences: null
    sequences_per_group: 3
    group_by: "subtype country year month"
```

| Key | Type | Default | Description |
|---|---|---|---|
| `qc.max_n_content` | Float | `5.0` | Maximum percentage of ambiguous `N` bases tolerated per sequence. |
| `qc.dedup` | Boolean | `true` | Enable MD5-based exact-duplicate removal. |
| `subsampling.group_by` | String | `"country year month"` | Space-separated metadata columns used to stratify the balanced sample. |
| `subsampling.sequences_per_group` | Integer | — | Maximum sequences per stratum. Mutually exclusive with `max_sequences`. |
| `subsampling.max_sequences` | Integer/null | `null` | Hard upper bound on total sequences. Applied after stratification. |
| `subsampling.min_length` | Integer | `800` | Minimum sequence length; shorter sequences are discarded by `augur filter`. |
| `subsampling.seed` | Integer | `42` | Random seed for reproducibility of subsampling. |
| `subsampling.query` | String/null | `null` | Pandas query string for pre-filtering before subsampling (e.g., restrict to SEA countries). |
| `subsampling.focal.enabled` | Boolean | `false` | Enable focal+contextual sampling mode. |
| `subsampling.focal.query` | String | — | Pandas query defining the focal (high-density) region. |
| `subsampling.focal.sequences_per_group` | Integer | — | Max sequences per group for focal region. |
| `subsampling.contextual.query` | String | — | Pandas query defining the contextual (background) sequences. |
| `subsampling.contextual.sequences_per_group` | Integer | — | Max sequences per group for contextual region. |

---

## M04: Phylogenetic Inference

```yaml
augur:
  threads: 16

timetree:
  clock_rate: 0.003
  clock_std_dev: 0.0015
```

| Key | Type | Default | Description |
|---|---|---|---|
| `augur.threads` | Integer | `8` | CPU threads for IQ-TREE2 (M04) and MAFFT (M03). |
| `timetree.clock_rate` | Float | — | Prior substitution rate (substitutions/site/year). **Pathogen-specific — do not use defaults.** |
| `timetree.clock_std_dev` | Float | — | Allowed deviation from `clock_rate`. Typical range: `10–30%` of `clock_rate`. |

---

## M05: HyPhy Selection Analysis

```yaml
hyphy:
  threads: 8
  p_value_threshold: 0.1
  methods:
    - "FUBAR"
    - "MEME"
  target_genes:
    - "ALL"
  subsample:
    enabled: true
    max_sequences: 2000
    min_sequences: 50
    strategy: "balanced"
    group_by:
      subtype: true
      country: true
      year: true
    seed: 42
```

| Key | Type | Default | Description |
|---|---|---|---|
| `hyphy.threads` | Integer | `32` | CPU threads per HyPhy method. Note: HyPhy scales poorly beyond 8–12 threads for short alignments; more threads allow more parallel jobs. |
| `hyphy.p_value_threshold` | Float | `0.1` | MEME p-value threshold below which a site is considered under episodic selection. |
| `hyphy.methods` | List | `["FUBAR"]` | HyPhy methods to run. Options: `FUBAR`, `MEME`, `GARD`. Add `GARD` only for non-segmented viruses. |
| `hyphy.target_genes` | List | `["ALL"]` | Genes to process. `ALL` uses all CDS genes detected from the GFF3. Specify gene names (e.g., `["M1", "M2"]`) to restrict. |
| `hyphy.subsample.enabled` | Boolean | `false` | Enable diversity-aware subsampling before HyPhy. Recommended for datasets > 1,000 sequences. |
| `hyphy.subsample.max_sequences` | Integer | `300` | Maximum sequences for HyPhy input. Recommended ≤ 500 for tractable runtime. |
| `hyphy.subsample.min_sequences` | Integer | `50` | Minimum sequences required. Pipeline aborts if fewer than this remain after subsampling. |
| `hyphy.subsample.strategy` | String | `"balanced"` | Subsampling strategy. `"balanced"` = stratified adaptive; `"random"` = random draw. |
| `hyphy.subsample.group_by.subtype` | Boolean | `true` | Include subtype as a stratification dimension. |
| `hyphy.subsample.group_by.country` | Boolean | `true` | Include country as a stratification dimension. |
| `hyphy.subsample.group_by.year` | Boolean | `true` | Include collection year as a stratification dimension. |
| `hyphy.subsample.seed` | Integer | `42` | Random seed for subsample reproducibility. |

---

## M07: Bridge Builder

```yaml
bridge:
  strategy: "full_segment"
  window_size: 300
  min_region_length: 150
  overlap_threshold: 0.75
  buffer_bp: 50
  entropy_threshold: 0.020
  gard_buffer_bp: 50
  include_utr: true
  design_regions:
    - "segment-MP"
  scoring:
    core:
      conservation_weight: 0.66
      purifying_weight: 0.34
    multipliers:
      meme:
        thresholds:
          - max_sites: 0
            multiplier: 1.0
          - max_sites: 1
            multiplier: 0.9
          - max_sites: 3
            multiplier: 0.7
        veto_above: 3
      recombination:
        near_penalty: 0.5
      gc_content:
        ideal_min: 35
        ideal_max: 60
        marginal_ranges:
          - {min: 25, max: 34}
          - {min: 61, max: 70}
        marginal_multiplier: 0.7
        veto_below: 25
        veto_above: 70
      conserved_block:
        entropy_threshold: 0.10
        ideal_min_len: 35
        penalty_min_len: 20
        penalty_multiplier: 0.8
      homopolymer:
        max_repeats: 5
        penalty_multiplier: 0.6
  tiers:
    strict:
      min_conservation: 0.95
      min_purifying_prob: 0.9
    moderate:
      min_conservation: 0.85
      min_purifying_prob: 0.7
    survival:
      min_conservation: 0.80
```

| Key | Type | Default | Description |
|---|---|---|---|
| `bridge.strategy` | String | `"windowed"` | Scanning mode. `"full_segment"` treats the entire segment as one design target; `"windowed"` uses a sliding window. |
| `bridge.window_size` | Integer | `300` | Sliding window size in bp. Ignored when `strategy: full_segment`. |
| `bridge.min_region_length` | Integer | `150` | Minimum length (bp) of a candidate region to be included in the manifest. |
| `bridge.entropy_threshold` | Float | `0.020` | Per-position Shannon entropy cap. Positions exceeding this are excluded from primer binding sites. |
| `bridge.include_utr` | Boolean | `true` | Include UTR (untranslated region) positions in the conservation scan. |
| `bridge.design_regions` | List | `null` | Restrict manifest output to specific named regions. `null` exports all regions. |
| `bridge.scoring.core.conservation_weight` | Float | `0.66` | Weight of physical sequence conservation ($C$) in the composite score. |
| `bridge.scoring.core.purifying_weight` | Float | `0.34` | Weight of FUBAR purifying selection probability ($P$) in the composite score. |
| `bridge.scoring.multipliers.meme.veto_above` | Integer | `3` | Veto (discard) positions where more than this many MEME-detected episodic selection sites fall within the window. |
| `bridge.scoring.multipliers.gc_content.ideal_min` | Integer | `35` | Lower bound of ideal GC content range for primer design. |
| `bridge.scoring.multipliers.gc_content.veto_below` | Integer | `25` | GC content below this value vetoes the candidate regardless of other scores. |
| `bridge.scoring.multipliers.homopolymer.max_repeats` | Integer | `5` | Maximum consecutive identical nucleotides before a homopolymer penalty is applied. |
| `bridge.tiers.strict.min_conservation` | Float | `0.95` | Minimum conservation score for a region to qualify as a Platinum (Grade A) target. |
| `bridge.tiers.moderate.min_conservation` | Float | `0.85` | Minimum conservation score for Gold (Grade B) classification. |
| `bridge.tiers.survival.min_conservation` | Float | `0.80` | Minimum conservation score for Silver (Grade C) — minimum viable target. |

---

## M08: PrimerLab Integration

```yaml
primerlab:
  repo_path: "git+https://github.com/CORE-Lab-Research/primerlab-genomic.git@main"
  template: "config/primerlab_template.yaml"
  design:
    num_candidates: 1000
    search_strategy: "standard"
    hard_qc_filter: "false"
    max_search_iterations: 5
    vienna_ranking_limit: 500
    window_size: null
    overlap: 200
    validate: true
    entropy_threshold: 0.50
    min_exclusion_block_length: 10
    degenerate:
      enabled: true
      maf_threshold: 0.10
      max_multiplier: 16
      thermo_eval_mode: "worst_case"
      vienna_temp: 37.0
      target_coverage_increase: 0.05
```

| Key | Type | Default | Description |
|---|---|---|---|
| `primerlab.repo_path` | String | — | Installation source for PrimerLab. Use a local path (e.g., `E:/Project/primerlab-genomic` as a development path example) for local testing, or the official remote git URL (e.g., `git+https://github.com/CORE-Lab-Research/primerlab-genomic.git@main`) for production. Set to `"conda"` if pre-installed. |
| `primerlab.template` | String | — | Path to the PrimerLab base template YAML (`config/primerlab_template.yaml`). This defines the reaction type (RAA/PCR) and thermodynamic parameters. |
| `primerlab.design.num_candidates` | Integer | `1000` | Total candidate primer sets to generate via Primer3 before filtering. |
| `primerlab.design.vienna_ranking_limit` | Integer | `500` | Maximum candidates to pass through ViennaRNA secondary structure analysis. Higher values = more thorough but slower. |
| `primerlab.design.window_size` | Integer/null | `null` | Sliding window size for Primer3 search. `null` = process the entire design region as one block (recommended for short segments ≤ 2000 bp). |
| `primerlab.design.validate` | Boolean | `true` | Run in-silico PCR and probe binding validation within PrimerLab during design. |
| `primerlab.design.entropy_threshold` | Float | `0.50` | Shannon entropy threshold (0–log₂(5) ≈ 2.32) above which positions are excluded from primer binding sites. `0.50` corresponds to ~30% MAF. |
| `primerlab.design.min_exclusion_block_length` | Integer | `10` | Minimum consecutive excluded positions to form an exclusion block. Prevents single wobble-base positions from blocking all primer placement. |
| `primerlab.design.degenerate.enabled` | Boolean | `true` | Enable degenerate primer design using IUPAC ambiguity codes. |
| `primerlab.design.degenerate.maf_threshold` | Float | `0.10` | MAF threshold for introducing a degenerate base at a variable position. |
| `primerlab.design.degenerate.max_multiplier` | Integer | `16` | Maximum degeneracy multiplier (total sequence combinations in the pool). Primers exceeding this limit are rejected. |
| `primerlab.design.degenerate.thermo_eval_mode` | String | `"worst_case"` | How to evaluate thermodynamics for a degenerate pool: `"worst_case"` = most stringent; `"average"` = pool mean; `"best_case"` = least stringent. |
| `primerlab.design.degenerate.vienna_temp` | Float | `37.0` | Temperature (°C) for ViennaRNA structure folding. Use 37°C for RAA/RPA, 55–60°C for PCR. |

---

## M09: BLAST Specificity Verification

```yaml
blast_verification:
  enabled: true
  num_candidates: 100
  mode: "local"
  database: "custom_respiratory_db"
  local_db_path: "resources/blast_db"
  qcov_hsp_perc: 60
  min_coverage_percent: 30.0
  three_prime_critical_bp: 8
  custom_db_build:
    enabled: true
    output_dir: "resources/blast_db"
    max_seqs_per_taxon: 500
    seq_length_min: 500
    cleanup_intermediates: true
    hosts:
      - name: "T2T_CHM13"
        fasta_url: "https://..."
        enabled: true
    virus_taxons:
      - name: "Influenza A virus"
        taxid: 11320
```

| Key | Type | Default | Description |
|---|---|---|---|
| `blast_verification.mode` | String | `"remote"` | `"remote"` queries NCBI `nt`; `"local"` queries custom local databases. |
| `blast_verification.num_candidates` | Integer | From `vienna_ranking_limit` | Number of top candidates to subject to BLAST verification. |
| `blast_verification.qcov_hsp_perc` | Integer | `60` | Minimum query coverage (%) for a BLAST HSP to be reported. Filters out short local alignments. |
| `blast_verification.three_prime_critical_bp` | Integer | `8` | Number of 3' terminal bases considered critical. Mismatches here are weighted more heavily. |
| `blast_verification.custom_db_build.max_seqs_per_taxon` | Integer | `500` | Maximum sequences downloaded from NCBI per virus taxon when building the local virus database. |
| `blast_verification.custom_db_build.cleanup_intermediates` | Boolean | `true` | Delete intermediate FASTA and .gz files after database indexing. Saves ~12–15 GB disk space. |

---

## M10: In-Silico PCR

```yaml
insilico_pcr:
  enabled: true
  max_mismatch: 3
  num_candidates: 20
  min_subtype_seqs: 30
  raa_3prime_bp: 8
  raa_max_3prime_mm: 1
  raa_max_internal_mm: 4
```

| Key | Type | Default | Description |
|---|---|---|---|
| `insilico_pcr.max_mismatch` | Integer | `3` | Maximum total mismatches tolerated between primer and template during `seqkit amplicon` simulation. |
| `insilico_pcr.num_candidates` | Integer | `10` | Number of top BLAST-verified candidates to simulate. |
| `insilico_pcr.min_subtype_seqs` | Integer | `30` | Minimum sequences per subtype in the validation database. Subtypes below this threshold are collapsed to `"Other (rare)"` in coverage reports. |
| `insilico_pcr.raa_3prime_bp` | Integer | `8` | Number of 3' terminal bases subject to the strict 3' mismatch rule (recombinase extension region). |
| `insilico_pcr.raa_max_3prime_mm` | Integer | `1` | Maximum mismatches tolerated in the 3' critical region (`raa_3prime_bp` bases). Exceeding this is a hard rejection for RAA assays. |
| `insilico_pcr.raa_max_internal_mm` | Integer | `4` | Maximum mismatches tolerated in the 5' internal region (primer length minus `raa_3prime_bp`). |
