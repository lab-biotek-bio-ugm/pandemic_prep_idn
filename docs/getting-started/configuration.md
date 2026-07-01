# Configuration Guide

This guide explains the parameters within the YAML configuration files (e.g., `config/config_influenza_a_mp.yaml`). Properly calibrating these values is essential for adapting the pipeline to new pathogens or surveillance targets.

---

## Configuration File Structure

The pipeline uses standard YAML profiles. Below is a comprehensive breakdown of the core configuration blocks:

### 1. Build Metadata & Paths

*   `build_name`: A unique string identifier for the run. All outputs are placed under `results/{build_name}/`.
*   `local_input`:
    *   `sequences`: Path to the input FASTA file containing surveillance sequences.
    *   `metadata`: Path to the matching CSV/TSV metadata file.
    *   `format`: Parser layout for metadata headers (e.g., `gisaid_epiflu`, `ncbi_virus_csv`).

### 2. Quality Control & Subsampling

```yaml
max_n_content: 0.05       # Discards sequences with >5% N bases
min_length: 900           # Excludes short fragments
dedup: true               # Collapses duplicate records (MD5 hash)
subsampling:
  max_sequences: 3000     # Maximum sequences allowed into tree building
  group_by: "country year month" # Columns used to balance the dataset
  sequences_per_group: 10 # Sample density per stratum
```

---

## 3. Evolutionary Selection Parameters

These parameters control the clock-fitting and selection pressure models:

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `clock_rate` | Float | `0.003` | Substitution rate (substitutions/site/year). Pathogen-dependent. |
| `clock_std_dev` | Float | `0.001` | Permitted deviation in molecular clock fitting. |
| `hyphy.methods` | List | `["FUBAR", "MEME"]` | HyPhy methods to run (`FUBAR` for purifying selection, `MEME` for episodic selection, `GARD` for recombination). |
| `hyphy.p_value_threshold` | Float | `0.1` | Statistical significance threshold for selection tests. |

---

## 4. Bridge Builder Scoring (M07)

The `composite_score` balances physical conservation against evolutionary pressure to rate target regions.

$$\text{Composite Score} = (w_c \cdot C) + (w_p \cdot P) - (w_e \cdot E) - (w_r \cdot R)$$

Where:
*   $C$: Conservation score (Shannon Entropy)
*   $P$: Purifying selection density (FUBAR posterior probability)
*   $E$: Episodic diversifying selection penalty (MEME)
*   $R$: Recombination breakpoint proximity penalty (GARD)

### Default Parameter Weights (50/25/8/17)
These weights are pre-calibrated via sensitivity parameter sweeps:
```yaml
scoring_weights:
  conservation: 0.50
  purifying_selection: 0.25
  episodic_penalty: 0.08
  recombination_penalty: 0.17
```

---

## 5. PrimerLab Design & Specificity (M08-M09)

```yaml
primer_design:
  target_type: "RAA"            # "RAA" (isothermal) or "PCR"
  num_candidates: 10            # Number of candidates to generate from Primer3
  max_multiplier: 16            # Max degeneracy pool size (e.g., 16x)
  vienna_ranking_limit: 5       # Max candidates sent to ViennaRNA secondary structure analysis
  optimal_temperature: 39       # Reaction temperature in °C (39°C for RAA, 55-60°C for PCR)
blast_verification:
  mode: "local"                 # "local" (custom DB) or "remote" (NCBI)
  local_db_path: "data/blast_db/pathogen_db"
  evalue_cutoff: 0.01
```

---

## Calibration Guidelines

### Setting `min_length`
*   Locate your target pathogen/segment in the [NCBI Virus Database](https://www.ncbi.nlm.nih.gov/labs/virus/).
*   Identify the full length of the segment (e.g., Influenza A HA segment is ~1700bp).
*   Set `min_length` to **80% - 90%** of the full length to prevent short sequencing fragments from introducing gap artifacts during alignment.

### Setting `clock_rate`
*   Refer to active [Nextstrain builds](https://nextstrain.org/) or literature.
*   Standard reference values:
    *   **Influenza A HA**: `0.0035 - 0.0045` (fast evolving)
    *   **Influenza A MP**: `0.0020 - 0.0030` (slow internal gene)
    *   **SARS-CoV-2**: `0.0008` (proofread, slow)
