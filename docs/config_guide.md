# Pipeline Configuration Guide

This guide provides a detailed explanation of the parameters used in the YAML configuration profiles (e.g., `config_influenza_a_mp.yaml`). Understanding these parameters is crucial for adapting the pipeline to different pathogens and surveillance contexts.

---

## 1. Target Identity & Data Input

- `build_name`: A unique identifier for the analysis run. This name is used to organize output directories in the `results/` folder.
- `local_input`:
  - `sequences`: Path to the input FASTA file.
  - `metadata`: Path to the associated CSV/TSV metadata.
  - `format`: The metadata schema (e.g., `gisaid_epiflu`, `ncbi_virus_csv`). The pipeline uses this to map columns correctly.

## 2. Quality Control (QC)

- `max_n_content`: Maximum percentage of ambiguous bases ('N') allowed. Sequences exceeding this are discarded to ensure high-quality alignments.
- `dedup`: If set to `true`, identical sequences (based on MD5 hash) are removed to reduce computational load without losing diversity signals.

## 3. Subsampling (Spatio-Temporal Filter)

This module prevents the analysis from being biased by regions or years with excessive over-sampling.

- `group_by`: Metadata columns used to balance the data (e.g., `subtype country year month`).
- `max_sequences`: The total target number of sequences allowed into the core analysis.
- `min_length`: Minimum sequence length requirement. For fragmented segments like Influenza HA (~1700bp), this should be adjusted to keep full-length or near-full-length reads.

## 4. Timetree Analysis (Augur Refine)

Parameters for the molecular clock used to estimate divergence times.

- `clock_rate`: The expected substitution rate (substitutions/site/year).
  - *Influenza A MP*: ~0.003
  - *Influenza A HA*: ~0.0035 (HA evolves faster due to antigenic pressure)
  - *SARS-CoV-2*: ~0.0008
- `clock_std_dev`: The allowed deviation for the molecular clock.

## 5. Evolutionary Selection Engine (HyPhy)

- `methods`: The statistical tests to run.
  - `FUBAR`: Fast Unconstrained Bayesian Approximation (Detects pervasive purifying selection).
  - `MEME`: Mixed Effects Model of Evolution (Detects episodic diversifying selection).
  - `GARD`: Genetic Algorithm for Recombination Detection (Identifies breakpoints to prevent false-positive selection signals).
- `p_value_threshold`: Significance level for mutation detection (default `0.1`).

## 6. Bridge Builder Scoring (M07)

This is the core decision-making module. It calculates a `composite_score` to rank potential diagnostic targets.

### Weighted Scoring (50/25/8/17)
Target stability is prioritized through a balanced weighting system. These weights reflect **Data-Driven Defaults** derived from a **Sensitivity Analysis (Parameter Sweep)** to ensure that identified targets are evolutionarily robust.

1. **Conservation (50%)**: Physical sequence identity across the dataset. Ensures the primer/probe has a high binding probability.
2. **Purifying Selection (25%)**: Evolutionary restriction. Indicates that mutations in this region are likely deleterious to the virus, ensuring long-term target stability.
3. **Episodic Penalty (8%)**: Mutational noise. Penalizes regions that show sudden, temporary bursts of mutation.
4. **Recombination Penalty (17%)**: Genomic exchange. Penalizes regions near recombination breakpoints to prevent primer failure due to genetic shifts.

### Grading Tiers

- **Strict (Tier 1)**: Recommended for diagnostic primers/probes (Min Conservation > 0.98).
- **Moderate (Tier 2)**: Suitable for research assays or screening.
- **Survival (Tier 3)**: Minimum baseline for emergency surveillance.

---

## 7. How to Determine Parameter Values

If you are adapting this pipeline to a new virus or segment, use the following methodology to determine the correct values:

### Determining `min_length`
1.  Search for the official reference genome/segment on [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/).
2.  Check the total sequence length (e.g., Influenza A H1N1 HA is ~1700bp).
3.  Set `min_length` to approximately **80% - 90%** of the full length. This ensures you exclude short fragments that might introduce gaps in the alignment while keeping most of the relevant data.

### Determining `clock_rate`
1.  Consult the [Nextstrain Documentation](https://docs.nextstrain.org/en/latest/guides/index.html) or their active builds (e.g., `nextstrain.org/flu`).
2.  Look for the "Substitution Rate" in the "Calibration" or "Refine" section of their published builds.
3.  If Nextstrain is unavailable, search literature for *"Substitution rate per site per year"* for your target pathogen.
4.  Standard reference values:
    - Highly variable (HA/NA): `0.003 - 0.005`
    - High-fidelity/Slow (SARS-CoV-2): `0.0008`
    - Conservative internal segments (MP/NP): `0.002 - 0.003`

### Calibrating `subsampling`
- **Total sequences**: For global surveillance with stable RAM (32GB+), `3000-5000` sequences is a good balance. For rapid screening, `500-1000` is sufficient.
- **Group counts**: Set `sequences_per_group` low (e.g., `5-10`) if you want maximum diversity across countries, or higher if you want to see local outbreak patterns.

---

> [!TIP]
> **Adapting to New Pathogens**
> When working with highly diverse viruses, consider slightly lowering the `entropy_threshold` and increasing the `recombination` weight to account for higher genetic volatility.
