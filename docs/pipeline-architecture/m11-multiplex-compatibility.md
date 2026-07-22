# Module 11: Multiplex Compatibility Screening

The Multiplex Compatibility Screening module evaluates the thermodynamic, kinetic, and spectral compatibility of candidate primers and probes designed for different pathogens or segments, allowing them to be combined into a single reaction tube (Multiplex Assay).

This module supports both **Pairwise (dual-target)** and **Triplex (three-way)** assay compatibility checks.

---

## Overview of Screening Criteria

To ensure optimal performance and prevent mutual interference of oligos in a multiplex reaction, the screening engine evaluates three key criteria:

### 1. Cross-Dimerization Thermodynamics ($\Delta G$)
When multiple primers and probes are mixed together, they can bind to one another instead of the target template, forming primer-dimers.
*   **Calculation**: The engine expands all degenerate sequences into their canonical variants (up to 16 variants per oligo) and computes the minimum free energy ($\Delta G$ in kcal/mol) for every pairwise combination using **Primer3**'s heterodimerization algorithm.
*   **Thresholds**:
    *   $\Delta G \ge -9.0$ kcal/mol: **PASS** (Safe, minimal risk of primer-dimer formation at working temperatures).
    *   $-12.0 < \Delta G < -9.0$ kcal/mol: **WARNING** (Moderate risk, marginal compatibility).
    *   $\Delta G \le -12.0$ kcal/mol: **VETO** (High risk of stable dimer formation; candidate combination is discarded).

### 2. Kinetic Temperature Alignment ($\Delta T_m$)
For uniform amplification rates across all multiplex targets, the melting temperatures ($T_m$) of all forward and reverse primers in the reaction must align.
*   **Calculation**: The engine calculates the average $T_m$ across all degenerate variants for each primer. It then computes the maximum difference ($\Delta T_m$) between any two primers.
*   **Thresholds**:
    *   $\Delta T_m \le 4.0^\circ\text{C}$: **PASS** (Optimal for uniform amplification).
    *   $\Delta T_m > 4.0^\circ\text{C}$: **WARNING** (Amplification rates might be unbalanced).

### 3. Spectral Separation (Dye Assignment)
Probes for different targets must be labeled with fluorophores that have non-overlapping emission spectra to prevent crosstalk.
*   **Target A (Influenza A)**: Labeled with **FAM** (emission ~520 nm).
*   **Target B (Influenza B)**: Labeled with **HEX** (emission ~556 nm).
*   **Target C (Measles)**: Labeled with **Cy5** (emission ~670 nm).

---

## Orchestrating a Multiplex Run

The screening is managed by the master batch orchestrator:

### 1. Pairwise Screening (2 Targets)
```bash
bash scripts/run_multiplex_batch.sh \
  config/config_influenza_a_mp.yaml \
  config/config_influenza_b_mp.yaml \
  --use-conda -c 64
```
Runs the Snakemake pipelines for both targets and calls `scripts/m11_multiplex_checker.py` to evaluate the top $10 \times 10 = 100$ combinations.

### 2. Triplex Screening (3 Targets)
```bash
python scripts/run_multiplex_batch.py \
  --configs config/config_influenza_a_mp.yaml config/config_influenza_b_mp.yaml config/config_measles_all_partial.yaml \
  --cores 64 \
  --skip-snakemake
```
Runs the three-way compatibility check using `scripts/m11_triplex_checker.py` across all $10 \times 10 \times 10 = 1,000$ combinations.

---

## Key Outputs

Results are written to target-specific directories:
*   **Pairwise**: `results/multiplex-<targetA>_x_<targetB>/`
    *   `multiplex_report.md`: Markdown summary report.
    *   `multiplex_candidates.json`: Raw thermodynamic details.
*   **Triplex**: `results/multiplex-triplex-all-3/`
    *   `triplex_report.md`: Markdown summary report with synthesis specifications.
    *   `triplex_candidates.json`: Sorted list of triplex candidates.
