# Module 11: Multiplex Compatibility Screening

The Multiplex Compatibility Screening module evaluates the thermodynamic, kinetic, and spectral compatibility of candidate primers and probes designed for different pathogens or segments, allowing them to be combined into a single reaction tube (Multiplex Assay).

This module supports both **Pairwise (dual-target)** and **Triplex (three-way)** assay compatibility checks.

---

All thresholds, weights and the score formula live in one place — `scripts/tools/iupac_utils.py` — and are imported by both the pairwise and the triplex checker. They previously existed as separate literals in each script and had already drifted apart (the $\Delta T_m$ limit was 3.5°C in one and 4.0°C in the other), so an identical oligo set could be reported `PASS` in one report and `WARNING` in the other.

## Overview of Screening Criteria

To ensure optimal performance and prevent mutual interference of oligos in a multiplex reaction, the screening engine evaluates four criteria:

### 1. Cross-Dimerization Thermodynamics ($\Delta G$)
When multiple primers and probes are mixed together, they can bind to one another instead of the target template, forming primer-dimers.
*   **Calculation**: The engine expands all degenerate sequences into their canonical variants (up to 16 variants per oligo) and computes the minimum free energy ($\Delta G$ in kcal/mol) for every pairwise combination using **Primer3**'s heterodimerization algorithm.
*   **Thresholds** (evaluated at the 39°C RAA reaction temperature):
    *   $\Delta G \ge -9.0$ kcal/mol: **PASS**
    *   $-12.0 < \Delta G < -9.0$ kcal/mol: **WARNING** (Moderate risk, marginal compatibility).
    *   $\Delta G \le -12.0$ kcal/mol: **VETO** (High risk of stable dimer formation; candidate combination is discarded).

!!! note "Why these look loose next to PCR practice"
    There is **no validated RPA/RAA-specific $\Delta G$ threshold** in the literature. The TwistAmp Assay Design Manual §2.2 states that *"there are no fixed rules to predict how well a given amplification primer will work"* and prescribes empirical screening of 8–10 primers per direction instead. The values above therefore derive from general IDT OligoAnalyzer guidance (dimer $\Delta G$ should stay more positive than −9 kcal/mol).

    Crucially, these are computed **at 39°C, not at a qPCR temperature**. Since $\Delta G = \Delta H - T\Delta S$, the same duplex is markedly more stable at the RAA temperature. Measured over 200 random 32-mer pairs with no designed complementarity:

    | | median | 5th pct | below −9 | below −3 |
    |---|---|---|---|---|
    | 39°C (RAA) | −4.12 | −6.78 | 0.5% | 88.0% |
    | 60°C (qPCR) | −1.33 | −3.07 | 0.0% | 6.5% |

    A "−3 kcal/mol" cutoff carried over from 60°C practice would reject 88% of entirely clean pairs. Conversely −9 at 39°C corresponds to roughly −4.7 at 60°C — already stricter than the guidance it came from.

### 2. 3'-End Anchoring
Only a 3' terminus that anneals can be extended by the polymerase, and those extension products are what amplify as **primer noise** — which the TwistAmp manual §2.6 identifies as the factor limiting a primer pair's sensitivity. A global $\Delta G$ can look acceptable while it is the 3' ends that are interacting, so this is screened separately.

*   **Criterion**: the length of the longest perfectly complementary run anchored at either oligo's 3' terminus. Degenerate oligos are expanded and the worst variant decides — a synthesised pool contains all of them.
*   **Threshold**: more than **5 bp** → **WARNING**.

!!! note "Why a length, not a $\Delta G$"
    A $\Delta G$ over a short terminal window was implemented first and rejected on measurement. A 5 nt window **saturates** near −6.3 kcal/mol however much complementarity exists, and widening it to 8 nt still overlaps the clean background badly. Run length separates cleanly against the same random-pair baseline:

    | longest 3' anchor | share of clean pairs |
    |---|---|
    | ≥4 bp | 16.8% |
    | ≥5 bp | 4.2% |
    | ≥6 bp | 0.5% |

    Deliberate 5, 6 and 8 bp anchors are detected in 100% of cases. Because a multiplex evaluates many oligo pairs at once (9 pairwise, 27 triplex), a per-pair false-positive rate compounds; tolerating up to 5 bp keeps the whole-combination false-warning rate near 4% pairwise / 13% triplex. The $\Delta G$ figure is still reported as context for a reader comparing two flagged pairs, but nothing is gated on it.

### 3. Kinetic Temperature Alignment ($\Delta T_m$)
*   **Calculation**: The engine calculates the average $T_m$ across all degenerate variants for each primer. It then computes the maximum difference ($\Delta T_m$) between any two primers.
*   **Thresholds**:
    *   $\Delta T_m \le 4.0^\circ\text{C}$: **PASS**
    *   $\Delta T_m > 4.0^\circ\text{C}$: **WARNING** (Amplification rates might be unbalanced).

This is a **soft** criterion for the primary RAA mode. At 39°C strand invasion is recombinase-driven rather than annealing-driven, so inter-primer $T_m$ matching matters far less than it does in PCR; the check is retained mainly for the 55°C qPCR cross-check.

### 4. Spectral Separation (Dye Assignment)
Probes for different targets must be labeled with fluorophores that have non-overlapping emission spectra to prevent crosstalk. Assignments are read from `config.yaml` (`multiplex.dyes`); the defaults are:
*   **Target A (Influenza A)**: **FAM** (green channel, emission ~520 nm).
*   **Target B (Influenza B)**: **HEX** (yellow channel, emission ~556 nm).
*   **Target C (Measles)**: **Cy5** (far-red channel, emission ~670 nm).

Two targets mapping to the same detection channel are optically inseparable, so a collision is a hard **VETO** and zeroes the score regardless of thermodynamic quality. A dye that is missing or absent from the channel table reports `UNKNOWN` and does **not** pass.

---

## Compatibility Score

$$\text{score} = 100 + \Delta G_{\min} \times 4.0 - \Delta T_m \times 3.0 \qquad (0 \text{ if the dye channels collide})$$

A heuristic ranking aid on a 0–100 scale — **not a percentage and not a probability**. It has no experimental calibration; its job is to decide which combinations are worth synthesising, not to predict assay performance. Reports render it as a bare number for that reason.

$\Delta G$ is weighted above $\Delta T_m$ because that is where the reaction chemistry points: at 39°C strand invasion is enzyme-driven, so $T_m$ matching matters relatively little, while primer-dimer artefacts compete directly for enzyme, primer and energy. The pairwise checker previously used the opposite weighting ($\Delta G \times 3.5$, $\Delta T_m \times 4.0$); both now share one implementation, so pairwise and triplex scores are directly comparable.

!!! tip "Before ordering oligos"
    Treat the ranking as a shortlist, not a verdict. Following the vendor's own protocol, order **several** of the top-ranked combinations rather than rank 1 alone — no in-silico calculation predicts primer noise reliably, and the decisive comparison is sensitivity on the bench.

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
