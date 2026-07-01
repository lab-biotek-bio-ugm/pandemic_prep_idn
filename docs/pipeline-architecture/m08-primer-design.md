# Module M08: Thermodynamic Primer Design

Module M08 integrates with the **PrimerLab** engine to design primers and probes within the conserved coordinates identified by Module M07. It performs thermodynamic folding validations to prevent assay failure.

---

## The Primer3 Design Logic

The module feeds the IUPAC-aware consensus sequences of candidate regions into Primer3, applying specific parameters based on the reaction type:

*   **Isothermal (RAA/RPA)**:
    *   Primer Length: 30 - 35 bp (optimal 32 bp)
    *   Optimal Temperature: 37 - 42 °C (optimal 39 °C)
    *   GC Content: 30% - 70%
    *   Probe Length: 45 - 53 bp (requires an internal tetrahydrofuran [THF] site, flanked by a fluorophore and quencher)
*   **PCR**:
    *   Primer Length: 18 - 22 bp
    *   Optimal Temperature: 55 - 60 °C

---

## Degeneracy Multiplier Validation

To ensure that designed primers cover viral variants without losing analytical sensitivity, the pipeline checks the **Degeneracy Multiplier** of each candidate primer.

The multiplier ($M$) represents the total number of unique sequence combinations in the degenerate pool:

$$M = \prod_{i=1}^{L} \text{options}(c_i)$$

Where standard bases have 1 option, two-base degeneracies (R, Y, W, S, K, M) have 2 options, three-base degeneracies (B, D, H, V) have 3 options, and four-base degeneracies (N) have 4 options.

*   **Pipeline Constraint**: The pipeline enforces a maximum degeneracy limit of **$16\text{x}$** (defined via `max_multiplier` in config).
*   **Scientific Rationale**: Limiting degeneracy prevents the *primer dilution effect* (where matching primer concentration drops too low) and reduces *non-specific amplification* or *primer-dimer* competition [@chassalevris2020].

---

## ViennaRNA Secondary Structure Validation

Every candidate primer and probe is evaluated for secondary structure thermodynamics at the target reaction temperature (e.g., 39°C for RAA) using the **ViennaRNA package (RNAfold/RNAduplex)**:

1.  **Hairpins**: Calculates the minimum free energy ($\Delta G$) of self-folding. Hairpins with $\Delta G < -3.0\text{ kcal/mol}$ are flagged as high-risk.
2.  **Self-Dimers (Homodimers)**: Simulates binding between two molecules of the same primer. Excludes candidates with homodimerization $\Delta G < -5.0\text{ kcal/mol}$.
3.  **Cross-Dimers (Heterodimers)**: Simulates binding between the forward primer, reverse primer, and probe. Any combinations with cross-hybridization $\Delta G < -6.0\text{ kcal/mol}$ are filtered out.
