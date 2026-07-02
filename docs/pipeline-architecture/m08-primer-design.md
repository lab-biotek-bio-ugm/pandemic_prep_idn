---
title: "M08: Thermodynamic Primer Design"
description: "Executes PrimerLab integration to design degenerate primers and probes within M07-consensuses. Validates secondary structures via ViennaRNA."
---

# Module M08: Thermodynamic Primer Design

Module M08 integrates with the **PrimerLab** engine to design primers and probes within the target regions identified by Module M07. It applies thermodynamic folding validations to prevent primer-dimer and hairpin structures at reaction temperatures.

---

## The Primer3 Design Logic

The module feeds the IUPAC-aware consensus sequences from M07 into Primer3, applying specific parameters based on the reaction type:

*   **Isothermal (RAA/RPA)**:
    *   Primer Length: 30–35 bp (optimal 32 bp)
    *   Optimal Temperature: 37–42 °C (optimal 39 °C)
    *   GC Content: 30%–70%
    *   Probe Length: 45–53 bp (requires an internal tetrahydrofuran [THF] site, flanked by a fluorophore and quencher)
*   **PCR**:
    *   Primer Length: 18–22 bp
    *   Optimal Temperature: 55–60 °C

---

## Degeneracy & Multiplier Bounds

To ensure designed primers cover circulating variants without losing analytical sensitivity, the pipeline supports **degenerate base design**:
- Variable positions exceeding `degenerate.maf_threshold` (default `0.10` or 10% Minor Allele Frequency) are replaced with the corresponding IUPAC degenerate code (e.g., `R`, `Y`, `W`).
- The **degeneracy multiplier** ($M$) represents the total number of unique sequence combinations in the degenerate pool:

$$M = \prod_{i=1}^{L} \text{options}(c_i)$$

Where standard bases have 1 option, two-base degeneracies (R, Y, W, S, K, M) have 2 options, three-base degeneracies (B, D, H, V) have 3 options, and four-base degeneracies (N) have 4 options.

- **Pipeline Constraint**: The pipeline enforces a maximum degeneracy limit of **$16\text{x}$** (defined via `max_multiplier` in config).
- **Scientific Rationale**: Limiting degeneracy prevents the *primer dilution effect* (where matching primer concentration drops too low) and reduces *non-specific amplification* or *primer-dimer* competition.

---

## ViennaRNA Secondary Structure Validation

Every candidate primer and probe is evaluated for secondary structure thermodynamics at the target reaction temperature (e.g., 39°C for RAA) using the **ViennaRNA package (RNAfold/RNAduplex)**:

1.  **Hairpins**: Calculates the minimum free energy ($\Delta G$) of self-folding. Hairpins with $\Delta G < -3.0\text{ kcal/mol}$ are flagged as high-risk.
2.  **Self-Dimers (Homodimers)**: Simulates binding between two molecules of the same primer. Excludes candidates with homodimerization $\Delta G < -5.0\text{ kcal/mol}$.
3.  **Cross-Dimers (Heterodimers)**: Simulates binding between the forward primer, reverse primer, and probe. Any combinations with cross-hybridization $\Delta G < -6.0\text{ kcal/mol}$ are filtered out.

---

## How to Read the Output Files

Module M08 generates two primary JSON outputs under `results/{BUILD}/m08/{mode}/{region}/`:

### 1. `summary.json`
Contains a ranked list of the top designed candidate primer sets that passed all thermodynamic and degeneracy filters.

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
      "amplicon_length": 187,
      "vienna_score": 0.942,
      "rank": 1
    }
  ]
}
```

### 2. `alternatives.json`
Contains the complete details of all evaluated candidate sets, including their thermodynamic parameters, GC contents, and folding energies. Useful for auditing rejected candidates.

```json
[
  {
    "candidate_id": "flu_mp_raa_001",
    "region_id": "region_0",
    "forward": {
      "seq": "CTGACTAARGGRATTTTRGGRTTTGTATTCA",
      "gc": 41.9,
      "multiplier": 4,
      "hairpin_dg": -1.2,
      "homodimer_dg": -3.1
    },
    "reverse": {
      "seq": "GCAGTTAAACTATATARGAARCTGAARAGRG",
      "gc": 38.7,
      "multiplier": 2,
      "hairpin_dg": -0.8,
      "homodimer_dg": -2.5
    },
    "probe": {
      "seq": "CTTCTTGTTCTTGCTCGAAGACATTTCGTTCTG",
      "gc": 45.4,
      "multiplier": 1,
      "hairpin_dg": -0.5,
      "homodimer_dg": -1.2
    },
    "heterodimer_forward_reverse_dg": -4.2,
    "heterodimer_forward_probe_dg": -3.5,
    "heterodimer_reverse_probe_dg": -4.8,
    "amplicon_length": 187,
    "vienna_score": 0.942,
    "pass_thermo": true
  }
]
```

*Key fields to inspect in `alternatives.json`*:
- `multiplier`: The product of the forward, reverse, and probe degeneracy multipliers. Here, $4 \times 2 \times 1 = 8$, which is well below the `max_multiplier: 16` cap.
- `heterodimer_*_dg`: Ensure none of the heterodimer binding energies exceed the threshold (e.g., must be $> -6.0$ kcal/mol). Values like `-4.2` are stable.
- `vienna_score`: A normalized value representing the thermodynamic suitability. Candidates are sorted in descending order of this score.
