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
    *   Probe Length: 45–53 bp (requires an internal abasic site, flanked by a fluorophore and quencher — see below)
*   **PCR**:
    *   Primer Length: 18–22 bp
    *   Optimal Temperature: 55–60 °C

---

## RAA Probe Geometry (Abasic Site Placement)

An RAA probe is only useful if the kit's repair enzyme can actually cleave it. PrimerLab places the abasic residue according to the TwistAmp Assay Design Manual, and screens probe candidates for a legal placement **before** ranking them thermodynamically:

| Rule | Value | Why |
|---|---|---|
| Bases 5' of the abasic residue | ≥ 30 (working range 30–38) | The cleaved 5' fragment must still function as a primer. |
| Bases 3' of the abasic residue | ≥ 15 | Below this the exonuclease does not cut efficiently. |
| Bases between a dT label and the abasic residue | 0, 1 or 2 | Keeps fluorophore–quencher separation ≤ 5, preserving quenching. |
| Flanking label positions | must be **T** | The fluorophore and quencher reagents exist only as dT couplings. |
| The base replaced by the abasic residue | any | The manual states there is no sequence requirement here. |

The same geometry applies whether the kit uses Exonuclease III (`exo` probes) or Endonuclease IV (`nfo`/`fpg` probes); the probe type changes how the signal is read, not where the abasic residue goes.

!!! warning "Corrected in this revision"
    Two defects in the previous placement logic were found by checking it against the manual's worked example and have been fixed:

    1. The abasic residue could land **28–29** bases from the 5' end — a case the manual lists verbatim as a *"Poor probe design"*. There was no upper bound either.
    2. **The label logic was inverted.** The old code anchored the abasic residue onto the nearest `T` — consuming the one base that could have carried a label — and let both dT labels fall on whatever happened to be adjacent. Every probe could therefore silently carry two unnecessary probe–target mismatches.

### When no ideal site exists

Placement degrades in the order the manual itself prescribes, and every relaxation is reported rather than applied silently:

1. Pair of thymines inside the 30–38 window → compliant, no mismatch.
2. Pair of thymines beyond 38 → used, with a warning.
3. Only one usable thymine → one label mismatch. The manual explicitly sanctions this: *"one can simply ignore the mismatch of **one** of the thymines... There may be a reduction in the efficiency of the probe, the extent to which cannot be predicted."*
4. No usable thymine → two label mismatches, warned.
5. No legal geometry at all → the probe is **rejected**, not annotated at a clamped position.

Cases 3–5 are exported as `probe_label_mismatches` alongside `thf_index` and `probe_compliant`, so downstream validation counts them as design properties rather than as target variation.

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

1.  **Hairpins**: Calculates the minimum free energy ($\Delta G$) of self-folding. Hairpins with $\Delta G \text{ } \lt -3.0\text{ kcal/mol}$ are flagged as high-risk.
2.  **Self-Dimers (Homodimers)**: Simulates binding between two molecules of the same primer. Excludes candidates with homodimerization $\Delta G \text{ } \lt -5.0\text{ kcal/mol}$.
3.  **Cross-Dimers (Heterodimers)**: Simulates binding between the forward primer, reverse primer, and probe. Any combinations with cross-hybridization $\Delta G \text{ } \lt -6.0\text{ kcal/mol}$ are filtered out.

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
