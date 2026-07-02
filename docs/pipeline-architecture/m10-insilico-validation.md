---
title: "M10: In-Silico PCR Coverage Analysis"
description: "Simulates PCR/RAA amplification against historical sequences using SeqKit amplicon and generates interactive HTML dashboards and mismatch breakdown reports."
---

# Module M10: In-Silico PCR Coverage Analysis

Module M10 is the validation checkpoint of the pipeline. It evaluates the diagnostic sensitivity (coverage) of the final candidate assays by simulating amplification against a large database of genomic sequences.

---

## Amplification Emulation (SeqKit)

The pipeline utilizes **SeqKit** (`seqkit amplicon`) to simulate physical PCR/RAA amplification across thousands of sequences:

```bash
seqkit amplicon \
  --forward-primer CTGACTAARGGRATTTTRGGRTTTGTATTCA \
  --reverse-primer GCAGTTAAACTATATARGAARCTGAARAGRG \
  --max-mismatch 3 \
  --threads 8 \
  --sequences data/validation_sequences.fasta \
  --out results/m10/amplicons.fasta
```

---

## Custom Binding Kinetics & Mismatch Scoring

A simple mismatch count is not sufficient for evaluating assay performance. Our custom script (`scripts/m10_insilico_pcr.py`) parses the alignment and applies strict molecular kinetics filters:

### 1. 3' Anchoring Penalty (Recombinase Extension Region)
For RAA and RPA isothermal assays, the recombinase enzyme binds to the primer and initiates strand displacement. Mismatches near the 3' end prevent efficient polymerase extension.
- **Critical Region (`raa_3prime_bp`)**: The last **8 bp** of the 3' end of either primer.
- **Tolerance (`raa_max_3prime_mm`)**: Maximum **1 mismatch** allowed in this region. Exceeding this triggers a hard rejection.

### 2. Internal Mismatch Tolerance
- **RAA/RPA**: Up to 4 mismatches are tolerated in the 5' internal region (primer length minus 8 bp) due to the long primer length (30–35 bp) which accommodates higher mismatch loads.
- **PCR**: Up to 2 mismatches are tolerated in the 5' half; mismatches in the 3' critical end must be zero.

### 3. Probe Cleavage Verification
For fluorogenic probes containing an internal Tetrahydrofuran (THF) modification flanked by a fluorophore and quencher, the endonuclease IV (EndoIV) enzyme must cleave the THF site to separate the dyes.
- **Cleavage Window**: Mismatches within **2 bp** of the THF site are prohibited. Any mismatch here prevents enzyme binding and signal release, leading to false-negatives even if amplification succeeds.

---

## How to Read the Output Reports

Module M10 produces three validation reports under `results/{BUILD}/m10/{mode}/`:

### 1. `coverage_summary.tsv`
This file records the raw coverage percentages broken down by subtype, collection year, and country.

| candidate_id | year | subtype | total_sequences | amplified | coverage_pct | 3prime_failures | internal_failures |
|---|---|---|---|---|---|---|---|
| flu_mp_raa_001 | 2024 | H3N2 | 850 | 832 | 97.88 | 12 | 6 |
| flu_mp_raa_001 | 2024 | H1N1 | 540 | 535 | 99.07 | 2 | 3 |
| flu_mp_raa_001 | 2025 | H3N2 | 120 | 114 | 95.00 | 5 | 1 |

### 2. `coverage_report.md`
A human-readable markdown report that summarizes overall performance, highlights any sudden drops in coverage, and lists the accessions of sequences that failed validation.

```markdown
## Overall Coverage: 98.24% (1482 / 1510 sequences)

### Subtype Breakdown:
- **H1N1**: 99.07% (535 / 540)
- **H3N2**: 97.88% (832 / 850)
- **Other (rare, n<30)**: 95.45% (115 / 120)

### Mismatch Analysis & Escaping Strains:
The following accessions failed amplification due to critical 3' mismatches:
- `EPI_ISOLATE_987123` (H3N2, Indonesia, 2025) - Mismatch: A->G at 3' position 2.
- `EPI_ISOLATE_987541` (H3N2, Singapore, 2025) - Mismatch: T->C at 3' position 1.
```

### 3. `coverage_heatmap.html`
An interactive HTML dashboard containing a Plotly heatmap. 

**How to interpret the Heatmap:**
- The **X-axis** displays the collection years.
- The **Y-axis** displays the viral subtypes.
- The **colour scale** ranges from Red (0% coverage) to Green (100% coverage).
- Hovering over a cell shows the exact sample count (`amplified / total_sequences`) and coverage percentage.
- **Actionable Insight**: A sudden shift from green to red in recent years indicates **antigenic/diagnostic drift**. The primer set is losing specificity to emerging strains and needs to be replaced.
- **Filter Rare Subtypes**: Subtypes with fewer than `insilico_pcr.min_subtype_seqs` (default 30) sequences are collapsed into the `"Other (rare, n<30)"` row to prevent statistically insignificant coverage fluctuations from skewing the heatmap.
