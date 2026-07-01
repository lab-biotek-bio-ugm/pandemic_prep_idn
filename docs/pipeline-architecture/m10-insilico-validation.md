# Module M10: In-Silico Validation

Module M10 performs large-scale in-silico PCR and RAA amplification simulations against historical surveillance databases (e.g., GISAID metadata, local clinical samples) to validate the diagnostic sensitivity and coverage of the candidate assays.

---

## Amplification Emulation (SeqKit)

The pipeline utilizes **SeqKit** (`seqkit amplicon`) to simulate physical PCR/RAA amplification across thousands of genomic sequences:

```bash
seqkit amplicon \
  --forward-primer CTGACTAARGGRATTTTRGGRTTTGTATTCA \
  --reverse-primer GCAGTTAAACTATATARGAARCTGAARAGRG \
  --max-mismatch 3 \
  --sequences data/validation_sequences.fasta \
  --out results/m10/amplicons.fasta
```

---

## Custom Binding Kinetics & Mismatch Scoring

While SeqKit checks for simple mismatches, our custom validation script applies a thermodynamic/kinetic penalty score to evaluate the impact of mismatches on hybridization:

1.  **3' Anchoring Penalty**: Mismatches within the last 5 nucleotides of the 3' end of either primer are heavily penalized. A single mismatch here can decrease binding efficiency by $>50\%$.
2.  **Mismatch Count Tolerance**:
    *   **RAA/RPA**: Up to 3 mismatches are tolerated if they are distributed across the middle or 5' region of the primer.
    *   **PCR**: Up to 2 mismatches are tolerated in the 5' half; mismatches in the 3' end must be zero.
3.  **Probe Validation**: Mismatches in the probe sequence are tolerated up to 2 bases, but no mismatch is allowed near the internal THF modification site to ensure efficient cleavage and signal release.

---

## Coverage Analytics & Reports

The module aggregates the in-silico validation results into distinct reports:

### 1. Temporal Coverage Report
Calculates the percentage of sequences successfully amplified year-by-year (e.g., 2005 to 2026). This identifies whether the assay is vulnerable to diagnostic drift over time.

### 2. Regional Coverage Report
Compares the diagnostic sensitivity across different geographic divisions (e.g., Indonesia only vs. Global isolates) to check for regional variations in performance.

### 3. Interactive HTML Dashboard
Generates a standalone, interactive HTML dashboard containing:
*   A line plot of sensitivity (%) over time.
*   A bar chart of regional coverage.
*   A summary table of mismatched isolates (to identify emerging strains that might escape detection).
