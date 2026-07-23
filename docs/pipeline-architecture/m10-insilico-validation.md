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
- **Critical Region (`raa_3prime_bp`)**: The last **4 bp** of the 3' end of either primer.
- **Tolerance (`raa_max_3prime_mm`)**: Maximum **1 mismatch** allowed in this region. Exceeding this triggers a hard rejection.

!!! note "Why 4 bp, and what this model still misses"
    [Higgins et al. 2022, *J Mol Diagn* 24:1207–1216](https://doi.org/10.1016/j.jmoldx.2022.08.005) is the only study that characterises primer–template mismatches specifically in RPA. It defines the primer **anchor region** as the *pre-ante-penultimate to 3'-terminal position* — four nucleotides — and reports that lone mismatches **outside** the terminal position do not significantly affect time to positivity. Even a lone terminal mismatch left the probability of reaction success at 0.877.

    The previous value of **8** had no traceable source and charged positions 5–8 against the strict 3' quota instead of the internal one, so coverage was **understated** for sequences carrying variation there.

    **What this model still does not capture**: the same study found the identity of the mismatched pair to matter as much as its position. A terminal cytosine–cytosine mismatch was the most damaging, followed by guanine–adenine, while adenine–cytosine and thymine–guanine were *"highly tolerated, rarely resulting in amplification failure"*. The pipeline's position-only criterion treats all four alike. Treat borderline candidates accordingly — a rejection driven by a single A–C terminal mismatch is likely over-conservative.

### 2. Internal Mismatch Tolerance
- **RAA/RPA**: Up to 4 mismatches are tolerated in the internal region (primer length minus `raa_3prime_bp`), reflecting the long primer length (30–35 bp) which accommodates higher mismatch loads. This is conservative relative to the TwistAmp manual §2.7, which reports RPA performance unaffected by 5–9 mismatches across primers and probe.
- **PCR**: Up to 2 mismatches are tolerated in the 5' half; mismatches in the 3' critical end must be zero.

### 3. Probe Cleavage — What Is and Is Not Checked

Fluorogenic RAA probes carry an internal abasic residue (tetrahydrofuran, THF) between a dT-fluorophore and a dT-quencher. The signal only appears when the repair enzyme in the kit — Exonuclease III for `exo` probes, Endonuclease IV for `nfo`/`fpg` probes — cleaves that residue, separating the two dyes. Both enzyme systems share the same probe geometry; the choice affects how the signal is read (real-time fluorescence vs lateral flow), not where the abasic residue goes.

!!! warning "Corrected: there is no 2 bp cleavage window"
    Earlier revisions of this page stated that *"mismatches within 2 bp of the THF site are prohibited... leading to false-negatives even if amplification succeeds"*. **That rule has no support and has been removed.** No literature or vendor documentation was found for it, and the primary source contradicts it:

    - TwistAmp Assay Design Manual §2.7: *"RPA is generally quite robust to mismatches. It has been shown that the performance of RPA assays is not affected when primers and probe carry 5–9 base pair mismatches"*.
    - §3.1.1 sanctions a **deliberate** mismatch inside that very window: *"It is also possible to tolerate a mismatch in the probe such that if two conveniently separated thymines cannot be located one can simply ignore the mismatch of one of the thymines"*. Those thymines are the dT-label positions, which §3.1.2 places **0, 1 or 2 bases** from the THF.

    A supporting-looking result in [NAR 24(22):4572](https://academic.oup.com/nar/article/24/22/4572/2386348) — that Exonuclease III does not cleave duplexes containing an A:A or C:A mismatch — was checked and does not apply: those mismatches were placed **at** the abasic position, as a substitute for it. That is evidence of specificity *for* abasic sites, not of a positional window around one.

    Implementing the rule would also have been self-defeating. Where a probe cannot find two suitably spaced thymines, PrimerLab now deliberately places a dT label on a non-T base — the vendor-sanctioned fallback above. A window check would have flagged the pipeline's own design as a defect. Those positions are exported as `probe_label_mismatches` so they can be discounted rather than counted as target variation.

### What M10 checks instead

A probe is judged on its **total mismatch load** against each surveillance sequence — `raa_probe_max_mm`, default **5**, the conservative end of the 5–9 range the manual describes as tolerated. Probe geometry itself (≥30 bases 5' of the abasic residue, ≥15 bases 3', dT-flanked labels) is enforced upstream at design time in PrimerLab and is not re-derived here.

!!! warning "Corrected: a probe was being judged as if it were a primer"
    Until this revision, M10 ran the probe through the **same** `assess_raa_binding()` routine as the primers, with `raa_3prime_bp=8` and `raa_max_3prime_mm=1`. That applies primer extension kinetics to an oligo that is chemically prevented from being extended: an RAA probe carries a 3' blocking group (C3-spacer, phosphate or biotin-TEG) precisely so the polymerase cannot extend it, so a mismatch at its 3' terminus has none of the consequence that makes 3' mismatches decisive for a primer.

    The effect was a false rejection. Measured on a 46 nt probe:

    | Mismatch pattern | Old (primer kinetics) | New (total load) |
    |---|---|---|
    | 2 mismatches in the last 8 bases | **rejected** | detected (2 mm) |
    | 2 mismatches mid-probe | detected | detected (2 mm) |
    | 5 mismatches scattered | detected | detected (5 mm) |
    | 7 mismatches scattered | rejected | rejected (7 mm) |

    Only probes with more than one mismatch clustered near the 3' end change verdict — but those were being counted as detection failures, so reported RAA coverage was **understated** for exactly the sequences carrying 3'-proximal probe variation.

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
- **Other (rare, n &lt; 30)**: 95.45% (115 / 120)

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
- **Filter Rare Subtypes**: Subtypes with fewer than `insilico_pcr.min_subtype_seqs` (default 30) sequences are collapsed into the `"Other (rare, n &lt; 30)"` row to prevent statistically insignificant coverage fluctuations from skewing the heatmap.
