---
title: "M07: Diagnostic Bridge Builder & SNP Finder"
description: "Combines Track A sequence conservation with Track B selection pressure scores to identify and rank stable target regions for diagnostic assay design."
---

# Module M07: Diagnostic Bridge Builder & SNP Finder

M07 is the central integration checkpoint of the pipeline. It merges the physical conservation data (Track A) and evolutionary selection pressure scores (Track B) to identify, grade, and export candidate genomic regions suitable for primer and probe design.

---

## The Integration Scoring Algorithm

For every position $i$ in the alignment, the Bridge Builder script (`scripts/m07_bridge_builder.py`) calculates a normalized composite score $S_i$ using the following weighted formula:

$$S_i = (w_c \cdot C_i) + (w_p \cdot P_i)$$

Where:
- $C_i$ is the normalized physical conservation score derived from Shannon Entropy: $C_i = 1 - \frac{H_i}{\log_2(5)}$.
- $P_i$ is the FUBAR posterior probability that the site is under purifying selection ($\alpha > \beta$).
- $w_c$ is the conservation weight (default `0.66`).
- $w_p$ is the purifying weight (default `0.34`), such that $w_c + w_p = 1.0$.

### Window Scoring and Multipliers

Once position-level scores are calculated, the script uses a sliding window (default `window_size: 150` bp or evaluates the full segment dynamically depending on `bridge.strategy`) to compute the average base score. This base score is then adjusted by several **multiplicative penalties** (rather than subtractive penalties) to avoid regions with unfavorable molecular properties:

$$\text{Final Window Score} = S_{\text{avg}} \cdot M_{\text{meme}} \cdot M_{\text{reco}} \cdot M_{\text{gc}} \cdot M_{\text{block}} \cdot M_{\text{poly}}$$

| Multiplier | Trigger Condition | Default Penalty | Purpose |
|---|---|---|---|
| **$M_{\text{meme}}$** | Presence of MEME episodic positive selection sites in the window. | `0.9` for 1 site, `0.7` for 3 sites. **Vetoes** if $>3$ sites. | Avoids emerging mutational hotspots. |
| **$M_{\text{reco}}$** | Proximity to a GARD recombination breakpoint. | `0.5` within `gard_buffer_bp` (50 bp). | Prevents designing across recombinant boundaries. |
| **$M_{\text{gc}}$** | GC content falls outside the ideal range (35%–60%). | `0.7` for marginal ranges (25%–34% or 61%–70%). **Vetoes** if $<25\%$ or $>70\%$. | Ensures primer design viability. |
| **$M_{\text{poly}}$** | Presence of homopolymer runs exceeding `max_repeats` (e.g., 5 identical consecutive bases). | `0.6` multiplier. | Prevents primer synthesis and sequencing errors. |
| **$M_{\text{block}}$** | Length of the contiguous conserved block (entropy $< 0.10$) is too short. | `0.8` if block length is between 20 bp and 35 bp. | Ensures sufficient space for primer binding. |

---

## Grading Tiers

Candidate regions are classified into four tiers based on their average conservation and purifying selection metrics:

| Grade | Tier Name | Conservation ($C$) | Purifying Selection ($P$) | Recommendation |
|---|---|---|---|---|
| **A** | **Strict (Platinum)** | $\geq 0.95$ | $\geq 0.90$ | Highly recommended. Extremely stable targets. |
| **B** | **Moderate (Gold)** | $\geq 0.85$ | $\geq 0.70$ | Safe for routine surveillance assays. |
| **C** | **Survival (Silver)** | $\geq 0.80$ | Any | Targetable, but may require degenerate bases or validation. |
| **D** | **Reject** | $< 0.80$ | — | High mutational risk. Do not target. |

---

## How to Read the Output Files

Module M07 generates three primary outputs in `results/{BUILD}/m07/`:

### 1. `bridge_report.json`
This is the master summary of all candidate regions. It is parsed by downstream design scripts and human researchers.

```json
{
  "build_name": "influenza-a-mp-sea",
  "regions": [
    {
      "region_id": "region_0",
      "start": 187,
      "end": 340,
      "length": 153,
      "gene": "M1",
      "avg_conservation": 0.985,
      "avg_purifying_prob": 0.921,
      "meme_sites": 0,
      "near_recombination": false,
      "gc_content": 44.5,
      "composite_score": 0.963,
      "grade": "A",
      "tier": "strict",
      "coding": true
    }
  ]
}
```

*Key fields to inspect:*
- `grade`: Target regions with Grade `A` or `B` are the primary targets for PrimerLab.
- `coding`: Indicates if the region lies within a protein-coding sequence (CDS). Non-coding regions (like UTRs) will have `coding: false` and will not have purifying selection ($P$) scores.

### 2. `primer_targets_manifest.yaml`
This file is the **direct input for Module M08 (PrimerLab)**. It defines the exact coordinates, IUPAC consensus sequences, and exclusion coordinates (variable sites to be avoided by Primer3).

```yaml
candidate_regions:
  - region_id: "region_0"
    start_bp: 187
    end_bp: 340
    consensus_sequence: "ATGGATGTCAAYCCGACTCAGCTACTC..."
    exclusion_zones:
      - 214  # Minor SNP site (MAF > 15%)
      - 289
```

### 3. `conserved_regions.bed`
A standard BED format file containing the coordinates of all Grade A, B, and C regions. This file can be loaded directly into genome browsers (like IGV or UCSC) for visual verification.

```
NC_002016.1	187	340	region_0_Grade_A	963	+
```

---

## Visualizing Output in IGV

To verify the candidates visually:
1. Open the [Integrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/).
2. Load the reference genome: `config/references/flu_a_mp_reference.fasta`.
3. Load the candidate BED file: `results/{BUILD}/m07/conserved_regions.bed`.
4. Load the GFF3 annotation: `config/references/flu_a_mp_annotation.gff`.

You will see candidate regions displayed as tracks aligned to the reference genome and coding gene frames. Choose regions that fall fully within your target gene (e.g., `M1`) and avoid regions that sit near splicing boundaries unless targeting spliced transcripts specifically.