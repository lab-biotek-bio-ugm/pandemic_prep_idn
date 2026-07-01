---
title: "m07-bridge-builder"
---

# Module M07: Bridge Builder Engine

The **Bridge Builder Engine** is the central integration module. It combines physical sequence conservation with evolutionary selection scores to identify and rank candidate diagnostic target regions.

---

## The Composite Scoring Formula

The module calculates a `composite_score` for every position in the segment:

$\text{Composite Score} = (w_c \cdot C) + (w_p \cdot P) - (w_e \cdot E) - (w_r \cdot R)$

Where the scoring factors are:

### 1. Conservation Score ($C$)

Physical nucleotide identity across the alignment, calculated via Shannon Entropy ($H$): $C = 1 - \frac{H}{\log_2(5)}$ Where $H = - \sum (p_i \log_2 p_i)$ over the 5 possible states (A, T, G, C, -). A score of $1.0$ represents 100% conservation.

### 2. Purifying Selection ($P$)

Determined by FUBAR. It is the posterior probability that a site is under negative selection. Higher values ($>0.90$) indicate strong evolutionary constraint.

### 3. Episodic Penalty ($E$)

Determined by MEME. If a site is identified as undergoing episodic diversifying selection ($p\text{-value} < 0.1$), a penalty is applied to avoid mutational hotspots.

### 4. Recombination Breakpoint Penalty ($R$)

Determined by GARD. Positions near recombination breakpoints are penalized to prevent designing primers that span crossing-over boundaries.

---

## Grading Tiers

Candidate regions are ranked and assigned a Grade based on their average `composite_score` across the target window:

| Grade | Score Range | Performance | Recommendation |
| :-: | :-- | :-- | :-- |
| **A** | **Platinum** | Outstanding stability. Highly recommended for first-line diagnostics. |  |
| **B** | **Gold** | Strong target with minor low-frequency SNPs. Safe for routine use. |  |
| **C** | **Silver** | Targetable, but requires degeneracy or careful validation. |  |
| **D** | **Reject** | High mutational risk or low conservation. Avoid. |  |

---

## Output Manifests

- `bridge_report.json`: Master JSON containing scores and metadata flags for every candidate region.
- `conserved_regions.bed`: Browser track file displaying candidate coordinates.
- `entropy_profile.tsv`: Per-nucleotide TSV recording the base frequencies, Shannon entropy, and IUPAC consensus bases.

```json
// Example bridge_report.json entry
{
  "region_id": "region_0",
  "start": 187,
  "end": 218,
  "length": 31,
  "avg_conservation": 0.985,
  "composite_score": 0.932,
  "grade": "A",
  "coding": true,
  "near_recombination": false
}
```