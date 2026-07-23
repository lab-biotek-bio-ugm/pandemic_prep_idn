---
title: "M07: Diagnostic Bridge Builder & SNP Finder"
description: "Combines Track A sequence conservation with Track B selection pressure scores to identify and rank stable target regions for diagnostic assay design."
---

# Module M07: Diagnostic Bridge Builder & SNP Finder

M07 is the central integration checkpoint of the pipeline. It merges the physical conservation data (Track A) and evolutionary selection pressure scores (Track B) to identify, grade, and export candidate genomic regions suitable for primer and probe design.

---

## The Integration Scoring Algorithm

The core score is **conservation alone**. It is the only quantity that determines whether a primer matches the strains circulating right now:

$$C_i = 1 - \frac{H_i}{\log_2(5)}$$

Everything else — evolutionary and biomechanical — modulates it multiplicatively:

$$\text{Final Window Score} = C_{\text{avg}} \cdot M_{\text{dur}} \cdot M_{\text{meme}} \cdot M_{\text{reco}} \cdot M_{\text{gc}} \cdot M_{\text{block}} \cdot M_{\text{poly}}$$

!!! warning "Changed: purifying selection was an additive term and is now a multiplier"
    Earlier revisions used $S_i = w_c C_i + w_p P_i$ with $w_c = 0.66$, $w_p = 0.34$, where $P_i$ was the FUBAR posterior $\text{Prob}[\alpha > \beta]$. Three problems, all measured:

    **1. Conservation became substitutable.** An additive form lets evolutionary evidence compensate for mismatches. Over the full range of the $P$ term it could overturn a **25.8 percentage-point** gap in conservation — ranking a site at ~14.5% minor-allele frequency above a perfectly conserved one. A primer cannot trade a mismatch today for evidence of constraint tomorrow; it either hybridises or it does not.

    | | $C$ | MAF today | $P$ | old core |
    |---|---|---|---|---|
    | A | 0.80 | 9.86% | 1.00 | **0.868** ← won |
    | B | 0.95 | 1.56% | 0.50 | 0.797 |

    **2. Having HyPhy data could hurt you.** The non-coding branch was `core = C × 0.95` while the coding branch carried the additive term, putting the two on incomparable scales. A perfectly conserved *coding* region with an unmeasured selection signal scored **0.83** against a non-coding region's **0.95**.

    **3. The weights had no source.** No published scheme combines Shannon entropy with dN/dS posteriors in a weighted sum, so `0.66 / 0.34` could not be justified by citation — the quantity does not exist in the literature. The fix was not to find a citation but to restructure so no arbitrary trade-off needs one. The model now has a **single** free parameter, and it has a direct reading (see $M_{\text{dur}}$).

    The **tier** system is unchanged and still uses $\text{Prob}[\alpha > \beta] \ge 0.9$ as an independent gate — the conventional FUBAR criterion, and the classification to cite.

### Multipliers

| Multiplier | Trigger Condition | Default Penalty | Purpose |
|---|---|---|---|
| **$M_{\text{dur}}$** | Nucleotide substitution rate above the gene's median. | Linear from `1.0` at the median down to `floor` (`0.95`) at `rate_cap` (3× median). `0.95` when unmeasured. | Discounts regions whose bases turn over quickly, so the assay stays valid longer. |
| **$M_{\text{meme}}$** | Presence of MEME episodic positive selection sites in the window. | `0.9` for 1 site, `0.7` for 3 sites. **Vetoes** if $>3$ sites. | Avoids emerging mutational hotspots. |
| **$M_{\text{reco}}$** | Proximity to a GARD recombination breakpoint. | `0.5` within `gard_buffer_bp` (50 bp). | Prevents designing across recombinant boundaries. |
| **$M_{\text{gc}}$** | GC content falls outside the ideal range (35%–60%). | `0.7` for marginal ranges (25%–34% or 61%–70%). **Vetoes** if `&lt;` 25% or `&gt;` 70%. | Ensures primer design viability. |
| **$M_{\text{poly}}$** | Presence of homopolymer runs exceeding `max_repeats` (e.g., 5 identical consecutive bases). | `0.6` multiplier. | Prevents primer synthesis and sequencing errors. |
| **$M_{\text{block}}$** | Length of the contiguous conserved block (entropy `&lt;` 0.10) is too short. | `0.8` if block length is between 20 bp and 35 bp. | Ensures sufficient space for primer binding. |


### Why $M_{\text{dur}}$ uses substitution rate, not dN/dS

The purpose of bringing selection analysis into primer design is to avoid targets that look conserved today but drift away tomorrow. Conservation alone cannot distinguish two situations that are identical in an alignment:

- **conserved by functional constraint** — the site cannot change without breaking the protein;
- **conserved by chance** — it has not changed yet, but nothing prevents it.

The obvious metric for that is dN/dS, and that is what the pipeline used. **It measures the wrong thing.** FUBAR's $\text{Prob}[\alpha > \beta]$ says the *protein* is constrained. A primer does not read protein. Under strong purifying selection the nonsynonymous channel is suppressed while the synonymous one — concentrated at the wobble position — is not, so a codon can score highly on $\text{Prob}[\alpha > \beta]$ and still turn its nucleotides over rapidly.

Two codons, scored **identically** by the old model:

| Codon | $\alpha$ (synonymous rate) | $\beta$ (nonsynonymous rate) | $\text{Prob}[\alpha>\beta]$ | nucleotide rate |
|---|---|---|---|---|
| X | 0.1 | 0.05 | high | 0.0625 → **good target** |
| Y | 5.0 | 1.0 | high | 2.0000 → **primer dies fast** |

Y turns over **32×** faster than X. It is exactly the case the evolutionary analysis was added to catch, and the old metric ranked it as favourably as X.

FUBAR reports $\alpha$ and $\beta$ per codon in the same JSON; they were parsed away. The module now computes an expected nucleotide substitution rate,

$$r = 0.25\,\alpha + 0.75\,\beta$$

weighting each channel by roughly how many of a codon's three nucleotide sites belong to it, and compares a region's mean $r$ against **its own gene's median**. FUBAR rates are relative to each alignment's scaling, so a self-calibrating ratio is the only figure that transfers between genes and builds.

### Choosing the floor

`floor` is the model's only free parameter, and it is chosen by stating the trade-off it permits rather than by assertion. With floor $f$, the multiplier's full range can overturn a conservation gap of at most $1 - f$:

| floor | max conservation gap overturned | equivalent MAF |
|---|---|---|
| 0.90 | 10.0 pp | 3.78% |
| **0.95** (default) | **5.0 pp** | **1.56%** |
| 0.97 | 3.0 pp | 0.84% |
| 0.99 | 1.0 pp | 0.23% |
| *(old additive form)* | *25.8 pp* | *14.5%* |

Setting `floor: 1.0` disables the durability term entirely, reducing the model to conservation-only — useful as a baseline for comparison.

The default 0.95 states: *lack of evidence for nucleotide stability costs at most 5%, which can overturn a conservation gap of about 1.6% MAF.* That is a defensible trade — a site at 1.5% MAF under low substitution pressure may well outlast one at 0% that is conserved only by chance in a small sample.

!!! note "Not yet empirically validated"
    This is a mechanistically motivated model, not a proven one. Nobody has demonstrated that rate-informed target selection actually extends assay lifetime. The pipeline already has what is needed to test it: M10 computes coverage **per year**, so a target designed only on data up to year $X$ can be scored against sequences from $X{+}1$ onwards. The hypothesis — regions with low $r$ retain coverage better than equally conserved regions with high $r$ — is directly testable, and a negative result would be just as informative.

### Reported fields

`avg_subst_rate`, `rel_subst_rate` and `durability_multiplier` are written to the region report and CSV so the discount applied to any candidate can be traced.

---

## Grading Tiers

Candidate regions are classified into four tiers based on their average conservation and purifying selection metrics:

Two **different** classifications are emitted, and they must not be confused:

### `tier_grade` — the A/B/C/D scheme (read this one)

Keyed on the **raw** conservation ($C$) and purifying selection ($P$) metrics, and
reported alongside the numeric `tier` (1/2/3).

| `tier_grade` | `tier` | Tier Name | Conservation ($C$) | Purifying Selection ($P$) | Recommendation |
|---|---|---|---|---|---|
| **A** | 1 | **Strict (Platinum)** | $\geq 0.95$ | $\geq 0.90$ | Highly recommended. Extremely stable targets. |
| **B** | 2 | **Moderate (Gold)** | $\geq 0.85$ | $\geq 0.70$ | Safe for routine surveillance assays. |
| **C** | 3 | **Survival (Silver)** | $\geq 0.80$ | Any | Targetable, but may require degenerate bases or validation. |
| **D** | 3 | **Reject** | $< 0.80$ | — | High mutational risk. Do not target. |

When "Grade A or B are the primary targets for PrimerLab", this is the field meant.

### `grade` — fine-grained ranking on the composite score

A separate S / A+ / A / A− / B+ / B / C / D ladder keyed on the **multiplied**
`final_score` (after MEME, recombination, GC, block and homopolymer multipliers).
It is a ranking aid only. Its letters do **not** correspond to the table above —
`grade: "A"` merely means `final_score >= 0.85`, which says nothing about $C$ or $P$
individually.

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