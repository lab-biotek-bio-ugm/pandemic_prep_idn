---
title: "Selection Pressure: FUBAR vs. MEME vs. GARD"
description: "A molecular biology guide comparing the three core HyPhy methods used in the pipeline. Explains their mathematical models, detection capabilities, and specific use cases in diagnostic design."
---

# Selection Pressure: FUBAR vs. MEME vs. GARD

The pipeline leverages the **HyPhy (Hypothesis Testing using Phylogenies)** software suite to map evolutionary selection pressure across the target genome. Rather than relying on simple alignment frequencies, HyPhy uses codon-substitution models and phylogenetic topologies to determine if individual positions are actively changing or evolutionarily constrained.

This page compares the three core HyPhy methods implemented in Module M05: **FUBAR**, **MEME**, and **GARD**.

---

## Overview Comparison

| Property | FUBAR | MEME | GARD |
|---|---|---|---|
| **Full Name** | Fast Unconstrained Bayesian AppRoximation | Mixed Effects Model of Evolution | Genetic Algorithm for Recombination Detection |
| **Primary Target** | Pervasive purifying or diversifying selection | Episodic diversifying selection | Recombination breakpoints |
| **Statistical Framework** | Bayesian grid approximation | Maximum Likelihood (Mixed Effects) | Genetic algorithm + c-AIC model comparison |
| **Data Scope** | Across the entire tree | In a subset of branches | Whole-alignment partitions |
| **Assay Design Action** | Positive signal ($P$ score) for purifying | Penalty ($E$ multiplier) to avoid hotspots | Penalty ($R$ multiplier) near breakpoints |

---

## 1. FUBAR: Pervasive Purifying Selection

FUBAR detects codon positions that are under **pervasive selection** — selection pressure that is consistent across all branches of the phylogenetic tree.

### How it works
For each codon site, FUBAR calculates the rate of non-synonymous substitutions ($\beta$) and synonymous substitutions ($\alpha$). It uses a fast Bayesian Markov Chain Monte Carlo (MCMC) grid solver to estimate the posterior probability distribution of $\alpha$ and $\beta$:
- **Purifying (negative) selection**: $\alpha > \beta$ (non-synonymous mutations are being purged).
- **Diversifying (positive) selection**: $\beta > \alpha$ (mutations are being selected for).

### Role in Diagnostic Design
We specifically look for **purifying selection** ($P(\alpha > \beta) \geq 0.90$). 
- High purifying probability indicates that amino acid changes at this position are deleterious to the virus (e.g., they disrupt protein folding, enzyme activity, or virion assembly).
- If the virus mutates this position, the resulting virions are less fit or non-viable, meaning the mutation is highly unlikely to fix in the circulating population.
- **Assay Action**: FUBAR purifying probability acts as the $P$ component of the Bridge composite score, **increasing** the target score of the region.

---

## 2. MEME: Episodic Diversifying Selection

MEME detects codon positions that undergo **episodic selection** — selection pressure that occurs in only a small subset of lineages (branches) while remaining conserved across the rest of the tree.

### How it works
Standard methods (like FEL or FUBAR) average the selection rate across the entire tree. If a position is conserved in 95% of the lineages but undergoes intense positive selection in 5% of the lineages (e.g., during a localized immune-escape event), the average $\beta$ rate remains below $\alpha$, and the positive selection goes undetected.

MEME models the non-synonymous rate ($\beta$) at each site using a two-category distribution:
- $\beta^-$ (sub-neutral rate, $\beta^- \leq \alpha$) with probability $q$.
- $\beta^+$ (positive rate, unrestricted) with probability $1-q$.

It uses a likelihood ratio test (LRT) to compare this mixed-effects model against a null model where no positive selection is allowed.

### Role in Diagnostic Design
Episodic selection sites represent **emerging mutational hotspots**. 
- Even if a site appears 100% conserved in historical sequences, a significant MEME p-value ($p < 0.10$) indicates the site is biologically capable of diversifying under selective sweeps (e.g., vaccine escape or drug resistance).
- Designing primers over these sites creates a high risk of future diagnostic escape as those minor lineages expand.
- **Assay Action**: MEME sites act as the $E$ component in the Bridge Builder, applying a **multiplier penalty** ($<1.0$) to the region's target score or triggering a veto if multiple MEME sites are clustered in the window.

---

## 3. GARD: Recombination Breakpoint Detection

GARD identifies **recombination breakpoints** along the sequence alignment by partitioning the alignment and testing if different segments evolved along different trees.

### How it works
Phylogenetic methods (including FUBAR and MEME) assume that all positions in the alignment share a single evolutionary history (one tree). If recombination has occurred, different regions of the genome will have different histories, causing phylogenetic tools to report false-positive selection signals due to tree-mismatch artifacts.

GARD searches for breakpoint locations using a genetic algorithm. For each candidate breakpoint set, it:
1. Splits the alignment into partitions.
2. Infers a maximum-likelihood tree for each partition.
3. Calculates the combined small sample Akaike Information Criterion (c-AIC).
4. Selects the partition model that minimizes c-AIC.

### Role in Diagnostic Design
- Primers or probes that span a recombination breakpoint are highly vulnerable to assay failure if a new recombinant lineage emerges. Recombination reshuffles entire chunks of genomic material, leading to massive target sequence changes.
- **Assay Action**: GARD output marks breakpoint coordinates. The Bridge Builder applies a **proximity penalty** ($R$) to any candidate regions within $N$ base pairs of a detected breakpoint.

<Warning>
**Pathogen Eligibility**: GARD is computationally expensive and is **disabled for all Influenza builds**. Influenza evolves via segment reassortment, not intra-segment recombination. Enable GARD only for non-segmented RNA viruses like SARS-CoV-2, Dengue, or RSV.
</Warning>

---

## Summary of Assay Design Integration

| Metric | Condition | Biological Interpretation | Bridge Builder Action |
|---|---|---|---|
| **FUBAR posterior** | $\geq 0.90$ | Position is functionally critical; mutations are lethal or reduce fitness. | **Promote** ($+w_p \cdot P$) |
| **MEME p-value** | $< 0.10$ | Position is under active selective pressure in some lineages; prone to mutation. | **Penalize / Veto** |
| **GARD Breakpoint** | Within 50 bp | Proximity to genetic crossover event; unstable across recombinants. | **Penalize** |
