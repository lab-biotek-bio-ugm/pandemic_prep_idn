---
title: "Evolutionary Diagnostics"
description: "The scientific rationale for combining physical sequence conservation with evolutionary selection pressure analysis to design mutation-resistant diagnostic assays for pandemic surveillance."
---

# Evolutionary Diagnostics: Scientific Rationale

This pipeline is built on a paradigm called **evolutionary diagnostics** — the use of population genomics and molecular evolution data to engineer assay components that are resistant to viral mutation. This page explains the problem, the two-signal framework, and why naive approaches fail.

---

## The Diagnostic Escape Problem

Standard PCR and isothermal assay design assumes that the target sequence is stable. For slowly-evolving pathogens (bacteria, DNA viruses), this assumption holds. For RNA viruses — which have error-prone RNA-dependent RNA polymerases (RdRp) with mutation rates of $10^{-3}$ to $10^{-6}$ substitutions/site/replication cycle — this assumption is routinely violated.

**Diagnostic escape** occurs when a virus accumulates sufficient mutations at the primer or probe binding site that the assay no longer amplifies efficiently. The clinical consequence is a **false-negative result**: an infected patient is classified as disease-free.

Historical examples include:
- A 2009 H1N1 outbreak in which a widely-distributed PCR kit missed the pandemic strain due to a single mismatch at the 3' end of the forward primer.
- SARS-CoV-2 Alpha (B.1.1.7) and Delta variants that reduced sensitivity of several commercial diagnostic kits targeting the N and ORF1a genes.

---

## Why Physical Conservation Is Insufficient

The naive approach to mutation-resistant design is to find **conserved regions** — positions where the nucleotide is identical across all known sequences. This is a necessary condition, but not sufficient for three reasons:

### 1. Conservation Reflects History, Not Future Constraint

A position that is identical in all currently-deposited sequences may be conserved by chance (founder effect, sampling bias) rather than by functional necessity. If that position is not under evolutionary constraint, a single selective sweep can fix a new variant globally within one influenza season.

### 2. Episodic Diversification Is Invisible to Sequence Counts

Some positions undergo **bursts of selection** in specific lineages or time windows — they appear highly conserved when looking at the global alignment but are actually under active diversifying pressure in a subset of circulating strains. Standard conservation metrics (proportion identical, Shannon entropy) cannot distinguish between these scenarios.

### 3. Recombination Creates Local Instability

In some RNA viruses (coronaviruses, picornaviruses), intragenic recombination can dramatically reshuffle local sequence composition even in regions that appear conserved in the absence of recombination. Designing primers that span recombination breakpoints risks false-negative results in recombinant lineages.

---

## The Two-Track Framework

This pipeline integrates two independent evidence tracks to address all three limitations:

```
Track A: Physical Conservation (Shannon Entropy)
          ↓
          Per-position conservation score (C)
          
Track B: Evolutionary Selection Pressure (HyPhy)
          ↓
          Purifying selection score (P) — FUBAR
          Episodic diversification penalty (E) — MEME
          Recombination proximity penalty (R) — GARD
          
          ↓
Bridge Builder (M07): Composite Score
          
C × wc + P × wp - E × we - R × wr
```

### Track A — Shannon Entropy

Shannon entropy $H$ quantifies the randomness (uncertainty) of nucleotide identity at a given alignment position:

$$H_i = -\sum_{b \in \{A,T,G,C,-\}} p_{ib} \log_2(p_{ib})$$

where $p_{ib}$ is the proportion of base $b$ at position $i$. A perfectly conserved position has $H = 0$; a position with equal proportions of all four bases and a gap has $H = \log_2(5) \approx 2.32$.

The **conservation score** normalises entropy to $[0, 1]$:

$$C_i = 1 - \frac{H_i}{\log_2(5)}$$

$C = 1.0$ means the position is 100% identical across all sequences. $C = 0$ means completely random.

**Limitation**: $C$ is a historical frequency metric. It cannot distinguish between positions constrained by function and those that happen to be uniform due to founder effects.

### Track B — Evolutionary Selection Pressure

HyPhy models the **ratio of non-synonymous to synonymous substitution rates** ($dN/dS$, also written $\omega$) at each codon position across the phylogeny:

- **$\omega < 1$** (purifying/negative selection): Non-synonymous mutations are deleterious and are purged from the population by natural selection. The encoded amino acid is functionally critical.
- **$\omega = 1$** (neutral evolution): Mutations have no fitness effect.
- **$\omega > 1$** (positive/diversifying selection): Non-synonymous mutations confer a fitness advantage (e.g., immune escape, receptor-binding adaptation).

**FUBAR** (Fast Unconstrained Bayesian AppRoximation) estimates the Bayesian posterior probability $P(\omega < 1)$ — the probability that a site is under **pervasive purifying selection** across the entire evolutionary history of the pathogen. A high $P(\omega < 1)$ means the site has been evolutionarily constrained across many millions of replications — a much stronger signal than simple sequence counting.

**MEME** (Mixed Effects Maximum Likelihood Evolution) identifies sites under **episodic positive selection** in specific lineages. Unlike FUBAR (which tests for pervasive effects across the entire tree), MEME can detect positions that are conserved globally but experience selection in a subset of branches. The MEME p-value for a site is the probability of observing the detected rate variation by chance.

**GARD** (Genetic Algorithm for Recombination Detection) locates **recombination breakpoints** within the alignment using a genetic algorithm that optimises the c-AIC across candidate breakpoint models. Positions near detected breakpoints receive a penalty because primers spanning recombination boundaries are unreliable across recombinant lineages.

---

## The Composite Score and Grading

The Bridge Builder (M07) computes a weighted composite score $S$ for each candidate region:

$$S = (w_c \cdot C) + (w_p \cdot P) - \text{multipliers}(E, R, \text{GC}, \text{HP})$$

Where:
- $C$ = conservation score from Track A entropy analysis
- $P$ = FUBAR purifying selection posterior probability
- Multipliers are applied (not subtracted linearly) for MEME episodic sites, recombination proximity, GC content deviation, and homopolymer runs

Candidate regions are classified into tiers:

| Grade | Tier | Conservation Threshold | Purifying Selection | Recommendation |
|---|---|---|---|---|
| A | Strict (Platinum) | ≥ 0.95 | ≥ 0.90 posterior | First-line diagnostic target |
| B | Moderate (Gold) | ≥ 0.85 | ≥ 0.70 posterior | Safe for routine use |
| C | Survival (Silver) | ≥ 0.80 | Any | May require degenerate bases; validate carefully |
| D | Reject | < 0.80 | — | Do not use; high mutation risk |

---

## Why Both Tracks Are Necessary

The synergy between the two tracks is best illustrated by contrasting three scenarios:

| Scenario | Track A (Conservation) | Track B (Selection) | Interpretation |
|---|---|---|---|
| **Optimal target** | High ($C \geq 0.95$) | Strong purifying ($P \geq 0.90$) | Evolutionarily constrained — mutations are lethal to the virus |
| **Historically conserved, evolutionarily neutral** | High ($C \geq 0.95$) | Low/no signal ($P \approx 0.5$) | Conserved by chance; may diversify suddenly under immune pressure |
| **Actively diversifying** | Low ($C < 0.80$) | High positive ($P(\beta>\alpha) > 0.9$) | Immune escape hotspot — avoid entirely |
| **Episodic diversification** | Moderate-high ($C \approx 0.88$) | MEME significant ($p < 0.05$) | Globally conserved but actively selected in some lineages; apply penalty |

A design based solely on Track A would miss the "historically conserved, evolutionarily neutral" risk. A design based solely on Track B would miss regions that are evolutionarily constrained at the molecular level but happen to have low current sequence diversity.

---

## Further Reading

For the statistical models underlying this framework:

- Murrell et al. (2013). FUBAR: A Fast, Unconstrained Bayesian AppRoximation for Inferring Selection. *Molecular Biology and Evolution*.
- Murrell et al. (2012). Detecting Individual Sites Subject to Episodic Diversifying Selection. *PLOS Genetics*.
- Pond et al. (2006). GARD: a genetic algorithm for recombination detection. *Bioinformatics*.
- Shannon, C.E. (1948). A Mathematical Theory of Communication. *Bell System Technical Journal*.
