---
title: "Shannon Entropy in Diagnostic Design"
description: "Detailed mathematical explanation of Shannon Entropy, how it is calculated across multiple sequence alignments, and its application as the baseline signal for sequence conservation."
---

# Shannon Entropy in Diagnostic Design

In information theory, **Shannon entropy** quantifies the amount of uncertainty or surprise in a set of possible outcomes. In bioinformatics, we use it to measure sequence variability at individual positions within a multiple sequence alignment (MSA). 

This page provides the mathematical formulation of Shannon entropy, a step-by-step calculation example, and explains how the pipeline translates entropy into a normalized physical conservation score ($C$).

---

## Mathematical Formulation

For a given position (column) $i$ in a multiple sequence alignment, the Shannon entropy $H_i$ is calculated as:

$$H_i = -\sum_{b \in \mathcal{A}} p_{ib} \log_2(p_{ib})$$

Where:
- $\mathcal{A}$ is the alphabet of possible states at that position. In our pipeline, we use a 5-letter alphabet: $\mathcal{A} = \{A, T, G, C, -\}$ where `$-$` represents an alignment gap.
- $p_{ib}$ is the observed frequency (proportion) of state $b$ in column $i$, such that $\sum_{b \in \mathcal{A}} p_{ib} = 1.0$.
- By definition, if $p_{ib} = 0$, we define $p_{ib} \log_2(p_{ib}) = 0$.

### Entropy Bounds

- **Minimum Entropy ($H = 0$)**: Occurs when a position is 100% conserved (only one base is present across all sequences). There is zero uncertainty or variation.
- **Maximum Entropy ($H = \log_2(5) \approx 2.322$)**: Occurs when all 5 states in the alphabet are present in equal proportions ($p_{ib} = 0.2$ for each state). There is maximum uncertainty.

---

## Normalised Conservation Score ($C$)

To integrate entropy with other signals in the Bridge Builder (M07), we normalise the entropy value to a scale of $[0, 1]$, where $1.0$ represents maximum conservation and $0.0$ represents maximum randomness.

$$C_i = 1 - \frac{H_i}{\log_2(5)} \approx 1 - \frac{H_i}{2.322}$$

---

## Step-by-Step Calculation Example

Suppose we have an alignment of **10 viral sequences** at a specific position $i$, and we observe the following base counts:
- **Adenine (A)**: 8 sequences
- **Guanine (G)**: 1 sequence
- **Gap (-)**: 1 sequence
- **Thymine (T)** / **Cytosine (C)**: 0 sequences

### Step 1: Calculate Proportions ($p_{ib}$)

$$p_{i\text{A}} = \frac{8}{10} = 0.8$$
$$p_{i\text{G}} = \frac{1}{10} = 0.1$$
$$p_{i\text{-}} = \frac{1}{10} = 0.1$$
$$p_{i\text{T}} = p_{i\text{C}} = 0.0$$

### Step 2: Compute $p \log_2(p)$ for each active state

- **For A**: $0.8 \cdot \log_2(0.8) \approx 0.8 \cdot (-0.3219) \approx -0.2575$
- **For G**: $0.1 \cdot \log_2(0.1) \approx 0.1 \cdot (-3.3219) \approx -0.3322$
- **For -**: $0.1 \cdot \log_2(0.1) \approx 0.1 \cdot (-3.3219) \approx -0.3322$

### Step 3: Sum the terms to find $H_i$

$$H_i = -(-0.2575 - 0.3322 - 0.3322) = 0.9219 \text{ bits}$$

### Step 4: Calculate Normalized Conservation ($C_i$)

$$C_i = 1 - \frac{0.9219}{2.3219} \approx 1 - 0.3970 \approx 0.6030$$

At this position, the conservation score is **0.603**. This is well below the pipeline's strict target threshold (usually $C \geq 0.85$), meaning this position would be avoided when designing primer binding sites.

---

## Integrating Gaps in the Alphabet

Many bioinformatics tools exclude gaps from entropy calculations or treat them as missing data. However, **treating gaps as a 5th character state is critical for diagnostic design**:

- A high gap frequency at a position indicates that a subset of circulating strains contains a deletion relative to the reference.
- Designing a primer over a deletion region will cause immediate hybridization failure in those strains.
- By incorporating gaps into the alphabet, any deletion results in a significant increase in entropy, driving the conservation score $C$ down and disqualifying the region from primer design.

---

## Shannon Entropy vs. Minor Allele Frequency (MAF)

While Minor Allele Frequency (MAF) measures the frequency of the second-most common allele, Shannon entropy provides a more robust signal for multiple-base variation:

| Metric | Formula | 3-Base Wobble Behavior | Sensitivity |
|---|---|---|---|
| **MAF** | $1 - p_{\text{major}}$ | Does not distinguish between one or multiple minor alleles. | Focuses only on the presence of variation. |
| **Shannon Entropy** | $-\sum p \log p$ | Higher entropy if variation is split across multiple alternative alleles (e.g., A/G/C vs. A/G). | Captures the complexity/diversity of variation. |

In the pipeline:
- **Shannon Entropy** is used as the continuous signal to calculate the base conservation score.
- **MAF** is used as a discrete threshold (via `global_params.iupac_threshold`) to decide when to introduce degenerate IUPAC bases into the consensus sequence for PrimerLab.
