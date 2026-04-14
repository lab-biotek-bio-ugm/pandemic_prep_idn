# Technical Manual: Influenza Evolutionary Diagnostic Pipeline

This document provides a step-by-step technical explanation of the pipeline's instrumentation, logic, and output validation. It is intended for laboratory bioinformaticians and researchers to ensure scientific reproducibility and ease of maintenance.

---

## Architecture Overview

The pipeline is structured into 7 core modules (**M01 – M07**), moving from raw sequence acquisition to a finalized list of ranked diagnostic targets. It leverages the **Nextstrain (Augur)** ecosystem for phylogenomics and **HyPhy** for selection pressure analysis.

---

## M01: Data Acquisition & Enrichment

**Objective**: Standardize raw data from heterogeneous sources (GISAID, NCBI, local labs) and perform early categorization.

### Technical Process
1.  **Metadata Mapping**: Maps source-specific headers to a unified schema (Nextstrain/Augur internal format).
2.  **Nextclade Categorization**: Every sequence is run through Nextclade to assign official Clades and Lineages. This step ensures our local data is synchronized with the global database.

### Key Scripts
- `rules/m01_acquisition.smk`
- `scripts/prepare_metadata_nextclade.py`

> [!TIP]
> **CLI Debugging**:
> To run Nextclade manually for a specific sample:
> ```bash
> nextclade run --input-fasta sample.fasta --dataset-name flu_h3n2_ha --output-all results/debug_m01/
> ```

---

## M02: Quality Control & Filtering

**Objective**: Eliminate low-quality data and balance the dataset to prevent spatio-temporal bias.

### Technical Process
1.  **Hard QC**: Removes sequences with >5% N-content (`max_n_content`) and duplicate records.
2.  **Adaptive Subsampling**: A critical step that uses **Spatio-Temporal Filtering**. It balances the number of sequences across countries, years, and months to ensure the evolution of the virus is not skewed by regions with excessive reporting.

### Key Instrumentation
- `augur filter`: The primary engine for subsampling logic.

> [!IMPORTANT]
> **Rationale**: If we have 10,000 sequences from the USA and 10 from Indonesia, the tree will look like a "USA Tree". Subsampling reduces the USA sequences while keeping all Indonesian ones to provide a balanced evolutionary view.

---

## M03: Multiple Sequence Alignment (MSA)

**Objective**: Align all filtered sequences against a standard reference to establish a unified coordinate system.

### Technical Process
- **Reference-Based Alignment**: Unlike *de novo* alignment, we use a fixed **RefSeq** (e.g., NC_002016.1 for Flu A MP). 
- **Tool**: MAFFT (invoked via `augur align`).
- **Gap Filling**: Bases corresponding to gaps in the reference are filled with `-` or `N` to maintain fixed alignment width.

### Key Instrumentation
- `augur align --fill-gaps`

> [!TIP]
> **CLI Debugging**:
> To check the alignment for a single file:
> ```bash
> augur align --sequences filtered.fasta --reference-sequence ref.fasta --output aligned.fasta --fill-gaps
> ```

---

## M04: Core Phylodynamics & Timetree

**Objective**: Estimate the phylogenetic relationships and time-scale of divergence.

### Technical Process
1.  **Phylogeny Inference**: Uses **IQ-TREE** (Maximum Likelihood) to build the `tree_raw.nwk`. It models the substitution process to find the most likely tree topology.
2.  **Refinement**: Uses `augur refine` to:
    -   Resolve polytomies stochastically.
    -   Apply a **Molecular Clock** (`clock_rate`) to calculate branch lengths in years.
    -   Root the tree (typically using root-to-tip regression or a specified outgroup).

### Key Instrumentation
- `augur tree` (IQ-TREE)
- `augur refine`

> [!NOTE]
> **Clock Rate Calibration**: Getting the `clock_rate` right is essential for M06 (Ancestral Reconstruction). A wrong clock rate will skew the estimation of when mutations occurred.

---

## M05: Evolutionary Selection Analysis (HyPhy)

**Objective**: Identify which genomic regions are under active mutational pressure.

### Technical Process
1.  **Sanitization**: Converts the nucleotide alignment into a codon-aware alignment (stripping stop codons and partial triplets).
2.  **FUBAR (Purifying Selection)**: Identifies sites where amino acid changes are rejected by evolution. These are "Stable Zones" ideal for primers.
3.  **MEME (Episodic Selection)**: Identifies sites under "Mutational Noise" (transient bursts of mutation). These are "Danger Zones" to be avoided.
4.  **GARD (Recombination)**: Detects if parts of the genome have different evolutionary histories, preventing false-positive selection scores.

### Key Instrumentation
- `scripts/run_hyphy_gard_safe.py` (Custom wrapper for standard HyPhy methods).

---

## M06: Ancestral State Reconstruction

**Objective**: Inferred sequences for all internal nodes in the tree.

### Technical Process
- **Reconstruction**: Uses Joint or Marginal Likelihood (Augur) to "back-calculate" what the ancestors of your current samples looked like. 
- **Internal vs Global**: 
    - This process is **INTERNAL**: It uses only the tree and alignment from your current build.
    - It is **NOT** a comparison with global databases (that was M01/Nextclade). Instead, it creates a historical record of mutational events *specific to your dataset*.
- **Translation**: Converts nucleotide mutations into Amino Acid changes based on the reference annotation (GFF/GFF3).

### Key Instrumentation
- `augur ancestral`
- `augur translate`

---

## M07: Bridge Builder Integration

**Objective**: Calculate the final "Composite Stability Score" to rank targets.

### Technical Process
1.  **Conservation Scanner**: Calculates **Shannon Entropy** ($H$) for every position using a sliding window.
    - $H = -\sum p_i \log_2 p_i$
    - A position with $H=0$ is perfectly conserved.
2.  **Score Integration**: Merges the entropy data with HyPhy selection signals (M05) and Ancestral history (M06).
3.  **Weighted Ranking**: 
    - **Score** = `Conservation (50%)` + `Purifying (25%)` - `Episodic (8%)` - `Recombination (17%)`.

### Key Scripts
- `scripts/m07_bridge_builder.py`

---

## Final Output Checklist
Before declaring a "Gold Zone" for primer design, verify:
1. [ ] **IQ-TREE Likelihood**: Check `logs/m04` for high bootstrap support.
2. [ ] **HyPhy Convergence**: Verify that FUBAR completed without errors.
3. [ ] **M07 Report**: Review `bridge_report.json` for Grade A candidates.
