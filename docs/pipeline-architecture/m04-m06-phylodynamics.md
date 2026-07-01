# Modules M04 - M06: Phylodynamics & Selection

These modules reconstruct the evolutionary history of the viral population, calculate molecular clock rates, and apply phylogenetic models to detect selection pressure at individual codon sites.

---

## M04: Tree Reconstruction & Timetree Calibration

**Objective**: Estimate the phylogenetic tree topology and calibrate branch lengths in calendar years using a molecular clock.

### Technical Details
*   **Phylogenetic Inference**: Uses **IQ-TREE** (Maximum Likelihood method) via `augur tree` to construct a raw phylogenetic tree based on the GTR substitution model.
*   **Molecular Clock Fitting**: Processes the tree through `augur refine`. It resolves polytomies (nodes with multiple branches) and fits a strict or relaxed molecular clock based on sequence collection dates.
*   **Output**: A timetree (`tree.nwk`) where branch lengths represent years, along with estimated substitution rates.

### Debugging Clock Mismatch
If `augur refine` fails, it is usually due to conflicting sequence dates or an incorrect prior `clock_rate`. Run refine manually to inspect:
```bash
augur refine \
  --tree results/m04/tree_raw.nwk \
  --alignment results/m03/aligned.fasta \
  --metadata results/m01/metadata.tsv \
  --output-tree results/m04/tree.nwk \
  --clock-rate 0.003 \
  --clock-std-dev 0.001
```

---

## M05: Evolutionary Selection Pressure (HyPhy)

**Objective**: Identify which codon sites are under evolutionary constraint (purifying selection) or undergoing active mutation (episodic selection).

### Technical Details
*   **Codon Alignment Sanitization**: The aligned nucleotide sequences are converted into a codon-aware alignment (stripping stop codons and correcting partial triplets).
*   **Purifying Selection (FUBAR)**: Runs FUBAR (Fast Unconstrained Bayesian Approximation) to calculate posterior probabilities for negative selection ($\beta < \alpha$). High posterior probability ($> 0.90$) indicates that amino acid changes are rejected by evolution, marking these positions as stable diagnostic targets.
*   **Episodic Selection (MEME)**: Runs MEME (Mixed Effects Model of Evolution) to identify sites subject to episodic diversifying selection (transient evolutionary changes). These are marked as mutational hotspots to be avoided.
*   **Recombination Partitioning (GARD)**: Recombination events swap genomic segments, which can distort selection calculations. GARD partitions the alignment into independent evolutionary blocks to ensure that selection scores are biologically accurate.

---

## M06: Ancestral Reconstruction & Clade Assignment

**Objective**: Reconstruct the genetic sequences of ancestral nodes and assign viral clades to all branches.

### Technical Details
*   **Ancestral Reconstruction**: Uses `augur ancestral` to infer the nucleotide sequence of internal parent nodes in the timetree. This step tracks how mutations arose over time.
*   **Translation**: Translates nucleotide reconstructions to amino acids for ORFs (`augur translate`).
*   **Clade Definitions**: Applies official clade definitions (e.g., GISAID clades or Nextclade lineages) to assign clade names to every node on the tree, allowing clade-specific primer audits.
