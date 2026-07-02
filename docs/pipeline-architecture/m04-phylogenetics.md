---
title: "M04: Phylogenetic Tree Inference & Time-Calibration"
description: "Constructs a maximum-likelihood phylogenetic tree from the M03 alignment using IQ-TREE2 (via augur tree) and calibrates it against sample collection dates to produce a time-resolved phylogeny (via augur refine)."
---

# Module M04: Phylogenetic Tree Inference & Time-Calibration

M04 builds the core evolutionary framework for the pipeline. It executes two sequential steps:

1. **ML tree inference** (`augur tree` → IQ-TREE2): Infers the maximum-likelihood topology from the aligned sequences.
2. **Time-tree calibration** (`augur refine` → TreeTime): Fits a strict molecular clock model to the ML topology, using sample collection dates from the enriched metadata to root and time-scale the tree.

The output time-tree (`tree.nwk`) and branch-length node data (`branch_lengths.json`) propagate downstream to M05 (HyPhy — requires a pruned tree), M06 (ancestral reconstruction), and ultimately M07 (Bridge Builder, which uses per-branch mutation mapping from M06).

---

## Input Requirements

| File | Source | Description |
|---|---|---|
| `results/{BUILD}/m03/aligned.fasta` | M03 | Reference-coordinate-normalised MSA |
| `results/{BUILD}/m01/metadata_enriched.tsv` | M01 | Enriched metadata with `strain` and `date` columns |

### Configuration Parameters

| Config Key | Type | Example | Description |
|---|---|---|---|
| `timetree.clock_rate` | Float | `0.0025` | Prior molecular clock rate (substitutions/site/year) |
| `timetree.clock_std_dev` | Float | `0.0005` | Permitted deviation from the clock-rate prior |
| `augur.threads` | Integer | `8` | CPU threads for IQ-TREE2 |

#### Reference Clock Rates by Pathogen

| Target | Segment | `clock_rate` Range |
|---|---|---|
| Influenza A (H3N2) | HA | `0.0035 – 0.0045` |
| Influenza A (H1N1) | HA | `0.0032 – 0.0042` |
| Influenza A | MP (M1/M2) | `0.0020 – 0.0030` |
| Influenza A | NA | `0.0030 – 0.0040` |
| SARS-CoV-2 | Full genome | `0.0007 – 0.0010` |

<Info>
`clock_rate` is the single most sensitive parameter in the pipeline. Incorrect values cause `augur refine` to fail with negative branch lengths or produce a tree with implausible root date. Consult current Nextstrain builds for your target pathogen/segment before setting this value.
</Info>

---

## Processing Stages

### Stage 1 — ML Tree Inference (`augur_tree`)

**Tool:** IQ-TREE2 (called internally by `augur tree`)

```bash
augur tree \
    --alignment results/{BUILD}/m03/aligned.fasta \
    --output results/{BUILD}/m04/tree_raw.nwk \
    --method iqtree \
    --nthreads 8
```

`augur tree` wraps IQ-TREE2 and automatically selects the GTR+G substitution model (the standard for viral nucleotide alignments). The raw Newick output (`tree_raw.nwk`) is an unrooted, uncalibrated maximum-likelihood topology.

**IQ-TREE2 model selection note:** IQ-TREE2 supports ModelFinder (`-m TEST`) for automatic model selection, but `augur tree` uses GTR+G by default for speed. For publication-grade analyses, consider running IQ-TREE2 directly with `-m MFP` and supplying the resulting tree to `augur refine`.

### Stage 2 — Time-Tree Calibration (`augur_refine`)

**Tool:** TreeTime (called internally by `augur refine`)

```bash
augur refine \
    --tree results/{BUILD}/m04/tree_raw.nwk \
    --alignment results/{BUILD}/m03/aligned.fasta \
    --metadata results/{BUILD}/m01/metadata_enriched.tsv \
    --output-tree results/{BUILD}/m04/tree.nwk \
    --output-node-data results/{BUILD}/m04/branch_lengths.json \
    --timetree \
    --clock-rate 0.0025 \
    --clock-std-dev 0.0005 \
    --stochastic-resolve
```

**Key flags:**

| Flag | Effect |
|---|---|
| `--timetree` | Activates clock-model fitting. Without this, `augur refine` only re-roots the tree without time-calibration. |
| `--clock-rate` | Sets the prior substitution rate. TreeTime uses this as a constraint during regression against sample dates. |
| `--clock-std-dev` | Allowable deviation from `clock_rate`. A tighter value (`0.0002`) enforces stricter clock behaviour; a looser value (`0.001`) accommodates more evolutionary rate variation. |
| `--stochastic-resolve` | Randomly resolves polytomies (multifurcating nodes) into bifurcating clades. Required for many downstream analyses that assume a binary tree. Generates a "Polytomies resolved stochastically" warning that is safe to ignore. |

**TreeTime algorithm overview:**  
TreeTime performs a regression of root-to-tip divergence against sample collection dates. For each taxon, it expects a linear relationship between evolutionary distance from the root and collection date. Outlier sequences (those that violate the clock model beyond `clock_std_dev` tolerance) are flagged in the output logs with `Pruned outliers:`.

---

## Outputs

| File | Path | Description |
|---|---|---|
| Raw ML tree | `results/{BUILD}/m04/tree_raw.nwk` | Unrooted ML topology from IQ-TREE2 |
| **Time-calibrated tree** | `results/{BUILD}/m04/tree.nwk` | **Primary M04 output** — rooted, time-resolved Newick |
| Branch-length node data | `results/{BUILD}/m04/branch_lengths.json` | Per-node clock data for Auspice visualisation |
| IQ-TREE log | `logs/{BUILD}/m04_augur_tree.log` | IQ-TREE2 stdout/stderr |
| TreeTime log | `logs/{BUILD}/m04_augur_refine.log` | augur refine / TreeTime output |

### Reading the Branch Lengths JSON

```json
// results/{BUILD}/m04/branch_lengths.json (excerpt)
{
  "nodes": {
    "A/Indonesia/PDDV0001/2024": {
      "numdate": 2024.673,
      "clock_length": 0.00341,
      "mutation_length": 0.00289,
      "raw_date": "2024-09-02"
    }
  }
}
```

`numdate` is the decimal year (2024.673 ≈ 2 September 2024). `clock_length` is the branch length in time units; `mutation_length` is the raw ML branch length.

---

## Conda Environment

```
envs/nextstrain.yaml
```

Key packages: `augur`, `iqtree`, `treetime`

---

## Debugging Clock Failures

**Symptom:** `augur refine` exits with `ERROR: Negative branch length detected.`

**Diagnosis steps:**

```bash
# 1. Check for date format errors in metadata
awk -F'\t' 'NR>1 {print $NF, $1}' results/{BUILD}/m01/metadata_enriched.tsv | \
    grep -vE '^[0-9]{4}(-[0-9]{2}(-[0-9]{2})?)?$'

# 2. Plot root-to-tip divergence vs date to visualise outliers
# (use Tempest or augur's built-in clock filter)
augur clock-filter \
    --tree results/{BUILD}/m04/tree_raw.nwk \
    --alignment results/{BUILD}/m03/aligned.fasta \
    --metadata results/{BUILD}/m01/metadata_enriched.tsv \
    --clock-rate 0.0025 \
    --output-tree results/{BUILD}/m04/tree_filtered.nwk
```

**Resolution:** If root-to-tip regression fails, try incrementally adjusting `clock_std_dev` upward (e.g., to `0.001`) or removing sequences identified as outliers by `augur clock-filter`.

---

## Known Limitations & Edge Cases

- **Polytomies warning is normal**: The `--stochastic-resolve` flag produces a "Polytomies resolved stochastically" warning when sequences share identical sequence content. This is expected and does not indicate an error.
- **Protein phylogenetics (deferred)**: M04 currently only builds a nucleotide-level phylogeny. Protein-level phylogenetic inference (using WAG/LG substitution models for antigen/vaccine strain selection) is planned for a future phase when serological objectives are in scope.
- **Large dataset memory**: IQ-TREE2 tree inference for >3,000 sequences requires approximately 16–32 GB RAM. If the system runs out of memory, reduce `max_sequences` in M02 config.
- **Reassortment (Influenza)**: Influenza A evolves primarily via reassortment (whole-segment exchange), not intra-segment recombination. The single-segment phylogeny built by M04 represents the evolutionary history of that segment only. Distinct segment-level trees are expected to have different topologies, which is biologically correct.
