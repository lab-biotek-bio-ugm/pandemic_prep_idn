---
title: "M05: Evolutionary Selection Pressure Analysis"
description: "Runs HyPhy's FUBAR, MEME, and GARD methods on gene-sliced codon alignments to quantify pervasive purifying selection, episodic diversifying selection, and recombination breakpoints. These signals are consumed by the Bridge Builder (M07) as scoring components."
---

# Module M05: Evolutionary Selection Pressure Analysis

M05 is the evolutionary bioinformatics core of the pipeline. It generates **per-site selection pressure scores** that the Bridge Builder (M07) uses to distinguish diagnostically stable genomic regions from those undergoing active mutational pressure.

The module executes five sequential stages for each target gene:
1. Diversity-aware subsampling of the MSA (for HyPhy computational tractability)
2. Gene extraction and codon-level sanitisation
3. Phylogenetic tree pruning
4. FUBAR — pervasive purifying and positive selection
5. MEME — episodic diversifying selection
6. GARD — recombination breakpoint detection (conditional)

---

## Input Requirements

| File | Source | Description |
|---|---|---|
| `results/{BUILD}/m03/aligned.fasta` | M03 | Full reference-coordinate alignment (Track A, entropy) |
| `results/{BUILD}/m04/tree.nwk` | M04 | Time-calibrated phylogenetic tree |
| `results/{BUILD}/m01/metadata_enriched.tsv` | M01 | Clade/subtype metadata for subsample stratification |
| `config.references.annotation` | Config | GFF3 annotation file for gene coordinate extraction |

### Configuration Parameters

| Config Key | Type | Default | Description |
|---|---|---|---|
| `hyphy.methods` | List | `["FUBAR", "MEME"]` | HyPhy methods to run. Add `"GARD"` for non-segmented viruses. |
| `hyphy.threads` | Integer | `32` | CPU threads for all HyPhy methods |
| `hyphy.subsample.enabled` | Boolean | `false` | Enable diversity-aware subsampling before HyPhy |
| `hyphy.subsample.max_sequences` | Integer | `300` | Maximum sequences for HyPhy input (recommended ≤500) |
| `hyphy.subsample.min_sequences` | Integer | `50` | Minimum sequences required; pipeline aborts if below this |
| `hyphy.subsample.strategy` | String | `"balanced"` | Sampling strategy: `"balanced"` or `"random"` |
| `hyphy.subsample.seed` | Integer | `42` | Random seed for reproducibility |

---

## Processing Stages

### Stage 1 — Diversity-Aware Subsampling (`m05_hyphy_subsample`)

**Script:** `scripts/hyphy_subsample.py`  
**Environment:** `envs/hyphy_macse.yaml`

HyPhy's likelihood models become computationally intractable for alignments exceeding ~500 sequences. This stage reduces the input to a representative subset using adaptive stratified sampling:

**Stratification hierarchy (adaptive fallback):**
1. `subtype × clade × country × year` (most specific)
2. `subtype × clade × country` (if grouping 1 yields too few)
3. `subtype × clade` (further fallback)
4. Random draw (last resort)

<Note>
**Track A (M07 Bridge Builder entropy scan) uses the full `aligned.fasta` directly**, not this subsample. The subsample only feeds the HyPhy selection analysis (Track B).
</Note>

```bash
python scripts/hyphy_subsample.py \
    --alignment results/{BUILD}/m03/aligned.fasta \
    --metadata results/{BUILD}/m01/metadata_enriched.tsv \
    --output results/{BUILD}/m05/aligned_hyphy.fasta \
    --report results/{BUILD}/m05/hyphy_subsample_report.json \
    --max-sequences 300 \
    --min-sequences 50 \
    --strategy balanced \
    --seed 42
```

If `hyphy.subsample.enabled: false`, the full `aligned.fasta` is copied as-is and a bypass report is written.

### Stage 2 — Gene Extraction & Codon Sanitisation (`m05_extract_genes`)

**Script:** `scripts/extract_genes.py`  
**Environment:** `envs/hyphy_macse.yaml`

Slices coding sequences from the subsampled alignment using GFF3 coordinates. Gene names are **dynamically parsed** from the GFF3 at workflow load time — no hardcoded gene lists. The GFF3 parser reads `CDS` features and extracts gene names from `gene=`, `Name=`, or `locus_tag=` attributes in that priority order.

For each gene (e.g., `M1`, `M2` for Influenza A MP segment):
- Nucleotide positions are extracted from the alignment according to GFF3 CDS coordinates
- A **harsh masking** step removes:
  - Internal stop codons (which would cause HyPhy to crash)
  - Frame-shift insertions relative to the reference
  - Columns that are >95% gaps

Outputs: `results/{BUILD}/m05/cleaned_{GENE}_NT.fasta` for each gene.

### Stage 3 — Tree Pruning (`m05_prune_tree`)

**Script:** `scripts/prune_tree.py`

The full M04 time-tree may contain taxa that were removed during HyPhy subsampling. HyPhy requires that the tree tips exactly match the alignment sequence IDs. This rule prunes the time-tree to retain only the tips present in each gene-specific FASTA.

Output: `results/{BUILD}/m05/tree_{GENE}.nwk`

### Stage 4 — FUBAR: Pervasive Selection (`m05_hyphy_fubar`)

**Tool:** HyPhy FUBAR (Fast Unconstrained Bayesian AppRoximation)  
**Environment:** `envs/hyphy_macse.yaml`

FUBAR uses a Bayesian grid-approximation approach to compute the **posterior probability** that each codon site is under:
- **Pervasive negative (purifying) selection**: $P(\alpha > \beta)$ — evolutionarily constrained sites
- **Pervasive positive selection**: $P(\beta > \alpha)$ — sites under adaptive pressure across the entire tree

```bash
# Stale cache is removed before every run to prevent dimension mismatch errors
rm -f results/{BUILD}/m05/cleaned_{GENE}_NT.fasta.FUBAR.cache

hyphy CPU=32 fubar \
    --alignment results/{BUILD}/m05/cleaned_{GENE}_NT.fasta \
    --tree results/{BUILD}/m05/tree_{GENE}.nwk \
    --output results/{BUILD}/m05/fubar_{GENE}.json
```

<Info>
The cache-clearing step before each FUBAR run is critical. HyPhy persists its MCMC cache to disk. If the alignment dimensions change between runs (e.g., after subsampling parameter change), a stale cache causes a fatal dimension mismatch error. Always clear the cache when re-running with modified inputs.
</Info>

**Output schema (`fubar_{GENE}.json` — key fields):**

```json
{
  "MLE": {
    "headers": ["alpha", "beta", "beta-alpha", "Prob[alpha>beta]", "Prob[beta>alpha]"],
    "content": {
      "0": [[0.123, 0.001, -0.122, 0.97, 0.03], ...]
    }
  }
}
```

The Bridge Builder (M07) reads `Prob[alpha>beta]` as the **purifying selection score** $P$ for each codon position.

### Stage 5 — MEME: Episodic Selection (`m05_hyphy_meme`)

**Tool:** HyPhy MEME (Mixed Effects Maximum Likelihood Evolution)  
**Environment:** `envs/hyphy_macse.yaml`

MEME detects **episodic positive selection** — sites that experienced positive selection in only a subset of branches (lineages), while being conserved on most branches. This is the dominant mode of diversifying selection in RNA viruses.

```bash
rm -f results/{BUILD}/m05/cleaned_{GENE}_NT.fasta.MEME.cache

hyphy CPU=32 meme \
    --alignment results/{BUILD}/m05/cleaned_{GENE}_NT.fasta \
    --tree results/{BUILD}/m05/tree_{GENE}.nwk \
    --output results/{BUILD}/m05/meme_{GENE}.json
```

MEME sites with p-value < `hyphy.p_value_threshold` (default `0.1`) are flagged as episodic selection hotspots. M07 applies a **penalty score** ($E$) to these positions to disfavour them as primer binding sites.

**Output schema (`meme_{GENE}.json` — key fields):**

```json
{
  "MLE": {
    "headers": ["alpha", "beta-", "p-", "beta+", "p+", "LRT", "p-value", "# branches"],
    "content": {
      "0": [[0.12, 0.0, 1.0, 0.89, 0.23, 4.67, 0.031, 2], ...]
    }
  }
}
```

M07 reads the `p-value` column to determine episodic penalty positions.

### Stage 6 — GARD: Recombination Detection (`m05_hyphy_gard`)

**Tool:** HyPhy GARD (Genetic Algorithm for Recombination Detection)  
**Script wrapper:** `scripts/run_hyphy_gard_safe.py`  
**Environment:** `envs/hyphy_macse.yaml`

GARD detects intra-sequence recombination breakpoints using a genetic algorithm that optimises c-AIC across candidate breakpoint models.

<Warning>
**GARD is disabled for all Influenza A/B configurations by default.** Influenza evolves primarily via reassortment (whole-segment swapping), not intra-segment recombination. Running GARD on Influenza produces computationally expensive analyses that detect noise rather than true biological signal (Boni et al. 2010; Worobey & Holmes 1999).

**Enable GARD only for non-segmented RNA viruses** (e.g., SARS-CoV-2, Dengue, Zika) by adding `"GARD"` to `hyphy.methods`.
</Warning>

The `run_hyphy_gard_safe.py` wrapper performs pre-flight validation:
- If `alignment_length < 4 × num_sequences` (c-AIC requirement), writes an empty JSON and exits 0 (soft-fail — pipeline continues).
- If GARD crashes at runtime, catches the exception and soft-fails.

```bash
python scripts/run_hyphy_gard_safe.py \
    results/{BUILD}/m05/cleaned_{GENE}_NT.fasta \
    results/{BUILD}/m05/gard_{GENE}.json \
    --threads 32 \
    --rv GDD
```

---

## Outputs

| File | Path | Description |
|---|---|---|
| HyPhy input MSA | `results/{BUILD}/m05/aligned_hyphy.fasta` | Diversity-subsampled alignment |
| Subsample report | `results/{BUILD}/m05/hyphy_subsample_report.json` | Subsampling statistics |
| Gene FASTA per gene | `results/{BUILD}/m05/cleaned_{GENE}_NT.fasta` | Codon-sanitised per-gene alignment |
| Pruned tree per gene | `results/{BUILD}/m05/tree_{GENE}.nwk` | Tip-pruned phylogeny for each gene |
| **FUBAR results** | `results/{BUILD}/m05/fubar_{GENE}.json` | **Per-site purifying/positive selection probabilities** |
| **MEME results** | `results/{BUILD}/m05/meme_{GENE}.json` | **Per-site episodic selection p-values** |
| GARD results | `results/{BUILD}/m05/gard_{GENE}.json` | Recombination breakpoint positions (if enabled) |

---

## Conda Environment

```
envs/hyphy_macse.yaml
```

Key packages: `hyphy`, `macse`, `biopython`, `pandas`

---

## Interpreting Selection Results

| Score | Interpretation | M07 Action |
|---|---|---|
| FUBAR `Prob[alpha>beta]` > 0.90 | Strong pervasive purifying selection | Positive contribution to composite score ($P$) |
| FUBAR `Prob[beta>alpha]` > 0.90 | Pervasive positive selection | Negative signal — site likely escaping immune pressure |
| MEME p-value < 0.1 | Episodic diversifying selection detected | Penalty applied to composite score ($E$) |
| GARD breakpoint | Recombination boundary detected | Recombination proximity penalty ($R$) applied |

---

## Known Limitations & Edge Cases

- **Dropped: HyPhy FEL** — FEL (Fixed Effects Likelihood) detects the same pervasive selection signal as FUBAR but with a frequentist approach. It is only advantageous for very small datasets (<20 sequences). FUBAR was retained because it outperforms FEL for the typical dataset size range of this pipeline (50–1000 sequences).
- **Deferred: Clade-defining SNP finder (M05-1)** — Planned feature to identify mutations unique to each clade from `nt_muts.json`. Required for Route B (clade-specific primer targets). Priority: Phase 2.
- **GARD computational cost**: On a 1,000 nt gene alignment with 300 sequences, GARD can take 12–48 hours on a single machine. Always submit GARD runs to a compute cluster (HPC/Slurm).
- **Codon frame requirement**: FUBAR and MEME operate on codon models and require that the input alignment length is divisible by 3. The `extract_genes.py` sanitiser enforces this, but verify by checking that `len(sequence) % 3 == 0` for all sequences in `cleaned_{GENE}_NT.fasta`.
