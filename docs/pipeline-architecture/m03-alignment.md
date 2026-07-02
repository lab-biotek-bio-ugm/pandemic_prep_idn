---
title: "M03: Multiple Sequence Alignment"
description: "Performs reference-based multiple sequence alignment using MAFFT (via augur align), producing a coordinate-normalised alignment that forms the shared input for both the phylogenetic track (M04) and the entropy conservation track (M07)."
---

# Module M03: Multiple Sequence Alignment

M03 aligns the subsampled sequences from M02 against a fixed reference genome using MAFFT in reference-guided mode (invoked via `augur align`). The resulting alignment establishes a unified coordinate system where every column corresponds to a specific position in the reference sequence. This coordinate normalisation is a prerequisite for both:
- **Track A** (M07 Bridge Builder): per-position Shannon entropy calculation requires aligned coordinates.
- **Track B** (M04 → M05 phylogenetics): IQ-TREE and HyPhy require gap-aware codon-position alignments.

> [!NOTE]
> trimAl (alignment column trimming) was explicitly **dropped** from the pipeline design. Reference-guided alignment via `--fill-gaps` inherently produces cleaner columns than de novo alignment, and trimming would risk removing conserved regions critical for primer design. Codon-level sanitisation for HyPhy is handled separately in M05 by `clean_alignment_for_hyphy.py`.

---

## Input Requirements

| File | Source | Description |
|---|---|---|
| `results/{BUILD}/m02/filtered_sequences.fasta` | M02 | Subsampled, deduplicated sequence set |
| `config.references.genome` | Config | Reference GenBank or FASTA file (e.g., `config/reference_NC_002016.gb`) |

### Configuration Parameters

| Config Key | Type | Default | Description |
|---|---|---|---|
| `references.genome` | String | — | Path to reference genome file (GenBank `.gb` or FASTA) |
| `augur.threads` | Integer | `8` | Number of threads passed to MAFFT via `augur align --nthreads` |

---

## Processing

### Rule: `augur_align`

**Environment:** `envs/nextstrain.yaml`  
**Tool:** MAFFT (called internally by `augur align`)

```bash
augur align \
    --sequences results/{BUILD}/m02/filtered_sequences.fasta \
    --reference-sequence config/reference_NC_002016.gb \
    --output results/{BUILD}/m03/aligned.fasta \
    --fill-gaps \
    --nthreads 8
```

**Key flags:**

| Flag | Effect |
|---|---|
| `--reference-sequence` | Forces MAFFT into reference-guided mode. All sequences are aligned to this reference, not to each other. |
| `--fill-gaps` | After alignment, gaps at positions relative to the reference that represent missing sequence (e.g., assembly gaps at sequence termini) are filled with `N`. This prevents gaps from being misinterpreted as true deletions by downstream tools. |
| `--nthreads` | Parallelises MAFFT across available CPU cores. Derived from `config.augur.threads`. |

**What reference-guided alignment means in practice:**  
Rather than building a full multiple sequence alignment from scratch (which is $O(N^2)$ in time), MAFFT first aligns each query sequence pairwise to the reference, then reconciles insertions relative to reference coordinates. Insertions that are present in query sequences but absent in the reference are removed from the output alignment. This preserves a consistent coordinate system (position 1 in the output always corresponds to position 1 of the reference).

---

## Downstream Usage of `aligned.fasta`

The `aligned.fasta` file produced by M03 is a **shared, non-modified input** used by multiple downstream modules:

| Consumer | Purpose |
|---|---|
| M04 (`augur tree`) | Provides the sequence data for IQ-TREE phylogenetic inference |
| M04 (`augur refine`) | Used alongside the raw tree for time-tree calibration |
| M05 (`m05_hyphy_subsample`) | Subsampled and sliced into gene-specific FASTA for HyPhy |
| M06 (`augur ancestral`) | Source alignment for ancestral state reconstruction |
| M07 (`m07_bridge_builder.py`) | Scanned at per-nucleotide resolution for Shannon entropy |

---

## Outputs

| File | Path | Description |
|---|---|---|
| **Aligned FASTA** | `results/{BUILD}/m03/aligned.fasta` | **Primary M03 output** — reference-coordinate-normalised MSA |
| Alignment log | `logs/{BUILD}/m03_augur_align.log` | MAFFT stdout/stderr |

### Sanity-Checking the Alignment

After alignment, verify that all sequences have the same length as the reference:

```bash
# Check that all sequences have identical length (should return 1 unique value)
awk '/^>/{next} {print length($0)}' results/{BUILD}/m03/aligned.fasta | sort -u

# Count sequences in alignment
grep -c '^>' results/{BUILD}/m03/aligned.fasta
```

Expected output: one single length value equal to the reference genome length, and a sequence count matching the M02 output count.

---

## Conda Environment

```
envs/nextstrain.yaml
```

Key packages: `augur`, `mafft`

---

## Known Limitations & Edge Cases

- **Sequences diverged from reference**: Highly divergent sequences (e.g., inter-subtype Influenza A variants) may produce poorly-aligned regions when forced into a single reference coordinate system. Review the alignment manually in a tool like AliView or UGENE if M07 reports unexpectedly high entropy across the entire alignment.
- **Reference must match segment**: Ensure `references.genome` corresponds to the correct viral segment (e.g., `NC_002016.1` for Influenza A M segment). Aligning HA sequences against an NA reference will produce a meaningless alignment without error.
- **Gaps vs. missing data**: The `--fill-gaps` flag converts terminal gaps to `N`. This means 3' and 5' truncated sequences are padded rather than having real data. The N-content filter in M02 (`max_n_content`) helps exclude sequences that are too fragmented, but some padding remains unavoidable for partially complete assemblies.
- **Memory scaling**: MAFFT reference-guided mode scales approximately linearly in memory with sequence length × count. For whole-genome alignments (>10 kb) with >5,000 sequences, ensure at least 16 GB RAM is available.
