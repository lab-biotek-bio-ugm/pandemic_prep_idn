---
title: "M06: Ancestral Sequence Reconstruction & Export"
description: "Reconstructs ancestral nucleotide and amino acid mutations across the phylogenetic tree using joint maximum likelihood inference, then exports the complete annotated dataset as an Auspice-compatible JSON for Nextstrain visualisation."
---

# Module M06: Ancestral Sequence Reconstruction & Export

M06 completes the phylogenetic analysis track. It infers ancestral states at internal tree nodes — determining which nucleotide and amino acid mutations arose on which branches — and packages the entire annotated phylogenetic dataset into the Auspice v2 JSON format for Nextstrain visualisation. The ancestral mutation maps feed into M07 (Bridge Builder) to provide phylogenetic context for conservation scoring.

---

## Input Requirements

| File | Source | Description |
|---|---|---|
| `results/{BUILD}/m04/tree.nwk` | M04 | Time-calibrated phylogenetic tree |
| `results/{BUILD}/m03/aligned.fasta` | M03 | MSA — used for ancestral nucleotide inference |
| `results/{BUILD}/m04/branch_lengths.json` | M04 | Branch-length node data for Auspice export |
| `results/{BUILD}/m01/metadata_enriched.tsv` | M01 | Sample metadata for Auspice colour tracks |
| `config.references.annotation` | Config | GFF3 annotation file — used for CDS gene coordinate extraction and `augur translate` |

### Configuration Parameters

| Config Key | Type | Notes |
|---|---|---|
| `references.annotation` | String | Path to GFF3 annotation (e.g., `config/reference_NC_002016.gff`) |
| `nextclade.enabled` | Boolean | If `true`, Auspice colour-by tracks include `clade` and `lineage` |

---

## Processing Stages

### Stage 1 — Nucleotide Ancestral Reconstruction (`augur_ancestral`)

**Tool:** IQ-TREE2 ancestral reconstruction (called internally by `augur ancestral`)

```bash
augur ancestral \
    --tree results/{BUILD}/m04/tree.nwk \
    --alignment results/{BUILD}/m03/aligned.fasta \
    --output-node-data results/{BUILD}/m06/nt_muts.json \
    --inference joint
```

**`--inference joint`** uses **joint maximum likelihood** reconstruction, which simultaneously optimises all ancestral states across the entire tree. This is more accurate than marginal (independent per-node) reconstruction but requires more memory.

The output `nt_muts.json` records, for each internal node and tip, the list of nucleotide mutations relative to its parent node. Format: `A187T` means adenine at position 187 changed to thymine.

### Stage 2 — Amino Acid Translation (`augur_translate`)

**Tool:** `augur translate`

Translates per-branch nucleotide mutations (`nt_muts.json`) into amino acid changes for each annotated CDS gene.

**Zero hardcoding of gene names:** Gene names are dynamically extracted from the GFF3 annotation at workflow load time via `_parse_cds_genes_from_gff()`. Swapping the GFF3 file automatically changes which genes appear in the output without any config modification.

```bash
# Example for Influenza A MP (genes: M1, M2)
augur translate \
    --tree results/{BUILD}/m04/tree.nwk \
    --ancestral-sequences results/{BUILD}/m06/nt_muts.json \
    --reference-sequence config/reference_NC_002016.gff \
    --genes M1 M2 \
    --validation-mode warn \
    --output-node-data results/{BUILD}/m06/aa_muts.json
```

**`--validation-mode warn`** prevents stop codons in minority sequences from halting the pipeline. Stop codons are logged as warnings but do not cause a failure.

### Stage 3 — Auspice Dataset Export (`augur_export`)

**Tool:** `augur export v2`

Assembles all node data files and metadata into a single Auspice v2 JSON for Nextstrain visualisation:

```bash
augur export v2 \
    --tree results/{BUILD}/m04/tree.nwk \
    --metadata results/{BUILD}/m01/metadata_enriched.tsv \
    --node-data \
        results/{BUILD}/m04/branch_lengths.json \
        results/{BUILD}/m06/nt_muts.json \
        results/{BUILD}/m06/aa_muts.json \
    --color-by-metadata clade lineage country region date \
    --output results/{BUILD}/nextstrain/auspice_dataset.json \
    --title "Pandemic Preparedness: {BUILD}" \
    --build-url "https://github.com/engkinandatama/pandemic-preparedness"
```

**Colour-by metadata fields:**
- If `nextclade.enabled: true` → `clade lineage country region date`
- If `nextclade.enabled: false` → `subtype clade lineage country region host date`

---

## Outputs

| File | Path | Description |
|---|---|---|
| **NT mutation node data** | `results/{BUILD}/m06/nt_muts.json` | Per-branch nucleotide mutation list |
| **AA mutation node data** | `results/{BUILD}/m06/aa_muts.json` | Per-branch amino acid mutation list per gene |
| Ancestral log | `logs/{BUILD}/m06_augur_ancestral.log` | augur ancestral stderr |
| Translate log | `logs/{BUILD}/m06_augur_translate.log` | augur translate stderr |
| **Auspice JSON** | `results/{BUILD}/nextstrain/auspice_dataset.json` | Complete visualisation dataset for Nextstrain |
| Export log | `logs/{BUILD}/m06_augur_export.log` | augur export stderr |

### Reading `nt_muts.json`

```json
// results/{BUILD}/m06/nt_muts.json (excerpt)
{
  "nodes": {
    "A/Indonesia/PDDV0001/2024": {
      "muts": ["A187T", "G452A"],
      "aa_muts": {}
    },
    "NODE_0001234": {
      "muts": ["C891T"],
      "aa_muts": {
        "M1": ["I297V"]
      }
    }
  }
}
```

`muts` contains nucleotide mutations relative to the parent node; `aa_muts` contains amino acid changes per gene.

---

## Visualising with Nextstrain Auspice

To load the exported dataset in a local Auspice instance:

```bash
# Install Auspice globally (one time)
npm install --global auspice

# Serve the dataset locally
auspice view --datasetDir results/{BUILD}/nextstrain/
```

Open `http://localhost:4000` in a browser. The tree will be interactive with colour-by options for country, clade, lineage, and date.

---

## Conda Environment

```
envs/nextstrain.yaml
```

Key packages: `augur`, `auspice`, `iqtree`

---

## Known Limitations & Edge Cases

- **Deferred: Protein phylogeny (M04-1)** — A separate protein-level phylogenetic tree using IQ-TREE2 with protein substitution models (WAG, LG) for antigen/serological applications is planned for Phase 2. Currently, protein sequences are used only for HyPhy dN/dS analysis (M05), not for building protein phylogenies.
- **Stop codons in `augur translate`**: Some sequences may contain in-frame stop codons due to sequencing errors. `--validation-mode warn` allows the pipeline to continue, but these sequences may appear with truncated amino acid annotations in Auspice. Review `logs/{BUILD}/m06_augur_translate.log` for warnings.
- **Auspice colour-by fields require matching metadata columns**: The `--color-by-metadata` fields must exist as columns in `metadata_enriched.tsv`. If Nextclade was disabled (`nextclade.enabled: false`) but `clade` is listed in the colour-by fields, `augur export` will produce a warning and omit that field.
- **`nt_muts.json` is used by M07**: The Bridge Builder script (`m07_bridge_builder.py`) reads `nt_muts.json` to map branch-specific mutation events back to reference positions, providing phylogenetic context beyond simple column-level entropy. Ensure M06 completes successfully before M07 runs.
