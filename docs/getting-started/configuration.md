# Configuration Guide

This guide explains how to configure the pipeline using the YAML configuration profiles (e.g., `config/config_influenza_a_mp.yaml`). Properly calibrating these values is essential for adapting the pipeline to new pathogens, segments, or geographic regions.

For a complete list of every available configuration key and its default value, see the [Configuration Reference](/reference/configuration-reference).

---

## 1. Quality Control & Subsampling Calibration

The subsampling parameters control how sequences are filtered and balanced before tree inference. Properly setting these values avoids computational bottlenecks in downstream modules.

```yaml
qc:
  max_n_content: 5.0
subsampling:
  min_length: 800
  max_sequences: 3000
  group_by: "subtype country year month"
  sequences_per_group: 10
```

### Setting `min_length`
- **Rule of thumb**: Set `min_length` to **80%–90%** of the full segment size.
- **Example**: The Influenza A MP segment is 1,027 bp. We set `min_length: 800` to discard short sequencing reads or partial PCR fragments that would introduce massive gap columns in the alignment.
- For whole-genome targets (like SARS-CoV-2, 29.9 kb), set `min_length` to `25000` or higher to select only complete genome assemblies.

### Managing Subsampling Density
- If your input dataset is very large (>10,000 sequences), setting `max_sequences: null` with a dense `sequences_per_group` will cause the phylogenetic tree inference (M04) to crash due to running out of memory.
- Enforce a hard cap using `max_sequences: 3000` (or `2000` for faster runs) to ensure the alignment and tree building finish in a reasonable time.

---

## 2. Molecular Clock Rate Calibration (M04)

`timetree` parameters control how TreeTime fits a strict molecular clock to the Maximum Likelihood tree topology.

```yaml
timetree:
  clock_rate: 0.003
  clock_std_dev: 0.0015
```

### Reference Clock Rates by Pathogen
The molecular clock rate prior (`clock_rate`) is the most sensitive parameter in the pipeline. Using an incorrect rate will cause time-calibration to fail with negative branch length errors.

| Pathogen / Segment | `clock_rate` Prior | `clock_std_dev` | Rationale |
|---|---|---|---|
| Influenza A MP | `0.003` | `0.0015` | Conserved internal segment, slower rate. |
| Influenza A NP | `0.003` | `0.0015` | Conserved internal nucleoprotein. |
| Influenza A HA | `0.004` | `0.0010` | Highly dynamic surface antigen, fast rate. |
| SARS-CoV-2 | `0.0008` | `0.0002` | Proofread RNA virus, very slow rate. |

If `augur refine` crashes with clock fitting errors:
1. Verify that all dates in the metadata are valid ISO dates (`YYYY-MM-DD`). Look for placeholder dates like `2024-XX-XX` and verify they are handled correctly.
2. Try relaxing `clock_std_dev` to a larger value (e.g., `0.0020`) to accommodate higher rate variation among lineages.

---

## 3. Bridge Builder Scoring Weights (M07)

The `bridge.scoring` weights balance the importance of historical sequence conservation against current evolutionary selection pressure.

```yaml
bridge:
  scoring:
    core:
      conservation_weight: 0.66
      purifying_weight: 0.34
```

### Adjusting the Balance
- **Default (0.66 / 0.34)**: Prioritizes physical conservation while utilizing purifying selection as a verification signal. This is optimal for routine universal diagnostic design.
- **Strict Evolutionary Focus (0.50 / 0.50)**: Increases the weight of purifying selection. Recommended when designing assays for highly variable regions (like HA or S genes) where physical conservation is low, forcing the design to rely on functional constraint.
- **Pure Conservation (1.00 / 0.00)**: Disables selection pressure analysis. The composite score is determined purely by Shannon Entropy. Use only if HyPhy is disabled or fails to run.

---

## 4. PrimerLab Degeneracy Bounds (M08)

Degenerate bases allow primers to bind to multiple variant sequences, but excessive degeneracy dilutes the active primer concentration in the reaction tube.

```yaml
primerlab:
  design:
    degenerate:
      enabled: true
      maf_threshold: 0.10
      max_multiplier: 16
```

### Calibration Rules
- **`maf_threshold`**: Set to `0.10` (10% Minor Allele Frequency). Any position where a minor allele is present in at least 10% of the population will receive a degenerate base.
- **`max_multiplier`**: Set to `16` or `32` maximum. A multiplier of 16 means there are 16 unique sequence combinations in the primer pool. Going beyond 32 will lead to a **primer dilution effect** and cause amplification failure due to non-specific binding and primer-dimers.
