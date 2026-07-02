---
title: "Supported Pathogens"
description: "Reference profiles for all pathogens currently validated in the pipeline, including reference accessions, segment configurations, recommended clock rates, and Nextclade dataset identifiers."
---

# Supported Pathogens

This page documents all pathogen configurations validated and tested in the pipeline. Each profile includes the reference genome, annotation, molecular clock parameters, and key considerations for that target.

---

## Influenza A Virus — MP Segment (Matrix Protein)

**Use case:** Universal influenza A detection across all subtypes. The MP (matrix protein) segment is the most conserved segment of the influenza A genome, encoding the M1 (matrix protein) and M2 (ion channel / proton pump) proteins.

| Parameter | Value |
|---|---|
| Reference accession | `NC_002016.1` |
| Reference file | `config/references/flu_a_mp_reference.fasta` |
| Annotation | `config/references/flu_a_mp_annotation.gff` |
| Genes (from GFF3) | `M1`, `M2` |
| Segment length | 1,027 bp |
| Clock rate | `0.003` subst/site/year |
| Clock std dev | `0.0015` |
| Nextclade dataset | Not typically used (nextclade: disabled) |
| Config profile | `config/config_influenza_a_mp.yaml` |
| Build name convention | `influenza-a-mp-{region}` |

**Why MP?** M1 contains the majority of the segment (nt 26–784) and is the most sequence-conserved coding region in the influenza A genome due to its structural role in virion assembly. M2 (nt 26–51, 740–1007; spliced product) is critical for viral uncoating and is the target of adamantane antivirals — making it a diagnostically and clinically important locus. The full segment is typically used as a single design target (`strategy: full_segment`).

**Recommended Bridge settings:**

```yaml
bridge:
  strategy: "full_segment"
  entropy_threshold: 0.020     # MP is highly conserved — strict threshold appropriate
  design_regions:
    - "segment-MP"
  tiers:
    strict:
      min_conservation: 0.95
      min_purifying_prob: 0.9
```

---

## Influenza A Virus — NP Segment (Nucleoprotein)

**Use case:** Alternative universal influenza A detection, particularly useful when cross-reactivity with M2 ion channel inhibitor sequences is a concern.

| Parameter | Value |
|---|---|
| Reference accession | `NC_002019.1` |
| Reference file | `config/references/flu_a_np_reference.fasta` |
| Annotation | `config/references/flu_a_np_annotation.gff` |
| Genes (from GFF3) | `NP` |
| Segment length | 1,565 bp |
| Clock rate | `0.003` subst/site/year |
| Clock std dev | `0.0015` |
| Nextclade dataset | Not typically used |
| Config profile | `config/config_influenza_a_np.yaml` |
| Build name convention | `influenza-a-np-{region}` |
| `subsampling.min_length` | `1200` (longer than MP due to segment size) |

**Why NP?** The nucleoprotein is the primary target for type-specific detection (influenza A vs. B) due to strong type-level conservation with lower subtype-level variability than the surface proteins (HA, NA). At 1,565 bp, it offers more primer design real-estate than MP. The single `NP` CDS spans the full coding sequence. Recommended `bridge.strategy: "windowed"` with `window_size: 150` for NP to identify the most stable sub-regions.

**Recommended Bridge settings:**

```yaml
bridge:
  strategy: "windowed"
  window_size: 150
  entropy_threshold: 0.020
  min_region_length: 150
```

---

## Influenza A Virus — HA Segment (Haemagglutinin)

**Use case:** Subtype-specific detection and phylogenetic surveillance (H1N1, H3N2, H5N1, etc.). HA is the primary antigenic protein and the most rapidly evolving surface protein.

| Parameter | Value |
|---|---|
| Reference accession | `CY121680.1` (H3N2), `MK626062.1` (H1N1pdm09) |
| Segment length | ~1,701 bp (H3), ~1,704 bp (H1) |
| Clock rate | `0.004 – 0.005` subst/site/year |
| Clock std dev | `0.001` |
| Nextclade dataset | `nextstrain/flu/h3n2/ha`, `nextstrain/flu/h1n1pdm/ha/MW626062` |
| Build name convention | `influenza-a-h3n2-ha-{region}` |

<Warning>
HA is **highly variable**. Expect significantly lower composite Bridge scores compared to MP/NP. `bridge.entropy_threshold` should be relaxed to `0.05–0.10`. Most high-scoring regions for universal assay design in HA are concentrated in the stem (HA2 domain, stalk region), not the head (HA1 domain). The head is the primary antibody target and under intense positive selection.
</Warning>

**Config notes:**
- Set `nextclade.enabled: true` and configure `nextclade.dataset_name` to assign H3N2 or H1N1 clade/lineage labels.
- Use `subsampling.query` to restrict to a single subtype to avoid mixing H3N2 and H1N1 in the same alignment.
- GARD (`hyphy.methods: ["FUBAR", "MEME", "GARD"]`) is **NOT recommended** for segmented viruses including HA — reassortment, not intra-segment recombination, is the dominant evolutionary mechanism.

---

## SARS-CoV-2

**Use case:** Whole-genome or gene-specific diagnostic design for the SARS-CoV-2 RNA genome. Larger genome (≈29.9 kb) requires specific scaling considerations.

| Parameter | Value |
|---|---|
| Reference accession | `NC_045512.2` (Wuhan-Hu-1) |
| Segment length | 29,903 bp |
| Clock rate | `0.0008` subst/site/year |
| Clock std dev | `0.0002` |
| Nextclade dataset | `sars-cov-2` |
| Build name convention | `sars-cov-2-{gene}-{region}` |

<Info>
For SARS-CoV-2, always restrict analysis to a **single gene** (e.g., `N`, `ORF1ab`, `E`) rather than the full genome. Aligning 29.9 kb sequences with >3,000 sequences requires ≥32 GB RAM and generates very large intermediate files. Use `hyphy.target_genes` to select only the target gene rather than running FUBAR/MEME on all 10+ CDS genes simultaneously.
</Info>

**Key differences from Influenza:**

| Property | Influenza A | SARS-CoV-2 |
|---|---|---|
| Evolution mechanism | Reassortment (segment-level) | Mutation + recombination (intra-genomic) |
| GARD | Disabled (segmented) | **Recommended** (`hyphy.methods: ["FUBAR", "MEME", "GARD"]`) |
| Clock rate | ~0.003 subst/site/yr | ~0.0008 subst/site/yr |
| Genome size | 0.9–2.3 kb per segment | 29.9 kb |
| `min_length` | 800–1200 bp | 25,000 bp (for complete genomes) |
| Recommended strategy | `full_segment` | `windowed` (window_size: 300) |

**Recommended config snippet for SARS-CoV-2 N gene:**

```yaml
build_name: "sars-cov-2-n-global"
references:
  reference_id: "NC_045512.2"
  genome: "config/references/sars2_reference.gb"
  annotation: "config/references/sars2_annotation.gff"
hyphy:
  methods:
    - "FUBAR"
    - "MEME"
    - "GARD"       # Enabled for non-segmented virus
  target_genes:
    - "N"          # Process only N gene, not all 10+ ORFs
  subsample:
    max_sequences: 500   # Cap for tractable GARD runtime
timetree:
  clock_rate: 0.0008
  clock_std_dev: 0.0002
bridge:
  strategy: "windowed"
  window_size: 300
  entropy_threshold: 0.050   # More permissive than Influenza
nextclade:
  enabled: true
  dataset_name: "sars-cov-2"
```

---

## Adding a New Pathogen

To configure the pipeline for a pathogen not listed above:

1. **Identify the reference genome** — use a complete, high-quality annotated GenBank record (NCBI RefSeq preferred).
2. **Download reference files:**
   ```bash
   # Download GenBank record
   efetch -db nuccore -id NC_XXXXXX.X -format gb \
       > config/references/target_reference.gb
   
   # Convert GenBank CDS annotations to GFF3
   # (use Biopython or bp_seqfeature_load.pl)
   python scripts/gb_to_gff3.py \
       config/references/target_reference.gb \
       > config/references/target_annotation.gff
   ```
3. **Determine the molecular clock rate** — consult published Nextstrain builds or primary literature for the target pathogen. Key parameters: RdRp error rate, generation time, transmission bottleneck size.
4. **Create a new config profile** from the closest existing template (`config_influenza_a_mp.yaml` for segmented RNA viruses, SARS-CoV-2 snippet above for non-segmented).
5. **Determine GARD eligibility** — enable GARD only for non-segmented RNA viruses where intra-genomic recombination is documented.
6. **Run a dry-run to verify the DAG:**
   ```bash
   snakemake --use-conda \
       --configfile config/config_NEW_TARGET.yaml \
       --dryrun --quiet
   ```

---

## Molecular Clock Rate Reference Table

| Pathogen | Segment/Gene | Published Rate (subst/site/yr) | Pipeline `clock_rate` |
|---|---|---|---|
| Influenza A H3N2 | HA | 0.0037 | `0.004` |
| Influenza A H1N1pdm09 | HA | 0.0035 | `0.004` |
| Influenza A | MP | 0.0024 | `0.003` |
| Influenza A | NP | 0.0022 | `0.003` |
| Influenza A | NA | 0.0033 | `0.003` |
| Influenza B | HA | 0.0018 | `0.002` |
| SARS-CoV-2 | Full genome | 0.00074 | `0.0008` |
| SARS-CoV-2 | N gene | 0.00085 | `0.0009` |
| Dengue virus 1–4 | Full genome | 0.00079 | `0.0008` |
| RSV-A | Full genome | 0.00124 | `0.0013` |

*Sources: Bedford et al. (2015), Rambaut et al. (2020), Hadfield et al. (2018), Nextstrain community builds.*
