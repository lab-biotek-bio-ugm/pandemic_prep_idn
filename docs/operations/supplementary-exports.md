# Supplementary Exports Operational Guide

`scripts/export_supplementary.py` is the automated exporter tool for generating the complete supplementary materials package (Executive Summary and detailed S1–S5 supporting workbooks).

---

## 🚀 Basic Usage

```bash
# Auto-detect all target directories in results/ and export:
python scripts/export_supplementary.py \
  --results-dir results \
  --outdir supplementary \
  --only-degenerate
```

---

## 🎯 Target Filtering Options

### 1. Excluding Specific Targets (`--exclude-targets`)
To exclude specific pathogen targets (e.g. excluding Measles when presenting a dual-target Influenza assay package):

```bash
python scripts/export_supplementary.py \
  --results-dir results \
  --outdir supplementary \
  --only-degenerate \
  --exclude-targets measles
```

When `--exclude-targets` is supplied:
- **`00_Executive_Overview`**: Target is omitted from the executive retention table.
- **`01_Synthesis_Oligo_Orders`**: Oligos for the excluded target are omitted from ordering lists.
- **`S1_Primer_Candidates.xlsx`**: Worksheets for the excluded target are skipped.
- **`S2_BLAST_Specificity.xlsx`**: Off-target BLAST data for the excluded target are skipped.
- **`S4_QC_Summary.xlsx`** & **`S5_Conservation_Scoring.xlsx`**: Filtered out completely.

### 2. Config File Target Filtering
You can pass target configuration YAML files directly to restrict outputs to only those active targets:

```bash
python scripts/export_supplementary.py \
  --results-dir results \
  --outdir supplementary \
  --only-degenerate \
  config/config_influenza_a_mp.yaml config/config_influenza_b_mp.yaml
```

---

## 📊 Output Package Structure

| Output Workbook / Folder | Description |
|---|---|
| **`Executive_Summary_Diagnostic_Assays.xlsx`** | Executive Overview, Oligo Synthesis Orders (Rank 1 Best + 3 Alts), Pairwise/N-Plex multiplex screening, and BLAST summary. |
| **`S1_Primer_Candidates.xlsx`** | Complete candidate table (positions, $T_m$, GC%, scores, IUPAC sequences) for each target cohort. |
| **`S2_BLAST_Specificity.xlsx`** | Detailed BLAST off-target specificity logs. |
| **`S3_Multiplex_Compatibility.xlsx`** | Thermodynamic cross-dimer details ($\Delta G$ at 39°C and 55°C). |
| **`S4_QC_Summary.xlsx`** | Sequence acquisition and retention statistics. |
| **`S5_Conservation_Scoring.xlsx`** | Per-target entropy, conservation scores, and grading. |
| **`supplementary/csv/`** | CSV backups of every worksheet for automated downstream processing. |
