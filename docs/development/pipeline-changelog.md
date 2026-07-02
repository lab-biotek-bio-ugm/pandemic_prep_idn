# Changelog

All notable changes to the Evolutionary Diagnostic Pipeline will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<Note>
**Semantic Versioning & Git Commits**: Since this repository does not have formal tags yet, this changelog acts as a detailed development record. Each minor version group contains **dated patch logs** mapping exact git commits and dates to their functional updates.
</Note>

---

## [1.0.0] - 2026-07-02 (Official Stable Production Release)

The first stable production release of the 10-Module Evolutionary Diagnostic Pipeline.

### v1.0.0-patch.1 - 2026-07-02 (Final Documentation & Parser Polish)
- **Fixed**: MDX parser syntax exceptions; escaped all raw `<` signs in markdown files to `&lt;` and equations to `\lt` to prevent Mintlify JSX compilation errors.
- **Fixed**: Theme outlines for Mermaid flowcharts; removed hardcoded styling to enable dynamic high-contrast adaptiveness in light/dark themes.
- **Changed**: Renamed changelog file to `pipeline-changelog.md` to bypass Mintlify's auto-ignore filter.

### v1.0.0-rc.1 - 2026-07-02 (Modular Documentation Upgrade)
- **Added**: Complete technical manuals and pages for M01-M10 (split from old combined guides).
- **Added**: Scientific concept pages for Shannon Entropy calculations and Selection Pressure models (FUBAR vs MEME vs GARD).
- **Added**: Reference sections detailing configuration parameters, expected output files, and default pathogen attributes (Flu A, SARS-CoV-2).
- **Added**: Operations and Development guides covering WSL2 setup, GISAID DAA compliance, Slurm HPC monitoring, and integration testing frameworks.

---

## [0.9.0] - 2026-06-30 (Quarto Report Compilation & Mintlify Ingestion)

This minor version compiles the final Quarto scientific reports, sets up automatic pdf compilation scripts, and migrates user-facing manuals to the Mintlify tab-based web theme.

### v0.9.0-patch.2 - 2026-07-02
- **Changed**: Migrated documentation to Mintlify; setup UGM Biology favicon, Maple theme styling, and horizontal Mermaid diagram flows.
- **Fixed**: Cleaned up leading slashes in `docs.json` paths and resolved index page hydration crashes.

### v0.9.0-patch.1 - 2026-07-01
- **Added**: Supplementary materials section discussing degeneracy multiplier limitations and raw IUPAC TSV downloads.
- **Added**: Automatic PDF compiler script converting all 20 pages of the Quarto report into a single printable handbook.
- **Fixed**: Pandoc citation rendering and list indentations in Quarto results chapters.

### v0.9.0-rc.1 - 2026-06-30
- **Added**: Quarto diagnostic development report site with live embedded IGV.js track views.
- **Changed**: Documented real phylogenetic tree composition counts, clock deviations, and FUBAR purifying score plots for HA/MP segments.

---

## [0.8.0] - 2026-06-24 (In-Silico PCR Verification & HTML Heatmap Dashboard)

This minor version completes Module 10 (In-Silico Validation) with parallelized amplification checking and dark-mode temporal × regional coverage heatmaps.

### v0.8.0-patch.2 - 2026-06-29
- **Added**: Curated HA and NP reference genome templates.
- **Added**: Toggle support for host genome exclusions in the local BLAST database builder.

### v0.8.0-patch.1 - 2026-06-28
- **Added**: Dual-coverage statistics panel on the in-silico PCR dashboard (comparing design database vs surveillance database coverage).
- **Fixed**: Coverage heatmap orientation bug causing ~19% false-negative coverage reports.

### v0.8.0-rc.1 - 2026-06-24
- **Added**: `m10_insilico_pcr.py` with multi-core ProcessPoolExecutor execution.
- **Added**: Plotly-based premium dark-mode HTML coverage dashboard displaying temporal and regional sensitivity.

---

## [0.7.0] - 2026-06-11 (PrimerLab thermodynamic integration & BLAST M08-M09)

This minor version integrates PrimerLab RT-RAA degenerate design and builds custom local BLAST databases for specificity screening.

### v0.7.0-patch.2 - 2026-06-18
- **Fixed**: Entropy-based coordinate exclusion rules (IUPAC threshold 15%) to prevent tiny wobble-position blocks from obstructing primer design.
- **Changed**: Switched to a full-segment design strategy (always use full segment without UTR coordinate restriction).

### v0.7.0-patch.1 - 2026-06-14
- **Added**: ViennaRNA secondary structure checks (hairpin, homodimer, heterodimer $\Delta G$ filters) at reaction temperatures.
- **Fixed**: Shell-unsafe single quotes in UTR coordinate names (3'UTR, 5'UTR).

### v0.7.0-rc.1 - 2026-06-11
- **Added**: Custom BLAST database builder (`build_db.smk`) compiling human host and pathogen GenBank sequences.
- **Added**: M08 degenerate primer payload generator and automated CLI executor.

---

## [0.6.0] - 2026-06-02 (Positional Scoring & IGV Session Integrations)

This minor version implements the first version of the Bridge Builder (M07) scoring algorithm and exports XML session configurations for local IGV rendering.

### v0.6.0-patch.1 - 2026-06-08
- **Added**: Positional color-coded tracks for IGV (separating entropy, purifying selection, and recombination bounds).
- **Fixed**: KeyErrors during selection mapping for genes lacking selection pressure records.

### v0.6.0-rc.1 - 2026-06-02
- **Added**: M07 Bridge Builder composite score formulas ($w_c \cdot C + w_p \cdot P$).
- **Added**: Holistic single-region consensus generation mode for whole-segment strategies.

---

## [0.5.0] - 2026-05-05 (Pathogen Specificity Support & Nextclade Enhancements)

This minor version calibrates the pipeline for multiple pathogen segments (Influenza A HA, NP, MP) and implements Nextclade bypass mode.

### v0.5.0-patch.2 - 2026-05-21
- **Added**: Curated HA reference FASTA templates and Nextclade bypass switches.
- **Added**: Metadata sanitizer supporting GISAID EpiFlu `.xls` spreadsheets.

### v0.5.0-patch.1 - 2026-05-12
- **Added**: Automatic GFF translation feature extracting CDS gene names dynamically.
- **Fixed**: Spliced M2 gene parsing and translation alignments.

### v0.5.0-rc.1 - 2026-05-05
- **Added**: Pathogen Taxon ID mappings for custom database building.
- **Added**: Early query filter to handle mixed-subtype FASTA datasets.

---

## [0.4.0] - 2026-04-10 (HPC Parallelism & GARD Safe Wrapper)

This minor version prepares the pipeline for cluster execution, pins dependencies, and creates safe wrappers for recombination checks.

### v0.4.0-patch.2 - 2026-04-13
- **Added**: Thread count parameterization for `augur tree` and HyPhy rules.
- **Changed**: Disabled GARD for Influenza MP configs to optimize execution runtime.

### v0.4.0-patch.1 - 2026-04-11
- **Added**: Pre-flight checks wrapper `run_hyphy_gard_safe.py` bypassing c-AIC violations.
- **Fixed**: Conda package pins for `iqtree` 2.x and `augur`.

### v0.4.0-rc.1 - 2026-04-10
- **Added**: SLURM scheduler integration.
- **Added**: Multi-core configuration parameters in default profiles.

---

## [0.3.0] - 2026-04-06 (Selection Pressure & Alignment Sanitization)

This minor version integrates codon-level selection modeling (FUBAR/MEME/GARD) and alignment cleaning.

### v0.3.0-patch.2 - 2026-04-07
- **Fixed**: Reading-frame masking of stop codons and frame-shift insertions.
- **Fixed**: Duplicated pandas metadata columns (`_x`/`_y`).

### v0.3.0-patch.1 - 2026-04-06
- **Added**: Direct slicing of codon regions from GFF annotations.
- **Fixed**: Gap-sequence column filtering.

### v0.3.0-rc.1 - 2026-04-06
- **Added**: Automated codon alignment using MACSE.
- **Added**: Local HyPhy selection engine calls (FUBAR, MEME, GARD).

---

## [0.2.0] - 2026-04-02 (Core Augur Orchestration - M01-M06)

This minor version builds the basic Nextstrain Augur orchestration rules running through temporal tree refinement.

### v0.2.0-patch.1 - 2026-04-02
- **Added**: Nextclade clade-assignment rule injecting metadata fields.
- **Fixed**: Ambiguous date filtering causing silent sequence drops.

### v0.2.0-rc.1 - 2026-04-02
- **Added**: Rules for M01 sequence/metadata acquisition.
- **Added**: Rules for `augur filter`, `augur align` (MAFFT), `augur tree` (IQ-TREE2), `augur ancestral`, and `augur export`.

---

## [0.1.0] - 2026-03-18 (Inception & Master Plan)

The inception milestone establishing requirements and diagram flows.

### v0.1.0-patch.1 - 2026-04-01
- **Added**: Initial repository structures, `.gitignore`, and environment specifications.

### v0.1.0-rc.1 - 2026-03-18
- **Added**: Brainstorming document and master diagram flow designs.
