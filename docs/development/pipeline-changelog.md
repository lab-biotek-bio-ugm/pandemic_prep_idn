# Changelog

All notable changes to the Evolutionary Diagnostic Pipeline will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<Note>
**Semantic Versioning & Git Commits**: Since this repository does not have formal tags yet, this changelog acts as a detailed development record. Each minor version group contains **dated patch logs** mapping exact git commits and dates to their functional updates.
</Note>

---

## [1.1.0] - 2026-07-23 (Correctness Audit: Literature Validation of Design Rules)

A full audit of the pipeline's scripts and algorithms, followed by validation of the assay-design rules against primary sources. Several rules turned out to be unsupported or inverted; those are corrected here with the evidence recorded, not silently changed.

### Corrected — assay design rules

- **Removed the "2 bp THF cleavage window" rule (M10).** The documentation stated that mismatches within 2 bp of the abasic site were prohibited and caused false negatives. No supporting source exists, and the TwistAmp Assay Design Manual contradicts it directly: §2.7 reports RPA performance unaffected by 5–9 mismatches in primers and probe, and §3.1.1 sanctions a *deliberate* mismatch on a dT-label position — which §3.1.2 places 0–2 bases from the abasic residue. A supporting-looking Exonuclease III result ([NAR 24(22):4572](https://academic.oup.com/nar/article/24/22/4572/2386348)) was checked and does not apply: those mismatches sat **at** the abasic position, not near it.
- **Fixed abasic-site placement in PrimerLab (M08).** The residue could land 28–29 bases from the probe's 5' end, which the manual lists verbatim as a "Poor probe design", and the 30–38 working range had no upper bound. More seriously, the label logic was **inverted**: the code anchored the abasic residue onto the nearest `T` and let the dT-fluorophore and dT-quencher fall on whatever bases were adjacent, so every probe could silently carry two unnecessary probe–target mismatches. The flanking positions are what must be `T`; the abasic position has no sequence requirement. Placement now degrades in the order the manual prescribes and reports every relaxation.
- **Verified against the manual's worked example**, which the corrected implementation now reproduces exactly (30 bases 5', 15 bases 3', both labels on real thymines).

### Changed — multiplex screening (M11)

- **Unified all thresholds and the score formula** into `scripts/tools/iupac_utils.py`. Pairwise and triplex previously held separate literals that had already drifted: $\Delta T_m$ was 3.5°C in one and 4.0°C in the other, so an identical oligo set could be `PASS` in one report and `WARNING` in the other.
- **Unified the score weights** to $\Delta G \times 4.0$, $\Delta T_m \times 3.0$ (pairwise previously used 3.5/4.0, weighting them the opposite way). $\Delta G$ dominates because at 39°C strand invasion is recombinase-driven, so inter-primer $T_m$ matching matters far less than dimer formation.
- **Stopped rendering the compatibility score as a percentage.** It is an uncalibrated 0–100 heuristic, not a probability.
- **Added a 3'-end anchoring check.** Only an annealed 3' terminus can be extended, and those products amplify as primer noise — the factor TwistAmp §2.6 names as limiting sensitivity. The criterion is the longest 3'-anchored complementary run (> 5 bp warns), *not* a $\Delta G$: a $\Delta G$ over a short terminal window was implemented first and rejected on measurement, because a 5 nt window saturates near −6.3 kcal/mol regardless of how much complementarity exists and an 8 nt window overlaps the clean background.
- **Documented the $\Delta G$ thresholds' provenance.** No validated RPA/RAA-specific threshold exists — TwistAmp §2.2 says there are "no fixed rules" and prescribes bench screening. The −9/−12 kcal/mol values derive from general IDT guidance, and because they are evaluated at 39°C rather than 60°C they are stricter than they appear: over 200 random 32-mer pairs, 88% fall below −3 at 39°C versus 6.5% at 60°C.
- **Every report now carries a methodology note** stating these caveats, so a bench reader sees them without opening the docs.

### Added — PrimerLab data flow

- `Primer.to_dict()` now exports `thf_index`, `probe_type`, `probe_label_mismatches` and `probe_compliant`. The abasic-site position was computed and then discarded before serialisation, so no downstream module could locate the cleavage site.

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
