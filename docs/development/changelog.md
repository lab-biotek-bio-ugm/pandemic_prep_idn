# Changelog

All notable changes to this project will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.1.0] - 2026-07-02

### Added
- **Separate Module Documentation**: Split the old combined module documentation into individual pages: `m01-acquisition.md`, `m02-qc-subsampling.md`, `m03-alignment.md`, `m04-phylogenetics.md`, `m05-selection-analysis.md`, `m06-ancestral-mapping.md`.
- **Shannon Entropy Concept Page**: Math, formulas, step-by-step calculations, and gap handling logic.
- **Selection Pressure Concept Page**: Comparison between FUBAR, MEME, and GARD with molecular biology rationale.
- **Data Preparation Guide**: Detailed GISAID, NCBI Virus, GenBank download processes and schemas.
- **Output Files Catalog**: Complete schemas for `metadata_enriched.tsv`, `fubar_{gene}.json`, `meme_{gene}.json`, and more.
- **Supported Pathogens**: Configuration profiles, clock rates, reference IDs, and parameter defaults for Flu A MP, NP, HA, and SARS-CoV-2.
- **Monitoring Operations Page**: CPU/HPC/Slurm configs, checkpoint behaviors, and disk estimates.
- **Data Management & GISAID Compliance**: Guidelines for database compliance, DAA terms, and workspace conventions.
- **Developer & Testing Guides**: Complete repository layouts, rule/script coding standards, and unit/integration testing suites.

### Changed
- **Navigation Structure (`docs.json`)**: Restructured Mintlify navigation to separate Getting Started, Concepts, Pipeline Architecture, Reference, Operations, Development.
- **Overview.md Rewrite**: Added full module Input/Output reference table, environment maps, and deferred features tracker.
- **Configuration.md Rewrite**: Expansion on molecular clock rate priors, subsampling limits, and degeneracy multipliers.

---

## [1.0.0] - 2026-06-30

### Added
- **Initial Snakemake Diagnostic Pipeline**: Complete 10-module workflow orchestration.
- **Augur Integration**: Subsampling (`augur filter`), alignment (`mafft`), tree building (`iqtree`), ancestral reconstruction (`augur ancestral/translate`).
- **Selection Pressure Engine**: Local integration with HyPhy FUBAR, MEME, and GARD.
- **Diagnostic Bridge Builder (M07)**: Merging physical conservation and evolutionary constraint into a composite score.
- **PrimerLab (M08)**: Degenerate primer design, ViennaRNA hairpin and dimer thermodynamic constraints.
- **BLAST Specificity (M09)**: Cross-reactivity screening (human + target genomes) using remote NCBI nt or custom local indexed databases.
- **In-silico PCR (M10)**: seqkit amplicon-based validation with custom RAA 3' kinetic extension and probe THF cleavage windows.
- **Interactive HTML Dashboard**: Plotly-based HTML heatmaps displaying temporal × subtype diagnostic sensitivity.
- **Conda Environments**: Isolated bioconda/conda-forge environments for every step.
- **Initial Mintlify Documentation**: Overview, quickstart, troubleshooting, and developer guide stubs.
