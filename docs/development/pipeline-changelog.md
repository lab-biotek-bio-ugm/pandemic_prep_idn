# Changelog

All notable changes to the Evolutionary Diagnostic Pipeline will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<Note>
**Semantic Milestone Versioning**: Since this repository does not yet have official GitHub release tags, the versions below document **actual development milestones** reconstructed directly from the Git commit history.
</Note>

---

## [1.1.1] - 2026-07-02

This patch release fixes build compilation errors on the Mintlify documentation platform, improves contrast styling in light/dark themes, and resolves file-naming conflicts.

### Fixed
- **MDX Parser Exceptions**: Escaped all raw less-than (`<`) signs outside code spans across all `.md` files to `&lt;` to prevent the MDX parser from throwing JSX compilation errors (e.g. `Unexpected character 2 (U+0032) before name`).
- **LaTeX Math Support**: Converted standard math inequality operators in equations to LaTeX-safe expressions (e.g. `\lt` instead of `<`) to permanently bypass MDX parser conflicts in thermodynamic design equations.
- **Adaptive Flowchart Colors**: Removed hardcoded inline CSS fills (`fill:#e3f2fd`, etc.) from the Mermaid flowcharts in `index.md` and `overview.md` to allow theme-switching between light and dark modes to render high-contrast outlines and text natively.
- **Changelog 404 Resolution**: Renamed `changelog.md` to `pipeline-changelog.md` and updated the navigation index to bypass Mintlify's default auto-ignore rules for `changelog.md` files.

---

## [1.1.0] - 2026-07-01

This minor release introduces configuration profiles for new viral segments, parameterization of computing resource allocations, and stability fixes for evolutionary models.

### Added
- **Multi-Pathogen Configurations**: Added curated reference genomes, annotations, and default parameters for:
  - **Influenza A NP Segment** (`config_influenza_a_np.yaml` / reference `NC_002019.1`)
  - **Influenza A HA Segment** (`CY121680.1` and `MK626062.1` templates)
  - **SARS-CoV-2 Whole-Genome / Segmented ORFs** (`NC_045512.2` template)
- **HPC Parallelism & Threading Control**:
  - Parameterized CPU thread limits for `augur align`, `augur tree` (`augur.threads`), and HyPhy rules (`hyphy.threads`).
  - Added SLURM job submission script templates (`sbatch slurm-run.sh`) to support background job execution on compute clusters.
- **GARD Pre-Flight Validation**: Implemented `run_hyphy_gard_safe.py` wrapper to automatically bypass GARD and soft-fail (writing empty outputs and exiting code 0) if the sequence alignment length is too short to satisfy the c-AIC model criteria (`alignment_length < 4 * num_sequences`).

### Changed
- **Modular Documentation Migration**: Dismantled the bundled module documentation (`m01-m03-qc-alignment.md` and `m04-m06-phylodynamics.md`) and refactored them into individual, technically deep documentation pages.
- **M07 Scoring Optimization**: Refactored scoring weights, implemented Non-Maximum Suppression (NMS) for overlapping candidate regions, and added support for a sliding 150bp probe search window.
- **M02 Subsampling Refinement**: Standardized subsampling to use `--subsample-max-sequences` (mutually exclusive with `--sequences-per-group`) to prevent Augur execution errors.

---

## [1.0.0] - 2026-06-30

This major release marks the first stable, end-to-end version of the 10-Module Evolutionary Diagnostic Pipeline.

### Added
- **M01 Data Ingestion**: Added a metadata sanitizer with 3-tier validation and auto-format detection (supporting GISAID EpiFlu `.xls` spreadsheets, local GenBank CSV files, and NCBI Datasets API).
- **M02 Quality Control**: Implemented sequence deduplication (MD5 hashing), N-content filtering, and early query filtering.
- **M03 Alignment**: Integrated MAFFT reference-guided alignment with `--fill-gaps` to establish a uniform coordinate system.
- **M04 Phylogenetics**: Configured maximum-likelihood tree construction (IQ-TREE2) and temporal tree reconstruction (TreeTime).
- **M05 Selection Analysis**: Integrated HyPhy selection engine (FUBAR for pervasive selection, MEME for episodic hotspots, and GARD for recombination boundaries) with spliced-gene support.
- **M06 Ancestral Mapping**: Configured ancestral state reconstruction and Auspice JSON exports for Nextstrain visualization.
- **M07 Bridge Builder**: Integrated entropy scores (Track A) and selection pressure posterior probabilities (Track B) into a unified target ranking manifest.
- **M08 PrimerLab Integration**: Configured thermodynamic primer and probe design with degeneracy constraints (max 16x multiplier) and secondary structure validations (ViennaRNA).
- **M09 Specificity Check**: Built custom respiratory database generator (`build_db.smk`) and local BLAST specificity checks (including critical 3' end mismatch weighting).
- **M10 In-Silico PCR**: Simulated diagnostic coverage using `seqkit amplicon` with custom RAA 3' kinetic extension and probe THF cleavage validation.
