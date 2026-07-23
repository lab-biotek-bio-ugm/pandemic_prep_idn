---
title: "Testing & Validation Suite"
description: "Guidelines for testing the pipeline: dry-run validation, unit testing local Python scripts, and running integration verification on subset datasets."
---

# Testing & Validation Suite

This guide explains how to execute tests and validate changes in the pipeline codebase. It covers dry-runs, script unit tests, and full integration tests on mock data subsets.

---

## 1. Dry-Run DAG Validation

Before committing any modification to rules (`rules/*.smk`) or configurations, always run a dry-run to verify that Snakemake can successfully compile the Directed Acyclic Graph (DAG) for all supported profiles:

```bash
# Test Influenza A MP profile
snakemake --use-conda --configfile config/config_influenza_a_mp.yaml -n -q

# Test Influenza A NP profile
snakemake --use-conda --configfile config/config_influenza_a_np.yaml -n -q
```

A successful dry-run returns a job count table with zero warnings or errors. Common failure causes include syntax errors in rule files, duplicate rule definitions, or missing input declarations.

---

## 2. Unit Testing Python Scripts

### Automated pytest suite

A `pytest` suite in `tests/` covers the correctness-critical scientific functions
(entropy scoring, RAA primer binding, dye-channel separation, QC column detection).
These are regression guards for previously-fixed bugs — run them before every commit:

```bash
pytest -q
```

The suite and a Snakemake dry-run also run automatically in CI
(`.github/workflows/tests.yml`) on every push/PR that touches `scripts/`, `rules/`,
`tests/`, or the `Snakefile`. The dry-run step fails on rule-definition errors
(e.g. a wildcard that appears only in `input`), catching DAG breakage without data.

Current coverage:

| Test file | Guards |
|---|---|
| `tests/test_m07_entropy.py` | Gap-inclusive Shannon entropy; deletion hotspots score low, not "perfectly conserved" |
| `tests/test_m10_raa_binding.py` | Reverse-primer 3'-terminus orientation in RAA mismatch tolerance |
| `tests/test_m11_dye_and_utils.py` | Spectral dye-channel collision detection; export active-target filter |
| `tests/test_qc_prefilter.py` | M02 subtype-scarcity filter runs on the post-sanitize `strain`/`subtype` schema |

### Running scripts directly

The scripts in `scripts/` perform data sanitisation, filtering, scoring, and output generation. They are designed as standard modular Python files that can be run independently of Snakemake.

### Setting up a Test Environment
Activate the environment containing the dependencies of the script you want to test (e.g., `hyphy_macse` for Bridge scripts, `primerlab` for design verification).

```bash
# Activate appropriate conda environment
conda activate hyphy_macse
```

### Running Test Scripts
We provide mock input files in `data/test_subset/` to allow developers to verify script logic quickly:

```bash
# Test the metadata sanitizer
python scripts/sanitize_metadata.py \
    --input dataset/EpiFlu/MP/EpiFlu_A_MP_Combined_Metadata.csv \
    --output results/test_debug/sanitized_metadata.tsv \
    --format gisaid_epiflu

# Test the Bridge Builder logic
python scripts/m07_bridge_builder.py \
    --msa results/influenza-a-mp-sea/m03/aligned.fasta \
    --gff config/references/flu_a_mp_annotation.gff \
    --hyphy-dir results/influenza-a-mp-sea/m05 \
    --reference config/references/flu_a_mp_reference.fasta \
    --metadata results/influenza-a-mp-sea/m01/metadata_enriched.tsv \
    --outdir results/test_debug \
    --config results/influenza-a-mp-sea/m07/active_config.yaml \
    --build-name test-run \
    --strategy full_segment
```

Verify that the output files are generated without exceptions, and check that the resulting candidate logs match expected values.

---

## 3. Integration Testing with Mock Datasets

We recommend maintaining a small "mock dataset" (e.g., 50 sequences) to run full end-to-end integration tests in a fraction of the time required for complete surveillance datasets.

### Creating a Mock Run Profile
Create a configuration profile specifically for integration testing: `config/config_test_mock.yaml`.

```yaml
# config/config_test_mock.yaml
build_name: "test-mock-run"
database: "local"
references:
  reference_id: "NC_002016.1"
  genome: "config/references/flu_a_mp_reference.fasta"
  annotation: "config/references/flu_a_mp_annotation.gff"
local_input:
  sequences: "data/test_subset/mock_sequences.fasta"
  metadata: "data/test_subset/mock_metadata.csv"
  format: "tsv"
qc:
  max_n_content: 10.0
subsampling:
  max_sequences: 50
  group_by: "country"
  sequences_per_group: 50
hyphy:
  methods: ["FUBAR"]
  subsample:
    enabled: false  # Avoid down-sampling the mock data further
bridge:
  strategy: "full_segment"
  entropy_threshold: 0.20 # Relaxed for small mock dataset
```

### Executing the Integration Test
Run the full Snakemake workflow using the mock profile:

```bash
snakemake \
    --use-conda \
    --configfile config/config_test_mock.yaml \
    --cores 4 \
    --reason \
    2>&1 | tee logs/mock_test_run.log
```

**What to verify after completion:**
- All outputs are successfully written to `results/test-mock-run/`.
- No Snakemake rules exited with non-zero status.
- Inspect `results/test-mock-run/m10/coverage_summary.tsv` to ensure the final diagnostic candidates were validated and received non-zero coverage scores.
