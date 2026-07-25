---
title: "Testing & Validation Suite"
description: "Guidelines for testing the pipeline: dry-run validation, unit testing local Python scripts, and running integration verification on subset datasets."
---

# Testing & Validation Suite

This guide explains how to execute tests and validate changes in the pipeline codebase. It covers dry-runs, script unit tests, and full integration tests on mock data subsets.

---

## 0. The Local Gate

Everything CI runs, you can run in one line. Install the tooling once:

```bash
pip install -r requirements/dev.txt
pre-commit install          # optional: runs the same checks on git commit
```

Then, before every commit:

```bash
ruff check .            # undefined names, bare except, mutable defaults
pytest                  # the whole suite, about four seconds
```

`pre-commit install` wires the lint and the config contracts into `git commit`.
Every hook is check-only; none of them rewrites a file. That rule exists because
a formatter would renormalise line endings and bury the real change in a
whole-file diff, which is also why `.gitattributes` pins `* -text`.

---

## 1. Dry-Run DAG Validation

Before committing any modification to rules (`rules/*.smk`) or configurations, verify that Snakemake can compile the Directed Acyclic Graph for **every** pathogen profile, not just the one you are working on. This is the loop CI runs:

```bash
for cfg in config/config*.yaml; do
    echo "== $cfg"
    snakemake --dry-run --quiet --configfile "$cfg"
done
```

A profile passes when the output reaches `Building DAG of jobs`. A `MissingInputException` after that point is expected on a fresh clone, because no raw sequence data is committed; it means the rules are sound and only the cohort is absent. Anything else — `WildcardError`, `RuleException`, `KeyError`, an ambiguous rule — is a real defect.

One thing a dry-run cannot tell you: `Snakefile` declares `configfile: "config/config.yaml"`, and Snakemake **deep-merges** any `--configfile` on top of that rather than replacing it. A per-pathogen config that omits a key therefore does not fail; it silently inherits the influenza-A-MP value. Section 2 is what catches those.

---

## 2. Unit Testing Python Scripts

### Automated pytest suite

The `pytest` suite in `tests/` covers the correctness-critical scientific functions and, in `tests/contracts/`, the agreement between files. Run it before every commit:

```bash
pytest                      # everything
pytest tests/contracts/     # config and CI contracts only, sub-second
```

Current coverage:

| Test file | Guards |
|---|---|
| `tests/contracts/test_config_contract.py` | Required keys, types and bounds in every pathogen config; thresholds that must not drift between them; keys that no rule or script reads; duplicate keys; numbers folded into strings |
| `tests/contracts/test_ci_contract.py` | Every workflow job carries the dormancy guard; no YAML anchors; actions pinned |
| `tests/test_regression_fixes.py` | The 1.1.0 correctness audit: triplex export leak, dye emptiness, M08 sentinels, entropy normaliser, M11 threshold unification, M07 core-score model, M09 verdict labelling |
| `tests/test_candidate_selection.py` | Multiplex candidate selection and gating |
| `tests/test_scripts_execute.py` | M07 end-to-end as a subprocess against a synthetic fixture, including the negative path; benchmark instrumentation |
| `tests/test_multiplex_orchestrator.py` | Batch orchestrator and its bash wrapper; real dye configs are spectrally separable |
| `tests/test_m07_entropy.py` | Gap-inclusive Shannon entropy; deletion hotspots score low, not "perfectly conserved" |
| `tests/test_m10_raa_binding.py` | Reverse-primer 3'-terminus orientation in RAA mismatch tolerance |
| `tests/test_m11_dye_and_utils.py` | Spectral dye-channel collision detection; export active-target filter |
| `tests/test_proc_utils.py` | Subprocess helpers and binary lookup |
| `tests/test_qc_prefilter.py` | M02 subtype-scarcity filter runs on the post-sanitize `strain`/`subtype` schema |

### The config contract

`tests/contracts/config_contract.yaml` is the single declaration of what a
pathogen config must contain: required keys with their types and bounds,
thresholds that are pipeline-wide policy and must be identical everywhere, and
keys that are read by nothing and must not appear. Each file is loaded on its
own, without the Snakemake merge, which is the only way the holes are visible.

Violations that exist today are recorded under `known_gaps` rather than hidden.
The tests assert the violation set matches those lists **exactly**, so a new
violation still fails, and repairing a recorded one also fails until its entry
is deleted. The lists can only shrink.

### Continuous integration

`.github/workflows/ci.yml` runs three jobs on every push and pull request that
touches `scripts/`, `rules/`, `tests/`, `config/`, `envs/`, `docs/`,
`requirements/` or the `Snakefile`:

| Job | What it does | Cost |
|---|---|---|
| Contracts and lint | `ruff check` plus `pytest tests/contracts/`, on two dependencies | about a minute |
| Test suite | the full suite against pinned scientific libraries | two to three minutes |
| DAG validation | byte-compiles `scripts/`, then dry-runs every `config/config*.yaml` | about three minutes |

Both workflows are gated on `github.repository`, so they stay dormant on the
upstream lab repository, whose Actions quota they would otherwise spend.

The contract job selects tests by **path**, never by `-m contract`. Marker
deselection happens after collection, and collecting `tests/` imports pipeline
scripts that need pandas, Biopython and primer3 — one of them calls `sys.exit(1)`
on `ImportError`, which crashes collection outright. Path selection is what keeps
that job down to two dependencies.

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

This section is a **recipe, not a description of the repository**. Neither
`config/config_test_mock.yaml` nor `data/test_subset/` is shipped, and `data/`
is gitignored, so a mock run has to be assembled locally before it can be used.
CI does not run one: a real pipeline pass takes hours and needs conda
environments, HyPhy, IQ-TREE and a multi-gigabyte BLAST database.

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
