---
title: "Developer Guide"
description: "Guide for contributors and developers: repository structure, adding new rules and scripts, coding standards, environment management, and testing protocols."
---

# Developer Guide

This guide is for bioinformaticians and software engineers who want to extend, modify, or contribute to the pipeline codebase.

---

## Repository Structure

```
pandemic-preparedness/
├── Snakefile                    # Main entry point — imports all rules, sets TARGET_GENES
├── config/
│   ├── config_influenza_a_mp.yaml   # Flu A MP segment profile
│   ├── config_influenza_a_np.yaml   # Flu A NP segment profile
│   ├── primerlab_template.yaml      # PrimerLab RAA design template
│   └── references/
│       ├── flu_a_mp_reference.fasta
│       ├── flu_a_mp_annotation.gff
│       └── flu_a_np_reference.fasta
├── rules/
│   ├── m01_acquisition.smk
│   ├── m02_m06_augur_core.smk  # Contains M02–M03 (QC, alignment) + M04–M06 (tree, ancestral)
│   ├── m05_hyphy.smk
│   ├── m07_bridge.smk
│   ├── m08_primerlab.smk
│   ├── m09_blast.smk
│   ├── m10_insilico_pcr.smk
│   └── build_db.smk             # BLAST local database builder
├── scripts/
│   ├── sanitize_metadata.py
│   ├── enrich_metadata.py
│   ├── qc_prefilter.py
│   ├── m02_filter_wrapper.py
│   ├── extract_genes.py
│   ├── hyphy_subsample.py
│   ├── prune_tree.py
│   ├── run_hyphy_gard_safe.py
│   ├── m07_bridge_builder.py
│   ├── m08_payload_generator.py
│   ├── m09_blast_verification.py
│   └── m10_insilico_pcr.py
├── envs/
│   ├── nextstrain.yaml
│   ├── hyphy_macse.yaml
│   ├── primerlab.yaml
│   └── m10_insilico.yaml
├── docs/                        # Mintlify documentation
│   ├── docs.json
│   ├── index.md
│   ├── getting-started/
│   ├── concepts/
│   ├── pipeline-architecture/
│   ├── reference/
│   ├── operations/
│   └── development/
└── results/                     # Generated outputs (gitignored)
```

### Key Design Principle: Configuration-Driven, Zero Hardcoding

The pipeline follows a strict **zero-hardcoding** principle:
- Gene names are parsed dynamically from the GFF3 annotation at workflow load time via `_parse_cds_genes_from_gff()` in `m02_m06_augur_core.smk`.
- `TARGET_GENES` in the Snakefile is populated from the active config + GFF3 parse result.
- Adding support for a new virus means adding a GFF3 annotation file — not modifying any rule or script code.
- All BLAST database taxonomy IDs and host genomes are specified in config, not hardcoded in build scripts.

---

## Conda Environment Management

### Environment Files

| Environment | File | Used By |
|---|---|---|
| `nextstrain` | `envs/nextstrain.yaml` | M01 (Nextclade, NCBI CLI), M02 (augur), M03 (MAFFT), M04–M06 (augur) |
| `hyphy_macse` | `envs/hyphy_macse.yaml` | M05 (HyPhy), M07 (Bridge Builder) |
| `primerlab` | `envs/primerlab.yaml` | M08 (PrimerLab), M09 (BLAST+) |
| `m10_insilico` | `envs/m10_insilico.yaml` | M10 (seqkit) |

### Creating a New Environment

If a new module requires packages not in any existing environment:

1. Create `envs/new_env.yaml` following this template:

```yaml
name: new_env
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.10
  - biopython=1.81
  - pandas>=1.5
  - pip:
    - some-pypi-package>=1.0
```

2. Pin specific versions for all critical tools (never use unpinned `latest`).
3. Test the environment: `conda env create -f envs/new_env.yaml && conda activate new_env`
4. Reference the environment from your rule: `conda: "../envs/new_env.yaml"`

> [!IMPORTANT]
> Do not add packages to `nextstrain.yaml` unless they are required by a nextstrain/augur tool. The nextstrain environment is the largest and most fragile — unnecessary additions risk breaking conda's dependency solver.

---

## Adding a New Snakemake Rule

### Rule Template

All rules must follow this structure:

```python
rule new_rule_name:
    """One-line docstring: what this rule does, why it exists."""
    input:
        required_file=f"results/{BUILD}/mXX/input.fasta",
        config_key=config["some"]["key"]
    output:
        main_output=f"results/{BUILD}/mYY/output.tsv",
        marker=touch(f"results/{BUILD}/mYY/.rule_complete")  # Optional sentinel
    params:
        some_param=config.get("some_section", {}).get("param", default_value)
    threads:
        config.get("section", {}).get("threads", 8)
    conda:
        "../envs/appropriate_env.yaml"
    log:
        f"logs/{BUILD}/mYY_new_rule.log"
    shell:
        """
        set -euo pipefail
        mkdir -p results/{BUILD}/mYY
        
        python scripts/new_script.py \
            --input {input.required_file} \
            --output {output.main_output} \
            --param {params.some_param} \
            2>&1 | tee {log}
        """
```

**Mandatory conventions:**
- `set -euo pipefail` in all shell blocks — causes immediate exit on any error, undefined variable, or pipe failure.
- `mkdir -p` the output directory at the top of every shell block.
- Log file must capture both stdout and stderr: `2>&1 | tee {log}`.
- All parameters must come from `config` — never hardcode values inside a rule.
- Add the rule file import to `Snakefile` if creating a new `.smk` file.

### Adding to an Existing `.smk` File

When extending an existing rule file (e.g., adding a new QC step to `m02_m06_augur_core.smk`), follow the file's existing conventions for variable naming and config access patterns.

---

## Writing a New Python Script

All scripts called by Snakemake rules must follow these conventions:

### Script Template

```python
#!/usr/bin/env python3
"""
Script: scripts/new_script.py
Purpose: Brief description of what this script does.
Called by: Rule new_rule_name in rules/mXY_module.smk

Input:
    --input: Path to input file
Output:
    --output: Path to output file
"""

import argparse
import logging
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr
)
log = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="Path to input file")
    p.add_argument("--output", required=True, help="Path to output file")
    p.add_argument("--param", type=float, default=5.0, help="A parameter")
    return p.parse_args()


def main():
    args = parse_args()
    
    log.info(f"Starting processing: {args.input}")
    # ... logic here ...
    log.info(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
```

**Required conventions:**
- Use `argparse` for CLI — never use `sys.argv` directly.
- Log to `stderr` with `logging` — Snakemake captures stderr to the log file.
- Use `pathlib.Path` for file operations, not raw string concatenation.
- Exit with non-zero status on any unrecoverable error — Snakemake treats non-zero exit as rule failure.
- Include docstring with: purpose, calling rule, inputs, outputs.

### Snakemake `script:` Directive

For scripts that need access to Snakemake's `input`, `output`, `params`, `log`, and `config` objects directly (without `argparse`), use the `script:` directive instead of `shell:`:

```python
# In rule:
#   script:
#       "../scripts/snakemake_script.py"

# In script (note: no if __name__ == "__main__"):
import logging
logging.basicConfig(stream=open(snakemake.log[0], "w"), level=logging.INFO)
log = logging.getLogger()

# Access Snakemake objects directly
input_fasta = snakemake.input.fasta
output_tsv = snakemake.output.tsv
max_n = snakemake.params.max_n_pct

log.info(f"Processing {input_fasta}")
# ... logic ...
```

---

## Dry-Run Testing Protocol

Before pushing any rule or script change, always verify with a dry-run across all supported config profiles:

```bash
# Test against MP config
snakemake \
    --use-conda \
    --configfile config/config_influenza_a_mp.yaml \
    --dryrun --quiet 2>&1 | tail -5

# Test against NP config
snakemake \
    --use-conda \
    --configfile config/config_influenza_a_np.yaml \
    --dryrun --quiet 2>&1 | tail -5

# Validate that no new wildcards are unresolved
snakemake \
    --use-conda \
    --configfile config/config_influenza_a_mp.yaml \
    --dryrun 2>&1 | grep -i "error\|exception\|missing"
```

A clean dry-run shows only job counts with no `Error`, `Exception`, or `Missing` messages.

---

## PrimerLab Integration Development Notes

The pipeline integrates with `primerlab-genomic` as an external dependency. For local development of PrimerLab itself alongside this pipeline:

1. Set `primerlab.repo_path` to an absolute local path:
   ```yaml
   primerlab:
     repo_path: "E:/Project/primerlab-genomic"
   ```
2. M08's `m08_install_primerlab` rule installs PrimerLab in editable mode (`pip install -e`) when a local path is detected.
3. Changes to the local PrimerLab repository are reflected immediately on the next M08 run without re-installation.

**Cross-platform path handling:** The M08 install rule contains path-handling logic for Windows drive letters (e.g., `E:/Project/...`). When running inside WSL2, the rule converts Windows paths to WSL mount paths (`/mnt/e/Project/...`) automatically.

---

## Documentation Conventions

All documentation lives in `docs/` and is rendered by Mintlify. Conventions:

- Page frontmatter requires `title` and `description` fields.
- Use `> [!NOTE]`, `> [!IMPORTANT]`, `> [!WARNING]`, `> [!CAUTION]` for callout blocks.
- Code blocks must include a language identifier (` ```bash`, ` ```yaml`, ` ```python`, ` ```json`).
- All file paths in docs should be relative to the project root.
- Navigation is managed in `docs/docs.json` — add any new page to the appropriate group.

To preview documentation locally:

```bash
cd docs
npx mintlify dev
```
