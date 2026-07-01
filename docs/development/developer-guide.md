# Developer & Contributor Guide

This guide details the coding standards, workflow conventions, and architectural patterns required to maintain and extend the pipeline.

---

## 1. Extending the Pipeline (Adding a New Module)

To add a new processing step (e.g., Module M11) to the workflow:

1.  **Define Rules**: Create a new rule file in `rules/m11_new_feature.smk`.
2.  **Declare Conda Env**: Create a matching conda environment file in `envs/m11_env.yaml` with specific package pins.
3.  **Write Rule Logic**:
    *   Always use relative paths configured via `config.yaml`.
    *   Explicitly define `conda:` and `log:` parameters.
    *   Do not write processing logic directly in the rule; call external scripts under `scripts/`.
4.  **Include Rule**: Include the rule file in the main `Snakefile` at the root directory:
    ```python
    include: "rules/m11_new_feature.smk"
    ```
5.  **Declare Targets**: Add the final output files of your rule to the `all` rule in the `Snakefile` to ensure it executes.

---

## 2. Snakemake Rule Design Standards

All Snakemake rules must adhere to the following template:

```python
rule m11_example_step:
    input:
        sequences = "results/{build_name}/m03/aligned.fasta",
        metadata = "results/{build_name}/m01/metadata.tsv"
    output:
        processed = "results/{build_name}/m11/processed.tsv"
    log:
        "results/{build_name}/logs/m11_example_step.log"
    conda:
        "../envs/m11_env.yaml"
    params:
        threshold = config["new_feature"]["threshold"]
    shell:
        "python scripts/run_feature.py "
        "--input {input.sequences} "
        "--metadata {input.metadata} "
        "--threshold {params.threshold} "
        "--output {output.processed} > {log} 2>&1"
```

---

## 3. Python Coding Standards

To ensure readability and maintenance:

*   **Version**: Use Python 3.10+ features (e.g., structural pattern matching, union type hinting `int | None`).
*   **Style**: Adhere strictly to **PEP8** (4-space indentation, descriptive lowercase variable names with underscores).
*   **Logging**: Avoid using raw `print()` statements in workflow scripts. Use the built-in `logging` module to output standard levels:
    ```python
    import logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    logging.info("Starting analysis...")
    ```
*   **Error Handling**: Implement defensive try-except blocks. If a script fails, it must write a clear error description to the log and exit with a non-zero code (`sys.exit(1)`) to notify Snakemake of the failure.

---

## 4. Git & Commit Workflow

We enforce **Conventional Commits** to maintain a clean git history. Commit messages must use the following prefixes:

*   `feat(scope): ...` — A new rule, script, or feature.
*   `fix(scope): ...` — A bug fix in rules or scripts.
*   `docs(scope): ...` — Updates to documentation or comments.
*   `refactor(scope): ...` — Code reorganizations without logic changes.
*   `chore(scope): ...` — Dependency updates or build system changes.

### Pull Requests
Before opening a pull request, run a full dry-run on both H5N1 and SARS-CoV-2 config profiles to verify that the DAG constructs cleanly.
