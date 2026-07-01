# Troubleshooting Guide

This guide list common issues encountered while running or deploying the pipeline, along with step-by-step resolution instructions.

---

## 1. Conda & Environment Errors

### Error: CondaHTTPError / Connection Timeouts
Occurs during the environment build phase when Snakemake downloads bioconda dependencies.
*   **Cause**: Heavy traffic on NCBI/Anaconda servers, or proxy/firewall blockages.
*   **Resolution**: 
    1. Force conda/mamba to use clean channels and increase timeout limits:
       ```bash
       conda config --set conda_build:max_log_lines 50000
       conda config --set remote_read_timeout_secs 600.0
       ```
    2. Try running the build again. If it continues to fail, build the environment manually using `mamba env create -f envs/m01_acquisition.yaml`.

### Error: "vienna-rna" or "hyphy" installation failed
*   **Cause**: Mismatch in system architecture (e.g., M1/M2 Mac ARM vs. Linux x86_64).
*   **Resolution**: Ensure you are running the pipeline inside a Linux-64 environment (such as WSL2 on Windows) where all bioconda binary builds are fully supported.

---

## 2. Nextstrain / Augur Refine Errors

### Error: `augur refine` fails with negative branch lengths
*   **Cause**: The molecular clock rate specified (`clock_rate`) is incompatible with the temporal divergence of the input sequences, or some sequence dates are incorrect in the metadata.
*   **Resolution**:
    1. Check for metadata format errors. Look for sequences with future dates or invalid formats (e.g., `2026-02-30` or `2026-XX-XX`).
    2. Slightly increase or decrease the `clock_rate` prior (e.g., adjust from `0.003` to `0.002` or `0.004` based on the pathogen) to allow the clock fitting algorithm to converge.

### Warning: "Polytomies resolved stochastically"
*   **Cause**: Many sequences in the dataset are identical or highly similar, resulting in nodes with multiple branches.
*   **Resolution**: This is a normal warning and can be safely ignored. If you wish to reduce polytomies, increase the subsampling filters in `config.yaml` to lower the sequence count.

---

## 3. HyPhy selection engine issues

### Error: FUBAR or MEME crashes with stop codon errors
*   **Cause**: The nucleotide alignment contains internal stop codons or frame-shift mutations, which are invalid in a codon-aware alignment.
*   **Resolution**:
    1. The pipeline automatically runs a sanitization script (`scripts/sanitize_alignment.py`) to strip stop codons.
    2. If the crash persists, check the input sequences for sequencing errors (insertions/deletions) that introduce frame-shifts. Discard the affected sequences from the raw database.

---

## 4. Local BLAST Database Errors

### Error: "Database not found" during M09 specificity check
*   **Cause**: The local database directory path configured in `config.yaml` is invalid, or the BLAST database files are missing.
*   **Resolution**:
    1. Verify that database files (`.nhr`, `.nin`, `.nsq`) exist under the target directory.
    2. Re-create the database using the NCBI command line utility:
       ```bash
       makeblastdb -in path/to/reference.fasta -dbtype nucl -out data/blast_db/pathogen_db
       ```
    3. Update the `local_db_path` parameter in `config.yaml` to point to the exact database prefix (without file extensions).

---

## 5. Out of Memory (OOM) Errors

### Error: Snakemake rules killed by system (OOM)
*   **Cause**: Building maximum likelihood trees (IQ-TREE) or running HyPhy FUBAR/MEME on datasets with $>3000$ sequences exceeds system RAM (often on 8GB or 16GB workstations).
*   **Resolution**:
    1. Decrease the `max_sequences` parameter in `config.yaml` to `1000` or `1500`.
    2. Limit the number of CPU cores used by Snakemake to free up RAM per thread:
       ```bash
       snakemake --cores 4 --use-conda ...
       ```
