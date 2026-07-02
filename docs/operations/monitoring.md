---
title: "Pipeline Monitoring & Run Management"
description: "Best practices for monitoring long-running pipeline executions, interpreting progress, managing Snakemake checkpoints, and handling partial failures on HPC and local workstations."
---

# Pipeline Monitoring & Run Management

This page covers practical operational guidance for managing pipeline runs from launch through completion, including progress monitoring, checkpoint behaviour, partial-failure recovery, and resource tracking.

---

## Standard Run Command

The canonical command for a full pipeline run with reproducible environments:

```bash
snakemake \
    --use-conda \
    --configfile config/config_influenza_a_mp.yaml \
    --cores 16 \
    --rerun-incomplete \
    --printshellcmds \
    2>&1 | tee logs/snakemake_$(date +%Y%m%d_%H%M%S).log
```

| Flag | Purpose |
|---|---|
| `--use-conda` | Activates per-rule conda environments defined in `envs/` |
| `--cores N` | Maximum parallel CPU cores; Snakemake schedules jobs to stay within this limit |
| `--rerun-incomplete` | Re-executes any rules with incomplete output files from a previous interrupted run |
| `--printshellcmds` | Prints the exact shell command for each rule before executing — critical for debugging |
| `2>&1 \| tee` | Captures both stdout and stderr to a timestamped log file without suppressing terminal output |

---

## HPC / Slurm Cluster Execution

For large datasets or computationally intensive steps (M05 HyPhy, M04 IQ-TREE2), submit to a Slurm cluster using a Snakemake profile:

```bash
# Install the Snakemake Slurm profile (one-time setup)
pip install snakemake-executor-plugin-slurm

# Run on the cluster
snakemake \
    --use-conda \
    --configfile config/config_influenza_a_mp.yaml \
    --executor slurm \
    --default-resources slurm_partition=standard mem_mb=16000 \
    --set-resources m05_hyphy_fubar:mem_mb=32000 \
    --set-resources augur_tree:mem_mb=32000 \
    --jobs 50 \
    2>&1 | tee logs/snakemake_cluster.log
```

**Per-rule resource overrides** are needed for memory-intensive steps:

| Rule | Minimum RAM | Notes |
|---|---|---|
| `augur_align` | 8 GB | Scales with sequence count × length |
| `augur_tree` | 16 GB | IQ-TREE2 full alignment tree inference |
| `augur_refine` | 4 GB | TreeTime with >3,000 taxa |
| `m05_hyphy_fubar` | 8 GB | Per gene; run multiple genes in parallel |
| `m05_hyphy_gard` | 32 GB | GARD is memory-intensive for long alignments |
| `m10_insilico_pcr` | 4 GB | seqkit amplicon on large sequence database |

---

## Monitoring Progress

### Real-time Log Tailing

Monitor the main Snakemake log and any specific rule log simultaneously:

```bash
# Follow the main Snakemake log
tail -f logs/snakemake_*.log

# Follow a specific rule log (e.g., HyPhy FUBAR for M1 gene)
tail -f logs/influenza-a-mp-sea/m05_hyphy_fubar_M1.log

# Watch all M05 logs in real time (requires tmux or GNU screen)
tail -f logs/influenza-a-mp-sea/m05_hyphy_*.log
```

### Checking DAG Progress

```bash
# Count completed rules (files exist)
snakemake \
    --use-conda \
    --configfile config/config_influenza_a_mp.yaml \
    --dryrun --quiet 2>&1 | grep "remaining"

# Generate a visual DAG of the remaining jobs
snakemake \
    --use-conda \
    --configfile config/config_influenza_a_mp.yaml \
    --dryrun \
    --dag 2>/dev/null | dot -Tpng > logs/dag_remaining.png
```

### Expected Runtime by Module (Reference: 16-core workstation, 3,000 sequences, Flu A MP)

| Module | Typical Runtime | Bottleneck |
|---|---|---|
| M01 Download / Copy | < 1 min | I/O |
| M01 Nextclade | 2–5 min | CPU |
| M02 QC + Filter | 1–3 min | I/O + pandas |
| M03 Alignment (MAFFT) | 5–20 min | CPU (16 threads) |
| M04 IQ-TREE2 | 15–60 min | CPU (16 threads) |
| M04 TreeTime | 10–30 min | Single-threaded |
| M05 HyPhy FUBAR (per gene) | 10–60 min | CPU (8 threads) |
| M05 HyPhy MEME (per gene) | 20–120 min | CPU (8 threads) |
| M05 HyPhy GARD | 12–48 hr | CPU — HPC recommended |
| M06 Ancestral + Export | 5–15 min | CPU |
| M07 Bridge Builder | 2–10 min | Python I/O |
| M08 PrimerLab (per region) | 10–60 min | CPU + ViennaRNA |
| M09 BLAST local (per region) | 5–30 min | Disk I/O |
| M09 BLAST remote (per region) | 30–120 min | Network + NCBI API |
| M10 In-silico PCR | 5–20 min | seqkit + Python |

<Note>
M05 HyPhy GARD is the single most compute-intensive step. For short segments (< 1,200 bp) with < 300 sequences, GARD may complete in 1–2 hours. For whole-genome targets or large datasets, always submit GARD to an HPC cluster.
</Note>

---

## Checkpoint Behaviour (M07 Bridge Builder)

M07 is declared as a Snakemake **`checkpoint`** rather than a regular `rule`. This is a critical architectural detail:

- A checkpoint forces Snakemake to **re-evaluate the DAG** after it completes, because the number of output regions (and therefore the number of downstream M08/M09/M10 jobs) is not known until M07 finishes and writes the manifest.
- After M07 completes, Snakemake reads `results/{BUILD}/m07/primer_targets_manifest.yaml` to discover the list of `region_id` values, then dynamically schedules the appropriate number of M08, M09, and M10 instances.

**Practical implication:** If M07 fails, all downstream M08–M10 jobs are blocked. Fix M07 before attempting to re-run:

```bash
# Re-run only M07 (forces re-evaluation of entire downstream DAG)
snakemake \
    --use-conda \
    --configfile config/config_influenza_a_mp.yaml \
    --forcerun m07_bridge \
    --cores 16
```

---

## Partial Failure Recovery

### Resuming After Interruption

Snakemake automatically tracks which output files exist. Simply re-run the same command with `--rerun-incomplete`:

```bash
snakemake \
    --use-conda \
    --configfile config/config_influenza_a_mp.yaml \
    --cores 16 \
    --rerun-incomplete
```

If the pipeline was killed mid-rule (e.g., by a cluster timeout), the output file may be partially written. `--rerun-incomplete` detects these and marks them for re-execution.

### Unlocking a Locked Working Directory

If Snakemake was killed without a clean shutdown, it may leave a `.snakemake/locks/` directory that prevents new runs:

```bash
snakemake --unlock --configfile config/config_influenza_a_mp.yaml
```

### Forcing a Specific Module to Re-run

```bash
# Force re-run of M02 and all downstream rules
snakemake \
    --use-conda \
    --configfile config/config_influenza_a_mp.yaml \
    --forcerun qc_prefilter \
    --cores 16

# Force re-run of a specific HyPhy gene rule
snakemake \
    --use-conda \
    --configfile config/config_influenza_a_mp.yaml \
    --forcerun m05_hyphy_fubar \
    --cores 8
```

---

## Dry-Run Validation

Always run a dry-run before any full pipeline execution to validate the DAG without executing anything:

```bash
snakemake \
    --use-conda \
    --configfile config/config_influenza_a_mp.yaml \
    --dryrun \
    --printshellcmds \
    2>&1 | tee logs/dryrun_check.log
```

Check the output for:
1. Expected number of jobs (verify it matches your target)
2. Any `MissingInputException` errors (missing config files or input data)
3. Any `AmbiguousRuleException` errors (rule ordering conflicts)
4. Resource estimates in the `--dryrun` summary

---

## Disk Space Management

### Estimating Storage Requirements

| Data | Size per 3,000 sequences |
|---|---|
| Input FASTA | 5–15 MB |
| Aligned FASTA (`m03/aligned.fasta`) | 5–20 MB |
| Conda environments (`envs/`) | 2–8 GB (one-time) |
| HyPhy FUBAR/MEME outputs | 50–200 MB per gene |
| BLAST local database (`resources/blast_db/`) | 3–15 GB |
| Total results per build | 200–800 MB (excluding BLAST DB) |

### Cleaning Intermediate Files

After a successful run, large intermediate files can be safely removed:

```bash
# Remove Nextclade intermediate directory (keeps nextclade_output.tsv)
rm -rf results/{BUILD}/m01/nextclade/

# Remove HyPhy per-gene FASTA copies (original in m03 is preserved)
rm results/{BUILD}/m05/cleaned_*_NT.fasta

# Remove raw ML tree (time-calibrated tree.nwk is the authoritative output)
rm results/{BUILD}/m04/tree_raw.nwk
```

<Warning>
Do **not** delete `results/{BUILD}/m07/primer_targets_manifest.yaml` after M07 completes. This file is read at pipeline startup by M10's DAG evaluation function. Deleting it requires a full M07 re-run to regenerate.
</Warning>
