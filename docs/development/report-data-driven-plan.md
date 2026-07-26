---
title: "Making the Quarto report data-driven"
description: "Follow-up plan to stop report numbers drifting from the pipeline outputs."
---

# Making the Quarto report data-driven

## Why

Every number in `report/*.qmd` is currently **hand-transcribed** — there are no
executable code chunks. An audit found the report had drifted badly: tables were
spliced from two or three different runs, `sea-only`/`all-partial` figures had no
backing data, and dozens of values disagreed with the latest run. The values were
re-synced by hand once (see the git history), but hand-transcription guarantees
the same drift after the next run. This note is the plan to fix the root cause.

## Principle

A results page must **read the run's own artifacts at render time**, so a number
can never disagree with the data it claims to come from. Keep prose hand-written;
make every figure/table computed.

## Approach (incremental, one page at a time)

1. **Enable execution** on the results pages: add `engine: jupyter` (or knitr)
   and `execute: {freeze: auto}` so rendered numbers are frozen with the commit
   but recomputed when the data changes.
2. **A tiny shared loader** (`report/_load.py`): functions that read
   `results/<build>/m10/<mode>/coverage_summary.json`,
   `.../multiplex_candidates.json`, `blast_specificity_summary.csv`,
   `m02/qc_report.json`, `m07/candidate_targets_summary.csv`, and the M12
   `coverage_attribution.json`, returning tidy dataframes. One source of truth for
   paths, reused by every page.
3. **Replace transcribed tables with code chunks** that render from the loader:
   - `results_performance.qmd`: the coverage table + per-subtype table (global vs
     design/focal columns) straight from `coverage_summary.json` `breakdowns` /
     `design_db_breakdowns`.
   - `results_multiplex.qmd`: Top-N and the 10×10 matrix from
     `multiplex_candidates.json` (drop the static matrix).
   - `results_candidates.qmd`: candidate cards + BLAST table from
     `final_candidates.json` / `blast_specificity_summary.csv`.
   - `results_selection.qmd`: FUBAR/MEME site counts from the m05 JSONs.
   - `results_phylodynamics.qmd`: counts from `m02/qc_report.json` and the tree.
   - `supplementary.qmd`: the per-position frequency tables from
     `entropy_profile.tsv`.
4. **Scenario coverage**: drive the `indo-only` / `sea-only` / `all-partial`
   columns from whichever builds actually exist under `results/`; a missing build
   renders an explicit "awaiting re-run" cell instead of a stale number.
5. **CI guard** (optional): a check that every numeric claim in a page has a
   backing artifact, so a hand-typed number in a code-driven page fails review.

## Constraints

- Rendering then needs the run's `results/` present (it is a Quarto `resource`
  already). The publish workflow must have the data, or render locally/on HPC and
  push `_site` — decide per the deploy story.
- `results/` is gitignored (derived); the freeze cache (`_freeze/`) captures the
  rendered numbers for the committed site without committing the raw data.

## Interim state (already done, manual)

- Module pages (m02–m10) config snippets synced to the current `config.yaml`.
- Coverage two-lens (PCR-strict vs RAA-viable) + the PCR-cap finding documented
  in `m10.qmd` §1.3 and `results_performance.qmd`.
- Measles/triplex content removed; `sea-only`/`all-partial` marked "awaiting
  re-run"; broken `file:///` links made relative; discrete number mismatches
  fixed; voluminous transcribed tables (10×10 matrix, per-position frequencies)
  carry a "from a prior run" caveat pending this data-driven conversion.
