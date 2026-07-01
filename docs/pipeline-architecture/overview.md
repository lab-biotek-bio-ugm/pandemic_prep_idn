# Pipeline Architecture Overview

The pipeline is organized into ten distinct, modular steps (M01 to M10), orchestrated via **Snakemake**. Each module has dedicated Conda environment files and outputs structured logs.

---

## High-Level Data Flow

The diagram below outlines how genomic data and selection scores flow from raw sequences to validated diagnostic assays:

```mermaid
flowchart LR
    subgraph Prep [Data Prep]
        M01[M01: Acquisition] --> M02[M02: QC & Filter] --> M03[M03: Alignment]
    end

    subgraph Evo [Evolutionary Analysis]
        M03 --> M04[M04: IQ-TREE & Clock]
        M04 --> M05[M05: HyPhy Selection]
        M04 --> M06[M06: Ancestral Mapping]
    end

    subgraph Design [Design & Verification]
        M03 --> M07[M07: Bridge Builder]
        M05 --> M07
        M06 --> M07
        M07 --> M08[M08: Primer Design]
        M08 --> M09[M09: Specificity Check]
        M09 --> M10[M10: In-Silico Validation]
    end
```

---

## Directory Conventions

A standardized folder structure is maintained to ensure that the workflow remains modular and reproducible:

*   [`rules/`](file:///e:/Project/pandemic-preparedness/rules/): Snakemake rule files containing logic definitions (e.g., `rules/m01_acquisition.smk`).
*   [`scripts/`](file:///e:/Project/pandemic-preparedness/scripts/): Executable Python and bash scripts called by rules.
*   [`envs/`](file:///e:/Project/pandemic-preparedness/envs/): Conda environment definition YAML files, specifying package versions for replication.
*   [`config/`](file:///e:/Project/pandemic-preparedness/config/): Pathogen YAML profiles.
*   [`results/`](file:///e:/Project/pandemic-preparedness/results/): Output directory, segregated by `build_name` to prevent collisions.
