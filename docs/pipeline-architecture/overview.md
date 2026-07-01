# Pipeline Architecture Overview

The pipeline is organized into ten distinct, modular steps (M01 to M10), orchestrated via **Snakemake**. Each module has dedicated Conda environment files and outputs structured logs.

---

## High-Level Data Flow

The diagram below outlines how genomic data and selection scores flow from raw sequences to validated diagnostic assays:

```mermaid
flowchart TD
    subgraph Data Prep
        M01[M01: Acquisition] --> M02[M02: Quality Control]
        M02 --> M03[M03: MAFFT Alignment]
    end

    subgraph Evolution
        M03 --> M04[M04: IQ-TREE & Timetree]
        M04 --> M05[M05: HyPhy Selection Engine]
        M04 --> M06[M06: Ancestral Mapping]
    end

    subgraph Design & Verification
        M03 & M05 & M06 --> M07[M07: Bridge Builder Scoring]
        M07 --> M08[M08: PrimerLab Thermodynamic Design]
        M08 --> M09[M09: BLAST Specificity Check]
        M09 --> M10[M10: In-Silico PCR/RAA Emulation]
    end

    style M07 fill:#1a365d,stroke:#333,stroke-width:2px,color:#fff
    style M08 fill:#2b6cb0,stroke:#333,stroke-width:2px,color:#fff
    style M10 fill:#2b6cb0,stroke:#333,stroke-width:2px,color:#fff
```

---

## Directory Conventions

A standardized folder structure is maintained to ensure that the workflow remains modular and reproducible:

*   [`rules/`](file:///e:/Project/pandemic-preparedness/rules/): Snakemake rule files containing logic definitions (e.g., `rules/m01_acquisition.smk`).
*   [`scripts/`](file:///e:/Project/pandemic-preparedness/scripts/): Executable Python and bash scripts called by rules.
*   [`envs/`](file:///e:/Project/pandemic-preparedness/envs/): Conda environment definition YAML files, specifying package versions for replication.
*   [`config/`](file:///e:/Project/pandemic-preparedness/config/): Pathogen YAML profiles.
*   [`results/`](file:///e:/Project/pandemic-preparedness/results/): Output directory, segregated by `build_name` to prevent collisions.
