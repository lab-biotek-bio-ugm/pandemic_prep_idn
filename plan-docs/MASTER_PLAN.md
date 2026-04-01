# 🧬 Master Plan: Pandemic Preparedness Diagnostic Pipeline

> **Project**: Pandemic Preparedness — Development of Mux RDT POCT (Orthomyxovirus & Orthocoronavirinae)
> **Mandate**: Kementerian Kesehatan + WHO, kolaborasi BRIN-UGM
> **Role**: Bioinformatics Research Assistant
> **Goal**: Fully automated in-silico pipeline generating Multiplex RT-LAMP LFA panels (w/ PCR Benchmark) for SARS-CoV-2, Flu A, and Flu B.

---

## Pipeline Flow Diagram

![Pipeline Flow — 10 modules from data acquisition to report generation with nucleotide and protein analysis branches](./docs/pipeline_flow.png)

---

## Architecture Overview

Pipeline ini terdiri dari **10 modul** yang dikelompokkan menjadi 5 stage:

| Stage | Modules | Responsibility | Tech |
|-------|---------|---------------|------|
| **Data Preparation** | M1–M2 | Download & QC | BioPython, SeqKit |
| **Analysis** | M3–M6 | MSA, Phylogeny, Clade ID | MAFFT, IQ-TREE2, MACSE |
| **Bridge** | M7 | Translate phylo → primer targets | Custom Python |
| **PrimerLab** | M8–M9 | Primer design & validation | PrimerLab Genomic |
| **Output** | M10 | Comprehensive report | Quarto / Jinja2 |

---

## Technical Decisions (Locked-in Defaults)

Untuk memastikan pipeline siap dieksekusi pada fase *Proof of Concept*, keputusan tools/stack berikut telah dikunci sebagai default awal (dapat diubah nanti jika diminta oleh PI):

1. **Environment Management**: `Mamba` (lebih cepat dari Conda, terintegrasi mulus dengan Snakemake).
2. **Codon-Aware MSA**: `MACSE` (Menangani frameshift dan stop codon dengan menyisipkan gap `!` tanpa merusak reading frame).
3. **Phylogenetic Tree**: `IQ-TREE2` (Dilengkapi ModelFinder dan Ultrafast Bootstrap yang jauh lebih modern & cepat).
4. **Clade Assignment**: `Nextclade CLI` (Mendukung banyak virus respiratori sekaligus seperti SARS-CoV-2, Flu, RSV).
5. **dN/dS & Recombination**: `HyPhy` (Mencakup algoritma GARD untuk rekombinasi dan FEL/MEME untuk dN/dS).

---

## Module Details

### M1: Data Acquisition 📥

```
Input  : config.yaml (virus, date_range, country)
Output : raw_sequences.fasta, metadata.csv
Tools  : BioPython Entrez API, GISAID API, BV-BRC API
```

**Multi-Virus Config Templates (Gap 10):** Proses akuisisi diatur oleh template config spesifik per virus (misal: `influenza_a_h1n1.yaml`) yang otomatis memuat reference genome, metadata *segmented genome* (jika ada, *Gap 2*), dan parameter default lainnya agar *plug-and-play*.

| Step | Deskripsi |
|------|-----------|
| 1 | Query database dengan filter (organism, date, country, host) |
| 2 | Download sequences (FASTA format) |
| 3 | Download metadata (collection date, location, host, clade info) |
| 4 | Merge multi-source & deduplicate |
| 5 | Generate acquisition log |

**Database Sources:**

| Database | Virus Coverage | Access |
|----------|---------------|--------|
| NCBI GenBank / NCBI Virus | Semua virus | Open, Entrez API |
| GISAID | Influenza, SARS-CoV-2 | Registrasi, data agreement |
| BV-BRC (ViPR) | Semua virus | Open, API |

---

### M2: QC & Preprocessing 🧹

```
Input  : raw_sequences.fasta, metadata.csv
Output : cleaned_sequences.fasta, qc_report.json
Tools  : SeqKit, custom Python
```

| Step | Deskripsi | Threshold |
|------|-----------|-----------|
| 1 | Filter by sequence length | Within ±20% expected genome size |
| 2 | Filter by ambiguous bases (N%) | < 5% N content |
| 3 | Remove duplicate sequences | 100% identity dedup |
| 4 | Standardize FASTA headers | `>Accession|Virus|Country|Date|Clade` |
| 5 | **Subsampling Strategy (Gap 1)** | Batasi max sequences (temporal/geographic) agar MSA tidak hang jika input > 10.000 sekuens |
| 6 | Generate QC statistics summary | — |

---

### M3: Multiple Sequence Alignment 📐 — FORK POINT

```
Input  : cleaned_sequences.fasta
Output : nucleotide_alignment.aln, codon_alignment.aln, protein_alignment.aln
Tools  : MAFFT, MACSE/PAL2NAL, trimAl
```

Di modul ini, pipeline **bercabang** menjadi 2 track paralel:

```
cleaned_sequences.fasta
         │
    ┌────┴─────┐
    ▼          ▼
Track A     Track B
Nucleotide  Codon-Aware
MSA         MSA
(MAFFT)     (MACSE)
    │          │
    ▼          ├──► codon_alignment.aln  (nucleotide, codon-respected)
trimAl       ├──► protein_alignment.aln (auto-translate)
    │         └──► trimAl (codon-respect)
    ▼
nucleotide_alignment.aln
```

**Track A — Nucleotide MSA (Standard)**

- MAFFT default alignment
- trimAl untuk remove poorly aligned regions
- Output utama untuk primer design (conserved region detection)

**Track B — Codon-Aware MSA**

- Alignment per-codon (gap selalu kelipatan 3)
- Auto-translate → protein alignment (gratis!)
- Basis untuk dN/dS analysis dan protein phylogeny

---

### M4: Phylogenetic Analysis 🌳

```
Input  : nucleotide_alignment.aln, protein_alignment.aln
Output : nucleotide_tree.nwk, protein_tree.nwk, tree_visualization.html
Tools  : IQ-TREE2, (optional) BEAST2, ETE3/ggtree
```

| Analysis | Input | Tool | Output |
|----------|-------|------|--------|
| Nucleotide phylogeny | nucleotide_alignment.aln | IQ-TREE2 (ModelFinder + bootstrap) | nucleotide_tree.nwk |
| Protein phylogeny | protein_alignment.aln | IQ-TREE2 | protein_tree.nwk |
| Time-scaled tree (opt) | nucleotide_alignment.aln + dates | BEAST2 | timetree.nwk |
| Visualization | *.nwk | ETE3 / ggtree / iTOL | tree_viz.html |

---

### M5: Sequence Analysis 🔍

```
Input  : nucleotide_alignment.aln, codon_alignment.aln, metadata.csv
Output : conserved_regions.bed, variable_regions.bed, snp_matrix.csv,
         dnds_results.csv, epitope_candidates.json (optional)
Tools  : Custom Python, SNP-sites, PAML/HyPhy, BepiPred (optional)
```

**Nucleotide Track:**

| Step | Deskripsi | Output |
|------|-----------|--------|
| 1 | Detect conserved regions (>95% identity) | conserved_regions.bed |
| 2 | Detect variable regions (clade-defining SNPs) | variable_regions.bed |
| 3 | Extract full SNP matrix | snp_matrix.csv |
| 4 | Map mutations to gene annotations | annotated_mutations.json |

**Protein Track (dari codon-aware MSA):**

| Step | Deskripsi | Output |
|------|-----------|--------|
| 5 | dN/dS analysis (selection pressure) | dnds_results.csv |
| 6 | Conserved protein region detection | conserved_proteins.json |
| 7 | **Recombination Detection (Gap 3)** | Deteksi breakpoint dengan HyPhy GARD untuk flag sequence rekombinan (mencegah *misleading phylogeny*) |
| 8 | (Optional) Epitope prediction | epitope_candidates.json |

> **dN/dS & Recombination insight**: Region dengan dN/dS < 1 = purifying selection = functionally conserved = target terbaik. Rekombinasi wajib dideteksi di awal agar tidak merusak interpretasi kekerabatan.

---

### M6: Clade Classification 🏷️

```
Input  : nucleotide_tree.nwk, snp_matrix.csv, metadata.csv
Output : clade_assignments.csv, clade_mutations.json,
         clade_X_consensus.fasta (per clade), clade_distribution.json
Tools  : Augur/Nextstrain, Pangolin (SARS-CoV-2), custom
```

| Step | Deskripsi |
|------|-----------|
| 1 | Assign clades (existing nomenclature: Nextstrain, Pango, HA subtype) |
| 2 | Identify clade-defining mutations (synapomorphies) |
| 3 | Generate consensus sequence per clade |
| 4 | Map clade distribution (temporal + geographic) |

---

### M7: Primer Target Preparation 🔀 — BRIDGE MODULE

```
Input  : conserved_regions.bed, variable_regions.bed, clade_consensus.fasta,
         dnds_results.csv, epitope_candidates.json
Output : *_target.fasta + primerlab_config.yaml per target,
         batch_config.yaml, target_manifest.json, epitope_report.json
Tools  : Custom Python (BioPython)
```

> **Ini modul "translator"** antara dunia phylogenetics dan PrimerLab.
> PrimerLab tidak perlu tau context virus — modul ini yang menyiapkan semuanya.

```
┌─────────────────────────────────────────────────────────────────┐
│                    MODULE 7 — 4 Sub-workflows                   │
│                                                                 │
│  ┌──────────────────┐  ┌──────────────────┐                     │
│  │ 7A: UNIVERSAL    │  │ 7B: CLADE-       │                     │
│  │ PRIMER TARGETS   │  │ SPECIFIC TARGETS │                     │
│  │                  │  │                  │                     │
│  │ conserved region │  │ per clade:       │                     │
│  │ → consensus seq  │  │ mutation region  │                     │
│  │ → FASTA + config │  │ + flanking 200-  │                     │
│  │                  │  │ 500bp → FASTA +  │                     │
│  │                  │  │ config           │                     │
│  └────────┬─────────┘  └────────┬─────────┘                     │
│           │                     │                               │
│  ┌────────┴─────────┐  ┌───────┴──────────┐                    │
│  │ 7C: MULTIPLEX    │  │ 7D: EPITOPE      │                    │
│  │ PREPARATION      │  │ REPORT (opt)     │                    │
│  │                  │  │                  │                    │
│  │ batch config +   │  │ ranked epitope   │                    │
│  │ compatible Tm    │  │ candidates →     │                    │
│  │ range            │  │ handoff wet lab  │                    │
│  └──────────────────┘  └──────────────────┘                    │
└─────────────────────────────────────────────────────────────────┘
```

---

### M8: Primer/Probe Design 🎯 — PrimerLab Genomic

```
Input  : *_target.fasta + primerlab_config.yaml (dari M7)
Output : primers_*.json, probes.json, design_reports.html
Tools  : PrimerLab Genomic (primerlab-genomic)
```

**Integrasi via Python API:**

```python
from primerlab.api.public import design_pcr_primers, design_qpcr_assays

# Per target
result = design_pcr_primers(sequence, config)
result = design_qpcr_assays(sequence, config)  # TaqMan
```

| Design Type | PrimerLab Command | Status |
|-------------|------------------|--------|
| PCR primers | `design_pcr_primers()` | ✅ Ready |
| qPCR + TaqMan probe | `design_qpcr_assays()` | ✅ Ready |
| Nested PCR | `nested-design` | ✅ Ready |
| Semi-Nested PCR | `seminested-design` | ✅ Ready |
| Batch processing | Built-in + SQLite cache | ✅ Ready |

---

### M9: In-Silico Validation ✅ — PrimerLab Genomic

```
Input  : primers from M8 + cleaned_sequences.fasta (M2)
Output : validation_report.json, ranked_primers.json, cross_reactivity.csv
Tools  : PrimerLab Genomic
```

| Validation Step | PrimerLab Feature | Status |
|----------------|------------------|--------|
| In-silico PCR simulation | `primerlab insilico` | ✅ Ready |
| Off-target BLAST **(Gap 6)** | `core/offtarget/` vs *Human Genome (GRCh38) & Respiratory Flora* | ✅ Ready (needs DB setup) |
| Species cross-reactivity | `core/species/` | ✅ Ready |
| Multiplex compatibility | `primerlab check-compat` | ✅ Ready |
| Dimer matrix | `primerlab dimer-matrix` | ✅ Ready |
| Primer ranking | Custom scoring script | 🆕 Build in pipeline |

---

### M10: Report Generation 📊

```
Input  : ALL outputs from M1–M9
Output : report.html, report.pdf, figures/
Tools  : Quarto / Jinja2, plotly, matplotlib
```

**Report Contents:**

1. Executive Summary
2. Data Overview & QC Statistics
3. Phylogenetic Analysis (interactive tree)
4. Clade Distribution & Temporal Trends
5. Primer Target Rationale (conserved/variable regions, dN/dS)
6. Primer Design Results (embedded PrimerLab output)
7. Validation Results (specificity, sensitivity, cross-reactivity)
8. Recommended Diagnostic Panel
9. **MIQE Compliance Checklist (Gap 7)** (Tm, amplicon size, probe chemistry, dll)
10. Methods & Parameters (full reproducibility)
11. Appendices

---

## Scope & Output Summary

### ✅ IN SCOPE (Primary Nucleic Acid Diagnostics)

| Output | Deskripsi | Module |
|--------|-----------|--------|
| **Multiplex RT-LAMP LFA** | **TARGET UTAMA:** Isothermal amplification (6 primers) untuk POCT | M8-M9 |
| **Universal PCR Benchmark** | **GOLD STANDARD:** Assay TaqMan qPCR untuk komparasi / validasi klinis | M8-M9 |
| Phylogenetic tree | Nuc + Protein | M4 |
| Clade classification | Assignment + consensus | M6 |
| In-silico validation | isPCR, BLAST exclusivitas, cross-reactivity matrix | M9 |
| Comprehensive report | HTML + PDF | M10 |

### ✅ IN SCOPE (Protein-Level — dari Codon-Aware MSA)

| Output | Deskripsi | Module |
|--------|-----------|--------|
| Protein phylogeny | Dari auto-translate codon MSA | M4 |
| dN/dS analysis | Selection pressure per region | M5 |

### ❌ OUT OF SCOPE (Future/Other Teams)

| Item                      | Alasan                           |
|---------------------------|----------------------------------|
| CRISPR guide RNA (Cas12/13) | Terlalu mahal untuk RDT massal, butuh amplifikasi ganda |
| Epitope/Antibody production | Domain wet-lab immunology, bukan in-silico primer design |
| Metabolite-based (GeNose) | Bukan molecular biology          |

---

## PrimerLab Enhancement Roadmap

Agar integrasi M7→M8→M9 *seamless*, PrimerLab perlu enhancement berikut:

### Prioritas 1 — Harus ada sebelum integrasi

| Enhancement                 | Deskripsi                                                                                                          | Effort |
|-----------------------------|--------------------------------------------------------------------------------------------------------------------|--------|
| **Multiplex Panel Builder** | Dari `check-compat` (cek saja) → Panel Builder (seleksi kombinasi terbaik: ΔTm < 2°C + no cross-dimer)             | Medium |
| **Force 3'-End Constraint** | Config option untuk memaksa Primer3 meletakkan 3'-end primer tepat di posisi mutasi (wajib untuk clade-specific)     | Small  |
| **Batch Manifest CLI**      | `primerlab run batch-manifest --manifest target_manifest.json` → proses puluhan target sekaligus                   | Medium |

### Prioritas 2 — Nice to have

| Enhancement                  | Deskripsi                                              | Effort |
|------------------------------|--------------------------------------------------------|--------|
| **Degenerate base support**  | Handle IUPAC codes (R, Y, S, W, dll) di input sequence | Small  |
| **Primer coverage analysis** | Hitung berapa % sequences di dataset yang match primer | Medium |

### Prioritas 3 — Future work (JANGAN sekarang)

| Enhancement          | Deskripsi                                      | Effort |
|----------------------|------------------------------------------------|--------|
| **RT-LAMP module**   | Logic 6-primer design (F3, B3, FIP, BIP, LF, LB) | Large  |
| **CRISPR guide RNA** | PAM detection + crRNA design + RPA primers       | Large  |

---

## Tech Stack

| Component | Tool | Bahasa |
|-----------|------|--------|
| **Orchestration** | Snakemake | Python |
| **Sequence handling** | BioPython, SeqKit | Python, Go |
| **Nucleotide MSA** | MAFFT | C |
| **Codon-aware MSA** | MACSE / PAL2NAL | Java / Perl |
| **Alignment trimming** | trimAl | C++ |
| **Phylogenetics** | IQ-TREE2 | C++ |
| **Tree visualization** | ETE3, ggtree, iTOL | Python, R, Web |
| **dN/dS analysis** | PAML / HyPhy | C / C++ |
| **Epitope prediction** | BepiPred | Web/Python |
| **Primer design** | **PrimerLab Genomic** | Python |
| **Validation** | **PrimerLab Genomic** | Python |
| **Reporting** | Quarto / Jinja2, plotly | Python |
| **Package management** | Conda / Mamba | — |
| **Containerization** | Docker / Singularity | — |

---

## Development Phases (In-Silico Timeline: April - Juni)


### Phase 1: PrimerLab Engine Refactoring
- [ ] Implement Batch Manifest CLI & `SEQUENCE_INCLUDED_REGION`
- [ ] Integrasi algoritma desain RT-LAMP (6-primer mode)
- [ ] Bangun N-Plex Heterodimer Matrix (ViennaRNA Cross-Check)
- [ ] Unit Test PrimerLab dengan dummy sequence

### Phase 2: Pipeline Backbone M1-M3
- [ ] Setup Snakemake project skeleton & Conda Envs
- [ ] Implement M1 & M2 (Target: SARS-CoV-2, Flu A, Flu B)
- [ ] Implement M3 (MAFFT + Codon-aware)
- [ ] Dry-run akuisisi data 10.000 sekuens

### Phase 3: Evolutionary Engine M4-M6
- [ ] Implement M4 (IQ-TREE2) & M6 (Clade Mapping)
- [ ] Implement M5 (Conserved logic + HyPhy FUBAR/MEME/GARD)
- [ ] Validasi output BED file dan CSV dN/dS

### Phase 4: Bridge & In-Silico Validation M7-M10
- [ ] Implement M7 (Bridge to Manifest YAML)
- [ ] Run End-to-End M1 -> M9 (Desain PCR Benchmark + Mux RT-LAMP)
- [ ] Validasi Off-target BLAST vs Human/Resp Flora
- [ ] Polish Quarto Report generation (M10)

---

## Project File Structure

```
pandemic-diagnostic-pipeline/
├── Snakefile                          # Main workflow
├── config/
│   ├── config.yaml                    # User parameters
│   └── defaults.yaml                  # Default values
├── rules/
│   ├── m01_acquisition.smk
│   ├── m02_qc.smk
│   ├── m03_alignment.smk              # nuc + codon-aware fork
│   ├── m04_phylogeny.smk
│   ├── m05_analysis.smk               # conserved/variable + dN/dS
│   ├── m06_classification.smk
│   ├── m07_target_prep.smk            # bridge to PrimerLab
│   ├── m08_primer_design.smk          # PrimerLab API calls
│   ├── m09_validation.smk             # PrimerLab validation
│   └── m10_report.smk
├── scripts/
│   ├── download_ncbi.py
│   ├── download_gisaid.py
│   ├── qc_filter.py
│   ├── find_conserved.py
│   ├── find_variable.py
│   ├── run_dnds.py
│   ├── predict_epitopes.py
│   ├── assign_clades.py
│   ├── make_consensus.py
│   ├── prepare_targets.py
│   ├── run_primerlab.py
│   ├── validate_primers.py
│   ├── rank_primers.py
│   └── generate_report.py
├── templates/
│   └── report_template.html           # Jinja2 / Quarto template
├── envs/
│   ├── base.yaml                      # Core Python deps
│   ├── alignment.yaml                 # MAFFT, trimAl, MACSE
│   ├── phylogeny.yaml                 # IQ-TREE2
│   ├── analysis.yaml                  # PAML, HyPhy, SNP-sites
│   └── primerlab.yaml                 # primerlab-genomic
├── docs/
│   ├── pipeline_flow.png              # Flow diagram
│   └── user_guide.md
├── results/                           # Output per module
│   ├── m01_acquisition/
│   ├── m02_qc/
│   ├── m03_alignment/
│   ├── m04_phylogeny/
│   ├── m05_analysis/
│   ├── m06_classification/
│   ├── m07_targets/
│   ├── m08_primers/
│   ├── m09_validation/
│   └── m10_report/
├── environment.yml                    # Master conda environment
├── Dockerfile
├── README.md
└── .gitignore
```

---

## Questions for Supervisor Meeting

### Critical (HARUS ditanya)

1. Virus apa saja yang masuk scope? (Flu A/B, SARS-CoV-2, RSV, hMPV, Adenovirus?)
2. Fokus di RT-qPCR saja atau juga conventional PCR?
3. Apakah data publik (GISAID + GenBank) acceptable jika data Kemenkes belum tersedia?
4. Apakah protein-level analysis (epitope prediction) termasuk scope?
5. Apakah multiplex panel (multiple virus dalam 1 tube) diinginkan?

### Important

1. Timeline — kapan deliverable pertama?
2. Compute resources — akses ke HPC BRIN/UGM?
3. Apakah primer akan divalidasi di wet lab? Siapa?
4. Siapa yang akan menggunakan pipeline ini selain kamu?
5. Apakah ada rencana publikasi?

---

## Key Talking Points for Meeting

> 1. **Pipeline ini fully automated dan reproducible** — dari download data sampai recommended diagnostic panel, hanya butuh 1 config file.
>
> 2. **Dual analysis track** — nucleotide MSA untuk primer design + codon-aware MSA untuk protein insight dan selection pressure analysis, dari satu input.
>
> 3. **PrimerLab Genomic** sebagai core engine primer design yang sudah battle-tested (1286 tests) — bukan build from scratch.
>
> 4. **Update cepat** — ketika varian baru muncul, pipeline tinggal re-run. Primer recommendation baru dalam hitungan jam, bukan minggu.
>
> 5. **Extensible** — arsitektur siap untuk RT-LAMP dan CRISPR-based diagnostic di masa depan jika ada kebutuhan point-of-care testing.
