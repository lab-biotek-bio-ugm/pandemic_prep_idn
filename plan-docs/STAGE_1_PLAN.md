# 📥 Stage 1: Data Preparation — Detailed Plan

> **Scope:** Module 1 (Data Acquisition) & Module 2 (QC & Preprocessing)
> **Objective:** Mendapatkan raw sequences, melakukan standardisasi, QC, dan subsampling agar siap masuk MSA.

---

## 1. Config Architecture

### 1.1 File Structure

```text
config/
├── config.yaml                     # Main user config (EDITABLE)
└── templates/
    ├── sars-cov-2.yaml             # Single segment
    ├── influenza_a_h1n1_ha.yaml    # Segmented (HA)
    ├── influenza_a_h1n1_m.yaml     # Segmented (M)
    └── rsv_a.yaml                  # Single segment
data/
└── local/                          # Tempat user meletakkan data lokal
    ├── *.fasta                     # File FASTA (boleh 1 atau banyak)
    └── metadata.tsv                # Metadata (WAJIB mengikuti template)
```

### 1.2 `config.yaml` — User Config

```yaml
# ──── API Credentials ────
# PENTING: Jangan pernah commit file ini ke Git jika berisi API key!
api_keys:
  ncbi_email: "your.email@brin.go.id"
  ncbi_api_key: "your_api_key_here"   # Daftar di: https://www.ncbi.nlm.nih.gov/account/settings/

# ──── Targets (Batch — 3 Virus Utama Sesuai TOR) ────
targets:
  - name: "sarscov2_panel"
    virus_template: "templates/sars-cov-2.yaml"
    database: "ncbi"
    date_start: "2023-01-01"
    date_end: "2024-03-31"
    country: "Indonesia"

  - name: "flu_a_panel"
    virus_template: "templates/influenza_a_h1n1_ha.yaml"
    database: "ncbi"
    date_start: "2023-01-01"
    date_end: "2024-03-31"
    country: "Indonesia"

  - name: "flu_b_panel"
    virus_template: "templates/influenza_b.yaml"
    database: "ncbi"
    date_start: "2023-01-01"
    date_end: "2024-03-31"
    country: "Indonesia"

# ──── QC & Subsampling (M2) ────
qc:
  max_n_content: 5.0                    # Max % basa 'N' yang diizinkan
  length_filter: true                   # On/off filter berdasar panjang
  dedup: true                           # On/off deduplication
  subsample:
    enabled: true
    strategy: "random"                  # random | temporal | geographic
    max_sequences: 1000                 # Batas atas untuk cegah komputasi jebol
    seed: 42                            # Reproducible random sampling
```

### 1.3 Virus Template — `templates/sars-cov-2.yaml`

```yaml
# ──── Identitas Biologis ────
species: "Severe acute respiratory syndrome coronavirus 2"
short_name: "SARS-CoV-2"

# ──── Database-Specific Aliases (Solusi Naming Issue) ────
# Setiap database punya cara pencarian berbeda.
# Script M1 akan membaca alias sesuai database yang dipilih user.
aliases:
  ncbi_taxid: 2697049                   # TaxID = ID numerik absolut, kebal terhadap perbedaan penamaan
  gisaid_query: "hCoV-19"               # Nomenklatur resmi GISAID
  bvbrc_query: "SARS-CoV-2"             # Nomenklatur BV-BRC

# ──── Genome Properties ────
segmented: false
expected_length: 29903                  # bp (basis pairs) — panjang reference genome
length_tolerance_percent: 10.0          # Toleransi ± 10% → accept 26913–32893 bp

# ──── Reference Sequence ────
# Bukan yang di-download massal. Ini SATU sequence rujukan klasik
# yang digunakan sebagai "Titik Nol" untuk koordinat mutasi (M5) dan posisi primer (M8).
reference_accession: "NC_045512.2"
```

### 1.4 Virus Template — `templates/influenza_a_h1n1_ha.yaml`

```yaml
species: "Influenza A virus (A/H1N1)"
short_name: "Flu-A-H1N1-HA"
aliases:
  ncbi_taxid: 114727
  gisaid_query: "A / H1N1"
  bvbrc_query: "Influenza A virus"
segmented: true
target_segment: "HA"                    # Nama segmen target
ncbi_segment_filter: "hemagglutinin"    # Kata kunci tambahan untuk filter ESearch
expected_length: 1700
length_tolerance_percent: 15.0
reference_accession: "NC_026433.1"
```

---

## 2. Module 1: Data Acquisition (M1)

### 2.1 Goal

Menghasilkan 1 file FASTA gabungan (raw) dan 1 file TSV metadata **per target** yang terstandardisasi.

### 2.2 Input & Output

| Direction | File | Keterangan |
|-----------|------|------------|
| **Input** | `config.yaml` + `template.yaml` | Parameter pencarian / path lokal |
| **Output** | `results/m01_acquisition/{target_name}/raw_sequences.fasta` | FASTA gabungan, header sudah baku |
| **Output** | `results/m01_acquisition/{target_name}/metadata.tsv` | Metadata terstandarisasi |
| **Output** | `results/m01_acquisition/{target_name}/acquisition_log.json` | Log proses (jumlah downloaded, errors, timestamp) |

### 2.3 Output Schema

**`metadata.tsv` — Kolom Standar (output M1, berlaku untuk SEMUA sumber data)**

| Kolom | Tipe | Wajib | Contoh | Keterangan |
|-------|------|-------|--------|------------|
| `accession` | string | ✅ | `MT123456` | ID unik sequence |
| `date` | string | ✅ | `2024-03-15` | Format YYYY-MM-DD |
| `country` | string | ✅ | `Indonesia` | Negara asal |
| `length` | int | ✅ | `29850` | Panjang sequence (bp) |
| `source` | string | ✅ | `ncbi` / `local` | Asal data |
| `original_header` | string | ✅ | `>MT123456.1 SARS-CoV-2/human/...` | Header asli (backup) |

**`raw_sequences.fasta` — Format Header Baku**

```text
>MT123456|Indonesia|2024-03-15
ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAAT...
>ON789012|Indonesia|2024-02-10
ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAAT...
```

### 2.4 Script: `scripts/m01_acquire_data.py`

Script ini adalah **Router + Validator + Standardizer**.

**Logic Flow:**

```text
START
  │
  ├─ 1. Load config.yaml + virus template
  │
  ├─ 2. Validate config (cek template file exists, cek format date, dll)
  │     └─ FAIL? → Raise error dengan pesan jelas, pipeline berhenti
  │
  ├─ 3. Route berdasarkan `database`:
  │     │
  │     ├─ [NCBI] ──────────────────────────────────────────────────
  │     │   a. Set Entrez.email + Entrez.api_key
  │     │   b. Build query: txid{taxid}[Organism]
  │     │      + "{country}"[Country]
  │     │      + ("{date_start}"[PDAT] : "{date_end}"[PDAT])
  │     │      + (jika segmented: AND "{ncbi_segment_filter}"[Gene])
  │     │   c. ESearch → dapatkan daftar ID (bisa ribuan)
  │     │   d. EFetch FASTA dalam batch (chunk 500, retry 3x per chunk)
  │     │   e. EFetch metadata (GenBank XML) → parse ke TSV
  │     │   f. Auto-download reference_accession jika belum ada
  │     │
  │     ├─ [LOCAL] ─────────────────────────────────────────────────
  │     │   a. Scan semua *.fasta / *.fa di local_fasta_dir
  │     │   b. Load local_metadata (TSV)
  │     │   c. Validate: cek kolom wajib (accession, date, country)
  │     │   d. Validate: cek format date (YYYY-MM-DD)
  │     │   e. Validate: cross-check accession di TSV vs ID di FASTA
  │     │   f. FAIL? → error di baris/file mana yang bermasalah
  │     │
  │     └─ [GISAID] → raise NotImplementedError("TODO: GISAID API logic (Post-PoC)")
  │
  ├─ 4. Standardize Header FASTA → >{accession}|{country}|{date}
  │     └─ Simpan header asli di kolom 'original_header' metadata
  │
  ├─ 5. Generate acquisition_log.json (timestamp, counts, errors)
  │
  └─ 6. Write output files
END
```

**Error Handling & Resilience:**

| Skenario | Handling |
|----------|---------|
| NCBI timeout / network error | Retry 3x per batch (exponential backoff) |
| 0 sequences found | Raise error, jangan lanjut ke M2 |
| Local FASTA kosong | Raise error dengan path file yang bermasalah |
| Accession di TSV tidak ada di FASTA | Raise error dengan list accession yang hilang |
| Tanggal format salah | Raise error dengan baris TSV yang salah |
| Template file tidak ditemukan | Raise error di awal sebelum proses apapun |

### 2.5 Local Data Input Requirements

User yang membawa data sendiri **WAJIB** menyiapkan:

1. **File FASTA** (`.fasta` atau `.fa`) — boleh 1 file gabungan atau banyak file terpisah
2. **File `metadata.tsv`** — Tab-separated, minimal 3 kolom wajib

**Contoh `metadata.tsv` (Local):**

```text
accession	date	country
SAMPLE_001	2024-01-15	Indonesia
SAMPLE_002	2024-02-20	Indonesia
SAMPLE_003	2024-03-10	Malaysia
```

> **Kenapa TSV bukan CSV?** Karena koma (`,`) sering muncul di nama lokasi, isolate, dll. Tab sebagai separator jauh lebih aman dan merupakan standar di bioinformatics.

**Mapping FASTA ↔ Metadata:**
Script M1 membaca kata pertama setelah `>` di FASTA sebagai ID, lalu mencocokkan dengan kolom `accession` di TSV. Contoh:

```text
File FASTA:                               metadata.tsv:
>SAMPLE_001 some description here   →     accession = SAMPLE_001 ✅ match
>SAMPLE_002                         →     accession = SAMPLE_002 ✅ match
>UNKNOWN_ID                         →     ❌ tidak ada di TSV → ERROR
```

### 2.6 Snakemake Rule (Draft)

```python
rule m01_acquisition:
    output:
        fasta = "results/m01_acquisition/{target}/raw_sequences.fasta",
        meta = "results/m01_acquisition/{target}/metadata.tsv",
        log_json = "results/m01_acquisition/{target}/acquisition_log.json"
    params:
        target_config = lambda wc: get_target_config(config, wc.target)
    conda:
        "envs/acquisition.yaml"
    log:
        "logs/m01_acquisition_{target}.log"
    benchmark:
        "benchmarks/m01_acquisition_{target}.txt"
    script:
        "scripts/m01_acquire_data.py"
```

> **Catatan:** `{target}` adalah wildcard Snakemake. Jika config punya 2 target (`sarscov2_indonesia` dan `flu_h1n1_local`), Snakemake otomatis menjalankan M1 dua kali secara paralel.

---

## 3. Module 2: QC & Preprocessing (M2)

### 3.1 Goal

Membuang sekuens berkualitas rendah, menghapus duplikat, dan melakukan subsampling. Output = dataset bersih siap MSA.

### 3.2 Input & Output

| Direction | File | Keterangan |
|-----------|------|------------|
| **Input** | `results/m01_acquisition/{target}/raw_sequences.fasta` | FASTA mentah dari M1 |
| **Input** | `results/m01_acquisition/{target}/metadata.tsv` | Metadata dari M1 |
| **Output** | `results/m02_qc/{target}/cleaned_sequences.fasta` | FASTA bersih siap MSA |
| **Output** | `results/m02_qc/{target}/qc_metadata.tsv` | Metadata yang sudah difilter |
| **Output** | `results/m02_qc/{target}/qc_report.json` | Statistik QC lengkap |

### 3.3 Script: `scripts/m02_qc_filter.py`

**Logic Flow (berurutan, pipeline filter):**

```text
raw_sequences.fasta (N sequences)
        │
        ▼
  ┌─ 1. Filter N-content ──────────────────────────────────────┐
  │    Hitung % basa 'N' per sekuens.                          │
  │    Buang jika > max_n_content (default: 5%)                │
  └──────────────────────────────────────── dropped: X ────────┘
        │
        ▼
  ┌─ 2. Filter Length ─────────────────────────────────────────┐
  │    Dari template: expected_length ± length_tolerance_%     │
  │    SARS-CoV-2: accept 26913–32893 bp                       │
  │    Buang jika di luar range                                │
  └──────────────────────────────────────── dropped: Y ────────┘
        │
        ▼
  ┌─ 3. Deduplication ─────────────────────────────────────────┐
  │    Hash (MD5) tiap sequence string (case-insensitive).     │
  │    Simpan 1 representasi, catat ID duplikat di metadata.   │
  └──────────────────────────────────────── dropped: Z ────────┘
        │
        ▼
  ┌─ 4. Subsampling (jika enabled) ────────────────────────────┐
  │    random: uniform sampling dengan seed (reproducible)     │
  │    temporal: N/bulan agar distribusi waktu adil             │
  │    geographic: N/negara jika multi-country                  │
  └──────────────────────────────────── sampled to: max_seq ───┘
        │
        ▼
  cleaned_sequences.fasta (M sequences, M ≤ max_sequences)
```

### 3.4 QC Report Schema (`qc_report.json`)

```json
{
  "target": "sarscov2_indonesia",
  "timestamp": "2024-03-25T09:00:00",
  "input_count": 5300,
  "filter_results": {
    "n_content_removed": 120,
    "length_removed": 45,
    "duplicates_removed": 830
  },
  "post_filter_count": 4305,
  "subsample": {
    "enabled": true,
    "strategy": "random",
    "seed": 42,
    "final_count": 1000
  },
  "output_count": 1000,
  "length_stats": {
    "min": 27105,
    "max": 30012,
    "mean": 29856,
    "median": 29880
  }
}
```

### 3.5 Snakemake Rule (Draft)

```python
rule m02_qc:
    input:
        fasta = "results/m01_acquisition/{target}/raw_sequences.fasta",
        meta = "results/m01_acquisition/{target}/metadata.tsv"
    output:
        fasta = "results/m02_qc/{target}/cleaned_sequences.fasta",
        meta = "results/m02_qc/{target}/qc_metadata.tsv",
        report = "results/m02_qc/{target}/qc_report.json"
    params:
        qc_config = config["qc"],
        template = lambda wc: get_template(config, wc.target)
    conda:
        "envs/qc.yaml"
    log:
        "logs/m02_qc_{target}.log"
    benchmark:
        "benchmarks/m02_qc_{target}.txt"
    script:
        "scripts/m02_qc_filter.py"
```

---

## 4. Execution Checklist

Jika plan Stage 1 ini disetujui, berikut langkah eksekusi yang akan dilakukan:

- [ ] Buat folder structure (`config/`, `config/templates/`, `scripts/`, `envs/`, `rules/`, `results/`, `logs/`, `benchmarks/`)
- [ ] Tulis `config/config.yaml` dan `config/templates/sars-cov-2.yaml`
- [ ] Tulis `envs/acquisition.yaml` dan `envs/qc.yaml` (conda env definitions)
- [ ] Tulis skeleton `Snakefile` + `rules/m01_acquisition.smk` + `rules/m02_qc.smk`
- [ ] Implementasi `scripts/m01_acquire_data.py`
- [ ] Implementasi `scripts/m02_qc_filter.py`
- [ ] Dry-run Snakemake (`snakemake -n`) untuk verifikasi DAG
- [ ] Test run dengan SARS-CoV-2 Indonesia (kecil, 50–100 sequences) sebagai PoC
- [ ] Inspeksi `qc_report.json` untuk validasi hasil
