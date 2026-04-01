# 🧬 Master Plan (V2): Pandemic Preparedness Pipeline — Nextstrain Hybrid Architecture

> **Project**: Pandemic Preparedness — Development of Mux RDT POCT (Orthomyxovirus & Orthocoronavirinae)
> **Mandate**: Kementerian Kesehatan + WHO, kolaborasi BRIN-UGM
> **Role**: Bioinformatics Research Assistant
> **Goal**: Fully automated in-silico pipeline generating Multiplex RT-LAMP LFA panels (w/ PCR Benchmark) for SARS-CoV-2, Flu A, dan Flu B menggunakan ekosistem Nextstrain (Augur) dan PrimerLab.

---

## 1. Architecture Overview (The Nextstrain-PrimerLab Hybrid)

Pipeline ini telah direvisi secara fundamental untuk mengadopsi standar WHO (**Nextstrain**). Beban kerja *data engineering* mentah (M1-M4) sepenuhnya didelegasikan ke *Augur toolkit*, sementara perhitungan evolusi dN/dS (HyPhy) dan desain diagnostik (PrimerLab) menjadi modul spesifik andalan kita.

| Stage | Modules | Responsibility | Tech / Tools |
|-------|---------|---------------|------|
| **Nextstrain Core** | M1–M4, M6 | Download, Spatial/Temporal QC, MSA, Tree, Clade ID | `augur` (filter, align, tree, traits) |
| **Evolution Analytics**| M5 | Selection Pressure (dN/dS) | `HyPhy` (FUBAR/MEME/GARD) |
| **Bridge** | M7 | Translate JSON/CSV → `manifest.yaml` | Custom Python |
| **PrimerLab Engine** | M8–M9 | Multiplex RDT design & validation (BLAST) | `primerlab-genomic` |
| **Visual Reporting** | M10 | Laporan HTML dengan Web Auspice bersarang | Quarto + Auspice |

---

## 2. Module Details (Revised for V2)

### The Nextstrain Ecosystem 🌍 (M1, M2, M3, M4, M6)

Kelima modul kustom yang sangat merepotkan di V1 (Script Python Download, MAFFT manual, Logika Subsampling naïf, IQ-TREE manual) **sepenuhnya dihapus dan diganti** dengan *rule* Snakemake standar dari repositori `nextstrain`.

| Nextstrain Command | Peran Modul Lama | Deskripsi Keunggulan |
|--------------------|------------------|----------------------|
| `ncov-ingest` / `augur fetch` | **M1** (Data Acquisition) | Mengambil seluruh metadata dan FASTA murni dari NCBI/GenBank secara real-time. (Tanpa butuh login GISAID). |
| `augur filter` | **M2** (QC & Subsample) | Tidak sekadar memotong ukuran/N-content, melainkan melakukan **Spatio-temporal equitable subsampling** (Misal: max 10 sekuens per negara per minggu). Ini wajib mencegah bias. |
| `augur align` | **M3** (MSA) | Mengorkestrasi MAFFT dengan sangat stabil untuk puluhan ribu data. |
| `augur tree` | **M4** (Phylogeny) | Mengorkestrasi IQ-TREE / FastTree secara otomatis + *time-resolved* refiner. |
| `augur clades` / `traits` | **M6** (Clade ID & Protein) | Mendefinisikan Clade (seperti JN.1, H1N1-victoria) dan secara instan menerjemahkan nukleotida yang sejajar menjadi petaan titik mutasi Asam Amino. Memproduksi `auspice.json`. |

---

### The Custom Analytics 🔍 (M5)

```text
Input  : alignment (dari augur align), tree (dari augur tree)
Output : dnds_results.csv (Selection Pressure Scores per kodon)
Tools  : PAML / HyPhy
```
- Nextstrain **tidak bisa** mengukur gaya seleksi Darwinian. 
- Di sinilah Modul 5 *Custom* kita tetap hidup. Ia akan memakan output setengah-jadi dari Nextstrain lalu memompa angka dN/dS (HyPhy FUBAR/MEME) agar kita bisa memprediksi mutasi pelarian genetik (*genetic escape*).

---

### The Bridge 🔀 (M7)

```text
Input  : auspice.json (Nextstrain), dnds_results.csv (HyPhy)
Output : primerlab_batch_manifest.yaml
```
M7 adalah penerjemah cerdas. Ia menarik "Daerah tak tersentuh mutasi" (dari dN/dS dan frekuensi clade *augur*) lalu menulis formulir pesanan/tugas mendetail untuk dicetak oleh PrimerLab.

---

### PrimerLab Diagnostic Engine & Validation 🎯 (M8 & M9)

```text
Input  : primerlab_batch_manifest.yaml + Reference Fasta
Output : File JSON berisi Kandidat Primer / LfA
Tools  : primerlab-genomic 
```
- **M8 (Design):** Mengeksekusi instruksi desain dari M7. Berjalan secara paralel (Batch Manifest CLI). Mampu merakit mode Multiplex RT-LAMP (6-primer fallback Core-4), RT-RPA, dan RT-qPCR benchmark.
- **M9 (Validation):** Secara otomatis menembakkan *In-Silico PCR* dan *Local BLAST* (-/+15 GB Database) untuk memastikan semua desain PrimerLab kebal terhadap Flora Pernapasan & Genom Manusia GRCh38. Kombinasi yang reaktif menyilang (Heterodimer/Off-target) dimusnahkan.

---

### The Comprehensive Report 📊 (M10)

Laporan akhir PDF/HTML via Quarto. Yang akan membuat PI-mu takjub adalah kita akan "menancapkan" / melakukan *embed* visualisasi Pohon Evolusi Interaktif khas wabah COVID (Auspice Web App) langsung berdampingan dengan tabel Reagen (Manifest Lab) dari PrimerLab.

---

## 3. Development Phases (V2 Timeline: April - Juni)

> Total Waktu Estimasi In-Silico: ~8-10 Minggu. (Lebih cepat 3 minggu berkat penghematan *Nextstrain*).

### Phase 1: PrimerLab Engine Refactoring
- [ ] Implement `manifest-run` CLI secara paralel di PrimerLab (Extensible, tanpa merusak basis kode lama).
- [ ] Integrasi profil thermodinamika parameter: RT-LAMP (Full 6 & Core 4) dan RT-RPA.
- [ ] Bangun N-Plex Heterodimer Matrix ViennaRNA *checking*.

### Phase 2: Sangkur Nextstrain (M1-M4, M6)
- [ ] Kloning *skeleton / template* rule Snakemake dari Nextstrain (e.g. `ncov` untuk SARS-CoV-2).
- [ ] Set parameter `augur filter` untuk NCBI / GenBank (Tanpa GISAID).

### Phase 3: Analytics & The Bridge (M5, M7)
- [ ] Buat rule Snakemake ekstensi untuk menyambungkan MAFFT & Tree Nextstrain ke *HyPhy*.
- [ ] Buat script Python `m07_bridge.py` pembaca kombinasi JSON+CSV -> `manifest.yaml`.

### Phase 4: Validation & Auspice Reporting (M9-M10)
- [ ] Konfigurasi BLAST Database lokal (Human GRCh38 + Virus ISPA lain).
- [ ] Merakit *Front-End* laporan Web (Jinja2/Quarto) di mana PrimerLab Lab-Manifest disajikan di samping Web Auspice.
