# 🛠️ PrimerLab Genomic — Universal Upgrades & Enhancement Notes

> **Dokumen Induk:** Kompilasi seluruh modifikasi teknis (Fitur Baru) yang wajib ditambahkan ke *engine* `primerlab-genomic` agar 100% kompatibel dengan *Pandemic Preparedness Pipeline* (khususnya untuk target Multiplex RT-LAMP RDT LFA).
> **Prinsip Dasar:** Semua modifikasi bersifat *Additive/Extensible*. Tidak mengubah fundamental *core* PrimerLab sehingga tetap bisa dipakai pengguna umum untuk kebutuhan standar (misal `batch-run` lama tetap utuh).

---

## 1. CLI Baru: `manifest-run` (Modul 7 ➡️ Modul 8)
**Konteks:** Fitur `batch-run` saat ini hanya difungsikan untuk "Satu Aturan (Config) di-apply ke Semua Sequence" (Homogeneous). Pipeline kita membutuhkan eksekusi puluhan target berbeda dengan aturan/mode berbeda secara paralel (Heterogeneous).
**Kebutuhan:**
- Membuat *Command Line Interface* (CLI) baru: `primerlab manifest-run --yaml <batch_manifest.yaml>`.
- CLI ini bertindak sebagai "Command Center" yang mem-parsing yaml dari M7 dan mendistribusikan *job* ke *worker* Primer3.
- Menambahkan kapabilitas pembacaan referensi multi-target dari folder berbeda dalam satu kali jalan.

## 2. Spatial Constraints Injection (`SEQUENCE_INCLUDED_REGION`)
**Konteks:** PrimerLab harus mematuhi perintah mutlak dari algoritma Evolusi (dN/dS) / Modul 7.
**Kebutuhan:**
- Mendukung pemetaan (mapping) parameter `SEQUENCE_INCLUDED_REGION` langsung ke dalam *core wrapper* P3 (`primer3_wrapper.py`).
- Ini memaksa PrimerLab *HANYA* boleh menaruh primer di dalam jendela *conserved* (misal jendela 300bp), tidak boleh keluar sedikitpun.
- Tetap mempertahankan `SEQUENCE_TARGET` untuk mendesain primer spesifik-clade (wajib jepit mutasi SNP).

## 3. Penambahan Isothermal Target Modes (RDT POCT)
**Konteks:** PrimerLab sejauh ini hanya fokus pada PCR / TaqMan. Kita harus menambahkan profil desain isothermal.
**Kebutuhan:**
- **Mode 1: `rt-lamp` (Target Utama)**
  - Mendukung sintesis profil 6 primer (F3, B3, FIP [F1c+F2], BIP [B1c+B2], LF, LB).
  - Mengimplementasikan aturan *Gap/Distance* antar primer LAMP.
  - **Fallback Cerdas:** Logika khusus di mana jika algoritma gagal membongkar ruang untuk 6 primer harmonis, ia akan *"Fallback"* mendesain versi **Core-4 Primer** (tanpa LF & LB) dan memberinya *flag peringatan* (warning: Reaksi butuh 60 menit).
- **Mode 2: `rt-rpa` (Plan B Isothermal)**
  - Mendukung desain Recombinase Polymerase Amplification.
  - Berbasis probe (mirip TaqMan) namun dengan parameter yang dibalik: Tm rendah (37-42°C), Primer sangat panjang (30-35 bp), dan Probe raksasa (46-52 bp).
- **Mode 3: `rt-qpcr` (Gold Standard Benchmark)**
  - Mengubah terminologi valid dari `qpcr` menjadi `rt-qpcr` secara saintifik karena kita bekerja dengan Virus RNA (Orthomyxovirus & Coronavirus).

## 4. Multiplex Compatibility Matrix (N-Plex Engine)
**Konteks:** Membuat 3 set RT-LAMP (SARS-CoV-2, Flu A, Flu B) berarti mencampur 18 primer secara bersamaan (Multiplex LFA). Jika 1 primer salah berpasangan dengan primer virus lain, tes LFA akan *False Positive*.
**Kebutuhan:**
- Membangun kelas/modul baru: **N-Plex Heterodimer Engine**.
- Modul ini menerima *Top N Candidates* (misalnya top 3 atau top 5 dari M7) untuk setiap penyakit.
- Melakukan *Cross-Combination Looping* (Mencampur Kandidat 1 SC2 + Kandidat 2 Flu A + Kandidat 1 Flu B).
- Menghancurkan (*Reject*) kombinasi yang membentuk *cross-dimer* kuat (misal $\Delta G < -6.0 \text{ kcal/mol}$) pada 153 titik interaksi.
- Opsi untuk konfigurasi: `advanced.max_candidates` untuk dibatasi ke Top 100 agar simulasi matriks ViennaRNA tidak memakan RAM berlebihan di HPC.

## 5. Output Standarisasi untuk M9 (Validation Hook)
**Konteks:** Modul 9 (In-Silico Validation) yang berisi *Local BLAST* dan *isPCR* dikerjakan oleh **Snakemake**, bukan PrimerLab.
**Kebutuhan:**
- PrimerLab `manifest-run` di M8 harus me-return JSON berformat rapi per interaksi multiplex (thermodynamics data, koordinat primer, Fallback flag).
- JSON ini tidak akan dilihat oleh user (Dokter/Peneliti akhir), melainkan "dimakan" oleh M9 untuk diuji silang dengan *Human Genome* (GRCh38) & *Respiratory Flora* (15GB local BLAST db).

---

*Kesimpulan Arsitektur: Modifikasi di atas 100% menggunakan prinsip Independent Extension. PrimerLab tetap utuh sebagai Open-Source Engine publik, sementara fitur-fitur "Sultan" ini meningkatkan nilai jual repositorinya sebagai "The Ultimate Multiplex Assay Designer".*
