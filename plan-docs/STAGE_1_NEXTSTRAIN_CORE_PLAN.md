# 🌍 Stage 1: Nextstrain Core Engine (Replacing M1-M4, M6)

> **Objektif**: Merinci eksekusi *Augur Toolkit* di dalam `Snakefile` kita. Ini adalah dokumen cetak biru yang menggantikan dokumen lama `STAGE_1_PLAN.md` dan `STAGE_2_PLAN.md`.
> **Fokus**: Data Acquisition (NCBI), Spatio-temporal QC, Alignment, Fast-Phylogeny, dan Translasi Kodon.
> **Output Utama**: `auspice.json` (Tree & Metadata 100% matang) dan `aligned.fasta` (yang akan disetor ke M5 HyPhy).

---

## 1. Alur Rantai "Rules" Snakemake (Augur Chain)

Di dalam Snakemake kita nanti, modul-modul ini tidak lagi berupa script Python manual bergaris ratusan kode, melainkan pemanggilan CLI bawaan dari *Nextstrain Augur*.

### Step 1.1: Data Acquisition (M1) 📥
- **Konsep**: Kita menyediakan **DUA jalur database** (opsional di `config.yaml`: `database: "ncbi" | "gisaid"`). Ini memastikan pipeline kebal masa depan (*future-proof*) jika kamu sewaktu-waktu harus pindah atau Kemenkes mewajibkan data GISAID.
- **Rule Name**: `download_sequences` (dinamis)
- **Opsi A (NCBI GenBank/Virus)**: Menggunakan `datasets` CLI resmi milik NCBI atau modul `ncov-ingest` *open-data*. Sifatnya bebas hambatan birokrasi.
- **Opsi B (GISAID EpiCoV/EpiFlu)**: Menggunakan GISAID Provisioning CLI (`epicov`) yang dapat ditandemkan dengan Nextstrain. Sifatnya membutuhkan *Client-ID* dan kredensial akses yang disuntikkan secara aman via file `.env`.
- **Output Seragam**: Apapun opsinya, tujuannya menghasilkan file bernama `sequences.fasta` dan `metadata.tsv`!

### Step 1.2: Spatio-Temporal Filter & QC (M2) 🧹
- **Rule Name**: `augur_filter`
- **Logika Subsampling Anti-Bias**: Alih-alih merandom manual, kita me-*limit* data dengan cerdas.
- **Command Details**:
  ```bash
  augur filter \
    --sequences sequences.fasta \
    --metadata metadata.tsv \
    --min-length 29000 \                    # QC: Buang sekuens yang terpenggal (untuk SC2)
    --query "N_content < 0.05" \            # QC: Buang FASTA yang ambigunya > 5%
    --group-by country year month \         # SUBSAMPLE: Bagi keranjang berdasarkan Geografis-Waktu
    --sequences-per-group 10 \              # SUBSAMPLE: Cukup ambil 10 sekuens terbaik dari tiap keranjang
    --output filtered.fasta
  ```

### Step 1.3: Multiple Sequence Alignment (M3) 📐
- **Rule Name**: `augur_align`
- **Command Details**:
  Secara transparan memanggil MAFFT di belakang layar menggunakan CPU threading agar super cepat.
  `augur align --sequences filtered.fasta --reference-sequence ref.fasta --fill-gaps --output aligned.fasta`

### Step 1.4: Fast Phylogeny & Time-Scaling (M4) 🌳
- **Rule Name**: `augur_tree` & `augur_refine`
- **Command Details**:
  - `augur tree`: Membangun *Draft Tree* raksasa memanggil `IQ-TREE` (atau `FastTree` untuk kecepatan).
  - `augur refine`: Menskalir *tree* secara temporal (memetakan sumbu-x berdasarkan Waktu pengumpulan metadata pasien).

### Step 1.5: Mutational Traits & Cladistics (M6) 🧬
- **Rule Name**: `augur_ancestral`, `augur_translate`, `augur_clades`
- **Command Details**:
  - `ancestral`: Menjejaki nukleotida balik ke *Tree Node* moyangnya.
  - `translate`: Merubah (Coding Sequence / CDS) menjadi *Amino Acid* agar kita tahu efek mutasinya ke struktur protein (Misal: mutasi `K417N`).
  - `clades`: Melabeli ranting pohon dengan "Nama Varian" (Mendeteksi Clade WHO/Nextstrain).

### Step 1.6: The Grand Export 📦
- **Rule Name**: `augur_export`
- **Command Details**:
  Membungkus semua file Pohon, Petaan Asam Amino, dan Metadata Geografis dari Step 1.1 - 1.5 menjadi **SATU FILE** bernama `auspice_dataset.json`. 

---

## 2. Jembatan (Bridge) Emas ke Pipeline Modul Selanjutnya
Tahap *Nextstrain Core* ini secara absolut harus mengekspor 2 artefak esensial:
1. **`aligned.fasta`** (dari proses M3): Ini akan dicekik (diambil alih) oleh **Modul 5** kita secara paralel. Di M5, file ini akan dijejalkan ke **HyPhy** untuk dikonversi menjadi angka tekanan evolusi FUBAR/MEME (dN/dS).
2. **`auspice_dataset.json`** (dari proses export): Ini akan disedot oleh **Modul 7 (The Bridge Script)**. File ini berisi 100% cetak biru kekerabatan virus dan koordinat mutasinya. M7 akan mem-parsing JSON raksasa ini bersama CSV dN/dS HyPhy untuk menentukan titik *"Conserved (Dingin) vs Clade-Specific (Panas)"*.

---

## 3. Dependensi (Conda `envs/nextstrain.yaml`)
Tidak perlu pusing mengkonfigurasi biopython versi berapa; cukup gunakan environment rilis utama Nextstrain.
- `nextstrain-cli`
- `augur`
- `auspice`
- `ncbi-datasets-cli` (GenBank Fetcher)
- `mafft`
- `iqtree`

*Dengan lahirnya dokumen cetak biru pengganti ini, kita secara resmi "mematikan" (obsolete) dokumen STAGE_1_PLAN.md dan STAGE_2_PLAN.md lama, dan menjadikan dokumen ini panduan utama untuk menulis kerangka Snakefile.*
