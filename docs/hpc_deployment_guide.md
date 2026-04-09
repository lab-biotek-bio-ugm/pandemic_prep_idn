# Panduan HPC: Pandemic Preparedness Pipeline
## Server: `bio6000.igf-sequencing.com`

---

## 🏗️ 1. Persiapan Awal (One-Time Setup)

> [!IMPORTANT]
> **Home vs Datadrive:**  
> Home (`/home/engkinandatama`) punya quota kecil. Semua data berat WAJIB disimpan di `datadrive`.  
> Dari screenshot-mu, `datadrive` sudah berupa symlink — artinya admin HPC sudah menyiapkan jalur ini untukmu.

### A. Cek lokasi datadrive
Buka terminal VS Code (Ctrl+`) lalu jalankan:
```bash
ls -la ~/datadrive
```
Catat path asli yang ditunjukkan symlink tersebut (misalnya `/data/engkinandatama` atau `/scratch/engkinandatama`).

### B. Pindahkan cache Conda ke Datadrive
Miniforge3 sudah terinstal di home, tapi cache packages conda (`~/.conda`) akan membengkak. Pindahkan:
```bash
# Pindahkan folder .conda yang sudah ada ke datadrive
mv ~/.conda ~/datadrive/.conda

# Buat symlink agar conda tetap menemukan cache-nya
ln -s ~/datadrive/.conda ~/.conda
```

### C. Pastikan Conda aktif
```bash
# Cek apakah conda sudah di PATH
conda --version

# Jika belum, aktifkan dulu
source ~/miniforge3/etc/profile.d/conda.sh
conda activate base
```

### D. Install Snakemake (jika belum ada)
```bash
conda install -c conda-forge -c bioconda snakemake -y
```

---

## 📥 2. Clone Repo Private ke Datadrive

### A. Buat Personal Access Token (PAT) di GitHub
1. Buka https://github.com/settings/tokens
2. Klik **"Generate new token (classic)"**
3. Beri nama: `HPC bio6000`
4. Centang scope: `repo` (full control)
5. Klik Generate → **Salin token-nya** (hanya muncul sekali!)

### B. Clone ke datadrive
```bash
cd ~/datadrive
git clone https://github.com/engkinandatama/pandemic-preparedness.git

# Saat diminta password, paste-kan PAT token (bukan password GitHub)
# Username: engkinandatama
# Password: ghp_xxxxxxxxxxxx (PAT token)
```

### C. Kunci permission folder (agar PI/rekan lab tidak bisa intip)
```bash
chmod 700 ~/datadrive/pandemic-preparedness
```

### D. Masuk ke project
```bash
cd ~/datadrive/pandemic-preparedness
```

---

## 🧪 3. Mode Development (Testing Interaktif)

> [!TIP]
> Gunakan mode ini saat masih building/testing pipeline.  
> Kamu dapat melihat output langsung di terminal VS Code.

### A. Minta Compute Node secara interaktif
```bash
# Minta 4 core, 12 GB RAM, maksimal 2 jam
salloc --partition=normal --cpus-per-task=4 --mem=12G --time=02:00:00

# Jika HPC-mu tidak pakai SLURM, coba PBS:
# qsub -I -l nodes=1:ppn=4,mem=12gb,walltime=02:00:00

# Jika tidak ada scheduler sama sekali (kecil kemungkinan), langsung saja
```

> [!NOTE]
> Jika perintah `salloc` error, tanyakan admin HPC kamu: *"Apakah ada job scheduler di sini? SLURM atau PBS?"*  
> Kemungkinan nama partisi juga berbeda (bukan `normal` tapi `batch`, `compute`, `short`, dll).

### B. Aktifkan conda (jika belum)
```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate base
```

### C. Dry Run (validasi DAG tanpa eksekusi)
```bash
cd ~/datadrive/pandemic-preparedness
snakemake -np --use-conda
```

### D. Live Run (eksekusi sungguhan)
```bash
snakemake --use-conda --cores 4
```

### E. Selesai testing → keluar dari compute node
```bash
exit
```

---

## 🚀 4. Mode Production (Nanti, setelah pipeline final)

> [!CAUTION]
> Gunakan mode ini HANYA jika pipeline sudah selesai dan kamu ingin memproses >1.000 genome semalaman.

### A. Submit job
```bash
cd ~/datadrive/pandemic-preparedness
sbatch slurm-run.sh
```

### B. Monitor
```bash
# Cek status antrean
squeue -u engkinandatama

# Lihat live output
tail -f logs/slurm_<JOB_ID>.out
```

---

## 🔄 5. Sinkronisasi Perubahan (Laptop ↔ HPC)

### Dari Laptop → HPC (setelah edit di laptop)
```bash
# Di laptop:
git add -A && git commit -m "update pipeline" && git push

# Di HPC:
cd ~/datadrive/pandemic-preparedness
git pull
```

### Dari HPC → Laptop (jika edit langsung di HPC via VS Code)
```bash
# Di HPC:
git add -A && git commit -m "fix from hpc testing" && git push

# Di laptop:
git pull
```

---

## 📋 Checklist Ringkas

| Langkah | Perintah | Status |
|---------|----------|--------|
| Cek datadrive | `ls -la ~/datadrive` | ⬜ |
| Pindah .conda | `mv ~/.conda ~/datadrive/.conda && ln -s ...` | ⬜ |
| Cek conda version | `conda --version` | ⬜ |
| Install snakemake | `conda install snakemake` | ⬜ |
| Clone repo | `git clone ... ~/datadrive/pandemic-preparedness` | ⬜ |
| Chmod 700 | `chmod 700 ~/datadrive/pandemic-preparedness` | ⬜ |
| Dry run | `snakemake -np --use-conda` | ⬜ |
| Live run | `snakemake --use-conda --cores 4` | ⬜ |
