# Panduan Lengkap HPC: Development vs Production

Panduan ini mengatur Standar Operasional Prosedur (SOP) untuk memindahkan dan menjalankan pipeline **Pandemic Preparedness** di kluster HPC institusimu. Ada dua mode eksekusi yang bisa kamu gunakan: **Mode Development (Testing)** dan **Mode Production**.

---

## 🏗️ 1. Persiapan Awal (Wajib untuk Semua Mode)

> [!WARNING]  
> **Home Directory vs Datadrive**  
> Jangan melakukan komputasi berat atau instalasi environment `conda` di folder `/home/username` karena limitnya sangat kecil. Letakkan data di direktori komputasi (`/data`, `/scratch`, atau `/work`).

### A. Clone lewat Datadrive
Masuk pakai SSH ke HPC, lalu *clone* via GitHub:
```bash
# Ganti "/data/your_username" sesuai lokasi yang ditentukan admin HPC
cd /data/your_username/
git clone https://github.com/your-username/pandemic-preparedness.git
cd pandemic-preparedness
```

### B. Isolasi Cache Snakemake (Conda)
Agar Snakemake tidak memakan disk space *home directory* kamu saat men-download tools (Auspice/Nextclade), jalankan dua baris ini sekali saja di awal:
```bash
mkdir -p /data/your_username/.conda
ln -s /data/your_username/.conda ~/.conda
```

---

## 🧪 2. Mode Development (Testing Interaktif dengan Bash)

> [!TIP]  
> Gunakan mode ini saat skrip belum 100% jadi, ketika kamu ingin *dry-run*, *troubleshooting* error, atau menjalankan *subset* kecil data (mis. 10 sequences) sambil memantau output langsung dari layar *console/bash*.

Secara teknis, kamu **tidak boleh** langsung menjalankan beban program berat di *Login Node* (halaman awal saat login ke HPC). Kamu harus "menyewa" *Compute Node* secara interaktif.

### Langkah-Langkah Development:

1. **Meminta Compute Node secara Interaktif (Salloc)**  
   ```bash
   salloc --partition=normal --cpus-per-task=4 --mem=12G --time=02:00:00 bash
   ```  
   *Penjelasan:* Kamu meminta akses interaktif untuk 4 Core, 12 GB RAM, selama maksimal 2 jam. Layar bash-mu seketika akan terlempar dari komputer "Login" ke komputer "Pekerja".

2. **Test Dry Run Snakemake**  
   Untuk memverifikasi bahwa *rules* dan integrasi GitHub telah bekerja.
   ```bash
   snakemake -np --use-conda
   ```

3. **Live Testing**  
   Jalankan eksekusi uji coba seperti lazimnya di laptop.
   ```bash
   snakemake --use-conda --cores 4
   ```

4. **Keluar dari Mode Development**  
   Saat trial & error selesai, kembalikan core yang kamu sewa.
   ```bash
   exit
   ```

---

## 🚀 3. Mode Production (Otomatis & Skala Besar via SLURM)

> [!CAUTION]  
> Gunakan mode ini jika kamu menargetkan pengolahan >1.000 data genom SARS-CoV-2/RSV, yang bisa memakan waktu pengerjaan di atas 8 jam. Kamu hanya *submit job*, tutup laptop-mu, dan biarkan HPC bekerja di balik layar.

Pada mode ini, kita tidak menggunakan layar interaktif `bash`, melainkan *menitipkan* selembar kerja instruksi (*script*.sh) kepada satpam HPC (SLURM Job Scheduler). 

### Langkah-Langkah Production:

1. **Pastikan File `slurm-run.sh` Siap di Root Folder Project**  
   ```bash
   #!/bin/bash
   #SBATCH --job-name=pandemic_run
   #SBATCH --output=logs/slurm_%j.out
   #SBATCH --error=logs/slurm_%j.err
   #SBATCH --partition=normal         
   #SBATCH --nodes=1
   #SBATCH --ntasks=1                 
   #SBATCH --cpus-per-task=16         # Alokasi Core Besar
   #SBATCH --mem=32GB                 # Alokasi RAM Besar untuk Evolusi/MSA
   #SBATCH --time=12:00:00            # Batas waktu maksimal
   
   echo "Starting Snakemake on Node: $(hostname)"
   snakemake --use-conda --cores $SLURM_CPUS_PER_TASK --keep-going --printshellcmds
   ```

2. **Daftarkan Job (Submission)**  
   ```bash
   sbatch slurm-run.sh
   ```
   *(Terminal akan merespons misalnya dengan: `Submitted batch job 100234`)*

3. **Pantau Status Job (Tanpa Takut Terputus)**  
   Melihat antrean:
   ```bash
   squeue -u your_username
   ```  
   Membaca *live output* progress pipeline (bahkan kalau laptopmu habis dimatikan):
   ```bash
   tail -f logs/slurm_100234.out
   ```
