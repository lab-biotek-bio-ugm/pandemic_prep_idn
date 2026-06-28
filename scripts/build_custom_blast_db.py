"""
Build Custom BLAST Database for Respiratory Pathogen Specificity Verification
==============================================================================
Script ini membaca konfigurasi `blast_verification.custom_db_build` dari file
YAML Snakemake dan membangun dua database BLAST lokal:
  - resources/blast_db/human_db  : GRCh38 + T2T-CHM13
  - resources/blast_db/virus_db  : Respiratory viruses (dari NCBI Datasets CLI)

Fitur:
  - Hash-based cache: skip download jika config tidak berubah.
  - Configurable max_seqs_per_taxon untuk mengontrol ukuran database.
  - Validasi ketersediaan tools (makeblastdb, datasets).
  - Fallback error yang informatif jika tools tidak tersedia.
"""

import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

import yaml


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def log(msg: str):
    print(f"[build_db] {msg}", flush=True)


def error(msg: str):
    print(f"[build_db] ❌ ERROR: {msg}", file=sys.stderr, flush=True)


def compute_config_hash(db_config: dict) -> str:
    """Compute SHA256 hash from the db config dict for cache invalidation."""
    canonical = json.dumps(db_config, sort_keys=True, ensure_ascii=True)
    return hashlib.sha256(canonical.encode()).hexdigest()


def check_tool(name: str) -> str:
    """Return path to tool or exit with helpful message if missing."""
    path = shutil.which(name)
    if not path:
        error(
            f"Tool '{name}' tidak ditemukan di PATH.\n"
            f"  → Untuk makeblastdb: conda install -c bioconda blast\n"
            f"  → Untuk datasets:    conda install -c conda-forge ncbi-datasets-cli"
        )
        sys.exit(1)
    return path


def download_file(url: str, dest: Path) -> Path:
    """Download a file from URL to dest, skip if already exists."""
    if dest.exists():
        log(f"  ✓ File sudah ada, skip download: {dest.name}")
        return dest

    log(f"  ⬇ Mengunduh: {url}")
    log(f"     → Tujuan: {dest}")
    try:
        with urllib.request.urlopen(url, timeout=600) as resp, open(dest, "wb") as out:
            total = int(resp.headers.get("Content-Length", 0))
            downloaded = 0
            chunk_size = 1024 * 1024  # 1 MB
            while True:
                chunk = resp.read(chunk_size)
                if not chunk:
                    break
                out.write(chunk)
                downloaded += len(chunk)
                if total > 0:
                    pct = downloaded / total * 100
                    print(f"\r     {pct:.1f}% ({downloaded/1e6:.0f}/{total/1e6:.0f} MB)", end="", flush=True)
        print()
        log(f"  ✓ Download selesai: {dest.name}")
    except Exception as exc:
        error(f"Gagal mengunduh {url}: {exc}")
        dest.unlink(missing_ok=True)
        sys.exit(1)
    return dest


def extract_gz(gz_path: Path, out_path: Path) -> Path:
    """Extract a .gz file to out_path, skip if already exists."""
    if out_path.exists():
        log(f"  ✓ File sudah diekstrak, skip: {out_path.name}")
        return out_path

    log(f"  📦 Mengekstrak: {gz_path.name} → {out_path.name}")
    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    log(f"  ✓ Ekstrak selesai: {out_path.name}")
    return out_path


def filter_fasta(in_path: Path, out_path: Path, min_len: int) -> int:
    """
    Copy FASTA from in_path to out_path, skipping sequences shorter than min_len.
    Returns number of sequences written.
    """
    if out_path.exists():
        log(f"  ✓ File filtered sudah ada, skip: {out_path.name}")
        # Count sequences for reporting
        count = sum(1 for line in open(out_path) if line.startswith(">"))
        return count

    log(f"  🔎 Filtering sekuens < {min_len} bp dari {in_path.name}...")
    written = 0
    skipped = 0
    header = None
    seq_chunks = []

    with open(in_path) as f_in, open(out_path, "w") as f_out:
        for line in f_in:
            line = line.rstrip()
            if line.startswith(">"):
                # Flush previous
                if header and seq_chunks:
                    seq = "".join(seq_chunks)
                    if len(seq) >= min_len:
                        f_out.write(header + "\n" + seq + "\n")
                        written += 1
                    else:
                        skipped += 1
                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line)
        # Last record
        if header and seq_chunks:
            seq = "".join(seq_chunks)
            if len(seq) >= min_len:
                f_out.write(header + "\n" + seq + "\n")
                written += 1
            else:
                skipped += 1

    log(f"  ✓ Ditulis: {written} sekuens, dilewati (terlalu pendek): {skipped}")
    return written


def subsample_fasta(in_path: Path, out_path: Path, max_seqs: int) -> int:
    """
    Subsample FASTA to at most max_seqs sequences (take first N).
    Returns number of sequences written.
    """
    if out_path.exists():
        log(f"  ✓ File subsampled sudah ada, skip: {out_path.name}")
        count = sum(1 for line in open(out_path) if line.startswith(">"))
        return count

    log(f"  ✂ Subsample ke maks {max_seqs} sekuens dari {in_path.name}...")
    written = 0
    in_seq = False

    with open(in_path) as f_in, open(out_path, "w") as f_out:
        for line in f_in:
            if line.startswith(">"):
                if written >= max_seqs:
                    break
                written += 1
                in_seq = True
                f_out.write(line)
            elif in_seq:
                f_out.write(line)

    log(f"  ✓ Subsampled: {written} sekuens")
    return written


def make_blast_db(makeblastdb_path: str, fasta_path: Path, db_out_prefix: Path, db_type: str = "nucl"):
    """Run makeblastdb to create BLAST index from FASTA file."""
    log(f"  🔨 Membangun BLAST database: {db_out_prefix.name} ...")
    cmd = [
        makeblastdb_path,
        "-in", str(fasta_path),
        "-dbtype", db_type,
        "-out", str(db_out_prefix),
        "-parse_seqids",
        "-title", db_out_prefix.name
    ]
    log(f"     Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        error(f"makeblastdb gagal:\n{result.stderr}")
        sys.exit(1)
    log(f"  ✓ Database siap: {db_out_prefix}")


# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------

def cleanup_intermediates(output_dir: Path, hosts: list, virus_taxons: list):
    """
    Remove intermediate files that are no longer needed after BLAST DB is built.

    Removed:
      - Downloaded genome .gz files (large, can be re-downloaded if needed)
      - Extracted individual .fna genome files
      - Individual per-taxon FASTA files (taxon_*.fasta)
      - human_combined_raw.fasta  (pre-filter combined)
      - virus_combined_raw.fasta  (pre-filter combined)
      - virus_combined_filtered.fasta (post-filter, pre-subsample)

    Kept (needed by BLAST at query time):
      - human_combined_filtered.fasta → makeblastdb input for human_db
      - virus_combined_subsampled.fasta → makeblastdb input for virus_db
      - human_db.* / virus_db.* (BLAST index files)
      - build_complete.txt (cache marker)
    """
    log("\n🧹 Membersihkan file intermediate...")
    removed = 0
    freed_bytes = 0

    def try_remove(path: Path):
        nonlocal removed, freed_bytes
        if path.exists():
            size = path.stat().st_size
            path.unlink()
            freed_bytes += size
            removed += 1
            log(f"   🗑 Dihapus: {path.name} ({size/1e6:.1f} MB)")
        else:
            log(f"   ✓ Sudah tidak ada: {path.name}")

    # Remove downloaded .gz and extracted .fna files per host
    for host in hosts:
        name = host.get("name", "unknown")
        try_remove(output_dir / f"{name}_genome.fna.gz")
        try_remove(output_dir / f"{name}_genome.fna")

    # Remove per-taxon individual FASTA files
    for taxon in virus_taxons:
        name = taxon.get("name", "unknown")
        taxid = taxon.get("taxid", "")
        try_remove(output_dir / f"taxon_{taxid}_{name.replace(' ', '_')}.fasta")

    # Remove combined raw and intermediate FASTA files
    try_remove(output_dir / "human_combined_raw.fasta")
    try_remove(output_dir / "virus_combined_raw.fasta")
    try_remove(output_dir / "virus_combined_filtered.fasta")

    log(f"   ✅ Selesai: {removed} file dihapus, dibebaskan ~{freed_bytes/1e9:.2f} GB")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Build custom respiratory BLAST database for primer specificity verification."
    )
    parser.add_argument(
        "--config-file", required=True,
        help="Path ke file YAML konfigurasi Snakemake (config/config_influenza_a_mp.yaml)"
    )
    args = parser.parse_args()

    # Load full config
    config_path = Path(args.config_file)
    if not config_path.exists():
        error(f"File konfigurasi tidak ditemukan: {config_path}")
        sys.exit(1)

    with open(config_path) as f:
        config = yaml.safe_load(f)

    blast_config = config.get("blast_verification", {})
    db_build = blast_config.get("custom_db_build", {})

    if not db_build.get("enabled", False):
        log("custom_db_build.enabled = false, tidak ada yang dilakukan.")
        sys.exit(0)

    output_dir = Path(db_build.get("output_dir", "resources/blast_db"))
    output_dir.mkdir(parents=True, exist_ok=True)

    max_seqs = db_build.get("max_seqs_per_taxon", 500)
    min_len = db_build.get("seq_length_min", 500)
    do_cleanup = db_build.get("cleanup_intermediates", True)  # Default: hapus file sampah

    hosts = db_build.get("hosts", [])
    virus_taxons = db_build.get("virus_taxons", [])

    if not hosts and not virus_taxons:
        error("Tidak ada host atau virus taxon yang dikonfigurasi di custom_db_build.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Cache invalidation — compare config hash
    # ------------------------------------------------------------------
    marker_file = output_dir / "build_complete.txt"
    config_hash = compute_config_hash(db_build)

    if marker_file.exists():
        stored = marker_file.read_text().strip().split("\n")
        stored_hash = stored[0].replace("config_hash:", "").strip() if stored else ""
        if stored_hash == config_hash:
            log("✅ Database sudah up-to-date (config hash cocok). Skip build.")
            sys.exit(0)
        else:
            log("⚠ Config berubah (hash tidak cocok). Membangun ulang database...")

    # ------------------------------------------------------------------
    # Check required tools
    # ------------------------------------------------------------------
    makeblastdb_path = check_tool("makeblastdb")
    datasets_path = check_tool("datasets")
    log(f"✓ makeblastdb: {makeblastdb_path}")
    log(f"✓ datasets CLI: {datasets_path}")

    # ------------------------------------------------------------------
    # Build Human Database (human_db)
    # ------------------------------------------------------------------
    log("\n=== Membangun Human Database (human_db) ===")
    human_combined_fasta = output_dir / "human_combined_raw.fasta"
    human_filtered_fasta = output_dir / "human_combined_filtered.fasta"

    if not human_combined_fasta.exists():
        with open(human_combined_fasta, "w") as out_f:
            for host in hosts:
                name = host.get("name", "unknown")
                if not host.get("enabled", True):
                    log(f"  ✓ Host '{name}' dinonaktifkan di konfigurasi (enabled: false), dilewati.")
                    continue
                url = host.get("fasta_url", "")
                if not url:
                    log(f"  ⚠ Host '{name}' tidak memiliki fasta_url, dilewati.")
                    continue

                log(f"\n  → Memproses host: {name}")
                gz_path = output_dir / f"{name}_genome.fna.gz"
                fna_path = output_dir / f"{name}_genome.fna"

                download_file(url, gz_path)
                extract_gz(gz_path, fna_path)

                # Append to combined (no length filter for chromosomes)
                log(f"  📎 Menggabungkan {name} ke human_combined_raw.fasta ...")
                with open(fna_path) as f:
                    shutil.copyfileobj(f, out_f)
    else:
        log("  ✓ human_combined_raw.fasta sudah ada, skip download.")

    # Filter short contigs/scaffolds from human genome
    filter_fasta(human_combined_fasta, human_filtered_fasta, min_len=1000)

    # Build BLAST DB
    human_db_prefix = output_dir / "human_db"
    make_blast_db(makeblastdb_path, human_filtered_fasta, human_db_prefix)

    # ------------------------------------------------------------------
    # Build Virus Database (virus_db)
    # ------------------------------------------------------------------
    log("\n=== Membangun Virus Database (virus_db) ===")
    virus_combined_fasta = output_dir / "virus_combined_raw.fasta"
    virus_filtered_fasta = output_dir / "virus_combined_filtered.fasta"
    virus_subsampled_fasta = output_dir / "virus_combined_subsampled.fasta"

    if not virus_combined_fasta.exists():
        with open(virus_combined_fasta, "w") as out_f:
            for taxon in virus_taxons:
                name = taxon.get("name", "unknown")
                taxid = taxon.get("taxid")
                if not taxid:
                    log(f"  ⚠ Taxon '{name}' tidak memiliki taxid, dilewati.")
                    continue

                log(f"\n  → Mengunduh: {name} (taxid: {taxid})")
                taxon_fasta = output_dir / f"taxon_{taxid}_{name.replace(' ', '_')}.fasta"

                if taxon_fasta.exists():
                    log(f"  ✓ FASTA sudah ada, skip download: {taxon_fasta.name}")
                else:
                    # Use NCBI datasets CLI to download nucleotide sequences by taxon
                    with tempfile.TemporaryDirectory() as tmpdir:
                        zip_path = Path(tmpdir) / "ncbi_dataset.zip"
                        cmd = [
                            datasets_path, "download", "virus", "genome",
                            "taxon", str(taxid),
                            "--filename", str(zip_path),
                            "--no-progressbar"
                        ]
                        log(f"     Command: {' '.join(cmd)}")
                        result = subprocess.run(cmd, capture_output=True, text=True)
                        if result.returncode != 0:
                            log(f"  ⚠ Gagal mengunduh taxon {taxid} ({name}): {result.stderr[:200]}")
                            log(f"     Melanjutkan ke taxon berikutnya...")
                            continue

                        # Unzip and find FASTA
                        unzip_dir = Path(tmpdir) / "unzipped"
                        import zipfile
                        with zipfile.ZipFile(zip_path) as zf:
                            zf.extractall(unzip_dir)

                        # Find all .fna files in the extracted directory
                        fna_files = list(unzip_dir.rglob("*.fna"))
                        if not fna_files:
                            log(f"  ⚠ Tidak ada file FASTA ditemukan untuk taxon {taxid}")
                            continue

                        with open(taxon_fasta, "w") as tf:
                            for fna in fna_files:
                                with open(fna) as ff:
                                    shutil.copyfileobj(ff, tf)

                log(f"  📎 Menggabungkan {name} ke virus_combined_raw.fasta ...")
                with open(taxon_fasta) as f:
                    shutil.copyfileobj(f, out_f)
    else:
        log("  ✓ virus_combined_raw.fasta sudah ada, skip download.")

    # Filter & subsample virus sequences
    filter_fasta(virus_combined_fasta, virus_filtered_fasta, min_len=min_len)
    subsample_fasta(virus_filtered_fasta, virus_subsampled_fasta, max_seqs=max_seqs * len(virus_taxons))

    # Build BLAST DB
    virus_db_prefix = output_dir / "virus_db"
    make_blast_db(makeblastdb_path, virus_subsampled_fasta, virus_db_prefix)

    # ------------------------------------------------------------------
    # Write marker file with config hash
    # ------------------------------------------------------------------
    marker_file.write_text(
        f"config_hash:{config_hash}\n"
        f"human_db:{human_db_prefix}\n"
        f"virus_db:{virus_db_prefix}\n"
    )
    log(f"\n✅ Semua database berhasil dibangun!")
    log(f"   Human DB : {human_db_prefix}")
    log(f"   Virus  DB: {virus_db_prefix}")
    log(f"   Marker   : {marker_file}")

    # ------------------------------------------------------------------
    # Cleanup intermediate files (optional, default: enabled)
    # ------------------------------------------------------------------
    if do_cleanup:
        cleanup_intermediates(output_dir, hosts, virus_taxons)
    else:
        log("\nℹ cleanup_intermediates = false, file intermediate TIDAK dihapus.")
        log("  Untuk menghemat storage, set cleanup_intermediates: true di config.")


if __name__ == "__main__":
    main()
