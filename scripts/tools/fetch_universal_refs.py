import sys
import argparse
from pathlib import Path
from Bio import Entrez, SeqIO

Entrez.email = "lab.biotek@ugm.ac.id"
Entrez.tool = "PandemicPrep_UGM"

TARGETS = {
    "MP": "NC_002016.1",
    "NP": "NC_002019.1"
}

def fetch_and_curate(segment, accession, outdir):
    print(f"\n[{segment}] Fetching GenBank record for {accession} via E-Utilities...")
    try:
        with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
            record = SeqIO.read(handle, "genbank")
    except Exception as e:
        print(f"[{segment}] FAILED to fetch {accession}: {e}")
        return

    # 1. Tulis FASTA
    fasta_path = outdir / f"flu_a_{segment.lower()}_reference.fasta"
    with open(fasta_path, "w") as f:
        f.write(f">{record.id} Influenza A {segment} Reference\n")
        f.write(f"{str(record.seq)}\n")
    print(f"[{segment}] Saved FASTA -> {fasta_path.name}")

    # 2. Kurasi GFF3 Manual dan Indah
    gff_path = outdir / f"flu_a_{segment.lower()}_annotation.gff"
    
    gff_lines = ["##gff-version 3", f"##sequence-region {record.id} 1 {len(record.seq)}"]
    seq_len = len(record.seq)
    
    # Create structural parent without gene= to prevent Augur modulo-3 errors on UTRs
    gff_lines.append(f"{record.id}\tNCBI\tgene\t1\t{seq_len}\t.\t+\t.\tID=gene-{segment}")

    min_cds_start = seq_len
    max_cds_end = 0
    cds_found = False

    for feature in record.features:
        if feature.type == "CDS":
            gene_name = feature.qualifiers.get("gene", ["unknown"])[0]
            
            for part in feature.location.parts:
                cds_found = True
                part_start = int(part.start) + 1
                part_end = int(part.end)
                
                # Update global min/max for UTRs
                min_cds_start = min(min_cds_start, part_start)
                max_cds_end = max(max_cds_end, part_end)
                
                phase = "0"
                # Multi-part CDS (like M2) must share the SAME ID for augur translate to stitch them
                gff_lines.append(f"{record.id}\tNCBI\tCDS\t{part_start}\t{part_end}\t.\t+\t{phase}\tID=cds-{gene_name};Parent=gene-{segment};Name={gene_name}_protein;gene={gene_name}")

    if cds_found:
        if min_cds_start > 1:
            gff_lines.append(f"{record.id}\tNCBI\tfive_prime_UTR\t1\t{min_cds_start - 1}\t.\t+\t.\tID=utr-5;Parent=gene-{segment}")
        if max_cds_end < seq_len:
            gff_lines.append(f"{record.id}\tNCBI\tthree_prime_UTR\t{max_cds_end + 1}\t{seq_len}\t.\t+\t.\tID=utr-3;Parent=gene-{segment}")

    with open(gff_path, "w") as f:
        f.write("\n".join(gff_lines) + "\n")
    print(f"[{segment}] Saved Curated GFF3 -> {gff_path.name}")

def main():
    outdir = Path("config/references")
    outdir.mkdir(parents=True, exist_ok=True)
    
    print("Mendownload dan mengkurasi referensi MP dan NP dengan dukungan Splicing (M1 & M2)...")
    for segment, accession in TARGETS.items():
        fetch_and_curate(segment, accession, outdir)
    print("\nProses selesai!")

if __name__ == "__main__":
    main()
