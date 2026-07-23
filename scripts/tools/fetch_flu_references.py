import argparse
import hashlib
import os
import sys
from pathlib import Path
from Bio import Entrez, SeqIO

# Konfigurasi Entrez. Override via env so other labs need not edit source.
Entrez.email = os.environ.get("ENTREZ_EMAIL", "lab.biotek@ugm.ac.id")
Entrez.tool = os.environ.get("ENTREZ_TOOL", "PandemicPrep_UGM")


def validate_record(record, subtype, accession):
    """Sanity-check a downloaded reference before it can poison every downstream step.

    A retired/replaced accession or an API mix-up would otherwise write the WRONG
    genome into config/references/ and silently corrupt alignment, entropy scoring
    and every primer coordinate, with no point of detection.
    """
    problems = []

    got = (record.id or "").split(".")[0].upper()
    want = accession.split(".")[0].upper()
    if got != want:
        problems.append(f"accession mismatch: requested {accession}, received {record.id}")

    if len(record.seq) == 0:
        problems.append("empty sequence")

    # The subtype label (e.g. H3N2) should appear in the record description.
    desc = f"{record.description} {record.annotations.get('organism', '')}".lower()
    if subtype and subtype.lower() not in desc.replace(" ", ""):
        problems.append(
            f"subtype '{subtype}' not found in record description "
            f"({record.description[:80]!r}) — possible wrong reference"
        )

    return problems


def write_checksum(path: Path, subtype: str, accession: str):
    """Record a SHA256 of the written FASTA for reproducibility auditing."""
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    (path.parent / f"{path.stem}.sha256").write_text(
        f"{digest}  {path.name}\naccession: {accession}\nsubtype: {subtype}\n"
    )
    print(f"[{subtype}] SHA256 {digest[:16]}... recorded")

TARGETS = {
    "H3N2": "NC_007366.1",
    "H1N1": "NC_026433.1",
    "H5N1": "NC_007362.1"
}

def fetch_and_curate(subtype, accession, outdir):
    print(f"\n[{subtype}] Fetching GenBank record for {accession} via E-Utilities...")
    try:
        with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
            record = SeqIO.read(handle, "genbank")
    except Exception as e:
        raise SystemExit(f"[{subtype}] FAILED to fetch {accession}: {e}")

    problems = validate_record(record, subtype, accession)
    if problems:
        raise SystemExit(
            f"[{subtype}] Reference validation FAILED for {accession}:\n  - "
            + "\n  - ".join(problems)
            + "\nRefusing to write a possibly-wrong reference genome."
        )

    # 1. Tulis FASTA
    fasta_path = outdir / f"flu_a_{subtype.lower()}_ha.fasta"
    with open(fasta_path, "w") as f:
        f.write(f">{record.id} Influenza A {subtype} HA Reference\n")
        f.write(f"{str(record.seq)}\n")
    print(f"[{subtype}] Saved FASTA -> {fasta_path.name} ({len(record.seq)} bp)")
    write_checksum(fasta_path, subtype, accession)

    # 2. Kurasi GFF3 Manual dan Indah
    gff_path = outdir / f"flu_a_{subtype.lower()}_ha.gff"
    
    # Kita akan parse dari feature GenBank -> GFF3
    gff_lines = ["##gff-version 3", f"##sequence-region {record.id} 1 {len(record.seq)}"]
    
    seq_len = len(record.seq)
    cds_found = False
    cds_start, cds_end = 0, 0
    
    # Gene record — structural parent only (NO gene= attribute!)
    # Augur translate & extract_genes.py scan for gene= or Name= on CDS lines.
    # If we put gene=HA here, augur will validate this 1760bp span as coding → fail mod3.
    gff_lines.append(f"{record.id}\tNCBI\tgene\t1\t{seq_len}\t.\t+\t.\tID=gene-HA")

    for feature in record.features:
        if feature.type == "CDS":
            cds_found = True
            # GenBank feature coordinates are 0-based half-open (in Biopython)
            # GFF3 is 1-based closed
            cds_start = feature.location.start + 1
            cds_end = feature.location.end
            phase = "0"
            gff_lines.append(f"{record.id}\tNCBI\tCDS\t{cds_start}\t{cds_end}\t.\t+\t{phase}\tID=cds-HA;Parent=gene-HA;Name=HA_protein;gene=HA")
            
        elif feature.type == "mat_peptide":
            mp_start = feature.location.start + 1
            mp_end = feature.location.end
            product = feature.qualifiers.get("product", ["mat_peptide"])[0]
            # Replace spaces with underscores for GFF attributes
            product_clean = product.replace(" ", "_")
            gff_lines.append(f"{record.id}\tNCBI\tmat_peptide\t{mp_start}\t{mp_end}\t.\t+\t.\tID=mat-{product_clean};Parent=cds-HA;Name={product_clean}")

    # Kurasi 5'UTR dan 3'UTR (dihitung berdasarkan letak CDS dari ujung segmen)
    if cds_found:
        if cds_start > 1:
            gff_lines.append(f"{record.id}\tNCBI\tfive_prime_UTR\t1\t{cds_start - 1}\t.\t+\t.\tID=utr-5;Parent=gene-HA")
        if cds_end < seq_len:
            gff_lines.append(f"{record.id}\tNCBI\tthree_prime_UTR\t{cds_end + 1}\t{seq_len}\t.\t+\t.\tID=utr-3;Parent=gene-HA")

    with open(gff_path, "w") as f:
        f.write("\n".join(gff_lines) + "\n")
    print(f"[{subtype}] Saved Curated GFF3 -> {gff_path.name}")


def main():
    parser = argparse.ArgumentParser(description="Fetch and curate viral reference genomes from NCBI")
    parser.add_argument("--subtype", help="Virus subtype label (e.g. H5N1, H3N2, MPXV)")
    parser.add_argument("--accession", help="NCBI GenBank Accession ID (e.g. NC_007362.1)")
    parser.add_argument("--outdir", default="config/references", help="Output directory path")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    print("Mendownload dan mengkurasi referensi dari NCBI secara presisi menggunakan Biopython (Entrez)...")

    if args.subtype and args.accession:
        # Run specific target from CLI
        fetch_and_curate(args.subtype, args.accession, outdir)
    else:
        # Run default targets
        for subtype, accession in TARGETS.items():
            fetch_and_curate(subtype, accession, outdir)
            
    print("\nProses pengambilan dan kurasi GFF selesai!")

if __name__ == "__main__":
    main()
