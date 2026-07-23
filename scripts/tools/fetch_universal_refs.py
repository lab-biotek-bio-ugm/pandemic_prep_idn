import argparse
import os
import sys
from pathlib import Path
from Bio import Entrez, SeqIO

# Shared reference validation / checksum helpers (single implementation).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fetch_flu_references import validate_record, write_checksum

Entrez.email = os.environ.get("ENTREZ_EMAIL", "lab.biotek@ugm.ac.id")
Entrez.tool = os.environ.get("ENTREZ_TOOL", "PandemicPrep_UGM")

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
        raise SystemExit(f"[{segment}] FAILED to fetch {accession}: {e}")

    # Segment labels (MP/NP) are not subtype strings, so only the accession and
    # non-empty checks apply here.
    problems = validate_record(record, None, accession)
    if problems:
        raise SystemExit(
            f"[{segment}] Reference validation FAILED for {accession}:\n  - "
            + "\n  - ".join(problems)
            + "\nRefusing to write a possibly-wrong reference genome."
        )

    # 1. Tulis FASTA
    fasta_path = outdir / f"flu_a_{segment.lower()}_reference.fasta"
    with open(fasta_path, "w") as f:
        f.write(f">{record.id} Influenza A {segment} Reference\n")
        f.write(f"{str(record.seq)}\n")
    print(f"[{segment}] Saved FASTA -> {fasta_path.name} ({len(record.seq)} bp)")
    write_checksum(fasta_path, segment, accession)

    # 2. Kurasi GFF3 Manual dan Indah
    gff_path = outdir / f"flu_a_{segment.lower()}_annotation.gff"
    
    gff_lines = ["##gff-version 3", f"##sequence-region {record.id} 1 {len(record.seq)}"]
    seq_len = len(record.seq)
    
    # Structural parent ID — use unique suffix to avoid colliding with per-gene IDs
    root_id = f"segment-{segment}"
    gff_lines.append(f"{record.id}\tNCBI\tgene\t1\t{seq_len}\t.\t+\t.\tID={root_id}")

    min_cds_start = seq_len
    max_cds_end = 0
    cds_found = False

    # Group CDS parts by gene name
    genes_data = {}

    for feature in record.features:
        if feature.type == "CDS":
            gene_name = feature.qualifiers.get("gene", ["unknown"])[0]
            if gene_name not in genes_data:
                genes_data[gene_name] = []

            for i, part in enumerate(feature.location.parts):
                part_start = int(part.start) + 1
                part_end = int(part.end)
                genes_data[gene_name].append({
                    "start": part_start,
                    "end": part_end,
                    "part_idx": i + 1
                })

    # Write gene and CDS features per detected gene.
    # Strategy:
    #   - Non-spliced genes (1 CDS part): emit a 'gene' feature so augur translate
    #     can detect and translate them → they appear in auspice (AA mutations).
    #   - Spliced genes (>1 CDS part, e.g. M2): do NOT emit a 'gene' feature because
    #     the gene span includes the intron and is NOT a multiple of 3, which causes
    #     augur translate to crash. For these, we only emit CDS features.
    #     m07_bridge_builder.py's CDS-fallback auto-detection will synthesize the
    #     gene region from the CDS entries, so primer design still works.
    for gene_name, parts in genes_data.items():
        gene_start = min(p["start"] for p in parts)
        gene_end = max(p["end"] for p in parts)
        is_spliced = len(parts) > 1
        gene_span = gene_end - gene_start + 1

        if not is_spliced and gene_span % 3 == 0:
            # Non-spliced gene: safe for augur translate
            gene_id = f"gene-{gene_name}"
            gff_lines.append(f"{record.id}\tNCBI\tgene\t{gene_start}\t{gene_end}\t.\t+\t.\tID={gene_id};gene={gene_name}")
            parent_id = gene_id
            print(f"  [{segment}] Gene '{gene_name}': {gene_span} bp (single exon, multiple of 3) -> gene feature added")
        else:
            # Spliced gene or non-multiple-of-3: skip gene feature, CDS-only
            parent_id = root_id
            reason = "spliced" if is_spliced else f"{gene_span} bp not multiple of 3"
            print(f"  [{segment}] Gene '{gene_name}': {gene_span} bp ({reason}) -> CDS-only (m07 fallback)")

        for p in parts:
            cds_id = f"cds-{gene_name}_p{p['part_idx']}" if is_spliced else f"cds-{gene_name}"
            gff_lines.append(f"{record.id}\tNCBI\tCDS\t{p['start']}\t{p['end']}\t.\t+\t0\tID={cds_id};Parent={parent_id};gene={gene_name}")

            min_cds_start = min(min_cds_start, p["start"])
            max_cds_end = max(max_cds_end, p["end"])
            cds_found = True

    if cds_found:
        if min_cds_start > 1:
            gff_lines.append(f"{record.id}\tNCBI\tfive_prime_UTR\t1\t{min_cds_start - 1}\t.\t+\t.\tID=utr-5;Parent={root_id}")
        if max_cds_end < seq_len:
            gff_lines.append(f"{record.id}\tNCBI\tthree_prime_UTR\t{max_cds_end + 1}\t{seq_len}\t.\t+\t.\tID=utr-3;Parent={root_id}")

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
