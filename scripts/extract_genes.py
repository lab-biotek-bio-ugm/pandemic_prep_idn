#!/usr/bin/env python3
"""
M05a: Gene Extractor
Extracts unaligned gene nucleotide sequences from a whole-genome MSA based on a reference GFF.

Workflow:
1. Load multiple sequence alignment (MSA).
2. Locate the reference sequence within the MSA.
3. Map the reference sequence's physical coordinates to the MSA column indices.
4. For each coding sequence (CDS/gene) in the GFF:
     a) Identify MSA column start and end.
     b) Slice the MSA for all sequences in those columns.
     c) Strip gap characters ('-') to produce unaligned, gene-specific FASTA variants.
5. Save as separate FASTA files ready for Codon-Aware Alignment (MACSE).
"""

import argparse
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_gff(gff_path):
    """Parse GFF to return a dict of {gene_name: (start, end, strand)} in 0-based coords."""
    genes = {}
    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 9: continue
            
            # Only process 'gene' feature type to avoid duplicate CDS entries
            feat_type = parts[2]
            if feat_type != "gene": continue
            
            start = int(parts[3]) - 1 # 1-based to 0-based
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]
            
            # Extract gene name strictly from 'gene=' attribute only
            name = None
            for attr in attributes.split(";"):
                attr = attr.strip()
                if attr.startswith("gene="):
                    name = attr.split("=")[1].strip('"')
                    break
            
            if name and name not in genes:
                genes[name] = (start, end, strand)
    return genes

def build_coordinate_map(ref_seq_str):
    """
    Returns a list `col_map` where `col_map[i]` is the MSA column index
    corresponding to the i-th ungapped nucleotide of the reference.
    """
    col_map = []
    for col_idx, char in enumerate(ref_seq_str):
        if char != "-":
            col_map.append(col_idx)
    return col_map

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--alignment", required=True, help="MSA from augur align (FASTA)")
    parser.add_argument("--gff", required=True, help="Reference GFF annotation")
    parser.add_argument("--reference-id", required=True, help="ID of the reference sequence in the MSA (e.g., 'MN908947.3')")
    parser.add_argument("--outdir", required=True, help="Directory to save the extracted unaligned gene FASTAs")
    args = parser.parse_args()

    # 1. Parse GFF
    print(f"Parsing GFF: {args.gff}")
    genes = parse_gff(args.gff)
    if not genes:
        print("ERROR: No genes found in GFF", file=sys.stderr)
        sys.exit(1)
    print(f"Found {len(genes)} genes: {list(genes.keys())}")

    # 2. Load FASTA
    print(f"Loading multiple sequence alignment: {args.alignment}")
    records = list(SeqIO.parse(args.alignment, "fasta"))
    print(f"Loaded {len(records)} sequences.")

    # 3. Find reference and map coordinates
    ref_record = None
    for r in records:
        if args.reference_id in r.id:
            ref_record = r
            break
            
    if not ref_record:
        print(f"ERROR: Reference ID '{args.reference_id}' not found in alignment.", file=sys.stderr)
        sys.exit(1)
        
    print(f"Mapping coordinates against reference: {ref_record.id}")
    col_map = build_coordinate_map(str(ref_record.seq))

    os.makedirs(args.outdir, exist_ok=True)

    # 4 & 5. Extract and save unaligned genes
    for gene_name, (start, end, strand) in genes.items():
        if end > len(col_map):
            print(f"WARNING: Gene {gene_name} exceeds reference length. Skipping.")
            continue
            
        msa_start = col_map[start]
        # map 'end' index (exclusive). If end is exactly len(col_map), it maps to the very end of MSA.
        msa_end = col_map[end-1] + 1 if end <= len(col_map) else len(str(ref_record.seq))
        
        out_records = []
        for r in records:
            # We don't want the reference included in MACSE (or we do, as a frameshift anchor!)
            # But usually it's fine.
            sliced_msa = str(r.seq)[msa_start:msa_end]
            
            # Remove gaps to get unaligned raw sequences
            unaligned_seq = sliced_msa.replace("-", "")
            
            # Reverse complement if on negative strand
            if strand == "-":
                unaligned_seq = str(Seq(unaligned_seq).reverse_complement())
                
            if len(unaligned_seq) == 0:
                continue # Sequence has no data for this gene
                
            out_records.append(SeqRecord(Seq(unaligned_seq), id=r.id, description=""))
            
        out_path = os.path.join(args.outdir, f"unaligned_{gene_name}.fasta")
        SeqIO.write(out_records, out_path, "fasta")
        print(f"Extracted {len(out_records)} unaligned sequences for {gene_name} -> {out_path}")

if __name__ == "__main__":
    main()
