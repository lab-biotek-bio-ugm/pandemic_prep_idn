#!/usr/bin/env python3
"""
M05a: M03 Direct Slicer & Harsh Masker
Extracts gene nucleotide sequences directly from the M03 whole-genome MSA.
Bypasses MACSE by strictly enforcing reading-frame integrity directly on the global alignment:
1. Pads the extracted alignment slice to a multiple of 3.
2. Identifies frameshifting indels (codons with 1 or 2 gaps) and masks them with 'NNN'.
3. Identifies in-frame stop codons and masks them with '---'.
4. Outputs the cleaned alignment directly for HyPhy (FUBAR/MEME/GARD).
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
            
            feat_type = parts[2]
            if feat_type != "gene": continue
            
            start = int(parts[3]) - 1 # 1-based to 0-based
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]
            
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

def apply_harsh_masking(seq_str, strand):
    """
    Applies Harsh Masking to ensure reading frame and HyPhy compatibility:
    - Masks partial gaps (frameshifts) with 'NNN'
    - Masks stop codons with '---'
    Returns (cleaned_str, n_frameshift_masked, n_stop_masked).
    NOTE: seq_str must be a multiple of 3 in length.
    """
    # Reverse complement if on negative strand BEFORE checking codons
    if strand == "-":
        seq_str = str(Seq(seq_str).reverse_complement())

    stop_codons = {"TAA", "TAG", "TGA"}
    cleaned_seq = []
    n_frameshift = 0
    n_stop = 0

    for i in range(0, len(seq_str), 3):
        codon = seq_str[i:i+3]

        # Rule 1: Partial gaps (frameshift artifacts)
        if "-" in codon and codon != "---":
            cleaned_seq.append("NNN")
            n_frameshift += 1
            continue

        # Rule 2: In-frame Stop Codons
        if codon in stop_codons:
            cleaned_seq.append("---")
            n_stop += 1
            continue

        # Standard valid codon
        cleaned_seq.append(codon)

    return "".join(cleaned_seq), n_frameshift, n_stop

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--alignment", required=True, help="MSA from augur align (FASTA)")
    parser.add_argument("--gff", required=True, help="Reference GFF annotation")
    parser.add_argument("--reference-id", required=True, help="ID of the reference sequence in the MSA")
    parser.add_argument("--outdir", required=True, help="Directory to save the cleaned gene alignments")
    args = parser.parse_args()

    print(f"Parsing GFF: {args.gff}")
    genes = parse_gff(args.gff)
    if not genes:
        sys.exit("ERROR: No genes found in GFF")
        
    print(f"Loading MSA: {args.alignment}")
    records = list(SeqIO.parse(args.alignment, "fasta"))
    print(f"Loaded {len(records)} sequences.")

    ref_record = next((r for r in records if r.id == args.reference_id or r.id.split('.')[0] == args.reference_id), None)
    if not ref_record:
        sys.exit(f"ERROR: Reference '{args.reference_id}' not found in MSA.")
        
    print(f"Mapping coordinates against reference: {ref_record.id}")
    col_map = build_coordinate_map(str(ref_record.seq))

    os.makedirs(args.outdir, exist_ok=True)

    for gene_name, (start, end, strand) in genes.items():
        if start >= len(col_map):
            continue

        if end > len(col_map):
            end = len(col_map)

        msa_start = col_map[start]
        # In python slicing, end is exclusive. We want to include the base at `end-1`.
        msa_end = col_map[end - 1] + 1

        slice_width = msa_end - msa_start
        pad_len = (3 - (slice_width % 3)) % 3
        n_codons = (slice_width + pad_len) // 3

        out_records = []
        total_frameshift = 0
        total_stop = 0
        skipped_gappy = 0

        for r in records:
            sliced_msa = str(r.seq)[msa_start:msa_end].upper()
            sliced_msa += "N" * pad_len

            # Filter: skip sequences that are >90% gaps in this gene region
            gap_frac = sliced_msa.count("-") / len(sliced_msa)
            if gap_frac > 0.9:
                skipped_gappy += 1
                continue

            cleaned_seq_str, n_fs, n_st = apply_harsh_masking(sliced_msa, strand)
            total_frameshift += n_fs
            total_stop += n_st
            out_records.append(SeqRecord(Seq(cleaned_seq_str), id=r.id, description=""))

        out_path = os.path.join(args.outdir, f"cleaned_{gene_name}_NT.fasta")
        SeqIO.write(out_records, out_path, "fasta")

        # QC statistics
        print(f"[{gene_name}] {len(out_records)} seqs written, "
              f"{skipped_gappy} skipped (>90%% gap), "
              f"{n_codons} codons/seq, "
              f"masked: {total_frameshift} frameshifts + {total_stop} stop codons -> {out_path}")

if __name__ == "__main__":
    main()
