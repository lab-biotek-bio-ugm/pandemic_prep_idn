#!/usr/bin/env python3
"""
Cleans MACSE nucleotide alignment for HyPhy compatibility.
1. Replaces MACSE frameshift characters ('!') with standard gaps ('-').
2. Replaces in-frame stop codons (TAA, TAG, TGA) with gaps ('---') because HyPhy strictly crashes on stop codons.
3. Ensures alignment length is a multiple of 3.
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(input_fasta, output_fasta):
    records = list(SeqIO.parse(input_fasta, "fasta"))
    if not records:
        return

    stop_codons = {"TAA", "TAG", "TGA"}
    cleaned_records = []

    for record in records:
        seq_str = str(record.seq).upper()
        
        # 1. Replace MACSE exclamation marks (frameshifts) with gaps
        seq_str = seq_str.replace("!", "-")
        
        # 2. Ensure length is multiple of 3 (pad with gaps if needed, though MACSE should guarantee this)
        remainder = len(seq_str) % 3
        if remainder != 0:
            seq_str += "-" * (3 - remainder)
            
        # 3. Replace in-frame stop codons
        cleaned_seq = []
        for i in range(0, len(seq_str), 3):
            codon = seq_str[i:i+3]
            # If the codon is a stop codon, mask it with gaps
            if codon in stop_codons:
                cleaned_seq.append("---")
            else:
                cleaned_seq.append(codon)
                
        cleaned_seq_str = "".join(cleaned_seq)
        
        cleaned_records.append(SeqRecord(Seq(cleaned_seq_str), id=record.id, description=""))

    SeqIO.write(cleaned_records, output_fasta, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python clean_alignment_for_hyphy.py <input.fasta> <output.fasta>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
