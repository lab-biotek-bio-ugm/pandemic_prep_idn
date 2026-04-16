#!/usr/bin/env python3
"""
M05a: M03 Smart Slicer & Harsh Masker

Extracts gene nucleotide sequences directly from the M03 whole-genome MSA.
Bypasses MACSE by strictly enforcing reading-frame integrity directly on the
global alignment. Supports both:

  - Single-CDS genes    (SARS-CoV-2, Flu M1, HA, NP, etc.)
  - Multi-CDS / spliced genes  (Flu M2, NS2/NEP) — segments are slice-and-
    stitched before masking to reproduce the mature mRNA reading frame.

Processing per gene:
  1. Collect all CDS segments (sorted by start position).
  2. Map each segment to MSA columns using the reference coordinate map.
  3. Concatenate all segments (slice-and-stitch = exon joining).
  4. Pad concatenated sequence to a multiple of 3.
  5. Apply Harsh Masking:
       - Partial gap codons (frameshifts)  → 'NNN'
       - In-frame stop codons              → '---'
  6. Filter sequences with >90% gaps in the gene region.
  7. Output cleaned FASTA for HyPhy (FUBAR/MEME/GARD).
"""

import argparse
import sys
import os
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_gff(gff_path):
    """Parse GFF3 to return CDS segments grouped by gene name.

    Returns:
        dict: {gene_name: [(start0, end, strand), ...]}
              Coordinates are 0-based half-open [start, end).
              Each tuple is one CDS segment. Sorted by start position.
    """
    gene_segments = defaultdict(list)

    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            feat_type = parts[2]
            if feat_type not in ("CDS", "gene"):
                continue

            start = int(parts[3]) - 1  # 1-based → 0-based
            end = int(parts[4])         # end stays as-is (exclusive in 0-based)
            strand = parts[6]
            attributes = parts[8]

            # Coba ekstrak berbagai variasi atribut penamaan secara hierarkis
            name = None
            for key in ["gene=", "Name=", "locus_tag="]:
                for attr in attributes.split(";"):
                    attr = attr.strip()
                    if attr.startswith(key):
                        name = attr.split("=", 1)[1].strip('"\n\r')
                        break
                if name:
                    break

            if not name:
                continue

            if feat_type == "CDS":
                gene_segments[name].append((start, end, strand))
            elif feat_type == "gene":
                # Register gene with placeholder — ensures entry exists
                # even for genes that only have 'gene' feature (legacy GFFs)
                if name not in gene_segments:
                    gene_segments[name] = [(start, end, strand)]

    # Sort each gene's segments by start position
    for name in gene_segments:
        gene_segments[name].sort(key=lambda x: x[0])

    return dict(gene_segments)


def build_coordinate_map(ref_seq_str):
    """Build MSA column → reference position map.

    Returns:
        list: col_map[i] is the MSA column index of the i-th ungapped
              reference nucleotide.
    """
    col_map = []
    for col_idx, char in enumerate(ref_seq_str):
        if char != "-":
            col_map.append(col_idx)
    return col_map


def stitch_segments(seq_str, segments, col_map):
    """Slice and stitch MSA segments for a gene (exon joining).

    For each CDS segment, maps reference coords → MSA columns, slices
    the sequence, and concatenates. Works identically for single-segment
    genes (trivial case).

    Args:
        seq_str:  Full aligned sequence string for one record.
        segments: List of (start0, end, strand) tuples.
        col_map:  Reference-position → MSA-column mapping.

    Returns:
        tuple: (stitched_seq_str, strand)
    """
    strand = segments[0][2]  # All segments share the same strand
    pieces = []

    for seg_start, seg_end, _ in segments:
        seg_start = max(0, seg_start)
        seg_end = min(seg_end, len(col_map))
        if seg_start >= seg_end:
            continue

        msa_start = col_map[seg_start]
        msa_end = col_map[seg_end - 1] + 1
        pieces.append(seq_str[msa_start:msa_end].upper())

    return "".join(pieces), strand


def apply_harsh_masking(seq_str, strand):
    """Apply Harsh Masking for reading-frame integrity and HyPhy compatibility.

    Assumptions: seq_str length is a multiple of 3 (caller ensures this).

    Operations:
      - Reverse-complement if strand == '-'.
      - Partial gap codons (frameshifts) → 'NNN'.
      - In-frame stop codons            → '---'.

    Returns:
        tuple: (cleaned_str, n_frameshift_masked, n_stop_masked)
    """
    if strand == "-":
        seq_str = str(Seq(seq_str).reverse_complement())

    stop_codons = {"TAA", "TAG", "TGA"}
    cleaned_seq = []
    n_frameshift = 0
    n_stop = 0

    for i in range(0, len(seq_str), 3):
        codon = seq_str[i:i + 3]

        # Partial gap → frameshift artifact
        if "-" in codon and codon != "---":
            cleaned_seq.append("NNN")
            n_frameshift += 1
            continue

        # In-frame stop codon
        if codon in stop_codons:
            cleaned_seq.append("---")
            n_stop += 1
            continue

        cleaned_seq.append(codon)

    return "".join(cleaned_seq), n_frameshift, n_stop


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--alignment", required=True,
                        help="MSA from augur align (FASTA)")
    parser.add_argument("--gff", required=True,
                        help="Reference GFF3 annotation")
    parser.add_argument("--reference-id", required=True,
                        help="ID of the reference sequence in the MSA")
    parser.add_argument("--outdir", required=True,
                        help="Directory to save the cleaned gene alignments")
    args = parser.parse_args()

    print(f"Parsing GFF: {args.gff}")
    gene_segments = parse_gff(args.gff)
    if not gene_segments:
        sys.exit("ERROR: No genes found in GFF")

    for gene, segs in gene_segments.items():
        n_seg = len(segs)
        seg_info = ", ".join(f"{s}-{e}" for s, e, _ in segs)
        label = "spliced" if n_seg > 1 else "single-CDS"
        print(f"  [{label}] {gene}: {n_seg} segment(s) → {seg_info}")

    print(f"\nLoading MSA: {args.alignment}")
    records = list(SeqIO.parse(args.alignment, "fasta"))
    print(f"Loaded {len(records)} sequences.")

    ref_record = next(
        (r for r in records
         if r.id == args.reference_id
         or r.id.split(".")[0] == args.reference_id.split(".")[0]),
        None
    )
    if not ref_record:
        sys.exit(f"ERROR: Reference '{args.reference_id}' not found in MSA.")

    print(f"Mapping coordinates against reference: {ref_record.id}")
    col_map = build_coordinate_map(str(ref_record.seq))

    os.makedirs(args.outdir, exist_ok=True)

    for gene_name, segments in gene_segments.items():
        is_spliced = len(segments) > 1

        # Calculate total expected stitched length for QC
        total_ref_len = sum(
            min(e, len(col_map)) - max(0, s)
            for s, e, _ in segments
            if max(0, s) < min(e, len(col_map))
        )
        pad_len = (3 - (total_ref_len % 3)) % 3

        out_records = []
        total_frameshift = 0
        total_stop = 0
        skipped_gappy = 0

        for r in records:
            seq_str = str(r.seq)
            stitched, strand = stitch_segments(seq_str, segments, col_map)

            # Pad to multiple of 3
            stitched += "N" * pad_len

            # Filter: skip sequences with >90% gaps across the gene region
            gap_frac = stitched.count("-") / len(stitched) if stitched else 1.0
            if gap_frac > 0.9:
                skipped_gappy += 1
                continue

            cleaned_seq_str, n_fs, n_st = apply_harsh_masking(stitched, strand)
            total_frameshift += n_fs
            total_stop += n_st
            out_records.append(
                SeqRecord(Seq(cleaned_seq_str), id=r.id, description="")
            )

        out_path = os.path.join(args.outdir, f"cleaned_{gene_name}_NT.fasta")
        SeqIO.write(out_records, out_path, "fasta")

        n_codons = (total_ref_len + pad_len) // 3
        splice_tag = f" [SPLICED, {len(segments)} exons]" if is_spliced else ""
        print(
            f"[{gene_name}]{splice_tag} {len(out_records)} seqs written, "
            f"{skipped_gappy} skipped (>90% gap), "
            f"{n_codons} codons/seq, "
            f"masked: {total_frameshift} frameshifts + {total_stop} stop codons "
            f"→ {out_path}"
        )


if __name__ == "__main__":
    main()
