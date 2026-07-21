#!/usr/bin/env python3
"""
NCBI Measles Dataset Standardizer & Reference Fetcher
Standardizes NCBI Virus CSV metadata and FASTA headers for 100% Augur/Nextstrain compatibility.
Also fetches RefSeq accession NC_001498.1 for Measles reference FASTA & GFF annotation.
"""

import os
import sys
import re
import pandas as pd
from Bio import Entrez, SeqIO

# Set Entrez email
Entrez.email = "bioinformatics@ugm.ac.id"

# Country -> Region mapping
COUNTRY_TO_REGION = {
    "India": "Asia", "China": "Asia", "Japan": "Asia", "Malaysia": "Asia", "Indonesia": "Asia",
    "Thailand": "Asia", "Vietnam": "Asia", "Philippines": "Asia", "Singapore": "Asia", "South Korea": "Asia",
    "Australia": "Oceania", "New Zealand": "Oceania",
    "USA": "North America", "Canada": "North America", "Mexico": "North America",
    "UK": "Europe", "Germany": "Europe", "France": "Europe", "Italy": "Europe", "Spain": "Europe",
    "South Africa": "Africa", "Nigeria": "Africa", "Kenya": "Africa", "Egypt": "Africa"
}

def clean_dataset():
    csv_path = "dataset/NCBI/Measles/sequences.csv"
    fasta_path = "dataset/NCBI/Measles/sequences.fasta"
    
    out_fasta = "dataset/NCBI/Measles/sequences_clean.fasta"
    out_tsv = "dataset/NCBI/Measles/metadata_clean.tsv"
    
    print(f"[Measles Prep] Reading metadata CSV: {csv_path}")
    df = pd.read_csv(csv_path)
    
    # Fill missing country from Geo_Location
    def parse_country(row):
        c = row.get("Country")
        if pd.notna(c) and str(c).strip() != "":
            return str(c).strip()
        geo = row.get("Geo_Location")
        if pd.notna(geo):
            return str(geo).split(":")[0].strip()
        return "Unknown"
        
    df["country"] = df.apply(parse_country, axis=1)
    df["region"] = df["country"].map(lambda c: COUNTRY_TO_REGION.get(c, "Unknown"))
    
    # Map host
    df["host"] = df["Host"].fillna("Human").map(lambda h: "Human" if "Homo" in str(h) or "human" in str(h).lower() else str(h))
    
    # Map standard columns
    meta_df = pd.DataFrame()
    meta_df["strain"] = df["Accession"].astype(str)
    meta_df["virus"] = "measles"
    meta_df["gisaid_epi_isl"] = df["Accession"].astype(str)
    meta_df["date"] = df["Collection_Date"].fillna("XXXX-XX-XX")
    meta_df["date_submitted"] = df["Release_Date"].fillna("XXXX-XX-XX")
    meta_df["country"] = df["country"]
    meta_df["region"] = df["region"]
    meta_df["host"] = df["host"]
    meta_df["genotype"] = df["Genotype"].fillna("Unknown")
    meta_df["isolate"] = df["Isolate"].fillna("")
    meta_df["length"] = df["Length"]
    
    meta_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[Measles Prep] Wrote sanitized metadata TSV: {out_tsv} ({len(meta_df)} rows)")
    
    # Clean FASTA headers to match Accession ID exactly
    print(f"[Measles Prep] Processing FASTA: {fasta_path}")
    written_records = 0
    with open(fasta_path, "r") as fin, open(out_fasta, "w") as fout:
        for record in SeqIO.parse(fin, "fasta"):
            # record.id e.g. PZ660920.1
            raw_id = record.id.split(".")[0]  # PZ660920
            fout.write(f">{raw_id}\n{str(record.seq)}\n")
            written_records += 1
            
    print(f"[Measles Prep] Wrote cleaned FASTA: {out_fasta} ({written_records} records)")

def fetch_measles_reference():
    ref_accession = "NC_001498.1"
    ref_fasta_out = "config/references/measles_reference.fasta"
    ref_gff_out = "config/references/measles_annotation.gff"
    
    os.makedirs("config/references", exist_ok=True)
    
    print(f"[Measles Ref] Fetching NCBI RefSeq {ref_accession}...")
    handle = Entrez.efetch(db="nuccore", id=ref_accession, rettype="gbwithparts", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    
    with open(ref_fasta_out, "w") as f:
        f.write(f">{record.id}\n{str(record.seq)}\n")
    print(f"[Measles Ref] Saved RefSeq FASTA: {ref_fasta_out}")
    
    # Generate GFF3 annotation
    gff_lines = ["##gff-version 3"]
    gff_lines.append(f"{record.id}\tNCBI\tregion\t1\t{len(record.seq)}\t.\t+\t.\tID=region-1;Name=Measles_genome")
    
    for feat in record.features:
        if feat.type in ["gene", "CDS"]:
            start = int(feat.location.start) + 1
            end = int(feat.location.end)
            strand = "+" if feat.location.strand >= 0 else "-"
            gene_name = feat.qualifiers.get("gene", ["unknown"])[0]
            product = feat.qualifiers.get("product", [""])[0]
            
            if feat.type == "gene":
                gff_lines.append(f"{record.id}\tNCBI\tgene\t{start}\t{end}\t.\t{strand}\t.\tID=gene-{gene_name};gene={gene_name}")
            elif feat.type == "CDS":
                gff_lines.append(f"{record.id}\tNCBI\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds-{gene_name};Parent=gene-{gene_name};gene={gene_name};product={product}")
                
    with open(ref_gff_out, "w") as f:
        f.write("\n".join(gff_lines) + "\n")
    print(f"[Measles Ref] Saved GFF3 annotation: {ref_gff_out}")

def main():
    clean_dataset()
    fetch_measles_reference()

if __name__ == "__main__":
    main()
