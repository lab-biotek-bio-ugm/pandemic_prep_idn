import urllib.request
import os
from Bio import SeqIO

# Create directory
os.makedirs("config/references", exist_ok=True)

# 1. Download FASTA
fasta_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_002016.1&rettype=fasta&retmode=text"
print("Downloading FASTA for NC_002016...")
urllib.request.urlretrieve(fasta_url, "config/references/flu_a_mp_reference.fasta")

# 2. Download GenBank to extract features and convert to GFF3
gb_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_002016.1&rettype=gb&retmode=text"
print("Downloading GenBank for NC_002016...")
urllib.request.urlretrieve(gb_url, "config/references/temp.gb")

print("Processing to GFF3...")
with open("config/references/temp.gb") as in_handle, \
     open("config/references/flu_a_mp_annotation.gff", "w", encoding="utf-8") as out_handle:
    
    out_handle.write("##gff-version 3\n")
    
    for record in SeqIO.parse(in_handle, "genbank"):
        seq_id = record.id
        length = len(record.seq)
        out_handle.write(f"##sequence-region {seq_id} 1 {length}\n")
        
        for feature in record.features:
            if feature.type == "source":
                continue
                
            start = feature.location.start + 1
            end = feature.location.end
            strand = "+" if feature.location.strand == 1 else "-"
            phase = "0" if feature.type == "CDS" else "."
            
            # Simple attribute string
            qualifiers = []
            for k, v in feature.qualifiers.items():
                val = "%2C".join(v) if isinstance(v, list) else str(v)
                # Keep it simple, Biopython parses some as lists
                qualifiers.append(f"{k}={val}")
                
            attr_str = ";".join(qualifiers)
            if not attr_str:
                attr_str = "ID=unknown"
                
            out_handle.write(f"{seq_id}\tRefSeq\t{feature.type}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attr_str}\n")

os.remove("config/references/temp.gb")
print("Done! Downloaded FASTA and generated GFF3.")
