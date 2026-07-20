import os
from Bio import Entrez, SeqIO

def main():
    # Set email for NCBI Entrez
    Entrez.email = "example@example.com"
    
    accession = "NC_002210.1"
    print(f"Fetching {accession} from NCBI...")
    
    # Fetch GenBank record
    try:
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        
        print(f"Successfully fetched GenBank record: {record.id}")
        print(f"Length: {len(record)} bp")
        print(f"Description: {record.description}")
        
        # Save FASTA
        ref_dir = "/home/nanda/projects/pandemic-preparedness/config/references"
        os.makedirs(ref_dir, exist_ok=True)
        fasta_path = os.path.join(ref_dir, "flu_b_mp_reference.fasta")
        SeqIO.write(record, fasta_path, "fasta")
        print(f"Saved reference FASTA to {fasta_path}")
        
        # Parse features and write GFF3
        gff_path = os.path.join(ref_dir, "flu_b_mp_annotation.gff")
        with open(gff_path, "w") as f:
            f.write("##gff-version 3\n")
            f.write(f"##sequence-region {accession} 1 {len(record)}\n")
            # Write segment feature
            f.write(f"{accession}\tNCBI\tgene\t1\t{len(record)}\t.\t+\t.\tID=segment-MP\n")
            
            for feature in record.features:
                if feature.type == "CDS":
                    gene_name = feature.qualifiers.get("gene", ["Unknown"])[0]
                    # GenBank is 0-indexed half-open, GFF is 1-indexed closed
                    start = int(feature.location.start) + 1
                    end = int(feature.location.end)
                    strand = "+" if feature.location.strand >= 0 else "-"
                    
                    # Write gene feature
                    f.write(f"{accession}\tNCBI\tgene\t{start}\t{end}\t.\t{strand}\t.\tID=gene-{gene_name};gene={gene_name}\n")
                    # Write CDS feature
                    f.write(f"{accession}\tNCBI\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds-{gene_name};Parent=gene-{gene_name};gene={gene_name}\n")
                    
                    print(f"Feature: {gene_name} CDS -> {start} to {end}")
                    
        print(f"Saved GFF3 annotation to {gff_path}")
        
    except Exception as e:
        print(f"ERROR: {e}")

if __name__ == "__main__":
    main()
