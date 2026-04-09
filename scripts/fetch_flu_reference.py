import sys
import os
import argparse
from Bio import Entrez, SeqIO

Entrez.email = "example@example.com"

# Dictionary of standard Influenza A (PR8) Reference Segments
REF_SEGMENTS = {
    "PB2": {"id": "NC_002023.1", "name": "flu_a_pb2"},
    "PB1": {"id": "NC_002021.1", "name": "flu_a_pb1"},
    "PA":  {"id": "NC_002022.1", "name": "flu_a_pa"},
    "HA":  {"id": "NC_002017.1", "name": "flu_a_ha"},
    "NP":  {"id": "NC_002019.1", "name": "flu_a_np"},
    "NA":  {"id": "NC_002018.1", "name": "flu_a_na"},
    "MP":  {"id": "NC_002016.1", "name": "flu_a_mp"},
    "NS":  {"id": "NC_002020.1", "name": "flu_a_ns"}
}

def fetch_refseq(segment_name, out_dir):
    if segment_name not in REF_SEGMENTS:
        print(f"Error: Unknown segment {segment_name}. Choose from {list(REF_SEGMENTS.keys())}")
        sys.exit(1)
        
    ref_info = REF_SEGMENTS[segment_name]
    ref_id = ref_info["id"]
    base_name = ref_info["name"]
    
    os.makedirs(out_dir, exist_ok=True)
    fasta_path = os.path.join(out_dir, f"{base_name}_reference.fasta")
    gff_path = os.path.join(out_dir, f"{base_name}_annotation.gff")
    
    print(f"Fetching {segment_name} ({ref_id}) from NCBI...")
    
    try:
        # Fetch GenBank record
        with Entrez.efetch(db="nucleotide", id=ref_id, rettype="gb", retmode="text") as handle:
            record = SeqIO.read(handle, "genbank")
            
        # 1. Write FASTA
        with open(fasta_path, "w") as f:
            f.write(f">{ref_id} {record.description}\n{record.seq}\n")
        print(f"Saved FASTA: {fasta_path}")
        
        # 2. Write Simplified GFF
        # Note: We skip spliced genes (like M2 or NS2) that have CompoundLocations 
        # to ensure augur translate mathematically passes the multi-of-3 gap test.
        with open(gff_path, "w") as f:
            f.write("##gff-version 3\n")
            for feature in record.features:
                if feature.type == "CDS":
                    gene_name = feature.qualifiers.get("gene", ["unknown"])[0]
                    # If the location type has "join" (Biopython returns a CompoundLocation), skip it for stability
                    if "CompoundLocation" in str(type(feature.location)):
                        print(f"  Skipping spliced CDS '{gene_name}' to avoid augur translation gap errors.")
                        continue
                        
                    start_1b = feature.location.start + 1
                    end_1b = feature.location.end
                    
                    if (end_1b - start_1b + 1) % 3 != 0:
                        print(f"  Warning: CDS '{gene_name}' length is not divisible by 3. Skipping.")
                        continue
                        
                    gff_line = f"{ref_id}\tRefSeq\tCDS\t{start_1b}\t{end_1b}\t.\t+\t0\tgene={gene_name};\n"
                    f.write(gff_line)
                    print(f"  Added GFF Feature: {gene_name} ({start_1b}-{end_1b})")
                    
        print(f"Saved GFF: {gff_path}")
        
    except Exception as e:
        print(f"Failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch Influenza A Reference")
    parser.add_argument("--segment", required=True, choices=REF_SEGMENTS.keys(), help="Segment to fetch (e.g. MP, HA, NP)")
    parser.add_argument("--outdir", default="config/references", help="Output directory")
    args = parser.parse_args()
    
    fetch_refseq(args.segment, args.outdir)
