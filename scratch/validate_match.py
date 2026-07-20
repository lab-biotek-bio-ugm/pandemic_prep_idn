import pandas as pd
from Bio import SeqIO

def main():
    fasta_path = "/home/nanda/projects/pandemic-preparedness/dataset/EpiFlu/B_MP/EpiFlu_B_MP_Combined_Sequence.fasta"
    metadata_path = "/home/nanda/projects/pandemic-preparedness/dataset/EpiFlu/B_MP/EpiFlu_B_MP_Combined_Metadata.csv"
    
    print("Reading FASTA sequences...")
    fasta_ids = []
    fasta_duplicates = 0
    seen_fasta = set()
    for record in SeqIO.parse(fasta_path, "fasta"):
        # Extract the lookup ID (first part of header before |)
        fid = record.id.split("|")[0].strip()
        if fid in seen_fasta:
            fasta_duplicates += 1
        else:
            seen_fasta.add(fid)
            fasta_ids.append(fid)
            
    print(f"Total FASTA records: {len(seen_fasta) + fasta_duplicates}")
    print(f"Unique FASTA IDs: {len(seen_fasta)}")
    print(f"Duplicate FASTA IDs: {fasta_duplicates}")
    
    print("\nReading Metadata CSV...")
    df = pd.read_csv(metadata_path, dtype=str)
    
    # Check if 'Isolate_Id' exists
    if 'Isolate_Id' not in df.columns:
        print("ERROR: Isolate_Id column not found in metadata!")
        print(f"Columns are: {list(df.columns)}")
        return
        
    metadata_ids = df['Isolate_Id'].dropna().str.strip().tolist()
    unique_meta_ids = set(metadata_ids)
    
    print(f"Total rows in metadata: {len(df)}")
    print(f"Unique Isolate_Ids in metadata: {len(unique_meta_ids)}")
    print(f"Duplicate Isolate_Ids in metadata: {len(metadata_ids) - len(unique_meta_ids)}")
    
    # Check intersection
    matched = seen_fasta.intersection(unique_meta_ids)
    only_fasta = seen_fasta - unique_meta_ids
    only_meta = unique_meta_ids - seen_fasta
    
    print(f"\nIntersection Analysis:")
    print(f"  - Matched IDs (present in both): {len(matched)}")
    print(f"  - IDs only in FASTA (no metadata): {len(only_fasta)}")
    print(f"  - IDs only in Metadata (no FASTA): {len(only_meta)}")
    
    if len(only_fasta) > 0:
        print(f"  - Examples of IDs only in FASTA: {list(only_fasta)[:5]}")
    if len(only_meta) > 0:
        print(f"  - Examples of IDs only in Metadata: {list(only_meta)[:5]}")

if __name__ == "__main__":
    main()
