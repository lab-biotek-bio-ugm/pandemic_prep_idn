from Bio import SeqIO

def main():
    fasta_path = "/home/nanda/projects/pandemic-preparedness/dataset/EpiFlu/B_MP/EpiFlu_B_MP_Combined_Sequence.fasta"
    cleaned_path = "/home/nanda/projects/pandemic-preparedness/dataset/EpiFlu/B_MP/EpiFlu_B_MP_Combined_Sequence_Cleaned.fasta"
    
    seen_ids = set()
    records_to_keep = []
    duplicate_count = 0
    
    for record in SeqIO.parse(fasta_path, "fasta"):
        fid = record.id.split("|")[0].strip()
        if fid in seen_ids:
            duplicate_count += 1
        else:
            seen_ids.add(fid)
            records_to_keep.append(record)
            
    SeqIO.write(records_to_keep, cleaned_path, "fasta")
    print(f"Removed {duplicate_count} duplicate sequences.")
    print(f"Saved {len(records_to_keep)} unique sequences to {cleaned_path}")
    
    # Overwrite the original
    import os
    os.replace(cleaned_path, fasta_path)
    print("Overwrote original combined FASTA with the cleaned version.")

if __name__ == "__main__":
    main()
