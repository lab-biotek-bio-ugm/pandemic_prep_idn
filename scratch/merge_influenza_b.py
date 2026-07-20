import os
import glob
import pandas as pd

def main():
    input_dir = "/mnt/d/Dataset/B_MP"
    output_dir = "/home/nanda/projects/pandemic-preparedness/dataset/EpiFlu/B_MP"
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Merge FASTA files
    fasta_files = sorted(glob.glob(os.path.join(input_dir, "*_Sequence.fasta")))
    combined_fasta_path = os.path.join(output_dir, "EpiFlu_B_MP_Combined_Sequence.fasta")
    
    print(f"Merging FASTA files:")
    total_seq_count = 0
    with open(combined_fasta_path, "w") as outfile:
        for fpath in fasta_files:
            country = os.path.basename(fpath).split("_")[0]
            print(f"  - {fpath} ({country})")
            with open(fpath, "r") as infile:
                for line in infile:
                    if line.startswith(">"):
                        total_seq_count += 1
                    outfile.write(line)
    print(f"Combined FASTA written to {combined_fasta_path} (Total sequences: {total_seq_count})")
    
    # 2. Merge Metadata files (which are actually .xls binary Excel files)
    metadata_files = sorted(glob.glob(os.path.join(input_dir, "*_Metadata.csv")))
    combined_metadata_path = os.path.join(output_dir, "EpiFlu_B_MP_Combined_Metadata.csv")
    
    print(f"\nMerging Metadata files:")
    df_list = []
    for fpath in metadata_files:
        country = os.path.basename(fpath).split("_")[0]
        print(f"  - {fpath} ({country})")
        # Read Excel file
        df = pd.read_excel(fpath)
        print(f"    Loaded {len(df)} rows. Columns: {list(df.columns[:3])}...")
        df_list.append(df)
        
    combined_df = pd.concat(df_list, ignore_index=True)
    # Save as true CSV
    combined_df.to_csv(combined_metadata_path, index=False)
    print(f"Combined Metadata CSV written to {combined_metadata_path} (Total rows: {len(combined_df)})")

if __name__ == "__main__":
    main()
