import os
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Paths
REF_PATH = "config/references/sars_cov_2_reference.fasta"
OUT_FASTA = "results/sars-cov-2-trial/m01/raw_sequences.fasta"
OUT_TSV = "results/sars-cov-2-trial/m01/metadata.tsv"
OUT_JSONL = "results/sars-cov-2-trial/m01/raw_metadata.jsonl" # dummy to satisfy snakemake

os.makedirs("results/sars-cov-2-trial/m01", exist_ok=True)

# Read reference
ref_record = next(SeqIO.parse(REF_PATH, "fasta"))
ref_seq = str(ref_record.seq)

records = []
metadata_lines = ["strain\tdate\tregion\tcountry\tdivision\tlength\tdate_submitted"]

# Generate 20 dummy Indonesian strains with random mutations
for i in range(1, 21):
    strain_name = f"Indonesia/Dummy_Strain_{i:02d}/2023"
    
    # Introduce 1 to 5 random mutations
    mutated_seq = list(ref_seq)
    for _ in range(random.randint(1, 5)):
        pos = random.randint(100, len(mutated_seq) - 100)
        mutated_seq[pos] = random.choice(['A', 'C', 'G', 'T'])
        
    seq_str = "".join(mutated_seq)
    rec = SeqRecord(Seq(seq_str), id=strain_name, description="")
    records.append(rec)
    
    # Metadata
    month = random.randint(1, 12)
    day = random.randint(1, 28)
    date = f"2023-{month:02d}-{day:02d}"
    
    metadata_lines.append(f"{strain_name}\t{date}\tAsia\tIndonesia\tJakarta\t{len(seq_str)}\t{date}")

# Write to disk
SeqIO.write(records, OUT_FASTA, "fasta")

with open(OUT_TSV, "w") as f:
    f.write("\n".join(metadata_lines))
    
with open(OUT_JSONL, "w") as f:
    f.write("{}")

print("✅ Dummy SARS-CoV-2 dataset created successfully in results/sars-cov-2-trial/m01/")
print(f"Generated {len(records)} sequences and metadata entries.")
