import urllib.request
import zipfile
import io
import sys
from pathlib import Path

def download_h3n2_ref():
    url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/virus/accession/NC_007366.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF"
    out_dir = Path("config/references")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    fasta_out = out_dir / "flu_a_h3n2_ha.fasta"
    gff_out = out_dir / "flu_a_h3n2_ha.gff"
    
    print(f"Downloading H3N2 reference data from NCBI Datasets API...", file=sys.stderr)
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    
    try:
        with urllib.request.urlopen(req) as response:
            zip_data = response.read()
            
        print("Extracting ZIP archive...", file=sys.stderr)
        with zipfile.ZipFile(io.BytesIO(zip_data)) as z:
            for file_info in z.infolist():
                if file_info.filename.endswith(".fna"):
                    print(f"Extracting FASTA...", file=sys.stderr)
                    with z.open(file_info) as f:
                        fasta_out.write_bytes(f.read())
                elif file_info.filename.endswith(".gff"):
                    print(f"Extracting GFF...", file=sys.stderr)
                    with z.open(file_info) as f:
                        gff_out.write_bytes(f.read())
                        
        print(f"Success! Saved to:\n- {fasta_out}\n- {gff_out}", file=sys.stderr)
    except Exception as e:
        print(f"Failed to download or extract: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    download_h3n2_ref()
