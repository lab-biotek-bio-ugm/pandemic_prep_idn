import argparse
import sys
import subprocess
from pathlib import Path
from logging import getLogger, basicConfig, INFO

basicConfig(level=INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(description="Early Filter using Augur")
    parser.add_argument("--sequences", required=True, help="Input FASTA")
    parser.add_argument("--metadata", required=True, help="Input metadata TSV")
    parser.add_argument("--output-sequences", required=True, help="Output FASTA")
    parser.add_argument("--output-metadata", required=True, help="Output metadata TSV")
    parser.add_argument("--query", help="Pandas query string")
    parser.add_argument("--min-date", help="Minimum date")
    parser.add_argument("--max-date", help="Maximum date")
    
    args = parser.parse_args()
    
    cmd = [
        "augur", "filter",
        "--sequences", args.sequences,
        "--metadata", args.metadata,
        "--output-sequences", args.output_sequences,
        "--output-metadata", args.output_metadata
    ]
    
    # Optional flags
    if args.query and args.query != "None" and args.query.strip():
        cmd.extend(["--query", args.query])
        
    if args.min_date and args.min_date != "None" and args.min_date.strip():
        cmd.extend(["--min-date", args.min_date])
        
    if args.max_date and args.max_date != "None" and args.max_date.strip():
        cmd.extend(["--max-date", args.max_date])
        
    logger.info(f"Executing: {' '.join(cmd)}")
    
    try:
        subprocess.run(cmd, check=True, stderr=subprocess.PIPE, text=True)
        logger.info("Early filter completed successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Augur filter failed:\n{e.stderr}")
        sys.exit(1)

if __name__ == "__main__":
    main()
