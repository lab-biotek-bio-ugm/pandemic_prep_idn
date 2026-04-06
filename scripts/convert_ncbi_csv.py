#!/usr/bin/env python3
"""
Convert NCBI Virus web-portal CSV metadata to the TSV format
expected by the Augur pipeline.

Input : CSV from https://www.ncbi.nlm.nih.gov/labs/virus/
Output: TSV with standardised column names matching Augur expectations.
"""
import argparse
import csv
import sys


# Mapping: NCBI Virus CSV column → pipeline TSV column
COLUMN_MAP = {
    "Accession": "accession",
    "Organism_Name": "virus_name",
    "Species": "species",
    "Pangolin": "pangolin_lineage",
    "Isolate": "isolate",
    "Length": "length",
    "Nuc_Completeness": "completeness",
    "Geo_Location": "geo_location",
    "Host": "host",
    "Collection_Date": "date",
    "Release_Date": "release_date",
    "Molecule_type": "molecule_type",
}


def parse_geo(geo_location: str):
    """Split 'Thailand:Chiang Mai' → country='Thailand', region='Chiang Mai'."""
    parts = geo_location.split(":", 1)
    country = parts[0].strip() if parts else ""
    region = parts[1].strip() if len(parts) > 1 else ""
    return country, region


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="NCBI Virus CSV file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    with open(args.input, "r", encoding="utf-8") as fin, \
         open(args.output, "w", newline="", encoding="utf-8") as fout:

        reader = csv.DictReader(fin)
        out_fields = list(COLUMN_MAP.values()) + ["country", "region"]
        writer = csv.DictWriter(fout, fieldnames=out_fields, delimiter="\t")
        writer.writeheader()

        count = 0
        for row in reader:
            out_row = {}
            for csv_col, tsv_col in COLUMN_MAP.items():
                out_row[tsv_col] = row.get(csv_col, "")
            country, region = parse_geo(row.get("Geo_Location", ""))
            out_row["country"] = country
            out_row["region"] = region
            writer.writerow(out_row)
            count += 1

        print(f"Converted {count} records → {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
