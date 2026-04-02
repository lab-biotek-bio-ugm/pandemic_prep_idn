"""
Module 1: Data Acquisition
Handles fetching raw viral sequences from NCBI GenBank and formatting
the metadata into the standard TSV structure required by Nextstrain Augur.
"""

rule download_ncbi_data:
    output:
        fasta=f"results/{BUILD}/m01/raw_sequences.fasta",
        jsonl=f"results/{BUILD}/m01/raw_metadata.jsonl"
    params:
        organism=config["ncbi_query"]["organism"],
        location=config["ncbi_query"]["location"],
        zip_file=f"results/{BUILD}/m01/ncbi_dataset.zip"
    conda:
        "../envs/nextstrain.yaml"
    log:
        f"logs/{BUILD}/m01_download.log"
    shell:
        """
        mkdir -p results/{BUILD}/m01 logs/{BUILD}

        datasets download virus genome taxon "{params.organism}" \
            --geo-location "{params.location}" \
            --filename {params.zip_file} 2> {log}

        unzip -o -p {params.zip_file} ncbi_dataset/data/genomic.fna > {output.fasta}
        unzip -o -p {params.zip_file} ncbi_dataset/data/data_report.jsonl > {output.jsonl}

        rm -f {params.zip_file}
        """

rule convert_metadata_tsv:
    input:
        jsonl=f"results/{BUILD}/m01/raw_metadata.jsonl"
    output:
        tsv=f"results/{BUILD}/m01/metadata.tsv"
    run:
        import json
        
        with open(input.jsonl, 'r') as f_in, open(output.tsv, 'w') as f_out:
            # Enforce Nextstrain standard TSV headers
            f_out.write("strain\tdate\tregion\tcountry\tdivision\thost\tlength\n")
            
            for line in f_in:
                if not line.strip(): 
                    continue
                    
                data = json.loads(line)
                strain = data.get("accession", "Unknown")
                
                isolate = data.get("isolate", {})
                date = isolate.get("collection_date", "Unknown")
                
                # Impute missing temporal data for Augur compatibility
                if len(date) == 4: 
                    date += "-XX-XX"
                elif len(date) == 7: 
                    date += "-XX"
                
                # Parse NCBI nested geographic location format (e.g., Continent: Country: Division)
                location = data.get("location", {})
                geo_string = location.get("geographic_location", "")
                parts = [p.strip() for p in geo_string.split(":")]
                
                region = parts[0] if len(parts) > 0 and parts[0] else location.get("geographic_region", "Unknown")
                country = parts[1] if len(parts) > 1 else "Unknown"
                division = parts[2] if len(parts) > 2 else "Unknown"
                
                host = data.get("host", {}).get("common_name", "Unknown")
                length = str(data.get("sequence_length", "0"))
                
                f_out.write(f"{strain}\t{date}\t{region}\t{country}\t{division}\t{host}\t{length}\n")
