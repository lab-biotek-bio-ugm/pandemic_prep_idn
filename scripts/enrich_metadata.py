import pandas as pd
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Enrich base metadata with Nextclade TSV output")
    parser.add_argument("--base-metadata", required=True, help="Original metadata TSV")
    parser.add_argument("--nextclade-tsv", required=True, help="Nextclade output TSV")
    parser.add_argument("--output", required=True, help="Output enriched metadata TSV")
    
    args = parser.parse_args()
    
    try:
        # Load base metadata
        # Nextstrain generally expects 'strain' as the index column.
        df_base = pd.read_csv(args.base_metadata, sep='\t')
        
        # Load Nextclade TSV
        # Nextclade outputs 'seqName' which corresponds to 'strain'
        df_nc = pd.read_csv(args.nextclade_tsv, sep='\t')
        
        # Rename Nextclade's key to match base metadata
        if 'seqName' in df_nc.columns:
            df_nc = df_nc.rename(columns={'seqName': 'strain'})
            
        # Select important columns from Nextclade to avoid TSV bloat
        # Assuming we care about clades, lineages, and QC metrics.
        nc_cols_to_keep = ['strain', 'clade', 'Nextclade_pango', 'qc.overallStatus', 'totalSubstitutions']
        # only keep columns that actually exist in the nextclade output
        nc_cols = [c for c in nc_cols_to_keep if c in df_nc.columns]
        
        df_nc_subset = df_nc[nc_cols]
        
        # Rename so they look pretty in Auspice
        df_nc_subset = df_nc_subset.rename(columns={
            'Nextclade_pango': 'lineage',
            'qc.overallStatus': 'qc_status',
            'totalSubstitutions': 'total_mutations'
        })
        
        # Merge by strain (Left merge ensures we don't drop rows if nextclade failed on some)
        df_merged = pd.merge(df_base, df_nc_subset, on='strain', how='left')
        
        # Write out
        df_merged.to_csv(args.output, sep='\t', index=False)
        print(f"Successfully merged Nextclade data into {args.output}")
        
    except Exception as e:
        print(f"Error enriching metadata: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
