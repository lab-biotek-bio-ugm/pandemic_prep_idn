"""
Enrich base metadata TSV with Nextclade classification and QC results.

Merges clade (Nextstrain), lineage (Pango), and quality metrics
from Nextclade output into the pipeline's main metadata file.
"""

import pandas as pd
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Enrich base metadata with Nextclade TSV output"
    )
    parser.add_argument("--base-metadata", required=True, help="Original metadata TSV")
    parser.add_argument("--nextclade-tsv", required=True, help="Nextclade output TSV")
    parser.add_argument("--output", required=True, help="Output enriched metadata TSV")
    
    args = parser.parse_args()
    
    try:
        # Load base metadata
        df_base = pd.read_csv(args.base_metadata, sep='\t', dtype=str)
        
        # Load Nextclade TSV
        df_nc = pd.read_csv(args.nextclade_tsv, sep='\t', dtype=str)
        
        # Dynamic seqName matching to handle NCBI vs GISAID FASTA headers
        if 'seqName' in df_nc.columns:
            base_strains = set(df_base['strain'])
            matched_seqnames = []
            
            for seq_name in df_nc['seqName']:
                if seq_name in base_strains:
                    matched_seqnames.append(seq_name) # Exact match
                elif seq_name.split()[0] in base_strains:
                    matched_seqnames.append(seq_name.split()[0]) # NCBI style (space-separated)
                elif seq_name.split('|')[0] in base_strains:
                    matched_seqnames.append(seq_name.split('|')[0]) # GISAID style (pipe-separated)
                else:
                    matched_seqnames.append(seq_name) # Fallback
                    
            df_nc['strain'] = matched_seqnames
            
        # Select important columns from Nextclade to avoid TSV bloat
        nc_cols_wanted = [
            'strain',
            'clade',
            'Nextclade_pango',
            'qc.overallStatus',
            'qc.overallScore',
            'totalSubstitutions',
            'totalDeletions',
            'totalInsertions',
            'totalFrameShifts',
            'totalMissing',
        ]
        # Only keep columns that actually exist in the Nextclade output
        nc_cols = [c for c in nc_cols_wanted if c in df_nc.columns]
        
        if not nc_cols or 'strain' not in nc_cols:
            print("WARNING: Nextclade output does not contain expected columns. "
                  "Copying base metadata as-is.", file=sys.stderr)
            df_base.to_csv(args.output, sep='\t', index=False)
            return
        
        df_nc_subset = df_nc[nc_cols].copy()
        
        # Rename columns to Augur/Auspice-friendly names
        rename_map = {
            'Nextclade_pango': 'lineage',
            'qc.overallStatus': 'qc_status',
            'qc.overallScore': 'qc_score',
            'totalSubstitutions': 'total_substitutions',
            'totalDeletions': 'total_deletions',
            'totalInsertions': 'total_insertions',
            'totalFrameShifts': 'total_frameshifts',
            'totalMissing': 'total_missing',
        }
        df_nc_subset = df_nc_subset.rename(columns=rename_map)
        
        # Merge (left join so we keep all base rows even if Nextclade missed some)
        # Drop overlapping columns from base to prevent _x/_y duplication
        overlap_cols = set(df_nc_subset.columns) & set(df_base.columns) - {'strain'}
        if overlap_cols:
            df_base = df_base.drop(columns=list(overlap_cols))
        df_merged = pd.merge(df_base, df_nc_subset, on='strain', how='left')
        
        # Write out
        df_merged.to_csv(args.output, sep='\t', index=False)
        
        n_total = len(df_merged)
        n_claded = df_merged['clade'].notna().sum() if 'clade' in df_merged.columns else 0
        print(f"Successfully enriched metadata: {n_claded}/{n_total} sequences classified.")
        
    except Exception as e:
        print(f"Error enriching metadata: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
