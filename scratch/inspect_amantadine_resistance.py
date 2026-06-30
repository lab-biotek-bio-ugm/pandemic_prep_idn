import pandas as pd

def main():
    print("Reading metadata_enriched.tsv...")
    df = pd.read_csv('ssh-output/results/influenza-a-mp-indo-only/m01/metadata_enriched.tsv', sep='\t')
    
    # Check the columns related to Amantadine resistance
    print("\nColumns in metadata:")
    print(df.columns.tolist()[:15]) # print first 15 columns
    
    if 'adamantanes_resistance_geno' in df.columns:
        print("\nValue counts for adamantanes_resistance_geno:")
        print(df['adamantanes_resistance_geno'].value_counts(dropna=False))
        
    if 'adamantanes_resistance_pheno' in df.columns:
        print("\nValue counts for adamantanes_resistance_pheno:")
        print(df['adamantanes_resistance_pheno'].value_counts(dropna=False))

if __name__ == '__main__':
    main()
