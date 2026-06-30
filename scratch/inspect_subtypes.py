import pandas as pd

def main():
    print("Reading metadata_enriched.tsv to count subtypes...")
    df = pd.read_csv('ssh-output/results/influenza-a-mp-indo-only/m01/metadata_enriched.tsv', sep='\t')
    
    if 'subtype' in df.columns:
        print("\nSubtype value counts:")
        print(df['subtype'].value_counts(dropna=False))

if __name__ == '__main__':
    main()
