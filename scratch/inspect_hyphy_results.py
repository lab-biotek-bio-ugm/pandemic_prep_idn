import json

def parse_fubar(filepath):
    with open(filepath, 'r') as f:
        data = json.load(f)
    content = data['MLE']['content']['0']
    
    purifying_sites = []
    for idx, site_row in enumerate(content):
        # Index 3 is Prob[alpha>beta] (negative selection probability)
        prob = site_row[3]
        if prob >= 0.90:
            purifying_sites.append((idx + 1, prob))
            
    return len(content), purifying_sites

def parse_meme(filepath):
    with open(filepath, 'r') as f:
        data = json.load(f)
    content = data['MLE']['content']['0']
    
    positive_sites = []
    for idx, site_row in enumerate(content):
        # Index 6 is p-value
        pval = site_row[6]
        if pval <= 0.05:
            positive_sites.append((idx + 1, pval))
            
    return len(content), positive_sites

def main():
    print("=== M1 (Matrix Protein 1) ===")
    m1_len, m1_fubar = parse_fubar('ssh-output/results/influenza-a-mp-indo-only/m05/fubar_M1.json')
    print(f"FUBAR Purifying Sites (prob >= 0.90): {len(m1_fubar)} / {m1_len} sites ({len(m1_fubar)/m1_len*100:.1f}%)")
    
    m1_meme_len, m1_meme = parse_meme('ssh-output/results/influenza-a-mp-indo-only/m05/meme_M1.json')
    print(f"MEME Positive Sites (p <= 0.05): {len(m1_meme)} sites:")
    for site, pval in m1_meme:
        print(f"  Codon {site}: p-value = {pval:.4f}")
        
    print("\n=== M2 (Matrix Protein 2) ===")
    m2_len, m2_fubar = parse_fubar('ssh-output/results/influenza-a-mp-indo-only/m05/fubar_M2.json')
    print(f"FUBAR Purifying Sites (prob >= 0.90): {len(m2_fubar)} / {m2_len} sites ({len(m2_fubar)/m2_len*100:.1f}%)")
    
    m2_meme_len, m2_meme = parse_meme('ssh-output/results/influenza-a-mp-indo-only/m05/meme_M2.json')
    print(f"MEME Positive Sites (p <= 0.05): {len(m2_meme)} sites:")
    for site, pval in m2_meme:
        print(f"  Codon {site}: p-value = {pval:.4f}")

if __name__ == '__main__':
    main()
