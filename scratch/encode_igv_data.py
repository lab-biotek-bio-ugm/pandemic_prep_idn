import base64

def get_base64_data_url(filepath):
    with open(filepath, 'rb') as f:
        data = f.read()
    b64 = base64.b64encode(data).decode('utf-8')
    return f"data:text/plain;base64,{b64}"

def main():
    print("Encoding files to base64...")
    fasta_url = get_base64_data_url('config/references/flu_a_mp_reference.fasta')
    fai_url = get_base64_data_url('config/references/flu_a_mp_reference.fasta.fai')
    entropy_url = get_base64_data_url('ssh-output/results/influenza-a-mp-indo-only/m09/degenerate/segment-MP_region_1/high_entropy.gff3')
    candidates_url = get_base64_data_url('ssh-output/results/influenza-a-mp-indo-only/m09/degenerate/segment-MP_region_1/final_candidates.gff3')

    qmd_path = 'report/m09.qmd'
    print(f"Reading {qmd_path}...")
    with open(qmd_path, 'r', encoding='utf-8') as f:
        content = f.read()

    print("Injecting base64 URLs into placeholders...")
    content = content.replace('"INJECT_FASTA_URL"', f'"{fasta_url}"')
    content = content.replace('"INJECT_FAI_URL"', f'"{fai_url}"')
    content = content.replace('"INJECT_HIGH_ENTROPY_URL"', f'"{entropy_url}"')
    content = content.replace('"INJECT_CANDIDATES_URL"', f'"{candidates_url}"')

    print(f"Writing updated content back to {qmd_path}...")
    with open(qmd_path, 'w', encoding='utf-8') as f:
        f.write(content)

    print("Done! Injection completed successfully.")

if __name__ == '__main__':
    main()
