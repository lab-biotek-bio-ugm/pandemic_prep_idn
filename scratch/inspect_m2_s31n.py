import json

def main():
    print("Reading aa_muts.json to verify M2 codon 31...")
    with open('ssh-output/results/influenza-a-mp-indo-only/m06/aa_muts.json', 'r') as f:
        data = json.load(f)

    # Search for any mutation at position 31 of M2
    m2_mut_counts = {}
    nodes = data.get("nodes", {})
    for node_id, node_data in nodes.items():
        aa_muts = node_data.get("aa_muts", {})
        # Note: let's check what genes are annotated in the JSON first
        for gene_name, muts in aa_muts.items():
            for mut in muts:
                if "31" in mut:
                    print(f"Node {node_id} has mutation {mut} in gene {gene_name}")

if __name__ == '__main__':
    main()
