import json

def main():
    print("Reading aa_muts.json...")
    with open('ssh-output/results/influenza-a-mp-indo-only/m06/aa_muts.json', 'r') as f:
        data = json.load(f)

    # Let's count how many times each M1 amino acid position appears as mutated in the tree
    m1_counts = {}
    nodes = data.get("nodes", {})
    for node_id, node_data in nodes.items():
        aa_muts = node_data.get("aa_muts", {})
        m1_muts = aa_muts.get("M1", [])
        for mut in m1_muts:
            # mut is like 'N82T' or 'R210K'
            # let's extract the position
            pos_str = "".join([c for c in mut if c.isdigit()])
            if pos_str:
                pos = int(pos_str)
                m1_counts[pos] = m1_counts.get(pos, 0) + 1

    print("\nMost mutated positions in M1 across the phylogenetic tree:")
    sorted_muts = sorted(m1_counts.items(), key=lambda x: x[1], reverse=True)
    for pos, count in sorted_muts[:15]:
        print(f"AA Position {pos}: mutated {count} times in the tree")

if __name__ == '__main__':
    main()
