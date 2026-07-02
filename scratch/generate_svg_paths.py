import json

def get_svg_path(filepath, total_width=500, x_offset=60, y_min=180, y_max=40):
    with open(filepath, 'r') as f:
        data = json.load(f)
    content = data['MLE']['content']['0']
    total_codons = len(content)
    
    path_coords = []
    for idx, site_row in enumerate(content):
        # Index 3 is Prob[alpha>beta] (negative selection probability)
        prob = site_row[3]
        # Map X: from 0 to total_width
        x = x_offset + (idx / (total_codons - 1)) * total_width
        # Map Y: from y_min (0.0) to y_max (1.0)
        y = y_min - prob * (y_min - y_max)
        path_coords.append(f"{x:.1f},{y:.1f}")
        
    path_str = "M " + " L ".join(path_coords)
    return path_str, total_codons

def main():
    print("=== Generating M1 FUBAR Path ===")
    m1_path, m1_len = get_svg_path('ssh-output/results/influenza-a-mp-indo-only/m05/fubar_M1.json')
    print(f"M1 Codons: {m1_len}")
    print(m1_path[:200] + " ... " + m1_path[-200:])
    
    print("\n=== Generating M2 FUBAR Path ===")
    m2_path, m2_len = get_svg_path('ssh-output/results/influenza-a-mp-indo-only/m05/fubar_M2.json')
    print(f"M2 Codons: {m2_len}")
    print(m2_path[:200] + " ... " + m2_path[-200:])

    # Let's save them to a temporary file so we can copy them
    with open('scratch/svg_paths.txt', 'w') as f:
        f.write("M1 PATH:\n")
        f.write(m1_path + "\n\n")
        f.write("M2 PATH:\n")
        f.write(m2_path + "\n")

if __name__ == '__main__':
    main()
