import argparse
import yaml
from pathlib import Path

def generate_payload(manifest_path, template_path, out_dir):
    # 1. Baca manifest M07
    with open(manifest_path, 'r') as f:
        manifest = yaml.safe_load(f)
        
    # 2. Baca template PrimerLab sebagai TEKS MURNI untuk menjaga komentar dan format
    with open(template_path, 'r') as f:
        template_lines = f.readlines()
        
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # 3. Ekstrak setiap region kandidat
    for idx, region in enumerate(manifest.get('candidate_regions', [])):
        region_id = region['region_id']
        meme_pvals = region['positional_scores']['meme_pval']
        consensus_seq = region['sequences']['consensus_sequence']
        
        # Hitung excluded_regions
        excluded_list = []
        in_exclusion = False
        start_idx = 0
        
        for i, pval in enumerate(meme_pvals):
            if pval < 0.10:
                if not in_exclusion:
                    in_exclusion = True
                    start_idx = i
            else:
                if in_exclusion:
                    excluded_list.append(f"{{start: {start_idx}, length: {i - start_idx}}}")
                    in_exclusion = False
                    
        if in_exclusion:
            excluded_list.append(f"{{start: {start_idx}, length: {len(meme_pvals) - start_idx}}}")
            
        # Format menjadi string array YAML satu baris: [{start: X, length: Y}, ...]
        excluded_str = "[" + ", ".join(excluded_list) + "]"
        
        # 4. Cari dan Ganti Teks (Text Replacement)
        out_lines = []
        for line in template_lines:
            # Ganti sequence: null
            if line.startswith("  sequence: null"):
                out_lines.append(f"  sequence: {consensus_seq}             # [M08 INJECTED]\n")
            # Ganti type: auto
            elif line.startswith("  type: auto"):
                out_lines.append("  type: raw                 # [M08 INJECTED]\n")
            # Ganti excluded_regions: []
            elif line.startswith("  excluded_regions: []"):
                # Pertahankan komentar asli jika ada di sebelah kanan "[]"
                if "[]" in line:
                    left, right = line.split("[]", 1)
                    out_lines.append(f"{left}{excluded_str}{right}")
                else:
                    out_lines.append(f"  excluded_regions: {excluded_str}\n")
            else:
                out_lines.append(line)
                
        # 5. Tambahkan Metadata M07 di bagian paling bawah
        out_lines.append("\n# ==============================================================================\n")
        out_lines.append("# --- M07 EVOLUTIONARY METADATA (DARI PANDEMIC PIPELINE) ---\n")
        out_lines.append("# ==============================================================================\n")
        out_lines.append("metadata:\n")
        out_lines.append(f"  project_name: {region_id}\n")
        out_lines.append("m07_metadata:\n")
        out_lines.append(f"  region_id: {region_id}\n")
        out_lines.append(f"  composite_score: {region['composite_score']}\n")
        out_lines.append(f"  iupac_density: {region['sequences']['iupac_density']}\n")
        out_lines.append(f"  tier: {region['tier']}\n")
        
        # 6. Simpan Payload YAML
        out_path = out_dir / f"primerlab_config_{region_id}.yaml"
        with open(out_path, 'w') as f:
            f.writelines(out_lines)
            
        print(f"[M08] Generated PrimerLab YAML config: {out_path}")
        print(f"      - Target Region : {region_id}")
        print(f"      - Excluded Zones: {len(excluded_list)} blocks (MEME p-val < 0.1)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="M08 Payload Generator for PrimerLab")
    parser.add_argument("--manifest", required=True, help="Path to M07 primer_targets_manifest.yaml")
    parser.add_argument("--template", required=True, help="Path to PrimerLab raa_default.yaml")
    parser.add_argument("--outdir", required=True, help="Output directory for generated YAMLs")
    args = parser.parse_args()
    
    generate_payload(args.manifest, args.template, args.outdir)
