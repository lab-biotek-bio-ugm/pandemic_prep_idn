import argparse
import yaml
from pathlib import Path

def generate_payload(manifest_path, template_path, out_dir, args=None):
    # 1. Baca manifest M07
    with open(manifest_path, 'r') as f:
        manifest = yaml.safe_load(f)
        
    # 2. Baca template PrimerLab sebagai TEKS MURNI untuk menjaga komentar dan format
    with open(template_path, 'r') as f:
        template_lines = f.readlines()
        
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Optional parameters from CLI (args is passed explicitly to avoid NameError)
    num_candidates = getattr(args, "num_candidates", None) if args else None
    vienna_limit = getattr(args, "vienna_limit", None) if args else None
    window_size = getattr(args, "window_size", None) if args else None
    overlap = getattr(args, "overlap", None) if args else None
    search_strategy = getattr(args, "search_strategy", None) if args else None
    hard_qc_filter = getattr(args, "hard_qc_filter", None) if args else None
    max_search_iterations = getattr(args, "max_search_iterations", None) if args else None
    validate = getattr(args, "validate", None) if args else None
    blast_check = getattr(args, "blast_check", None) if args else None
    blast_limit = getattr(args, "blast_limit", None) if args else None
    blast_db = getattr(args, "blast_db", None) if args else None
    blast_mode = getattr(args, "blast_mode", None) if args else None
    blast_local_db_path = getattr(args, "blast_local_db_path", None) if args else None
    blast_targeted_check = getattr(args, "blast_targeted_check", None) if args else None
    blast_targeted_db_path = getattr(args, "blast_targeted_db_path", None) if args else None
    entropy_threshold = getattr(args, "entropy_threshold", 0.50) if args else 0.50
    # Minimum consecutive excluded positions to form a block — prevents hundreds of tiny
    # 1-2bp wobble-position blocks from blocking all primer placement in Primer3.
    # Default 10: only continuous variable stretches ≥10bp (like UTRs) are excluded.
    min_exclusion_block_length = getattr(args, "min_exclusion_block_length", 10) if args else 10

    candidate_regions = manifest.get('candidate_regions', [])
    if not candidate_regions:
        print("[M08] WARNING: No candidate regions found in manifest.")
        return

    # 3. Analisis strategi manifest
    lengths = [len(r['sequences']['consensus_sequence']) for r in candidate_regions]
    max_len = max(lengths)
    min_len = min(lengths)
    is_full_segment_strategy = (max_len - min_len > 100) or (max_len > 500)
    
    if is_full_segment_strategy:
        # Cari sekuens konsensus terpanjang sebagai representasi full segment
        full_segment_region = max(candidate_regions, key=lambda r: len(r['sequences']['consensus_sequence']))
        full_segment_seq = full_segment_region['sequences']['consensus_sequence']
        print(f"[M08] Detected full-segment strategy. Base template sequence length: {len(full_segment_seq)} bp")
    else:
        full_segment_seq = None
        print(f"[M08] Detected windowed strategy. Slices will be designed independently.")

    # 4. Ekstrak setiap region kandidat
    for idx, region in enumerate(candidate_regions):
        region_id = region['region_id']
        meme_pvals = region['positional_scores']['meme_pval']
        cz = region.get('conserved_zone', {})
        cz_start = cz.get('start', 0)
        
        # Tentukan template sekuens
        if is_full_segment_strategy and full_segment_seq:
            consensus_seq = full_segment_seq
        else:
            consensus_seq = region['sequences']['consensus_sequence']
            
        target_region_str = "null"
        included_region_str = "null"
            
        entropy_scores = region.get('positional_scores', {}).get('entropy', [])
            
        # Hitung excluded_regions (relatif terhadap template sekuens)
        # Gabungkan meme_pval < 0.10 dan entropy >= entropy_threshold
        excluded_list = []
        in_exclusion = False
        start_idx = 0
        
        def check_exclude(pos):
            is_meme_excl = (pos < len(meme_pvals)) and (meme_pvals[pos] is not None) and (meme_pvals[pos] < 0.10)
            is_ent_excl = (pos < len(entropy_scores)) and (entropy_scores[pos] is not None) and (entropy_scores[pos] >= entropy_threshold)
            return is_meme_excl or is_ent_excl
            
        total_len = max(len(meme_pvals), len(entropy_scores))
        
        for i in range(total_len):
            if check_exclude(i):
                if not in_exclusion:
                    in_exclusion = True
                    start_idx = i
            else:
                if in_exclusion:
                    block_len = i - start_idx
                    if block_len >= min_exclusion_block_length:
                        offset = cz_start if is_full_segment_strategy and full_segment_seq else 0
                        excluded_list.append(f"{{start: {start_idx + offset}, length: {block_len}}}")
                    in_exclusion = False
                    
        if in_exclusion:
            block_len = total_len - start_idx
            if block_len >= min_exclusion_block_length:
                offset = cz_start if is_full_segment_strategy and full_segment_seq else 0
                excluded_list.append(f"{{start: {start_idx + offset}, length: {block_len}}}")
            
        # Format menjadi string array YAML satu baris: [{start: X, length: Y}, ...]
        excluded_str = "[" + ", ".join(excluded_list) + "]"
        
        # 5. Cari dan Ganti Teks (Text Replacement)
        out_lines = []
        for line in template_lines:
            # Ganti sequence: null
            if line.startswith("  sequence: null"):
                out_lines.append(f"  sequence: {consensus_seq}             # [M08 INJECTED]\n")
            # Ganti type: auto
            elif line.startswith("  type: auto"):
                out_lines.append("  type: raw                 # [M08 INJECTED]\n")
            # Ganti target_region: null
            elif line.strip().startswith("target_region:"):
                comment = line.split("#", 1)[1] if "#" in line else ""
                indent = line[:line.find("target_region")]
                out_lines.append(f"{indent}target_region: {target_region_str}       # [M08 INJECTED]{' #' + comment if comment else ''}\n")
            # Ganti included_region: null
            elif line.strip().startswith("included_region:"):
                comment = line.split("#", 1)[1] if "#" in line else ""
                indent = line[:line.find("included_region")]
                out_lines.append(f"{indent}included_region: {included_region_str}       # [M08 INJECTED]{' #' + comment if comment else ''}\n")
            # Ganti excluded_regions
            elif line.strip().startswith("excluded_regions:"):
                # Pertahankan komentar asli jika ada
                comment = line.split("#", 1)[1] if "#" in line else ""
                indent = line[:line.find("excluded_regions")]
                out_lines.append(f"{indent}excluded_regions: {excluded_str}       # [M08 INJECTED]{' #' + comment if comment else ''}\n")
            
            # --- New Overrides ---
            elif num_candidates is not None and line.strip().startswith("num_candidates:"):
                comment = line.split("#", 1)[1] if "#" in line else ""
                out_lines.append(f"  num_candidates: {num_candidates}       # [M08 OVERRIDE]{'#' + comment if comment else ''}\n")
            elif search_strategy is not None and line.strip().startswith("search_strategy:"):
                comment = line.split("#", 1)[1] if "#" in line else ""
                out_lines.append(f"  search_strategy: {search_strategy}  # [M08 OVERRIDE]{'#' + comment if comment else ''}\n")
            elif hard_qc_filter is not None and line.strip().startswith("hard_qc_filter:"):
                comment = line.split("#", 1)[1] if "#" in line else ""
                out_lines.append(f"  hard_qc_filter: {str(hard_qc_filter).lower()}       # [M08 OVERRIDE]{'#' + comment if comment else ''}\n")
            elif max_search_iterations is not None and line.strip().startswith("max_search_iterations:"):
                comment = line.split("#", 1)[1] if "#" in line else ""
                out_lines.append(f"  max_search_iterations: {max_search_iterations} # [M08 OVERRIDE]{'#' + comment if comment else ''}\n")
            elif vienna_limit is not None and line.strip().startswith("vienna_ranking_limit:"):
                comment = line.split("#", 1)[1] if "#" in line else ""
                out_lines.append(f"  vienna_ranking_limit: {vienna_limit}   # [M08 OVERRIDE]{'#' + comment if comment else ''}\n")
            elif vienna_limit is not None and line.strip().startswith("num_results_to_export:"):
                comment = line.split("#", 1)[1] if "#" in line else ""
                out_lines.append(f"  num_results_to_export: {vienna_limit}   # [M08 OVERRIDE]{'#' + comment if comment else ''}\n")
            elif blast_check is not None and line.strip().startswith("blast_offtarget_check:"):
                comment = line.split("#", 1)[1] if "#" in line else ""
                out_lines.append(f"  blast_offtarget_check: {str(blast_check).lower()} # [M08 OVERRIDE]{'#' + comment if comment else ''}\n")
            elif blast_limit is not None and line.strip().startswith("blast_ranking_limit:"):
                comment = line.split("#", 1)[1] if "#" in line else ""
                out_lines.append(f"  blast_ranking_limit: {blast_limit}     # [M08 OVERRIDE]{'#' + comment if comment else ''}\n")
            elif blast_db is not None and line.strip().startswith("blast_database:"):
                comment = line.split("#", 1)[1] if "#" in line else ""
                out_lines.append(f"  blast_database: \"{blast_db}\"         # [M08 OVERRIDE]{'#' + comment if comment else ''}\n")
                if blast_mode is not None:
                    out_lines.append(f"  blast_mode: \"{blast_mode}\"\n")
                if blast_local_db_path is not None:
                    out_lines.append(f"  blast_local_db_path: \"{blast_local_db_path}\"\n")
                if blast_targeted_check is not None:
                    out_lines.append(f"  blast_targeted_check: {str(blast_targeted_check).lower()}\n")
                if blast_targeted_db_path is not None:
                    out_lines.append(f"  blast_targeted_db_path: \"{blast_targeted_db_path}\"\n")
            elif validate is not None and line.strip().startswith("auto_validate:"):
                comment = line.split("#", 1)[1] if "#" in line else ""
                out_lines.append(f"  auto_validate: {str(validate).lower()} # [M08 OVERRIDE]{'#' + comment if comment else ''}\n")
            elif window_size is not None and line.strip().startswith("# window_size:"):
                # Uncomment and override
                out_lines.append(f"  window_size: {window_size}             # [M08 OVERRIDE]\n")
            elif overlap is not None and line.strip().startswith("# overlap:"):
                # Uncomment and override
                out_lines.append(f"  overlap: {overlap}                 # [M08 OVERRIDE]\n")
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
        print(f"      - Excluded Zones: {len(excluded_list)} blocks (entropy>={entropy_threshold} OR meme_pval<0.1, min_block_len={min_exclusion_block_length})")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="M08 Payload Generator for PrimerLab")
    parser.add_argument("--manifest", required=True, help="Path to M07 primer_targets_manifest.yaml")
    parser.add_argument("--template", required=True, help="Path to PrimerLab raa_default.yaml")
    parser.add_argument("--outdir", required=True, help="Output directory for generated YAMLs")
    
    # Optional Overrides
    parser.add_argument("--num-candidates", type=int, help="Override number of candidates")
    parser.add_argument("--vienna-limit", type=int, help="Override ViennaRNA ranking limit")
    parser.add_argument("--window-size", type=int, help="Enable and override window size")
    parser.add_argument("--overlap", type=int, help="Override window overlap")
    parser.add_argument("--search-strategy", type=str, help="Override search strategy (e.g. standard, diversified)")
    parser.add_argument("--hard-qc-filter", type=str, help="Override hard qc filter (true/false)")
    parser.add_argument("--max-search-iterations", type=int, help="Override max search iterations")
    parser.add_argument("--validate", type=str, help="Override in-silico validation (true/false)")
    parser.add_argument("--blast-check", type=str, help="Override blast off-target check (true/false)")
    parser.add_argument("--blast-limit", type=int, help="Override blast ranking limit")
    parser.add_argument("--blast-db", type=str, help="Override blast database")
    parser.add_argument("--blast-mode", type=str, help="Override blast mode (remote/local)")
    parser.add_argument("--blast-local-db-path", type=str, help="Override local blast database path")
    parser.add_argument("--blast-targeted-check", type=str, help="Override blast targeted check (true/false)")
    parser.add_argument("--blast-targeted-db-path", type=str, help="Override targeted blast database path")
    parser.add_argument("--entropy-threshold", type=float, default=0.50, help="Entropy threshold for exclusion (Shannon entropy)")
    parser.add_argument("--min-exclusion-block-length", type=int, default=10, help="Minimum consecutive excluded positions to form a block (prevents tiny wobble-position blocks from blocking all primers)")
    
    args = parser.parse_args()
    
    generate_payload(args.manifest, args.template, args.outdir, args)
