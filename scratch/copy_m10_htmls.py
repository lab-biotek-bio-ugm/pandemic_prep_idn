import os
import shutil

def main():
    src_base = "ssh-output/results"
    dst_base = "report/results"
    
    scenarios = [
        "influenza-a-mp-indo-only",
        "influenza-a-mp-sea-only",
        "influenza-a-mp-all-partial"
    ]
    
    files_to_copy = [
        "m10/comparison_report.html",
        "m10/degenerate/coverage_heatmap.html",
        "m10/non_degenerate/coverage_heatmap.html"
    ]
    
    for scenario in scenarios:
        for f in files_to_copy:
            src_path = os.path.join(src_base, scenario, f)
            dst_path = os.path.join(dst_base, scenario, f)
            
            if os.path.exists(src_path):
                # Ensure destination directory exists
                os.makedirs(os.path.dirname(dst_path), exist_ok=True)
                shutil.copy2(src_path, dst_path)
                print(f"Copied: {src_path} -> {dst_path}")
            else:
                print(f"Warning: Source not found: {src_path}")

if __name__ == "__main__":
    main()
