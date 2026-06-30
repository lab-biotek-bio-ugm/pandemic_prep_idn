import os

def parse_coverage_report(report_path):
    if not os.path.exists(report_path):
        return "Not found"
    
    with open(report_path, 'r') as f:
        content = f.read()
        
    # We want to extract the first table or key stats.
    # Usually has a line like "| Overall |" or mentions coverage percentage.
    lines = content.split('\n')
    stats = []
    for line in lines[:30]:
        if '|' in line:
            stats.append(line)
    return "\n".join(stats)

def main():
    base_dir = 'ssh-output/results'
    runs = ['influenza-a-mp-indo-only', 'influenza-a-mp-sea-only', 'influenza-a-mp-all-partial']
    for r in runs:
        print(f"\n==========================================")
        print(f"RUN: {r}")
        print(f"==========================================")
        
        # Degenerate
        degen_path = os.path.join(base_dir, r, 'm10', 'degenerate', 'coverage_report.md')
        print("--- DEGENERATE COVERAGE ---")
        print(parse_coverage_report(degen_path))
        
        # Non-Degenerate
        non_degen_path = os.path.join(base_dir, r, 'm10', 'non_degenerate', 'coverage_report.md')
        print("--- NON-DEGENERATE COVERAGE ---")
        print(parse_coverage_report(non_degen_path))

if __name__ == '__main__':
    main()
