import os

def print_first_lines(label, path):
    if not os.path.exists(path):
        print(f"{label}: NOT FOUND")
        return
    with open(path, 'r') as f:
        lines = f.readlines()
    print(f"\n--- {label} ({path}) ---")
    # print up to 15 lines of the table
    count = 0
    for line in lines:
        if '|' in line:
            print(line.strip())
            count += 1
            if count >= 8:
                break

def main():
    runs = ['influenza-a-mp-indo-only', 'influenza-a-mp-sea-only', 'influenza-a-mp-all-partial']
    for r in runs:
        degen = f"ssh-output/results/{r}/m10/degenerate/coverage_report.md"
        non_degen = f"ssh-output/results/{r}/m10/non_degenerate/coverage_report.md"
        print_first_lines(f"{r} - DEGENERATE", degen)
        print_first_lines(f"{r} - NON-DEGENERATE", non_degen)

if __name__ == '__main__':
    main()
