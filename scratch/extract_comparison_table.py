import os
import re

def extract_table_rows(html_path):
    if not os.path.exists(html_path):
        return "Not found"
    
    with open(html_path, 'r', encoding='utf-8') as f:
        html = f.read()
        
    # Find tables
    tables = re.findall(r'<table[^>]*>.*?</table>', html, re.DOTALL)
    if not tables:
        return "No tables found"
        
    # Let's get the first table
    table = tables[0]
    
    # Extract rows (tr)
    rows = re.findall(r'<tr[^>]*>.*?</tr>', table, re.DOTALL)
    
    clean_rows = []
    for r in rows:
        # Find cells (th or td)
        cells = re.findall(r'<t[d|h][^>]*>(.*?)</t[d|h]>', r, re.DOTALL)
        clean_cells = [re.sub(r'<[^>]+>', '', c).strip() for c in cells]
        if clean_cells:
            clean_rows.append(clean_cells)
            
    return clean_rows

def main():
    runs = ['influenza-a-mp-indo-only', 'influenza-a-mp-sea-only', 'influenza-a-mp-all-partial']
    for r in runs:
        path = f"ssh-output/results/{r}/m10/comparison_report.html"
        print(f"\n=== Table in {path} ===")
        rows = extract_table_rows(path)
        if isinstance(rows, str):
            print(rows)
        else:
            for r_idx, row in enumerate(rows[:6]): # Print first 6 rows (header + top 5 rows)
                print(f"Row {r_idx}: {row}")

if __name__ == '__main__':
    main()
