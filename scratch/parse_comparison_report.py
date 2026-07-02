import os
import re

def parse_html_report(path):
    if not os.path.exists(path):
        return "Not found"
    
    with open(path, 'r', encoding='utf-8') as f:
        content = f.read()
        
    # Let's extract heading tags (h1, h2, h3) and key table headers/labels
    headings = re.findall(r'<h[1-3][^>]*>(.*?)</h[1-3]>', content)
    tables = re.findall(r'<table[^>]*>.*?</table>', content, re.DOTALL)
    
    print(f"\nReport Path: {path}")
    print("Headings:")
    for h in headings[:10]:
        # Strip HTML tags inside heading
        clean_h = re.sub(r'<[^>]+>', '', h)
        print(f"  - {clean_h.strip()}")
        
    print(f"Total tables found: {len(tables)}")
    if tables:
        # Let's print some text from the first table
        clean_t = re.sub(r'<[^>]+>', ' ', tables[0])
        words = clean_t.split()
        print("First table preview:")
        print(" ".join(words[:50]) + "...")

def main():
    runs = ['influenza-a-mp-indo-only', 'influenza-a-mp-sea-only', 'influenza-a-mp-all-partial']
    for r in runs:
        path = f"ssh-output/results/{r}/m10/comparison_report.html"
        parse_html_report(path)

if __name__ == '__main__':
    main()
