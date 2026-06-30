import os
import sys

# Reconfigure stdout to use utf-8 to prevent encoding errors on Windows terminal
if sys.platform.startswith('win'):
    sys.stdout.reconfigure(encoding='utf-8')

def main():
    path = "ssh-output/results/influenza-a-mp-indo-only/m10/comparison_report.html"
    if not os.path.exists(path):
        print("Not found")
        return
        
    with open(path, 'r', encoding='utf-8') as f:
        html = f.read()
        
    print(f"Length of HTML: {len(html)} characters")
    
    # Let's search for table or JSON inside
    # Some Quarto/Pandas reports embed JSON data or script tags.
    # Let's see if there is table or if it uses javascript to render tables.
    table_count = html.count('<table')
    print(f"Total '<table' instances: {table_count}")
    
    # Find all table matches and show start
    import re
    tables = re.findall(r'<table[^>]*>.*?</table>', html, re.DOTALL)
    for i, t in enumerate(tables):
        print(f"\n--- TABLE {i+1} PREVIEW ---")
        # print first 500 characters of table
        print(t[:500].strip())

if __name__ == '__main__':
    main()
