import os
import re
import sys
import json

if sys.platform.startswith('win'):
    sys.stdout.reconfigure(encoding='utf-8')

def main():
    path = "ssh-output/results/influenza-a-mp-indo-only/m10/comparison_report.html"
    if not os.path.exists(path):
        print("Not found")
        return
        
    with open(path, 'r', encoding='utf-8') as f:
        html = f.read()
        
    # Search for script tags
    scripts = re.findall(r'<script[^>]*>(.*?)</script>', html, re.DOTALL)
    print(f"Found {len(scripts)} script tags")
    
    for idx, s in enumerate(scripts):
        if 'const candidates =' in s or 'var candidates =' in s or 'candidates =' in s:
            print(f"\n--- Script {idx+1} contains candidates data! ---")
            # Let's print the variable definition
            match = re.search(r'(const|var)?\s*candidates\s*=\s*(.*?);', s, re.DOTALL)
            if match:
                data_str = match.group(2).strip()
                print(f"Data string length: {len(data_str)}")
                # Try to print first 500 characters of the JSON/array
                print("Data preview:")
                print(data_str[:1500])
        
        if 'const subtypeData =' in s or 'subtypeData =' in s:
            print(f"\n--- Script {idx+1} contains subtypeData! ---")
            match = re.search(r'(const|var)?\s*subtypeData\s*=\s*(.*?);', s, re.DOTALL)
            if match:
                data_str = match.group(2).strip()
                print(data_str[:1000])

if __name__ == '__main__':
    main()
