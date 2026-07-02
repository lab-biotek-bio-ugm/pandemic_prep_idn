import os
import re
import sys

if sys.platform.startswith('win'):
    sys.stdout.reconfigure(encoding='utf-8')

def main():
    path = "ssh-output/results/influenza-a-mp-indo-only/m10/comparison_report.html"
    if not os.path.exists(path):
        print("Not found")
        return
        
    with open(path, 'r', encoding='utf-8') as f:
        html = f.read()
        
    scripts = re.findall(r'<script[^>]*>(.*?)</script>', html, re.DOTALL)
    if scripts:
        print("Script tag content preview (first 3000 chars):")
        print(scripts[0][:3000])
    else:
        print("No script tags found")

if __name__ == '__main__':
    main()
