import os

def main():
    path = "ssh-output/results/influenza-a-mp-indo-only/m10/comparison_report.html"
    if os.path.exists(path):
        with open(path, 'r', encoding='utf-8') as f:
            content = f.read(15000) # Read first 15KB
        print("HTML HEAD PREVIEW:")
        print(content[:3000])
        print("\nSEARCHING FOR TABLES:")
        # Look for table tag
        idx = content.find('<table')
        if idx != -1:
            print(content[idx:idx+2500])
        else:
            print("No <table tag found in first 15KB")

if __name__ == '__main__':
    main()
