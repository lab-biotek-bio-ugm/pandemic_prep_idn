import os
import re

def check_file(path):
    with open(path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        
    errors = []
    for i in range(1, len(lines)):
        line = lines[i].strip()
        # check if line starts with a list bullet (* or - followed by spaces)
        if re.match(r'^[\*\-]\s+', line):
            # check the previous line
            prev_line = lines[i-1].strip()
            if prev_line:
                # skip if previous line is also a list item
                if re.match(r'^[\*\-]\s+', prev_line):
                    continue
                # skip if previous line starts with markdown block markers
                if prev_line.startswith('#') or prev_line.startswith('|') or prev_line.startswith(':::') or prev_line.startswith('---'):
                    continue
                # We have a list following a paragraph without a blank line!
                errors.append((i+1, prev_line, line))
    return errors

def main():
    report_dir = 'report'
    for f in os.listdir(report_dir):
        if f.endswith('.qmd'):
            path = os.path.join(report_dir, f)
            errors = check_file(path)
            if errors:
                print(f"\nIn file: {path}")
                for line_num, prev, curr in errors:
                    # Avoid emoji printing issues by encoding/decoding or clean printing
                    print(f"  Line {line_num}:")
                    print(f"    Prev: {prev.encode('ascii', errors='replace').decode('ascii')}")
                    print(f"    List: {curr.encode('ascii', errors='replace').decode('ascii')}")

if __name__ == '__main__':
    main()
