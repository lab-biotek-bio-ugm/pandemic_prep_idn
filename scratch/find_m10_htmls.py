import os

def main():
    base_dir = 'ssh-output/results'
    for root, dirs, files in os.walk(base_dir):
        for f in files:
            if f.endswith('.html') and 'm10' in root:
                full_path = os.path.join(root, f)
                print(f"Found HTML: {full_path} (size: {os.path.getsize(full_path)} bytes)")

if __name__ == '__main__':
    main()
