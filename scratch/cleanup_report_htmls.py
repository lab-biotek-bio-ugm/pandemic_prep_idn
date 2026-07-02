import os
import glob
import shutil

def main():
    # Remove untracked html files directly in report/
    html_files = glob.glob('report/*.html')
    for f in html_files:
        if os.path.exists(f):
            os.remove(f)
            print(f"Removed untracked file: {f}")
            
    # Also remove report/site_libs if it exists
    libs_dir = 'report/site_libs'
    if os.path.exists(libs_dir):
        shutil.rmtree(libs_dir)
        print(f"Removed untracked directory: {libs_dir}")

if __name__ == '__main__':
    main()
