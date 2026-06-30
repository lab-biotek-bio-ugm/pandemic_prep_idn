import os
import subprocess
import sys

def main():
    chrome_path = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
    if not os.path.exists(chrome_path):
        print(f"Error: Chrome not found at {chrome_path}")
        sys.exit(1)
        
    dst_dir = "pdf-report"
    os.makedirs(dst_dir, exist_ok=True)
    
    pages = [
        ("results_phylodynamics.html", "results_phylodynamics.pdf"),
        ("results_selection.html", "results_selection.pdf"),
        ("results_candidates.html", "results_candidates.pdf"),
        ("results_performance.html", "results_performance.pdf")
    ]
    
    current_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    site_dir = os.path.join(current_dir, "report", "_site")
    
    for html_file, pdf_file in pages:
        html_path = os.path.join(site_dir, html_file)
        pdf_path = os.path.join(current_dir, dst_dir, pdf_file)
        
        if not os.path.exists(html_path):
            print(f"Warning: HTML source not found: {html_path}")
            continue
            
        url = f"file:///{html_path.replace(os.sep, '/')}"
        
        cmd = [
            chrome_path,
            "--headless",
            "--disable-gpu",
            "--no-sandbox",
            f"--print-to-pdf={pdf_path}",
            url
        ]
        
        print(f"Printing: {html_file} -> {pdf_file}...")
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            print(f"Successfully generated: {pdf_path}")
        except subprocess.CalledProcessError as e:
            print(f"Error printing {html_file}: {e}")
            print(f"Stdout: {e.stdout}")
            print(f"Stderr: {e.stderr}")

if __name__ == "__main__":
    main()
