import os
import subprocess
import sys
from pypdf import PdfWriter

def main():
    chrome_path = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
    if not os.path.exists(chrome_path):
        print(f"Error: Chrome not found at {chrome_path}")
        sys.exit(1)
        
    dst_dir = "pdf-report"
    os.makedirs(dst_dir, exist_ok=True)
    
    # List of all 20 pages in the order they appear in the website sidebar
    pages = [
        "index.html",
        "literature.html",
        "eda_justifikasi.html",
        "eda_statistik.html",
        "config_global.html",
        "m01.html",
        "m02.html",
        "m03.html",
        "m04.html",
        "m05.html",
        "m06.html",
        "m07.html",
        "m08.html",
        "m09.html",
        "m10.html",
        "results_phylodynamics.html",
        "results_selection.html",
        "results_candidates.html",
        "results_performance.html",
        "references.html"
    ]
    
    current_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    site_dir = os.path.join(current_dir, "report", "_site")
    
    temp_pdfs = []
    
    # 1. Print all pages to temporary PDFs
    for idx, html_file in enumerate(pages):
        html_path = os.path.join(site_dir, html_file)
        temp_pdf_name = f"temp_{idx:02d}_{html_file.replace('.html', '.pdf')}"
        temp_pdf_path = os.path.join(current_dir, dst_dir, temp_pdf_name)
        
        if not os.path.exists(html_path):
            print(f"Warning: HTML source not found: {html_path}")
            continue
            
        url = f"file:///{html_path.replace(os.sep, '/')}"
        
        cmd = [
            chrome_path,
            "--headless",
            "--disable-gpu",
            "--no-sandbox",
            f"--print-to-pdf={temp_pdf_path}",
            url
        ]
        
        print(f"[{idx+1}/{len(pages)}] Printing: {html_file}...")
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            temp_pdfs.append(temp_pdf_path)
        except subprocess.CalledProcessError as e:
            print(f"Error printing {html_file}: {e}")
            
    # 2. Merge all temporary PDFs into a single file
    if not temp_pdfs:
        print("Error: No PDF pages were generated. Exiting.")
        sys.exit(1)
        
    output_pdf_path = os.path.join(current_dir, dst_dir, "pandemic_preparedness_report.pdf")
    print(f"\nMerging {len(temp_pdfs)} PDF pages into single file: {output_pdf_path}...")
    
    try:
        merger = PdfWriter()
        for pdf_path in temp_pdfs:
            merger.append(pdf_path)
        merger.write(output_pdf_path)
        merger.close()
        print("Successfully merged PDF!")
    except Exception as e:
        print(f"Error merging PDFs: {e}")
        sys.exit(1)
        
    # 3. Clean up temporary PDF files
    print("Cleaning up temporary page PDF files...")
    for pdf_path in temp_pdfs:
        try:
            os.remove(pdf_path)
        except OSError as e:
            print(f"Warning: Could not remove temporary file {pdf_path}: {e}")
            
    print("\nPDF Report generation completed successfully!")
    print(f"Final report saved to: {output_pdf_path}")

if __name__ == "__main__":
    main()

