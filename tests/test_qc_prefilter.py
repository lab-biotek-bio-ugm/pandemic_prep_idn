"""Regression test for M02 qc_prefilter column detection (scripts/qc_prefilter.py).

Guards the fix for the dead subtype-scarcity filter: after sanitize_metadata.py the
ID column is 'strain' and the subtype column is lowercase 'subtype'. The pre-fix code
looked for capitalised 'Subtype'/'Segment_Id' and silently matched nothing, so the
documented rare-subtype filter never removed anything.
"""
import json
import subprocess
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parent.parent


def test_subtype_scarcity_filter_runs_on_sanitized_schema(tmp_path):
    meta = tmp_path / "meta.tsv"
    meta.write_text(
        "strain\tsubtype\tdate\n"
        "EPI1\tH3N2\t2020\n"
        "EPI2\tH3N2\t2020\n"
        "EPI3\tH1N1\t2021\n"  # H1N1 is scarce (n=1)
    )
    fasta = tmp_path / "seqs.fasta"
    fasta.write_text(
        ">EPI1\nACGTACGTACGT\n"
        ">EPI2\nACGTACGTACGA\n"
        ">EPI3\nTTTTGGGGCCCC\n"
    )
    out = tmp_path / "out.fasta"
    report = tmp_path / "rep.json"

    res = subprocess.run(
        [sys.executable, str(_REPO / "scripts" / "qc_prefilter.py"),
         "--input-fasta", str(fasta), "--output-fasta", str(out),
         "--report", str(report), "--metadata", str(meta),
         "--min-subtype-count", "2", "--max-n-pct", "5"],
        capture_output=True, text=True,
    )
    assert res.returncode == 0, res.stderr

    rep = json.loads(report.read_text())
    # The scarce H1N1 sequence (EPI3) must be removed by the subtype filter.
    assert rep["subtype_filtering"]["removed_count"] == 1
    assert rep["output_count"] == 2
    # Column auto-detection should have engaged (proves the filter isn't dead code).
    assert "Using ID column 'strain'" in res.stdout
