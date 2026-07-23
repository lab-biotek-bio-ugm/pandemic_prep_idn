"""Guards against bugs that only appear when a script actually RUNS.

The rest of the suite imports modules and calls individual functions. That leaves
two blind spots, and both bit in production:

  * `m07_bridge_builder.main()` referenced an undefined name. Importing the module
    never executes main(), and CI's `compileall` only catches syntax errors, so it
    passed everything and then failed on the cluster after ~15 minutes of HyPhy.
  * `_generate_individual_forbidden_tracks()` in M09 used an undefined parameter
    inside a `try/except Exception`, so the NameError was swallowed and the
    entropy/MEME/FUBAR/GARD tracks were silently never written.
"""
import json
import random
import subprocess
import sys

import pytest

yaml = pytest.importorskip("yaml")


def test_no_undefined_names_in_any_script():
    """Static undefined-name check across every script.

    `snakemake` is excluded: Snakemake injects it as a global into scripts used
    via the `script:` directive, so pyflakes cannot see it.
    """
    pyflakes = pytest.importorskip("pyflakes")  # noqa: F841
    import glob

    files = sorted(glob.glob("scripts/*.py") + glob.glob("scripts/tools/*.py"))
    out = subprocess.run([sys.executable, "-m", "pyflakes", *files],
                         capture_output=True, text=True).stdout
    offenders = [l for l in out.splitlines()
                 if "undefined name" in l and "'snakemake'" not in l]
    assert not offenders, "undefined names found:\n" + "\n".join(offenders)


@pytest.fixture
def m07_inputs(tmp_path):
    """Minimal but structurally real inputs for a full M07 run."""
    random.seed(3)
    gene = "".join(random.choice("ACGT") for _ in range(90))
    ref = "A" * 30 + gene + "T" * 30

    (tmp_path / "ref.fasta").write_text(">NC_TEST\n" + ref + "\n")
    recs = []
    for i in range(12):
        s = list(ref)
        for p in random.sample(range(len(ref)), 3):
            s[p] = random.choice("ACGT")
        recs.append(f">seq{i}|2024-01-01\n" + "".join(s))
    (tmp_path / "aln.fasta").write_text("\n".join(recs) + "\n")

    # GFF3 is 1-based inclusive; the CDS spans the 90 nt gene.
    (tmp_path / "ann.gff").write_text(
        "##gff-version 3\n"
        "NC_TEST\ttest\tCDS\t31\t120\t.\t+\t0\tID=cds-M1;gene=M1\n")
    (tmp_path / "meta.tsv").write_text(
        "strain\tdate\tsubtype\tcountry\n"
        + "\n".join(f"seq{i}|2024-01-01\t2024-01-01\tH3N2\tIndonesia"
                    for i in range(12)) + "\n")

    hyphy = tmp_path / "hyphy"
    hyphy.mkdir()
    headers = [["alpha", ""], ["beta", ""], ["alpha-beta", ""],
               ["Prob[alpha>beta]", ""], ["Prob[alpha<beta]", ""],
               ["BayesFactor[alpha<beta]", ""]]
    rows = [[round(random.uniform(0.05, 3), 3), round(random.uniform(0.01, 1), 3),
             0.0, round(random.uniform(0.5, 0.99), 3), 0.01, 0.2]
            for _ in range(30)]
    (hyphy / "fubar_M1.json").write_text(
        json.dumps({"MLE": {"headers": headers, "content": {"0": rows}}}))
    (hyphy / "meme_M1.json").write_text(json.dumps({"MLE": {
        "headers": [["alpha", ""], ["beta-", ""], ["p-value", ""]],
        "content": {"0": [[1, 1, round(random.uniform(0.01, 0.9), 3)]
                          for _ in range(30)]}}}))

    cfg = yaml.safe_load(open("config/config.yaml"))
    (tmp_path / "active.yaml").write_text(yaml.safe_dump(cfg))
    (tmp_path / "out").mkdir()
    return tmp_path


def _run_m07(d):
    return subprocess.run(
        [sys.executable, "scripts/m07_bridge_builder.py",
         "--msa", str(d / "aln.fasta"), "--gff", str(d / "ann.gff"),
         "--hyphy-dir", str(d / "hyphy"), "--reference", str(d / "ref.fasta"),
         "--metadata", str(d / "meta.tsv"), "--outdir", str(d / "out"),
         "--config", str(d / "active.yaml"), "--build-name", "test",
         "--strategy", "full_segment", "--min-subtype-count", "30",
         "--iupac-threshold", "0.15"],
        capture_output=True, text=True)


def test_m07_runs_end_to_end(m07_inputs):
    r = _run_m07(m07_inputs)
    assert r.returncode == 0, f"M07 failed:\n{r.stdout}\n{r.stderr}"
    for name in ("conserved_regions.bed", "bridge_report.json",
                 "primer_targets_manifest.yaml",
                 "candidate_targets_summary.csv"):
        assert (m07_inputs / "out" / name).exists(), f"{name} was not written"


def test_m07_populates_the_durability_fields(m07_inputs):
    """The whole point of parsing FUBAR's alpha/beta is that the rate reaches the
    report. If it silently stays empty the durability multiplier is a no-op and
    every region gets the same discount."""
    assert _run_m07(m07_inputs).returncode == 0
    region = json.loads(
        (m07_inputs / "out" / "bridge_report.json").read_text())["regions"][0]
    assert region["avg_subst_rate"] != "", "FUBAR rates never reached the report"
    assert region["rel_subst_rate"] != ""
    assert 0 < region["durability_multiplier"] <= 1.0

    header = (m07_inputs / "out" / "candidate_targets_summary.csv"
              ).read_text().splitlines()[0]
    for col in ("avg_subst_rate", "rel_subst_rate", "durability_multiplier"):
        assert col in header, f"{col} missing from the CSV summary"


def test_m07_refuses_a_reference_of_the_wrong_length(m07_inputs):
    """The coordinate invariant check must fire, not crash. It is the guard that
    stops M07 emitting silently wrong coordinates when the reference passed to it
    is not the one M03 aligned against."""
    (m07_inputs / "ref.fasta").write_text(">NC_TEST\n" + "A" * 50 + "\n")
    r = _run_m07(m07_inputs)
    assert r.returncode == 1
    assert "!= alignment length" in r.stderr, r.stderr
