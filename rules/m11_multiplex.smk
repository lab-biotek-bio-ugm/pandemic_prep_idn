# ==============================================================================
# MODULE 11: MULTIPLEX COMPATIBILITY SCREENING (PAIRWISE & TRIPLEX)
# ==============================================================================
# Evaluates cross-dimer thermodynamics, spectral fluorophore separation,
# and kinetic Tm alignment for:
#   - Pairwise (dual-target): m11_multiplex_compatibility — 2 targets, 6 oligos, 9 inter-assay pairs
#   - Triplex (three-target): m11_triplex_compatibility   — 3 targets, 9 oligos, 27 inter-assay pairs
#
# Both rules use the shared IUPAC expansion engine in scripts/tools/iupac_utils.py.
# The IUPAC expansion cap (max_variants) must match config.yaml:
#   degenerate.max_iupac_variants (default: 16)
#
# These rules are ORCHESTRATOR-DRIVEN: scripts/run_multiplex_batch.py runs each
# build's full pipeline (M01-M09) first, THEN screens compatibility. Because the
# per-build M09 outputs already exist at that point, the target segment is
# discovered from disk via an input function rather than carried as a wildcard.
# Encoding segment as a wildcard that only appears in `input` (not `output`) was a
# bug — Snakemake cannot infer it and raised WildcardError.
# ==============================================================================

import glob as _glob


# Build names never contain '/' or the pairwise '_x_' separator; constrain the
# wildcards so the pairwise and triplex output patterns cannot ambiguously match.
wildcard_constraints:
    build_a = r"[^/]+",
    build_b = r"[^/]+",
    build_c = r"[^/]+"


# Per-build detection fluorophore / display label for the multiplex rules.
# The main Snakefile loads ONE config, so cross-build dyes must be declared
# explicitly, e.g.:
#   multiplex:
#     dyes:   {influenza-a-mp-indo-only: FAM, influenza-b-mp-indo-only: HEX}
#     labels: {influenza-a-mp-indo-only: "Influenza A MP"}
# Unknown builds resolve to UNRESOLVED so check_dye_separation reports UNKNOWN
# rather than silently assuming the targets are spectrally separated.
_MUX_CFG = config.get("multiplex", {}) or {}
_DYE_MAP = _MUX_CFG.get("dyes", {}) or {}
_LABEL_MAP = _MUX_CFG.get("labels", {}) or {}


def _dye_for(build):
    return _DYE_MAP.get(build, "UNRESOLVED")


def _label_for(build):
    return _LABEL_MAP.get(build, build)


def _find_m09_candidates(build):
    """Resolve a build's degenerate M09 final_candidates.json without hardcoding the
    segment. Returns the on-disk path if present, else a deterministic expected path
    (whose absence yields a clear MissingInputException instead of a WildcardError)."""
    base = f"results/{build}/m09/degenerate"
    hits = sorted(_glob.glob(f"{base}/segment-*_region_1/final_candidates.json"))
    return hits[0] if hits else f"{base}/segment-UNKNOWN_region_1/final_candidates.json"


# ------------------------------------------------------------------------------
# PAIRWISE: Dual-Target Screening (e.g. Influenza A + Influenza B)
# ------------------------------------------------------------------------------
rule m11_multiplex_compatibility:
    input:
        target_a = lambda wc: _find_m09_candidates(wc.build_a),
        target_b = lambda wc: _find_m09_candidates(wc.build_b)
    output:
        # NOTE: multiplex_matrix.html was declared here but m11_multiplex_checker.py
        # never writes an HTML file, so the rule always died with
        # MissingOutputException. Removed rather than faked.
        json   = "results/multiplex-{build_a}_x_{build_b}/multiplex_candidates.json",
        report = "results/multiplex-{build_a}_x_{build_b}/multiplex_report.md"
    params:
        outdir         = "results/multiplex-{build_a}_x_{build_b}",
        num_candidates = 10,
        max_variants   = config.get("primerlab", {}).get("design", {}).get("degenerate", {}).get("max_iupac_variants", 16),
        dye_a          = lambda wc: _dye_for(wc.build_a),
        dye_b          = lambda wc: _dye_for(wc.build_b),
        label_a        = lambda wc: _label_for(wc.build_a),
        label_b        = lambda wc: _label_for(wc.build_b)
    log:
        "logs/m11_multiplex_{build_a}_x_{build_b}.log"
    shell:
        """
        python scripts/m11_multiplex_checker.py \
            --file-a {input.target_a} \
            --file-b {input.target_b} \
            --outdir {params.outdir} \
            --num-candidates {params.num_candidates} \
            --max-variants {params.max_variants} \
            --dye-a {params.dye_a} --dye-b {params.dye_b} \
            --label-a "{params.label_a}" --label-b "{params.label_b}" > {log} 2>&1
        """

# ------------------------------------------------------------------------------
# TRIPLEX: Three-Target Screening (e.g. Influenza A + Influenza B + Measles)
# ------------------------------------------------------------------------------
rule m11_triplex_compatibility:
    input:
        target_a = lambda wc: _find_m09_candidates(wc.build_a),
        target_b = lambda wc: _find_m09_candidates(wc.build_b),
        target_c = lambda wc: _find_m09_candidates(wc.build_c)
    output:
        json   = "results/multiplex-triplex-{build_a}_{build_b}_{build_c}/triplex_candidates.json",
        report = "results/multiplex-triplex-{build_a}_{build_b}_{build_c}/triplex_report.md"
    params:
        outdir         = "results/multiplex-triplex-{build_a}_{build_b}_{build_c}",
        num_candidates = 10,
        max_variants   = config.get("primerlab", {}).get("design", {}).get("degenerate", {}).get("max_iupac_variants", 16),
        dye_a          = lambda wc: _dye_for(wc.build_a),
        dye_b          = lambda wc: _dye_for(wc.build_b),
        dye_c          = lambda wc: _dye_for(wc.build_c),
        label_a        = lambda wc: _label_for(wc.build_a),
        label_b        = lambda wc: _label_for(wc.build_b),
        label_c        = lambda wc: _label_for(wc.build_c)
    log:
        "logs/m11_triplex_{build_a}_{build_b}_{build_c}.log"
    shell:
        """
        python scripts/m11_triplex_checker.py \
            --file-a {input.target_a} \
            --file-b {input.target_b} \
            --file-c {input.target_c} \
            --outdir {params.outdir} \
            --num-candidates {params.num_candidates} \
            --max-variants {params.max_variants} \
            --dye-a {params.dye_a} --dye-b {params.dye_b} --dye-c {params.dye_c} \
            --label-a "{params.label_a}" --label-b "{params.label_b}" --label-c "{params.label_c}" > {log} 2>&1
        """

# ------------------------------------------------------------------------------
# EXPORT: Consolidate all pipeline outputs for supervisor review
# ------------------------------------------------------------------------------
def _existing_multiplex_outputs(wildcards):
    """Discover any multiplex/triplex candidate JSONs already produced. The export
    rule has a fully-literal output, so it carries NO wildcards — its trigger cannot
    reference build wildcards. Returning whatever exists keeps the rule valid and
    orders it after multiplex screening when those outputs are present."""
    return sorted(_glob.glob("results/multiplex*/*_candidates.json"))


rule m11_export_supplementary:
    input:
        # Discovery trigger (may be empty): the script itself scans all of results/.
        trigger = _existing_multiplex_outputs
    output:
        xlsx_candidates = "supplementary/S1_Primer_Candidates.xlsx",
        xlsx_specificity = "supplementary/S2_BLAST_Specificity.xlsx",
        xlsx_multiplex  = "supplementary/S3_Multiplex_Compatibility.xlsx",
        xlsx_qc         = "supplementary/S4_QC_Summary.xlsx",
        xlsx_conserved  = "supplementary/S5_Conservation_Scoring.xlsx",
        qc_dossier      = "supplementary/qc_dossier/index.html",
        readme          = "supplementary/README.md"
    params:
        results_dir     = "results",
        outdir          = "supplementary",
        only_degenerate = "--only-degenerate" if config.get("primerlab", {}).get("design", {}).get("degenerate", {}).get("enabled", True) else "",
        # Restrict exports to the builds active in this run so excluded viruses do
        # not leak into the supplementary workbooks.
        active_flag     = ("--active-targets " + " ".join(config["active_targets"])) if config.get("active_targets") else "",
        mirror_dir      = config.get("supplementary_mirror_dir", "")
    log:
        "logs/m11_export_supplementary.log"
    shell:
        """
        python scripts/export_supplementary.py \
            --results-dir {params.results_dir} \
            --outdir {params.outdir} \
            {params.only_degenerate} \
            {params.active_flag} > {log} 2>&1

        # Per-assay QC decision dossiers (HTML + markdown), from the same results/.
        # index.html is always written (even with zero assays), so it anchors the
        # rule; a run with no multiplex results yet is tolerated.
        python scripts/qc_dossier.py \
            --results-dir {params.results_dir} \
            --out {params.outdir}/qc_dossier >> {log} 2>&1 || true

        # Package manifest (this output was previously declared but never written).
        printf '%s\n' \
            "Supplementary package generated by rule m11_export_supplementary." \
            "Executive_Summary_Diagnostic_Assays.xlsx: sheets 00-06 incl. the human" \
            "off-target audit (05) and the cited QC decision audit (06)." \
            "qc_dossier/index.html: per-assay QC decision dossiers (HTML + markdown)." \
            > {params.outdir}/README.md

        # Optional mirror to an external folder (config: supplementary_mirror_dir).
        if [ -n "{params.mirror_dir}" ] && [ -d "$(dirname '{params.mirror_dir}')" ]; then
            echo "[Snakemake] Mirroring supplementary folder to {params.mirror_dir}..." >> {log}
            mkdir -p "{params.mirror_dir}"
            cp -r {params.outdir}/* "{params.mirror_dir}/" >> {log} 2>&1
        fi
        """
