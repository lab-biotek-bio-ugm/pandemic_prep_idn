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
# ==============================================================================

# ------------------------------------------------------------------------------
# PAIRWISE: Dual-Target Screening (e.g. Influenza A + Influenza B)
# ------------------------------------------------------------------------------
rule m11_multiplex_compatibility:
    input:
        target_a = "results/{build_a}/m09/degenerate/segment-{segment}_region_1/final_candidates.json",
        target_b = "results/{build_b}/m09/degenerate/segment-{segment}_region_1/final_candidates.json"
    output:
        json   = "results/multiplex-{build_a}_x_{build_b}/multiplex_candidates.json",
        report = "results/multiplex-{build_a}_x_{build_b}/multiplex_report.md",
        html   = "results/multiplex-{build_a}_x_{build_b}/multiplex_matrix.html"
    params:
        outdir         = "results/multiplex-{build_a}_x_{build_b}",
        num_candidates = 10,
        max_variants   = config.get("primerlab", {}).get("design", {}).get("degenerate", {}).get("max_iupac_variants", 16)
    log:
        "logs/m11_multiplex_{build_a}_x_{build_b}.log"
    shell:
        """
        python scripts/m11_multiplex_checker.py \
            --file-a {input.target_a} \
            --file-b {input.target_b} \
            --outdir {params.outdir} \
            --num-candidates {params.num_candidates} \
            --max-variants {params.max_variants} > {log} 2>&1
        """

# ------------------------------------------------------------------------------
# TRIPLEX: Three-Target Screening (e.g. Influenza A + Influenza B + Measles)
# ------------------------------------------------------------------------------
rule m11_triplex_compatibility:
    input:
        target_a = "results/{build_a}/m09/degenerate/segment-{seg_a}_region_1/final_candidates.json",
        target_b = "results/{build_b}/m09/degenerate/segment-{seg_b}_region_1/final_candidates.json",
        target_c = "results/{build_c}/m09/degenerate/segment-{seg_c}_region_1/final_candidates.json"
    output:
        json   = "results/multiplex-triplex-{build_a}_{build_b}_{build_c}/triplex_candidates.json",
        report = "results/multiplex-triplex-{build_a}_{build_b}_{build_c}/triplex_report.md"
    params:
        outdir         = "results/multiplex-triplex-{build_a}_{build_b}_{build_c}",
        num_candidates = 10,
        max_variants   = config.get("primerlab", {}).get("design", {}).get("degenerate", {}).get("max_iupac_variants", 16)
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
            --max-variants {params.max_variants} > {log} 2>&1
        """

# ------------------------------------------------------------------------------
# EXPORT: Consolidate all pipeline outputs for supervisor review
# ------------------------------------------------------------------------------
rule m11_export_supplementary:
    input:
        # Requires multiplexing results to be generated first
        # (This is a pseudo-trigger; it discovers all contents of results/)
        trigger = "results/multiplex-triplex-{build_a}_{build_b}_{build_c}/triplex_candidates.json"
    output:
        xlsx_candidates = "supplementary/S1_Primer_Candidates.xlsx",
        xlsx_specificity = "supplementary/S2_BLAST_Specificity.xlsx",
        xlsx_multiplex  = "supplementary/S3_Multiplex_Compatibility.xlsx",
        xlsx_qc         = "supplementary/S4_QC_Summary.xlsx",
        xlsx_conserved  = "supplementary/S5_Conservation_Scoring.xlsx",
        readme          = "supplementary/README.md"
    params:
        results_dir     = "results",
        outdir          = "supplementary",
        only_degenerate = "--only-degenerate" if config.get("primerlab", {}).get("design", {}).get("degenerate", {}).get("enabled", True) else ""
    log:
        "logs/m11_export_supplementary.log"
    shell:
        """
        python scripts/export_supplementary.py \
            --results-dir {params.results_dir} \
            --outdir {params.outdir} \
            {params.only_degenerate} > {log} 2>&1
        """
