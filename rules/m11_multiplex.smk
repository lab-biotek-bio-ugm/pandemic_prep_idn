# ==============================================================================
# MODULE 11: DUAL-TARGET MULTIPLEX COMPATIBILITY SCREENING
# ==============================================================================
# Evaluates pairwise cross-dimer thermodynamics, spectral fluorophore separation,
# and kinetic Tm alignment between two target diagnostic candidate sets.

rule m11_multiplex_compatibility:
    input:
        target_a = "results/{build_a}/m09/degenerate/segment-{segment}_region_1/final_candidates.json",
        target_b = "results/{build_b}/m09/degenerate/segment-{segment}_region_1/final_candidates.json"
    output:
        json = "results/multiplex-{build_a}_x_{build_b}/multiplex_candidates.json",
        report = "results/multiplex-{build_a}_x_{build_b}/multiplex_report.md",
        html = "results/multiplex-{build_a}_x_{build_b}/multiplex_matrix.html"
    params:
        outdir = "results/multiplex-{build_a}_x_{build_b}",
        num_candidates = 10
    log:
        "logs/m11_multiplex_{build_a}_x_{build_b}.log"
    shell:
        """
        python scripts/m11_multiplex_checker.py \
            --file-a {input.target_a} \
            --file-b {input.target_b} \
            --outdir {params.outdir} \
            --num-candidates {params.num_candidates} > {log} 2>&1
        """
