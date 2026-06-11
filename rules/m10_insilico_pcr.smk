"""
Module 10: In-silico PCR Coverage Analysis
Validates candidate primers and probes against a large sequence database using seqkit amplicon
and calculates coverage metrics broken down by subtype, collection year, and country.
"""

def get_m10_input_candidates(wildcards):
    import yaml
    # Force Snakemake to evaluate checkpoint m07_bridge
    checkpoint_output = checkpoints.m07_bridge.get(BUILD=wildcards.BUILD)
    manifest_path = checkpoint_output.output.manifest
    
    with open(manifest_path, "r") as f:
        manifest = yaml.safe_load(f)
    
    regions = [r["region_id"] for r in manifest.get("candidate_regions", [])]
    return [
        f"results/{wildcards.BUILD}/m09/{region}/final_candidates.json"
        for region in regions
    ]

rule m10_insilico_pcr:
    """
    Executes in-silico PCR simulation using seqkit amplicon and custom Python scripts
    to validate candidate coverage against a large sequence database.
    """
    input:
        candidates=get_m10_input_candidates,
        sequences=config["local_input"]["sequences"],
        metadata=config["local_input"]["metadata"],
        gff=config["references"]["annotation"]
    output:
        summary_tsv="results/{BUILD}/m10/coverage_summary.tsv",
        summary_json="results/{BUILD}/m10/coverage_summary.json",
        report_md="results/{BUILD}/m10/coverage_report.md",
        heatmap_html="results/{BUILD}/m10/coverage_heatmap.html",
        gff3_track="results/{BUILD}/m10/m10_candidates.gff3"
    params:
        max_mismatch=config["insilico_pcr"].get("max_mismatch", 3),
        num_candidates=config["insilico_pcr"].get("num_candidates", 10),
        outdir="results/{BUILD}/m10"
    conda:
        "../envs/m10_insilico.yaml"
    log:
        "logs/{BUILD}/m10/insilico_pcr.log"
    shell:
        """
        python scripts/m10_insilico_pcr.py \
            --candidates {input.candidates} \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --gff {input.gff} \
            --outdir {params.outdir} \
            --max-mismatch {params.max_mismatch} \
            --num-candidates {params.num_candidates} > {log} 2>&1
        """
