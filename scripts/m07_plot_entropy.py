"""
M07 Entropy Profile Visualisation

Generates a publication-quality PNG/SVG plot of per-position Shannon entropy
across the reference genome, with primer candidate zones highlighted.

Inputs:
  --entropy-tsv : Path to entropy_profile.tsv from M07 Bridge Builder
  --bridge-json : Path to bridge_report.json (for region annotations)
  --candidates  : (Optional) Path to final_candidates.json from M09 for primer overlay
  --outdir      : Output directory for generated plots
  --format      : Output format (png, svg, pdf). Default: png

Outputs:
  entropy_profile_plot.{format}  — Full-genome entropy profile
"""

import argparse
import csv
import json
import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")  # Non-interactive backend for server/CI environments
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.collections import BrokenBarHCollection
except ImportError:
    print("❌ matplotlib is required. Install via: pip install matplotlib")
    sys.exit(1)


def load_entropy_profile(tsv_path: str) -> list:
    """Load entropy_profile.tsv into a list of dicts."""
    profile = []
    with open(tsv_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            profile.append({
                "position": int(row["position"]),
                "entropy": float(row["entropy"]),
                "consensus_base": row.get("consensus_base", ""),
                "iupac_base": row.get("iupac_base", ""),
            })
    return profile


def load_bridge_report(json_path: str) -> dict:
    """Load bridge_report.json."""
    with open(json_path, "r") as f:
        return json.load(f)


def load_candidates(json_path: str) -> list:
    """Load final_candidates.json from M09."""
    with open(json_path, "r") as f:
        return json.load(f)


def plot_entropy(profile, bridge_data, candidates, outdir, fmt="png"):
    """Generate entropy profile plot with annotations."""
    positions = [p["position"] for p in profile]
    entropies = [p["entropy"] for p in profile]

    fig, ax = plt.subplots(figsize=(16, 5), dpi=150)

    # Main entropy line
    ax.fill_between(positions, entropies, alpha=0.3, color="#2196F3", linewidth=0)
    ax.plot(positions, entropies, color="#1565C0", linewidth=0.5, alpha=0.8)

    # Highlight conserved regions from bridge report
    regions = bridge_data.get("regions", [])
    for region in regions:
        start = region["start"]
        end = region["end"]
        grade = region.get("grade", "D")
        # Color by grade
        color_map = {
            "S": "#4CAF50",
            "A+": "#66BB6A",
            "A": "#81C784",
            "A-": "#A5D6A7",
            "B+": "#FFB74D",
            "B": "#FFA726",
            "C": "#EF5350",
            "D": "#BDBDBD",
        }
        color = color_map.get(grade, "#BDBDBD")
        ax.axvspan(start, end, alpha=0.15, color=color, zorder=0)

    # Overlay primer candidates (from M09 final_candidates.json)
    if candidates:
        for cand in candidates:
            rank = cand.get("rank", "?")
            amp_start = cand.get("amplicon_start")
            amp_end = cand.get("amplicon_end")
            status = cand.get("status", "PASS")

            if amp_start is not None and amp_end is not None:
                edge_color = "#4CAF50" if status == "PASS" else "#FF9800"
                ax.axvspan(amp_start, amp_end, alpha=0.08, color=edge_color, zorder=1)
                ax.axvline(amp_start, color=edge_color, linewidth=0.8, linestyle="--", alpha=0.6)
                ax.axvline(amp_end, color=edge_color, linewidth=0.8, linestyle="--", alpha=0.6)
                mid = (amp_start + amp_end) / 2
                ax.text(mid, max(entropies) * 0.95, f"C{rank}",
                        ha="center", va="top", fontsize=7, fontweight="bold",
                        color=edge_color, alpha=0.8,
                        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor=edge_color, alpha=0.7))

    # Highlight GARD breakpoints
    gard_bps = bridge_data.get("gard_breakpoints", [])
    for bp in gard_bps:
        gp = bp.get("genome_pos", -1)
        if gp >= 0:
            ax.axvline(gp, color="#F44336", linewidth=1.2, linestyle=":", alpha=0.7, zorder=2)

    # Gene annotations from bridge report
    gene_regions = {}
    for region in regions:
        gene = region.get("gene", "")
        if gene and gene not in gene_regions:
            gene_regions[gene] = {"start": region["start"], "end": region["end"]}
        elif gene:
            gene_regions[gene]["start"] = min(gene_regions[gene]["start"], region["start"])
            gene_regions[gene]["end"] = max(gene_regions[gene]["end"], region["end"])

    # Entropy threshold line
    ax.axhline(y=0.02, color="#FF5722", linewidth=0.8, linestyle="-.", alpha=0.5,
               label="Conservation threshold (H=0.02)")

    # Formatting
    ax.set_xlabel("Genome Position (bp)", fontsize=11)
    ax.set_ylabel("Shannon Entropy (bits)", fontsize=11)
    ax.set_title("Per-Position Shannon Entropy Profile with Primer Candidate Zones", fontsize=13, fontweight="bold")
    ax.set_xlim(min(positions), max(positions))
    ax.set_ylim(0, min(max(entropies) * 1.1, 2.0))

    # Legend
    legend_patches = [
        mpatches.Patch(color="#1565C0", alpha=0.5, label="Entropy"),
        mpatches.Patch(color="#4CAF50", alpha=0.3, label="Grade S/A Regions"),
        mpatches.Patch(color="#FFB74D", alpha=0.3, label="Grade B Regions"),
        mpatches.Patch(color="#4CAF50", alpha=0.15, label="PASS Primer Amplicon"),
        mpatches.Patch(color="#FF9800", alpha=0.15, label="BORDERLINE Primer Amplicon"),
        plt.Line2D([0], [0], color="#F44336", linewidth=1.2, linestyle=":", label="GARD Breakpoint"),
        plt.Line2D([0], [0], color="#FF5722", linewidth=0.8, linestyle="-.", label="Threshold (H=0.02)"),
    ]
    ax.legend(handles=legend_patches, loc="upper right", fontsize=7, framealpha=0.9)

    ax.grid(axis="y", alpha=0.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    outpath = Path(outdir) / f"entropy_profile_plot.{fmt}"
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(outpath), format=fmt, bbox_inches="tight")
    plt.close()
    print(f"📊 Entropy profile plot saved: {outpath}")
    return str(outpath)


def main():
    parser = argparse.ArgumentParser(description="M07 Entropy Profile Visualisation")
    parser.add_argument("--entropy-tsv", required=True, help="Path to entropy_profile.tsv")
    parser.add_argument("--bridge-json", required=True, help="Path to bridge_report.json")
    parser.add_argument("--candidates", default=None, help="Path to final_candidates.json (optional)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--format", default="png", choices=["png", "svg", "pdf"], help="Output format")

    args = parser.parse_args()

    print(f"📈 Loading entropy profile from: {args.entropy_tsv}")
    profile = load_entropy_profile(args.entropy_tsv)
    print(f"   → {len(profile)} positions loaded.")

    print(f"📋 Loading bridge report from: {args.bridge_json}")
    bridge = load_bridge_report(args.bridge_json)

    candidates = None
    if args.candidates and Path(args.candidates).exists():
        print(f"🎯 Loading primer candidates from: {args.candidates}")
        candidates = load_candidates(args.candidates)
        print(f"   → {len(candidates)} candidates loaded.")

    plot_entropy(profile, bridge, candidates, args.outdir, fmt=args.format)
    print("✅ Done.")


if __name__ == "__main__":
    main()
