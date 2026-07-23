"""
M07 Entropy Profile Visualisation  (v2.0 — Audit Fixes ISU-1/2/3/6)

Generates a publication-quality PNG/SVG plot of per-position Shannon entropy
across the reference genome, with:
  - Normalised entropy (H_norm = H / 2) on Y-axis  [ISU-2]
  - Gene annotation track drawn from bridge_report  [ISU-1]
  - 3' terminus highlight for each primer candidate [ISU-3]
  - MEME and FUBAR selection site markers           [ISU-6]
  - GARD recombination breakpoints (unchanged)
  - Conserved region shading from bridge_report

Inputs:
  --entropy-tsv : Path to entropy_profile.tsv from M07 Bridge Builder
  --bridge-json : Path to bridge_report.json (for region annotations)
  --candidates  : (Optional) Path to final_candidates.json from M09 for primer overlay
  --outdir      : Output directory for generated plots
  --format      : Output format (png, svg, pdf). Default: png
  --prime3-len  : Number of bp at the 3' end of each primer to highlight (default: 8)

Outputs:
  entropy_profile_plot.{format}  — Full-genome entropy profile with all annotations
"""

import argparse
import csv
import json
import math
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

# ---------------------------------------------------------------------------
# H_MAX: Maximum Shannon entropy for the 5-symbol alphabet {A,T,G,C,-} (log2(5)).
# Must match MAX_ENTROPY_BITS in m07_bridge_builder.py: entropy now includes gaps,
# so raw values reach ~2.322 and a hardcoded 2.0 would push deletion hotspots
# off the top of the plot — exactly the columns this view exists to expose.
# Used to normalise raw entropy to [0, 1]:  H_norm = H / H_MAX
# ---------------------------------------------------------------------------
H_MAX = math.log2(5)


def load_entropy_profile(tsv_path: str) -> list:
    """Load entropy_profile.tsv into a list of dicts.

    Adds a normalised 'entropy_norm' key: entropy_norm = entropy / H_MAX
    so downstream code always works with the dimensionless [0, 1] scale.
    """
    profile = []
    with open(tsv_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            raw_h = float(row["entropy"])
            profile.append({
                "position": int(row["position"]),
                "entropy": raw_h,
                "entropy_norm": raw_h / H_MAX,          # [ISU-2] normalised
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


def _gene_spans_from_bridge(bridge_data: dict) -> dict:
    """Extract gene-level spans from bridge_report regions.

    Returns {gene_name: {"start": int, "end": int}}.
    """
    gene_regions: dict = {}
    for region in bridge_data.get("regions", []):
        gene = region.get("gene", "")
        if not gene:
            continue
        if gene not in gene_regions:
            gene_regions[gene] = {"start": region["start"], "end": region["end"]}
        else:
            gene_regions[gene]["start"] = min(gene_regions[gene]["start"], region["start"])
            gene_regions[gene]["end"]   = max(gene_regions[gene]["end"],   region["end"])
    return gene_regions


def plot_entropy(
    profile: list,
    bridge_data: dict,
    candidates: list,
    outdir: str,
    fmt: str = "png",
    prime3_len: int = 8,
) -> str:
    """Generate entropy profile plot with all annotations.

    Parameters
    ----------
    profile     : list of per-position dicts (with entropy_norm key)
    bridge_data : parsed bridge_report.json
    candidates  : list of candidate dicts from M09 final_candidates.json, or None
    outdir      : output directory
    fmt         : image format ('png', 'svg', 'pdf')
    prime3_len  : number of bp at 3' end to highlight with a critical marker [ISU-3]
    """
    positions  = [p["position"] for p in profile]
    # ISU-2: use normalised entropy, not raw H
    entropies  = [p["entropy_norm"] for p in profile]
    max_h_norm = max(entropies) if entropies else 1.0

    # Build figure with two sub-axes:
    #  ax_gene  — narrow gene annotation track  [ISU-1]
    #  ax_main  — entropy + selections + primers
    fig, (ax_gene, ax_main) = plt.subplots(
        2, 1,
        figsize=(16, 6.5),
        dpi=150,
        gridspec_kw={"height_ratios": [1, 6], "hspace": 0.05},
        sharex=True,
    )

    # ─────────────────────────────────────────────────────────────────────────
    # TRACK 1 — Gene Annotation  [ISU-1]
    # ─────────────────────────────────────────────────────────────────────────
    gene_colors = [
        "#3B82F6", "#10B981", "#F59E0B", "#EF4444",
        "#8B5CF6", "#EC4899", "#14B8A6", "#F97316",
    ]
    gene_regions = _gene_spans_from_bridge(bridge_data)
    ax_gene.set_ylim(0, 1)
    ax_gene.set_ylabel("Genes", fontsize=8, labelpad=4)
    ax_gene.tick_params(axis="both", which="both", length=0, labelleft=False)
    ax_gene.spines["top"].set_visible(False)
    ax_gene.spines["right"].set_visible(False)
    ax_gene.spines["bottom"].set_visible(False)
    ax_gene.set_facecolor("#F8FAFC")

    if gene_regions:
        for idx, (gene, span) in enumerate(sorted(gene_regions.items(), key=lambda x: x[1]["start"])):
            color = gene_colors[idx % len(gene_colors)]
            bar_width = span["end"] - span["start"]
            ax_gene.barh(
                0.5, bar_width, left=span["start"], height=0.55,
                color=color, alpha=0.80, edgecolor="white", linewidth=1.2,
            )
            mid = span["start"] + bar_width / 2
            ax_gene.text(
                mid, 0.5, gene,
                ha="center", va="center", fontsize=7.5, fontweight="bold",
                color="white",
                clip_on=True,
            )
    else:
        # Fallback: label as "Genomic Segment" if no gene info available
        xmin, xmax = min(positions), max(positions)
        ax_gene.barh(0.5, xmax - xmin, left=xmin, height=0.55,
                     color="#6B7280", alpha=0.6, edgecolor="white", linewidth=1.2)
        ax_gene.text((xmin + xmax) / 2, 0.5, "Genomic Segment",
                     ha="center", va="center", fontsize=7.5, fontweight="bold",
                     color="white")

    # ─────────────────────────────────────────────────────────────────────────
    # TRACK 2 — Conserved region shading (from bridge_report)
    # ─────────────────────────────────────────────────────────────────────────
    regions = bridge_data.get("regions", [])
    color_map = {
        "S":  "#4CAF50", "A+": "#66BB6A", "A":  "#81C784", "A-": "#A5D6A7",
        "B+": "#FFB74D", "B":  "#FFA726", "C":  "#EF5350", "D":  "#BDBDBD",
    }
    for region in regions:
        start = region["start"]
        end   = region["end"]
        grade = region.get("grade", "D")
        color = color_map.get(grade, "#BDBDBD")
        ax_main.axvspan(start, end, alpha=0.10, color=color, zorder=0)

    # ─────────────────────────────────────────────────────────────────────────
    # TRACK 2 — Normalised entropy line  [ISU-2]
    # ─────────────────────────────────────────────────────────────────────────
    ax_main.fill_between(positions, entropies, alpha=0.25, color="#2196F3", linewidth=0)
    ax_main.plot(positions, entropies, color="#1565C0", linewidth=0.6, alpha=0.85)

    # Conservation threshold (normalised): 0.02 raw → 0.01 normalised  [ISU-2]
    THRESHOLD_NORM = 0.01
    ax_main.axhline(
        y=THRESHOLD_NORM, color="#FF5722", linewidth=0.9, linestyle="-.",
        alpha=0.6, label=f"Conservation threshold (H_norm = {THRESHOLD_NORM})",
    )

    # ─────────────────────────────────────────────────────────────────────────
    # TRACK 2 — GARD recombination breakpoints
    # ─────────────────────────────────────────────────────────────────────────
    gard_bps = bridge_data.get("gard_breakpoints", [])
    for bp in gard_bps:
        gp = bp.get("genome_pos", -1)
        if gp >= 0:
            ax_main.axvline(
                gp, color="#F44336", linewidth=1.2, linestyle=":", alpha=0.75, zorder=2,
                label="_GARD_bp",
            )

    # ─────────────────────────────────────────────────────────────────────────
    # ISU-6 — MEME episodic selection sites (blue dashed)
    # ─────────────────────────────────────────────────────────────────────────
    meme_sites = bridge_data.get("meme_sites", [])
    _meme_labeled = False
    for site in meme_sites:
        pos = site.get("genome_pos", site.get("pos", -1))
        if pos < 0:
            continue
        lbl = "_MEME_site" if _meme_labeled else "MEME episodic selection site"
        ax_main.axvline(
            pos, color="#3F51B5", linewidth=0.9, linestyle="--", alpha=0.55,
            zorder=2, label=lbl,
        )
        _meme_labeled = True

    # ─────────────────────────────────────────────────────────────────────────
    # ISU-6 — FUBAR constitutive positive-selection sites (orange dotted)
    # ─────────────────────────────────────────────────────────────────────────
    fubar_pos_sites = bridge_data.get("fubar_positive_sites", [])
    _fubar_labeled = False
    for site in fubar_pos_sites:
        pos = site.get("genome_pos", site.get("pos", -1))
        if pos < 0:
            continue
        lbl = "_FUBAR_site" if _fubar_labeled else "FUBAR positive selection site"
        ax_main.axvline(
            pos, color="#FF9800", linewidth=0.8, linestyle=":", alpha=0.55,
            zorder=2, label=lbl,
        )
        _fubar_labeled = True

    # ─────────────────────────────────────────────────────────────────────────
    # TRACK 2 — Primer candidate amplicons + 3' end markers  [ISU-3]
    # ─────────────────────────────────────────────────────────────────────────
    if candidates:
        for cand in candidates:
            rank      = cand.get("rank", "?")
            amp_start = cand.get("amplicon_start")
            amp_end   = cand.get("amplicon_end")
            fwd_end   = cand.get("fwd_end")    # 3' end of forward primer
            rev_start = cand.get("rev_start")  # 3' end of reverse primer (reversed)
            status    = cand.get("status", "PASS")

            if amp_start is None or amp_end is None:
                continue

            edge_color = "#4CAF50" if status == "PASS" else "#FF9800"

            # Amplicon shading
            ax_main.axvspan(amp_start, amp_end, alpha=0.07, color=edge_color, zorder=1)
            ax_main.axvline(amp_start, color=edge_color, linewidth=0.9, linestyle="--", alpha=0.6)
            ax_main.axvline(amp_end,   color=edge_color, linewidth=0.9, linestyle="--", alpha=0.6)

            # Candidate rank label
            mid = (amp_start + amp_end) / 2
            ax_main.text(
                mid, max_h_norm * 0.93, f"C{rank}",
                ha="center", va="top", fontsize=7, fontweight="bold",
                color=edge_color, alpha=0.85,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          edgecolor=edge_color, alpha=0.75),
            )

            # ISU-3: 3' end critical zone highlights
            # Forward primer 3' end = last `prime3_len` bp at fwd_end
            if fwd_end is not None:
                p3_fwd_start = max(amp_start, fwd_end - prime3_len)
                ax_main.axvspan(
                    p3_fwd_start, fwd_end,
                    alpha=0.25, color="#EF4444", zorder=2,
                    label="3ʹ Primer terminus (critical zone)" if rank == 1 else "_3prime",
                )
            # Reverse primer 3' end = first `prime3_len` bp from rev_start (reversed orientation)
            if rev_start is not None:
                p3_rev_end = min(amp_end, rev_start + prime3_len)
                ax_main.axvspan(
                    rev_start, p3_rev_end,
                    alpha=0.25, color="#EF4444", zorder=2,
                    label="_3prime",
                )

    # ─────────────────────────────────────────────────────────────────────────
    # Formatting — main axis
    # ─────────────────────────────────────────────────────────────────────────
    ax_main.set_xlabel("Genome Position (bp)", fontsize=11)
    ax_main.set_ylabel("Normalised Shannon Entropy  (H / H_max)", fontsize=10)  # [ISU-2]
    ax_main.set_xlim(min(positions), max(positions))
    ax_main.set_ylim(0, min(max_h_norm * 1.12, 1.05))
    ax_main.grid(axis="y", alpha=0.18)
    ax_main.spines["top"].set_visible(False)
    ax_main.spines["right"].set_visible(False)

    # ─────────────────────────────────────────────────────────────────────────
    # Legend  [ISU-1: Gene track now separate; ISU-2: threshold label updated]
    # ─────────────────────────────────────────────────────────────────────────
    legend_patches = [
        mpatches.Patch(color="#1565C0", alpha=0.55, label="Normalised Entropy (H/H_max)"),
        mpatches.Patch(color="#4CAF50", alpha=0.25, label="Grade S/A Conserved Regions"),
        mpatches.Patch(color="#FFB74D", alpha=0.25, label="Grade B Regions"),
        mpatches.Patch(color="#4CAF50", alpha=0.12, label="PASS Primer Amplicon"),
        mpatches.Patch(color="#FF9800", alpha=0.12, label="BORDERLINE Primer Amplicon"),
        mpatches.Patch(color="#EF4444", alpha=0.30, label="3ʹ Primer Terminus (Critical)"),
        plt.Line2D([0], [0], color="#F44336", linewidth=1.2, linestyle=":", label="GARD Breakpoint"),
        plt.Line2D([0], [0], color="#3F51B5", linewidth=0.9, linestyle="--", label="MEME Episodic Selection"),
        plt.Line2D([0], [0], color="#FF9800", linewidth=0.8, linestyle=":", label="FUBAR Positive Selection"),
        plt.Line2D([0], [0], color="#FF5722", linewidth=0.9, linestyle="-.",
                   label=f"Conservation threshold (H_norm ≤ {THRESHOLD_NORM})"),
    ]
    ax_main.legend(handles=legend_patches, loc="upper right", fontsize=7, framealpha=0.92)

    fig.suptitle(
        "Per-Position Normalised Shannon Entropy Profile with Primer Candidate Zones",
        fontsize=12, fontweight="bold", y=1.01,
    )

    plt.tight_layout()

    outpath = Path(outdir) / f"entropy_profile_plot.{fmt}"
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(outpath), format=fmt, bbox_inches="tight")
    plt.close()
    print(f"📊 Entropy profile plot saved: {outpath}")
    return str(outpath)


def main():
    parser = argparse.ArgumentParser(description="M07 Entropy Profile Visualisation (v2.0)")
    parser.add_argument("--entropy-tsv", required=True, help="Path to entropy_profile.tsv")
    parser.add_argument("--bridge-json", required=True, help="Path to bridge_report.json")
    parser.add_argument("--candidates",  default=None,  help="Path to final_candidates.json (optional)")
    parser.add_argument("--outdir",      required=True, help="Output directory")
    parser.add_argument(
        "--format", default="png", choices=["png", "svg", "pdf"], help="Output format",
    )
    parser.add_argument(
        "--prime3-len", type=int, default=8,
        help="Number of bp at 3' end of each primer to mark as 'critical zone' (default: 8)",
    )

    args = parser.parse_args()

    print(f"📈 Loading entropy profile from: {args.entropy_tsv}")
    profile = load_entropy_profile(args.entropy_tsv)
    print(f"   → {len(profile)} positions loaded (entropy normalised by /{H_MAX}).")

    print(f"📋 Loading bridge report from: {args.bridge_json}")
    bridge = load_bridge_report(args.bridge_json)

    candidates = None
    if args.candidates and Path(args.candidates).exists():
        print(f"🎯 Loading primer candidates from: {args.candidates}")
        candidates = load_candidates(args.candidates)
        print(f"   → {len(candidates)} candidates loaded.")

    plot_entropy(
        profile, bridge, candidates, args.outdir,
        fmt=args.format, prime3_len=args.prime3_len,
    )
    print("✅ Done.")


if __name__ == "__main__":
    main()
