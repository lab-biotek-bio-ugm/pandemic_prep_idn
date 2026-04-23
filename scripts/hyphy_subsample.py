#!/usr/bin/env python3
"""
M05-pre: Diversity-Aware Subsampler for HyPhy (Track B)

Reduces the full MSA to a representative subset before gene extraction
and HyPhy analysis. Runs ONCE on aligned.fasta — all genes share the
same consistent subset.

Adaptive grouping strategy (tries levels in order until quota >= 2 per group):
  Level 1: clade x country x year
  Level 2: clade x country
  Level 3: clade only
  Level 4: random (final fallback)

Bypass conditions:
  - total sequences <= min_sequences: use all (no subsampling needed)
  - total sequences <= max_sequences: use all (already within budget)
  - Nextclade disabled (no clade column): falls back to country x year
"""

import argparse
import csv
import json
import random
import sys
from pathlib import Path

from Bio import SeqIO


def load_metadata(metadata_path, id_col="strain"):
    """Load metadata TSV. Returns dict {strain_id: {clade, country, year}}."""
    records = {}
    with open(metadata_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            strain = row.get(id_col, "").strip()
            if not strain:
                continue
            date = row.get("date", "") or ""
            year = date[:4] if len(date) >= 4 else "unknown"
            records[strain] = {
                "subtype": (row.get("subtype") or "unknown").strip() or "unknown",
                "clade": (row.get("clade") or "unknown").strip() or "unknown",
                "country": (row.get("country") or "unknown").strip() or "unknown",
                "year": year,
            }
    return records


def build_groups(alignment_ids, metadata, dims):
    """Group alignment IDs by the given metadata dimensions.
    Returns dict {group_key: [list of strain IDs]}.
    Strains with no metadata entry are placed into '_ungrouped'.
    """
    groups = {}
    ungrouped = []

    for sid in alignment_ids:
        meta = metadata.get(sid)
        if meta is None:
            ungrouped.append(sid)
            continue
        key = "|".join(meta.get(d, "unknown") for d in dims)
        groups.setdefault(key, []).append(sid)

    if ungrouped:
        groups.setdefault("_ungrouped", []).extend(ungrouped)

    return groups


def adaptive_subsample(alignment_ids, metadata, max_sequences, seed, group_config):
    """Try grouping levels from most to least specific until quota >= 2.
    Returns (selected_ids, level_used, group_report).
    """
    use_subtype = group_config.get("subtype", True)
    use_clade = group_config.get("clade", False)
    use_country = group_config.get("country", True)
    use_year = group_config.get("year", True)

    # Build candidate levels, most specific first
    levels = []
    dims = []
    if use_subtype:
        dims.append("subtype")
    if use_clade:
        dims.append("clade")
    if use_country:
        dims.append("country")
    if use_year:
        dims.append("year")

    # Generate progressively simpler combinations
    for drop_count in range(len(dims) + 1):
        candidate = dims[:len(dims) - drop_count]
        if candidate:
            label = " x ".join(candidate)
            levels.append((label, candidate))

    # Always include pure random as final fallback
    levels.append(("random", []))

    rng = random.Random(seed)

    for level_name, level_dims in levels:
        if not level_dims:
            # Pure random fallback
            selected = rng.sample(alignment_ids, min(max_sequences, len(alignment_ids)))
            return selected, "random", {"note": "pure random fallback"}

        groups = build_groups(alignment_ids, metadata, level_dims)
        n_groups = len(groups)
        quota = max_sequences // n_groups

        if quota < 2 and n_groups > 1:
            print(f"[SUBSAMPLE]   Level '{level_name}': {n_groups} groups, "
                  f"quota={quota} — too sparse, trying simpler grouping",
                  file=sys.stderr)
            continue

        # Sample from each group
        selected = []
        report = {}
        for group_key, members in sorted(groups.items()):
            take = min(len(members), quota)
            chosen = rng.sample(members, take) if len(members) > take else list(members)
            selected.extend(chosen)
            report[group_key] = {"available": len(members), "selected": take}

        # Top-up to max_sequences with remaining sequences
        if len(selected) < max_sequences:
            selected_set = set(selected)
            remainder = [s for s in alignment_ids if s not in selected_set]
            extra = min(max_sequences - len(selected), len(remainder))
            if extra > 0:
                selected.extend(rng.sample(remainder, extra))

        print(f"[SUBSAMPLE] Level used: '{level_name}' "
              f"({n_groups} groups, ~{quota} seqs/group)")
        return selected, level_name, report

    # Should never reach here
    selected = rng.sample(alignment_ids, min(max_sequences, len(alignment_ids)))
    return selected, "random", {}


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--alignment", required=True,
                        help="Full MSA FASTA from M03 (aligned.fasta)")
    parser.add_argument("--metadata", required=True,
                        help="Enriched metadata TSV (metadata_enriched.tsv)")
    parser.add_argument("--output", required=True,
                        help="Output subsampled FASTA")
    parser.add_argument("--report", required=True,
                        help="Output JSON report")
    parser.add_argument("--max-sequences", type=int, default=300,
                        help="Maximum sequences to retain (default: 300)")
    parser.add_argument("--min-sequences", type=int, default=50,
                        help="If total <= this, skip subsampling (default: 50)")
    parser.add_argument("--strategy", default="balanced",
                        choices=["balanced", "random"],
                        help="Sampling strategy (default: balanced)")
    parser.add_argument("--group-by-subtype", default=True, action="store_true")
    parser.add_argument("--no-group-by-subtype", dest="group_by_subtype",
                        action="store_false")
    parser.add_argument("--group-by-clade", default=True, action="store_true")
    parser.add_argument("--no-group-by-clade", dest="group_by_clade",
                        action="store_false")
    parser.add_argument("--group-by-country", default=True, action="store_true")
    parser.add_argument("--no-group-by-country", dest="group_by_country",
                        action="store_false")
    parser.add_argument("--group-by-year", default=True, action="store_true")
    parser.add_argument("--no-group-by-year", dest="group_by_year",
                        action="store_false")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility (default: 42)")
    args = parser.parse_args()

    print(f"[SUBSAMPLE] Loading alignment: {args.alignment}")
    records = list(SeqIO.parse(args.alignment, "fasta"))
    alignment_ids = [r.id for r in records]
    n_total = len(alignment_ids)
    print(f"[SUBSAMPLE] Total sequences in MSA: {n_total}")
    print(f"[SUBSAMPLE] Budget: max={args.max_sequences}, min_bypass={args.min_sequences}")

    # Bypass conditions
    if n_total <= args.min_sequences:
        print(f"[SUBSAMPLE] BYPASS: total ({n_total}) <= min_sequences ({args.min_sequences}). "
              f"Using all sequences.")
        selected_ids = alignment_ids
        level = "bypass (total <= min_sequences)"
        group_report = {}

    elif n_total <= args.max_sequences:
        print(f"[SUBSAMPLE] BYPASS: total ({n_total}) <= max_sequences ({args.max_sequences}). "
              f"Using all sequences.")
        selected_ids = alignment_ids
        level = "bypass (total <= max_sequences)"
        group_report = {}

    elif args.strategy == "random":
        rng = random.Random(args.seed)
        selected_ids = rng.sample(alignment_ids, args.max_sequences)
        level = "random (strategy override)"
        group_report = {}
        print(f"[SUBSAMPLE] Random strategy: selected {len(selected_ids)} sequences")

    else:
        # Balanced / adaptive grouping
        print(f"[SUBSAMPLE] Loading metadata: {args.metadata}")
        metadata = load_metadata(args.metadata)
        matched = sum(1 for sid in alignment_ids if sid in metadata)
        print(f"[SUBSAMPLE] Metadata matches: {matched}/{n_total} sequences")

        if matched == 0:
            print("[SUBSAMPLE] WARNING: No metadata matches found. "
                  "Falling back to random sampling.", file=sys.stderr)
            rng = random.Random(args.seed)
            selected_ids = rng.sample(alignment_ids, args.max_sequences)
            level = "random (no metadata matches)"
            group_report = {}
        else:
            group_config = {
                "subtype": args.group_by_subtype,
                "clade": args.group_by_clade,
                "country": args.group_by_country,
                "year": args.group_by_year,
            }
            selected_ids, level, group_report = adaptive_subsample(
                alignment_ids, metadata, args.max_sequences, args.seed, group_config
            )

    # Write subsampled FASTA
    id_set = set(selected_ids)
    out_records = [r for r in records if r.id in id_set]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(out_records, args.output, "fasta")

    n_retained = len(out_records)
    n_skipped = n_total - n_retained

    print(f"[SUBSAMPLE] Retained: {n_retained} / {n_total} "
          f"(skipped: {n_skipped})")
    print(f"[SUBSAMPLE] Written to: {args.output}")

    # Write JSON report
    report = {
        "original_sequences": n_total,
        "retained_sequences": n_retained,
        "skipped_sequences": n_skipped,
        "level_used": level,
        "max_sequences": args.max_sequences,
        "min_sequences": args.min_sequences,
        "seed": args.seed,
        "strategy": args.strategy,
        "group_by": {
            "clade": args.group_by_clade,
            "country": args.group_by_country,
            "year": args.group_by_year,
        },
        "groups": group_report,
    }
    with open(args.report, "w") as f:
        json.dump(report, f, indent=2)
    print(f"[SUBSAMPLE] Report written: {args.report}")


if __name__ == "__main__":
    main()
