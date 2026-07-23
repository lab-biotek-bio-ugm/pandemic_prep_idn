#!/usr/bin/env python3
"""
Prune a Newick tree to only contain tips that are present in the FASTA alignment.
Called as a Snakemake `script:` directive.
"""

from Bio import Phylo, SeqIO
import sys
import io

# Fraction of tree tips that may legitimately lack an alignment sequence.
# Normal `augur refine --timetree` filtering removes only a small share; anything
# beyond this points at a systematic ID-format mismatch, not real filtering.
MAX_PRUNED_FRACTION = 0.5

def main(tree_path, alignment_path, output_path, log_path):
    with open(log_path, "w") as logf:
        try:
            # Get tip names from alignment
            alignment_ids = set()
            for record in SeqIO.parse(alignment_path, "fasta"):
                alignment_ids.add(record.id)
            logf.write(f"Alignment contains {len(alignment_ids)} sequences.\n")

            # Read tree
            tree = Phylo.read(tree_path, "newick")
            all_tips = [tip.name for tip in tree.get_terminals()]
            logf.write(f"Tree contains {len(all_tips)} tips.\n")

            # Determine tips to prune (in tree but NOT in alignment)
            tips_to_prune = [tip for tip in all_tips if tip not in alignment_ids]
            logf.write(f"Pruning {len(tips_to_prune)} tips not found in alignment.\n")

            for tip_name in tips_to_prune:
                tree.prune(tip_name)

            remaining_tips = {tip.name for tip in tree.get_terminals()}
            remaining = len(remaining_tips)
            logf.write(f"Pruned tree has {remaining} tips.\n")

            if remaining < 3:
                logf.write("ERROR: Pruned tree has fewer than 3 tips. HyPhy requires at least 3.\n")
                sys.exit(1)

            # Guard against tree<->alignment desync, WITHOUT breaking legitimate runs.
            #
            # The tree comes from `augur refine --timetree`, which deliberately drops
            # tips it cannot date and clock outliers. So alignment sequences with no
            # matching tree tip are EXPECTED and must not be fatal — they are only
            # reported. HyPhy takes the pruned tree plus the alignment; extra
            # sequences in the alignment are tolerated.
            #
            # The real failure mode is a systematic ID-format mismatch (e.g. Newick
            # escaping of strain names containing '/'), where almost nothing matches.
            # That shows up as a huge pruned fraction, which we DO fail on.
            aln_only = alignment_ids - remaining_tips
            if aln_only:
                preview = ", ".join(sorted(aln_only)[:10])
                logf.write(
                    f"NOTE: {len(aln_only)} alignment sequence(s) have no tree tip "
                    f"(expected — augur refine drops undated/clock-outlier tips). "
                    f"Examples: {preview}\n"
                )

            pruned_frac = len(tips_to_prune) / len(all_tips) if all_tips else 0.0
            if pruned_frac > MAX_PRUNED_FRACTION:
                logf.write(
                    f"ERROR: {pruned_frac:.1%} of tree tips ({len(tips_to_prune)}/{len(all_tips)}) "
                    f"had no matching alignment sequence, exceeding the "
                    f"{MAX_PRUNED_FRACTION:.0%} limit. This indicates a systematic "
                    f"tree<->alignment ID mismatch rather than normal filtering — "
                    f"HyPhy would run on a wrong taxon set. Aborting.\n"
                )
                sys.exit(1)

            # Write pruned tree
            Phylo.write(tree, output_path, "newick")
            logf.write(f"Pruned tree written to {output_path}\n")

        except Exception as e:
            logf.write(f"ERROR: {e}\n")
            raise

# Snakemake script interface
main(
    tree_path=snakemake.input.tree,
    alignment_path=snakemake.input.alignment,
    output_path=str(snakemake.output.tree),
    log_path=str(snakemake.log[0])
)
