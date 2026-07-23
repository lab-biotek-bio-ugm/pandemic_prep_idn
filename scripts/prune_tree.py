#!/usr/bin/env python3
"""
Prune a Newick tree to only contain tips that are present in the FASTA alignment.
Called as a Snakemake `script:` directive.
"""

from Bio import Phylo, SeqIO
import sys
import io

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

            # Fail loudly on tree<->alignment desync. Pruning only removes tips that
            # are in the tree but not the alignment; if any alignment sequence has NO
            # matching tree tip, HyPhy would silently run on a mismatched taxon set and
            # produce misleading dN/dS. The tip set MUST be a subset of alignment IDs.
            aln_only = alignment_ids - remaining_tips
            if aln_only:
                preview = ", ".join(sorted(aln_only)[:10])
                logf.write(
                    f"ERROR: {len(aln_only)} alignment sequence(s) have no matching tree tip "
                    f"(tree<->alignment desync). Examples: {preview}\n"
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
