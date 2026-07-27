#!/usr/bin/env python3
"""
Prune a Newick tree to only contain tips that are present in the FASTA alignment.
Called as a Snakemake `script:` directive.
"""

from Bio import Phylo, SeqIO
import sys
import io

# Fraction of the gene ALIGNMENT that may legitimately lack a tree tip. The gene
# alignment comes from the HyPhy subsample (m05, capped at ~2000), while the tree
# (m04) is built on the full m02 cohort, so pruning the tree down to the subsample
# routinely removes most tips -- that is expected and NOT an error. A real
# tree<->alignment ID mismatch (e.g. Newick escaping of '/') instead shows up as
# alignment sequences that match no tree tip; that fraction is what we bound here.
MAX_UNMATCHED_ALIGN_FRACTION = 0.5

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
            # A large pruned fraction is EXPECTED here: the tree is built on the full
            # m02 cohort while this gene alignment is the HyPhy subsample (~2000), so
            # for big builds (all-partial, large sea-only) most tips are legitimately
            # pruned. What must never happen is the alignment sequences failing to map
            # onto the tree at all — the systematic ID-format mismatch (e.g. Newick
            # escaping of strain names containing '/'). We fail on THAT fraction.
            #
            # A few unmatched alignment sequences are also normal, because
            # `augur refine --timetree` drops tips it cannot date / clock outliers, so
            # a subsampled sequence can be absent from the refined tree.
            aln_only = alignment_ids - remaining_tips
            unmatched_frac = len(aln_only) / len(alignment_ids) if alignment_ids else 0.0
            pruned_frac = len(tips_to_prune) / len(all_tips) if all_tips else 0.0
            logf.write(
                f"Pruned {pruned_frac:.1%} of tree tips (expected when the tree is "
                f"larger than the HyPhy subsample); {unmatched_frac:.1%} of alignment "
                f"sequences have no tree tip.\n"
            )
            if aln_only:
                preview = ", ".join(sorted(aln_only)[:10])
                logf.write(
                    f"NOTE: {len(aln_only)} alignment sequence(s) have no tree tip "
                    f"(some expected — augur refine drops undated/clock-outlier tips). "
                    f"Examples: {preview}\n"
                )

            if unmatched_frac > MAX_UNMATCHED_ALIGN_FRACTION:
                logf.write(
                    f"ERROR: {unmatched_frac:.1%} of alignment sequences "
                    f"({len(aln_only)}/{len(alignment_ids)}) had no matching tree tip, "
                    f"exceeding the {MAX_UNMATCHED_ALIGN_FRACTION:.0%} limit. This "
                    f"indicates a systematic tree<->alignment ID mismatch — HyPhy "
                    f"would run on a wrong taxon set. Aborting.\n"
                )
                sys.exit(1)

            # Write pruned tree
            Phylo.write(tree, output_path, "newick")
            logf.write(f"Pruned tree written to {output_path}\n")

        except Exception as e:
            logf.write(f"ERROR: {e}\n")
            raise

# Snakemake script interface. Guarded so the module can be imported standalone
# (e.g. by unit tests); Snakemake injects `snakemake` into globals at run time.
if "snakemake" in globals():
    main(
        tree_path=snakemake.input.tree,
        alignment_path=snakemake.input.alignment,
        output_path=str(snakemake.output.tree),
        log_path=str(snakemake.log[0])
    )
