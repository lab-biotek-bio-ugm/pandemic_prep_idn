"""Unit tests for scripts/prune_tree.py.

The prune step restricts the full m04 tree (built on the whole m02 cohort) down to
the HyPhy subsample (m05, ~2000 seqs). For large builds (all-partial, big
sea-only) that legitimately removes most tips, so the fatal guard must key off
alignment sequences that fail to map onto the tree -- NOT the pruned fraction.
This regression pins that behaviour.
"""
import importlib.util

import pytest

from conftest import REPO

pytestmark = pytest.mark.unit

Phylo = pytest.importorskip("Bio.Phylo")


def _load_prune_tree():
    spec = importlib.util.spec_from_file_location("prune_tree", REPO / "scripts" / "prune_tree.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # safe: the snakemake call is guarded
    return module


def _write_star_tree(path, names):
    path.write_text("(" + ",".join(f"{n}:0.1" for n in names) + ");\n", encoding="utf-8")


def _write_fasta(path, names):
    path.write_text("".join(f">{n}\nACGT\n" for n in names), encoding="utf-8")


def _count_tips(path):
    return len(Phylo.read(str(path), "newick").get_terminals())


def test_large_tree_pruned_to_subsample_succeeds(tmp_path):
    # Tree has the full cohort (10 tips); alignment is the subsample (3 seqs), all
    # present in the tree. 70% of tips are pruned -- expected, must NOT be fatal.
    tree = tmp_path / "tree.nwk"
    aln = tmp_path / "aln.fasta"
    out = tmp_path / "tree_M1.nwk"
    log = tmp_path / "prune.log"
    _write_star_tree(tree, [f"t{i}" for i in range(1, 11)])
    _write_fasta(aln, ["t1", "t2", "t3"])

    _load_prune_tree().main(str(tree), str(aln), str(out), str(log))

    assert out.exists(), "pruned tree should be written for a legitimate subsample prune"
    assert _count_tips(out) == 3
    assert "Pruned 70.0% of tree tips" in log.read_text(encoding="utf-8")


def test_systematic_id_mismatch_aborts(tmp_path):
    # 3 alignment ids match the tree, but 7 do not -> 70% of the alignment has no
    # tree tip. That is the real desync the guard must still catch.
    tree = tmp_path / "tree.nwk"
    aln = tmp_path / "aln.fasta"
    out = tmp_path / "tree_M1.nwk"
    log = tmp_path / "prune.log"
    _write_star_tree(tree, [f"t{i}" for i in range(1, 11)])
    _write_fasta(aln, ["t1", "t2", "t3"] + [f"x{i}" for i in range(1, 8)])

    with pytest.raises(SystemExit):
        _load_prune_tree().main(str(tree), str(aln), str(out), str(log))

    assert not out.exists()
    assert "systematic tree<->alignment ID mismatch" in log.read_text(encoding="utf-8")


def test_too_few_tips_aborts(tmp_path):
    # Fewer than 3 surviving tips is unusable for HyPhy regardless of matching.
    tree = tmp_path / "tree.nwk"
    aln = tmp_path / "aln.fasta"
    out = tmp_path / "tree_M1.nwk"
    log = tmp_path / "prune.log"
    _write_star_tree(tree, ["t1", "t2", "t3", "t4", "t5"])
    _write_fasta(aln, ["t1", "t2"])

    with pytest.raises(SystemExit):
        _load_prune_tree().main(str(tree), str(aln), str(out), str(log))
    assert not out.exists()
