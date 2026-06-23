#!/usr/bin/env python3
import json
import sys
import subprocess
import tempfile
import os
import itertools
from pathlib import Path

CANDIDATES_JSON = "ssh-output/influenza-a-mp-indonesia-fix-1/m09/degenerate/segment-MP_region_1/final_candidates.json"
SURVEILLANCE_FASTA = "dataset/EpiFlu/MP/EpiFlu_A_MP_Combined_Sequence.fasta"
DESIGN_FASTA = "ssh-output/influenza-a-mp-indonesia-fix-1/m02/filtered_sequences.fasta"
MAX_MISMATCH = 3

def run_seqkit_single(fwd, rev, fasta, label=""):
    cmd = ["seqkit", "amplicon", "-F", fwd, "-R", rev, "-m", str(MAX_MISMATCH), "-j", "4", "--bed", fasta]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    lines = [l for l in r.stdout.splitlines() if l.strip() and not l.startswith("#")]
    seqs = set(l.split("\t")[0] for l in lines if "\t" in l)
    print(f"  [seqkit single {label}] found: {len(seqs)} sequences | rc={r.returncode}")
    return seqs

def run_seqkit_pool(fwd_list, rev_list, fasta, label=""):
    combos = list(itertools.product(fwd_list, rev_list))
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False, encoding="utf-8") as tmp:
        for i, (f, r) in enumerate(combos):
            tmp.write(f"comb_{i}\t{f}\t{r}\n")
        tmp_path = tmp.name

    cmd = ["seqkit", "amplicon", "-p", tmp_path, "-m", str(MAX_MISMATCH), "-j", "4", "--bed", fasta]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    os.remove(tmp_path)
    
    lines = [l for l in r.stdout.splitlines() if l.strip() and not l.startswith("#")]
    seqs = set(l.split("\t")[0] for l in lines if "\t" in l)
    print(f"  [seqkit pool ({len(combos)} combos) {label}] found: {len(seqs)} sequences | rc={r.returncode}")
    return seqs

def main():
    if not os.path.exists(CANDIDATES_JSON):
        print(f"Candidates file not found: {CANDIDATES_JSON}")
        return

    with open(CANDIDATES_JSON) as f:
        candidates = json.load(f)

    cand = candidates[0]  # Rank 1
    print(f"Candidate Rank 1: {cand.get('cand_id')}")
    fwd_maj = cand["fwd_seq"]
    rev_maj = cand["rev_seq"]
    fwd_pool = cand.get("fwd_possible_seqs", [])
    rev_pool = cand.get("rev_possible_seqs", [])

    print(f"Fwd majority: {fwd_maj} (len {len(fwd_maj)})")
    print(f"Rev majority: {rev_maj} (len {len(rev_maj)})")
    print(f"Fwd pool: {len(fwd_pool)} seqs, Rev pool: {len(rev_pool)} seqs")

    print("\nRunning Seqkit on Design FASTA:")
    if os.path.exists(DESIGN_FASTA):
        run_seqkit_single(fwd_maj, rev_maj, DESIGN_FASTA, "majority")
        if fwd_pool and rev_pool:
            run_seqkit_pool(fwd_pool, rev_pool, DESIGN_FASTA, "pool")
            fwd_pool_plus = list(set(fwd_pool + [fwd_maj]))
            rev_pool_plus = list(set(rev_pool + [rev_maj]))
            run_seqkit_pool(fwd_pool_plus, rev_pool_plus, DESIGN_FASTA, "pool+majority")
    else:
        print(f"Design FASTA not found: {DESIGN_FASTA}")

if __name__ == "__main__":
    main()
