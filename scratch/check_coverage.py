#!/usr/bin/env python3
import re
import os

# IUPAC table
IUPAC_GROUPS = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'}, 'U': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'}, 'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'}, '-': {'-'}
}

def match_base(p_base: str, a_base: str) -> bool:
    p_base = p_base.upper()
    a_base = a_base.upper()
    p_set = IUPAC_GROUPS.get(p_base, {p_base})
    a_set = IUPAC_GROUPS.get(a_base, {a_base})
    return not p_set.isdisjoint(a_set)

def reverse_complement(seq: str) -> str:
    complement = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
        'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N', '-': '-'
    }
    return "".join(complement.get(base.upper(), base) for base in reversed(seq))

def min_mismatches_sliding_window(probe_seq: str, amplicon_seq: str) -> int:
    probe_len = len(probe_seq)
    amp_len = len(amplicon_seq)
    if amp_len < probe_len:
        return probe_len
    
    best_mismatch = probe_len
    for p_seq in [probe_seq, reverse_complement(probe_seq)]:
        for i in range(amp_len - probe_len + 1):
            window = amplicon_seq[i:i + probe_len]
            mismatches = 0
            for j in range(probe_len):
                if not match_base(p_seq[j], window[j]):
                    mismatches += 1
                    if mismatches >= best_mismatch:
                        break
            if mismatches < best_mismatch:
                best_mismatch = mismatches
            if best_mismatch == 0:
                return 0
    return best_mismatch

def main():
    fasta_path = "ssh-output/influenza-a-mp-indonesia-fix-1/m02/filtered_sequences.fasta"
    if not os.path.exists(fasta_path):
        print(f"File not found: {fasta_path}")
        return

    # Load candidates degenerate (Cand 1)
    # Fwd: CTGACTAARGGRATTTTRGGRTTTGTATTCA
    # Rev: GCAGTTAAACTATATARGAARCTGAARAGRG
    # Let's check both as majority (ATGC) and degenerate (IUPAC)
    fwd_deg = "CTGACTAARGGRATTTTRGGRTTTGTATTCA"
    rev_deg = "GCAGTTAAACTATATARGAARCTGAARAGRG" # Reverse: CYCTYTTCAGYTTYCTATATAGTTTAACTGC
    
    fwd_maj = "CTGACTAAGGGGATTTTAGGATTTGTATTCAC"
    rev_maj = "TCTCTTTTCAGCTTCTTATATAGTTTAACTGC"

    # Read FASTA
    sequences = []
    current_seq = []
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith(">"):
                if current_seq:
                    sequences.append("".join(current_seq))
                current_seq = []
            else:
                current_seq.append(line.strip().upper())
        if current_seq:
            sequences.append("".join(current_seq))

    print(f"Loaded {len(sequences)} sequences from Design FASTA.")

    # Let's count coverage for degenerate
    fwd_deg_ok = 0
    rev_deg_ok = 0
    both_deg_ok = 0

    fwd_maj_ok = 0
    rev_maj_ok = 0
    both_maj_ok = 0

    for i, seq in enumerate(sequences):
        # Degenerate
        f_deg_mm = min_mismatches_sliding_window(fwd_deg, seq)
        r_deg_mm = min_mismatches_sliding_window(reverse_complement(rev_deg), seq)
        
        # Majority
        f_maj_mm = min_mismatches_sliding_window(fwd_maj, seq)
        r_maj_mm = min_mismatches_sliding_window(reverse_complement(rev_maj), seq)

        if f_deg_mm <= 3:
            fwd_deg_ok += 1
        if r_deg_mm <= 3:
            rev_deg_ok += 1
        if f_deg_mm <= 3 and r_deg_mm <= 3:
            both_deg_ok += 1

        if f_maj_mm <= 3:
            fwd_maj_ok += 1
        if r_maj_mm <= 3:
            rev_maj_ok += 1
        if f_maj_mm <= 3 and r_maj_mm <= 3:
            both_maj_ok += 1

    print("\nDegenerate (IUPAC-aware, max 3 mm):")
    print(f"  Fwd degenerate covers: {fwd_deg_ok}/{len(sequences)} ({fwd_deg_ok/len(sequences)*100:.2f}%)")
    print(f"  Rev degenerate covers: {rev_deg_ok}/{len(sequences)} ({rev_deg_ok/len(sequences)*100:.2f}%)")
    print(f"  Both degenerate cover: {both_deg_ok}/{len(sequences)} ({both_deg_ok/len(sequences)*100:.2f}%)")

    print("\nMajority (Non-degenerate, max 3 mm):")
    print(f"  Fwd majority covers:   {fwd_maj_ok}/{len(sequences)} ({fwd_maj_ok/len(sequences)*100:.2f}%)")
    print(f"  Rev majority covers:   {rev_maj_ok}/{len(sequences)} ({rev_maj_ok/len(sequences)*100:.2f}%)")
    print(f"  Both majority cover:   {both_maj_ok}/{len(sequences)} ({both_maj_ok/len(sequences)*100:.2f}%)")

if __name__ == "__main__":
    main()
