import os
import json
import shutil
import pandas as pd
import yaml

def reverse_complement(seq: str) -> str:
    complement = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
        'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N', '-': '-'
    }
    return "".join(complement.get(base.upper(), base) for base in reversed(seq))

def match_base(p_base: str, a_base: str) -> bool:
    iupac_groups = {
        'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'}, 'U': {'T'},
        'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
        'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'},
        'D': {'A', 'G', 'T'}, 'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
        'N': {'A', 'C', 'G', 'T'}, '-': {'-'}
    }
    p_set = iupac_groups.get(p_base.upper(), {p_base.upper()})
    a_set = iupac_groups.get(a_base.upper(), {a_base.upper()})
    return not p_set.isdisjoint(a_set)

def audit_primer(profile, start, end, seq, maj_seq, is_rev=False):
    # Slice the profile dataframe
    slice_df = profile[(profile["position"] >= start) & (profile["position"] < end)].copy()
    
    rows = []
    # If the primer is reverse-oriented, let's reverse the primer index for display
    # but keep alignment positions.
    for i, (idx, row) in enumerate(slice_df.iterrows()):
        pos = int(row["position"])
        cons = row["consensus_base"]
        iupac = row["iupac_base"]
        
        freqs = {
            "A": row["freq_A"],
            "T": row["freq_T"],
            "G": row["freq_G"],
            "C": row["freq_C"]
        }
        
        # Determine primer relative position (5' -> 3')
        if is_rev:
            primer_pos = (end - 1) - pos
        else:
            primer_pos = pos - start
            
        is_degenerate = iupac != cons
        
        rows.append({
            "pos": pos,
            "primer_pos": primer_pos,
            "cons": cons,
            "iupac": iupac,
            "A": freqs["A"],
            "T": freqs["T"],
            "G": freqs["G"],
            "C": freqs["C"],
            "is_degenerate": is_degenerate
        })
        
    # Sort by primer position 5' -> 3'
    rows.sort(key=lambda x: x["primer_pos"])
    return rows

def generate_markdown_table(rows):
    lines = []
    lines.append("| Primer Pos (5' → 3') | Segment Position | Majority Base | IUPAC Base | A (%) | T (%) | G (%) | C (%) | Decision Rationale |")
    lines.append("| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |")
    
    for r in rows:
        rationale = ""
        if r["is_degenerate"]:
            cleared = []
            for base in ["A", "T", "G", "C"]:
                pct = r[base] * 100
                if pct >= 15.0:
                    cleared.append(f"**{base}** ({pct:.1f}%)")
            rationale = f"**Degenerate**: Cleared &ge;15% MAF threshold: " + ", ".join(cleared)
            iupac_str = f"**{r['iupac']}**"
        else:
            vars_present = []
            for base in ["A", "T", "G", "C"]:
                pct = r[base] * 100
                if 0.0 < pct < 15.0:
                    vars_present.append(f"{base} ({pct:.1f}%)")
            if vars_present:
                rationale = f"Conserved (minor variation &lt;15% ignored: " + ", ".join(vars_present) + ")"
            else:
                rationale = "100% Conserved"
            iupac_str = r["iupac"]
            
        lines.append(f"| {r['primer_pos']} | {r['pos']} | {r['cons']} | {iupac_str} | {r['A']*100:.1f}% | {r['T']*100:.1f}% | {r['G']*100:.1f}% | {r['C']*100:.1f}% | {rationale} |")
        
    return "\n".join(lines)

def main():
    # 1. Source and Dest Paths
    src_profile = "ssh-output/results/influenza-a-mp-indo-only/m07/entropy_profile.tsv"
    dest_dir = "report/results/influenza-a-mp-indo-only"
    dest_profile = os.path.join(dest_dir, "entropy_profile.tsv")
    
    primers_path = "ssh-output/results/influenza-a-mp-indo-only/m08/degenerate/segment-MP_region_1/best_primers.json"
    manifest_path = "ssh-output/results/influenza-a-mp-indo-only/m07/primer_targets_manifest.yaml"
    
    # Create dest folder if not exists
    os.makedirs(dest_dir, exist_ok=True)
    
    # Copy profile TSV for serving
    print(f"Copying {src_profile} to {dest_profile}...")
    shutil.copy(src_profile, dest_profile)
    
    # Load profile data
    profile = pd.read_csv(src_profile, sep="\t")
    
    # Load best primers
    with open(primers_path, "r") as f:
        primers = json.load(f)
        
    fwd = primers["forward"]
    rev = primers["reverse"]
    probe = primers["probe"]
    
    # Perform audits
    fwd_audit = audit_primer(profile, fwd["start"], fwd["end"], fwd["sequence"], fwd["majority_sequence"], is_rev=False)
    rev_audit = audit_primer(profile, rev["start"], rev["end"], rev["sequence"], rev["majority_sequence"], is_rev=True)
    
    # Audit Probe
    probe_audit = audit_primer(profile, probe["start"], probe["end"], probe["sequence"], probe["sequence"], is_rev=False)
    
    # Write supplementary.qmd
    qmd_content = f"""---
title: "Supplementary Material: IUPAC Degeneracy Audit"
subtitle: "A detailed biochemical validation appendix for Candidate 1 primers & probe"
---

::: {{.callout-note}}
## Purpose of this Appendix
This supplementary document provides a granular, nucleotide-by-nucleotide audit of the degenerate IUPAC positions within the designed primers and probe for **Candidate 1** (designed under the **Indonesia Only** evolutionary model). All figures are derived directly from the Multiple Sequence Alignment (MSA) of Indonesian H5N1, H3N2, and H1N1 isolates, and are recorded in the raw database:
[Download Raw Entropy Profile TSV (entropy_profile.tsv)](results/influenza-a-mp-indo-only/entropy_profile.tsv)
:::

## 1. IUPAC Degeneracy Decision Rules

In our pipeline, designed primers are back-mapped to the alignment consensus. Ambiguity codes (degeneracy) are selected using a strict minor-allele frequency (MAF) cutoff:

1. **Unambiguous Base (A/T/G/C)**: If a single base occurs at a frequency $> 85\\%$ (minor alleles $< 15\\%$) at a given position, it is written as a standard nucleotide. Minor variations under $15\\%$ are ignored to prevent excessive degeneracy multipliers.
2. **Degenerate Base (IUPAC)**: If two or more nucleotides at a single position each clear the **$15\\%$ threshold**, the position is mapped to the appropriate IUPAC degeneracy code (e.g., `R` for A/G, `Y` for C/T, `D` for A/G/T).

### 1.1. Understanding the Degeneracy Multiplier

A degenerate primer is not a single sequence, but a physical mixture of multiple unique primer sequences synthesized together in the same tube. The **Degeneracy Multiplier** (e.g., `16x`) represents the total number of unique sequence combinations resulting from the IUPAC positions in the primer.

It is calculated by multiplying the number of possible bases at each degenerate position:
$$M = \prod_{{i=1}}^{{L}} \\text{{options}}(c_i)$$

Where:
* Standard bases (`A`, `T`, `G`, `C`) have **1** option.
* Two-base degeneracies (`R`, `Y`, `W`, `S`, `K`, `M`) have **2** options.
* Three-base degeneracies (`B`, `D`, `H`, `V`) have **3** options.
* Four-base degeneracies (`N`) have **4** options.

For example, a forward primer with sequence `CTGACTAARGGRATTTTRGGRTTTGTATTCA` has four `R` positions (each allowing `A` or `G`). This results in a multiplier of $2 \\times 2 \\times 2 \\times 2 = 16\\text{{x}}$. Thus, the final tube contains exactly 16 distinct primer variations in equal proportions.

**Why limit the multiplier?**
* **Primer Dilution**: If the multiplier is too high (e.g., $>64\\text{{x}}$), the concentration of any single primer variant in the reaction buffer becomes too low, which can severely reduce amplification efficiency.
* **Non-Specific Amplification**: High degeneracy increases the risk of off-target binding and primer-dimer formation.
To balance universality and performance, the pipeline enforces a maximum degeneracy limit of **$16\\text{{x}}$** for RAA designs.

---

## 2. Forward Primer Audit

* **ID**: `{fwd['id']}`
* **Position on Segment**: `{fwd['start']} - {fwd['end']} bp` (0-based)
* **Degenerate Sequence (5' &rarr; 3')**: `{fwd['sequence']}`
* **Majority Sequence (5' &rarr; 3')**: `{fwd['majority_sequence']}`
* **Degeneracy Multiplier**: `{fwd['degeneracy_multiplier']}x`

::: {{.small-table}}
{generate_markdown_table(fwd_audit)}
:::

---

## 3. Reverse Primer Audit

* **ID**: `{rev['id']}`
* **Position on Segment**: `{rev['start']} - {rev['end']} bp` (0-based)
* **Degenerate Sequence (5' &rarr; 3')**: `{rev['sequence']}`
* **Majority Sequence (5' &rarr; 3')**: `{rev['majority_sequence']}`
* **Degeneracy Multiplier**: `{rev['degeneracy_multiplier']}x`

::: {{.small-table}}
{generate_markdown_table(rev_audit)}
:::

---

## 4. Probe Audit

* **ID**: `{probe['id']}`
* **Position on Segment**: `{probe['start']} - {probe['end']} bp` (0-based)
* **Degenerate Sequence (5' &rarr; 3')**: `{probe['sequence']}`
* **Degeneracy Multiplier**: `{probe['degeneracy_multiplier']}x`

::: {{.small-table}}
{generate_markdown_table(probe_audit)}
:::
"""
    
    with open("report/supplementary.qmd", "w", encoding="utf-8") as f:
        f.write(qmd_content)
    print("Successfully generated report/supplementary.qmd")

if __name__ == "__main__":
    main()
