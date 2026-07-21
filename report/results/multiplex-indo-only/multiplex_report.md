# Dual-Target Multiplex Assay Compatibility & Synthesis Report
## Target Pair: Influenza A MP (FAM) + Influenza B MP (HEX) — Indonesia Cohort

Screened **100 candidate combinations** ($10 \times 10$) across two temperature regimes ($39.0^\circ\text{C}$ RT-RAA & $55.0^\circ\text{C}$ qPCR).

### Top 5 Recommended Multiplex Combinations

| Rank | Candidate Pair | Compatibility Score | RAA ΔG (39°C) | qPCR ΔG (55°C) | ΔTm (°C) | Status (39°C) |
|---|---|---|---|---|---|---|
| 1 | **FluA_Rank7_x_FluB_Rank6** | **73.9%** | -7.13 kcal/mol | -5.01 kcal/mol | 0.29°C | `PASS` |
| 2 | **FluA_Rank7_x_FluB_Rank11** | **73.8%** | -7.13 kcal/mol | -5.01 kcal/mol | 0.31°C | `PASS` |
| 3 | **FluA_Rank7_x_FluB_Rank5** | **73.2%** | -7.13 kcal/mol | -5.01 kcal/mol | 0.45°C | `PASS` |
| 4 | **FluA_Rank4_x_FluB_Rank6** | **72.1%** | -7.13 kcal/mol | -5.01 kcal/mol | 0.74°C | `PASS` |
| 5 | **FluA_Rank6_x_FluB_Rank6** | **72.1%** | -7.13 kcal/mol | -5.01 kcal/mol | 0.74°C | `PASS` |


### Rank 1 Oligo Synthesis Specifications (`FluA_Rank7_x_FluB_Rank6`)

#### 🟩 Target A: Influenza A MP (FAM Channel)

| Oligo Name | Type | Degenerate Sequence (5' -> 3') | Length | Tm (°C) | GC (%) | Modification / Label |
|---|---|---|---|---|---|---|
| FluA_Fwd | Fwd Primer | `CTGACTAARGGRATTTTRGGRTTTGTATT` | 29 nt | 60.1°C | 27.6% | Unlabeled |
| FluA_Rev | Rev Primer | `AGCAGTTAAACTATATARGAARCTGAARAGR` | 31 nt | 60.4°C | 25.8% | Unlabeled |
| FluA_Probe | Exo-Probe | `ACCGTGCCCAGTGAGCGAGGACTGCAGCGTAGACGCTTTRTCCARAATGCC` | 51 nt | 81.6°C | 56.9% | 5'-FAM / 3'-BHQ1 / Internal-THF |

#### 🟨 Target B: Influenza B MP (HEX Channel)

| Oligo Name | Type | Degenerate Sequence (5' -> 3') | Length | Tm (°C) | GC (%) | Modification / Label |
|---|---|---|---|---|---|---|
| FluB_Fwd | Fwd Primer | `GAAATAYCTATAATGYTYGAACCATTYCAGA` | 31 nt | 60.4°C | 25.8% | Unlabeled |
| FluB_Rev | Rev Primer | `TTTATTGTCTCTTTRTTTGGMCCYTTTATTC` | 31 nt | 60.4°C | 25.8% | Unlabeled |
| FluB_Probe | Exo-Probe | `TYTTATCAGCTCTCCATTTCATGGCTTGGACAATAGGRCATTTRAATC` | 48 nt | 71.3°C | 35.4% | 5'-HEX / 3'-BHQ1 / Internal-THF |

### Cross-Dimer Inter-Assay Interactions

| Oligo Pair | Cross-Dimer ΔG (39°C) | Cross-Dimer ΔG (55°C) |
|---|---|---|
| Fwd_A_vs_Fwd_B | -4.14 kcal/mol | -2.05 kcal/mol |
| Fwd_A_vs_Rev_B | -7.13 kcal/mol | -5.01 kcal/mol |
| Fwd_A_vs_Probe_B | -3.09 kcal/mol | -1.53 kcal/mol |
| Rev_A_vs_Fwd_B | -5.95 kcal/mol | -3.06 kcal/mol |
| Rev_A_vs_Rev_B | -4.37 kcal/mol | -2.27 kcal/mol |
| Rev_A_vs_Probe_B | -5.69 kcal/mol | -3.49 kcal/mol |
| Probe_A_vs_Fwd_B | -4.77 kcal/mol | -2.83 kcal/mol |
| Probe_A_vs_Rev_B | -6.12 kcal/mol | -3.49 kcal/mol |
| Probe_A_vs_Probe_B | -6.9 kcal/mol | -4.22 kcal/mol |
