# Dual-Target Multiplex Assay Compatibility Report
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


### Cross-Dimer Interactions for Rank 1 Combination

**Pair ID**: `FluA_Rank7_x_FluB_Rank6` (Score: 73.9%)

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
| Probe_A_vs_Probe_B | -6.76 kcal/mol | -1.29 kcal/mol |
