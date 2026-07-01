# Module M09: Specificity Check & IGV Output

Module M09 verifies the diagnostic specificity of the designed primers and probe against host genomes and other non-target pathogens. It also generates visualization tracks for IGV.

---

## BLAST Specificity Audit

To ensure the assay does not generate false positives, candidate primers and probes are run through a BLAST search:

1.  **Database Modes**:
    *   `remote`: Queries the NCBI `nt` database via the cloud API.
    *   `local` (Recommended): Queries local BLAST databases (e.g., Human GRCh38 genome and a database of common respiratory viruses like Adenovirus, RSV, and Rhinoviruses).
2.  **Cross-Reactivity Evaluation**:
    *   Candidates with hits showing $>80\%$ identity and $>80\%$ coverage to non-target templates are flagged.
    *   Special attention is paid to the 3' end of the primers; mismatches in the last 3 bases of the 3' end are heavily weighted, as they effectively prevent Taq or rekombinase extension.
    *   Candidates that pass are labeled as `PASS`, while candidates with minor cross-reactivity risks are labeled as `BORDERLINE`. High-risk candidates are discarded.

---

## GARD Breakpoint Cross-Validation

The module overlays the coordinates of the designed primers with the recombination breakpoints identified by **HyPhy GARD** in Module 05:

*   **Rule**: Primers and probes must not cross a recombination breakpoint.
*   **Rationale**: Recombination breakpoints are sites of frequent genetic exchange. Designing an assay that spans a breakpoint increases the risk of diagnostic failure due to sudden genomic reassortments.

---

## Integrative Genomics Viewer (IGV) Export

To allow visual validation of the candidates, the module automatically generates standard genomic tracks:

*   **`final_candidates.gff3`**: Contains GFF3 records for the forward primer, reverse primer, and probe.
*   **`final_candidates.bed`**: BED file representing the amplified amplicon region.
*   **`igv_session.xml`**: A pre-configured XML session file. 

### Loading the Session in IGV
1.  Open **IGV**.
2.  Select **File > Open Session...**
3.  Choose the `results/{build_name}/m09/igv_session.xml` file.
4.  IGV will automatically load the reference sequence, the candidate primer tracks, and the BED amplicon tracks, aligned together for review.
