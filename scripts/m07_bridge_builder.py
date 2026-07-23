"""
M07 Bridge Builder v2.2 — Conservation & Selection Integrator

Merges outputs from two independent analysis tracks:
  Track A (Augur): MSA conservation via Shannon Entropy
  Track B (HyPhy): Evolutionary selection pressure (FUBAR/MEME/GARD)

Version 2.2 adds:
  - IUPAC-aware consensus sequence (iupac_consensus_sequence).
    Positions where a minor allele frequency >= IUPAC_THRESHOLD are encoded
    using standard IUPAC ambiguity codes (R, Y, W, S, K, M, B, D, H, V, N)
    instead of the majority-vote base. This is critical for RAA primer design
    where degenerate positions must be known and avoided.
  - Per-position base frequency data stored in entropy profile TSV.
  - Three sequence representations now exported in all manifests:
      consensus_sequence       — majority-vote per column (plain ATGC)
      iupac_consensus_sequence — IUPAC-encoded where ambiguous
      reference_sequence       — raw slice from reference genome

Version 2.1 adds:
  - 150 bp reference consensus sequence embedded in all output manifests.
  - Formal scientific grading annotations (replacing informal descriptors).
  - Automatic CSV export alongside JSON/YAML for spreadsheet-based review.

Inputs: MSA, GFF, HyPhy JSONs, Reference Fasta
Outputs:
  - conserved_regions.bed             (0-based half-open)
  - bridge_report.json                (0-based half-open)
  - primer_targets_manifest.yaml
  - entropy_profile.tsv               (now includes iupac_base + base freqs)
  - candidate_targets_summary.csv     (flat table for PI/Excel review)
"""

import argparse
import csv
import json
import logging
import math
import statistics
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import yaml
from Bio import SeqIO

# --- Step 1: Conservation & Profile Generation ---

# Minimum minor-allele frequency to trigger IUPAC ambiguity encoding.
# Set to 0.15 (15%) for RAA diagnostic sensitivity: a 15% variant means the
# primer would mismatch on ~15% of circulating strains, which is unacceptable
# for universal diagnostic assays.  Lower = more positions flagged = safer.
# Global parameters will be overridden by args
IUPAC_THRESHOLD = 0.15

# Maximum Shannon entropy for the documented 5-symbol alphabet {A,T,G,C,-}
# (see docs/concepts/shannon-entropy.md). Used to normalise entropy to a
# conservation score C = 1 - H / MAX_ENTROPY_BITS in [0, 1].
MAX_ENTROPY_BITS = math.log2(5)  # ~= 2.3219

# Significance threshold for MEME episodic-selection sites. HyPhy recommends
# p <= 0.1 rather than the conventional 0.05, because MEME is a conservative test
# on small alignments (hyphy.org selection tutorial). Overridable via
# config hyphy.p_value_threshold; FUBAR uses a posterior probability instead
# (>= 0.9, see bridge.tiers.*.min_purifying_prob).
MEME_PVALUE_THRESHOLD = 0.1

# Full IUPAC ambiguity lookup: frozenset of observed bases -> IUPAC code.
# Only 2-base and 3-base combos are listed; 4-base = N.
_IUPAC_MAP: dict = {
    frozenset("AG"): "R",  # puRine
    frozenset("CT"): "Y",  # pYrimidine
    frozenset("AT"): "W",  # Weak
    frozenset("CG"): "S",  # Strong
    frozenset("GT"): "K",  # Keto
    frozenset("AC"): "M",  # aMino
    frozenset("CGT"): "B", # not A
    frozenset("AGT"): "D", # not C
    frozenset("ACT"): "H", # not G
    frozenset("ACG"): "V", # not T
    frozenset("ACGT"): "N",
}


def sanitize_region_id(gene_name: str) -> str:
    """Return a shell-safe region_id by replacing characters that break bash
    quoting or filesystem paths.  Specifically:
      - single-quote (') -> p  (e.g. 3'UTR -> 3pUTR, 5'UTR -> 5pUTR)
      - double-quote (") -> removed
      - whitespace      -> underscore
      - any remaining non-alphanumeric except - and _ -> removed
    This is applied at manifest write time so all downstream tools
    (M08, M09, M10) always receive safe filesystem-friendly names.
    """
    import re
    s = gene_name.replace("'", "p").replace('"', "")
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_\-]", "", s)
    return s


def _to_iupac_base(freq: dict, threshold: float = IUPAC_THRESHOLD) -> str:
    """Return IUPAC ambiguity code for a position given per-base frequencies.

    If every base other than the majority is below *threshold*, return the
    majority base (plain A/T/G/C).  If two or more bases each clear the
    threshold, return the appropriate IUPAC code.  Falls back to 'N' for
    fully ambiguous positions.

    Args:
        freq: Dict mapping base (A/T/G/C) to fractional frequency (0-1).
        threshold: Minor-allele frequency cutoff for IUPAC encoding.

    Returns:
        Single-character string: plain base or IUPAC ambiguity code.
    """
    # Collect bases that are meaningfully present.
    present = {b for b, f in freq.items() if f >= threshold}
    if len(present) == 0:
        return "-"
    if len(present) == 1:
        return present.pop()  # unambiguous majority
    key = frozenset(present)
    return _IUPAC_MAP.get(key, "N")


def build_entropy_profile(records):
    """Calculate Shannon entropy, majority-vote consensus, and IUPAC base
    for every column in a multiple sequence alignment (MSA).

    Returns a list of per-position dicts with keys:
        position       — 0-based column index
        consensus_base — majority-vote base (plain A/T/G/C or '-')
        iupac_base     — IUPAC-encoded base if minor allele >= IUPAC_THRESHOLD,
                         else same as consensus_base
        entropy        — Shannon entropy in bits over the 5-symbol alphabet
                         {A,T,G,C,-} (0 = fully conserved, log2(5) ~= 2.322 = max);
                         gaps and other non-ACGT symbols count toward variability
                         so deletions correctly drive conservation down.
        freq_A, freq_T, freq_G, freq_C — per-base frequencies (gap-excluded,
                         used for consensus / IUPAC degenerate-base decisions)
    """
    if not records:
        return [], 0, 0

    n_sequences = len(records)
    aln_length  = len(records[0].seq)

    _empty_freq = {"freq_A": 0.0, "freq_T": 0.0, "freq_G": 0.0, "freq_C": 0.0}

    profile = []
    for i in range(aln_length):
        column = [str(rec.seq[i]).upper() for rec in records]
        bases  = [b for b in column if b in "ACGT"]
        n      = len(bases)                       # ACGT-only count
        gap_like = n_sequences - n                # gaps + N + other ambiguity

        if n == 0:
            # No usable base in ANY sequence (all gap/N): undesignable column.
            # Assign maximum entropy -> conservation score 0 so the region is
            # disqualified rather than mislabelled as a perfectly conserved gap.
            profile.append({
                "position": i,
                "consensus_base": "-",
                "iupac_base": "-",
                "entropy": MAX_ENTROPY_BITS,
                **_empty_freq,
            })
            continue

        counter = Counter(bases)
        consensus = counter.most_common(1)[0][0]

        # ACGT-only frequencies drive consensus / IUPAC degeneracy and the
        # reported freq_* columns (documented as gap-excluded).
        freq = {b: c / n for b, c in counter.items()}
        iupac = _to_iupac_base(freq)

        # Entropy uses the documented gap-inclusive 5-symbol alphabet: the
        # denominator is ALL sequences and every non-ACGT symbol (gap, N,
        # ambiguity) is folded into a single "gap" state, since like a deletion
        # it cannot be relied on for primer hybridisation.
        entropy = 0.0
        ent_counts = list(counter.values())
        if gap_like > 0:
            ent_counts.append(gap_like)
        for c in ent_counts:
            p = c / n_sequences
            if p > 0:
                entropy -= p * math.log2(p)

        profile.append({
            "position":      i,
            "consensus_base": consensus,
            "iupac_base":    iupac,
            "entropy":       entropy,
            "freq_A":        round(freq.get("A", 0.0), 4),
            "freq_T":        round(freq.get("T", 0.0), 4),
            "freq_G":        round(freq.get("G", 0.0), 4),
            "freq_C":        round(freq.get("C", 0.0), 4),
        })

    return profile, n_sequences, aln_length


def export_entropy_profile(profile, path):
    """Write the full per-position entropy profile to a TSV file.

    Columns: position, consensus_base, iupac_base, entropy,
             freq_A, freq_T, freq_G, freq_C.
    """
    fieldnames = [
        "position", "consensus_base", "iupac_base", "entropy",
        "freq_A", "freq_T", "freq_G", "freq_C",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(profile)
    print(f"[M07] Wrote entropy profile: {path}")

def scan_conservation(profile, window_size, entropy_threshold):
    """Scan pre-computed entropy profile for conserved regions using sliding-window."""
    aln_length = len(profile)
    entropies = [p["entropy"] for p in profile]
    regions = []
    
    for i in range(aln_length - window_size + 1):
        window = entropies[i:i + window_size]
        mean_e = sum(window) / len(window)
        if mean_e < entropy_threshold:
            regions.append({
                "start": i,
                "end": i + window_size,
                "mean_entropy": mean_e,
                "conservation_pct": round((1 - mean_e / MAX_ENTROPY_BITS) * 100, 2),
                "length": window_size
            })
            
    if not regions:
        print(f"[M07] WARNING: No windows with entropy < {entropy_threshold}")
        
    print(f"[M07] Found {len(regions)} conserved raw windows (threshold={entropy_threshold})")
    return regions

# --- Biomechanical Scanners ---

def extract_reference_subsequence(ref_seq, start, end):
    """Extract subsequence from the raw reference genome.
    Note: if ref lacks gaps relative to MSA, mapping should ideally map column to ref pos.
    For simplicity, assuming 1:1 loosely aligned or gaps stripped prior.
    """
    return str(ref_seq[start:end]).upper()

def calculate_gc_content(sequence):
    length = len(sequence)
    if length == 0:
        return 0.0
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / length) * 100.0

def detect_homopolymer(sequence, max_repeats=5):
    for base in "ACGT":
        if (base * (max_repeats + 1)) in sequence:
            return True
    return False

def find_max_contiguous_block(profile, start, end, threshold=0.10):
    max_len = 0
    current_len = 0
    for i in range(start, end):
        if profile[i]["entropy"] < threshold:
            current_len += 1
            if current_len > max_len:
                max_len = current_len
        else:
            current_len = 0
    return max_len

def get_mismatch_positions(profile, start, end, threshold=0.10):
    mismatches = []
    for i in range(start, end):
        if profile[i]["entropy"] > threshold:
            mismatches.append(i)
    return mismatches

def parse_gff(gff_path):
    """Parse GFF3 file to extract gene and UTR coordinates (0-based half-open).

    Supports spliced genes by grouping CDS segments by name under 'segments' coordinate lists.
    """
    features_by_gene = {}  # {gene_name: {"CDS": [], "gene": [], "chrom": chrom, "strand": strand}}
    utrs = []

    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            chrom = parts[0]
            feature_type = parts[2]
            start = int(parts[3]) - 1  # GFF is 1-based → convert to 0-based
            end = int(parts[4])        # Keep as exclusive end
            strand = parts[6]

            # Parse attributes into dict
            attrs = {}
            for attr in parts[8].split(";"):
                attr = attr.strip()
                if "=" in attr:
                    key, val = attr.split("=", 1)
                    attrs[key] = val.strip('"')

            if feature_type in ("five_prime_UTR", "three_prime_UTR"):
                utr_name = "5'UTR" if "five" in feature_type else "3'UTR"
                utrs.append({"name": utr_name, "start": start, "end": end, "chrom": chrom})
                continue

            if feature_type in ("gene", "CDS"):
                name = attrs.get("gene") or attrs.get("Name") or attrs.get("ID") or ""
                if not name:
                    continue
                if name not in features_by_gene:
                    features_by_gene[name] = {"CDS": [], "gene": [], "chrom": chrom, "strand": strand}
                features_by_gene[name][feature_type].append((start, end))

    genes = []
    for name, data in features_by_gene.items():
        segs = data["CDS"] if data["CDS"] else data["gene"]
        if not segs:
            continue
        # Sort segments by start coordinate
        segs.sort(key=lambda x: x[0])
        
        gene_start = min(s[0] for s in segs)
        gene_end = max(s[1] for s in segs)
        
        genes.append({
            "name": name,
            "start": gene_start,
            "end": gene_end,
            "chrom": data["chrom"],
            "strand": data["strand"],
            "segments": segs
        })

    return genes, utrs

def find_gene_for_position(pos, genes, utrs):
    # First check if the position is inside any CDS segment of a gene.
    # Prioritize genes where the position is actually in a CDS segment,
    # and if multiple, sort by length (smallest first) to handle overlapping genes.
    coding_genes = []
    for g in genes:
        for seg_start, seg_end in g.get("segments", [(g["start"], g["end"])]):
            if seg_start <= pos < seg_end:
                coding_genes.append(g)
                break
    if coding_genes:
        coding_genes.sort(key=lambda x: x["end"] - x["start"])
        return coding_genes[0]["name"], True

    # Fallback to general span (e.g. if position is in UTR or intron)
    for g in sorted(genes, key=lambda x: x["end"] - x["start"]):
        if g["start"] <= pos < g["end"]:
            return g["name"], True

    for u in sorted(utrs, key=lambda x: x["end"] - x["start"]):
        if u["start"] <= pos < u["end"]:
            return u["name"], False

    return "intergenic", False

# --- HyPhy Score Mapper ---

def parse_fubar_json(filepath):
    """Extract per-codon FUBAR estimates.

    Returns entries with three fields, which serve different purposes:

      fubar_purifying_prob — Prob[alpha>beta], the posterior probability that the
        codon is under purifying selection AT THE PROTEIN LEVEL. Used by the tier
        system, where the conventional >= 0.9 cut-off applies.

      alpha — mean posterior SYNONYMOUS substitution rate.
      beta  — mean posterior NONSYNONYMOUS substitution rate.

    The rates were previously discarded, and that mattered. A primer does not read
    protein: what breaks it is nucleotide change at any of the three codon
    positions. Under strong purifying selection the nonsynonymous channel is
    suppressed while the synonymous one — concentrated at the wobble position — is
    not, so a codon can have a high Prob[alpha>beta] and still turn over its
    nucleotides rapidly. Ranking primer targets on Prob[alpha>beta] alone therefore
    scores a fast-drifting site identically to a genuinely static one.
    """
    with open(filepath, "r") as f:
        data = json.load(f)
    mle = data.get("MLE", {})
    headers = mle.get("headers", [])
    rows = mle.get("content", {}).get("0", [])

    def _find(*predicates, fallback=None):
        for i, h in enumerate(headers):
            label = (h[0] if isinstance(h, list) else str(h)).lower().strip()
            if any(pred(label) for pred in predicates):
                return i
        return fallback

    # Match Prob[alpha>beta] before the bare "alpha" column, or the substring test
    # would claim the wrong index.
    purifying_idx = _find(lambda l: "alpha>beta" in l.replace(" ", ""), fallback=5)
    alpha_idx = _find(lambda l: l == "alpha", lambda l: l.startswith("alpha "),
                      fallback=0)
    beta_idx = _find(lambda l: l == "beta", lambda l: l.startswith("beta "),
                     fallback=1)

    results = []
    for i, row in enumerate(rows):
        results.append({
            "codon_idx": i,
            "fubar_purifying_prob": row[purifying_idx] if purifying_idx < len(row) else 0.5,
            "alpha": row[alpha_idx] if alpha_idx is not None and alpha_idx < len(row) else None,
            "beta": row[beta_idx] if beta_idx is not None and beta_idx < len(row) else None,
        })
    return results


# Average number of synonymous / nonsynonymous nucleotide sites in a codon. A
# codon has 3 positions; roughly three quarters of possible single-nucleotide
# changes are nonsynonymous. Used to convert FUBAR's two rates into one expected
# nucleotide substitution rate, which is the quantity a primer actually cares
# about.
SYN_SITE_FRACTION = 0.25


def codon_substitution_rate(alpha, beta):
    """Expected nucleotide substitution rate for a codon, or None if unmeasured.

    alpha and beta are per-site rates for the synonymous and nonsynonymous
    channels; weighting them by how many of the codon's nucleotide sites belong to
    each channel gives the rate at which ANY of its three bases changes -- which is
    what determines whether a primer overlapping it keeps matching.
    """
    if alpha is None or beta is None:
        return None
    return SYN_SITE_FRACTION * alpha + (1.0 - SYN_SITE_FRACTION) * beta

def parse_meme_json(filepath):
    with open(filepath, "r") as f:
        data = json.load(f)
    mle = data.get("MLE", {})
    headers = mle.get("headers", [])
    rows = mle.get("content", {}).get("0", [])

    pval_idx = None
    for i, h in enumerate(headers):
        label = h[0] if isinstance(h, list) else str(h)
        if "p-value" in label.lower() or "p value" in label.lower():
            pval_idx = i
            break
    if pval_idx is None:
        pval_idx = len(headers) - 1

    results = []
    for i, row in enumerate(rows):
        pval = row[pval_idx] if pval_idx < len(row) else 1.0
        results.append({"codon_idx": i, "meme_pvalue": pval})
    return results

def parse_gard_json(filepath):
    with open(filepath, "r") as f:
        data = json.load(f)
    if data.get("skipped", False):
        return []
    partitions = data.get("data partitions", data.get("dataPartitions", {}))
    if not partitions or len(partitions) <= 1:
        return []

    boundaries = []
    for key in sorted(partitions.keys(), key=lambda x: int(x) if x.isdigit() else x):
        coverage = partitions[key].get("coverage", [])
        for span in coverage:
            if isinstance(span, list) and len(span) == 2:
                boundaries.append({"start": span[0], "end": span[1]})

    breakpoints = []
    for i in range(1, len(boundaries)):
        breakpoints.append(boundaries[i]["start"])
    return breakpoints

def map_hyphy_to_genome(genes, hyphy_dir, target_genes):
    codon_data = {}       # {gene_name: {genomic_pos: entry}}
    gene_median_rate = {} # {gene_name: median codon substitution rate or None}
    gard_breakpoints = []

    for gene in genes:
        name = gene["name"]
        if name not in target_genes:
            continue

        # Construct cds_coords: consecutive nucleotide coordinates across all CDS segments
        segments = gene.get("segments", [(gene["start"], gene["end"])])
        cds_coords = []
        for start, end in segments:
            cds_coords.extend(range(start, end))

        if gene.get("strand") == "-":
            cds_coords.reverse()

        fubar_path = hyphy_dir / f"fubar_{name}.json"
        fubar_data = []
        if fubar_path.exists():
            try: fubar_data = parse_fubar_json(str(fubar_path))
            except Exception: pass

        meme_path = hyphy_dir / f"meme_{name}.json"
        meme_data = []
        if meme_path.exists():
            try: meme_data = parse_meme_json(str(meme_path))
            except Exception: pass

        if fubar_data or meme_data:
            codon_data[name] = {}
            n_codons = max(len(fubar_data), len(meme_data))
            for ci in range(n_codons):
                # Stop mapping if codon coordinates exceed GFF segment lengths (prevent M1/M2 overflow)
                if ci * 3 + 2 >= len(cds_coords):
                    break

                entry = {
                    "gene": name,
                    "codon_idx": ci,
                    "fubar_purifying_prob": 0.5,
                    "meme_pvalue": 1.0,
                    "subst_rate": None,
                }
                if ci < len(fubar_data):
                    entry["fubar_purifying_prob"] = fubar_data[ci]["fubar_purifying_prob"]
                    entry["subst_rate"] = codon_substitution_rate(
                        fubar_data[ci].get("alpha"), fubar_data[ci].get("beta"))
                if ci < len(meme_data):
                    entry["meme_pvalue"] = meme_data[ci]["meme_pvalue"]

                # Map the 3 nucleotides of this codon to genomic coordinates
                for offset in range(3):
                    gpos = cds_coords[ci * 3 + offset]
                    codon_data[name][gpos] = entry

            # FUBAR rates are relative to this alignment's own scaling, so an
            # absolute cut-off would not transfer between genes or builds. Record
            # the gene's median so a region can be judged against its own gene.
            rates = sorted({id(e): e for e in codon_data[name].values()}.values(),
                           key=lambda e: e["codon_idx"])
            measured = [e["subst_rate"] for e in rates if e["subst_rate"] is not None]
            gene_median_rate[name] = (statistics.median(measured) if measured else None)
            if measured:
                print(f"[M07] {name}: median codon substitution rate = "
                      f"{gene_median_rate[name]:.4f} over {len(measured)} codons")
            else:
                print(f"[M07] {name}: no FUBAR rate data — durability multiplier "
                      f"will fall back to the no-evidence value")

        gard_path = hyphy_dir / f"gard_{name}.json"
        if gard_path.exists():
            try:
                bp_codons = parse_gard_json(str(gard_path))
                for bp_codon in bp_codons:
                    if bp_codon * 3 < len(cds_coords):
                        gard_breakpoints.append({
                            "genome_pos": cds_coords[bp_codon * 3],
                            "gene": name,
                            "codon_idx": bp_codon,
                        })
            except Exception: pass

    return codon_data, gard_breakpoints, gene_median_rate

# --- Multipliers Computation ---

def compute_meme_multiplier(sites, config):
    if sites >= config.get("veto_above", 4):
        return 0.0
    for t in config.get("thresholds", []):
        if sites <= t.get("max_sites", 0):
            return t.get("multiplier", 1.0)
    return 1.0

def compute_gc_multiplier(gc_pct, config):
    if gc_pct < config.get("veto_below", 25) or gc_pct > config.get("veto_above", 70):
        return 0.0
    if config.get("ideal_min", 35) <= gc_pct <= config.get("ideal_max", 60):
        return 1.0
    marginal = config.get("marginal_ranges", [])
    for m in marginal:
        if m.get("min", 0) <= gc_pct <= m.get("max", 100):
            return config.get("marginal_multiplier", 0.7)
    return 0.0

def compute_block_multiplier(block_len, config):
    if block_len < config.get("penalty_min_len", 20):
        return 0.0 # Veto
    if block_len >= config.get("ideal_min_len", 35):
        return 1.0
    return config.get("penalty_multiplier", 0.8)

def compute_homo_multiplier(has_homo, config):
    return config.get("penalty_multiplier", 0.6) if has_homo else 1.0

# --- Step 3: Intersection & Composite Scoring ---

GRADE_TEXT = {
    "S": "Optimal Target — Composite score ≥0.95; no episodic selection detected; high confidence for universal primer design.",
    "A+": "High-Priority Target — Composite score ≥0.90; near-optimal conservation; negligible evolutionary risk.",
    "A": "Recommended Target — Composite score ≥0.85; strong conservation profile; minor multiplier penalties applied.",
    "A-": "Acceptable Target — Composite score ≥0.80; conserved with minor biomechanical caveats; suitable for validation.",
    "B+": "Moderate Priority — Composite score ≥0.75; satisfactory conservation; one or more moderate penalties detected.",
    "B": "Low Priority — Composite score ≥0.70; approaching lower tolerance bounds; requires careful thermodynamic validation.",
    "C": "At-Risk Target — Composite score ≥0.60; significant penalties (e.g., episodic selection, GC imbalance); not recommended without further review.",
    "D": "Excluded — Composite score <0.60; failed one or more critical biomechanical or evolutionary constraints.",
}

def compute_durability_multiplier(rel_rate, cfg):
    """Discount a region by how fast its nucleotides turn over.

    `rel_rate` is the region's mean codon substitution rate divided by the gene's
    median, so it is self-calibrating: 1.0 means "typical for this gene". None
    means the rate could not be measured (non-coding region, or no FUBAR data).

    This replaces an additive purifying-selection term in the core score. The
    additive form let evolutionary evidence SUBSTITUTE for conservation: over the
    full range of the old term it could overturn a 25.8 percentage-point gap in
    conservation, i.e. rank a site with ~14.5% minor-allele frequency above a
    perfectly conserved one. A primer cannot trade a mismatch today for evidence
    of constraint tomorrow -- it either hybridises or it does not. Conservation is
    therefore the core, and durability only modulates it, exactly as every other
    evolutionary and biomechanical factor in this module already does.

    The floor is the single free parameter and has a direct reading: with floor f,
    the multiplier's full range can overturn a conservation gap of at most (1 - f).
    At the default 0.95 that is 5 percentage points, about the difference between
    0% and 1.6% minor-allele frequency -- a trade the pipeline is willing to make,
    since a site at 1.5% MAF under low substitution pressure may well outlast one
    at 0% that is conserved only by chance in a small sample.
    """
    floor = cfg.get("floor", 0.95)
    no_data = cfg.get("no_evidence", 0.95)
    cap = cfg.get("rate_cap", 3.0)   # rel_rate at or above which the floor applies

    if rel_rate is None:
        return no_data
    if rel_rate <= 1.0:
        return 1.0
    if rel_rate >= cap:
        return floor
    # Linear between the gene median and the cap.
    frac = (rel_rate - 1.0) / (cap - 1.0)
    return 1.0 - frac * (1.0 - floor)


def score_regions(regions, profile, ref_seq, codon_data, gard_breakpoints, genes, utrs,
                  scoring_config, tiers_config, gard_buffer_bp,
                  meme_pvalue_threshold=MEME_PVALUE_THRESHOLD,
                  gene_median_rate=None):
    """Score each conserved region using Core x Multipliers."""
    
    # The core is conservation alone. Purifying selection used to be a second
    # additive term here (conservation_weight / purifying_weight); it is now the
    # `durability` multiplier below. See compute_durability_multiplier().
    gene_median_rate = gene_median_rate or {}
    
    # Multiplier Configs
    mult_cfgs = scoring_config.get("multipliers", {})
    meme_cfg = mult_cfgs.get("meme", {})
    rec_cfg = mult_cfgs.get("recombination", {})
    gc_cfg = mult_cfgs.get("gc_content", {})
    block_cfg = mult_cfgs.get("conserved_block", {})
    homo_cfg = mult_cfgs.get("homopolymer", {})
    dur_cfg = mult_cfgs.get("durability", {})

    scored = []

    for region in regions:
        start, end = region["start"], region["end"]
        # Always compute mid — used by both gene lookup and GARD check below
        mid = (start + end) // 2
        
        if region.get("gene"):
            gene_name = region["gene"]
            is_coding = True
        else:
            gene_name, is_coding = find_gene_for_position(mid, genes, utrs)

        # 1. Gather HyPhy Data
        fubar_probs, meme_pvals, subst_rates = [], [], []
        has_hyphy = False
        gene_codon_data = codon_data.get(gene_name, {})
        for pos in range(start, end):
            if pos in gene_codon_data:
                has_hyphy = True
                fubar_probs.append(gene_codon_data[pos]["fubar_purifying_prob"])
                meme_pvals.append(gene_codon_data[pos]["meme_pvalue"])
                r = gene_codon_data[pos].get("subst_rate")
                if r is not None:
                    subst_rates.append(r)

        avg_fubar = sum(fubar_probs) / len(fubar_probs) if fubar_probs else 0.0

        # Nucleotide turnover relative to this gene's own median, so the figure is
        # comparable across genes and builds despite FUBAR's per-alignment scaling.
        median_rate = gene_median_rate.get(gene_name) if gene_name else None
        avg_rate = (sum(subst_rates) / len(subst_rates)) if subst_rates else None
        if avg_rate is not None and median_rate:
            rel_rate = avg_rate / median_rate
        else:
            rel_rate = None
        # HyPhy recommends p <= 0.1 for MEME, deliberately looser than the usual
        # 0.05 because MEME is conservative on small alignments. This was hardcoded
        # to 0.05 while config carried an unread hyphy.p_value_threshold.
        meme_danger_sites = sum(1 for p in meme_pvals if p < meme_pvalue_threshold)

        # Check GARD
        near_recomb = False
        half_len = (end - start) // 2
        for bp in gard_breakpoints:
            if abs(bp["genome_pos"] - mid) <= gard_buffer_bp + half_len:
                near_recomb = True
                break

        # 2. Biomechanical features
        ref_subseq = extract_reference_subsequence(ref_seq, start, end)
        gc_pct = calculate_gc_content(ref_subseq)
        
        block_t_hr = block_cfg.get("entropy_threshold", 0.10)
        max_block_len = find_max_contiguous_block(profile, start, end, threshold=block_t_hr)
        
        has_homo = detect_homopolymer(ref_subseq, homo_cfg.get("max_repeats", 5))
        mismatches = get_mismatch_positions(profile, start, end, threshold=block_t_hr)

        # 3. Core Score Calculation
        # Normalized conservation over the 5-symbol alphabet {A,T,G,C,-}:
        # max entropy = log2(5); clamp guards against float drift / all-gap columns.
        cons_norm = 1.0 - (region["mean_entropy"] / MAX_ENTROPY_BITS)
        cons_norm = max(0.0, min(1.0, cons_norm))
        
        # The core is conservation alone: it is the only quantity that determines
        # whether the primer matches the strains circulating right now. Everything
        # else modulates it. The old code capped non-coding regions at 0.95 here
        # while coding regions carried an additive purifying term, which put the
        # two branches on incomparable scales -- a perfectly conserved coding
        # region with a weak selection signal scored 0.83 against a non-coding
        # region's 0.95. Both now go through the same durability multiplier.
        core = cons_norm

        # 4. Multipliers
        m_dur = compute_durability_multiplier(rel_rate, dur_cfg)
        m_meme = compute_meme_multiplier(meme_danger_sites, meme_cfg)
        m_rec = rec_cfg.get("near_penalty", 0.5) if near_recomb else 1.0
        m_gc = compute_gc_multiplier(gc_pct, gc_cfg)
        m_block = compute_block_multiplier(max_block_len, block_cfg)
        m_homo = compute_homo_multiplier(has_homo, homo_cfg)

        final_score = core * m_dur * m_meme * m_rec * m_gc * m_block * m_homo
        final_score = max(0.0, min(1.0, final_score))

        # --- Tier System (Fail-safe checks on raw metrics) ---
        strict = tiers_config.get("strict", {})
        moderate = tiers_config.get("moderate", {})
        cons_frac = cons_norm # Using normalized score directly
        tier = 3

        if has_hyphy:
            if (cons_frac >= strict.get("min_conservation", 0.95)
                    and avg_fubar >= strict.get("min_purifying_prob", 0.9)
                    and meme_danger_sites == 0):
                tier = 1
            elif (cons_frac >= moderate.get("min_conservation", 0.85)
                  and avg_fubar >= moderate.get("min_purifying_prob", 0.7)):
                tier = 2
        else:
            if cons_frac >= strict.get("min_conservation", 0.95):
                tier = 1
            elif cons_frac >= moderate.get("min_conservation", 0.85):
                tier = 2

        # --- Grade Classification ---
        # NOTE: two distinct classifications are emitted, previously easy to confuse
        # because the code's letters clashed with the documentation's:
        #   grade      — fine-grained ranking on the MULTIPLIED final_score (S..D).
        #                Code-internal; "A" here means final_score >= 0.85, which is
        #                NOT the docs' Grade A.
        #   tier_grade — the A/B/C/D scheme documented in
        #                docs/pipeline-architecture/m07-bridge-builder.md, keyed on
        #                the RAW conservation (C) and purifying (P) metrics via tier.
        #                This is the one to read as "Grade A / B are primary targets".
        if final_score >= 0.95: grade = "S"
        elif final_score >= 0.90: grade = "A+"
        elif final_score >= 0.85: grade = "A"
        elif final_score >= 0.80: grade = "A-"
        elif final_score >= 0.75: grade = "B+"
        elif final_score >= 0.70: grade = "B"
        elif final_score >= 0.60: grade = "C"
        else: grade = "D"

        # Documented tier -> grade mapping (Strict=A, Moderate=B, Survival=C, Reject=D).
        survival_min = tiers_config.get("survival", {}).get("min_conservation", 0.80)
        if tier == 1:
            tier_grade = "A"
        elif tier == 2:
            tier_grade = "B"
        elif cons_frac >= survival_min:
            tier_grade = "C"
        else:
            tier_grade = "D"

        # Build the three sequence representations for this window:
        #   consensus_sequence       — majority-vote per column (plain ATGC)
        #   iupac_consensus_sequence — IUPAC-encoded where minor allele >= threshold
        #   reference_sequence       — raw slice from reference genome
        consensus_seq = "".join(
            profile[i]["consensus_base"] for i in range(start, end)
        )
        iupac_consensus_seq = "".join(
            profile[i]["iupac_base"] for i in range(start, end)
        )
        # Count IUPAC-ambiguous positions for downstream awareness
        n_iupac_positions = sum(
            1 for i in range(start, end)
            if profile[i]["iupac_base"] not in "ACGT-"
        )
        iupac_density = round(n_iupac_positions / max(end - start, 1), 4)

        scored.append({
            "gene": gene_name,
            "start": start,
            "end": end,
            "length": end - start,
            "conservation_pct": round(region["conservation_pct"], 2),
            "mean_entropy": round(region["mean_entropy"], 6),
            "avg_fubar_purifying_prob": round(avg_fubar, 4),
            # Nucleotide turnover: what actually determines how long a primer keeps
            # matching. avg_fubar measures protein-level constraint, which is not
            # the same thing — see parse_fubar_json().
            "avg_subst_rate": round(avg_rate, 5) if avg_rate is not None else "",
            "rel_subst_rate": round(rel_rate, 3) if rel_rate is not None else "",
            "durability_multiplier": round(m_dur, 3),
            "meme_episodic_sites": meme_danger_sites,
            "gc_pct": round(gc_pct, 2),
            "max_block_len": max_block_len,
            "has_homopolymer": has_homo,
            "composite_score": round(final_score, 4),
            "tier": tier,
            "grade": grade,
            "tier_grade": tier_grade,
            "recommendation": GRADE_TEXT[grade],
            # Three sequence representations
            "consensus_sequence": consensus_seq,
            "iupac_consensus_sequence": iupac_consensus_seq,
            "iupac_positions": n_iupac_positions,
            "iupac_density": iupac_density,
            "reference_sequence": ref_subseq,
            "flags": {
                "coding": is_coding,
                "has_hyphy_data": has_hyphy,
                "near_recombination": near_recomb,
            },
            "multipliers": {
                "meme": m_meme,
                "recombination": m_rec,
                "gc_content": m_gc,
                "conserved_block": m_block,
                "homopolymer": m_homo
            },
            "mismatch_positions": mismatches
        })

    scored.sort(key=lambda r: r["composite_score"], reverse=True)
    for i, r in enumerate(scored):
        r["rank"] = i + 1

    return scored

def apply_nms(scored_regions, overlap_threshold=0.5):
    scored_regions.sort(key=lambda x: x["composite_score"], reverse=True)
    kept = []
    
    for current in scored_regions:
        is_suppressed = False
        for winner in kept:
            if current["gene"] == winner["gene"]:
                overlap_start = max(current["start"], winner["start"])
                overlap_end = min(current["end"], winner["end"])
                if overlap_end > overlap_start:
                    overlap_len = overlap_end - overlap_start
                    overlap_ratio = overlap_len / current["length"]
                    if overlap_ratio > overlap_threshold:
                        is_suppressed = True
                        break
        if not is_suppressed:
            kept.append(current)
            
    for idx, r in enumerate(kept):
        r["rank"] = idx + 1
        
    return kept

def write_bed(regions, ref_id, path):
    with open(path, "w") as f:
        for r in regions:
            rid = f"{r['gene']}_region_{r['rank']}"
            f.write(f"{ref_id}\t{r['start']}\t{r['end']}\t{rid}\t"
                    f"{r['composite_score']}\t+\n")
    print(f"[M07] Wrote BED: {path}")

def write_report(regions, meta, path):
    report = {
        "pipeline_version": "2.0",
        "scoring_method": "core_x_multiplier",
        "coordinate_system": "0-based half-open",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "track_name": meta.get("track_name", "Universal"),
        "msa_sequences": meta.get("n_sequences", 0),
        "msa_length": meta.get("aln_length", 0),
        "genes_analyzed": meta.get("genes_analyzed", []),
        "hyphy_genes": meta.get("hyphy_genes", []),
        "total_conserved_regions": len(regions),
        "grade_summary": {g: sum(1 for r in regions if r["grade"] == g)
                          for g in ["S", "A+", "A", "A-", "B+", "B", "C", "D"]},
        "tier_summary": {f"tier_{t}": sum(1 for r in regions if r["tier"] == t)
                         for t in (1, 2, 3)},
        "regions": regions,
        "gard_breakpoints": meta.get("gard_breakpoints", []),
        "warnings": meta.get("warnings", []),
    }
    with open(path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"[M07] Wrote report: {path}")

def write_manifest(regions, build_name, ref_genome, buffer_bp, path, profile, codon_data, design_regions=None):
    """Write the YAML manifest. If design_regions is set (e.g. ['segment-MP']),
    only export regions whose gene label appears in that list."""
    if design_regions:
        filtered = [r for r in regions if r.get("gene", "") in design_regions]
        if not filtered:
            print(f"[M07] WARNING: design_regions filter {design_regions} matched no regions. "
                  f"Exporting all {len(regions)} regions instead.")
            filtered = regions
        else:
            print(f"[M07] design_regions filter: {len(regions)} → {len(filtered)} regions "
                  f"(kept: {design_regions})")
        regions = filtered

    manifest = {
        "pipeline_run_id": f"{build_name}_bridge",
        "source": "pandemic-preparedness-m07",
        "pipeline_version": "2.0",
        "coordinate_system": "0-based half-open",
        "reference_genome": str(ref_genome),
        "buffer_bp": buffer_bp,
        "total_candidates": len(regions),
        "candidate_regions": [],
    }

    for r in regions:
        inc_start = max(0, r["start"] - buffer_bp)
        inc_end = r["end"] + buffer_bp

        manifest["candidate_regions"].append({
            "region_id": f"{sanitize_region_id(r['gene'])}_region_{r['rank']}",
            "rank": r["rank"],
            "gene": r["gene"],
            "tier": r["tier"],
            "grade": r["grade"],
            "tier_grade": r.get("tier_grade", ""),
            "composite_score": r["composite_score"],
            "conserved_zone": {
                "start": r["start"],
                "end": r["end"],
                "length": r["length"],
            },
            "included_region": {
                "start": inc_start,
                "end": inc_end,
                "length": inc_end - inc_start,
            },
            "biomechanics": {
                "conservation_pct": r["conservation_pct"],
                "mean_entropy": r["mean_entropy"],
                "avg_fubar_purifying_prob": r["avg_fubar_purifying_prob"],
                "rel_subst_rate": r.get("rel_subst_rate", ""),
                "durability_multiplier": r.get("durability_multiplier", ""),
                "meme_episodic_sites": r["meme_episodic_sites"],
                "gc_pct": r["gc_pct"],
                "max_block_len": r["max_block_len"],
                "has_homopolymer": r["has_homopolymer"],
                "mismatch_positions": r["mismatch_positions"],
                "multipliers": r["multipliers"],
            },
            # --- Three sequence representations ---
            # consensus_sequence       : majority-vote per MSA column (plain ATGC)
            # iupac_consensus_sequence : IUPAC-encoded where minor allele >= 20%
            #                           Use this for: understanding variability,
            #                           choosing excluded_regions in primer design.
            #                           Do NOT use directly as RAA primer template.
            # reference_sequence       : raw slice from reference genome (NC_*)
            "sequences": {
                "consensus_sequence":        r["consensus_sequence"],
                "iupac_consensus_sequence":  r["iupac_consensus_sequence"],
                "reference_sequence":        r["reference_sequence"],
                "iupac_positions":           r["iupac_positions"],
                "iupac_density":             r["iupac_density"],
                "iupac_threshold":           IUPAC_THRESHOLD,
            },
            "positional_scores": {
                "entropy": [round(profile[i]["entropy"], 4) for i in range(r["start"], r["end"])],
                "fubar_prob": [round(codon_data.get(r["gene"], {}).get(i, {}).get("fubar_purifying_prob", 0.5), 4) for i in range(r["start"], r["end"])],
                "meme_pval": [round(codon_data.get(r["gene"], {}).get(i, {}).get("meme_pvalue", 1.0), 4) for i in range(r["start"], r["end"])],
            },
            "flags": r["flags"],
            "recommendation": r["recommendation"],
        })

    with open(path, "w") as f:
        yaml.dump(manifest, f, default_flow_style=False,
                  sort_keys=False, allow_unicode=True)
    print(f"[M07] Wrote manifest: {path}")


def write_csv(regions, path):
    """Export a flat CSV table of scored regions for spreadsheet-based review (PI).

    Each row corresponds to one candidate primer target region. Nested fields
    (multipliers, flags) are flattened into individual columns.
    """
    fieldnames = [
        "rank",
        "region_id",
        "gene",
        "tier",
        "grade",
        "tier_grade",
        "composite_score",
        "start_0based",
        "end_0based",
        "length_bp",
        "conservation_pct",
        "mean_entropy",
        "avg_fubar_purifying_prob",
        "avg_subst_rate",
        "rel_subst_rate",
        "durability_multiplier",
        "meme_episodic_sites",
        "iupac_positions",
        "iupac_density",
        "gc_pct",
        "max_conserved_block_len",
        "has_homopolymer",
        "mult_meme",
        "mult_recombination",
        "mult_gc_content",
        "mult_conserved_block",
        "mult_homopolymer",
        "flag_coding_region",
        "flag_has_hyphy_data",
        "flag_near_recombination",
        "consensus_sequence",
        "iupac_consensus_sequence",
        "reference_sequence",
        "recommendation",
    ]

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in regions:
            writer.writerow({
                "rank": r.get("rank", ""),
                "region_id": f"{sanitize_region_id(r['gene'])}_region_{r.get('rank', '')}",
                "gene": r.get("gene", ""),
                "tier": r.get("tier", ""),
                "grade": r.get("grade", ""),
                "tier_grade": r.get("tier_grade", ""),
                "composite_score": r.get("composite_score", ""),
                "start_0based": r.get("start", ""),
                "end_0based": r.get("end", ""),
                "length_bp": r.get("length", ""),
                "conservation_pct": r.get("conservation_pct", ""),
                "mean_entropy": r.get("mean_entropy", ""),
                "avg_fubar_purifying_prob": r.get("avg_fubar_purifying_prob", ""),
                "avg_subst_rate": r.get("avg_subst_rate", ""),
                "rel_subst_rate": r.get("rel_subst_rate", ""),
                "durability_multiplier": r.get("durability_multiplier", ""),
                "meme_episodic_sites": r.get("meme_episodic_sites", ""),
                "iupac_positions": r.get("iupac_positions", ""),
                "iupac_density": r.get("iupac_density", ""),
                "gc_pct": r.get("gc_pct", ""),
                "max_conserved_block_len": r.get("max_block_len", ""),
                "has_homopolymer": r.get("has_homopolymer", ""),
                "mult_meme": r.get("multipliers", {}).get("meme", ""),
                "mult_recombination": r.get("multipliers", {}).get("recombination", ""),
                "mult_gc_content": r.get("multipliers", {}).get("gc_content", ""),
                "mult_conserved_block": r.get("multipliers", {}).get("conserved_block", ""),
                "mult_homopolymer": r.get("multipliers", {}).get("homopolymer", ""),
                "flag_coding_region": r.get("flags", {}).get("coding", ""),
                "flag_has_hyphy_data": r.get("flags", {}).get("has_hyphy_data", ""),
                "flag_near_recombination": r.get("flags", {}).get("near_recombination", ""),
                "consensus_sequence": r.get("consensus_sequence", ""),
                "iupac_consensus_sequence": r.get("iupac_consensus_sequence", ""),
                "reference_sequence": r.get("reference_sequence", ""),
                "recommendation": r.get("recommendation", ""),
            })
    print(f"[M07] Wrote CSV summary: {path}")

def execute_track(track_name, records, genes, utrs, target_genes, hyphy_dir, 
                  window_size, entropy_threshold, overlap_threshold, gard_buffer_bp, 
                  scoring_config, tiers_config, ref_id, ref_seq, outdir, build_name, buffer_bp, args_reference,
                  strategy="windowed", design_regions=None, prefix="",
                  meme_pvalue_threshold=MEME_PVALUE_THRESHOLD):
    print(f"\n[M07] ============================================")
    print(f"[M07]  Track: {track_name}")
    print(f"[M07] ============================================")

    print(f"\n[M07] [{track_name}] Step 1: Generating Entropy Profile")
    profile, n_seq, aln_len = build_entropy_profile(records)
    
    sanitized_name = track_name.replace('/', '_').replace(' ', '_')
    prefix = "" if track_name == "Universal" else f"subtype_{sanitized_name}_"
    
    export_entropy_profile(profile, outdir / f"{prefix}entropy_profile.tsv")

    print(f"\n[M07] [{track_name}] Step 2: Conservation Scanner (Strategy: {strategy})")
    if strategy == "full_segment":
        # Strategy: Full Segment — Create one region per gene/UTR (preserving gene labels
        # so codon_data lookups in score_regions work correctly for FUBAR/MEME data).
        regions = []

        # Build gene-level regions from GFF annotations
        all_features = []
        for g in genes:
            all_features.append({
                "start": g["start"],
                "end": min(g["end"], aln_len),
                "gene": g["name"],
            })
        for u in utrs:
            all_features.append({
                "start": u["start"],
                "end": min(u["end"], aln_len),
                "gene": u["name"],
            })

        if all_features:
            all_features.sort(key=lambda x: x["start"])
            for feat in all_features:
                s, e = feat["start"], feat["end"]
                if e <= s:
                    continue
                pos_entropies = [profile[i]["entropy"] for i in range(s, e) if i < len(profile)]
                mean_ent = sum(pos_entropies) / max(1, len(pos_entropies))
                regions.append({
                    "start": s,
                    "end": e,
                    "mean_entropy": mean_ent,
                    "conservation_pct": round((1 - mean_ent / MAX_ENTROPY_BITS) * 100, 2),
                    "length": e - s,
                    "gene": feat["gene"],
                })
            print(f"[M07] [{track_name}] Full-segment: created {len(regions)} gene-level regions "
                  f"({', '.join(r['gene'] for r in regions)})")
        else:
            # Fallback: no GFF features found → single holistic region
            mean_ent = sum(p["entropy"] for p in profile) / max(1, aln_len)
            regions = [{
                "start": 0,
                "end": aln_len,
                "mean_entropy": mean_ent,
                "conservation_pct": round((1 - mean_ent / MAX_ENTROPY_BITS) * 100, 2),
                "length": aln_len,
                "gene": "Whole_Segment",
            }]
            print(f"[M07] [{track_name}] WARNING: No GFF genes found — fell back to single "
                  f"Whole_Segment region ({aln_len} bp)")
    else:
        # Strategy: Windowed (Original)
        regions = scan_conservation(profile, window_size, entropy_threshold)

        if not regions:
            print(f"[M07] [{track_name}] WARNING: Relaxing threshold significantly to find candidates (survival mode)")
            # 0.45 was calibrated on the old 4-symbol entropy scale; rescale by
            # log2(5)/log2(4) so the survival fallback keeps the same strictness.
            regions = scan_conservation(profile, window_size, 0.45 * (MAX_ENTROPY_BITS / 2.0))

        if not regions:
             print(f"[M07] [{track_name}] No conserved regions found even after max relaxation. Skipping.")
             return

    print(f"\n[M07] [{track_name}] Step 3: HyPhy Score Mapper")
    codon_data, gard_bps, gene_median_rate = map_hyphy_to_genome(genes, hyphy_dir, target_genes)
    total_positions = sum(len(codon_data[g]) for g in codon_data)
    print(f"[M07] [{track_name}] Mapped {total_positions} positions across {len(codon_data)} genes, {len(gard_bps)} GARD breakpoints")

    print(f"\n[M07] [{track_name}] Step 4: Biomechanical Analysis & Core × Multipliers Scoring")
    scored = score_regions(regions, profile, ref_seq, codon_data, gard_bps, genes, utrs,
                           scoring_config, tiers_config, gard_buffer_bp,
                           meme_pvalue_threshold=meme_pvalue_threshold,
                           gene_median_rate=gene_median_rate)

    print(f"\n[M07] [{track_name}] Step 5: Non-Maximum Suppression (NMS)")
    final_scored = apply_nms(scored, overlap_threshold=overlap_threshold)

    grade_keys = ["S", "A+", "A", "A-", "B+", "B", "C", "D"]
    grades = {g: sum(1 for r in final_scored if r["grade"] == g) for g in grade_keys}
    grade_str = " ".join([f"{g}={grades[g]}" for g in grade_keys])
    print(f"[M07] [{track_name}] Scored {len(final_scored)} disjoint/optimal regions: {grade_str}")

    warnings = []
    if not codon_data:
        warnings.append("No HyPhy data loaded. All scores are conservation-only.")

    meta = {
        "track_name": track_name,
        "n_sequences": n_seq,
        "aln_length": aln_len,
        "genes_analyzed": [g["name"] for g in genes] + [u["name"] for u in utrs],
        "hyphy_genes": [g for g in target_genes if (hyphy_dir / f"fubar_{g}.json").exists()],
        "gard_breakpoints": [{"genome_pos": bp["genome_pos"], "gene": bp["gene"]} for bp in gard_bps],
        "warnings": warnings,
    }

    print(f"\n[M07] [{track_name}] Step 6: Output Generation")
    write_bed(final_scored, ref_id, outdir / f"{prefix}conserved_regions.bed")
    write_report(final_scored, meta, outdir / f"{prefix}bridge_report.json")
    write_manifest(final_scored, build_name, args_reference, buffer_bp,
                   outdir / f"{prefix}primer_targets_manifest.yaml", profile, codon_data,
                   design_regions=design_regions)
    write_csv(final_scored, outdir / f"{prefix}candidate_targets_summary.csv")

    print(f"\n[M07] [{track_name}] Complete — {len(final_scored)} candidate regions")
    if final_scored:
        top = final_scored[0]
        print(f"[M07] [{track_name}] #1: {top['gene']} ({top['start']}-{top['end']}) Grade {top['grade']}, Score {top['composite_score']}")


def main():
    parser = argparse.ArgumentParser(
        description="M07 Bridge Builder v2.0 — Core x Multipliers Integrator")
    parser.add_argument("--msa", required=True,
                        help="Aligned FASTA (MSA from augur align / M03)")
    parser.add_argument("--gff", required=True,
                        help="Gene annotation GFF3 file")
    parser.add_argument("--hyphy-dir", required=True,
                        help="Directory containing fubar_*.json, meme_*.json, gard_*.json")
    parser.add_argument("--reference", required=True,
                        help="Reference genome FASTA")
    parser.add_argument("--metadata", required=False,
                        help="Metadata TSV file for subtype/clade splitting (Route B)")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for BED, JSON, YAML")
    parser.add_argument("--config", required=True,
                        help="Pipeline config.yaml path")
    parser.add_argument("--build-name", default="unknown",
                        help="Build name for manifest identification")
    parser.add_argument("--strategy", default="windowed", choices=["windowed", "full_segment"],
                        help="Analysis strategy: windowed (slices) or full_segment (global)")
    parser.add_argument("--min-subtype-count", type=int, default=30,
                        help="Minimum sequences required to analyze a subtype (default: 30)")
    parser.add_argument("--iupac-threshold", type=float, default=0.15, help="IUPAC minor frequency threshold")
    args = parser.parse_args()

    # Load config
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    bridge = config.get("bridge", {})
    window_size = bridge.get("window_size", 150)
    entropy_threshold = bridge.get("entropy_threshold", 0.015)
    overlap_threshold = bridge.get("overlap_threshold", 0.75)
    buffer_bp = bridge.get("buffer_bp", 100)
    gard_buffer_bp = bridge.get("gard_buffer_bp", 50)
    
    scoring_config = bridge.get("scoring", {})
    if not scoring_config and "weights" in bridge:
        print("[M07] WARNING: Falling back to v1.0 weights config. Please update to 'scoring' config block.", file=sys.stderr)
        # Rudimentary fallback mapping. The composite score S = w_c*C + w_p*P assumes
        # w_c + w_p == 1.0; the legacy defaults (0.5 / 0.25) sum to 0.75, which
        # deflates EVERY score by 25% and makes the documented 0.95/0.85 tier
        # thresholds unreachable even for a perfect region. Renormalise so the
        # legacy ratio is preserved but the invariant holds.
        _w_c = bridge["weights"].get("conservation", 0.5)
        _w_p = bridge["weights"].get("purifying", 0.25)
        _total = _w_c + _w_p
        if _total <= 0:
            print("[M07] ERROR: legacy bridge.weights sum to zero.", file=sys.stderr)
            sys.exit(1)
        if abs(_total - 1.0) > 1e-9:
            print(f"[M07] WARNING: legacy weights sum to {_total:.3f}, not 1.0 — "
                  f"renormalising to conservation={_w_c/_total:.3f}, "
                  f"purifying={_w_p/_total:.3f}.", file=sys.stderr)
        scoring_config = {
            "core": {
                "conservation_weight": _w_c / _total,
                "purifying_weight": _w_p / _total
            },
            "multipliers": {}
        }
        
    tiers_config = bridge.get("tiers", {})
    target_genes = config.get("hyphy", {}).get("target_genes", [])
    # hyphy.p_value_threshold existed in every config but was read by nothing;
    # M07 counted MEME danger sites against a hardcoded 0.05 instead. Wire it, and
    # keep HyPhy's own recommendation (0.1) as the fallback.
    meme_p = config.get("hyphy", {}).get("p_value_threshold", MEME_PVALUE_THRESHOLD)
    print(f"[M07] MEME episodic-selection threshold: p < {meme_p}")
    design_regions = bridge.get("design_regions", None)   # e.g. ["segment-MP"] or null
    if design_regions:
        print(f"[M07] design_regions filter active: {design_regions}")
    hyphy_dir = Path(args.hyphy_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    global IUPAC_THRESHOLD
    IUPAC_THRESHOLD = args.iupac_threshold

    all_records = list(SeqIO.parse(args.msa, "fasta"))
    if not all_records:
        print("ERROR: MSA file is empty", file=sys.stderr)
        sys.exit(1)
        
    ref_records = list(SeqIO.parse(args.reference, "fasta"))
    if not ref_records:
         print("ERROR: Reference file is empty", file=sys.stderr)
         sys.exit(1)
         
    ref_seq = str(ref_records[0].seq).upper()
    ref_id = ref_records[0].id

    # COORDINATE INVARIANT: this module treats an alignment column index as a
    # reference genome index (extract_reference_subsequence, HyPhy codon mapping,
    # conserved_regions.bed, the manifest that M08/M09/IGV consume as ground truth).
    # That only holds because M03 runs `augur align --reference-sequence ... --fill-gaps`,
    # which pins the alignment length to the reference via mafft --keeplength.
    # The invariant is enforced entirely OUTSIDE this script, so verify it here:
    # if a different reference is ever passed, ref_seq[start:end] slices the wrong
    # bases and Python raises nothing — every downstream coordinate is silently wrong.
    aln_len_check = len(all_records[0].seq)
    if len(ref_seq) != aln_len_check:
        print(
            f"ERROR: reference length ({len(ref_seq)} bp, {ref_id}) != alignment length "
            f"({aln_len_check} bp).\n"
            f"       M07 maps alignment columns directly onto reference coordinates, so "
            f"these MUST match.\n"
            f"       Likely cause: the --reference given to M07 is not the one used for "
            f"`augur align` in M03.\n"
            f"       Refusing to emit coordinates that would be silently wrong.",
            file=sys.stderr,
        )
        sys.exit(1)

    print("[M07] ============================================")
    print("[M07]  Bridge Builder v2.0 — Biomechanical Scoring")
    print("[M07] ============================================")
    print(f"[M07] Reference : {ref_id} ({len(ref_seq)} bp, matches alignment length)")
    print(f"[M07] Strategy  : {args.strategy}")
    if args.strategy == "windowed":
        print(f"[M07] Window    : {window_size} bp")
    print(f"[M07] Threshold : {entropy_threshold}")

    # ----- Parse GFF -----
    genes, utrs = parse_gff(args.gff)
    print(f"[M07] GFF: {len(genes)} genes, {len(utrs)} UTRs")
    if "ALL" in target_genes:
        target_genes = [g["name"] for g in genes]

    # Execute Universal Track
    execute_track("Universal", all_records, genes, utrs, target_genes, hyphy_dir,
                  window_size, entropy_threshold, overlap_threshold, gard_buffer_bp,
                  scoring_config, tiers_config, ref_id, ref_seq, outdir, args.build_name, buffer_bp, args.reference,
                  strategy=args.strategy, design_regions=design_regions,
                  meme_pvalue_threshold=meme_p)

    # Route B: Subtype-Specific Tracks
    if args.metadata and Path(args.metadata).exists():
        print(f"\n[M07] ============================================")
        print(f"[M07]  ROUTE B: Subtype Auto-Detection")
        print(f"[M07] ============================================")
        strain_to_subtype = {}
        with open(args.metadata, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                strain = row.get("strain")
                subtype = row.get("subtype")
                if strain and subtype and subtype != "unknown" and subtype != "?":
                    strain_to_subtype[strain] = subtype
                    
        from collections import defaultdict
        subtype_groups = defaultdict(list)
        for rec in all_records:
             rec_id = rec.id.split('|')[0]
             subtype = strain_to_subtype.get(rec_id)
             if subtype:
                 subtype_groups[subtype].append(rec)
                 
        for subtype, recs in subtype_groups.items():
            if len(recs) >= args.min_subtype_count:
                execute_track(subtype, recs, genes, utrs, target_genes, hyphy_dir,
                              window_size, entropy_threshold, overlap_threshold, gard_buffer_bp,
                              scoring_config, tiers_config, ref_id, ref_seq, outdir, args.build_name, buffer_bp, args.reference,
                              strategy=args.strategy, design_regions=design_regions,
                              meme_pvalue_threshold=meme_p)
            else:
                print(f"[M07] Skipping subtype {subtype} (only {len(recs)} sequences, needs >={args.min_subtype_count})")


if __name__ == "__main__":
    main()
