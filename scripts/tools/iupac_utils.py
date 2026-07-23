#!/usr/bin/env python3
"""
Shared IUPAC Degenerate Base Utilities
=======================================
Centralised expansion and thermodynamic utility functions shared by all
Module 11 screening scripts (m11_multiplex_checker.py, m11_triplex_checker.py).

Extracted from duplicated logic so that:
  - The IUPAC cap (max_variants) is configurable from CLI / caller.
  - A single implementation is tested and maintained.
  - Downstream scripts simply import from here.

Usage
-----
from scripts.tools.iupac_utils import expand_iupac, calculate_gc, calculate_tm

Or from the same directory:
  import sys; sys.path.insert(0, str(Path(__file__).parent.parent.parent))
  from scripts.tools.iupac_utils import expand_iupac
"""

import itertools
from functools import lru_cache
from typing import List, Tuple

# ---------------------------------------------------------------------------
# IUPAC Base Expansion Table
# ---------------------------------------------------------------------------
IUPAC_MAP = {
    'A': ('A',), 'C': ('C',), 'G': ('G',), 'T': ('T',),
    'R': ('A', 'G'),  # puRine
    'Y': ('C', 'T'),  # pYrimidine
    'S': ('G', 'C'),  # Strong
    'W': ('A', 'T'),  # Weak
    'K': ('G', 'T'),  # Keto
    'M': ('A', 'C'),  # aMino
    'B': ('C', 'G', 'T'),  # not A
    'D': ('A', 'G', 'T'),  # not C
    'H': ('A', 'C', 'T'),  # not G
    'V': ('A', 'C', 'G'),  # not T
    'N': ('A', 'C', 'G', 'T'),  # aNy
}

# Default cap — matches config.yaml degenerate.max_iupac_variants = 16
DEFAULT_MAX_VARIANTS = 16


def expand_iupac(seq: str, max_variants: int = DEFAULT_MAX_VARIANTS) -> List[str]:
    """Expand a degenerate IUPAC sequence into all canonical combinations.

    Parameters
    ----------
    seq : str
        IUPAC-encoded sequence (e.g. "CTGACTAARGGR").
    max_variants : int
        Maximum number of canonical variants to return. Variants beyond this
        cap are silently truncated. Default: 16 (matches config.yaml setting).
        Set to 0 or None to disable the cap (use with caution on long seqs).

    Returns
    -------
    List[str]
        List of fully canonical (A/C/G/T-only) sequences.
    """
    seq = seq.upper()
    options = [IUPAC_MAP.get(base, (base,)) for base in seq]
    # Use itertools.islice for memory efficiency when cap is set
    if max_variants:
        combos = []
        for combo in itertools.product(*options):
            combos.append(''.join(combo))
            if len(combos) >= max_variants:
                break
        return combos
    return [''.join(p) for p in itertools.product(*options)]


def count_iupac_variants(seq: str) -> int:
    """Return the total number of canonical variants (before capping).

    Parameters
    ----------
    seq : str
        IUPAC-encoded sequence.

    Returns
    -------
    int
        Product of the number of bases at each position.
    """
    seq = seq.upper()
    total = 1
    for base in seq:
        total *= len(IUPAC_MAP.get(base, (base,)))
    return total


def calculate_gc(seq: str) -> float:
    """Calculate expected GC content percentage from an IUPAC sequence.

    Every degenerate code contributes the FRACTION of its possible bases that are
    G or C, which is the expected GC content over the primer pool it represents.
    Counting only G/C/S (the previous behaviour) systematically under-estimates a
    degenerate oligo: B (C/G/T) is 2/3 GC and M (A/C) is 1/2 GC, but both scored 0.
    That can push a borderline design across a GC cutoff without anyone noticing.

    Parameters
    ----------
    seq : str
        IUPAC or canonical sequence.

    Returns
    -------
    float
        Expected GC percentage (0.0 – 100.0), rounded to 1 decimal place.
        Gaps and unknown characters are excluded from the denominator.
    """
    seq = (seq or "").upper()
    gc = 0.0
    counted = 0
    for base in seq:
        bases = IUPAC_MAP.get(base)
        if not bases:
            continue  # gaps / unknown symbols are not part of the oligo
        counted += 1
        gc += sum(1 for b in bases if b in ("G", "C")) / len(bases)
    return round((gc / counted) * 100.0, 1) if counted else 0.0


def calculate_tm(seq: str, max_variants: int = DEFAULT_MAX_VARIANTS) -> float:
    """Calculate the mean melting temperature across all IUPAC variants.

    Uses primer3.calc_tm on each canonical expansion and averages the result.

    Parameters
    ----------
    seq : str
        IUPAC-encoded primer or probe sequence.
    max_variants : int
        Cap for expand_iupac; also limits the Tm averaging set.

    Returns
    -------
    float
        Mean Tm in °C.
    """
    try:
        import primer3
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "primer3-py is required for Tm calculations. "
            "Install via: pip install primer3-py"
        ) from exc
    variants = expand_iupac(seq, max_variants=max_variants)
    tms = [primer3.calc_tm(v) for v in variants]
    return sum(tms) / len(tms)


# ---------------------------------------------------------------------------
# Heterodimerisation with Dual-Level Caching
# ---------------------------------------------------------------------------

# Module-level shared caches — populated across successive calls.
# Callers sharing the same process automatically reuse cached results.
_heterodimer_cache: dict = {}
_raw_p3_cache: dict = {}


def evaluate_pair_heterodimer(
    seq1: str,
    seq2: str,
    temp_c: float,
    max_variants: int = DEFAULT_MAX_VARIANTS,
) -> Tuple[float, bool]:
    """Calculate worst-case heterodimerisation ΔG between two IUPAC sequences.

    Expands both sequences into all canonical variants and evaluates every
    pairwise combination via primer3.calc_heterodimer, returning the *minimum*
    ΔG (most stable / most dangerous dimer).

    Parameters
    ----------
    seq1, seq2 : str
        IUPAC-encoded sequences.
    temp_c : float
        Reaction temperature in °C (e.g. 39.0 for RAA, 55.0 for qPCR).
    max_variants : int
        IUPAC expansion cap applied to both sequences.

    Returns
    -------
    (min_dg_kcal, structure_found) : Tuple[float, bool]
        min_dg_kcal — most negative ΔG found (kcal/mol); 0.0 if no structure.
        structure_found — True if at least one variant pair formed a dimer.
    """
    try:
        import primer3
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "primer3-py is required. Install via: pip install primer3-py"
        ) from exc

    key = (seq1, seq2, temp_c)
    if key in _heterodimer_cache:
        return _heterodimer_cache[key]
    key_rev = (seq2, seq1, temp_c)
    if key_rev in _heterodimer_cache:
        return _heterodimer_cache[key_rev]

    variants1 = expand_iupac(seq1, max_variants=max_variants)
    variants2 = expand_iupac(seq2, max_variants=max_variants)

    min_dg = 99.0
    struct_found = False

    for v1 in variants1:
        for v2 in variants2:
            p3_key = (v1, v2, temp_c)
            if p3_key in _raw_p3_cache:
                dg_kcal, sf = _raw_p3_cache[p3_key]
            else:
                res = primer3.calc_heterodimer(v1, v2, temp_c=temp_c)
                dg_kcal = res.dg / 1000.0
                sf = res.structure_found
                _raw_p3_cache[p3_key] = (dg_kcal, sf)
                _raw_p3_cache[(v2, v1, temp_c)] = (dg_kcal, sf)

            if sf:
                struct_found = True
                if dg_kcal < min_dg:
                    min_dg = dg_kcal
            else:
                if 0.0 < min_dg:
                    min_dg = 0.0

    if min_dg == 99.0:
        min_dg = 0.0

    result = (min_dg, struct_found)
    _heterodimer_cache[key] = result
    return result


def clear_caches() -> None:
    """Clear all module-level primer3 result caches.

    Call between independent screening runs if memory is a concern.
    """
    _heterodimer_cache.clear()
    _raw_p3_cache.clear()


# ---------------------------------------------------------------------------
# Multiplex Compatibility Thresholds
# ---------------------------------------------------------------------------
# Single source of truth for pairwise AND triplex screening. These lived as
# separate literals in each checker and had already drifted: the pairwise ΔTm
# limit was 3.5 while the triplex and the documentation used 4.0, so an
# identical pair was PASS in one report and WARNING in the other.
#
# ΔTm tolerance is a soft criterion for the 39°C RAA mode (isothermal strand
# invasion is enzyme-driven, not annealing-driven); it matters more for the 55°C
# qPCR benchmark. 4.0°C matches docs/pipeline-architecture/m11-multiplex-compatibility.md.
MAX_DELTA_TM_C = 4.0

# Cross-dimer ΔG limits, evaluated at the RAA reaction temperature of 39°C.
#
# There is NO validated RPA/RAA-specific ΔG threshold in the literature — the
# TwistAmp Assay Design Manual §2.2 states outright that "there are no fixed
# rules to predict how well a given amplification primer will work" and
# prescribes empirical screening of 8–10 primers per direction instead. These
# values therefore derive from the general IDT OligoAnalyzer guidance (self- and
# hetero-dimer ΔG should stay more positive than −9 kcal/mol).
#
# Read them at 39°C, not with PCR intuition. ΔG = ΔH − TΔS, so the same duplex is
# markedly more stable at the RAA temperature than at a qPCR temperature. Measured
# over 200 random 32-mer pairs with no designed complementarity:
#
#     39°C: median −4.12, 5th pct −6.78, 0.5% below −9, 88% below −3
#     60°C: median −1.33, 5th pct −3.07, 0.0% below −9,  6.5% below −3
#
# So a "−3 kcal/mol" cutoff borrowed from 60°C practice would reject 88% of
# entirely clean pairs here, while −9 at 39°C corresponds to roughly −4.7 at 60°C
# — already stricter than the general IDT guidance it came from.
MIN_CROSS_DIMER_DG = -9.0    # warn below this
VETO_CROSS_DIMER_DG = -12.0  # hard veto below this

# 3'-end complementarity. Only a 3' terminus that anneals can be extended by the
# polymerase, and it is those extension products that seed "primer noise" — which
# the TwistAmp manual §2.6 identifies as the factor that limits the sensitivity of
# a given primer pair. A global ΔG can look acceptable while the 3' ends are the
# part actually interacting, so this is screened separately.
#
# The criterion is the LENGTH of the longest perfectly complementary run anchored
# at a 3' terminus, not a ΔG. A ΔG computed over a short terminal window was tried
# first and rejected: a 5 nt window saturates near −6.3 kcal/mol regardless of how
# much complementarity exists, and widening it to 8 nt still overlaps the clean
# background badly (400 random 32-mer pairs reached −7.57 at worst, while a
# deliberate 8 bp 3' complement sat at −6.65 on average — the metric cannot
# separate them).
#
# The run length does separate, measured the same way:
#
#     background (400 random pairs)   >=4 bp: 16.8%   >=5 bp: 4.2%   >=6 bp: 0.5%
#     designed 5 bp complementarity   detected at >=5 bp in 100% of cases
#     designed 6/8 bp                 detected in 100% of cases
#
# A multiplex evaluates many oligo pairs at once (9 pairwise, 27 triplex), so a
# per-pair false-positive rate compounds. Tolerating up to 5 bp keeps the
# whole-combination false-warning rate near 4% pairwise / 13% triplex while still
# catching every deliberate 6 bp anchor.
MAX_THREE_PRIME_ANCHOR_BP = 5

# Terminal window used for the *reported* 3'-end ΔG. This is context for a reader
# comparing two flagged pairs; it does not gate anything, for the reasons above.
THREE_PRIME_WINDOW_NT = 8

# Compatibility score weights. Both checkers compute
#     score = 100 + min_dg * SCORE_WEIGHT_CROSS_DIMER - delta_tm * SCORE_WEIGHT_DELTA_TM
# and previously used different weights (pairwise 3.5/4.0, triplex 4.0/3.0), so
# the two reports' scores were not comparable despite sharing a 0–100 scale.
#
# The weighting favours ΔG over ΔTm because that is the direction the reaction
# chemistry points: at 39°C strand invasion is recombinase-driven rather than
# annealing-driven, so inter-primer Tm matching matters relatively little, while
# primer-dimer artefacts compete directly for enzyme, primer and energy.
#
# The result is a heuristic ranking aid with no experimental calibration. It picks
# which combinations are worth synthesising; it is not a probability and must not
# be reported as a percentage.
SCORE_WEIGHT_CROSS_DIMER = 4.0
SCORE_WEIGHT_DELTA_TM = 3.0


def compatibility_score(min_dg: float, delta_tm: float, dye_ok: bool = True) -> float:
    """Heuristic 0–100 multiplex compatibility score (higher is better).

    Shared by the pairwise and triplex checkers so their numbers mean the same
    thing. A dye-channel collision zeroes the score: the targets are optically
    inseparable and no thermodynamic quality can compensate.
    """
    if not dye_ok:
        return 0.0
    score = 100.0
    if min_dg < 0:
        score += min_dg * SCORE_WEIGHT_CROSS_DIMER
    score -= delta_tm * SCORE_WEIGHT_DELTA_TM
    return max(0.0, round(score, 1))


def screening_methodology_note() -> str:
    """Markdown block explaining what the thresholds mean and where they come from.

    Appended verbatim to both the pairwise and triplex reports so a bench reader
    gets the same caveats either way, and so the justification lives next to the
    constants it describes rather than being restated in two report writers.
    """
    return f"""
---

## How to Read These Numbers

### Compatibility score
A heuristic ranking aid on a 0–100 scale, **not a percentage and not a
probability**:

```
score = 100 + min_dG x {SCORE_WEIGHT_CROSS_DIMER} - dTm x {SCORE_WEIGHT_DELTA_TM}   (0 if the dye channels collide)
```

It has no experimental calibration. Its job is to decide which combinations are
worth synthesising, not to predict assay performance. Pairwise and triplex
reports use identical weights, so their scores are directly comparable.

### Why the ΔG thresholds look loose next to PCR practice
All ΔG values here are computed at the **39°C RAA reaction temperature**, not at
a qPCR temperature, and ΔG = ΔH − TΔS means the same duplex is markedly more
stable at 39°C. Measured over 200 random 32-mer pairs with no designed
complementarity:

| | median | 5th pct | below −9 | below −3 |
|---|---|---|---|---|
| 39°C (RAA) | −4.12 | −6.78 | 0.5% | 88.0% |
| 60°C (qPCR) | −1.33 | −3.07 | 0.0% | 6.5% |

A "−3 kcal/mol" cutoff carried over from 60°C practice would reject 88% of
entirely clean pairs. The {MIN_CROSS_DIMER_DG} kcal/mol threshold used here
corresponds to roughly −4.7 at 60°C — already stricter than the general IDT
OligoAnalyzer guidance it derives from.

There is **no validated RPA/RAA-specific ΔG threshold** in the literature. The
TwistAmp Assay Design Manual §2.2 states that "there are no fixed rules to
predict how well a given amplification primer will work" and prescribes
empirical screening of 8–10 primers per direction instead.

### 3'-end anchoring
Only a 3' terminus that anneals can be extended by the polymerase, and those
extension products are what amplify as *primer noise* — which the TwistAmp manual
§2.6 identifies as the factor limiting a primer pair's sensitivity.

The criterion is the **length of the longest perfectly complementary run anchored
at a 3' terminus**, not a ΔG. More than {MAX_THREE_PRIME_ANCHOR_BP} bp is flagged
`WARNING`. Measured against the same random-pair baseline:

| longest 3' anchor | share of clean pairs |
|---|---|
| ≥4 bp | 16.8% |
| ≥5 bp | 4.2% |
| ≥6 bp | 0.5% |

Deliberate 5, 6 and 8 bp anchors are detected in 100% of cases. A ΔG over a short
terminal window was tried first and rejected — a 5 nt window saturates near −6.3
kcal/mol however much complementarity is present, and an 8 nt window overlaps the
clean background — so the ΔG column is reported for context only and nothing is
gated on it.

### ΔTm
Tolerance is {MAX_DELTA_TM_C}°C. This is a soft criterion for RAA: at 39°C strand
invasion is recombinase-driven rather than annealing-driven, so inter-primer Tm
matching matters far less than it does in PCR. It is retained mainly for the 55°C
qPCR cross-check.

### Before ordering oligos
Treat this ranking as a shortlist, not a verdict. Following the vendor's own
protocol, order **several** of the top-ranked combinations rather than rank 1
alone — no in-silico calculation predicts primer noise reliably, and the decisive
comparison is sensitivity on the bench.
"""


_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    return (seq or "").upper().translate(_COMPLEMENT)[::-1]


def three_prime_anchor_length(
    seq1: str,
    seq2: str,
    max_len: int = 12,
    max_variants: int = DEFAULT_MAX_VARIANTS,
) -> int:
    """Longest perfectly complementary run anchored at either oligo's 3' terminus.

    A run of length k means the last k bases of one oligo can anneal somewhere on
    the other — terminally or internally — giving the polymerase a 3'-OH to extend
    from. The extension product is a primer artefact that then amplifies on its
    own and competes with the real target for enzyme, primer and energy.

    Degenerate oligos are expanded, and the WORST variant decides: a synthesised
    pool contains all of them, so an anchor present in any one variant is present
    in the tube.
    """
    if not seq1 or not seq2:
        return 0
    best = 0
    for tail_src, partner in ((seq1, seq2), (seq2, seq1)):
        partners = expand_iupac(partner, max_variants=max_variants)
        limit = min(max_len, len(tail_src))
        for k in range(1, limit + 1):
            tails = expand_iupac(tail_src[-k:], max_variants=max_variants)
            if any(reverse_complement(t) in p for t in tails for p in partners):
                best = max(best, k)
    return best


def evaluate_three_prime_dimer(
    seq1: str,
    seq2: str,
    temp_c: float,
    window: int = THREE_PRIME_WINDOW_NT,
    max_variants: int = DEFAULT_MAX_VARIANTS,
) -> float:
    """ΔG for either oligo's 3' terminal `window` bases annealing to the other.

    Reported alongside three_prime_anchor_length() as context for a reader
    comparing two flagged pairs. Do NOT threshold on this value — a short terminal
    window saturates (a 5 nt window bottoms out near −6.3 kcal/mol however much
    complementarity exists) and an 8 nt window overlaps the clean background, so
    it cannot separate designed complementarity from chance. The anchor length is
    the criterion.

    Returns the most negative ΔG in kcal/mol; 0.0 if no structure.
    """
    try:
        import primer3
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "primer3-py is required. Install via: pip install primer3-py"
        ) from exc

    if not seq1 or not seq2:
        return 0.0

    min_dg = 0.0
    for tail_src, partner in ((seq1, seq2), (seq2, seq1)):
        tails = expand_iupac(tail_src[-window:], max_variants=max_variants)
        partners = expand_iupac(partner, max_variants=max_variants)
        for t in tails:
            for p in partners:
                key = (t, p, temp_c)
                if key in _raw_p3_cache:
                    dg_kcal, sf = _raw_p3_cache[key]
                else:
                    res = primer3.calc_heterodimer(t, p, temp_c=temp_c)
                    dg_kcal, sf = res.dg / 1000.0, res.structure_found
                    _raw_p3_cache[key] = (dg_kcal, sf)
                    _raw_p3_cache[(p, t, temp_c)] = (dg_kcal, sf)
                if sf and dg_kcal < min_dg:
                    min_dg = dg_kcal
    return min_dg


# ---------------------------------------------------------------------------
# Spectral Dye Channel Separation
# ---------------------------------------------------------------------------
# Maps common qPCR fluorophores to their detection channel. Two targets read on
# the SAME channel cannot be distinguished in a multiplex, so they collide.
DYE_CHANNELS = {
    "FAM": "green", "SYBR": "green", "SYTO9": "green", "ALEXA488": "green",
    "HEX": "yellow", "VIC": "yellow", "JOE": "yellow", "TET": "yellow", "YAKIMA": "yellow",
    "CY3": "orange", "NED": "orange", "TAMRA": "orange", "ATTO550": "orange",
    "ROX": "red", "TEXASRED": "red", "CALRED610": "red", "CALFLUORRED610": "red",
    "CY5": "farred", "QUASAR670": "farred", "ATTO647N": "farred",
    "CY5.5": "infrared", "QUASAR705": "infrared",
}


def _normalize_dye(dye: str) -> str:
    return (dye or "").upper().replace(" ", "").replace("-", "").replace("_", "")


def check_dye_separation(*dyes: str) -> Tuple[str, bool]:
    """Verify that every target in a multiplex is read on a distinct spectral channel.

    Returns (human_readable_status, ok). ``ok`` is False if any two dyes map to the
    same detection channel (a collision — targets become indistinguishable) or if a
    dye is not in DYE_CHANNELS (cannot be verified — reported, never silently passed).
    """
    if not dyes:
        return ("UNKNOWN (no dyes supplied; cannot verify spectral separation)", False)

    channels = []
    unknown = []
    for d in dyes:
        # An empty/None dye is MISSING information, not an absent target. Dropping
        # it would let check_dye_separation("FAM", "") report a clean separation
        # for a target whose channel nobody ever specified.
        if not d or not str(d).strip():
            unknown.append("<empty>")
            channels.append(("<empty>", None))
            continue
        ch = DYE_CHANNELS.get(_normalize_dye(d))
        if ch is None:
            unknown.append(d)
        channels.append((d, ch))

    if unknown:
        return (f"UNKNOWN (dye(s) not in channel table: {', '.join(unknown)}; "
                f"cannot verify spectral separation)", False)

    # Detect any shared channel between two distinct targets.
    seen = {}
    for d, ch in channels:
        if ch in seen:
            return (f"FAIL ({seen[ch]} & {d} share the '{ch}' channel — "
                    f"targets cannot be distinguished)", False)
        seen[ch] = d

    labels = ", ".join(f"{d}={ch}" for d, ch in channels)
    return (f"PASSED (distinct channels: {labels})", True)
