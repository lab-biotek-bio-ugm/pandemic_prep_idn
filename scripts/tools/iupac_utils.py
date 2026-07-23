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
    """Calculate GC content percentage from an IUPAC sequence.

    S (G/C) counts as a G/C base.

    Parameters
    ----------
    seq : str
        IUPAC or canonical sequence.

    Returns
    -------
    float
        GC percentage (0.0 – 100.0), rounded to 1 decimal place.
    """
    seq = seq.upper()
    gc = sum(1 for b in seq if b in ('G', 'C', 'S'))
    return round((gc / len(seq)) * 100.0, 1) if seq else 0.0


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
