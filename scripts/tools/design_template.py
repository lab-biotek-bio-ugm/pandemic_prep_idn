#!/usr/bin/env python3
"""
Shared M08 design-template resolution.

Under the "full segment" strategy the payload generator does NOT design each
candidate region against its own consensus: it substitutes the single longest
region's consensus as the Primer3 template for every region, so Primer3's
reported start/end are relative to THAT sequence.

The degenerate rescue step then back-maps those coordinates onto a consensus
slice. If it re-derives the sequence per region (as it used to), it slices a
different — usually shorter — sequence than the one Primer3 actually saw, and
silently produces primers with wrong positions or truncated content.

Both scripts now resolve the template through this single implementation so the
two can never drift apart.
"""
from __future__ import annotations

from typing import List, Tuple

# A manifest is treated as "full segment" when region lengths differ markedly or
# any region is long enough to be a whole segment rather than a window.
FULL_SEGMENT_LENGTH_SPREAD = 100
FULL_SEGMENT_MIN_LENGTH = 500


def _consensus_of(region: dict) -> str:
    return region.get("sequences", {}).get("consensus_sequence", "") or ""


def detect_full_segment_strategy(candidate_regions: List[dict]) -> bool:
    """True when the manifest describes one full segment rather than windows."""
    lengths = [len(_consensus_of(r)) for r in candidate_regions if _consensus_of(r)]
    if not lengths:
        return False
    return (max(lengths) - min(lengths) > FULL_SEGMENT_LENGTH_SPREAD) or (
        max(lengths) > FULL_SEGMENT_MIN_LENGTH
    )


def template_region(candidate_regions: List[dict]) -> dict:
    """The region whose consensus is used as the Primer3 design template."""
    return max(candidate_regions, key=lambda r: len(_consensus_of(r)))


def resolve_design_template(candidate_regions: List[dict], region_id: str) -> Tuple[str, str, str]:
    """Return (consensus_sequence, iupac_consensus_sequence, template_region_id).

    Under the full-segment strategy every region resolves to the longest region's
    sequences — matching what Primer3 was given — otherwise a region uses its own.
    """
    by_id = {r.get("region_id"): r for r in candidate_regions}
    own = by_id.get(region_id)
    if own is None:
        raise KeyError(f"region_id {region_id!r} not present in manifest candidate regions")

    if detect_full_segment_strategy(candidate_regions):
        tpl = template_region(candidate_regions)
    else:
        tpl = own

    seqs = tpl.get("sequences", {})
    return (
        seqs.get("consensus_sequence", ""),
        seqs.get("iupac_consensus_sequence", ""),
        tpl.get("region_id", region_id),
    )
