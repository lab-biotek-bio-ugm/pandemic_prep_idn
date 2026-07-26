#!/usr/bin/env python3
"""M12 — coverage-failure attribution.

Answers, for each subtype the universal assay covers poorly (<90%), the auditable
question WHY: which oligo, which position (3'-weighted), is it a dominant variant
or diffuse divergence, and where do the failing strains come from (clade/country).

It recomputes on the failing subset only (a few thousand sequences), so it is
cheap and never re-runs the hour-long M10. It reads only artifacts a run already
wrote, and uses the same IUPAC-aware alignment core as M10 (best_binding), so the
mismatch semantics match. Deterministic; statistics are numpy/stdlib only (no
scipy, no LLM/API) so the tool runs in any lab.

The comparison that carries the causal story: within a failing subtype, sequences
whose clade the design set has SEEN (focal) vs clades it has NOT (divergent).
A primer designed on Indonesian H5N1 (older clades) that fails on the globally
dominant 2.3.4.4b shows up here as a 3'-terminal base that is conserved in the
focal group and switched in the divergent group.
"""
import argparse
import csv
import json
import math
import os
import sys
from collections import Counter
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
from m10_insilico_pcr import best_binding, find_metadata_for_header, match_base  # noqa: E402

DETECTED_FLOOR = 90.0        # a subtype below this is "failing" and gets attributed
MIN_STRATUM_SEQS = 30        # below this, too sparse to attribute with confidence
DRIVER_DIVERGENT_RATE = 0.5  # a position where >=50% of divergent seqs mismatch...
DRIVER_FOCAL_RATE = 0.2      # ...while <=20% of focal seqs do, is a "driver"
_ID_COLS = ("strain", "isolate_id", "isolate_name", "gisaid_isolate_id", "accession")


# --------------------------------------------------------------------- loading

def load_metadata(path):
    """strain/isolate id -> {subtype, clade, country}. Keyed by every id column so
    a GISAID EPI_ISL_ FASTA header joins regardless of which column holds it."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        sample = fh.readline()
    delim = "\t" if sample.count("\t") >= sample.count(",") else ","
    out = {}
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        lower = {c: (c or "").strip().lower() for c in (reader.fieldnames or [])}
        for row in reader:
            r = {lower[k]: (v or "").strip() for k, v in row.items() if k in lower}
            meta = {
                "subtype": r.get("subtype") or "Unknown",
                "clade": r.get("clade") or "Unknown",
                "country": r.get("country") or r.get("source_country") or "Unknown",
            }
            for col in _ID_COLS:
                key = r.get(col)
                if key:
                    out.setdefault(key, meta)
    return out


def iter_fasta(path):
    header, chunks = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def _meta_of(header, strain_to_meta):
    return find_metadata_for_header(header, strain_to_meta) or {}


def load_group_sequences(fasta_path, strain_to_meta, wanted_subtypes):
    """subtype -> list of {seq, clade, country} for the wanted subtypes."""
    groups = {s: [] for s in wanted_subtypes}
    for header, seq in iter_fasta(fasta_path):
        m = _meta_of(header, strain_to_meta)
        st = m.get("subtype", "Unknown")
        if st in groups and seq:
            groups[st].append({"seq": seq.upper(),
                               "clade": m.get("clade", "Unknown"),
                               "country": m.get("country", "Unknown")})
    return groups


def design_clades_by_subtype(design_fasta, strain_to_meta):
    """subtype -> set of clades present in the design set (the 'seen' context)."""
    seen = {}
    for header, _seq in iter_fasta(design_fasta):
        m = _meta_of(header, strain_to_meta)
        st, cl = m.get("subtype", "Unknown"), m.get("clade", "Unknown")
        seen.setdefault(st, set()).add(cl)
    return seen


# ------------------------------------------------------------------ statistics

def g_test_2x2(a, b, c, d):
    """G-test (log-likelihood ratio) for a 2x2 table, 1 df. Returns (G, p).

    Chosen over Fisher exact because n here is large (hundreds-thousands), where
    Fisher is expensive and the G approximation is excellent. The 1-df chi-square
    survival function is exact via erfc, so no scipy is needed.
    """
    n = a + b + c + d
    if n == 0:
        return 0.0, 1.0
    row1, row2, col1, col2 = a + b, c + d, a + c, b + d
    g = 0.0
    for obs, (rr, cc) in zip((a, b, c, d),
                             ((row1, col1), (row1, col2), (row2, col1), (row2, col2))):
        exp = rr * cc / n
        if obs > 0 and exp > 0:
            g += obs * math.log(obs / exp)
    g *= 2.0
    p = math.erfc(math.sqrt(g / 2.0)) if g > 0 else 1.0
    return g, p


def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted p-values (stdlib)."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [1.0] * m
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        k = m - rank + 1
        prev = min(prev, pvals[idx] * m / k)
        adj[idx] = prev
    return adj


# ------------------------------------------------------------------ attribution

def _tvd(counter_a, na, counter_b, nb):
    """Total-variation distance between two base-frequency distributions."""
    if na == 0 or nb == 0:
        return 0.0
    bases = set(counter_a) | set(counter_b)
    return 0.5 * sum(abs(counter_a.get(x, 0) / na - counter_b.get(x, 0) / nb) for x in bases)


def _align(oligo, seq):
    """Return (aligned_bases | None, mismatch_positions) in the oligo's 5'->3' frame."""
    _c, offset, _s, positions, aligned = best_binding(oligo, seq)
    if offset is None:
        return None, []
    return aligned, positions


def _raa_detected(aln, crit, max3p, maxint, probe_max):
    """Mechanistic RAA detection gate from alignments (no seqkit).

    Both primers must bind with <=max3p mismatches in the 3'-critical window and
    <=maxint elsewhere; the probe on total load. This mirrors M10's RAA model
    (assess_raa_binding + assess_probe_binding) closely enough to split the
    subtype into 'detected' vs 'not detected' for attribution.
    """
    for name in ("forward", "reverse"):
        aligned, mm_pos, ln = aln[name]
        if aligned is None:
            return False
        n3 = sum(1 for j in mm_pos if j >= ln - crit)
        if n3 > max3p or (len(mm_pos) - n3) > maxint:
            return False
    aligned, mm_pos, _ln = aln["probe"]
    if aligned is None or len(mm_pos) > probe_max:
        return False
    return True


def _profile_from_aligned(oligo, aligned_list):
    """Per-position base counts + oligo-mismatch counts over aligned-base lists."""
    ln = len(oligo)
    base_counts = [Counter() for _ in range(ln)]
    mismatch = [0] * ln
    used = 0
    for aligned in aligned_list:
        if aligned is None:
            continue
        used += 1
        for j in range(ln):
            b = aligned[j]
            base_counts[j][b] += 1
            if not match_base(oligo[j], b):
                mismatch[j] += 1
    return base_counts, mismatch, used


def attribute_oligo(name, oligo, detected_aligned, failing_aligned, crit_3prime):
    """Contrast the DETECTED vs NOT-DETECTED sequences of a subtype at this oligo.

    The driver positions are where the failing group mismatches the oligo and the
    detected group does not — the mechanistic cause of the coverage gap.
    """
    ln = len(oligo)
    dc, dm, dn = _profile_from_aligned(oligo, detected_aligned)   # detected
    fc, fm, fn = _profile_from_aligned(oligo, failing_aligned)    # not detected
    positions, pvals = [], []
    for j in range(ln):
        det_rate = dm[j] / dn if dn else 0.0
        fail_rate = fm[j] / fn if fn else 0.0
        _g, p = g_test_2x2(dn - dm[j], dm[j], fn - fm[j], fm[j])
        pvals.append(p)
        positions.append({
            "pos": j, "three_prime": j >= ln - crit_3prime,
            "detected_mismatch_rate": round(det_rate, 4),
            "failing_mismatch_rate": round(fail_rate, 4),
            "tvd": round(_tvd(dc[j], dn, fc[j], fn), 4),
            "detected_bases": dict(dc[j].most_common(3)),
            "failing_bases": dict(fc[j].most_common(3)),
            "is_driver": (fail_rate >= DRIVER_DIVERGENT_RATE and det_rate <= DRIVER_FOCAL_RATE),
            "g_pvalue": p,
        })
    for pos, adj in zip(positions, bh_fdr(pvals)):
        pos["g_pvalue_fdr"] = round(adj, 6)
        pos["g_pvalue"] = round(pos["g_pvalue"], 6)
    drivers = [p for p in positions if p["is_driver"]]
    return {
        "oligo": name, "sequence": oligo, "length": ln,
        "detected_n": dn, "failing_n": fn,
        "binding_site_divergence_tvd": round(
            sum(p["tvd"] for p in positions) / ln if ln else 0.0, 4),
        "driver_positions": drivers,
        "has_3prime_driver": any(p["three_prime"] for p in drivers),
        "positions": positions,
    }


def classify(oligo_results):
    """dominant_snp if a few (esp. 3') positions drive it; diffuse otherwise."""
    all_drivers = [d for r in oligo_results for d in r["driver_positions"]]
    any_3p = any(r["has_3prime_driver"] for r in oligo_results)
    if any_3p and len(all_drivers) <= 3:
        return "dominant_snp"
    if all_drivers:
        return "diffuse_divergence" if len(all_drivers) > 5 else "dominant_snp"
    return "diffuse_divergence"


# ------------------------------------------------------------------- driver

def _pick_candidate(summary, rank):
    for c in summary:
        if c.get("rank") == rank:
            return c
    return summary[0] if summary else None


def attribute_build(results_dir, build, rank, mode="degenerate", floor=DETECTED_FLOOR,
                    max_mismatch=3):
    base = Path(results_dir) / build
    cov_path = base / "m10" / mode / "coverage_summary.json"
    summary = json.loads(cov_path.read_text(encoding="utf-8"))
    cand = _pick_candidate(summary, rank)
    if not cand:
        raise SystemExit(f"No candidate rank {rank} in {cov_path}")

    oligos = [
        ("forward", cand.get("fwd_degenerate_seq") or cand.get("fwd_seq") or ""),
        ("reverse", cand.get("rev_degenerate_seq") or cand.get("rev_seq") or ""),
        ("probe", cand.get("probe_seq") or ""),
    ]
    raa = cand.get("raa_params") or {}
    crit = int(raa.get("critical_3prime_bp", 4))

    failing = [b for b in (cand.get("breakdowns", {}).get("subtype") or [])
               if b.get("raa_detected_pct") is not None
               and b["raa_detected_pct"] < floor
               and not str(b.get("subtype", "")).startswith("Other")]

    max3p = int(raa.get("max_3prime_mm", 1))
    maxint = int(raa.get("max_internal_mm", 4))
    probe_max = int(raa.get("probe_max_mm", 5))

    strain_to_meta = load_metadata(base / "m01" / "metadata_enriched.tsv")
    design_seen = design_clades_by_subtype(base / "m02" / "filtered_sequences.fasta",
                                           strain_to_meta)
    wanted = {b["subtype"] for b in failing}
    seqs_by_subtype = load_group_sequences(base / "m01" / "raw_sequences.fasta",
                                           strain_to_meta, wanted)

    reports = []
    for b in sorted(failing, key=lambda x: x["raa_detected_pct"]):
        st = b["subtype"]
        seqs = seqs_by_subtype.get(st, [])
        seen = {c for c in design_seen.get(st, set()) if c != "Unknown"}
        rec = {"subtype": st, "coverage_pct": b["raa_detected_pct"],
               "n_global": b.get("total"), "n_analysed": len(seqs)}

        if not seen:
            rec["failure_mode"] = "no_focal_data"
            rec["explanation"] = ("This subtype/its clades are absent from the design "
                                  "set, so a universal primer built without it cannot "
                                  "be expected to cover it. Not a design defect.")
            reports.append(rec)
            continue
        if len(seqs) < MIN_STRATUM_SEQS:
            rec["failure_mode"] = "too_sparse"
            rec["explanation"] = f"Only {len(seqs)} sequences; too few to attribute confidently."
            reports.append(rec)
            continue

        # Align all three oligos to every sequence once; gate each into
        # detected / not-detected (the mechanistic RAA model, no seqkit); then
        # contrast the two groups. This asks the right question -- why do THESE
        # sequences fail -- rather than assuming an unseen clade is to blame.
        # Align all three oligos to every sequence once, then apply TWO lenses:
        #   strict-PCR: total primer mismatch <= max_mismatch (reproduces M10's
        #               seqkit gate, which is what caps the reported number);
        #   RAA-viable: the RAA tolerance model (<=1 3' + <=4 internal), which is
        #               the actual assay chemistry.
        # The gap between them is the story; the residue that fails even RAA is the
        # truly-uncoverable set, attributed by driver position below.
        detected = {nm: [] for nm, _ in oligos}
        failing_al = {nm: [] for nm, _ in oligos}
        n_raa, n_pcr = 0, 0
        fail_recs = []
        for s in seqs:
            # aln[name] = (aligned_bases|None, mismatch_positions, oligo_len)
            aln = {nm: (_align(seq, s["seq"]) + (len(seq),)) for nm, seq in oligos if seq}
            primer_mm = sum(len(aln[nm][1]) for nm in ("forward", "reverse") if nm in aln)
            if primer_mm <= max_mismatch:
                n_pcr += 1
            ok = _raa_detected(aln, crit, max3p, maxint, probe_max)
            bucket = detected if ok else failing_al
            for nm, _seq in oligos:
                if nm in aln:
                    bucket[nm].append(aln[nm][0])
            if ok:
                n_raa += 1
            else:
                fail_recs.append(s)

        n = len(seqs)
        rec["reported_pct"] = b["raa_detected_pct"]        # M10, PCR-capped
        rec["pcr_strict_pct"] = round(100.0 * n_pcr / n, 2)   # reproduces ~M10
        rec["raa_viable_pct"] = round(100.0 * n_raa / n, 2)   # RAA tolerance
        rec["failing_n"] = len(fail_recs)
        rec["detected_n"] = n_raa
        rec["seen_clades"] = sorted(seen)
        rec["failing_top_clades"] = Counter(s["clade"] for s in fail_recs).most_common(5)
        rec["failing_top_countries"] = Counter(s["country"] for s in fail_recs).most_common(5)

        if len(fail_recs) < MIN_STRATUM_SEQS or n_raa < 5:
            # Almost everything is RAA-viable: the low reported number is a
            # strict-PCR (-m{max}) artifact, not a real coverage hole.
            rec["failure_mode"] = "raa_rescued"
            rec["verdict"] = "strict_pcr_artifact"
            rec["explanation"] = _narrate(rec)
            reports.append(rec)
            continue

        oligo_results = [attribute_oligo(nm, seq, detected[nm], failing_al[nm], crit)
                         for nm, seq in oligos if seq]
        rec["failure_mode"] = "mechanistic"
        rec["verdict"] = classify(oligo_results)
        rec["oligos"] = oligo_results
        rec["explanation"] = _narrate(rec)
        reports.append(rec)

    return {"build": build, "candidate_rank": cand.get("rank"),
            "detected_floor": floor, "max_mismatch": max_mismatch,
            "finding_note": (
                "The reported coverage is gated by M10's seqkit in-silico PCR at "
                f"max_mismatch={max_mismatch}. For an RAA assay that is conservative: "
                "RPA/RAA tolerates more mismatches than PCR (TwistAmp: 5-9 bp), so a "
                "'raa_viable_pct' well above the reported number means the assay likely "
                "works despite a low PCR-strict score. This is flagged, not fixed; "
                "changing M10's coverage model is a separate decision."),
            "stats_note": ("Effect sizes (mismatch-rate difference, per-position TVD) "
                           "are primary. G-test p-values are FDR-corrected but "
                           "ANTICONSERVATIVE: viral sequences are phylogenetically "
                           "non-independent, so treat them as ranking aids, not proof."),
            "subtypes": reports}


def _narrate(rec):
    top_clade = rec["failing_top_clades"][0] if rec.get("failing_top_clades") else ("?", 0)
    top_country = rec["failing_top_countries"][0] if rec.get("failing_top_countries") else ("?", 0)
    lead = (f"{rec['subtype']}: reported {rec['reported_pct']}% (strict-PCR gate reproduces "
            f"{rec['pcr_strict_pct']}%), but RAA-viable {rec['raa_viable_pct']}% under RAA "
            f"mismatch tolerance. The low reported number is largely a strict-PCR artifact: "
            f"most strains carry more mismatches than PCR allows but within RAA's tolerance.")

    if rec.get("failure_mode") == "raa_rescued":
        return (lead + f" Only {rec['failing_n']} strains fail even the RAA model — too few "
                "to attribute a driver; the coverage 'gap' is essentially a PCR-strictness "
                "artifact, not a real RAA hole.")

    gap = (f" Of the ~{rec['failing_n']} strains that fail even RAA, they concentrate in "
           f"clade {top_clade[0]} (n={top_clade[1]}) and {top_country[0]} (n={top_country[1]}).")

    def _fmt(bases):
        tot = max(sum(bases.values()), 1)
        return "/".join(f"{k} {round(100 * v / tot)}%" for k, v in bases.items())

    drivers = [(o["oligo"], d) for o in rec.get("oligos", []) for d in o["driver_positions"]]
    if not drivers:
        return (lead + gap + " No single position separates them — divergence is diffuse, "
                "so recovering them needs a dedicated assay, not a base tweak.")
    parts = []
    for oligo, d in drivers[:4]:
        where = "3'-terminal" if d["three_prime"] else "internal"
        parts.append(f"{oligo} {where} position {d['pos']}: RAA-detected [{_fmt(d['detected_bases'])}] "
                     f"vs RAA-failing [{_fmt(d['failing_bases'])}] (mismatch "
                     f"{round(100 * d['failing_mismatch_rate'])}% vs "
                     f"{round(100 * d['detected_mismatch_rate'])}%)")
    return lead + gap + " Driver(s): " + "; ".join(parts) + "."


# ---------------------------------------------------------------------- output

def write_markdown(result, path):
    lines = [f"# Coverage-failure attribution — {result['build']} (candidate rank "
             f"{result['candidate_rank']})", "",
             f"> {result.get('finding_note', '')}", "",
             f"_{result['stats_note']}_", ""]
    for r in result["subtypes"]:
        lines.append(f"## {r['subtype']} — reported {r.get('reported_pct', r.get('coverage_pct'))}% "
                     f"({r['failure_mode']})")
        lines.append("")
        lines.append(r.get("explanation", ""))
        lines.append("")
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description="Attribute why the universal assay under-covers subtypes.")
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--build", required=True)
    ap.add_argument("--rank", type=int, default=1, help="Coverage-summary candidate rank to attribute.")
    ap.add_argument("--mode", default="degenerate")
    ap.add_argument("--floor", type=float, default=DETECTED_FLOOR)
    ap.add_argument("--max-mismatch", type=int, default=3,
                    help="M10's seqkit PCR mismatch cap; the strict-PCR lens reproduces it.")
    ap.add_argument("--out", required=True, help="Output directory.")
    args = ap.parse_args()

    result = attribute_build(args.results_dir, args.build, args.rank, args.mode,
                             args.floor, args.max_mismatch)
    os.makedirs(args.out, exist_ok=True)
    (Path(args.out) / "coverage_attribution.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8")
    write_markdown(result, Path(args.out) / "coverage_attribution.md")
    mech = sum(1 for r in result["subtypes"] if r["failure_mode"] == "mechanistic")
    print(f"Attributed {len(result['subtypes'])} failing subtype(s) "
          f"({mech} mechanistic) -> {args.out}")


if __name__ == "__main__":
    main()
