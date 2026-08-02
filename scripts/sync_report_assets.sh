#!/usr/bin/env bash
#
# Refresh the interactive run artifacts that the Quarto report links to, on the
# cf-pages deploy branch, in one command.
#
# Why this exists: results/** is gitignored derived data, so the publish workflow
# (which renders on a runner without results/) can never regenerate it. publish.yml
# uses keep_files:true so a publish no longer wipes these assets — but when a run
# changes you still need to push the new assets once. That is this script.
#
# It syncs, for every influenza-a-mp-* build found under the results dir:
#   results/<build>/m10/**                         (comparison_report.html, heatmaps)
# plus, for indo-only, the linked M09 BLAST summaries. It also relabels the
# dataset toggle in any comparison_report.html that still carries the old
# hardcoded "Design DB (Indonesia 645)" / "Global Surveillance DB (50k)" text,
# using the real per-build counts from coverage_summary.json.
#
# Usage:
#   scripts/sync_report_assets.sh [--results-dir DIR] [--remote NAME]
#     --results-dir  where the run outputs live (default: ./results)
#     --remote       git remote whose cf-pages feeds the site (default: origin)
#
set -euo pipefail

RESULTS_DIR="results"
REMOTE="origin"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-dir) RESULTS_DIR="$2"; shift 2 ;;
    --remote) REMOTE="$2"; shift 2 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

REPO="$(git rev-parse --show-toplevel)"
cd "$REPO"
[[ -d "$RESULTS_DIR" ]] || { echo "ERROR: results dir not found: $RESULTS_DIR" >&2; exit 1; }

WT="$(mktemp -d)/cf-pages"
cleanup() { git worktree remove "$WT" --force 2>/dev/null || true; }
trap cleanup EXIT

echo "Fetching ${REMOTE}/cf-pages ..."
git fetch "$REMOTE" cf-pages
git worktree add --detach "$WT" "${REMOTE}/cf-pages" >/dev/null

# --- copy linked assets + relabel toggles (label from coverage_summary data) ---
python3 - "$RESULTS_DIR" "$WT" <<'PY'
import json, glob, os, shutil, sys
results_dir, wt = sys.argv[1], sys.argv[2]

def counts(build):
    """Return (total_seqs, design_db_total) from the degenerate coverage summary."""
    p = glob.glob(os.path.join(results_dir, build, "m10", "degenerate", "coverage_summary.json"))
    if not p:
        return None, None
    c = json.load(open(p[0]))
    c = c[0] if isinstance(c, list) else c
    return c.get("total_seqs"), c.get("design_db_total")

def relabel(html_path, build):
    if not os.path.exists(html_path):
        return
    g, d = counts(build)
    t = open(html_path, encoding="utf-8").read()
    if g is not None:
        t = t.replace("🌍 Global Surveillance DB (50k)", f"🌍 Global Surveillance DB ({g:,})")
    if d is not None:
        t = t.replace("🎯 Design DB (Indonesia 645)", f"🎯 Design DB ({d:,} seqs)")
    open(html_path, "w", encoding="utf-8").write(t)

synced = []
for src in sorted(glob.glob(os.path.join(results_dir, "influenza-a-mp-*", "m10"))):
    build = os.path.basename(os.path.dirname(src))
    dst = os.path.join(wt, "results", build, "m10")
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(src, dst)
    relabel(os.path.join(dst, "comparison_report.html"), build)
    synced.append(build + "/m10")

# indo-only M09 BLAST summaries linked from the candidates page
d9 = glob.glob(os.path.join(results_dir, "influenza-a-mp-indo-only", "m09",
                            "degenerate", "*", "blast_specificity_summary.csv"))
if d9:
    srcdir = os.path.dirname(d9[0])
    dstdir = os.path.join(wt, os.path.relpath(srcdir, "."))  # mirror path under wt
    # rebuild the results/... path explicitly (relpath from results_dir)
    rel = os.path.relpath(srcdir, results_dir)
    dstdir = os.path.join(wt, "results", rel)
    os.makedirs(dstdir, exist_ok=True)
    for f in ("blast_specificity_summary.csv", "blast_specificity_summary.md",
              "final_candidates_report.md"):
        s = os.path.join(srcdir, f)
        if os.path.exists(s):
            shutil.copy2(s, os.path.join(dstdir, f))
    synced.append("indo-only/m09 blast summaries")

print("synced: " + ", ".join(synced) if synced else "nothing to sync")
PY

cd "$WT"
git add results/
if git diff --cached --quiet; then
  echo "No asset changes — cf-pages already up to date."
  exit 0
fi
git -c user.name="engkinandatama" -c user.email="kaysan0200@gmail.com" \
    commit -q -m "deploy: refresh report-linked M10/M09 assets on cf-pages"
git push "$REMOTE" HEAD:cf-pages
echo "Pushed refreshed assets to ${REMOTE}/cf-pages. Cloudflare will rebuild; hard-refresh the site."
