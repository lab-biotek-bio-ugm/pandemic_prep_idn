"""Checks on the workflow files themselves.

The repository guard is the one piece of CI configuration with a cost attached to
getting it wrong. Actions minutes on the private lab repository are metered, and a
job that runs on the wrong repository spends the wrong quota. Each workflow is
deliberately pinned to a single repository:

  - ci.yml (unit tests) runs on the LAB upstream, where the supervisor wants the
    tests visible; it is skipped on the personal fork.
  - publish.yml (Cloudflare deploy) runs on the personal FORK, where the
    Cloudflare Pages integration is wired; it must never publish to the lab's
    cf-pages branch.

GitHub Actions has no workflow-level `if`, so the guard has to be repeated on
every job, and a job added later is exactly the kind of thing that gets it wrong.
"""
import pytest
import yaml

from conftest import REPO

pytestmark = pytest.mark.contract

WORKFLOW_DIR = REPO / ".github" / "workflows"
WORKFLOWS = sorted(p for p in WORKFLOW_DIR.glob("*.y*ml"))
FORK = "engkinandatama/pandemic-preparedness"
LAB = "lab-biotek-bio-ugm/pandemic_prep_idn"
# Which repository each workflow is pinned to run on.
EXPECTED_GUARD = {
    "ci.yml": f"github.repository == '{LAB}'",
    "publish.yml": f"github.repository == '{FORK}'",
}


def _load(path):
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def test_there_is_at_least_one_workflow():
    assert WORKFLOWS, f"no workflow files found under {WORKFLOW_DIR}"


@pytest.mark.parametrize("wf", WORKFLOWS, ids=lambda p: p.name)
def test_every_job_carries_the_repository_guard(wf):
    jobs = _load(wf).get("jobs", {})
    assert jobs, f"{wf.name} defines no jobs"

    guard = EXPECTED_GUARD.get(wf.name)
    assert guard, (
        f"{wf.name} has no expected repository guard registered in EXPECTED_GUARD; "
        f"add one so a new workflow cannot run on the wrong repository unnoticed."
    )

    unguarded = [name for name, job in jobs.items() if guard not in str(job.get("if", ""))]
    assert not unguarded, (
        f"{wf.name}: these jobs are missing the pinning guard and would run on the "
        f"wrong repository: {unguarded}\n"
        f"Add to each one:\n    if: {guard}"
    )


@pytest.mark.parametrize("wf", WORKFLOWS, ids=lambda p: p.name)
def test_no_yaml_anchors(wf):
    """GitHub Actions rejects anchors and aliases outright, while PyYAML resolves
    them happily -- so a plain YAML check reports green on a workflow GitHub
    refuses to parse, and the failure mode is that nothing runs at all.

    Detected from the parser's event stream rather than by grepping the text: a
    substring search for " &" also matches the `&&` in any shell `run:` line.
    """
    offenders = []
    for event in yaml.parse(wf.read_text(encoding="utf-8")):
        anchor = getattr(event, "anchor", None)
        if isinstance(event, yaml.AliasEvent):
            offenders.append(f"  line {event.start_mark.line + 1}: alias *{anchor}")
        elif anchor:
            offenders.append(f"  line {event.start_mark.line + 1}: anchor &{anchor}")

    assert not offenders, (
        f"{wf.name} uses a YAML anchor or alias, which GitHub Actions does not "
        f"support. Write the value out in full instead:\n" + "\n".join(offenders)
    )


@pytest.mark.parametrize("wf", WORKFLOWS, ids=lambda p: p.name)
def test_actions_are_pinned_to_a_major_version(wf):
    """An unpinned `uses:` silently changes what CI runs."""
    unpinned = [
        step["uses"]
        for job in _load(wf).get("jobs", {}).values()
        for step in job.get("steps", [])
        if isinstance(step, dict) and "uses" in step and "@" not in step["uses"]
    ]
    assert not unpinned, f"{wf.name} uses unpinned actions: {unpinned}"
