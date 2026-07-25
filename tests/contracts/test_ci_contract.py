"""Checks on the workflow files themselves.

The dormancy guard is the one piece of CI configuration with a cost attached to
getting it wrong. The upstream lab repository is private, and a workflow running
there spends the lab's Actions quota -- including for pull requests, which bill
to the base repository rather than the fork. GitHub Actions has no
workflow-level `if`, so the guard has to be repeated on every job, and a job
added later is exactly the kind of thing that gets it wrong.
"""
import pytest
import yaml

from conftest import REPO

pytestmark = pytest.mark.contract

WORKFLOW_DIR = REPO / ".github" / "workflows"
WORKFLOWS = sorted(p for p in WORKFLOW_DIR.glob("*.y*ml"))
GUARD = "github.repository == 'engkinandatama/pandemic-preparedness'"


def _load(path):
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def test_there_is_at_least_one_workflow():
    assert WORKFLOWS, f"no workflow files found under {WORKFLOW_DIR}"


@pytest.mark.parametrize("wf", WORKFLOWS, ids=lambda p: p.name)
def test_every_job_carries_the_dormancy_guard(wf):
    jobs = _load(wf).get("jobs", {})
    assert jobs, f"{wf.name} defines no jobs"

    unguarded = [name for name, job in jobs.items() if GUARD not in str(job.get("if", ""))]
    assert not unguarded, (
        f"{wf.name}: these jobs would also run on the upstream lab repository and "
        f"spend its Actions quota: {unguarded}\n"
        f"Add to each one:\n    if: {GUARD}"
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
