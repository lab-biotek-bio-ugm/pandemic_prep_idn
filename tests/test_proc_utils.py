"""Tests for the centralised fail-closed subprocess helper (scripts/tools/proc_utils.py).

The whole point of this helper is that a tool which never ran (missing binary,
non-zero exit) raises ToolError, while a tool that ran successfully returns its
result even if the output is empty — so callers can no longer confuse "did not run"
with "ran and found nothing".
"""
import sys

import pytest

from conftest import load_module  # noqa: E402

proc = load_module("scripts/tools/proc_utils.py", "proc_utils")


def test_successful_run_returns_completed_process():
    r = proc.run_tool([sys.executable, "-c", "print('ok')"], tool="py")
    assert r.returncode == 0
    assert "ok" in r.stdout


def test_successful_run_with_empty_output_is_not_an_error():
    # Empty output on a zero exit is a legitimate result, NOT a failure.
    r = proc.run_tool([sys.executable, "-c", ""], tool="py")
    assert r.returncode == 0
    assert r.stdout == ""


def test_nonzero_exit_raises_tool_error():
    with pytest.raises(proc.ToolError) as exc:
        proc.run_tool([sys.executable, "-c", "import sys; sys.exit(3)"], tool="py")
    assert exc.value.returncode == 3


def test_missing_binary_raises_tool_error():
    with pytest.raises(proc.ToolError):
        proc.run_tool(["definitely-not-a-real-binary-xyzzy"])


def test_allow_nonzero_returns_instead_of_raising():
    r = proc.run_tool([sys.executable, "-c", "import sys; sys.exit(7)"],
                      tool="py", allow_nonzero=True)
    assert r.returncode == 7


def test_require_tool_finds_python():
    assert proc.require_tool(sys.executable.split("/")[-1] if "/" in sys.executable else "python3")


def test_require_tool_missing_raises():
    with pytest.raises(proc.ToolError):
        proc.require_tool("definitely-not-a-real-binary-xyzzy")
