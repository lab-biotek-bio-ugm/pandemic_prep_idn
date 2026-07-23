#!/usr/bin/env python3
"""
Centralised subprocess execution for pipeline scripts.

Purpose
-------
Across the pipeline, external tools (seqkit, hyphy, blast/primerlab, makeblastdb,
NCBI datasets) were each wrapped in ad-hoc ``try/except Exception`` blocks that
returned an empty/zero/neutral default on failure. That "fail-open" pattern makes a
tool that *never ran* (missing binary, crash, bad args) indistinguishable from a
tool that *ran and legitimately found nothing* — e.g. a seqkit crash reads as "0%
amplification" and a broken BLAST reads as "no off-target hits = specific".

This module provides ONE primitive that keeps those two cases distinct:

    - The tool could not be executed / exited non-zero  -> raise ToolError (fail-closed)
    - The tool ran successfully                          -> return its CompletedProcess

Callers then interpret a *successful* run's (possibly empty) output themselves, and
let ToolError propagate so a real failure aborts the step instead of masquerading as
a valid negative result.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from typing import Optional, Sequence


class ToolError(RuntimeError):
    """An external tool was missing, crashed, or exited non-zero.

    Distinct from a successful run that simply produced empty/negative output —
    catching this should mean "the analysis did not run", never "the analysis ran
    and found nothing".
    """

    def __init__(self, tool: str, message: str, returncode: Optional[int] = None,
                 stderr: str = ""):
        self.tool = tool
        self.returncode = returncode
        self.stderr = stderr
        super().__init__(message)


def require_tool(name: str) -> str:
    """Return the resolved path to an executable, or raise ToolError if absent.

    Use at the top of a step to fail loudly and early when a dependency is not on
    PATH (e.g. the conda env was not activated), rather than deep inside a loop.
    """
    path = shutil.which(name)
    if path is None:
        raise ToolError(name, f"Required tool '{name}' not found on PATH. "
                              f"Is the correct conda environment activated?")
    return path


def run_tool(
    cmd: Sequence[str],
    *,
    tool: Optional[str] = None,
    allow_nonzero: bool = False,
    capture_output: bool = True,
    text: bool = True,
    check_binary: bool = True,
    **kwargs,
) -> subprocess.CompletedProcess:
    """Run an external command, raising ToolError on any execution failure.

    Parameters
    ----------
    cmd : the argv list (first element is the executable).
    tool : display name for errors (defaults to ``cmd[0]``).
    allow_nonzero : if True, a non-zero exit is NOT an error and is returned to the
        caller (use only when a specific non-zero code is a meaningful, documented
        result — never as a blanket "ignore failures").
    check_binary : verify the executable exists first for a clearer error.

    Returns
    -------
    subprocess.CompletedProcess on success. On success the output may still be empty
    — that is the caller's to interpret; it is NOT a failure.

    Raises
    ------
    ToolError if the binary is missing, cannot be launched, or (unless
    ``allow_nonzero``) exits non-zero.
    """
    name = tool or (cmd[0] if cmd else "<empty>")
    if check_binary and cmd:
        require_tool(cmd[0])

    try:
        result = subprocess.run(cmd, capture_output=capture_output, text=text, **kwargs)
    except FileNotFoundError as e:
        raise ToolError(name, f"Could not execute '{name}': {e}") from e
    except OSError as e:
        raise ToolError(name, f"OS error launching '{name}': {e}") from e

    if result.returncode != 0 and not allow_nonzero:
        # Read stderr whenever it was actually captured. Callers may pipe it
        # explicitly via stderr=PIPE while passing capture_output=False; keying on
        # capture_output alone discarded the tool's own error message, leaving the
        # operator with nothing but an exit code.
        stderr = result.stderr or ""
        raise ToolError(
            name,
            f"'{name}' failed with exit code {result.returncode}."
            + (f" stderr: {stderr.strip()[:1000]}" if stderr.strip() else ""),
            returncode=result.returncode,
            stderr=stderr,
        )
    return result


def warn(msg: str) -> None:
    """Emit a visible warning to stderr (kept simple; callers may use logging)."""
    print(f"[proc_utils][WARN] {msg}", file=sys.stderr)
