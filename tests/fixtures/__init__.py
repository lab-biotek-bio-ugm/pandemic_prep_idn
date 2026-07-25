"""Reusable builders for synthetic pipeline inputs and fake external binaries.

The pipeline shells out to blast, seqkit, hyphy and primerlab, none of which are
present in CI. Rather than run them, the tests fake their OUTPUT: a shim binary
on PATH that writes the artefact the wrapper reads, then exits with a chosen
code. That is the only way to exercise the wrapper logic -- which is where the
correctness bugs have historically lived -- without the toolchain.
"""
import os
import stat
from pathlib import Path


def install_shim(bin_dir, name, body):
    """Write an executable Python shim `name` into bin_dir and return bin_dir.

    Prepend the returned directory to PATH so the wrapper picks up the shim
    instead of the real binary.
    """
    bin_dir = Path(bin_dir)
    bin_dir.mkdir(parents=True, exist_ok=True)
    shim = bin_dir / name
    shim.write_text("#!/usr/bin/env python3\n" + body, encoding="utf-8")
    shim.chmod(shim.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return bin_dir


def path_with(bin_dir):
    """A PATH string with bin_dir first."""
    return f"{bin_dir}{os.pathsep}{os.environ['PATH']}"
