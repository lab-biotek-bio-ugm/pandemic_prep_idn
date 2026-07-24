#!/usr/bin/env python3
"""
Per-module timing, appended to one machine-readable file.
=============================================================================
The first production run left no way to answer "how long did M10 take?". Not one
log recorded a duration, so sizing decisions -- how many candidates to evaluate,
whether the 16-worker cap is right -- had to be argued from a model of the code
instead of measurement, and the runtime cost of raising a limit could only be
expressed as a multiple rather than in minutes.

Every entry is appended as one JSON line to `results/<build>/benchmark.jsonl`,
which survives partial runs: a module that crashes still records the attempt and
its error, because a step that dies after 40 minutes is exactly the one worth
knowing about.

Usage
-----
    from tools.benchmark import timed

    with timed("m10_insilico_pcr", build, mode="degenerate", n_candidates=40):
        ...

Or around a whole script:

    from tools.benchmark import timed_main
    timed_main("m07_bridge", build)   # registers an atexit hook
"""

import atexit
import json
import os
import platform
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Optional

ENV_DISABLE = "PANDEMIC_BENCHMARK_OFF"


def _record_path(build: Optional[str]) -> Optional[Path]:
    if os.environ.get(ENV_DISABLE):
        return None
    if not build:
        return None
    return Path("results") / build / "benchmark.jsonl"


def _cpu_count() -> int:
    try:
        return len(os.sched_getaffinity(0))  # respects cgroup/slurm pinning
    except AttributeError:
        return os.cpu_count() or 0


def record(step: str, build: Optional[str], seconds: float,
           status: str = "ok", error: str = "", **fields: Any) -> None:
    """Append one timing entry. Never raises: instrumentation must not break a run."""
    path = _record_path(build)
    if path is None:
        return
    entry: Dict[str, Any] = {
        "step": step,
        "seconds": round(seconds, 2),
        "status": status,
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                     time.gmtime(time.time() - seconds)),
        "host": platform.node(),
        "cpus_available": _cpu_count(),
    }
    if error:
        entry["error"] = error[:300]
    entry.update({k: v for k, v in fields.items() if v is not None})
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(entry) + "\n")
    except Exception:
        pass


@contextmanager
def timed(step: str, build: Optional[str], **fields: Any):
    """Time a block and record it, including when it raises.

    A step that fails after forty minutes is the one most worth measuring, so the
    duration is written on the error path too and the exception is re-raised
    unchanged.
    """
    start = time.time()
    try:
        yield
    except BaseException as exc:
        record(step, build, time.time() - start, status="failed",
               error=f"{type(exc).__name__}: {exc}", **fields)
        raise
    else:
        record(step, build, time.time() - start, status="ok", **fields)


def timed_main(step: str, build: Optional[str], **fields: Any) -> None:
    """Time a whole script via an atexit hook.

    For scripts whose work is not confined to one block. Uses atexit rather than
    wrapping main() so it also fires on sys.exit().
    """
    start = time.time()
    atexit.register(lambda: record(step, build, time.time() - start, **fields))


def summarise(build: str) -> str:
    """Human-readable table of a build's recorded timings."""
    path = _record_path(build)
    if path is None or not path.exists():
        return f"No benchmark data for {build}."

    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        try:
            rows.append(json.loads(line))
        except Exception:
            continue
    if not rows:
        return f"No benchmark data for {build}."

    total = sum(r["seconds"] for r in rows if r.get("status") == "ok")
    width = max(len(r["step"]) for r in rows)
    out = [f"Timing for {build} (total {total / 60:.1f} min across {len(rows)} steps)", ""]
    out.append(f"{'step'.ljust(width)}  {'minutes':>9} {'share':>7}  detail")
    for r in sorted(rows, key=lambda x: -x["seconds"]):
        share = (r["seconds"] / total * 100) if total else 0.0
        detail = ", ".join(
            f"{k}={v}" for k, v in r.items()
            if k not in ("step", "seconds", "status", "started_utc", "host",
                         "cpus_available", "error"))
        flag = "" if r.get("status") == "ok" else f"  [{r['status']}]"
        out.append(f"{r['step'].ljust(width)}  {r['seconds'] / 60:9.2f} "
                   f"{share:6.1f}%  {detail}{flag}")
    return "\n".join(out)


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("usage: python scripts/tools/benchmark.py <build_name>")
        raise SystemExit(2)
    print(summarise(sys.argv[1]))
