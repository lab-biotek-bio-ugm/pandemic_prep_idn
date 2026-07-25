# Environment locks

Each `envs/*.yaml` declares dependencies with version *ranges* (`>=`), so solving
it produces different exact packages as the conda channels move. For a pipeline
whose output goes into a paper that is a reproducibility hole: the run that
generated a figure cannot be recreated a year later.

The `*.conda-lock.yml` files here are the fix. Each pins one env to exact builds
for `linux-64` (the HPC target), produced by [conda-lock](https://github.com/conda/conda-lock).
They are committed on purpose: a lock is a *source of truth* for reproducibility,
not a derived build artefact, and the pipeline's provenance depends on it.

## Install an env from its lock (reproducible)

```bash
conda-lock install --name pandemic_nextstrain envs/locks/nextstrain.conda-lock.yml
# or, without conda-lock:
conda create --name pandemic_nextstrain --file envs/locks/nextstrain.conda-lock.yml
```

Snakemake's `--use-conda` still reads `envs/*.yaml`. The locks are for exact
reproduction and for archiving alongside a submission; use them when a run must
match a previous one build-for-build.

## Regenerate after editing an env

Whenever you change an `envs/*.yaml`, re-lock it so the two do not drift:

```bash
conda-lock lock -f envs/<name>.yaml -p linux-64 \
    --lockfile envs/locks/<name>.conda-lock.yml
```

`tests/contracts/test_env_locks.py` checks that every env has a lock and that
each lock targets linux-64 and points back at its env file. It does not detect a
lock that is stale relative to its env (that needs conda-lock and a network
solve); re-lock by hand after an edit, or run the command above for all envs.
