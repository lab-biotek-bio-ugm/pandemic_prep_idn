"""Every conda env must have a committed lock, and each lock must match its env.

The locks (envs/locks/*.conda-lock.yml) pin each environment to exact builds for
reproducibility. This guards the cheap, per-PR failure: an env added or renamed
without a lock, a lock left behind for a deleted env, or a lock built for the
wrong platform. It does NOT detect a lock that is stale relative to an edited env
-- that needs conda-lock and a network solve, and is a manual re-lock (see
envs/locks/README.md). Only pyyaml is needed, so this runs in the contract job.
"""
import pytest
import yaml

from conftest import REPO

pytestmark = pytest.mark.contract

ENV_DIR = REPO / "envs"
LOCK_DIR = ENV_DIR / "locks"
TARGET_PLATFORM = "linux-64"

ENV_FILES = sorted(ENV_DIR.glob("*.yaml"))
LOCK_FILES = sorted(LOCK_DIR.glob("*.conda-lock.yml"))


def _lock_for(env_path):
    return LOCK_DIR / f"{env_path.stem}.conda-lock.yml"


def test_there_are_env_files():
    assert ENV_FILES, f"no env files under {ENV_DIR}"


@pytest.mark.parametrize("env_path", ENV_FILES, ids=lambda p: p.name)
def test_every_env_has_a_lock(env_path):
    lock = _lock_for(env_path)
    assert lock.exists(), (
        f"{env_path.name} has no lock at envs/locks/{lock.name}. Generate it with:\n"
        f"  conda-lock lock -f envs/{env_path.name} -p {TARGET_PLATFORM} "
        f"--lockfile envs/locks/{lock.name}"
    )


@pytest.mark.parametrize("lock_path", LOCK_FILES, ids=lambda p: p.name)
def test_lock_targets_linux64_and_points_back_at_its_env(lock_path):
    data = yaml.safe_load(lock_path.read_text(encoding="utf-8"))
    meta = data.get("metadata", {})

    assert TARGET_PLATFORM in meta.get("platforms", []), (
        f"{lock_path.name} is not locked for {TARGET_PLATFORM}: {meta.get('platforms')}"
    )

    env_stem = lock_path.name[: -len(".conda-lock.yml")]
    sources = meta.get("sources", [])
    assert any(f"{env_stem}.yaml" in str(s) for s in sources), (
        f"{lock_path.name} does not name {env_stem}.yaml in its sources ({sources}); "
        "it may have been copied from another env."
    )
    assert data.get("package"), f"{lock_path.name} lists no packages"


@pytest.mark.parametrize("lock_path", LOCK_FILES, ids=lambda p: p.name)
def test_no_orphan_lock(lock_path):
    env_stem = lock_path.name[: -len(".conda-lock.yml")]
    assert (ENV_DIR / f"{env_stem}.yaml").exists(), (
        f"{lock_path.name} has no matching envs/{env_stem}.yaml; delete the stale lock."
    )
