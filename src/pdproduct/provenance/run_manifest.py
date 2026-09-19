"""Run manifest creation: every execution leaves a replayable record.

Implements the run_manifest contract (contracts/run_manifest.schema.json).
A run directory without a schema-valid manifest is not a completed run.
"""
from __future__ import annotations

import datetime as _dt
import json
import os
import platform
import subprocess
from pathlib import Path

from ..contracts.loader import repo_root, sha256_file, validate

NO_GIT_COMMIT = "NO_GIT"


def _now_iso() -> str:
    return _dt.datetime.now(_dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def git_state() -> tuple[str, bool]:
    """Return (HEAD commit, dirty flag), or an explicit no-Git sentinel."""
    try:
        repo = str(repo_root())
        commit_result = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=repo, capture_output=True, text=True, check=True
        )
        porcelain_result = subprocess.run(
            ["git", "status", "--porcelain"], cwd=repo, capture_output=True, text=True, check=True
        )
    except (OSError, subprocess.SubprocessError, RuntimeError):
        return NO_GIT_COMMIT, False
    return commit_result.stdout.strip(), bool(porcelain_result.stdout.strip())


def _env_digest() -> str | None:
    """SHA256 of the frozen environment lock (G0: environment into the manifest)."""
    import hashlib
    try:
        lock = repo_root() / "env" / "environment.lock.txt"
    except RuntimeError:
        return None
    if not lock.exists():
        return None
    return hashlib.sha256(lock.read_bytes()).hexdigest()


def new_run_manifest(
    run_id: str,
    *,
    approved_scope: str,
    command: str,
    operator: str | None = None,
    parent_run_id: str | None = None,
    seed: int | None = None,
    environment_digest: str | None = None,
) -> dict:
    commit, dirty = git_state()
    manifest = {
        "run_id": run_id,
        "parent_run_id": parent_run_id,
        "operator": operator or os.environ.get("USER") or os.environ.get("USERNAME") or "unknown",
        "approved_scope": approved_scope,
        "started_at": _now_iso(),
        "ended_at": None,
        "host": platform.node(),
        "container": None,
        "gpu": None,
        "git_commit": commit,
        "git_dirty": dirty,
        "config_sha256": None,
        "dataset_manifest_sha256": None,
        "split_sha256": None,
        "environment_digest": environment_digest or _env_digest(),
        "command": command,
        "seed": seed,
        "exit_code": None,
        "peak_ram_gb": None,
        "peak_gpu_mem_mib": None,
        "outputs": [],
        "warnings": (["git metadata unavailable"] if commit == NO_GIT_COMMIT else []),
        "refusal": None,
        "human_review": None,
        "status": "not-started",
    }
    validate(manifest, "run_manifest")
    return manifest


def write_run_manifest(manifest: dict, run_dir: str | Path) -> Path:
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    path = run_dir / "run_manifest.json"
    if path.exists():
        raise FileExistsError(f"run manifest already exists, overwriting forbidden: {path}")
    validate(manifest, "run_manifest")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False, sort_keys=True)
    return path


def finalize_run_manifest(
    manifest: dict,
    *,
    status: str,
    exit_code: int | None = None,
    output_paths: list[str | Path] | None = None,
    run_dir: str | Path | None = None,
    warnings: list[str] | None = None,
    refusal: dict | None = None,
) -> dict:
    """Close out a run: stamp outputs with SHA256, set final status, persist."""
    manifest["ended_at"] = _now_iso()
    manifest["status"] = status
    manifest["exit_code"] = exit_code
    if warnings:
        manifest["warnings"] = list(manifest.get("warnings", [])) + list(warnings)
    if refusal is not None:
        manifest["refusal"] = refusal
    for p in output_paths or []:
        path = Path(p)
        if not path.is_file():
            raise FileNotFoundError(f"run output does not exist: {path}")
        if any(entry.get("path") == str(path) for entry in manifest["outputs"]):
            continue
        manifest["outputs"].append({
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    validate(manifest, "run_manifest")
    if run_dir is not None:
        out = Path(run_dir) / "run_manifest.json"
        with open(out, "w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=2, ensure_ascii=False, sort_keys=True)
    return manifest
