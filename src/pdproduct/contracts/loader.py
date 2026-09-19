"""Contract schema loading and machine validation.

The schemas live in <repo>/contracts/*.schema.json and are the single
authority for dataset manifests, run manifests and artifact contracts.
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import jsonschema

SCHEMA_FILES = {
    "dataset_manifest": "dataset_manifest.schema.json",
    "run_manifest": "run_manifest.schema.json",
    "artifact": "artifact.schema.json",
    "evidence_card": "evidence_card.schema.json",
    "plan_ir": "plan_ir.schema.json",
}


def repo_root() -> Path:
    """Locate the repository root by walking up to pyproject.toml."""
    for parent in Path(__file__).resolve().parents:
        if (parent / "pyproject.toml").exists():
            return parent
    raise RuntimeError("repo root not found: no pyproject.toml in package ancestors")


def contracts_dir() -> Path:
    return repo_root() / "contracts"


@lru_cache(maxsize=None)
def load_schema(name: str) -> dict:
    if name not in SCHEMA_FILES:
        raise KeyError(f"unknown contract schema {name!r}; known: {sorted(SCHEMA_FILES)}")
    with open(contracts_dir() / SCHEMA_FILES[name], encoding="utf-8") as f:
        return json.load(f)


def check_schema_self_consistency() -> None:
    """Every schema file must itself be a valid JSON Schema (Draft 2020-12)."""
    for name in SCHEMA_FILES:
        jsonschema.Draft202012Validator.check_schema(load_schema(name))


def validate(instance: object, schema_name: str) -> True:
    """Validate an in-memory instance against a named contract schema.

    Raises jsonschema.ValidationError on failure; returns True on success.
    """
    jsonschema.validate(instance, load_schema(schema_name))
    return True


def validate_file(path: str | Path, schema_name: str) -> True:
    with open(path, encoding="utf-8") as f:
        return validate(json.load(f), schema_name)


def dataset_manifest_path(dataset_id: str) -> Path:
    """Resolve a repository-local dataset manifest across supported layouts."""
    candidates = (
        repo_root() / "manifests_datasets" / f"{dataset_id}.json",
        repo_root() / "data" / "manifests" / f"{dataset_id}.json",
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    rendered = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(f"dataset manifest {dataset_id!r} not found; checked: {rendered}")


def sha256_file(path: str | Path) -> str:
    import hashlib

    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
