"""Generate a strict run ledger with manifest and output integrity checks."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from pdproduct.contracts.loader import repo_root, sha256_file, validate_file


def _output_checks(run_dir: Path, manifest: dict) -> tuple[str, str | None]:
    outputs = manifest.get("outputs", [])
    if not outputs:
        return "not-applicable", None
    for entry in outputs:
        raw_path = Path(entry["path"])
        if raw_path.is_absolute():
            candidates = (raw_path,)
        else:
            candidates = (run_dir / raw_path, repo_root() / raw_path)
        path = next((candidate for candidate in candidates if candidate.is_file()), candidates[0])
        if not path.is_file():
            return "missing", f"output missing: {path}"
        actual = sha256_file(path)
        if actual != entry["sha256"]:
            return "hash-mismatch", f"output hash mismatch: {path}"
        if "bytes" in entry and path.stat().st_size != entry["bytes"]:
            return "size-mismatch", f"output size mismatch: {path}"
    return "ok", None


def scan_run_directory(run_dir: Path) -> dict:
    mf = run_dir / "run_manifest.json"
    base = {"run_id": run_dir.name, "status": "unreadable", "started": "?", "commit": "?",
            "manifest_check": "missing", "outputs_check": "not-applicable", "error": None}
    if not mf.is_file():
        base["error"] = "run_manifest.json missing"
        return base
    try:
        manifest = json.loads(mf.read_text(encoding="utf-8"))
    except Exception as exc:
        base["error"] = str(exc)
        return base
    base["status"] = manifest.get("status", "?")
    started = manifest.get("started_at", "?")
    base["started"] = started[:16] if isinstance(started, str) else "?"
    commit = manifest.get("git_commit", "?")
    base["commit"] = commit[:8] if isinstance(commit, str) else "?"
    try:
        validate_file(mf, "run_manifest")
    except Exception as exc:
        base["manifest_check"] = "invalid"
        base["error"] = f"manifest: {exc}"
        return base
    base["manifest_check"] = "ok"
    check, error = _output_checks(run_dir, manifest)
    base["outputs_check"] = check
    base["error"] = error
    return base


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs-dir", required=True)
    ap.add_argument("--out-md", required=True)
    args = ap.parse_args()
    runs_dir = Path(args.runs_dir)
    rows = [scan_run_directory(run_dir) for run_dir in sorted(runs_dir.iterdir()) if run_dir.is_dir()]
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    invalid = sum(row["manifest_check"] != "ok" or row["outputs_check"] in {"missing", "hash-mismatch", "size-mismatch"} for row in rows)
    out = Path(args.out_md)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        f.write("# RUN LEDGER\n\n")
        f.write(f"总 run 数：{len(rows)}｜状态分布：{counts}｜校验失败：{invalid}\n\n")
        f.write("| run_id | status | started | commit | manifest | outputs | error |\n|---|---|---|---|---|---|---|\n")
        for row in rows:
            error = (row["error"] or "").replace("|", "\|").replace("\n", " ")
            f.write(f"| {row['run_id']} | {row['status']} | {row['started']} | {row['commit']} | {row['manifest_check']} | {row['outputs_check']} | {error} |\n")
    print(json.dumps({"n_runs": len(rows), "counts": counts, "validation_failures": invalid, "ledger": str(out)}, ensure_ascii=False))
    return 1 if invalid else 0


if __name__ == "__main__":
    sys.exit(main())
