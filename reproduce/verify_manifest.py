"""Verify the released file inventory without importing scientific dependencies."""
import hashlib
import json
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / "release_manifest.json").read_text())
    errors = []
    for name, expected in manifest["files"].items():
        path = root / name
        if not path.is_file():
            errors.append(f"missing: {name}")
            continue
        content = path.read_bytes()
        if len(content) != expected["bytes"] or hashlib.sha256(content).hexdigest() != expected["sha256"]:
            errors.append(f"changed: {name}")
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"PASS: {len(manifest['files'])} released files match their SHA-256 hashes and sizes.")


if __name__ == "__main__":
    main()
