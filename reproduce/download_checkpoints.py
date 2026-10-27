"""Download and verify the ten authoritative PeriScope checkpoints.

Example: python reproduce/download_checkpoints.py --output /path/to/models
The output must be a new directory. No GitHub token is needed for this release.
"""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile
import tempfile
import urllib.request


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--archive", type=Path, help="verify/extract an already downloaded archive")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / "reproduce/checkpoint_release.json").read_text())
    registry = json.loads((root / manifest["registry"]).read_text())
    expected = {f"seed_{seed}/model.pt": row["sha256"] for seed, row in registry["models"].items()}
    if manifest["members"] != expected:
        raise SystemExit("Release inventory does not match the scientific model registry")
    out = args.output.expanduser().resolve()
    if out.exists():
        parser.error("output already exists; choose a new directory")
    with tempfile.TemporaryDirectory(prefix="periscope-download-") as temp:
        archive = args.archive.expanduser().resolve() if args.archive else Path(temp) / "models.tar.gz"
        if not args.archive:
            print("Downloading", manifest["url"], flush=True)
            urllib.request.urlretrieve(manifest["url"], archive)
        if archive.stat().st_size != manifest["bytes"] or digest(archive) != manifest["sha256"]:
            raise SystemExit("Archive size or SHA-256 mismatch")
        with tarfile.open(archive, "r:gz") as tar:
            members = tar.getmembers()
            if len(members) != len(expected) or {m.name for m in members} != set(expected):
                raise SystemExit("Unexpected checkpoint archive members")
            for member in members:
                if not member.isfile() or member.size > 100 * 1024 * 1024:
                    raise SystemExit(f"Invalid checkpoint member: {member.name}")
                content = tar.extractfile(member).read()
                if hashlib.sha256(content).hexdigest() != expected[member.name]:
                    raise SystemExit(f"Checkpoint hash mismatch: {member.name}")
            out.mkdir(parents=True)
            for member in members:
                path = out / member.name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(tar.extractfile(member).read())
    print(f"PASS: {len(expected)} registry-verified checkpoints installed in {out}")


if __name__ == "__main__":
    main()
