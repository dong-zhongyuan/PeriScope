"""Thin entrypoint: make the donor-level split for GSE178265 SN donors.

Locked-test donors are reserved once; the split file + sha256 sidecar go to
the asset root and are referenced by run manifests from W3 onward.
"""
import argparse
import json
import sys
from pathlib import Path

from pdproduct.datasets.splits import make_donor_split, save_split
from pdproduct.provenance import finalize_run_manifest, new_run_manifest, write_run_manifest


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc-json", required=True)
    ap.add_argument("--out-split", required=True)
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    manifest = new_run_manifest(
        "split-gse178265-v1",
        approved_scope="phaseB-W1 donor split pre-registration (ADR-0001)",
        command="python scripts/make_split_gse178265.py",
        seed=args.seed,
    )
    write_run_manifest(manifest, args.run_dir)

    qc = json.load(open(args.qc_json, encoding="utf-8"))
    # SN donors only: the caudate (CN-) donors are a different tissue and stay
    # out of the SN discovery split (recorded as excluded).
    sn_facts = {d: f for d, f in qc["donors"].items() if d.startswith("SN-")}
    excluded = sorted(d for d in qc["donors"] if d.startswith("CN-"))

    split = make_donor_split(
        sn_facts,
        ratios={"train": 0.55, "calibration": 0.15, "model_selection": 0.15, "locked_test": 0.15},
        seed=args.seed,
    )
    split["meta"]["dataset"] = "GSE178265"
    split["meta"]["scope"] = "SN donors only; caudate donors excluded"
    split["meta"]["excluded_donors"] = excluded
    path = save_split(split, args.out_split)

    finalize_run_manifest(
        manifest, status="completed", exit_code=0,
        output_paths=[path, str(path).replace(".json", ".sha256")],
        run_dir=args.run_dir,
    )
    counts = {}
    for d, p in split["donors"].items():
        key = (p, sn_facts[d]["condition"])
        counts[key] = counts.get(key, 0) + 1
    print(json.dumps({"split_file": str(path), "partition_x_condition": {f"{k[0]}/{k[1]}": v for k, v in sorted(counts.items())}, "excluded": excluded}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
