"""Thin entrypoint: ingest Moquin-Beaudry 2025 PBMC (30 donors incl MSA/PSP) into interim h5ad."""
import argparse
import json
import sys
from pathlib import Path

from pdproduct.datasets.tenx import ingest_mb2025
from pdproduct.provenance import finalize_run_manifest, new_run_manifest, write_run_manifest


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--extract-dir", required=True)
    ap.add_argument("--out-h5ad", required=True)
    ap.add_argument("--report-json", required=True)
    ap.add_argument("--run-dir", required=True)
    args = ap.parse_args()

    manifest = new_run_manifest(
        "ingest-mb2025",
        approved_scope="pd_product direction-analysis prereg v2 (f64fc3e)",
        command="python scripts/preprocess_mb2025.py",
    )
    write_run_manifest(manifest, args.run_dir)

    report = ingest_mb2025(Path(args.extract_dir), Path(args.out_h5ad))
    with open(args.report_json, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False, default=str)

    finalize_run_manifest(
        manifest, status="completed", exit_code=0,
        output_paths=[args.out_h5ad, args.report_json], run_dir=args.run_dir,
    )
    summary = {k: report[k] for k in ("n_cells", "n_genes", "n_donors", "per_donor_cells", "out_h5ad_bytes") if k in report}
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
