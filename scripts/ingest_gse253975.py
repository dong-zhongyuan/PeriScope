"""Thin entrypoint: ingest GSE253975 GeoMx spatial data into AnnData."""
import argparse
import json
import sys
from pathlib import Path

from pdproduct.datasets.geomx import ingest_gse253975
from pdproduct.provenance import finalize_run_manifest, new_run_manifest, write_run_manifest


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--extract-dir", required=True)
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--out-h5ad", required=True)
    ap.add_argument("--report-json", required=True)
    ap.add_argument("--run-dir", required=True)
    args = ap.parse_args()

    manifest = new_run_manifest(
        "ingest-gse253975-v1",
        approved_scope="phaseB spatial anchor ingestion (ADR-0001)",
        command="python scripts/ingest_gse253975.py",
    )
    write_run_manifest(manifest, args.run_dir)

    report = ingest_gse253975(Path(args.extract_dir), Path(args.xlsx), Path(args.out_h5ad))
    with open(args.report_json, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False, default=str)

    finalize_run_manifest(manifest, status="completed", exit_code=0, output_paths=[args.out_h5ad, args.report_json], run_dir=args.run_dir)
    print(json.dumps({k: report[k] for k in ("n_donors", "n_rois_total", "n_genes", "condition_donor_counts", "median_roi_total_counts", "out_h5ad_bytes")}, ensure_ascii=False, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
