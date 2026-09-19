"""Thin entrypoint: ingest GSE223138 PBMC (6 donors) into interim h5ad."""
import argparse
import json
import sys
from pathlib import Path

from pdproduct.datasets.tenx import ingest_gse223138
from pdproduct.provenance import finalize_run_manifest, new_run_manifest, write_run_manifest


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-tar", required=True)
    ap.add_argument("--extract-dir", required=True)
    ap.add_argument("--meta-json", required=True)
    ap.add_argument("--out-h5ad", required=True)
    ap.add_argument("--report-json", required=True)
    ap.add_argument("--run-dir", required=True)
    args = ap.parse_args()

    manifest = new_run_manifest(
        "ingest-gse223138",
        approved_scope="phaseB-W1 data ingestion (ADR-0001)",
        command="python scripts/preprocess_gse223138.py",
    )
    write_run_manifest(manifest, args.run_dir)

    report = ingest_gse223138(Path(args.raw_tar), Path(args.extract_dir), Path(args.out_h5ad), Path(args.meta_json))
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
