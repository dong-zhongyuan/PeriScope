"""Thin entrypoint: ingest GSE178265 Homo matrix into interim h5ad with donor annotation.

Run under a run manifest (phaseB-W1); raw is never modified.
"""
import argparse
import json
import sys
from pathlib import Path

from pdproduct.datasets.tenx import ingest_gse178265
from pdproduct.provenance import finalize_run_manifest, new_run_manifest, write_run_manifest


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--meta-json", required=True)
    ap.add_argument("--out-h5ad", required=True)
    ap.add_argument("--report-json", required=True)
    ap.add_argument("--run-dir", required=True)
    args = ap.parse_args()

    run_id = Path(args.report_json).stem.replace("_qc", "")
    manifest = new_run_manifest(
        f"ingest-{run_id}",
        approved_scope="phaseB-W1 data ingestion (ADR-0001)",
        command="python scripts/preprocess_gse178265.py",
    )
    write_run_manifest(manifest, args.run_dir)

    report = ingest_gse178265(Path(args.raw_dir), Path(args.meta_json), Path(args.out_h5ad))
    report_path = Path(args.report_json)
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False, default=str)

    finalize_run_manifest(
        manifest,
        status="completed",
        exit_code=0,
        output_paths=[args.out_h5ad, args.report_json],
        run_dir=args.run_dir,
    )

    summary = {k: report[k] for k in (
        "n_cells", "n_genes", "n_donors_with_gsm", "n_sn_donors", "n_cn_donors",
        "sn_condition_counts", "n_unmatched_barcodes", "median_genes_per_cell",
        "median_pct_mt", "out_h5ad_bytes",
    ) if k in report}
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
