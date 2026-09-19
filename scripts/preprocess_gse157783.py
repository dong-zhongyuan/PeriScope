"""Thin entrypoint: ingest GSE157783 (Smajic 2022 midbrain snRNA, PD vs control) into interim h5ad."""
import argparse
import json
import sys
from pathlib import Path

from pdproduct.datasets.tenx import ingest_gse157783
from pdproduct.provenance import finalize_run_manifest, new_run_manifest, write_run_manifest


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--umi-tsv", required=True)
    ap.add_argument("--cell-tsv", required=True)
    ap.add_argument("--genes-tsv", required=True)
    ap.add_argument("--hgnc-map", default=None,
                    help="genenames.org custom tsv (ensg|symbol|prev) to attach approved symbols")
    ap.add_argument("--out-h5ad", required=True)
    ap.add_argument("--report-json", required=True)
    ap.add_argument("--run-dir", required=True)
    args = ap.parse_args()

    manifest = new_run_manifest(
        "ingest-gse157783",
        approved_scope="phaseB third-cohort external validation (MOBP oligodendrocyte, HANDOFF §20.2)",
        command="python scripts/preprocess_gse157783.py",
    )
    write_run_manifest(manifest, args.run_dir)

    report = ingest_gse157783(
        Path(args.umi_tsv), Path(args.cell_tsv), Path(args.genes_tsv), Path(args.out_h5ad),
        hgnc_map_tsv=Path(args.hgnc_map) if args.hgnc_map else None,
    )
    with open(args.report_json, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False, default=str)

    finalize_run_manifest(
        manifest, status="completed", exit_code=0,
        output_paths=[args.out_h5ad, args.report_json], run_dir=args.run_dir,
    )
    summary = {k: report[k] for k in (
        "n_cells", "n_genes", "n_donors", "condition_donor_counts", "out_h5ad_bytes",
    ) if k in report}
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
