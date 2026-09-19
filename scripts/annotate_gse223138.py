"""Annotate GSE223138 PBMC cells: pre-registered filter -> lognorm -> marker scoring.

Copy+edit of annotate_gse178265.py (pd_product rule 2026-08-31). Blood-side input for
the CCWM Block A stratified coupling: PBMC scRNA, 6 donors (Early PD x2 / Late PD x2 /
Normal x2). Quality signal to check: cell-type composition by condition. Output:
processed h5ad + annotation report.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import scanpy as sc

from pdproduct.datasets.annotate import score_cell_types
from pdproduct.provenance import finalize_run_manifest, new_run_manifest, write_run_manifest

# pre-registered filter (frozen for this run; changes require a new run id)
FILTER = {"min_genes": 500, "max_pct_mt": 10.0}

# PBMC marker panel (canonical PBMC conventions; shared genes like CD3D across T
# subsets are resolved by argmax of the mean-score, same convention as SN markers)
PBMC_MARKERS: dict[str, list[str]] = {
    "CD4_T": ["CD3D", "CD3E", "TRAC", "CD4", "IL7R", "CCR7"],
    "CD8_T": ["CD3D", "CD3E", "CD8A", "CD8B", "GZMK"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "TCL1A"],
    "NK": ["NKG7", "GNLY", "KLRD1", "NCR1"],
    "Monocyte_CD14": ["LYZ", "CD14", "S100A8", "S100A9", "FCN1"],
    "Monocyte_CD16": ["FCGR3A", "MS4A7", "LST1", "CDKN1C"],
    "DC": ["FCER1A", "CST3", "CLEC10A", "HLA-DQA1"],
    "pDC": ["LILRA4", "IL3RA", "CLEC4C"],
}


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-h5ad", required=True)
    ap.add_argument("--out-h5ad", required=True)
    ap.add_argument("--report-json", required=True)
    ap.add_argument("--run-dir", required=True)
    args = ap.parse_args()

    manifest = new_run_manifest(
        "annotate-gse223138-v1",
        approved_scope="pd_product PLAN step1 blood-side prep (user 2026-08-31 go)",
        command="python scripts/annotate_gse223138.py",
    )
    write_run_manifest(manifest, args.run_dir)

    adata = sc.read_h5ad(args.in_h5ad)
    n_before = adata.n_obs
    mask = (adata.obs["n_genes_by_counts"] >= FILTER["min_genes"]) & (adata.obs["pct_mt"] <= FILTER["max_pct_mt"])
    adata = adata[mask].copy()
    n_after = adata.n_obs

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()

    ann = score_cell_types(adata, markers=PBMC_MARKERS)

    ct_donor = (
        adata.obs.groupby(["donor", "cell_type"], observed=True).size().unstack(fill_value=0)
    )
    # blood-side check: cell-type composition by condition (fractions, donor-level then pooled)
    cond = adata.obs.groupby("donor", observed=True)["condition"].first()
    frac = ct_donor.div(ct_donor.sum(axis=1), axis=0)
    comp_by_cond = {
        c: {
            "mean_fraction": {k: float(v) for k, v in frac.loc[[d for d in frac.index if cond[d] == c]].mean().items() if v > 0},
            "donors": [d for d in frac.index if cond[d] == c],
        }
        for c in sorted(set(cond.tolist()))
    }

    report = {
        "filter": FILTER,
        "input_h5ad": args.in_h5ad,
        "input_sha256": sha256_file(args.in_h5ad),
        "n_cells_before": int(n_before),
        "n_cells_after": int(n_after),
        "n_per_type": {k: int(v) for k, v in ann["n_per_type"].items()},
        "markers_missing": ann["markers_missing"],
        "composition_by_condition": comp_by_cond,
        "per_donor_celltype_head": {d: {k: int(v) for k, v in ct_donor.loc[d].sort_values(ascending=False).head(4).items()} for d in ct_donor.index[:6]},
    }
    Path(args.out_h5ad).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out_h5ad, compression="gzip")
    with open(args.report_json, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False, default=str)

    finalize_run_manifest(manifest, status="completed", exit_code=0, output_paths=[args.out_h5ad, args.report_json], run_dir=args.run_dir)
    print(json.dumps({k: report[k] for k in ("n_cells_before", "n_cells_after", "n_per_type", "composition_by_condition")}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
