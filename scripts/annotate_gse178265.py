"""Annotate GSE178265 SN cells: pre-registered filter -> lognorm -> marker scoring.

Quality signal to check: DA_neuron fraction should be lower in Disease donors
(classic SN pathology). Output: processed h5ad + annotation report.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import scanpy as sc

from pdproduct.datasets.annotate import score_cell_types
from pdproduct.provenance import finalize_run_manifest, new_run_manifest, write_run_manifest

# pre-registered filter (frozen for this run; changes require a new run id)
FILTER = {"min_genes": 500, "max_pct_mt": 10.0}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-h5ad", required=True)
    ap.add_argument("--out-h5ad", required=True)
    ap.add_argument("--report-json", required=True)
    ap.add_argument("--run-dir", required=True)
    args = ap.parse_args()

    manifest = new_run_manifest(
        "annotate-gse178265-v1",
        approved_scope="phaseB-W3/4 PD-side annotation (ADR-0001)",
        command="python scripts/annotate_gse178265.py",
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

    ann = score_cell_types(adata)

    ct_donor = (
        adata.obs.groupby(["donor", "cell_type"]).size().unstack(fill_value=0)
    )
    # restrict to SN donors for the pathology check
    sn = ct_donor[ct_donor.index.str.startswith("SN-")]
    cond = adata.obs.groupby("donor")["condition"].first()
    da_frac = (sn.get("DA_neuron", 0).T / sn.sum(axis=1)).T
    da_by_cond = {
        c: {"mean_da_fraction": float(da_frac[[d for d in da_frac.index if cond[d] == c]].mean()),
            "donors": [d for d in da_frac.index if cond[d] == c]}
        for c in sorted(set(cond[sn.index].tolist()))
    }

    report = {
        "filter": FILTER,
        "n_cells_before": int(n_before),
        "n_cells_after": int(n_after),
        "n_per_type": {k: int(v) for k, v in ann["n_per_type"].items()},
        "markers_missing": ann["markers_missing"],
        "da_fraction_by_condition": da_by_cond,
        "per_donor_celltype_head": {d: {k: int(v) for k, v in ct_donor.loc[d].sort_values(ascending=False).head(4).items()} for d in ct_donor.index[:6]},
    }
    Path(args.out_h5ad).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out_h5ad, compression="gzip")
    with open(args.report_json, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False, default=str)

    finalize_run_manifest(manifest, status="completed", exit_code=0, output_paths=[args.out_h5ad, args.report_json], run_dir=args.run_dir)
    print(json.dumps({k: report[k] for k in ("n_cells_before", "n_cells_after", "n_per_type", "da_fraction_by_condition")}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
