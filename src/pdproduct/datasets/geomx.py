"""GeoMx DSP spatial data ingestion (GSE253975, ROI-level).

Counts files are gene x ROI whitespace matrices; ROI column names embed a
16bp spatial barcode (slide-position traceable via oligo map, retrofittable).
Condition is joined via the GSE253981 44-sample metadata (shared T-donor ids).
Spatial semantics: region/ROI-level measured — no spot/cell-level claims.
"""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

T_DONOR_TISSUE = ("Substantia Nigra",)


def load_t_donor_conditions(xlsx_path: Path) -> dict[str, str]:
    """T-donor -> condition from the 44-sample metadata sheet."""
    df = pd.read_excel(xlsx_path, header=None)
    data = df[~df[0].astype(str).str.startswith("#")].dropna(how="all")
    hdr = data.iloc[1].tolist()
    rows = data.iloc[2:]
    tcol, ccol, ticol = hdr.index("title"), hdr.index("Condition"), hdr.index("tissue")
    sub = rows[[tcol, ccol, ticol]].dropna()
    sub.columns = ["title", "condition", "tissue"]
    out = {}
    for r in sub.itertuples():
        t = str(r.title)
        if t.startswith("T-") and r.tissue in T_DONOR_TISSUE:
            out[t] = str(r.condition)
    return out


def ingest_gse253975(extract_dir: Path, xlsx_path: Path, out_h5ad: Path) -> dict:
    cond_map = load_t_donor_conditions(xlsx_path)
    parts = []
    donor_facts = {}
    for counts in sorted(extract_dir.glob("GSM*_T*.txt.gz")):
        gsm = counts.name.split("_")[0]
        donor = "T-" + counts.stem.split("_T")[1].split(".")[0]
        condition = cond_map.get(donor)
        if condition is None:
            raise ValueError(f"{donor} has no SN condition in 44-sample metadata")
        df = pd.read_csv(counts, sep=" ", index_col=0)
        X = ad.AnnData(
            X=sp_csr(df),
            obs=pd.DataFrame(
                {"sample_gsm": gsm, "donor": donor, "condition": condition},
                index=pd.Index([f"{donor}_{c}" for c in df.columns]),
            ),
            var=pd.DataFrame(index=pd.Index([str(g) for g in df.index])),
        )
        parts.append(X)
        donor_facts[donor] = {"condition": condition, "gsm": gsm, "n_rois": df.shape[1]}
    adata = ad.concat(parts, join="inner", index_unique=None)
    totals = np.asarray(adata.X.sum(axis=1)).ravel()
    n_genes_roi = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    adata.obs["total_counts"] = totals
    adata.obs["n_genes"] = n_genes_roi
    cond_counts = pd.Series({d: f["condition"] for d, f in donor_facts.items()}).value_counts().to_dict()
    report = {
        "n_donors": len(donor_facts),
        "n_rois_total": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "condition_donor_counts": cond_counts,
        "donors": donor_facts,
        "median_roi_total_counts": float(np.median(totals)),
        "spatial_semantics": "ROI/region-level measured (GeoMx DSP); no spot/cell-level claims",
    }
    out_h5ad = Path(out_h5ad)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad, compression="gzip")
    report["out_h5ad"] = str(out_h5ad)
    report["out_h5ad_bytes"] = out_h5ad.stat().st_size
    return report


def sp_csr(df: pd.DataFrame):
    import scipy.sparse as sp

    return sp.csr_matrix(df.T.values.astype(np.float32))
