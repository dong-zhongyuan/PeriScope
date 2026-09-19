"""10x matrix ingestion: read raw counts, attach donor/condition metadata, QC.

Ingestion is deliberately non-destructive: raw counts land in interim h5ad
unaltered; filtering happens later under pre-registered parameters only.
"""
from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

from .geo_meta import gse178265_donor_of


def _read_lines_gz(path: Path) -> list[str]:
    with gzip.open(path, "rt") as f:
        return [ln.rstrip("\n") for ln in f]


def read_10x_mtx(barcodes_tsvgz: Path, features_tsvgz: Path, matrix_mtxgz: Path) -> ad.AnnData:
    """Read one 10x mtx triplet (cells x genes, raw counts)."""
    barcodes = _read_lines_gz(barcodes_tsvgz)
    feat_lines = _read_lines_gz(features_tsvgz)
    feat = pd.DataFrame(
        [ln.split("\t") for ln in feat_lines],
        columns=["gene_id", "gene_symbol", "feature_type"],
    )
    with gzip.open(matrix_mtxgz, "rb") as f:
        m = sio.mmread(f)  # genes x cells in 10x convention
    X = sp.csr_matrix(m.T)  # -> cells x genes
    obs = pd.DataFrame(index=pd.Index(barcodes, name=None))
    var = feat.set_index("gene_id")
    adata = ad.AnnData(X=X, obs=obs, var=var)
    return adata


def basic_qc(adata: ad.AnnData) -> dict:
    """Per-cell QC metrics (computed, not applied)."""
    X = adata.X
    total_counts = np.asarray(X.sum(axis=1)).ravel()
    n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    if "gene_symbol" in adata.var.columns:
        symbols = adata.var["gene_symbol"].astype(str)
    else:
        symbols = pd.Series(adata.var_names.astype(str), index=adata.var.index)
    mito_mask = symbols.str.startswith("MT-").values | np.asarray(
        pd.Index([str(v) for v in adata.var_names]).str.startswith("MT-")
    )
    mito_counts = np.asarray(X[:, mito_mask].sum(axis=1)).ravel() if mito_mask.any() else np.zeros(X.shape[0])
    adata.obs["n_genes_by_counts"] = n_genes
    adata.obs["total_counts"] = total_counts
    adata.obs["pct_mt"] = np.where(total_counts > 0, 100.0 * mito_counts / total_counts, 0.0)
    return {
        "n_cells": int(X.shape[0]),
        "n_genes": int(X.shape[1]),
        "median_genes_per_cell": float(np.median(n_genes)),
        "median_counts_per_cell": float(np.median(total_counts)),
        "median_pct_mt": float(np.median(adata.obs["pct_mt"])),
        "n_mito_genes": int(mito_mask.sum()),
    }


def ingest_gse178265(raw_dir: Path, human_meta_json: Path, out_h5ad: Path) -> dict:
    """Ingest the merged Homo matrix; donors from barcode prefixes joined to GSM facts."""
    with open(human_meta_json, encoding="utf-8") as f:
        meta = json.load(f)

    donor_facts: dict[str, dict] = {}
    for s in meta["human_samples"]:
        donor = gse178265_donor_of(s["title"])
        if donor is None:
            continue
        ch = s["characteristics"]
        facts = donor_facts.setdefault(
            donor,
            {
                "condition": ch.get("disease", "?"),
                "sex": ch.get("Sex", ch.get("sex", "?")),
                "age": ch.get("age", "?"),
                "tissue": ("caudate" if donor.startswith("CN-") else "substantia_nigra"),
                "n_gsm": 0,
            },
        )
        facts["n_gsm"] += 1

    barcodes = _read_lines_gz(raw_dir / "GSE178265_Homo_bcd.tsv.gz")
    with gzip.open(raw_dir / "GSE178265_Homo_matrix.mtx.gz", "rb") as f:
        m = sio.mmread(f)
    X = sp.csr_matrix(m.T)
    feat_lines = _read_lines_gz(raw_dir / "GSE178265_Homo_features.tsv.gz")
    var = pd.DataFrame([ln.split("\t") for ln in feat_lines], columns=["gene_id", "gene_symbol", "feature_type"]).set_index("gene_id")

    donors = [gse178265_donor_of(bc) for bc in barcodes]
    obs = pd.DataFrame(
        {
            "donor": donors,
            "barcode": barcodes,
        },
        index=pd.Index([f"c{i}" for i in range(len(barcodes))]),
    )
    unmatched = int(obs["donor"].isna().sum())
    obs["donor"] = obs["donor"].fillna("UNMATCHED")
    obs["condition"] = obs["donor"].map(lambda d: donor_facts.get(d, {}).get("condition", "?"))
    obs["tissue"] = obs["donor"].map(lambda d: donor_facts.get(d, {}).get("tissue", "?"))
    obs["sex"] = obs["donor"].map(lambda d: donor_facts.get(d, {}).get("sex", "?"))

    adata = ad.AnnData(X=X, obs=obs, var=var)
    qc = basic_qc(adata)

    per_donor = obs.groupby("donor").size().sort_values(ascending=False)
    sn_donors = sorted(d for d in donor_facts if d.startswith("SN-"))
    cn_donors = sorted(d for d in donor_facts if d.startswith("CN-"))
    sn_cond = {d: donor_facts[d]["condition"] for d in sn_donors}
    report = {
        **qc,
        "n_unmatched_barcodes": unmatched,
        "n_donors_with_gsm": len(donor_facts),
        "n_sn_donors": len(sn_donors),
        "n_cn_donors": len(cn_donors),
        "sn_condition_counts": pd.Series(sn_cond).value_counts().to_dict(),
        "donors": donor_facts,
        "per_donor_cells_head": {k: int(v) for k, v in per_donor.head(12).items()},
    }
    out_h5ad = Path(out_h5ad)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad, compression="gzip")
    report["out_h5ad"] = str(out_h5ad)
    report["out_h5ad_bytes"] = out_h5ad.stat().st_size
    return report


# --- GSE223138 (PBMC, 6 donors: 2 Normal / 2 Early PD / 2 Advance PD) ---

_GSE223138_DONOR = {
    "Normal 1": ("BL-N1", "Normal"),
    "Normal 2": ("BL-N2", "Normal"),
    "Early PD 1": ("BL-E1", "Early PD"),
    "Early PD 2": ("BL-E2", "Early PD"),
    "Advance PD 1": ("BL-L1", "Late PD"),
    "Advance PD 2": ("BL-L2", "Late PD"),
}


def ingest_gse223138(raw_tar: Path, extract_dir: Path, out_h5ad: Path, meta_json: Path) -> dict:
    """Extract per-GSM 10x triplets, concatenate, tag donors/conditions.

    Hard fact recorded in DATA_EVIDENCE_LEDGER: only 6 donors — donor-level
    splits on this cohort alone cannot support generalization claims.
    """
    import tarfile

    extract_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(raw_tar) as tf:
        tf.extractall(extract_dir)

    import re as _re

    gsm_files: dict[str, dict[str, Path]] = {}
    for p in sorted(extract_dir.rglob("GSM*")):
        m = _re.match(r"(GSM\d+)_(\d+)_(barcodes|features|matrix)\.(tsv|mtx)\.gz", p.name)
        if not m:
            continue
        gsm, _, kind = m.group(1), m.group(2), m.group(3)
        gsm_files.setdefault(gsm, {})[kind] = p

    with open(meta_json, encoding="utf-8") as f:
        meta = json.load(f)
    gsm_to_title = {s["accession"]: s["title"].split(",")[0].strip() for s in meta["human_samples"]}

    parts = []
    donor_facts = {}
    for gsm, files in sorted(gsm_files.items()):
        a = read_10x_mtx(files["barcodes"], files["features"], files["matrix"])
        title_kind = gsm_to_title.get(gsm)
        if title_kind is None or title_kind not in _GSE223138_DONOR:
            raise ValueError(f"cannot map {gsm} to donor (title_kind={title_kind!r})")
        donor, cond = _GSE223138_DONOR[title_kind]
        a.obs["donor"] = donor
        a.obs["condition"] = cond
        donor_facts[donor] = {"condition": cond, "gsm": gsm}
        parts.append(a)

    adata = ad.concat(parts, join="outer", label="sample_batch", keys=[p.obs["donor"][0] for p in parts], index_unique="_")
    # anndata concat drops var columns; reattach from the first part (gene_id -> symbol map is stable)
    adata.var = parts[0].var.copy()
    qc = basic_qc(adata)
    per_donor = adata.obs.groupby("donor").size().to_dict()
    report = {
        **qc,
        "n_donors": len(donor_facts),
        "per_donor_cells": {k: int(v) for k, v in per_donor.items()},
        "donors": donor_facts,
        "cohort_warning": "6 donors only: no donor-level generalization claims from this cohort alone",
    }
    out_h5ad = Path(out_h5ad)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad, compression="gzip")
    report["out_h5ad"] = str(out_h5ad)
    report["out_h5ad_bytes"] = out_h5ad.stat().st_size
    return report


# --- GSE184950 (postmortem SN; Parkinson's Disease Dementia vs Unaffected Control) ---


def ingest_gse184950(raw_tar: Path, extract_dir: Path, out_h5ad: Path, meta_json: Path) -> dict:
    """Two-layer extraction (outer tar -> per-GSM tar.gz -> 10x triplet).

    donor = brain bank donor id from GEO characteristics; same-donor
    multi-sample (verified: HBHF x2, HBNX x2) is recorded at both sample
    and donor level — splits must bind by donor.
    """
    import tarfile

    extract_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(raw_tar) as tf:
        tf.extractall(extract_dir)

    with open(meta_json, encoding="utf-8") as f:
        meta = json.load(f)
    gsm_facts = {}
    for s in meta["human_samples"]:
        ch = s["characteristics"]
        donor = next((v for k, v in ch.items() if "donor" in k.lower()), None)
        cond = next((v for k, v in ch.items() if "disease" in k.lower() or "status" in k.lower()), "?")
        if donor:
            gsm_facts[s["accession"]] = {
                "donor": f"PD-{donor}",
                "condition": ("PDD" if "Parkinson" in cond else ("Control" if "Unaffected" in cond or "Control" in cond else cond)),
                "raw_condition": cond,
                "age": ch.get("age", "?"),
                "sample_title": s["title"],
            }

    import re as _re

    parts, sample_rows = [], []
    for inner in sorted(extract_dir.glob("GSM*_*.tar.gz")):
        gsm = _re.match(r"(GSM\d+)_", inner.name).group(1)
        if gsm not in gsm_facts:
            raise ValueError(f"{gsm} has no donor facts in metadata")
        sub = extract_dir / inner.stem
        with tarfile.open(inner) as tf:
            tf.extractall(sub)
        tenx_dir = sub / "filtered_feature_bc_matrix"
        a = read_10x_mtx(
            tenx_dir / "barcodes.tsv.gz", tenx_dir / "features.tsv.gz", tenx_dir / "matrix.mtx.gz"
        )
        facts = gsm_facts[gsm]
        a.obs["donor"] = facts["donor"]
        a.obs["condition"] = facts["condition"]
        a.obs["sample"] = facts["sample_title"]
        sample_rows.append({"gsm": gsm, **facts, "n_cells": a.n_obs})
        parts.append(a)

    adata = ad.concat(parts, join="outer", label="sample_batch", keys=[p.obs["sample"][0] for p in parts], index_unique="_")
    adata.var = parts[0].var.copy()
    qc = basic_qc(adata)
    donor_facts = {}
    for r in sample_rows:
        donor_facts.setdefault(r["donor"], {"condition": r["condition"], "n_samples": 0, "ages": set()})
        donor_facts[r["donor"]]["n_samples"] += 1
        donor_facts[r["donor"]]["ages"].add(str(r["age"]))
    for d in donor_facts:
        donor_facts[d]["ages"] = sorted(donor_facts[d]["ages"])
    report = {
        **qc,
        "n_samples": len(sample_rows),
        "n_donors": len(donor_facts),
        "condition_donor_counts": pd.Series({d: v["condition"] for d, v in donor_facts.items()}).value_counts().to_dict(),
        "multi_sample_donors": {d: v["n_samples"] for d, v in donor_facts.items() if v["n_samples"] > 1},
        "donors": donor_facts,
        "samples": sample_rows,
        "note": "cohort is Parkinson's Disease Dementia (PDD) vs control — use as staging/external-validation candidate, not classic-PD discovery",
    }
    out_h5ad = Path(out_h5ad)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad, compression="gzip")
    report["out_h5ad"] = str(out_h5ad)
    report["out_h5ad_bytes"] = out_h5ad.stat().st_size
    return report


def ingest_gse157783(umi_tsv: Path, cell_tsv: Path, genes_tsv: Path, out_h5ad: Path,
                     hgnc_map_tsv: Path | None = None) -> dict:
    """Ingest Smajic 2022 (GSE157783) midbrain/SN snRNA: PD vs control.

    Source layout (GEO suppl): three gzipped tars -> plain TSVs.
    UMI tsv: header row = barcodes only, data rows = integers only
    (genes x cells, no row labels; gene i = genes.tsv row-column i, verified
    by vst.mean vs per-gene total spearman == 1.0).
    cell tsv: barcode | cell_ontology | patient (PD*/C* prefix = condition).
    var_names are ENSG ids; if hgnc_map_tsv given (genenames.org custom
    download), attach approved-symbol column for cross-cohort matching.
    """
    cells = pd.read_csv(cell_tsv, sep="\t", dtype=str)
    genes = pd.read_csv(genes_tsv, sep="\t", dtype={0: str, "row": int}).rename(columns={"gene": "gene_id"})
    genes = genes.sort_values("row")  # UMI row i = genes row column i (verified: no row-label column in UMI tsv)
    # UMI tsv: header row = 41435 barcodes, data rows = integers only (genes x cells, no labels)
    df = pd.read_csv(umi_tsv, sep="\t", header=0, dtype=np.int32)
    if df.shape[1] != len(cells) or df.shape[0] != len(genes):
        raise ValueError(f"shape mismatch: UMI {df.shape} vs cells {len(cells)} x genes {len(genes)}")
    X = sp.csr_matrix(df.values.T)  # cells x genes
    obs = cells.set_index("barcode").loc[df.columns].copy()
    obs.index.name = None
    obs["donor"] = obs["patient"]
    obs["condition"] = np.where(obs["patient"].str.startswith("PD"), "PD", "Control")
    obs["cell_type"] = obs["cell_ontology"]
    var = pd.DataFrame(index=pd.Index(genes["gene_id"].values, name="gene_id"))
    n_mapped = 0
    if hgnc_map_tsv is not None:
        m = pd.read_csv(hgnc_map_tsv, sep="\t")
        m.columns = ["ensg", "symbol", "prev"]
        m = m[m["ensg"].astype(str).str.startswith("ENSG")].drop_duplicates("ensg")
        var["gene_symbol"] = var.index.map(dict(zip(m["ensg"], m["symbol"]))).fillna("")
        n_mapped = int((var["gene_symbol"] != "").sum())
    adata = ad.AnnData(X=X, obs=obs, var=var)
    qc = basic_qc(adata)
    donor_counts = obs.groupby(["condition", "donor"]).size().reset_index(name="n_cells")
    report = {
        **qc,
        "n_samples": int(obs["patient"].nunique()),
        "n_donors": int(obs["patient"].nunique()),
        "condition_donor_counts": obs.groupby("condition")["donor"].nunique().to_dict(),
        "donors": {
            r["donor"]: {"condition": r["condition"], "n_cells": int(r["n_cells"])}
            for _, r in donor_counts.iterrows()
        },
        "cell_type_counts": obs["cell_type"].value_counts().to_dict(),
        "n_symbol_mapped": n_mapped,
        "note": "Smajic 2022 Cell Reports midbrain snRNA; third-cohort external validation for oligodendrocyte MOBP (discovery: GSE178265, second: GSE184950-staging)",
    }
    out_h5ad = Path(out_h5ad)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad, compression="gzip")
    report["out_h5ad"] = str(out_h5ad)
    report["out_h5ad_bytes"] = out_h5ad.stat().st_size
    return report


def ingest_mb2025(extract_dir: Path, out_h5ad: Path) -> dict:
    """Ingest Moquin-Beaudry 2025 (Zenodo 14372434, CC-BY-4.0) Seurat export.

    Input: RDS converted via container ai-env Seurat (read-only) to MatrixMarket;
    matrix.mtx is the author LogNormalize 'data' layer (raw counts not exported).
    QC columns nFeature_RNA / nCount_RNA / percent.mt come from author metadata.
    obs schema aligned to ingest_gse223138 (donor/condition/n_genes_by_counts/
    total_counts/pct_mt); author annotations kept as extra columns.
    """
    meta = pd.read_csv(extract_dir / "metadata.csv", index_col=0)
    features = [ln.strip() for ln in open(extract_dir / "features.tsv", encoding="utf-8")]
    barcodes = [ln.strip() for ln in open(extract_dir / "barcodes.tsv", encoding="utf-8")]
    assert len(barcodes) == meta.shape[0], "barcodes/metadata row mismatch"
    m = sio.mmread(extract_dir / "matrix.mtx")  # genes x cells
    X = sp.csr_matrix(m.T)  # -> cells x genes (lognorm)
    obs = meta.loc[barcodes].copy()
    obs.index = pd.Index(obs["simpl_ID"].astype(str), dtype=str)
    obs["donor"] = obs["PaperDonor"].astype(str)
    obs["condition"] = obs["Group"].astype(str)
    obs["n_genes_by_counts"] = obs["nFeature_RNA"].astype(float)
    obs["total_counts"] = obs["nCount_RNA"].astype(float)
    obs["pct_mt"] = obs["percent.mt"].astype(float)
    keep = ["donor", "condition", "n_genes_by_counts", "total_counts", "pct_mt",
            "Donor", "Group", "Sex", "Age", "predicted.celltype.l1",
            "predicted.celltype.l2", "final_broad", "Simplified_Idents"]
    n_dup = len(features) - len(set(features))
    if n_dup:
        seen: dict[str, int] = {}
        fixed = []
        for g in features:
            if g in seen:
                seen[g] += 1
                fixed.append(f"{g}_dup{seen[g]}")
            else:
                seen[g] = 0
                fixed.append(g)
        features = fixed
    var = pd.DataFrame(index=pd.Index(features, dtype=str))
    var["gene_symbol"] = [g.split("_dup")[0] for g in features]
    adata = ad.AnnData(X=X, obs=obs[keep], var=var)
    adata.write_h5ad(out_h5ad, compression="gzip")
    per_donor = obs.groupby(["condition", "donor"]).size().to_dict()
    return {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_donors": int(obs["donor"].nunique()),
        "per_donor_cells": {f"{k[0]}_{k[1]}": int(v) for k, v in per_donor.items()},
        "duplicate_gene_symbols_suffixed": n_dup,
        "out_h5ad_bytes": int(Path(out_h5ad).stat().st_size),
    }
