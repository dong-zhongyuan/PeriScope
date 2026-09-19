"""TCR clonotype analysis (GSE253981, iRepertoire, 44 samples).

Independent donor-level test of the Ma 2025 claim direction (increased
clonal expansion in PD substantia nigra). Join key: library name (PDCxxxx.n)
between RAW.tar filenames and the 44-sample metadata xlsx.
Evidence class: C1 association, donor/sample-level only.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def parse_tcr_table(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as f:
        df = pd.read_csv(f, sep="\t")
    return df


def clonotype_metrics(df: pd.DataFrame, expanded_min_reads: int = 3) -> dict:
    reads = df["Reads"].astype(int).values
    total = reads.sum()
    if total == 0:
        return {"n_clonotypes": 0, "n_reads": 0}
    p = reads / total
    return {
        "n_clonotypes": int(len(df)),
        "n_reads": int(total),
        "top_clone_freq": float(reads.max() / total),
        "simpson_diversity": float(1.0 - np.sum(p**2)),
        "n_expanded_ge3": int((reads >= expanded_min_reads).sum()),
        "frac_expanded_ge3": float((reads >= expanded_min_reads).mean()),
    }


def load_metadata(xlsx_path: Path) -> pd.DataFrame:
    df = pd.read_excel(xlsx_path, header=None)
    data = df[~df[0].astype(str).str.startswith("#")].dropna(how="all")
    hdr = data.iloc[1].tolist()
    rows = data.iloc[2:]
    cols = {c: i for i, c in enumerate(hdr)}
    sub = rows[[cols["library name"], cols["title"], cols["Condition"], cols["tissue"]]].dropna()
    sub.columns = ["library", "donor", "condition", "tissue"]
    return sub


def analyze(extract_dir: Path, xlsx_path: Path) -> dict:
    meta = load_metadata(xlsx_path)
    lib2meta = {str(r.library): r for r in meta.itertuples()}

    rows = []
    for tsv in sorted(extract_dir.glob("GSM*_PDC*.filtJS.tsv.gz")):
        lib = tsv.name.split("_", 1)[1].replace(".filtJS.tsv.gz", "")
        if lib not in lib2meta:
            continue
        m = lib2meta[lib]
        df = parse_tcr_table(tsv)
        rows.append({"library": lib, "donor": str(m.donor), "condition": str(m.condition),
                     "tissue": str(m.tissue), **clonotype_metrics(df)})
    res = pd.DataFrame(rows)
    res.to_csv(extract_dir.parent / "tcr_per_sample_metrics.csv", index=False)

    sn = res[res["tissue"] == "Substantia Nigra"]
    out = {"n_samples_total": int(len(res)), "tissue_counts": res["tissue"].value_counts().to_dict()}
    tests = {}
    for metric in ["top_clone_freq", "simpson_diversity", "frac_expanded_ge3", "n_clonotypes"]:
        a = sn.loc[sn.condition == "Control", metric].values
        b = sn.loc[sn.condition == "PD", metric].values
        u, p = stats.mannwhitneyu(b, a, alternative="greater")
        tests[metric] = {
            "PD_mean": float(np.mean(b)), "Ctrl_mean": float(np.mean(a)),
            "n_PD": int(len(b)), "n_Ctrl": int(len(a)),
            "mw_p_one_sided_PD_greater": float(p),
        }
    out["sn_tests"] = tests
    return out
