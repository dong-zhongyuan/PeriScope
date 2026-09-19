"""Pre-registered gene sets and donor-level set-score testing.

Sets are fixed before looking at stage-2 results (pre-registration note in
run manifest). Set score = mean log1p-CPM of member genes present in >=50%
of donors per group; tested across donors (Welch t, two-sided), BH-adjusted
across sets x cell types.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

GENE_SETS: dict[str, list[str]] = {
    "PD_risk_core": ["SNCA", "LRRK2", "GBA", "PINK1", "PARK7", "VPS35", "MAPT", "SYNJ1", "FBXO7", "PLA2G6", "ATP13A2"],
    "neuroinflammation": ["AIF1", "CX3CR1", "CCL2", "CCL3", "CCL4", "IL1B", "IL6", "TNF", "TYROBP", "APOE", "LYZ", "CTSS", "CD68", "ITGAM"],
    "lysosome_autophagy": ["LAMP1", "LAMP2", "CTSD", "CTSB", "CTSL", "GBA", "SQSTM1", "ATG5", "ATG7", "BECN1", "TFEB", "VCP"],
    "antigen_presentation": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-A", "HLA-B", "HLA-C", "B2M", "CD74"],
    "t_cell_activity": ["CD3D", "CD3E", "CD8A", "CD8B", "GZMB", "PRF1", "IFNG", "CXCL9", "CXCL10"],
}
# mitochondrial respiratory chain by objective gene-name prefixes
MITO_PREFIXES = ("NDUF", "COX", "ATP5", "UQCR", "SDH")


def build_sets(gene_names) -> dict[str, list[str]]:
    genes = [str(g) for g in gene_names]
    sets = {k: [g for g in v if g in genes] for k, v in GENE_SETS.items()}
    sets["mito_respiratory_chain"] = [g for g in genes if g.startswith(MITO_PREFIXES)]
    return {k: v for k, v in sets.items() if len(v) >= 4}


def set_scores(records, genes, condition_of, cell_type_name, ctrl_label, disease_label, gene_index=None):
    """Per-donor set scores for one cell type. Returns DataFrame donors x sets.

    gene_index: optional {symbol: column} map when genes are ENSG ids;
    built externally from var_names + gene_symbol column.
    """
    group = [(d, s) for ct, d, s in records if ct == cell_type_name]
    donors = [d for d, _ in group]
    M = np.stack([s for _, s in group]).astype(np.float64)
    cpm = M / np.clip(M.sum(axis=1, keepdims=True), 1, None) * 1e6
    logcpm = np.log1p(cpm)
    gene_list = [str(g) for g in genes]
    gidx = gene_index if gene_index is not None else {g: i for i, g in enumerate(gene_list)}
    cond = np.array([condition_of[d] for d in donors])

    sets = build_sets(list(gidx.keys()))
    out = {}
    for name, members in sets.items():
        cols = [gidx[g] for g in members if g in gidx]
        sub = logcpm[:, cols]
        present = ((M[:, cols] > 0).mean(axis=0) >= 0.5)
        if present.sum() < 4:
            continue
        out[name] = sub[:, present].mean(axis=1)
    df = pd.DataFrame(out, index=donors)
    df["__cond__"] = cond
    return df


def test_sets(df, ctrl_label, disease_label) -> pd.DataFrame:
    rows = []
    a = df["__cond__"] == ctrl_label
    b = df["__cond__"] == disease_label
    for col in [c for c in df.columns if c != "__cond__"]:
        if a.sum() < 2 or b.sum() < 2:
            continue
        t, p = stats.ttest_ind(df.loc[b, col], df.loc[a, col], equal_var=False)
        rows.append({
            "set": col,
            "mean_disease": float(df.loc[b, col].mean()),
            "mean_ctrl": float(df.loc[a, col].mean()),
            "diff": float(df.loc[b, col].mean() - df.loc[a, col].mean()),
            "p": float(p),
            "n": f"{int(b.sum())}v{int(a.sum())}",
        })
    res = pd.DataFrame(rows)
    if not len(res):
        return res
    res = res.sort_values("p")
    if len(res):
        n = len(res)
        order = res["p"].rank(method="first")
        res["q"] = (res["p"] * n / order).clip(upper=1.0)
    return res.reset_index(drop=True)
