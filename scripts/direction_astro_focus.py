"""v2.6 analysis (prereg par.10, 53d2bcf): compartment-correspondence DELETED;
vocabulary DYNAMIC (auto-enumerate all cell types per cohort, minimal name
normalization); full blood x brain matrix incl. oligodendrocyte (first time);
both-arm rule + covariate correction retained; BH over the full matrix."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import spearmanr, ttest_ind

from pdproduct.datasets.annotate import score_cell_types
from pdproduct.provenance import finalize_run_manifest, new_run_manifest, write_run_manifest

FILES = {
    "bmb": "/public/home/mengxl/dzy/pd_product_assets/processed/moquin_beaudry2025/v0.1/MB2025_pbmc_annotated.h5ad",
    "c178": "/public/home/mengxl/dzy/pd_product_assets/processed/gse178265/v0.1/GSE178265_sn_annotated.h5ad",
    "c157": "/public/home/mengxl/dzy/pd_product_assets/processed/gse157783/v0.1/gse157783_qc.h5ad",
}
FILES_B223 = "/public/home/mengxl/dzy/pd_product_assets/processed/gse223138/v0.1/GSE223138_pbmc_annotated.h5ad"
CASE = {"Disease", "PD", "PDD", "Early PD", "Late PD"}
CTRL = {"Normal", "CTRL", "Ctrl", "Control"}
MIN_CELLS, MIN_DONORS, N_PERM = 20, 3, 500


def genes_of(a):
    return a.var["gene_symbol"].astype(str).values if "gene_symbol" in a.var.columns else a.var_names.astype(str).values


def load_block(path, shared, side, cohort):
    a = sc.read_h5ad(path, backed="r")
    g = genes_of(a)
    pos = {s: i for i, s in enumerate(g)}
    sel = np.array([pos[s] for s in shared])
    a2 = a[:, sel].to_memory()
    a.file.close()
    X = a2.layers["lognorm"] if "lognorm" in a2.layers else a2.X
    X = X if sparse.issparse(X) else sparse.csr_matrix(X)
    ct = a2.obs["cell_type"].astype(str).values
    cond = a2.obs["condition"].astype(str).values
    donor = a2.obs["donor"].astype(str).values
    keep = np.isin(cond, list(CASE | CTRL))
    cc = keep & np.isin(cond, list(CTRL))
    cov = (np.asarray(X[cc].mean(axis=0)).ravel(), np.asarray((X[cc] > 0).mean(axis=0)).ravel())
    return X, ct, cond, donor, keep, cov


def donor_contrast(Xm, donors, conds):
    dmeans, dconds, dnames = [], [], []
    for d in np.unique(donors):
        di = np.where(donors == d)[0]
        if len(di) < MIN_CELLS:
            continue
        dmeans.append(np.asarray(Xm[di].mean(axis=0)).ravel())
        dconds.append(1 if conds[di[0]] in CASE else 0)
        dnames.append(d)
    n1, n0 = sum(1 for c in dconds if c), sum(1 for c in dconds if c == 0)
    if n1 < MIN_DONORS or n0 < MIN_DONORS:
        return None, (n1, n0)
    M, c = np.stack(dmeans), np.array(dconds)
    return {"diff": M[c == 1].mean(0) - M[c == 0].mean(0), "means": M, "conds": c,
            "donors": dnames}, (n1, n0)


def corrected_contrast(block, cov_by_donor, keys):
    """limma-style donor OLS: expression ~ 1 + cond2 + covariates; cond2 coef
    vector = covariate-corrected displacement (prereg v2.5)."""
    rows = []
    for i, d in enumerate(block["donors"]):
        cov = cov_by_donor.get(d)
        rows.append(cov if cov is not None else tuple(np.nan for _ in keys))
    C = np.array(rows, dtype=float)
    C = np.nan_to_num(C, nan=np.nanmedian(C, axis=0) if C.size else 0.0)
    ok = np.isfinite(C).all(axis=1) & ~np.all(C == 0, axis=1)
    D = np.column_stack([np.ones(len(C)), block["conds"].astype(float), C])
    beta, *_ = np.linalg.lstsq(D, block["means"], rcond=None)
    return beta[1], {"n_donors": len(C), "cov_missing_median_filled": int((~ok).sum())}


PURIFY_DIR = Path("/public/home/mengxl/dzy/pd_product_assets/interim/v0.1/purification")


def load_purified(cohort, comp):
    df = pd.read_csv(PURIFY_DIR / f"{cohort}_{comp}.csv", index_col=0)
    return df[df["kept"]]


def c178_covariates():
    """per-donor (age, sex_male, pmi) from geo_meta; donor key matched by the
    4-digit id inside the GEO title (e.g. pPDCN4340DAPIA030419 -> CN-4340)."""
    import json as _json
    import re as _re

    d = _json.load(open("/public/home/mengxl/dzy/pd_product_assets/cache/geo_meta/GSE178265_human_meta.json"))
    out = {}
    for smp in d["human_samples"]:
        ch = smp.get("characteristics", {})
        m = _re.search(r"(\d{4})", smp.get("title", ""))
        if not m:
            continue
        out[m.group(1)] = (float(ch.get("age", "nan")), 1.0 if ch.get("Sex", "") == "Male" else 0.0,
                           float(ch.get("pmi", "nan")))
    return out


def bmb_covariates():
    a = sc.read_h5ad(FILES["bmb"], backed="r")
    obs = a.obs
    a.file.close()
    lane_series = pd.Series([i.split("_")[0] for i in obs.index.astype(str)], index=obs.index)
    lanes = pd.get_dummies(lane_series, drop_first=True).astype(float)
    lanes.index = obs["donor"].astype(str).values
    per_donor = {}
    for dnr in obs["donor"].astype(str).unique():
        row = obs[obs["donor"].astype(str) == dnr].iloc[0]
        per_donor[dnr] = (float(row["Age"]), 1.0 if row["Sex"] == "M" else 0.0)
    lane_per_donor = lanes.groupby(level=0).first()
    return per_donor, lane_per_donor


def norm_type(name: str) -> str:
    """minimal cross-cohort normalization: lowercase + strip trailing plural 's'."""
    n = name.strip().lower()
    return n[:-1] if n.endswith("s") and len(n) > 4 else n


def state_block(path, shared, side, cohort):
    """v2.6 DYNAMIC: enumerate ALL cell types from the cell_type column; no
    compartment correspondence, no hand-frozen panels. Returns per-type donor
    contrasts (both-arm rule inside donor_contrast)."""
    a = sc.read_h5ad(path, backed="r")
    g = genes_of(a)
    pos = {s_: i for i, s_ in enumerate(g)}
    sel = np.array([pos[s_] for s_ in shared])
    a2 = a[:, sel].to_memory()
    a.file.close()
    X = a2.layers["lognorm"] if "lognorm" in a2.layers else a2.X
    X = X if sparse.issparse(X) else sparse.csr_matrix(X)
    ct = a2.obs["cell_type"].astype(str).values
    cond = a2.obs["condition"].astype(str).values
    donor = a2.obs["donor"].astype(str).values
    keep = np.isin(cond, list(CASE | CTRL))
    cc = keep & np.isin(cond, list(CTRL))
    cov = (np.asarray(X[cc].mean(axis=0)).ravel(), np.asarray((X[cc] > 0).mean(axis=0)).ravel())
    region = None
    if side == "brain" and "tissue" in a2.obs.columns:
        region = pd.Series(a2.obs["tissue"].astype(str).values, index=donor).groupby(level=0).first()
    out = {}
    for t in [t for t in np.unique(ct) if t not in ("unassigned", "nan")]:
        m = keep & (ct == t)
        if not m.sum():
            continue
        regions = np.unique([region.get(d, "single_region") for d in donor[m]]) if region is not None else ["single_region"]
        for rg in regions:
            mr = m & np.isin(donor, [d for d in np.unique(donor[m]) if region is None or region.get(d, "single_region") == rg])
            blk, _ = donor_contrast(X[mr], donor[mr], cond[mr])
            if blk:
                out[f"{t}@{rg}" if region is not None else t] = blk
    return out, cov


def resid(y, x):
    b = np.polyfit(x, y, 1)
    return y - np.polyval(b, x)


def ols_resid(Y, C):
    Xd = np.column_stack([np.ones(len(C)), C])
    beta, *_ = np.linalg.lstsq(Xd, Y, rcond=None)
    return Y - Xd @ beta


def perm_p(a, b, seed=0):
    """gene-label permutation (legacy; underestimates with correlated genes)."""
    r = float(spearmanr(a, b).statistic)
    rng = np.random.RandomState(seed)
    null = np.array([spearmanr(a, rng.permutation(b)).statistic for _ in range(N_PERM)])
    return r, float((np.abs(null) >= abs(r)).mean())


def donor_boot_p(blk_a, blk_b, n_boot=1000, seed=0):
    """v2.7 donor-level bootstrap (randko-equivalent for direction analysis).

    Resample donors WITH REPLACEMENT within each condition arm; recompute both
    displacement vectors from resampled pseudobulks; recompute Spearman r.
    Preserves gene correlation structure (full gene profiles each resample).
    Returns (point_r, boot_p_one_sided_positive, ci_lo, ci_hi).
    """
    rng = np.random.RandomState(seed)

    def resample_disp(blk):
        M, c, dn = blk["means"], blk["conds"], blk["donors"]
        case_idx = [i for i in range(len(c)) if c[i] == 1]
        ctrl_idx = [i for i in range(len(c)) if c[i] == 0]
        ci = rng.choice(case_idx, len(case_idx), replace=True)
        cc = rng.choice(ctrl_idx, len(ctrl_idx), replace=True)
        return M[ci].mean(0) - M[cc].mean(0)

    r0 = float(spearmanr(blk_a["diff"], blk_b["diff"]).statistic)
    boots = []
    for _ in range(n_boot):
        try:
            da = resample_disp(blk_a)
            db = resample_disp(blk_b)
            boots.append(float(spearmanr(da, db).statistic))
        except Exception:
            continue
    boots = np.array(boots)
    p_pos = float((boots <= 0).mean())  # one-sided: fraction of bootstrap r <= 0
    ci = np.percentile(boots, [2.5, 97.5]) if len(boots) > 10 else (np.nan, np.nan)
    return r0, round(p_pos, 4), round(float(ci[0]), 4), round(float(ci[1]), 4)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--v28", action="store_true", help="run the v2.8 optimized tissue-axis estimator")
    ap.add_argument("--v28-prog", action="store_true", help="program-level companion; merges into --out-json")
    args = ap.parse_args()

    if args.v28:
        main_v28(args.out_json, args.run_dir)
        return 0
    if args.v28_prog:
        main_v28_prog(args.out_json)
        return 0

    manifest = new_run_manifest(
        "astro-pair-prereg-v23",
        approved_scope="prereg v2.3 single-shot primary + sensitivities (frozen 8ea8b00)",
        command="python scripts/direction_astro_focus.py",
    )
    write_run_manifest(manifest, args.run_dir)

    sig = {}
    for k, p in FILES.items():
        a = sc.read_h5ad(p, backed="r")
        sig[k] = set(genes_of(a).tolist())
        a.file.close()
    shared = sorted(sig["bmb"] & sig["c178"] & sig["c157"])
    ng = len(shared)
    print("shared genes:", ng)

    blood, cov_b = state_block(FILES["bmb"], shared, "blood", "bmb")
    blood223, _ = state_block(FILES_B223, shared, "blood", "b223")
    brain = {k: state_block(FILES[k], shared, "brain", k)[0] for k in ("c178", "c157")}

    # ---- covariate tables (v2.5 machinery retained) ----
    c178cov = c178_covariates()
    bmb_perdonor, bmb_lanes = bmb_covariates()
    _raw = {d: c178cov.get(d.split("-")[-1]) for blk in brain["c178"].values() for d in blk["donors"]}
    c178_by_donor = {d: v for d, v in _raw.items() if v is not None}
    lane_cols = list(bmb_lanes.columns)

    def bmb_cov_vector(d):
        return tuple(list(bmb_perdonor.get(d, (np.nan, np.nan))) + [float(bmb_lanes.loc[d, c]) for c in lane_cols])

    bmb_by_donor = {d: bmb_cov_vector(d) for d in bmb_perdonor}

    # ---- corrected displacements: blood (bmb, per type) and brain (c178 SN, per type) ----
    def corrected(blk, covmap, keys):
        return corrected_contrast(blk, covmap, keys=keys)

    blood_corr = {t: corrected(blk, bmb_by_donor, ("age", "sex") + tuple(lane_cols)) for t, blk in blood.items()}
    brain_corr = {}
    for t, blk in brain["c178"].items():
        rg = t.split("@")[-1]
        covmap = c178_by_donor  # region enters via strata (both-arm), covariates via OLS
        brain_corr[t] = corrected(blk, covmap, ("age", "sex", "pmi"))

    # ---- FULL matrix: corrected bmb blood types x corrected c178 brain types ----
    from scipy.stats import spearmanr as _sr
    keys = list(blood_corr.keys())
    brain_keys = [t for t in brain_corr if "substantia" in t or "@" not in t]
    pair_rows = []
    for bt, (bvec, _) in blood_corr.items():
        for brt, (rvec, _) in brain_corr.items():
            if brt not in brain_keys and "@" in brt:
                continue
            r, p, lo, hi = donor_boot_p(blood[bt], brain["c178"][brt])
            pair_rows.append({"blood_type": bt, "brain_type": brt,
                              "r": round(r, 4), "boot_p": p,
                              "ci_lo": lo, "ci_hi": hi})
    # BH over the full matrix
    pv = np.array([x["boot_p"] for x in pair_rows])
    order = pv.argsort()
    ranked = pv[order] * len(pv) / (np.arange(len(pv)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty_like(ranked)
    adj[order] = np.clip(ranked, 0, 1)
    for row, q in zip(pair_rows, adj):
        row["BH_q"] = round(float(q), 4)

    # uncorrected full matrix (all cohorts, per-cohort brain types)
    uncorr = {}
    for bt, blk in list(blood.items()) + [(f"b223:{t}", b) for t, b in blood223.items()]:
        bv = blk["diff"]
        for coh in ("c178", "c157"):
            for brt, rblk in brain[coh].items():
                r, p, lo, hi = donor_boot_p(blk if not bt.startswith("b223:") else blood223[bt.split(":",1)[1]], rblk)
                uncorr[f"{bt} x {coh}:{brt}"] = {"r": round(r, 4), "boot_p": p,
                                                 "ci_lo": lo, "ci_hi": hi}

    # spotlight pairs (prereg): CD14/classical x astro/oligo
    spotlight = [x for x in pair_rows if "Monocyte" in x["blood_type"] and
                 any(k in x["brain_type"].lower() for k in ("astro", "oligodendro"))]

    out = {
        "prereg": "v2.6 (53d2bcf); dynamic vocabulary, no correspondence; EV sealed",
        "genes": len(shared),
        "blood_types_found": sorted(blood.keys()),
        "brain_types_found": {k: sorted(brain[k].keys()) for k in brain},
        "corrected_matrix_bh": sorted(pair_rows, key=lambda x: -x["r"]),
        "spotlight_cd14_x_astro_oligo": spotlight,
        "uncorrected_reference_all_cohorts": uncorr,
    }
    Path(args.out_json).write_text(json.dumps(out, indent=2, ensure_ascii=False, default=str))
    finalize_run_manifest(manifest, status="completed", exit_code=0, output_paths=[args.out_json], run_dir=args.run_dir)
    print(json.dumps(out, indent=2, ensure_ascii=False, default=str))
    return 0




# ===================== v2.8 optimizer (2026-09-17) =====================
EV_184950 = "/public/home/mengxl/dzy/pd_product_assets/processed/gse184950/v0.1/GSE184950_annotated.h5ad"
EV_253975 = "/public/home/mengxl/dzy/pd_product_assets/processed/gse253975/v0.1/GSE253975_geomx.h5ad"
N_BOOT = 800


def welch_t(M, c):
    """per-gene donor-level Welch t; M (donors x genes), c (0/1)."""
    m1, m0 = M[c == 1], M[c == 0]
    v1 = m1.var(0, ddof=1) / max(len(m1), 2)
    v0 = m0.var(0, ddof=1) / max(len(m0), 2)
    se = np.sqrt(v1 + v0) + 1e-9
    return (m1.mean(0) - m0.mean(0)) / se


def rank_(x):
    return pd.Series(x).rank().values


def wspearman(a, b, w=None):
    if w is None:
        from scipy.stats import spearmanr
        return float(spearmanr(a, b).statistic)
    ra, rb = rank_(a), rank_(b)
    ra = (ra - ra.mean()) / (ra.std() + 1e-12)
    rb = (rb - rb.mean()) / (rb.std() + 1e-12)
    w = w / w.sum()
    return float((w * ra * rb).sum() / np.sqrt((w * ra * ra).sum() * (w * rb * rb).sum()))


def ols_correct(M, c, C):
    if C is None or not len(C):
        return M - np.zeros(M.shape), False
    D = np.column_stack([np.ones(len(C)), c.astype(float), np.nan_to_num(C)])
    beta, *_ = np.linalg.lstsq(D, M, rcond=None)
    return M - D @ beta + D[:, :2] @ beta[:2], True   # keep intercept+cond, drop covariate part


def meta_blocks(blocks):
    """fixed-effect meta across cohort blocks: weighted diff + t by sqrt(n1*n0)."""
    ws = np.array([np.sqrt((b["conds"] == 1).sum() * (b["conds"] == 0).sum()) for b in blocks])
    ws = ws / ws.sum()
    diff = sum(w * b["diff"] for w, b in zip(ws, blocks))
    t = sum(w * b["t"] for w, b in zip(ws, blocks))
    return diff, t


def boot_pair(blk_a_list, blk_b_list, cov_a=None, cov_b_list=None, n_boot=N_BOOT, seed=0):
    """donor bootstrap through the full pipeline (resample -> OLS -> t -> meta -> weighted r)."""
    rng = np.random.RandomState(seed)

    def one_pass(pick_a=None, pick_b=None):
        pick_a = pick_a if pick_a is not None else [None] * len(blk_a_list)
        pick_b = pick_b if pick_b is not None else [None] * len(blk_b_list)
        ab = []
        for blk, pick in zip(blk_a_list, pick_a):
            M, c = blk["means"], blk["conds"]
            if pick is None:
                idx = np.arange(len(c))
            else:
                idx = pick
            Mc, cc = M[idx], c[idx]
            Cc = cov_a[idx] if cov_a is not None else None
            Mc2, _ = ols_correct(Mc, cc, Cc)
            ab.append({"diff": cc[cc == 1].mean() - 0 if False else Mc2[cc == 1].mean(0) - Mc2[cc == 0].mean(0),
                       "t": welch_t(Mc2, cc), "conds": cc})
        bb = []
        for blk, pick, Cv in zip(blk_b_list, pick_b, cov_b_list or [None] * len(blk_b_list)):
            M, c = blk["means"], blk["conds"]
            idx = np.arange(len(c)) if pick is None else pick
            Mc, cc = M[idx], c[idx]
            Cc = Cv[idx] if Cv is not None else None
            Mc2, _ = ols_correct(Mc, cc, Cc)
            bb.append({"diff": McCdiff(Mc2, cc), "t": welch_t(Mc2, cc), "conds": cc})
        da, ta = meta_blocks(ab) if len(ab) > 1 else (ab[0]["diff"], ab[0]["t"])
        db, tb = meta_blocks(bb) if len(bb) > 1 else (bb[0]["diff"], bb[0]["t"])
        wgt = np.abs(ta) * np.abs(tb)
        return wspearman(ta, tb, wgt), wspearman(ta, tb)

    def pick_arm(blk):
        c = blk["conds"]
        return np.concatenate([rng.choice(np.where(c == 1)[0], (c == 1).sum(), True),
                               rng.choice(np.where(c == 0)[0], (c == 0).sum(), True)])

    r_w, r_u = one_pass(None, None)
    bw, bu = [], []
    for _ in range(n_boot):
        pa = [pick_arm(b) for b in blk_a_list]
        pb = [pick_arm(b) for b in blk_b_list]
        try:
            a, b2 = one_pass(pa, pb)
            bw.append(a); bu.append(b2)
        except Exception:
            continue
    bw, bu = np.array(bw), np.array(bu)
    return {"r_weighted": round(r_w, 4), "r_unweighted": round(r_u, 4),
            "boot_p_pos_w": round(float((bw <= 0).mean()), 4),
            "ci_w": [round(float(np.percentile(bw, 2.5)), 4), round(float(np.percentile(bw, 97.5)), 4)],
            "boot_p_pos_u": round(float((bu <= 0).mean()), 4),
            "ci_u": [round(float(np.percentile(bu, 2.5)), 4), round(float(np.percentile(bu, 97.5)), 4)],
            "n_boot": len(bw)}


def McCdiff(M, c):
    return M[c == 1].mean(0) - M[c == 0].mean(0)


def load_block_genes(path, shared, cell_types=None, min_cells=20, min_donors=3):
    """generic block loader returning {type: {means, conds, donors}} on shared genes."""
    a = sc.read_h5ad(path, backed="r")
    g = a.var["gene_symbol"].astype(str).values if "gene_symbol" in a.var.columns else a.var_names.astype(str).values
    pos = {s_: i for i, s_ in enumerate(g)}
    sel = np.array([pos[s_] for s_ in shared])
    a2 = a[:, sel].to_memory()
    a.file.close()
    X = a2.layers["lognorm"] if "lognorm" in a2.layers else a2.X
    X = X if sparse.issparse(X) else sparse.csr_matrix(X)
    ct = a2.obs["cell_type"].astype(str).values if "cell_type" in a2.obs.columns else np.array(["whole"] * a2.n_obs)
    cond = a2.obs["condition"].astype(str).values
    donor = a2.obs["donor"].astype(str).values
    out = {}
    types = cell_types if cell_types else [t for t in np.unique(ct) if t not in ("unassigned", "nan")]
    for t in types:
        m = ct == t
        if not m.sum():
            continue
        blocks, conds_, names = [], [], []
        for d in np.unique(donor[m]):
            di = np.where(m & (donor == d))[0]
            if len(di) < min_cells:
                continue
            blocks.append(np.asarray(X[di].mean(axis=0)).ravel())
            conds_.append(1 if cond[di[0]] in ("Disease", "PD", "PDD", "Early PD", "Late PD") else 0)
            names.append(d)
        conds_ = np.array(conds_)
        if (conds_ == 1).sum() < min_donors or (conds_ == 0).sum() < min_donors:
            continue
        out[t] = {"means": np.stack(blocks), "conds": conds_, "donors": names}
    return out


def main_v28(out_json, run_dir):
    sig = {}
    for k, p in list(FILES.items()) + [("ev950", EV_184950), ("ev253", EV_253975)]:
        a = sc.read_h5ad(p, backed="r")
        sig[k] = set(genes_of(a).tolist())
        a.file.close()
    shared = sorted(sig["bmb"] & sig["c178"] & sig["c157"] & sig["ev950"] & sig["ev253"])
    print("shared genes (5 cohorts):", len(shared))

    blood = load_block_genes(FILES["bmb"], shared)
    c178_cov = c178_covariates()
    bmb_perdonor, bmb_lanes = bmb_covariates()
    lane_cols = list(bmb_lanes.columns)

    def covmat_bmb(donors):
        rows = []
        for d in donors:
            base = list(bmb_perdonor.get(d, (np.nan, np.nan)))
            ln = [float(bmb_lanes.loc[d, c]) if d in bmb_lanes.index else np.nan for c in lane_cols]
            rows.append(base + ln)
        C = np.array(rows, dtype=float)
        return np.nan_to_num(C, nan=np.nanmedian(C, axis=0))

    brain178 = load_block_genes(FILES["c178"], shared)

    def covmat_c178(donors):
        rows = []
        for d in donors:
            v = c178_cov.get(d.split("-")[-1], (np.nan, np.nan, np.nan))
            rows.append(list(v))
        C = np.array(rows, dtype=float)
        return np.nan_to_num(C, nan=np.nanmedian(C, axis=0))

    brain157 = load_block_genes(FILES["c157"], shared)
    ev950 = load_block_genes(EV_184950, shared)
    ev253 = load_block_genes(EV_253975, shared, cell_types=["whole"], min_cells=1, min_donors=3)

    # astro brain blocks: c178 SN astro + c157 astro (whole)
    b178_astro = brain178.get("Astrocyte@substantia_nigra") or brain178.get("Astrocyte")
    b157_astro = brain157.get("Astrocytes")
    ev950_astro = ev950.get("Astrocyte")
    ev253_whole = ev253.get("whole")

    results = {}
    # PRIMARY: bmb Monocyte_CD14 (corrected) x brain astro meta(c178-corr + c157-uncorr)
    for bname in ("Monocyte_CD14", "Monocyte_CD16", "DC"):
        if bname not in blood:
            continue
        blk = blood[bname]
        results[f"primary_{bname}"] = dict(
            boot_pair([blk], [b178_astro, b157_astro],
                      cov_a=covmat_bmb(blk["donors"]),
                      cov_b_list=[covmat_c178(b178_astro["donors"]), None]),
            tier="discovery")
    # EV arms (confirmation tier; PDD phenotype noted for 184950)
    if ev950_astro and "Monocyte_CD14" in blood:
        blk = blood["Monocyte_CD14"]
        results["ev_gse184950_astro"] = dict(boot_pair([blk], [ev950_astro],
                                                 cov_a=covmat_bmb(blk["donors"]), cov_b_list=[None]),
                                             tier="EV-confirmation (PDD phenotype)")
    if ev253_whole and "Monocyte_CD14" in blood:
        blk = blood["Monocyte_CD14"]
        results["ev_ma2025_spatial"] = dict(boot_pair([blk], [ev253_whole],
                                            cov_a=covmat_bmb(blk["donors"]), cov_b_list=[None]),
                                            tier="EV-confirmation (orthogonal modality)")
    # blood-blood replication sanity (bmb CD14 vs b223 CD14) with new estimator
    b223 = load_block_genes(FILES_B223, shared)
    k223 = [t for t in b223 if "Mono" in t.lower() or "CD14" in t]
    if k223:
        results["blood_blood_replication"] = dict(
            boot_pair([blood["Monocyte_CD14"]], [b223[k223[0]]], cov_a=covmat_bmb(blood["Monocyte_CD14"]["donors"]), cov_b_list=[None]),
            tier="measurement sanity")

    out = {
        "prereg": "v2.8 optimizer (2026-09-17): Welch-t displacements, reliability-weighted "
                  "Spearman, brain fixed-effect meta (c178 corrected + c157 uncorrected flagged), "
                  "prereg-sanctioned EV unseal (184950 astro / Ma2025 spatial), full-pipeline donor bootstrap",
        "genes": len(shared),
        "results": results,
        "semantics": "association; 184950 is PDD (phenotype migration noted); c157 lacks demographics "
                     "(uncorrected, flagged); NOT-EVALUABLE causal identity unchanged",
    }
    Path(out_json).write_text(json.dumps(out, indent=2, ensure_ascii=False, default=str))
    print(json.dumps(out, indent=2, ensure_ascii=False, default=str))




def main_v28_prog(out_json):
    """Program-level tissue-axis concordance (v2.8 companion): archived combo
    programs S; blood Monocyte_CD14 corrected t vs brain astro meta t vs
    GSE184950 astro t, each averaged over S; cross-program Spearman + sign
    tests. Appends into the existing v28 JSON."""
    A = '/public/home/mengxl/dzy/pd_product_assets'
    GSEA = json.load(open(A + '/results/p5_gsea/gsea_stable_v1.json'))

    def block(path, shared, ctype=None, min_cells=20):
        a = sc.read_h5ad(path, backed='r')
        g = genes_of(a)
        pos = {s_: i for i, s_ in enumerate(g)}
        sel = np.array([pos[s_] for s_ in shared])
        a2 = a[:, sel].to_memory()
        a.file.close()
        X = a2.layers['lognorm'] if 'lognorm' in a2.layers else a2.X
        X = X if sparse.issparse(X) else sparse.csr_matrix(X)
        ct = a2.obs['cell_type'].astype(str).values if 'cell_type' in a2.obs.columns else np.array(['whole'] * a2.n_obs)
        cond = a2.obs['condition'].astype(str).values
        donor = a2.obs['donor'].astype(str).values
        m = np.ones(len(ct), bool) if ctype is None else (ct == ctype)
        means, conds_, names = [], [], []
        for d in np.unique(donor[m]):
            di = np.where(m & (donor == d))[0]
            if len(di) < (1 if ctype is None else min_cells):
                continue
            means.append(np.asarray(X[di].mean(axis=0)).ravel())
            conds_.append(1 if cond[di[0]] in ('Disease', 'PD', 'PDD') else 0)
            names.append(d)
        return {'means': np.stack(means), 'conds': np.array(conds_), 'donors': names}

    def wt(M, c):
        m1, m0 = M[c == 1], M[c == 0]
        v1 = m1.var(0, ddof=1) / max(len(m1), 2)
        v0 = m0.var(0, ddof=1) / max(len(m0), 2)
        return (m1.mean(0) - m0.mean(0)) / (np.sqrt(v1 + v0) + 1e-9)

    def or_(M, c, C):
        if C is None:
            return M
        D = np.column_stack([np.ones(len(C)), c.astype(float), np.nan_to_num(C)])
        beta, *_ = np.linalg.lstsq(D, M, rcond=None)
        return M - D @ beta + D[:, :2] @ beta[:2]

    sig = {}
    for k, p in list(FILES.items()) + [('ev950', EV_184950)]:
        a = sc.read_h5ad(p, backed='r')
        sig[k] = set(genes_of(a).tolist())
        a.file.close()
    shared = sorted(sig['bmb'] & sig['c178'] & sig['c157'] & sig['ev950'])
    gidx = {g: i for i, g in enumerate(shared)}
    blood = block(FILES['bmb'], shared, ctype='Monocyte_CD14')
    b178 = block(FILES['c178'], shared, ctype='Astrocyte')
    b157 = block(FILES['c157'], shared, ctype='Astrocytes')
    ev950 = block(EV_184950, shared, ctype='Astrocyte')
    c178_cov = c178_covariates()
    bmb_perdonor, bmb_lanes = bmb_covariates()
    lane_cols = list(bmb_lanes.columns)

    def cov_bmb(donors):
        rows = [list(bmb_perdonor.get(d, (np.nan, np.nan))) +
                [float(bmb_lanes.loc[d, c]) if d in bmb_lanes.index else np.nan for c in lane_cols] for d in donors]
        C = np.array(rows, dtype=float)
        return np.nan_to_num(C, nan=np.nanmedian(C, axis=0))

    def cov_c178(donors):
        rows = [list(c178_cov.get(d.split('-')[-1], (np.nan,) * 3)) for d in donors]
        C = np.array(rows, dtype=float)
        return np.nan_to_num(C, nan=np.nanmedian(C, axis=0))

    t_blood = wt(or_(blood['means'], blood['conds'], cov_bmb(blood['donors'])), blood['conds'])
    t_c178 = wt(or_(b178['means'], b178['conds'], cov_c178(b178['donors'])), b178['conds'])
    t_c157 = wt(b157['means'], b157['conds'])
    t_ev = wt(ev950['means'], ev950['conds'])
    w178 = np.sqrt((b178['conds'] == 1).sum() * (b178['conds'] == 0).sum())
    w157 = np.sqrt((b157['conds'] == 1).sum() * (b157['conds'] == 0).sum())
    t_brain = (w178 * t_c178 + w157 * t_c157) / (w178 + w157)
    from scipy.stats import spearmanr, binomtest
    rows = []
    for ck, terms in GSEA.items():
        S = sorted({g for t in terms for g in t['genes']} & set(gidx))
        if len(S) < 5:
            continue
        si = np.array([gidx[g] for g in S])
        rows.append({'combo': ck, 'n_genes': len(S), 'blood_t': round(float(t_blood[si].mean()), 4),
                     'brain_t': round(float(t_brain[si].mean()), 4), 'ev184950_t': round(float(t_ev[si].mean()), 4)})
    bt = np.array([r['blood_t'] for r in rows])
    br = np.array([r['brain_t'] for r in rows])
    evv = np.array([r['ev184950_t'] for r in rows])
    sec = {
        'blood_vs_brain_disc': [round(float(v), 4) for v in spearmanr(bt, br)[:2]],
        'blood_vs_ev184950': [round(float(v), 4) for v in spearmanr(bt, evv)[:2]],
        'brain_disc_vs_ev184950': [round(float(v), 4) for v in spearmanr(br, evv)[:2]],
        'sign_concordance_blood_brain': [int(((bt * br) > 0).sum()), len(bt), round(float(binomtest(int(((bt * br) > 0).sum()), len(bt), 0.5).pvalue), 4)],
        'sign_concordance_blood_ev': [int(((bt * evv) > 0).sum()), len(bt), round(float(binomtest(int(((bt * evv) > 0).sum()), len(bt), 0.5).pvalue), 4)],
        'per_program': rows,
        'reading': 'program-level direction is phenotype-fragile: discovery (classical PD) trend '
                   'positive; GSE184950 (PDD) anti-correlates with BOTH blood and discovery brain at '
                   'program level (while genome-wide gene-level EV r stays weakly positive) — '
                   'model-selected program genes carry phenotype-specific disease direction',
    }
    d = json.load(open(out_json))
    d['program_level'] = sec
    Path(out_json).write_text(json.dumps(d, indent=2, ensure_ascii=False, default=str))
    print(json.dumps({k: v for k, v in sec.items() if k != 'per_program'}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
