"""Spot-level spatial test of the priority programs in recipient-enriched Visium spots (GSE253975,
10 donors, 5 PD / 5 control), using RCTD RNA fractions to select spots and donor-level pseudobulk.
Competitive expression-matched null (5000 sets); directional p (down = PD decrease, as in matched brain).
usage: python spot_level_test.py OUTDIR
"""
import json, sys, glob, warnings
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets"
R = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results/optimization_20261003"
OUT = sys.argv[1]; rng = np.random.RandomState(42)
pd.set_option("display.width", 260); pd.set_option("display.max_columns", 40)
programs = json.load(open(R + "/joint_analysis/programs.json"))
pr = pd.read_csv(R + "/mr_parallel_gate_20261003/priority_all_programs.csv")
six = ["raw::CD294__classical_mono__x__astro__up", "raw::CD294__pDC__x__astro__up", "raw::CD30__classical_mono__x__astro__up",
       "raw::CD80__nonclassical_mono__x__astro__up", "raw::CD80__pDC__x__astro__up", "relative::CD274__pDC__x__astro__up"]
a = ad.read_h5ad(A + "/processed/gse253975/v0.1/GSE253975_geomx.h5ad")
print("spatial:", a.shape, "| X max", float(a.X.max()), "| integer counts:", bool(np.allclose(a.X.data[:1000], np.round(a.X.data[:1000]))))
genes = np.array(a.var_names.astype(str)); uniq = sorted(set(genes)); gidx = {g: i for i, g in enumerate(uniq)}
collapse = sparse.csr_matrix((np.ones(len(genes)), (np.arange(len(genes)), [gidx[g] for g in genes])), shape=(len(genes), len(uniq)))
X = a.X @ collapse
obs = a.obs.copy(); obs["donor"] = obs["donor"].astype(str); obs["condition"] = obs["condition"].astype(str)
fr = pd.concat([pd.read_csv(f, index_col=0) for f in sorted(glob.glob(R + "/spatial_subtype_author_20261003/samples/GSM*/cell_state_RNA_fractions.csv"))])
fr = fr.reindex(obs.index); print("fractions aligned:", int(fr.notna().all(1).sum()), "of", len(obs))
fr["astro_total"] = fr["astro_fibrous"] + fr["astro_protoplasmic"]
donors = sorted(set(obs.donor)); is_pd = np.array([obs.condition[obs.donor == d].iloc[0] == "PD" for d in donors])
print("donors:", dict(zip(donors, is_pd)))
def pseudobulk(mask):
    M = []
    for d in donors:
        mm = (mask & (obs.donor == d).values)
        s = np.asarray(X[mm].sum(0)).ravel(); M.append(np.log1p(s / max(s.sum(), 1) * 1e4))
    return np.stack(M)
def welch(M):
    a_, b_ = M[is_pd], M[~is_pd]
    return (a_.mean(0) - b_.mean(0)) / np.sqrt(a_.var(0, ddof=1) / len(a_) + b_.var(0, ddof=1) / len(b_) + 1e-12)
def selections(col, q):
    sel = np.zeros(len(obs), bool); counts = {}
    for d in donors:
        m = (obs.donor == d).values; thr = np.nanquantile(fr.loc[m, col], q); pick = m & (fr[col].values >= thr); sel |= pick; counts[d] = int(pick.sum())
    return sel, counts
designs = {}
for name, col, q in [("astro_top25", "astro_total", 0.75), ("astro_top50", "astro_total", 0.5), ("fibrous_top25", "astro_fibrous", 0.75), ("mhc2_top25", "microglia_mhc2", 0.75), ("mhc2_top50", "microglia_mhc2", 0.5), ("all_spots", None, None)]:
    if col is None: sel = np.ones(len(obs), bool); counts = {d: int((obs.donor == d).sum()) for d in donors}
    else: sel, counts = selections(col, q)
    M = pseudobulk(sel); t = welch(M); detected = (M > 0).sum(0) >= 3
    designs[name] = dict(t=pd.Series(np.where(detected, t, np.nan), index=uniq), expr=pd.Series(M.mean(0), index=uniq), spots=counts)
    print(f"{name:14s} spots per donor: {min(counts.values())}-{max(counts.values())} | genes detected: {int(detected.sum())}")
def matched_p(d, gset, B=5000):
    t = d["t"].dropna(); e = d["expr"].reindex(t.index); bins = pd.qcut(e.rank(method="first"), 20, labels=False)
    gi = [g for g in gset if g in t.index]
    if len(gi) < 4 or len(gi) / len(gset) < 0.5: return dict(stat=np.nan, p_down=np.nan, p_up=np.nan, n_used=len(gi))
    stat = t[gi].mean(); pools = {b: np.where(bins.values == b)[0] for b in np.unique(bins.values)}; gb = bins[gi].values; tv = t.values
    null = np.array([np.mean([tv[rng.choice(pools[b])] for b in gb]) for _ in range(B)])
    return dict(stat=round(float(stat), 3), p_down=round((1 + np.sum(null <= stat)) / (B + 1), 4), p_up=round((1 + np.sum(null >= stat)) / (B + 1), 4), n_used=len(gi))
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); run = 1.0
    for r in range(n, 0, -1):
        i = o[r - 1]; run = min(run, p[i] * n / r); q[i] = run
    return q
rows = []
for combo in list(pr.combo) + six:
    gs = programs[combo]; bs = combo.split("__x__")[1].split("__")[0]
    primary = "astro_top25" if bs == "astro" else "mhc2_top25"; secondary = "astro_top50" if bs == "astro" else "mhc2_top50"
    r = dict(program=combo.replace("__x__", " x ").replace("raw::", "").replace("relative::", "rel::"), priority=combo in set(pr.combo), n=len(gs), design=primary)
    for lab, dn in [("primary", primary), ("secondary", secondary), ("allspots", "all_spots")]:
        res = matched_p(designs[dn], gs); r[lab + "_stat"] = res["stat"]; r[lab + "_p_down"] = res["p_down"]; r[lab + "_n"] = res["n_used"]
    rows.append(r)
df = pd.DataFrame(rows)
for lab in ["primary", "secondary", "allspots"]:
    m = df.priority.values; q = np.full(len(df), np.nan); q[m] = bh(df.loc[m, lab + "_p_down"].fillna(1).values); df[lab + "_q16"] = np.round(q, 4)
df.to_csv(OUT + "/spot_level_test.csv", index=False)
print("\n=== spot-level recipient-enriched test (p_down = PD decrease; q16 = BH over the 16 priority programs) ===")
print(df.to_string())
print("\npriority programs with primary q16<=0.05:", int((df.loc[df.priority, "primary_q16"] <= 0.05).sum()), "| secondary:", int((df.loc[df.priority, "secondary_q16"] <= 0.05).sum()), "| all spots:", int((df.loc[df.priority, "allspots_q16"] <= 0.05).sum()))
