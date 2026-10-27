"""MG-AP / astro: every benchmarked method under one consistent aggregation rule (equal weight over the four
blood routes, then over the available seeds) - the same rule the pig evaluation uses as its primary endpoint.
usage: python mgap_allmethods.py OUTDIR
"""
import sys, os, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets/results/benchmark_20261002"; OUT = sys.argv[1]
pd.set_option("display.width", 250)
routes = ["cDC", "classical_mono", "nonclassical_mono", "pDC"]; receivers = ["microglia_mhc2", "astro"]
methods = ["periscope", "shared_label_mean", "shared_label_marginal", "ot_ridge", "cellot", "conditional_icnn", "conditional_mlp"]
obs = {}
for rc in receivers:
    k = np.load(A + f"/cache/cDC__{rc}.npz", allow_pickle=True); obs[rc] = k["profiles1"].mean(0) - k["profiles0"].mean(0)

rows = []
for rc in receivers:
    print(f"\n=== {rc} ===")
    for mth in methods:
        seeds = [s for s in range(42, 52) if os.path.exists(A + f"/results/seed_{s}/{mth}__cDC__{rc}.npz")]
        per_seed_agg, per_route = [], {}
        for s in seeds:
            ds = []
            for rt in routes:
                f = A + f"/results/seed_{s}/{mth}__{rt}__{rc}.npz"
                if not os.path.exists(f): continue
                z = np.load(f); d = z["map1_blood1"].mean(0) - z["map0_blood0"].mean(0); ds.append(d)
                per_route.setdefault(rt, []).append(d)
            if len(ds) == len(routes): per_seed_agg.append(np.mean(ds, 0))
        if not per_seed_agg: continue
        ens = np.mean(per_seed_agg, 0)
        rho_ens = stats.spearmanr(ens, obs[rc]).correlation if np.ptp(ens) > 1e-9 else np.nan
        rs = [stats.spearmanr(d, obs[rc]).correlation for d in per_seed_agg if np.ptp(d) > 1e-9]
        rroutes = {rt: stats.spearmanr(np.mean(v, 0), obs[rc]).correlation for rt, v in per_route.items() if np.ptp(np.mean(v, 0)) > 1e-9}
        print(f"  {mth:22s} seeds={len(seeds):2d} | equal-route ensemble rho {rho_ens:+.3f} | per-seed mean {np.mean(rs):+.3f} ({np.mean(np.array(rs)>0)*100:.0f}% positive) | per-route {' '.join(f'{k}:{v:+.2f}' for k,v in rroutes.items())}")
        rows.append(dict(receiver=rc, method=mth, n_seeds=len(seeds), rho_equal_route_ensemble=rho_ens, rho_per_seed_mean=np.mean(rs), frac_positive=np.mean(np.array(rs) > 0), **{f"rho_{k}": v for k, v in rroutes.items()}))
    # both receivers combined, concatenated contrast
print("\n=== both receivers pooled (concatenated gene contrasts, equal-route ten-seed ensemble) ===")
for mth in methods:
    cat_p, cat_o = [], []
    for rc in receivers:
        seeds = [s for s in range(42, 52) if os.path.exists(A + f"/results/seed_{s}/{mth}__cDC__{rc}.npz")]
        agg = []
        for s in seeds:
            ds = [np.load(A + f"/results/seed_{s}/{mth}__{rt}__{rc}.npz")["map1_blood1"].mean(0) - np.load(A + f"/results/seed_{s}/{mth}__{rt}__{rc}.npz")["map0_blood0"].mean(0) for rt in routes if os.path.exists(A + f"/results/seed_{s}/{mth}__{rt}__{rc}.npz")]
            if len(ds) == len(routes): agg.append(np.mean(ds, 0))
        if not agg: continue
        cat_p.append(np.mean(agg, 0)); cat_o.append(obs[rc])
    if len(cat_p) == 2:
        p = np.concatenate(cat_p); o = np.concatenate(cat_o)
        print(f"  {mth:22s} pooled rho {stats.spearmanr(p,o).correlation:+.3f}")
pd.DataFrame(rows).to_csv(OUT + "/mgap_allmethods.csv", index=False)
print("\nwrote", OUT + "/mgap_allmethods.csv")
