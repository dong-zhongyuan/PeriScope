"""MG-AP (and astro) benchmark: per-gene contrast decomposition.
(a) reliability ceiling of the observed test-donor contrast (donor jackknife split correlation);
(b) per-route and equal-route, ten-seed aggregate direction correlation - the same aggregation rule the
    pig evaluation uses as its primary - versus the per-route single-seed numbers now reported;
(c) where the label-mean baseline gets its advantage.
usage: python mgap_routes.py OUTDIR
"""
import sys, json, itertools, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets/results/benchmark_20261002"; OUT = sys.argv[1]
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
routes = ["cDC", "classical_mono", "nonclassical_mono", "pDC"]; receivers = ["microglia_mhc2", "astro"]
seeds = list(range(42, 52))

import os
print("=== available result dirs ===", sorted(os.listdir(A + "/results")))
print("=== methods present per seed dir ===")
for s in seeds:
    d = A + f"/results/seed_{s}"
    if os.path.isdir(d): print(" seed", s, sorted({f.split("__")[0] for f in os.listdir(d)}))

obs = {}; rel = {}
for rc in receivers:
    k = np.load(A + f"/cache/cDC__{rc}.npz", allow_pickle=True)
    p0, p1 = k["profiles0"], k["profiles1"]
    obs[rc] = p1.mean(0) - p0.mean(0)
    # donor jackknife: correlation between contrasts computed from disjoint-as-possible donor subsets
    pairs = []
    for i in range(len(p0)):
        for j in range(len(p1)):
            c_in = p1[j] - p0[i]
            oi = [x for x in range(len(p0)) if x != i]; oj = [x for x in range(len(p1)) if x != j]
            if not oi or not oj: continue
            c_out = p1[oj].mean(0) - p0[oi].mean(0)
            pairs.append(stats.spearmanr(c_in, c_out).correlation)
    rel[rc] = (np.mean(pairs) if pairs else np.nan, len(pairs), len(p0), len(p1))
    print(f"\n{rc}: {len(p0)} control / {len(p1)} PD test donors | observed contrast RMS {np.sqrt((obs[rc]**2).mean()):.4f}")
    print(f"  donor split-half reliability of the observed contrast: mean Spearman {rel[rc][0]:+.3f} over {rel[rc][1]} splits")

rows = []
for rc in receivers:
    ens = {}
    for rt in routes:
        per_seed = []
        for s in seeds:
            f = A + f"/results/seed_{s}/periscope__{rt}__{rc}.npz"
            if not os.path.exists(f): continue
            z = np.load(f); d = z["map1_blood1"].mean(0) - z["map0_blood0"].mean(0)
            per_seed.append(d)
        if not per_seed: continue
        per_seed = np.array(per_seed)
        rhos = [stats.spearmanr(d, obs[rc]).correlation for d in per_seed]
        ensd = per_seed.mean(0); ens[rt] = ensd
        rows.append(dict(receiver=rc, route=rt, n_seeds=len(per_seed), rho_seed_mean=np.mean(rhos), rho_seed_min=np.min(rhos), rho_seed_max=np.max(rhos), rho_ten_seed_ensemble=stats.spearmanr(ensd, obs[rc]).correlation))
        print(f"{rc:16s} {rt:20s} per-seed rho mean {np.mean(rhos):+.3f} (range {np.min(rhos):+.3f} to {np.max(rhos):+.3f}) | ten-seed ensemble {stats.spearmanr(ensd, obs[rc]).correlation:+.3f}")
    if len(ens) == len(routes):
        agg = np.mean([ens[r] for r in routes], 0)
        r_agg = stats.spearmanr(agg, obs[rc]).correlation
        # equal-route aggregate per seed, then across seeds
        per_seed_agg = []
        for s in seeds:
            ds = []
            for rt in routes:
                f = A + f"/results/seed_{s}/periscope__{rt}__{rc}.npz"
                if os.path.exists(f):
                    z = np.load(f); ds.append(z["map1_blood1"].mean(0) - z["map0_blood0"].mean(0))
            if len(ds) == len(routes): per_seed_agg.append(stats.spearmanr(np.mean(ds, 0), obs[rc]).correlation)
        print(f"  >>> {rc}: EQUAL-ROUTE ten-seed ensemble rho = {r_agg:+.3f} | equal-route single-seed rho mean {np.mean(per_seed_agg):+.3f} (range {np.min(per_seed_agg):+.3f} to {np.max(per_seed_agg):+.3f}, {np.mean(np.array(per_seed_agg)>0):.0%} positive)")
        rows.append(dict(receiver=rc, route="EQUAL-ROUTE AGGREGATE", n_seeds=len(per_seed_agg), rho_seed_mean=np.mean(per_seed_agg), rho_seed_min=np.min(per_seed_agg), rho_seed_max=np.max(per_seed_agg), rho_ten_seed_ensemble=r_agg))
        print(f"       reliability-adjusted (rho / donor split-half {rel[rc][0]:+.3f}) = {r_agg/rel[rc][0]:+.3f}")
pd.DataFrame(rows).to_csv(OUT + "/mgap_routes.csv", index=False)
print("\nwrote", OUT + "/mgap_routes.csv")
