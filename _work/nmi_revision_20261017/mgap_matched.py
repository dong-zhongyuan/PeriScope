"""MG-AP / astro: matched-seed comparison. Every method, seeds 42-44 only (the seeds every comparator has),
equal weight over the four blood routes. Also a donor-level jackknife of the observed contrast so the
comparison does not rest on one test-donor configuration.
usage: python mgap_matched.py OUTDIR
"""
import sys, os, warnings, itertools
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets/results/benchmark_20261002"; OUT = sys.argv[1]
pd.set_option("display.width", 250)
routes = ["cDC", "classical_mono", "nonclassical_mono", "pDC"]; receivers = ["microglia_mhc2", "astro"]
methods = ["periscope", "shared_label_mean", "shared_label_marginal", "ot_ridge", "cellot", "conditional_icnn", "conditional_mlp"]
LBL = {"periscope": "PeriScope", "shared_label_mean": "label/state brain mean", "shared_label_marginal": "label/state brain marginal",
       "ot_ridge": "OT-Ridge", "cellot": "CellOT", "conditional_icnn": "conditional ICNN", "conditional_mlp": "conditional MLP"}
prof = {rc: np.load(A + f"/cache/cDC__{rc}.npz", allow_pickle=True) for rc in receivers}
obs = {rc: prof[rc]["profiles1"].mean(0) - prof[rc]["profiles0"].mean(0) for rc in receivers}

def agg(mth, rc, seeds):
    out = []
    for s in seeds:
        ds = []
        for rt in routes:
            f = A + f"/results/seed_{s}/{mth}__{rt}__{rc}.npz"
            if not os.path.exists(f): return None
            z = np.load(f); ds.append(z["map1_blood1"].mean(0) - z["map0_blood0"].mean(0))
        out.append(np.mean(ds, 0))
    return np.mean(out, 0)

print("=== matched seeds 42-44, equal weight over the four blood routes ===")
rows = []
for mth in methods:
    r = {}
    for rc in receivers:
        a = agg(mth, rc, [42, 43, 44])
        r[rc] = stats.spearmanr(a, obs[rc]).correlation if a is not None and np.ptp(a) > 1e-9 else np.nan
    pooled_p = np.concatenate([agg(mth, rc, [42, 43, 44]) for rc in receivers]); pooled_o = np.concatenate([obs[rc] for rc in receivers])
    pr = stats.spearmanr(pooled_p, pooled_o).correlation
    print(f"  {LBL[mth]:28s} MG-AP {r['microglia_mhc2']:+.3f} | Astro {r['astro']:+.3f} | both receivers pooled {pr:+.3f} | lower of the two {min(r.values()):+.3f}")
    rows.append(dict(method=LBL[mth], key=mth, rho_MGAP=r["microglia_mhc2"], rho_Astro=r["astro"], rho_pooled=pr, rho_min=min(r.values())))
t = pd.DataFrame(rows).sort_values("rho_pooled", ascending=False)
print("\nranking by pooled correlation:"); print(t[["method", "rho_MGAP", "rho_Astro", "rho_pooled"]].round(3).to_string(index=False))

print("\n=== donor jackknife: drop one test brain donor at a time, recompute the observed contrast, matched seeds ===")
jl = []
for rc in receivers:
    p0, p1 = prof[rc]["profiles0"], prof[rc]["profiles1"]
    for which, idx in [("ctrl", range(len(p0))), ("pd", range(len(p1)))]:
        for i in idx:
            keep0 = [x for x in range(len(p0)) if not (which == "ctrl" and x == i)]
            keep1 = [x for x in range(len(p1)) if not (which == "pd" and x == i)]
            if not keep0 or not keep1: continue
            o = p1[keep1].mean(0) - p0[keep0].mean(0)
            for mth in methods:
                a = agg(mth, rc, [42, 43, 44])
                if a is None or np.ptp(a) < 1e-9: continue
                jl.append(dict(receiver=rc, dropped=f"{which}{i}", method=LBL[mth], rho=stats.spearmanr(a, o).correlation))
j = pd.DataFrame(jl)
pv = j.pivot_table(index="method", columns="receiver", values="rho", aggfunc=["min", "mean"])
print(pv.round(3).to_string())
print("\nfraction of donor-jackknife configurations in which each method is positive at BOTH receivers:")
for mth in methods:
    sub = j[j.method == LBL[mth]]
    print(f"  {LBL[mth]:28s} MG-AP positive in {int((sub[sub.receiver=='microglia_mhc2'].rho>0).sum())}/{len(sub[sub.receiver=='microglia_mhc2'])} | Astro positive in {int((sub[sub.receiver=='astro'].rho>0).sum())}/{len(sub[sub.receiver=='astro'])}")
t.to_csv(OUT + "/mgap_matched.csv", index=False); j.to_csv(OUT + "/mgap_jackknife.csv", index=False)
print("\nwrote", OUT + "/mgap_matched.csv and mgap_jackknife.csv")
