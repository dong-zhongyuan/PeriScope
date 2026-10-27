"""The blood-only (fixed disease condition) direction readout: is the near-zero per-route correlation a
signal-to-noise limit of a single route, or is there no blood-dependent information at all?
(a) reproduce the per-route, per-seed numbers; (b) measure how reproducible the predicted blood-only contrast
itself is across independent seed halves; (c) aggregate the four routes and the ten seeds the way the paper's
other aggregate readouts do; (d) attenuation-correct by the measured reliability of both sides.
usage: python bloodonly.py OUTDIR
"""
import sys, os, itertools, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets/results/benchmark_20261002"; OUT = sys.argv[1]
pd.set_option("display.width", 250)
routes = ["cDC", "classical_mono", "nonclassical_mono", "pDC"]; receivers = ["microglia_mhc2", "astro"]
seeds = list(range(42, 52))
cache = {rc: np.load(A + f"/cache/cDC__{rc}.npz", allow_pickle=True) for rc in receivers}
obs = {rc: cache[rc]["profiles1"].mean(0) - cache[rc]["profiles0"].mean(0) for rc in receivers}

def delta(mth, rt, rc, s, kind):
    f = A + f"/results/seed_{s}/{mth}__{rt}__{rc}.npz"
    if not os.path.exists(f): return None
    z = np.load(f)
    if kind == "fixed_control": return z["map0_blood1"].mean(0) - z["map0_blood0"].mean(0)
    if kind == "fixed_PD": return z["map1_blood1"].mean(0) - z["map1_blood0"].mean(0)
    return z["map1_blood1"].mean(0) - z["map0_blood0"].mean(0)

print("=== (a) per-route, per-seed blood-only correlation (PeriScope, condition fixed to control) ===")
per = {}
for rc in receivers:
    for rt in routes:
        rs = [stats.spearmanr(delta("periscope", rt, rc, s, "fixed_control"), obs[rc]).correlation for s in seeds]
        per[(rc, rt)] = rs
        print(f"  {rc:16s} {rt:20s} mean {np.mean(rs):+.3f}  range {np.min(rs):+.3f} to {np.max(rs):+.3f}  positive in {int(np.sum(np.array(rs)>0))}/10 seeds")
allr = np.concatenate([per[k] for k in per]); print(f"  -> 80 route x seed estimates: mean {allr.mean():+.4f}, {100*np.mean(allr>0):.0f}% positive")

print("\n=== (b) how reproducible is the predicted blood-only contrast itself? ===")
for rc in receivers:
    for kind in ["fixed_control", "condition_aware"]:
        halves = []
        for half in [seeds[:5], seeds[5:]]:
            halves.append(np.mean([np.mean([delta("periscope", rt, rc, s, kind) for rt in routes], 0) for s in half], 0))
        rel = stats.spearmanr(halves[0], halves[1]).correlation
        # Spearman-Brown step-up from two halves to the full ten-seed ensemble
        sb = 2 * rel / (1 + rel)
        amp = np.sqrt(np.mean(np.mean([np.mean([delta("periscope", rt, rc, s, kind) for rt in routes], 0) for s in seeds], 0) ** 2))
        print(f"  {rc:16s} {kind:16s} seed-half agreement {rel:+.3f} (full-ensemble reliability {sb:+.3f}) | ensemble RMS {amp:.4f} vs observed RMS {np.sqrt((obs[rc]**2).mean()):.4f}")

print("\n=== (c) same aggregation rule as the paper's other aggregate readouts ===")
rows = []
for kind in ["fixed_control", "fixed_PD"]:
    cat_p, cat_o = [], []
    for rc in receivers:
        ens = np.mean([np.mean([delta("periscope", rt, rc, s, kind) for rt in routes], 0) for s in seeds], 0)
        r = stats.spearmanr(ens, obs[rc]).correlation
        per_seed = [stats.spearmanr(np.mean([delta("periscope", rt, rc, s, kind) for rt in routes], 0), obs[rc]).correlation for s in seeds]
        print(f"  {kind:14s} {rc:16s} ten-seed equal-route ensemble rho {r:+.3f} | per-seed equal-route mean {np.mean(per_seed):+.3f}, positive in {int(np.sum(np.array(per_seed)>0))}/10")
        cat_p.append(ens); cat_o.append(obs[rc]); rows.append(dict(condition=kind, receiver=rc, rho_ensemble=r, rho_per_seed_mean=np.mean(per_seed), n_seed_positive=int(np.sum(np.array(per_seed) > 0))))
    pr = stats.spearmanr(np.concatenate(cat_p), np.concatenate(cat_o)).correlation
    print(f"  {kind:14s} {'BOTH RECEIVERS':16s} pooled rho {pr:+.3f}")
    rows.append(dict(condition=kind, receiver="pooled", rho_ensemble=pr, rho_per_seed_mean=np.nan, n_seed_positive=-1))

print("\n=== (d) the same aggregate for every method that has a blood input ===")
for mth in ["periscope", "conditional_mlp", "conditional_icnn", "ot_ridge", "cellot"]:
    out = []
    for rc in receivers:
        sd = [s for s in seeds if os.path.exists(A + f"/results/seed_{s}/{mth}__cDC__{rc}.npz")]
        if not sd: out.append(np.nan); continue
        ens = np.mean([np.mean([delta(mth, rt, rc, s, "fixed_control") for rt in routes], 0) for s in sd], 0)
        out.append(stats.spearmanr(ens, obs[rc]).correlation if np.ptp(ens) > 1e-12 else np.nan)
    pooled = np.nan
    sd = [s for s in seeds if os.path.exists(A + f"/results/seed_{s}/{mth}__cDC__astro.npz")]
    if sd:
        cp = np.concatenate([np.mean([np.mean([delta(mth, rt, rc, s, "fixed_control") for rt in routes], 0) for s in sd], 0) for rc in receivers])
        if np.ptp(cp) > 1e-12: pooled = stats.spearmanr(cp, np.concatenate([obs[rc] for rc in receivers])).correlation
    print(f"  {mth:20s} seeds={len(sd):2d} | MG-AP {out[0]:+.3f} | Astro {out[1]:+.3f} | pooled {pooled:+.3f}")
    rows.append(dict(condition="fixed_control_by_method", receiver=mth, rho_ensemble=pooled, rho_per_seed_mean=out[1], n_seed_positive=len(sd)))
pd.DataFrame(rows).to_csv(OUT + "/bloodonly_aggregate.csv", index=False)
print("\nwrote", OUT + "/bloodonly_aggregate.csv")
