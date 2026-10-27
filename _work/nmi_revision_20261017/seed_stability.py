"""(1) Do the two route-specific positive blood-only correlations quoted in the manuscript hold across the ten
released models, or only across the three benchmark runs?
(2) Distribution fit: is PeriScope's gap to the best method larger than the run-to-run variation?
usage: python seed_stability.py OUTDIR
"""
import sys, os, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets/results/benchmark_20261002"
F = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/figure_source_data"
OUT = sys.argv[1]; pd.set_option("display.width", 250)
routes = ["cDC", "classical_mono", "nonclassical_mono", "pDC"]; receivers = ["astro", "microglia_mhc2"]
cache = {rc: np.load(A + f"/cache/cDC__{rc}.npz", allow_pickle=True) for rc in receivers}
obs = {rc: cache[rc]["profiles1"].mean(0) - cache[rc]["profiles0"].mean(0) for rc in receivers}

print("=== (1) blood-only route correlations: three benchmark runs vs all ten released models ===")
rows = []
for rc in receivers:
    for rt in routes:
        for kind, key in [("fixed control", "map0_blood"), ("fixed PD", "map1_blood")]:
            per = {}
            for s in range(42, 52):
                f = A + f"/results/seed_{s}/periscope__{rt}__{rc}.npz"
                if not os.path.exists(f): continue
                z = np.load(f); per[s] = stats.spearmanr(z[key + "1"].mean(0) - z[key + "0"].mean(0), obs[rc]).correlation
            three = np.mean([per[s] for s in (42, 43, 44)]); ten = np.mean(list(per.values()))
            npos = int(np.sum(np.array(list(per.values())) > 0))
            rows.append(dict(receiver=rc, route=rt, readout=kind, mean_3_runs=three, mean_10_models=ten, n_positive_of_10=npos))
t = pd.DataFrame(rows)
print(t.round(3).to_string(index=False))
t.to_csv(OUT + "/bloodonly_seed_stability.csv", index=False)
print("\n  quoted in the manuscript: cMono-to-Astro 0.129 / nMono-to-Astro 0.078 under control, 0.103 / 0.070 under PD")
q = t[(t.receiver == "astro") & (t.route.isin(["classical_mono", "nonclassical_mono"]))]
print(q.round(3).to_string(index=False))

print("\n\n=== (2) distribution fit: PeriScope's gap to the best method vs run-to-run variation ===")
b = pd.read_csv(F + "/benchmark_raw.csv"); b = b[b.readout == "raw_shared"]
for col in ["mmd", "mean_mse", "energy_per_sqrt_gene"]:
    print(f"\n  -- {col}")
    g = b.groupby(["method", "seed"])[col].mean().unstack()
    print(g.round(4).to_string())
    per_route = b.pivot_table(index=["blood_state", "brain_state", "seed"], columns="method", values=col)
    best = [m for m in ["conditional_mlp", "conditional_icnn"] if m in per_route.columns]
    d = per_route["periscope"] - per_route[best].min(1)
    print(f"    PeriScope minus best conditional model, over {len(d)} route x run cells: mean {d.mean():+.4f}, sd {d.std():.4f}")
    print(f"    paired Wilcoxon: P = {stats.wilcoxon(per_route['periscope'], per_route[best].min(1)).pvalue:.4f}")
    spread = b[b.method == 'periscope'].groupby(['blood_state', 'brain_state'])[col].std().mean()
    print(f"    PeriScope run-to-run sd within a route: {spread:.4f} | gap/run-sd = {abs(d.mean())/spread:.2f}")
    lab = per_route[[c for c in ['shared_label_mean', 'shared_label_marginal', 'ot_ridge'] if c in per_route.columns]].min(1)
    print(f"    PeriScope vs best non-conditional baseline: mean difference {(per_route['periscope']-lab).mean():+.4f} (PeriScope better by {100*(1-per_route['periscope'].mean()/lab.mean()):.0f}%)")
