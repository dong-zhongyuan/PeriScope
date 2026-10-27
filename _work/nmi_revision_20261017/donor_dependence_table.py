"""Table S9 numbers: how much does the predicted brain profile move with the identity of the held-out blood
donor supplied, relative to how much it moves with that donor's disease label? Ten frozen models.
usage: python donor_dependence_table.py OUTDIR
"""
import sys, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets"; OUT = sys.argv[1]
pd.set_option("display.width", 250)
BS = ["cDC", "classical_mono", "nonclassical_mono", "pDC"]; RS = ["astro", "microglia_mhc2"]
LBL = {"astro": "Astro", "microglia_mhc2": "MG-AP"}
RT = {"cDC": "cDC", "classical_mono": "cMono", "nonclassical_mono": "nMono", "pDC": "pDC"}
b = dict(np.load(A + "/interim/rescreen_20261002/blood.npz")); dn = b["donor_names"].astype(str)
held = b["is_val"] | b["is_test"]
z = np.load(OUT + "/bloodonly_profiles.npz", allow_pickle=True)
prof = {}
for k in z.files:
    s, bs, rs, cb, d = k.split("|"); prof[(int(s), bs, rs, int(cb), int(d))] = z[k]
seeds = sorted({k[0] for k in prof}); print("seeds:", seeds)
ctrl = sorted({int(d) for d in np.unique(b["donor"][held & (b["cond2"] == 0)])})
pdd = sorted({int(d) for d in np.unique(b["donor"][held & (b["cond2"] == 1)])})
print("control donors:", [dn[d] for d in ctrl], "| PD donors:", [dn[d] for d in pdd])

rows = []
for rs in RS:
    for bs in BS:
        c0 = [np.mean([prof[(s, bs, rs, 0, d)] for s in seeds if (s, bs, rs, 0, d) in prof], 0) for d in ctrl if any((s, bs, rs, 0, d) in prof for s in seeds)]
        c1 = [np.mean([prof[(s, bs, rs, 0, d)] for s in seeds if (s, bs, rs, 0, d) in prof], 0) for d in pdd if any((s, bs, rs, 0, d) in prof for s in seeds)]
        if len(c0) + len(c1) < 3: 
            rows.append(dict(receiver=LBL[rs], route=RT[bs], n_control=len(c0), n_pd=len(c1), between=np.nan, within=np.nan, ratio=np.nan)); continue
        c0 = np.array(c0); c1 = np.array(c1); allp = np.concatenate([c0, c1])
        between = float(np.sqrt(((c1.mean(0) - c0.mean(0)) ** 2).mean()))
        within = float(np.sqrt(np.mean([((x - allp.mean(0)) ** 2).mean() for x in allp])))
        rows.append(dict(receiver=LBL[rs], route=RT[bs], n_control=len(c0), n_pd=len(c1), between=between, within=within, ratio=within / between))
t = pd.DataFrame(rows)
print("\n=== donor-driven vs label-driven variation of the predicted brain profile (ten models) ===")
print(t.round(4).to_string(index=False))
sub = t[t.route != "pDC"].dropna()
print(f"\nacross the cDC and monocyte routes: ratio {sub.ratio.min():.0f} to {sub.ratio.max():.0f} | all routes {t.ratio.min():.0f} to {t.ratio.max():.0f}")
t.to_csv(OUT + "/donor_dependence.csv", index=False)

print("\n=== LaTeX rows ===")
for _, r in t.iterrows():
    if np.isnan(r.between): print(f"{r.receiver} & {r.route} & {int(r.n_control)} / {int(r.n_pd)} & -- & -- & -- \\\\")
    else: print(f"{r.receiver} & {r.route} & {int(r.n_control)} / {int(r.n_pd)} & {r.between:.4f} & {r.within:.4f} & {r.ratio:.0f} \\\\")
