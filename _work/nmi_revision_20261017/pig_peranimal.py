"""Pig transfer, sixth pass: does the animal-level agreement survive adjustment for expression level?
The whole statistic (per-animal rank agreement, Mann-Whitney over animals) is recomputed under each of the
210 treatment labelings, so the test stays exact.
usage: python pig_peranimal.py OUTDIR
"""
import sys, itertools, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
P = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results/pig_external_validation_20261003"
OUT = sys.argv[1]
obs = pd.read_csv(P + "/observed_brain_log2_response.csv", index_col=0)
animals = list(obs.columns); lps = np.array([c.endswith("LPS") for c in animals]); O = obs.values
T = pd.read_csv(P + "/inputs/brain_one2one_sum_provided_transcript_TPM.csv", index_col=0)[animals].reindex(obs.index).values
lab = np.log10(T.mean(1) + 0.1); rz = stats.rankdata(lab); Z = np.c_[np.ones_like(rz), rz, rz ** 2]
Zp = Z @ np.linalg.pinv(Z)
methods = ["periscope", "conditional_mlp", "cellot", "conditional_icnn", "ot_ridge"]
preds = {}
for m in methods:
    try: preds[m] = pd.read_csv(P + f"/predicted_response_{m}.csv", index_col=0).reindex(obs.index)[animals].values
    except Exception: pass
labelings = [np.isin(np.arange(10), list(s)) for s in itertools.combinations(range(10), 4)]

def rhos_for(Pm, m, adjust):
    out = []
    for i in range(10):
        if m[i]: ref = ~m
        else: ref = (~m) & (np.arange(10) != i)
        od = O[:, i] - O[:, ref].mean(1); pv = Pm[:, i] - Pm[:, ref].mean(1)
        if adjust:
            a = stats.rankdata(od); b = stats.rankdata(pv)
            a = a - Zp @ a; b = b - Zp @ b
            out.append(float(np.corrcoef(a, b)[0, 1]))
        else:
            out.append(stats.spearmanr(pv, od).correlation)
    return np.array(out)

print("=== per-animal agreement, exact relabeling of the whole statistic ===")
rows = []
for adjust in [False, True]:
    tag = "expression-adjusted" if adjust else "unadjusted        "
    for name, Pm in preds.items():
        r = rhos_for(Pm, lps, adjust)
        U = stats.mannwhitneyu(r[lps], r[~lps], alternative="greater").statistic
        null = np.array([stats.mannwhitneyu(*(lambda x, mm: (x[mm], x[~mm]))(rhos_for(Pm, mm, adjust), mm), alternative="greater").statistic for mm in labelings])
        p = float(np.mean(null >= U - 1e-9))
        if name == "periscope":
            print(f"  [{tag}] {name:18s} U={U:.0f}/24  exact P={p:.4f}")
            print(f"      LPS animals  {np.round(r[lps],3)}  (n=4)")
            print(f"      controls     {np.round(r[~lps],3)}  (n=6)")
        else:
            print(f"  [{tag}] {name:18s} U={U:.0f}/24  exact P={p:.4f}")
        rows.append(dict(statistic="per-animal Mann-Whitney", adjusted=adjust, method=name, U=float(U), exact_p=p,
                         lps_rho=";".join(f"{x:.3f}" for x in r[lps]), control_rho=";".join(f"{x:.3f}" for x in r[~lps])))
pd.DataFrame(rows).to_csv(OUT + "/pig_per_animal_exact.csv", index=False)
print("\nwrote", OUT + "/pig_per_animal_exact.csv")
