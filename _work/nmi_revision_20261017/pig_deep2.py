"""Pig transfer, second pass: what is the dominant label-free brain axis (PC1), does the prediction capture the
LPS-specific axis beyond it, is the per-animal agreement valid under exact relabeling, and does sampling time explain
the animal-level pattern. usage: python pig_deep2.py OUTDIR
"""
import sys, itertools, warnings, json
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets"; R = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results"; P = R + "/pig_external_validation_20261003"; OUT = sys.argv[1]
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
obs = pd.read_csv(P + "/observed_brain_log2_response.csv", index_col=0); animals = list(obs.columns); lps = np.array([c.endswith("LPS") for c in animals])
pred = pd.read_csv(P + "/predicted_response_periscope.csv", index_col=0).reindex(obs.index); O = obs.values; Pm = pred[animals].values
pe = Pm[:, lps].mean(1) - Pm[:, ~lps].mean(1); oe = O[:, lps].mean(1) - O[:, ~lps].mean(1)
labelings = [np.isin(np.arange(len(animals)), list(sel)) for sel in itertools.combinations(range(len(animals)), int(lps.sum()))]
# ---- PC1 characterization
Oc = O - O.mean(1, keepdims=True); U, S, Vt = np.linalg.svd(Oc, full_matrices=False); load1 = U[:, 0] * S[0]; load2 = U[:, 1] * S[1]
rp = pd.read_csv(A + "/results/optimization_20261003/spatial_subtype_author_20261003/inputs/reference_profiles.csv.gz", index_col=0)
lin = {"astro_fibrous": "astro", "astro_protoplasmic": "astro", "microglia_homeostatic": "microglia", "microglia_mhc2": "microglia", "neuron_GABAergic": "neuron", "neuron_MEG3": "neuron", "neuron_dopaminergic": "neuron", "Oligo": "oligo", "OPC": "opc", "Endothelial": "endo", "BAM": "other", "Tcells": "other"}
share = rp.T.groupby(lambda c: lin[c]).sum().T; share = share.div(share.sum(1) + 1e-12, axis=0).reindex(obs.index)
print("=== PC1 loading vs human lineage share (Spearman) ===")
for c in ["neuron", "astro", "microglia", "oligo", "endo"]:
    ok = share[c].notna(); print(f"  {c:10s} PC1 {stats.spearmanr(load1[ok], share.loc[ok, c]).correlation:+.3f} | PC2 {stats.spearmanr(load2[ok], share.loc[ok, c]).correlation:+.3f} | predicted effect {stats.spearmanr(pe[ok], share.loc[ok, c]).correlation:+.3f} | observed LPS effect {stats.spearmanr(oe[ok], share.loc[ok, c]).correlation:+.3f}")
top = pd.Series(load1, index=obs.index).sort_values(); print("PC1 most negative genes:", list(top.index[:15])); print("PC1 most positive genes:", list(top.index[-15:]))
# known marker checks
for name, gs in {"immediate-early/inflammation": ["FOS", "JUNB", "EGR1", "IL1B", "TNF", "CCL2", "NFKBIA", "SOCS3", "CXCL10", "IRF1", "STAT1", "ISG15", "MX1", "OAS1"], "microglia": ["AIF1", "CX3CR1", "P2RY12", "TMEM119", "CSF1R", "C1QA", "C1QB"], "astro": ["GFAP", "AQP4", "SLC1A3", "ALDH1L1"], "neuron": ["SNAP25", "SYT1", "RBFOX3", "GAD1", "SLC17A7"], "oligo": ["MBP", "PLP1", "MOG"], "mito/ribo-like": [g for g in obs.index if g.startswith("MT-") or g.startswith("RPL") or g.startswith("RPS")][:40]}.items():
    gs = [g for g in gs if g in obs.index]
    if gs: print(f"  {name:28s} n={len(gs):2d} | mean PC1 loading {np.mean(pd.Series(load1, index=obs.index)[gs]):+.3f} | mean PC2 loading {np.mean(pd.Series(load2, index=obs.index)[gs]):+.3f} | mean observed LPS effect {np.mean(pd.Series(oe, index=obs.index)[gs]):+.3f} | mean predicted {np.mean(pd.Series(pe, index=obs.index)[gs]):+.4f}")
# ---- does the prediction capture LPS signal beyond PC1? residualize observed responses on PC1 (label-free), recompute
O_res = Oc - np.outer(load1, Vt[0]); oe_res = O_res[:, lps].mean(1) - O_res[:, ~lps].mean(1)
r_res = stats.spearmanr(pe, oe_res).correlation; null = [stats.spearmanr(pe, O_res[:, m].mean(1) - O_res[:, ~m].mean(1)).correlation for m in labelings]
print(f"\n=== after removing PC1 from observed data: rho = {r_res:.3f} | exact P = {np.mean(np.array(null) >= r_res - 1e-12):.4f}")
O_res2 = O_res - np.outer(load2, Vt[1]); oe_res2 = O_res2[:, lps].mean(1) - O_res2[:, ~lps].mean(1); print(f"after removing PC1 and PC2: rho = {stats.spearmanr(pe, oe_res2).correlation:.3f}")
print(f"prediction vs PC2 loading (the LPS-separating axis): rho = {stats.spearmanr(pe, load2).correlation:+.3f}")
# ---- per-animal MWU statistic under exact relabeling (statistic recomputed per labeling)
def per_animal_U(m):
    rhos = []
    for i in range(len(animals)):
        if m[i]: od = O[:, i] - O[:, ~m].mean(1); pdv = Pm[:, i] - Pm[:, ~m].mean(1)
        else:
            others = (~m) & (np.arange(len(animals)) != i); od = O[:, i] - O[:, others].mean(1); pdv = Pm[:, i] - Pm[:, others].mean(1)
        rhos.append(stats.spearmanr(pdv, od).correlation)
    rhos = np.array(rhos); return stats.mannwhitneyu(rhos[m], rhos[~m], alternative="greater").statistic, rhos
U_obs, rhos_obs = per_animal_U(lps); nullU = np.array([per_animal_U(m)[0] for m in labelings])
print(f"\n=== per-animal agreement, exact relabeling of the whole statistic: U = {U_obs:.0f} | exact P = {np.mean(nullU >= U_obs - 1e-9):.4f} | 90th pct of null U = {np.quantile(nullU, .9):.1f}")
# ---- sampling time
try:
    ma = pd.read_csv(P + "/inputs/matched_animals.csv"); print("\n=== matched animals ==="); print(ma.to_string(index=False))
except Exception as e: print("matched_animals:", e)
pa = pd.DataFrame(dict(animal=animals, LPS=lps, per_animal_rho=np.round(rhos_obs, 3), PC1=np.round(Vt[0], 3), PC2=np.round(Vt[1], 3))); print(pa.to_string(index=False))
pa.to_csv(OUT + "/pig_per_animal_axes.csv", index=False)
