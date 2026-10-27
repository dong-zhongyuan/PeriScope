"""Cell-type-matched evaluation of the pig cross-species transfer: the model predicts glial (astro / MG-AP)
responses, the observed readout is bulk prefrontal cortex. Partition the 3,730 orthologs by dominant lineage in
the human Ma reference (no pig data used for the partition) and evaluate rank agreement and the exact brain-label
permutation test within each partition, for PeriScope and all comparators. usage: python pig_lineage_eval.py OUTDIR
"""
import sys, itertools, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets"; R = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results"
P = R + "/pig_external_validation_20261003"; OUT = sys.argv[1]
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
obs = pd.read_csv(P + "/observed_brain_log2_response.csv", index_col=0)
animals = list(obs.columns); lps = np.array([c.endswith("LPS") for c in animals]); print("animals:", animals)
pred = {m: pd.read_csv(P + f"/predicted_response_{m}.csv", index_col=0).reindex(obs.index) for m in ["periscope", "conditional_mlp", "conditional_icnn", "ot_ridge", "cellot"]}
eff = pd.read_csv(P + "/gene_level_group_effects.csv", index_col=0).reindex(obs.index)
pred_eff = {m: pred[m][animals].values[:, lps].mean(1) - pred[m][animals].values[:, ~lps].mean(1) for m in pred}
assert np.allclose(pred_eff["periscope"], eff.periscope.values, atol=1e-9), "effect reconstruction mismatch"
rp = pd.read_csv(A + "/results/optimization_20261003/spatial_subtype_author_20261003/inputs/reference_profiles.csv.gz", index_col=0)
lin = {"astro_fibrous": "astro", "astro_protoplasmic": "astro", "microglia_homeostatic": "microglia", "microglia_mhc2": "microglia", "neuron_GABAergic": "neuron", "neuron_MEG3": "neuron", "neuron_dopaminergic": "neuron", "Oligo": "oligo", "OPC": "opc", "Endothelial": "endo", "BAM": "other", "Tcells": "other"}
share = rp.T.groupby(lambda c: lin[c]).sum().T; share = share.div(share.sum(1) + 1e-12, axis=0)
dom = share.idxmax(1).reindex(obs.index).fillna("unannotated")
glia_share = (share["astro"] + share["microglia"]).reindex(obs.index)
parts = {"all": np.ones(len(obs), bool), "glia_dominant": dom.isin(["astro", "microglia"]).values, "astro_dominant": (dom == "astro").values, "microglia_dominant": (dom == "microglia").values,
         "neuron_dominant": (dom == "neuron").values, "oligo_dominant": (dom == "oligo").values, "other/endo/opc": dom.isin(["other", "endo", "opc"]).values, "glia_share>=0.5": (glia_share >= 0.5).values, "glia_share<0.2": (glia_share < 0.2).values}
# receiver marker sets defined from human training data only (their protocol)
mk = pd.read_csv(P + "/human_training_marker_sets.csv"); parts["receiver_markers(human)"] = np.asarray(obs.index.isin(mk.gene))
labelings = list(itertools.combinations(range(len(animals)), int(lps.sum()))); O = obs.values
def exact_p(pe, mask):
    rho = stats.spearmanr(pe[mask], O[mask][:, lps].mean(1) - O[mask][:, ~lps].mean(1)).correlation; null = []
    for sel in labelings:
        m = np.zeros(len(animals), bool); m[list(sel)] = True
        null.append(stats.spearmanr(pe[mask], O[mask][:, m].mean(1) - O[mask][:, ~m].mean(1)).correlation)
    null = np.array(null); return rho, float(np.mean(null >= rho - 1e-12))
rows = []
for name, mask in parts.items():
    r = dict(partition=name, n_genes=int(mask.sum()))
    for m in pred_eff:
        rho, p = exact_p(pred_eff[m], mask); r[m + "_rho"] = round(rho, 3); r[m + "_exactP"] = round(p, 4)
    rows.append(r)
df = pd.DataFrame(rows); df.to_csv(OUT + "/pig_lineage_eval.csv", index=False)
print("=== Spearman agreement and exact brain-label permutation P (210 labelings), by human lineage partition ===")
print(df.to_string())
# animal-level leave-one-out for the glia partition (PeriScope)
mask = parts["glia_dominant"]; loo = []
for i, a in enumerate(animals):
    keep = np.arange(len(animals)) != i; l2 = lps[keep]
    pe = pred["periscope"][animals].values[:, keep]; pe = pe[:, l2].mean(1) - pe[:, ~l2].mean(1); oe = O[:, keep][:, l2].mean(1) - O[:, keep][:, ~l2].mean(1)
    loo.append((a, round(stats.spearmanr(pe[mask], oe[mask]).correlation, 3)))
print("\nglia-dominant partition, PeriScope leave-one-animal-out:", loo)
# sign agreement and magnitude within glia partition
pe = pred_eff["periscope"][mask]; oe = (O[:, lps].mean(1) - O[:, ~lps].mean(1))[mask]
print("glia partition: sign agreement =", round(float(np.mean(np.sign(pe) == np.sign(oe))), 3), "| magnitude ratio =", round(float(np.sqrt(np.mean(pe ** 2)) / np.sqrt(np.mean(oe ** 2))), 4))
# top-responding observed genes: agreement restricted to the top 20% |observed| within glia partition (reported, not selected)
thr = np.quantile(np.abs(oe), 0.8); sel = np.abs(oe) >= thr
print(f"glia partition, top-20% abs(observed) genes (n={int(sel.sum())}): rho = {stats.spearmanr(pe[sel], oe[sel]).correlation:.3f}")
