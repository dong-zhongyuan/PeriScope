"""Pig transfer, fifth pass: is the all-gene agreement carried by a shared dependence on expression level
(a compression artefact of the log-ratio readout) or by gene-specific response prediction?
Partial Spearman controlling mean log expression, with the full exact 210-labeling null recomputed the same way.
usage: python pig_confound.py OUTDIR
"""
import sys, itertools, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
P = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results/pig_external_validation_20261003"
OUT = sys.argv[1]
obs = pd.read_csv(P + "/observed_brain_log2_response.csv", index_col=0)
animals = list(obs.columns); lps = np.array([c.endswith("LPS") for c in animals]); O = obs.values
pred = pd.read_csv(P + "/predicted_response_periscope.csv", index_col=0).reindex(obs.index)[animals].values
T = pd.read_csv(P + "/inputs/brain_one2one_sum_provided_transcript_TPM.csv", index_col=0)[animals].reindex(obs.index).values
lab = np.log10(T.mean(1) + 0.1)
labelings = [np.isin(np.arange(10), list(s)) for s in itertools.combinations(range(10), 4)]
pe = pred[:, lps].mean(1) - pred[:, ~lps].mean(1); oe = O[:, lps].mean(1) - O[:, ~lps].mean(1)

print("=== dependence of each side on expression level ===")
print(f"Spearman(predicted contrast, mean log TPM)      = {stats.spearmanr(pe, lab).correlation:+.3f}")
print(f"Spearman(observed  contrast, mean log TPM)      = {stats.spearmanr(oe, lab).correlation:+.3f}")
print(f"Spearman(|observed contrast|, mean log TPM)     = {stats.spearmanr(np.abs(oe), lab).correlation:+.3f}")
print(f"Spearman(predicted, observed) unadjusted        = {stats.spearmanr(pe, oe).correlation:+.3f}")

def partial(x, y, z):
    rx = stats.rankdata(x); ry = stats.rankdata(y); rz = stats.rankdata(z)
    Z = np.c_[np.ones_like(rz), rz, rz ** 2]
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]; ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])

r_par = partial(pe, oe, lab)
null_par = np.array([partial(pe, O[:, m].mean(1) - O[:, ~m].mean(1), lab) for m in labelings])
null_raw = np.array([stats.spearmanr(pe, O[:, m].mean(1) - O[:, ~m].mean(1)).correlation for m in labelings])
print(f"\npartial Spearman controlling mean log TPM       = {r_par:+.3f} | exact P = {np.mean(null_par >= r_par - 1e-12):.4f}")
print(f"  (unadjusted exact P for comparison            = {np.mean(null_raw >= stats.spearmanr(pe, oe).correlation - 1e-12):.4f})")

print("\n=== agreement within expression strata (quintiles of mean TPM), exact test inside each stratum ===")
q = pd.qcut(lab, 5, labels=False, duplicates="drop")
for s in range(q.max() + 1):
    m = q == s
    r = stats.spearmanr(pe[m], oe[m]).correlation
    nl = np.array([stats.spearmanr(pe[m], O[m][:, b].mean(1) - O[m][:, ~b].mean(1)).correlation for b in labelings])
    print(f"  quintile {s+1} (n={m.sum():4d}, TPM {10**lab[m].min():.2f}-{10**lab[m].max():.1f}): rho={r:+.3f} exact P={np.mean(nl >= r - 1e-12):.4f}")

print("\n=== pooled within-stratum (Fisher-combined over the five quintiles, exact null recomputed identically) ===")
def pooled(obsm):
    zs = []
    for s in range(q.max() + 1):
        m = q == s; r = stats.spearmanr(pe[m], obsm[m]).correlation; zs.append(np.arctanh(np.clip(r, -.999, .999)) * np.sqrt(m.sum() - 3))
    return float(np.sum(zs) / np.sqrt(len(zs)))
z_obs = pooled(oe); z_null = np.array([pooled(O[:, b].mean(1) - O[:, ~b].mean(1)) for b in labelings])
print(f"  pooled within-stratum Stouffer z = {z_obs:+.2f} | exact P = {np.mean(z_null >= z_obs - 1e-12):.4f}")
pd.DataFrame(dict(gene=obs.index, predicted=pe, observed=oe, mean_log_tpm=lab, stratum=q)).to_csv(OUT + "/pig_confound.csv", index=False)
print("\nwrote", OUT + "/pig_confound.csv")
