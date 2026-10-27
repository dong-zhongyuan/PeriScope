"""Pig, new angle: how much of the brain response is predictable from blood in this cohort at all?
Compare the model against the direct blood-response reference (the animals' own measured blood response used as
if it were the brain response), and against an abundance-stratified gene permutation that holds the labels fixed
and asks whether gene identity adds anything beyond the abundance-driven component.
usage: python pig_upper_bound.py OUTDIR
"""
import sys, itertools, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
P = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results/pig_external_validation_20261003"
OUT = sys.argv[1]; pd.set_option("display.width", 250)
obs = pd.read_csv(P + "/observed_brain_log2_response.csv", index_col=0)
animals = list(obs.columns); lps = np.array([c.endswith("LPS") for c in animals]); O = obs.values
pred = pd.read_csv(P + "/predicted_response_periscope.csv", index_col=0).reindex(obs.index)[animals].values
bt = pd.read_csv(P + "/inputs/blood_one2one_sum_provided_transcript_TPM.csv", index_col=0)[animals].reindex(obs.index)
brt = pd.read_csv(P + "/inputs/brain_one2one_sum_provided_transcript_TPM.csv", index_col=0)[animals].reindex(obs.index)

def response(T):
    v = T.values.astype(float)
    v = v / v.sum(0, keepdims=True) * 1e4
    ctrl = v[:, ~lps].mean(1, keepdims=True)
    return np.log2((v + .1) / (ctrl + .1))

B = response(bt)
labelings = [np.isin(np.arange(10), list(s)) for s in itertools.combinations(range(10), 4)]
oe = O[:, lps].mean(1) - O[:, ~lps].mean(1)

def exact(p):
    r = stats.spearmanr(p, oe).correlation
    null = np.array([stats.spearmanr(p, O[:, m].mean(1) - O[:, ~m].mean(1)).correlation for m in labelings])
    return r, float(np.mean(null >= r - 1e-12))

print("=== how predictable is the pig brain response from blood at all? ===")
be = B[:, lps].mean(1) - B[:, ~lps].mean(1)
pe = pred[:, lps].mean(1) - pred[:, ~lps].mean(1)
for name, p in [("measured blood response used directly", be), ("PeriScope prediction", pe)]:
    r, pv = exact(p); print(f"  {name:40s} rho {r:+.3f} | exact label P {pv:.4f}")
print(f"  Spearman(PeriScope prediction, measured blood response) = {stats.spearmanr(pe, be).correlation:+.3f}")

print("\n=== abundance-stratified gene permutation: does gene identity add anything beyond abundance? ===")
lab = np.log10(brt.values.astype(float).mean(1) + 0.1)
q = pd.qcut(lab, 20, labels=False, duplicates="drop")
rng = np.random.default_rng(0)
def strat_null(p, n=5000):
    out = np.empty(n)
    idx = [np.where(q == s)[0] for s in range(q.max() + 1)]
    for k in range(n):
        perm = np.arange(len(p))
        for ii in idx: perm[ii] = rng.permutation(ii)
        out[k] = stats.spearmanr(p[perm], oe).correlation
    return out
for name, p in [("measured blood response", be), ("PeriScope prediction", pe)]:
    r = stats.spearmanr(p, oe).correlation
    nl = strat_null(p)
    print(f"  {name:28s} rho {r:+.3f} | abundance-stratified null mean {nl.mean():+.3f} sd {nl.std():.3f} | P {np.mean(nl>=r):.4f}")

print("\n=== does the model add anything over the measured blood response? ===")
rk = lambda x: stats.rankdata(x)
Z = np.c_[np.ones(len(be)), rk(be)]
resid_pred = rk(pe) - Z @ np.linalg.lstsq(Z, rk(pe), rcond=None)[0]
resid_obs = rk(oe) - Z @ np.linalg.lstsq(Z, rk(oe), rcond=None)[0]
print(f"  partial Spearman(prediction, observed | measured blood response) = {np.corrcoef(resid_pred, resid_obs)[0,1]:+.3f}")
pd.DataFrame(dict(gene=obs.index, predicted=pe, observed=oe, blood=be, abundance_stratum=q)).to_csv(OUT + "/pig_upper_bound.csv", index=False)
print("\nwrote", OUT + "/pig_upper_bound.csv")
