"""Deeper diagnosis of the pig cross-species transfer.
(1) Why is the exact label-permutation P only 0.10? Test whether the predicted LPS effect aligns with a label-free
    principal axis of between-animal brain variation.
(2) Animal-level agreement: per-animal Spearman between predicted and observed deviation from control means;
    LPS vs control animals (exact Mann-Whitney, 4 vs 6).
(3) Robust-response gene subset with a label-permutation-valid definition (genes whose observed response has a
    consistent sign across all LPS animals under each labeling).
(4) Per-route agreement (all_six_axis_metrics) and the amplitude mechanism relative to the human blood-only amplitude.
usage: python pig_deep_diag.py OUTDIR
"""
import sys, itertools, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
R = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results"; P = R + "/pig_external_validation_20261003"; OUT = sys.argv[1]
F = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/figure_source_data"
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
obs = pd.read_csv(P + "/observed_brain_log2_response.csv", index_col=0); animals = list(obs.columns); lps = np.array([c.endswith("LPS") for c in animals])
pred = pd.read_csv(P + "/predicted_response_periscope.csv", index_col=0).reindex(obs.index)
O = obs.values; Pm = pred[animals].values
pe = Pm[:, lps].mean(1) - Pm[:, ~lps].mean(1); oe = O[:, lps].mean(1) - O[:, ~lps].mean(1)
print("primary rho:", round(stats.spearmanr(pe, oe).correlation, 3))
# ---- (1) label-free principal axis of brain variation
Oc = O - O.mean(1, keepdims=True); U, S, Vt = np.linalg.svd(Oc, full_matrices=False)
var_expl = S ** 2 / (S ** 2).sum()
print("\n=== between-animal brain variation: PC variance explained ===", np.round(var_expl[:5], 3))
for k in range(3):
    load = U[:, k] * S[k]; scores = Vt[k]
    print(f"PC{k+1}: Spearman(predicted effect, PC loading) = {stats.spearmanr(pe, load).correlation:+.3f} | LPS vs control scores: {np.round(scores[lps], 2)} vs {np.round(scores[~lps], 2)} | MWU p = {stats.mannwhitneyu(scores[lps], scores[~lps]).pvalue:.3f}")
# how often does a random labeling produce a contrast aligned with PC1?
labelings = list(itertools.combinations(range(len(animals)), int(lps.sum())))
null = []
for sel in labelings:
    m = np.zeros(len(animals), bool); m[list(sel)] = True; null.append(stats.spearmanr(pe, O[:, m].mean(1) - O[:, ~m].mean(1)).correlation)
null = np.array(null); print("exact null: mean %.3f, sd %.3f, 90th pct %.3f, observed %.3f" % (null.mean(), null.std(), np.quantile(null, .9), stats.spearmanr(pe, oe).correlation))
# ---- (2) animal-level agreement
ctrl_mean_o = O[:, ~lps].mean(1, keepdims=True); ctrl_mean_p = Pm[:, ~lps].mean(1, keepdims=True)
per_animal = []
for i, a in enumerate(animals):
    # deviation from the other controls' mean for control animals (leave-self-out) to avoid self-reference
    if lps[i]: od = O[:, i] - ctrl_mean_o[:, 0]; pdv = Pm[:, i] - ctrl_mean_p[:, 0]
    else:
        others = (~lps) & (np.arange(len(animals)) != i); od = O[:, i] - O[:, others].mean(1); pdv = Pm[:, i] - Pm[:, others].mean(1)
    per_animal.append(dict(animal=a, LPS=bool(lps[i]), rho=round(stats.spearmanr(pdv, od).correlation, 3), pearson=round(np.corrcoef(pdv, od)[0, 1], 3)))
pa = pd.DataFrame(per_animal); print("\n=== per-animal agreement (predicted vs observed deviation from control mean) ===\n", pa.to_string(index=False))
u = stats.mannwhitneyu(pa.rho[pa.LPS], pa.rho[~pa.LPS], alternative="greater"); print("LPS > control (one-sided Mann-Whitney):", u)
# ---- (3) robust-response subset with permutation-valid definition
def robust_rho(m):
    o_lps = O[:, m]; o_ctl = O[:, ~m]; eff = o_lps.mean(1) - o_ctl.mean(1)
    consistent = (np.sign(o_lps) == np.sign(eff)[:, None]).all(1)  # every treated animal moves in the effect direction
    return stats.spearmanr(pe[consistent], eff[consistent]).correlation, int(consistent.sum())
r_obs, n_obs = robust_rho(lps); null_r = []
for sel in labelings:
    m = np.zeros(len(animals), bool); m[list(sel)] = True; null_r.append(robust_rho(m)[0])
null_r = np.array(null_r); print(f"\n=== robust-response subset (sign-consistent across treated animals, recomputed per labeling) ===\nobserved rho = {r_obs:.3f} on {n_obs} genes | exact P = {np.mean(null_r >= r_obs - 1e-12):.4f} | null mean {np.nanmean(null_r):.3f}")
# same with magnitude threshold: top-25% |effect| genes (recomputed per labeling)
def top_rho(m, q=0.75):
    eff = O[:, m].mean(1) - O[:, ~m].mean(1); sel = np.abs(eff) >= np.quantile(np.abs(eff), q); return stats.spearmanr(pe[sel], eff[sel]).correlation
t_obs = top_rho(lps); null_t = np.array([top_rho(np.isin(np.arange(len(animals)), list(sel))) for sel in labelings])
print(f"top-25% |observed effect| genes (recomputed per labeling): rho = {t_obs:.3f} | exact P = {np.mean(null_t >= t_obs - 1e-12):.4f}")
# ---- (4) per-route metrics and amplitude mechanism
ax = pd.read_csv(P + "/all_six_axis_metrics.csv"); print("\n=== per-route metrics (PeriScope, condition 0) ===")
cols = [c for c in ["axis", "method", "condition", "subset", "spearman", "exact_group_p_positive", "magnitude_ratio", "sign_agreement"] if c in ax.columns]
sub = ax[(ax.method == "periscope") & (ax.condition == 0)] if "condition" in ax.columns else ax[ax.method == "periscope"]
print(sub[cols].round(3).to_string(index=False)[:2500])
b = pd.read_csv(F + "/benchmark_raw.csv"); bp = b[b.method == "periscope"]
print("\n=== human benchmark amplitude: predicted/observed RMS of the contrast ===")
print("condition-aware:", round((bp.condition_aware_predicted_rms / bp.condition_aware_observed_rms).mean(), 3), "| fixed control (blood-only):", round((bp.fixed_control_predicted_rms / bp.fixed_control_observed_rms).mean(), 3), "| fixed PD:", round((bp.fixed_PD_predicted_rms / bp.fixed_PD_observed_rms).mean(), 3))
pa.to_csv(OUT + "/pig_per_animal.csv", index=False)
