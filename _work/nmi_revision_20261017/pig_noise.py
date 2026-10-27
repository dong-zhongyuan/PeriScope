"""Pig transfer, fourth pass: locate the noise in the OBSERVED contrast.
(a) how much of the 3,730-gene endpoint is carried by genes that are essentially unmeasured;
(b) label-free expression filters (fixed before any labelling) with the exact 210-labeling test;
(c) a variance-stabilised statistic (moderated t) recomputed in full under every labeling;
(d) the released high-confidence ortholog subset.
usage: python pig_noise.py OUTDIR
"""
import sys, itertools, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
P = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results/pig_external_validation_20261003"
OUT = sys.argv[1]
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)

obs = pd.read_csv(P + "/observed_brain_log2_response.csv", index_col=0)
animals = list(obs.columns); lps = np.array([c.endswith("LPS") for c in animals]); O = obs.values
pred = pd.read_csv(P + "/predicted_response_periscope.csv", index_col=0).reindex(obs.index)[animals].values
tpm = pd.read_csv(P + "/inputs/brain_one2one_sum_provided_transcript_TPM.csv", index_col=0)[animals].reindex(obs.index)
T = tpm.values
labelings = [np.isin(np.arange(10), list(s)) for s in itertools.combinations(range(10), 4)]

print("=== measured depth of the 3,730 endpoint genes (provided pig brain TPM) ===")
print(f"genes with TPM=0 in at least one animal: {(T==0).any(1).sum()} | in at least half: {((T==0).sum(1)>=5).sum()} | in all: {(T==0).all(1).sum()}")
for q in [0, 10, 25, 50, 75, 90]:
    print(f"  {q:2d}th percentile of mean TPM: {np.percentile(T.mean(1), q):.3f}")
lowm = T.mean(1) < 1
print(f"genes with mean TPM < 1: {lowm.sum()} ({lowm.mean():.1%})")
print(f"  their share of total variance of the observed response: {np.var(O[lowm])*lowm.sum()/(np.var(O)*len(O)):.3f}")
print(f"  sd of observed group effect: low-expression genes {(O[lowm][:,lps].mean(1)-O[lowm][:,~lps].mean(1)).std():.3f} vs rest {(O[~lowm][:,lps].mean(1)-O[~lowm][:,~lps].mean(1)).std():.3f}")

def exact(predv, obsm, sub=None, stat="mean"):
    if sub is None: sub = np.ones(len(obsm), bool)
    pe = predv[sub][:, lps].mean(1) - predv[sub][:, ~lps].mean(1)
    def S(m):
        X = obsm[sub]
        if stat == "mean":
            return X[:, m].mean(1) - X[:, ~m].mean(1)
        a, b = X[:, m], X[:, ~m]
        df = X.shape[1] - 2
        s2 = ((a.var(1, ddof=1) * (m.sum() - 1)) + (b.var(1, ddof=1) * ((~m).sum() - 1))) / df
        s2m = (df * s2 + df * np.median(s2)) / (2 * df)
        return (a.mean(1) - b.mean(1)) / np.sqrt(s2m * (1 / m.sum() + 1 / (~m).sum()))
    r = stats.spearmanr(pe, S(lps)).correlation
    null = np.array([stats.spearmanr(pe, S(m)).correlation for m in labelings])
    return r, float(np.mean(null >= r - 1e-12)), int(sub.sum())

print("\n=== label-free expression filters, group-mean statistic (exact 210-labeling test) ===")
rows = []
for name, sub in [("all mapped genes", None),
                  ("mean TPM >= 0.5", T.mean(1) >= 0.5),
                  ("mean TPM >= 1", T.mean(1) >= 1),
                  ("mean TPM >= 5", T.mean(1) >= 5),
                  ("mean TPM >= 10", T.mean(1) >= 10),
                  ("TPM > 0 in all 10", (T > 0).all(1)),
                  ("TPM >= 1 in all 10", (T >= 1).all(1))]:
    r, p, n = exact(pred, O, sub); print(f"  {name:22s} n={n:5d}  rho={r:+.3f}  exact P={p:.4f}"); rows.append(dict(subset=name, n=n, rho=r, p=p, stat="group mean"))

print("\n=== same filters, variance-stabilised (moderated t) statistic recomputed under every labeling ===")
for name, sub in [("all mapped genes", None), ("mean TPM >= 1", T.mean(1) >= 1), ("TPM >= 1 in all 10", (T >= 1).all(1)), ("mean TPM >= 5", T.mean(1) >= 5)]:
    r, p, n = exact(pred, O, sub, stat="modt"); print(f"  {name:22s} n={n:5d}  rho={r:+.3f}  exact P={p:.4f}"); rows.append(dict(subset=name, n=n, rho=r, p=p, stat="moderated t"))

aud = pd.read_csv(P + "/inputs/existing_TPM_feature_audit.csv")
hc = aud[aud.tissue == "Brain"].set_index("model_gene").high_confidence.reindex(obs.index).fillna(False).values.astype(bool)
print("\n=== released high-confidence ortholog subset ===")
for stat in ["mean", "modt"]:
    r, p, n = exact(pred, O, hc, stat=stat); print(f"  high-confidence orthologs ({stat}): n={n} rho={r:+.3f} exact P={p:.4f}"); rows.append(dict(subset="high-confidence orthologs", n=n, rho=r, p=p, stat=stat))
    s2 = hc & (T.mean(1) >= 1); r, p, n = exact(pred, O, s2, stat=stat); print(f"  high-confidence AND mean TPM>=1 ({stat}): n={n} rho={r:+.3f} exact P={p:.4f}"); rows.append(dict(subset="high-confidence and mean TPM>=1", n=n, rho=r, p=p, stat=stat))
pd.DataFrame(rows).to_csv(OUT + "/pig_noise_filters.csv", index=False)
print("\nwrote", OUT + "/pig_noise_filters.csv")
