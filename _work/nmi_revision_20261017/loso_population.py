"""Leave-one-seed-out retention of the 16 priority programs, placed in the population of all programs that
passed the full five-seed confirmation test. Hypothesis: retention is governed by how far a program's full-data
q lies below 0.05 (removing one of five confirmation seeds removes 20% of the data), and by recipient state.
usage: python loso_population.py OUTDIR
"""
import sys, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
R = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results/optimization_20261003/seed_robustness"
OUT = sys.argv[1]; pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
s = pd.read_csv(R + "/leave_one_seed_summary.csv")
print("leave_one_seed_summary.csv", s.shape); print(list(s.columns)); print(s.head(3).to_string(index=False))
pr = pd.read_csv(R + "/priority_program_retention.csv")
key = "combo" if "combo" in s.columns else s.columns[0]
qcol = [c for c in s.columns if c in ("full_q", "q_full", "full_FDR_q")]
rcol = [c for c in s.columns if c in ("retained_runs", "n_retained")]
print("\nkey", key, "| q col", qcol, "| retained col", rcol)
if not qcol or not rcol: sys.exit("columns not found")
qcol, rcol = qcol[0], rcol[0]
passed = s[s[qcol] <= 0.05].copy()
print(f"\nprograms passing the full five-seed confirmation test: {len(passed)} (of {len(s)})")
passed["state"] = np.where(passed[key].str.contains("microglia_mhc2"), "MG-AP", np.where(passed[key].str.contains("astro"), "Astro", np.where(passed[key].str.contains("homeostatic"), "MG-H", "other")))
print("\n=== retention distribution among all confirmation-passing programs ===")
print(passed.groupby("state")[rcol].describe().round(2).to_string())
print(f"fraction retained in all five omissions: overall {np.mean(passed[rcol]==5):.3f}")
for st, g in passed.groupby("state"): print(f"  {st:6s} n={len(g):5d} all-five {np.mean(g[rcol]==5):.3f} | >=3 {np.mean(g[rcol]>=3):.3f}")

print("\n=== retention as a function of full-data q (all confirmation-passing programs) ===")
bins = [0, 1e-3, 5e-3, 0.01, 0.02, 0.03, 0.04, 0.05]
passed["qbin"] = pd.cut(passed[qcol], bins, include_lowest=True)
tab = passed.groupby("qbin")[rcol].agg(["count", "mean", lambda x: np.mean(x == 5), lambda x: np.mean(x >= 3)]); tab.columns = ["n", "mean_retained", "frac_all5", "frac_ge3"]
print(tab.round(3).to_string())
rho = stats.spearmanr(passed[qcol], passed[rcol]).correlation
print(f"Spearman(full q, retained runs) = {rho:+.3f}")

print("\n=== the 16 priority programs against programs with similar full q and the same recipient state ===")
rows = []
for _, r in pr.iterrows():
    st = "MG-AP" if "microglia_mhc2" in r.combo else "Astro"
    lo, hi = r.full_q / 2, min(r.full_q * 2, 0.05)
    ref = passed[(passed.state == st) & (passed[qcol] >= lo) & (passed[qcol] <= hi)]
    exp = ref[rcol].mean() if len(ref) else np.nan
    pct = np.mean(ref[rcol] <= r.retained_runs) if len(ref) else np.nan
    rows.append(dict(program=r.combo.replace("relative::", "rel::").replace("raw::", ""), state=st, full_q=r.full_q, retained=int(r.retained_runs), expected_given_q_state=exp, n_ref=len(ref), percentile_in_ref=pct))
t = pd.DataFrame(rows).sort_values(["retained", "full_q"]); print(t.round(4).to_string(index=False))
print(f"\npriority programs: observed mean retention {t.retained.mean():.2f} vs expected {t.expected_given_q_state.mean():.2f} given full q and state")
print(f"programs at or above the expectation for their q and state: {int((t.retained >= t.expected_given_q_state - 1e-9).sum())}/16")
t.to_csv(OUT + "/loso_population.csv", index=False); tab.to_csv(OUT + "/loso_by_q.csv")
print("\nwrote", OUT + "/loso_population.csv")
