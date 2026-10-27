"""Robustness of the Ma 2025 third-cohort replication: donor-score exact permutation test, leave-one-donor-out,
minimum-cell sensitivity, and the Figure 6B donor-score effects (PD minus control) in the assigned Ma state.
usage: python ma_robustness.py DIAGDIR
"""
import json, sys, itertools, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
R = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results/optimization_20261003"
D = sys.argv[1]; rng = np.random.RandomState(42)
pd.set_option("display.width", 260); pd.set_option("display.max_columns", 40)
programs = json.load(open(R + "/joint_analysis/programs.json")); pr = pd.read_csv(R + "/mr_parallel_gate_20261003/priority_all_programs.csv")
pb = pd.read_csv(D + "/ma/ma_pseudobulk_counts.csv", index_col=0); cells = pd.read_csv(D + "/ma/ma_pseudobulk_cells.csv").set_index("group").n_cells
cond = pd.read_csv(D + "/ma/ma_donors.csv").set_index("Donor").Condition.astype(str).to_dict()
cols = pd.DataFrame([c.split("|") for c in pb.columns], columns=["donor", "lineage", "sub"], index=pb.columns)
def donor_logcp(lineage, subs, min_cells):
    sel = cols[(cols.lineage == lineage) & (cols["sub"].isin(subs))]; M = {}; ncell = {}
    for d, g in sel.groupby("donor"):
        n = int(cells[g.index].sum())
        if n < min_cells or d not in cond: continue
        s = pb[g.index].sum(1).values; M[d] = np.log1p(s / max(s.sum(), 1) * 1e4); ncell[d] = n
    M = pd.DataFrame(M, index=pb.index).T; is_pd = np.array([cond[d] == "PD" for d in M.index]); return M, is_pd, ncell
def welch(M, is_pd):
    a_, b_ = M.values[is_pd], M.values[~is_pd]
    t = (a_.mean(0) - b_.mean(0)) / np.sqrt(a_.var(0, ddof=1) / len(a_) + b_.var(0, ddof=1) / len(b_) + 1e-12); t[(M.values > 0).sum(0) < 3] = np.nan
    return pd.Series(t, index=M.columns), pd.Series(M.values.mean(0), index=M.columns)
def matched_p(t, e, gset, B=3000):
    t = t.dropna(); e = e.reindex(t.index); bins = pd.qcut(e.rank(method="first"), 20, labels=False).values; gi = [g for g in gset if g in t.index]
    if len(gi) < 4 or len(gi) / len(gset) < 0.5: return np.nan, np.nan, len(gi)
    tv = t.values; pos = {g: i for i, g in enumerate(t.index)}; gix = np.array([pos[g] for g in gi]); stat = tv[gix].mean()
    pools = {b: np.where(bins == b)[0] for b in np.unique(bins)}; gb = bins[gix]
    null = np.array([np.mean([tv[rng.choice(pools[b])] for b in gb]) for _ in range(B)])
    return round(float(stat), 3), round((1 + np.sum(null <= stat)) / (B + 1), 4), len(gi)
def donor_score(M, gset):
    gi = [g for g in gset if g in M.columns]; Z = M[gi]; Z = (Z - Z.mean(0)) / (Z.std(0, ddof=0) + 1e-9); Z = Z.loc[:, Z.std(0) > 0]
    return Z.mean(1)
def exact_p(score, is_pd):
    n = len(score); k = int(is_pd.sum()); obs = score[is_pd].mean() - score[~is_pd].mean(); diffs = []
    for sel in itertools.combinations(range(n), k):
        m = np.zeros(n, bool); m[list(sel)] = True; diffs.append(score[m].mean() - score[~m].mean())
    diffs = np.array(diffs); return float(obs), float(np.mean(diffs <= obs))  # one-sided (decrease)
states = {"mg_activated": ("Myeloid", ["Activated"]), "astro_fibrous": ("Astrocyte", ["Fibrous"])}
rows = []; effects = {}
for combo in pr.combo:
    gset = programs[combo]; bs = combo.split("__x__")[1].split("__")[0]; st = "astro_fibrous" if bs == "astro" else "mg_activated"; lineage, subs = states[st]
    r = dict(program=combo.replace("__x__", " x ").replace("raw::", "").replace("relative::", "rel::"), state=st)
    for mc in [25, 50]:
        M, is_pd, ncell = donor_logcp(lineage, subs, mc); t, e = welch(M, is_pd); stat, p, n = matched_p(t, e, gset)
        r[f"min{mc}_nPD"] = int(is_pd.sum()); r[f"min{mc}_nC"] = int((~is_pd).sum()); r[f"min{mc}_stat"] = stat; r[f"min{mc}_p_down"] = p
        if mc == 25:
            sc = donor_score(M, gset); eff, pex = exact_p(sc.values, is_pd); r["score_effect_PD_minus_C"] = round(eff, 3); r["exact_p_one_sided"] = round(pex, 4)
            effects[combo] = round(eff, 4)
            # leave-one-donor-out on the competitive statistic
            loo = []
            for d in M.index:
                keep = M.index != d
                if is_pd[keep].sum() < 2 or (~is_pd[keep]).sum() < 2: continue
                t2, e2 = welch(M[keep], is_pd[keep]); loo.append(matched_p(t2, e2, gset, B=600)[1])
            r["loo_max_p_down"] = round(float(np.nanmax(loo)), 4) if loo else np.nan; r["loo_n"] = len(loo)
    rows.append(r)
df = pd.DataFrame(rows); df.to_csv(D + "/ma_robustness.csv", index=False)
print(df.to_string())
pd.Series(effects, name="cell_state_sensitivity_ma_effect").to_csv(D + "/ma_fig6b_effects.csv")
print("\nMG-AP programs exact p<=0.05:", int((df[df.state == "mg_activated"].exact_p_one_sided <= 0.05).sum()), "/", int((df.state == "mg_activated").sum()), "| astro programs:", int((df[df.state == "astro_fibrous"].exact_p_one_sided <= 0.05).sum()), "/", int((df.state == "astro_fibrous").sum()))
