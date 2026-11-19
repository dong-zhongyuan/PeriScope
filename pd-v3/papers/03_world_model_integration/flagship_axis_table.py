"""Two separate ceiling-normalised tables on the discriminating axis.

The merge check failed on purpose-built evidence: the benchmark harness pins every
method to one frozen brain representation and decoder and swaps only the blood-side
map, while the ablation pipeline measures each model end to end through its own
decoder. On the same axis and seed the two differ by up to 0.0595, so they are not
interchangeable and are reported separately.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

B = Path("/public/home/mengxl/dzy/pd_product_assets/results/benchmark_20261002")
A = Path("/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003")
OUT = Path("/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/flagship_axis_20261118")
OUT.mkdir(parents=True, exist_ok=True)
CEILING = {"astro": 0.5537, "microglia_mhc2": 0.6116}
SEEDS = [42, 43, 44]


def norm(df):
    df = df.copy()
    df["brain_state"] = df.axis.str.split("__").str[1]
    df["ceiling"] = df.brain_state.map(CEILING)
    df["pct_of_ceiling"] = 100 * df.rho / df.ceiling
    return df


def summarise(df, label):
    s = df.groupby("method").agg(rho_mean=("rho", "mean"), rho_sd=("rho", "std"),
                                 pct_of_ceiling=("pct_of_ceiling", "mean"),
                                 cells=("rho", "size")).sort_values("rho_mean", ascending=False)
    print("=== %s ===" % label)
    print(s.round(4).to_string())
    return s


def paired(df, reference):
    base = df[df.method.eq(reference)].set_index(["axis", "seed"]).rho
    rows = []
    for m in df.method.unique():
        if m == reference:
            continue
        other = df[df.method.eq(m)].set_index(["axis", "seed"]).rho
        j = pd.concat([base.rename("ref"), other.rename("oth")], axis=1).dropna()
        d = j.ref - j.oth
        rows.append(dict(vs=m, n_pairs=len(d), mean_diff=d.mean(), sd_diff=d.std(),
                         reference_wins="%d/%d" % (int((d > 0).sum()), len(d))))
    return pd.DataFrame(rows).sort_values("mean_diff", ascending=False)


# ---------- Table A: benchmark harness, shared frozen representation ----------
p = pd.read_csv(B / "paired_comparisons.csv")
p = p[p.readout.eq("raw_shared")]
rows = []
for _, r in p.iterrows():
    rows.append(dict(method="periscope", axis=r["axis"], seed=int(r["seed"]), rho=r["condition_aware_rho_periscope"]))
    rows.append(dict(method=r["comparator"], axis=r["axis"], seed=int(r["seed"]), rho=r["condition_aware_rho_other"]))
tabA = norm(pd.DataFrame(rows).drop_duplicates(subset=["method", "axis", "seed"]))
sA = summarise(tabA, "A: benchmark harness (one frozen brain representation, blood-side swapped), seeds 42-44, 8 axes")
print()
print("paired vs periscope, same axis and seed:")
pA = paired(tabA, "periscope")
print(pA.round(4).to_string())

# cellot lives only in the evaluation summary, not in paired_comparisons
ev = pd.read_csv(B / "benchmark_evaluation_summary.csv")
ce = ev[ev.method.eq("cellot") & ev.metric.eq("condition_aware_rho")]
if len(ce):
    print()
    print("cellot condition_aware_rho: mean %.4f over %d axis rows, n_seeds=%s"
          % (ce.estimate.mean(), len(ce), sorted(ce.n_seeds.unique())))
    print("  cellot has no seed-matched paired test: it is absent from paired_comparisons.csv")


def per_seed_view(df, reference, label):
    """Average over axes within a seed first: axes are not independent."""
    s = df.groupby(["method", "seed"]).rho.mean().unstack("seed")
    s["mean"] = s.mean(axis=1)
    s["sd_across_seeds"] = df.groupby(["method", "seed"]).rho.mean().unstack("seed").std(axis=1)
    s["pct_of_ceiling"] = 100 * s["mean"] / np.mean(list(CEILING.values()))
    s = s.sort_values("mean", ascending=False)
    print()
    print("=== %s: axis-averaged within seed, then across seeds ===" % label)
    print(s.round(4).to_string())
    base = df[df.method.eq(reference)].groupby("seed").rho.mean()
    out = []
    for m in df.method.unique():
        if m == reference:
            continue
        o = df[df.method.eq(m)].groupby("seed").rho.mean()
        d = (base - o).dropna()
        out.append(dict(vs=m, n_seeds=len(d), mean_diff=d.mean(), min_diff=d.min(),
                        max_diff=d.max(), reference_wins="%d/%d" % (int((d > 0).sum()), len(d))))
    o = pd.DataFrame(out).sort_values("mean_diff", ascending=False)
    print()
    print("paired by seed vs %s (every seed is one paired observation):" % reference)
    print(o.round(4).to_string())
    return s, o

# ---------- Table B: end-to-end own decoder, the ablation pipeline ----------
SRC = {"full": A / "evaluation",
       "coupling_off": A / "ablations/coupling_off/evaluation",
       "random_pairing": A / "ablations/random_pairing/evaluation"}
rows = []
for label, where in SRC.items():
    for s in SEEDS:
        f = where / ("brain_direction_seed%d.csv" % s)
        if not f.exists():
            continue
        d = pd.read_csv(f)
        for _, r in d[d["mode"].eq("condition_aware")].iterrows():
            rows.append(dict(method=label, axis="%s__%s" % (r["blood_state"], r["brain_state"]),
                             seed=int(r["seed"]), rho=r["rho"]))
tabB = norm(pd.DataFrame(rows))
print()
sB = summarise(tabB, "B: end-to-end through each model's own decoder, seeds 42-44, 6 axes")
print()
print("paired vs full, same axis and seed:")
pB = paired(tabB, "full")
print(pB.round(4).to_string())

# fixed-mode diagnostic for table B
rows = []
for label, where in SRC.items():
    for s in SEEDS:
        f = where / ("brain_direction_seed%d.csv" % s)
        if f.exists():
            d = pd.read_csv(f)
            for _, r in d[d["mode"].ne("condition_aware")].iterrows():
                rows.append(dict(method=label, mode=r["mode"], rho=r["rho"]))
fx = pd.DataFrame(rows)
print()
print("=== blood-contribution diagnostic (label held fixed), table B models ===")
print(fx.pivot_table(index="mode", columns="method", values="rho", aggfunc="median").round(4).to_string())

sA2, pA2 = per_seed_view(tabA, "periscope", "A: benchmark harness")
sB2, pB2 = per_seed_view(tabB, "full", "B: end-to-end")
sA2.to_csv(OUT / "tableA_per_seed.csv")
pA2.to_csv(OUT / "tableA_per_seed_paired.csv", index=False)
sB2.to_csv(OUT / "tableB_per_seed.csv")
pB2.to_csv(OUT / "tableB_per_seed_paired.csv", index=False)
tabA.to_csv(OUT / "tableA_benchmark_harness_per_cell.csv", index=False)
sA.to_csv(OUT / "tableA_summary.csv")
pA.to_csv(OUT / "tableA_paired_vs_periscope.csv", index=False)
tabB.to_csv(OUT / "tableB_end_to_end_per_cell.csv", index=False)
sB.to_csv(OUT / "tableB_summary.csv")
pB.to_csv(OUT / "tableB_paired_vs_full.csv", index=False)
(OUT / "METHOD.json").write_text(json.dumps(dict(
    created_at_utc=pd.Timestamp.utcnow().isoformat(),
    metric="condition_aware_rho: Spearman across 4741 genes between predicted and observed PD-minus-control direction",
    seeds=SEEDS,
    ceiling=CEILING,
    ceiling_source="measured this session: Spearman-Brown reliability of the 3-PD-donor mean direction, then sqrt",
    why_two_tables="benchmark pins all methods to one frozen brain representation and decoder and swaps the blood-side map; "
                   "the ablation pipeline runs each model end to end through its own decoder. On matched axis and seed the "
                   "two differ by up to 0.0595, so they are not interchangeable.",
    tableA="benchmark harness, 8 axes, periscope vs 5 comparators, seed-matched paired",
    tableB="end-to-end, 6 axes (no pDC in the ablation run), full vs coupling_off vs random_pairing, seed-matched paired",
    open_item="to place coupling_off / random_pairing inside table A they must be run through the benchmark harness; "
              "their checkpoints exist for seeds 42-44 but the harness has not been run on them",
    caveats=["PeriScope is not best on MMD or mean_mse; the advantage is specific to condition_aware_rho",
             "with the condition label held fixed every method is near zero: blood alone carries no direction",
             "cellot has no seed-matched paired test in paired_comparisons.csv"],
), indent=2) + "\n")
print()
print("written:", OUT)
