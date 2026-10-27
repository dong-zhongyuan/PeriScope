"""Recompute NMI revision statistics from released inputs, without GPU inference.

Frozen results are never modified. --ma also reruns the 5,000-set Ma cohort
analysis from donor pseudobulk counts and may take several minutes.
"""
import argparse
import contextlib
import json
from pathlib import Path
import platform
import shutil
import sys

import numpy as np
import pandas as pd
import scipy
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "reproduce/inputs"
RESULTS = ROOT / "results/nmi_revision_20261017"
SCRIPTS = ROOT / "_work/nmi_revision_20261017"
METHODS = ["periscope", "shared_label_mean", "shared_label_marginal", "ot_ridge",
           "cellot", "conditional_icnn", "conditional_mlp"]
ROUTES = ["cDC", "classical_mono", "nonclassical_mono", "pDC"]
RECEIVERS = ["microglia_mhc2", "astro"]


def compare(path):
    actual = pd.read_csv(path)
    expected = pd.read_csv(RESULTS / path.name)
    pd.testing.assert_frame_equal(actual, expected, check_dtype=False,
                                  rtol=1e-7, atol=1e-10)
    return {"file": path.name, "rows": len(actual), "matches_archive": True}


def benchmark(out):
    rows, pooled, ensembles, matched = [], [], {}, []
    labels = dict(zip(METHODS, ["PeriScope", "label/state brain mean", "label/state brain marginal",
                              "OT-Ridge", "CellOT", "conditional ICNN", "conditional MLP"]))
    with np.load(INPUTS / "benchmark_contrasts.npz", allow_pickle=False) as z:
        for rc in RECEIVERS:
            observed = z[f"observed|{rc}"]
            for method in METHODS:
                seeds = [s for s in range(42, 52) if f"{method}|{s}|cDC|{rc}" in z]
                per_seed, per_route = [], {route: [] for route in ROUTES}
                for seed in seeds:
                    ds = [z[f"{method}|{seed}|{route}|{rc}"] for route in ROUTES]
                    per_seed.append(np.mean(ds, axis=0))
                    for route, d in zip(ROUTES, ds):
                        per_route[route].append(d)
                ens = np.mean(per_seed, axis=0)
                ensembles[method, rc] = ens
                rs = [stats.spearmanr(d, observed).statistic for d in per_seed
                      if np.ptp(d) > 1e-9]
                rows.append(dict(receiver=rc, method=method, n_seeds=len(seeds),
                    rho_equal_route_ensemble=stats.spearmanr(ens, observed).statistic,
                    rho_per_seed_mean=np.mean(rs), frac_positive=np.mean(np.array(rs) > 0),
                    **{f"rho_{route}": stats.spearmanr(np.mean(ds, axis=0), observed).statistic
                       for route, ds in per_route.items() if np.ptp(np.mean(ds, axis=0)) > 1e-9}))
        for method in METHODS:
            predicted = np.concatenate([ensembles[method, rc] for rc in RECEIVERS])
            observed = np.concatenate([z[f"observed|{rc}"] for rc in RECEIVERS])
            pooled.append(dict(method=method, n_seeds=10 if method == "periscope" else 3,
                               rho=stats.spearmanr(predicted, observed).statistic))
            primary = {rc: np.mean([np.mean([z[f"{method}|{s}|{rt}|{rc}"] for rt in ROUTES], axis=0)
                                    for s in [42, 43, 44]], axis=0) for rc in RECEIVERS}
            rs = {rc: stats.spearmanr(primary[rc], z[f"observed|{rc}"]).statistic for rc in RECEIVERS}
            matched.append(dict(method=labels[method], key=method, rho_MGAP=rs["microglia_mhc2"],
                                rho_Astro=rs["astro"], rho_pooled=stats.spearmanr(
                                    np.concatenate([primary[rc] for rc in RECEIVERS]), observed).statistic,
                                rho_min=min(rs.values())))
    pd.DataFrame(rows).to_csv(out / "mgap_allmethods.csv", index=False)
    pd.DataFrame(pooled).to_csv(out / "benchmark_pooled_all_available.csv", index=False)
    pd.DataFrame(matched).sort_values("rho_pooled", ascending=False).to_csv(out / "mgap_matched.csv", index=False)
    return [compare(out / "mgap_allmethods.csv"), compare(out / "mgap_matched.csv")]


def donor_dependence(out):
    metadata = pd.read_csv(INPUTS / "blood_donor_metadata.csv")
    held = metadata[metadata.split.isin(["val", "test"])]
    groups = [sorted(held.loc[held.blood_condition == c, "donor_id"].unique()) for c in [0, 1]]
    with np.load(RESULTS / "bloodonly_profiles.npz", allow_pickle=False) as z:
        prof = {}
        for key in z.files:
            s, bs, rs, cb, d = key.split("|")
            prof[int(s), bs, rs, int(cb), int(d)] = z[key]
    seeds = sorted({key[0] for key in prof})
    labels = {"astro": "Astro", "microglia_mhc2": "MG-AP"}
    route_labels = dict(zip(ROUTES, ["cDC", "cMono", "nMono", "pDC"]))
    rows = []
    for rc in ["astro", "microglia_mhc2"]:
        for route in ROUTES:
            profiles = [[np.mean([prof[s, route, rc, 0, d] for s in seeds
                                      if (s, route, rc, 0, d) in prof], axis=0)
                         for d in donors if any((s, route, rc, 0, d) in prof for s in seeds)]
                        for donors in groups]
            c0, c1 = profiles
            row = dict(receiver=labels[rc], route=route_labels[route], n_control=len(c0), n_pd=len(c1))
            if len(c0) + len(c1) < 3:
                row.update(between=np.nan, within=np.nan, ratio=np.nan)
            else:
                c0, c1 = np.array(c0), np.array(c1)
                allp = np.concatenate([c0, c1])
                between = float(np.sqrt(((c1.mean(0) - c0.mean(0)) ** 2).mean()))
                within = float(np.sqrt(np.mean([((x - allp.mean(0)) ** 2).mean() for x in allp])))
                row.update(between=between, within=within, ratio=within / between)
            rows.append(row)
    pd.DataFrame(rows).to_csv(out / "donor_dependence.csv", index=False)
    return [compare(out / "donor_dependence.csv")]


def run_archived(name, out):
    """Run the original scientific calculation; substitute file locations only."""
    path = SCRIPTS / name
    code = path.read_text()
    code = code.replace("/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product", str(ROOT))
    code = code.replace('W + "/diag/', repr(str(RESULTS)) + ' + "/')
    old_argv = sys.argv
    try:
        sys.argv = [str(path), str(out)]
        with (out / (path.stem + ".log")).open("w") as log, contextlib.redirect_stdout(log):
            exec(compile(code, str(path), "exec"), {"__name__": "__main__", "__file__": str(path)})
    finally:
        sys.argv = old_argv
    return compare(out / (path.stem + ".csv"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path, help="new output directory")
    parser.add_argument("--ma", action="store_true", help="also rerun the Ma donor-level null tests")
    args = parser.parse_args()
    out = args.output.expanduser().resolve()
    if out.exists():
        parser.error("output directory already exists; choose a new directory")
    if out == ROOT / "results" or ROOT / "results" in out.parents:
        parser.error("output must not be inside the frozen results directory")
    out.mkdir(parents=True)
    checks = []
    for fn in [benchmark, donor_dependence]:
        for check in fn(out):
            checks.append(check)
            print("PASS:", check, flush=True)
    for name in ["mr_steiger.py", "ma_spatial_complement.py"]:
        check = run_archived(name, out)
        checks.append(check)
        print("PASS:", check, flush=True)
    if args.ma:
        shutil.copytree(RESULTS / "ma", out / "ma")
        print("Recomputing Ma cohort from pseudobulk counts (5,000 matched null sets)...", flush=True)
        check = run_archived("ma_cohort_test.py", out)
        checks.append(check)
        print("PASS:", check, flush=True)
    report = dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                  scipy=scipy.__version__, checks=checks, ma_recomputed=args.ma,
                  scope="Statistics from frozen predictions, harmonized MR and donor pseudobulk inputs; no training or fresh spatial deconvolution.")
    (out / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print("Primary matched comparison (all methods: benchmark seeds 42–44):")
    print(pd.read_csv(out / "mgap_matched.csv").to_string(index=False))


if __name__ == "__main__":
    main()
