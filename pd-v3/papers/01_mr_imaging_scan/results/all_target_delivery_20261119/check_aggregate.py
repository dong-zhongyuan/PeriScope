#!/usr/bin/env python
"""Independent check of the all-target delivery tables.

Does not import aggregate_all_target_results.py. Row counts and a random sample of values are
recomputed from the raw per-region files, the plan, the audit result files and unit receipts.
A small sample of selection-null values is re-simulated with separately written code and
different random numbers, and compared within Monte Carlo error. Read-only except for one
report file next to this script.
"""
import argparse, glob, gzip, hashlib, json, os, sys
import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
P = "/public/home/mengxl/dzy/pd-v3/"
SCAN = P + "papers/01_mr_imaging_scan/"
COLOC = SCAN + "results/all_target_coloc_20261118/"
REVERSE = SCAN + "results/all_target_reverse_20261118/"
REGIONAL = P + "data/all_target_followup_20261118/regional_all/"
ACCEPT = P + "data/all_target_followup_20261118/acceptance_20261118/"
FAST = P + "data/all_target_followup_20261118/fast_range_20261118/"
RESEARCH = P + "papers/03_world_model_integration/results/research_20261118/"
GRID_KEYS = [(1e-5, "1e5"), (1e-6, "1e6")]
NULL_ORDER = ["s1p000", "s1p025", "s1p030", "s1p050", "s1p100", "lambda_median_family", "lambda_median_each", "lambda_meanchi2_each", "realvar_signflip", "realvar_winsor4p5"]
REF = ["s1p025", "s1p030", "lambda_median_family", "lambda_median_each", "lambda_meanchi2_each", "realvar_signflip", "realvar_winsor4p5"]
checks = []


def draws_for(name):
    a = RESEARCH + "verify-integration-INT-1-a/"
    r = RESEARCH + "verify-integration-INT-1-r/own_maxnull_draws.npz"
    if name.startswith("s1p"):
        return np.load(RESEARCH + "integration/atlas/maxnull_draws.npy").astype(float) * float(name[1:].replace("p", "."))
    if name.startswith("lambda"):
        return np.load(a + "null_" + name + ".npy").astype(float)
    return np.load(r)[{"realvar_signflip": "V1_realvar_signflip", "realvar_winsor4p5": "V2_realvar_winsor4.5"}[name]].astype(float)


def any_by(hit, codes):
    u, inv = np.unique(codes, return_inverse=True)
    out = np.zeros((len(u),) + hit.shape[1:], bool)
    np.logical_or.at(out, inv, hit)
    return out


def my_clump(au, pairs, drop):
    ld = {}
    for a, b, r in zip(pairs.variant_a, pairs.variant_b, pairs.r):
        ld[(a, b)] = r
        ld[(b, a)] = r
    ok = au[au.ld_selection_status.isin(["retained_by_exposure_only_LD_pruning", "LD_pruned", "LD_not_estimable"]) & ~au.variant_key.isin(drop)]
    ok = ok.assign(neg=-ok.discovery_neglog10_p).sort_values(["neg", "variant_key"])
    chosen = []
    for row in ok.itertuples(index=False):
        clash = False
        for other in chosen:
            if row.chromosome == other.chromosome and abs(row.position_grch37 - other.position_grch37) <= 10_000_000:
                r = ld.get((row.variant_key, other.variant_key), float("nan"))
                clash = not (r == r) or r * r > 0.001
                if clash:
                    break
        if not clash:
            chosen.append(row)
    return {c.variant_key for c in chosen}


def my_mr(v):
    bx, sx, by, sy = (v[c].values.astype(float) for c in ("beta_imaging", "se_imaging", "protein_beta_aligned", "protein_se"))
    if len(v) == 1:
        b, se = by[0] / bx[0], abs(sy[0] / bx[0])
    else:
        w = 1 / sy ** 2
        b = (w * bx * by).sum() / (w * bx ** 2).sum()
        se = (w * bx ** 2).sum() ** -.5
    return b, se, 2 * norm.sf(abs(b / se))


def rec(name, passed, detail=""):
    checks.append(dict(check=name, passed=bool(passed), detail=str(detail)))
    print(("PASS " if passed else "FAIL ") + name + ("  " + str(detail) if detail != "" else ""), flush=True)


def near(a, b, tol=1e-9):
    a, b = float(a), float(b)
    return (np.isnan(a) and np.isnan(b)) or abs(a - b) <= tol


def logsum(x, axis=-1):
    top = np.max(x, axis=axis, keepdims=True)
    return np.squeeze(top, axis=axis) + np.log(np.sum(np.exp(x - top), axis=axis))


def coloc_h4(labf1, labf2, p1, p12):
    # labf2: draws x SNPs
    t1 = logsum(labf1)
    t2 = logsum(labf2, axis=1)
    t12 = logsum(labf1 + labf2, axis=1)
    sumlog = t1 + t2
    gap = np.minimum(t12 - sumlog, -1e-12)
    h3 = 2 * np.log(p1) + sumlog + np.log1p(-np.exp(gap))
    terms = np.column_stack([np.zeros(len(t2)), np.log(p1) + t1 + np.zeros(len(t2)), np.log(p1) + t2, h3, np.log(p12) + t12])
    return np.exp(terms[:, 4] - logsum(terms, axis=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, default=60)
    ap.add_argument("--resim", type=int, default=6)
    ap.add_argument("--seed", type=int, default=20261119)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    T = pd.read_csv(HERE + "/all_target_region_table.csv.gz", low_memory=False, dtype={"idp": str, "chromosome": str})
    plan = json.load(open(REGIONAL + "PLAN.json"))["regions"]
    rec("row count equals planned regions", len(T) == len(plan) == 1375, f"table {len(T)} plan {len(plan)}")
    rec("one row per region id", T.region_id.is_unique and set(T.region_id) == set(plan))
    cov = pd.read_csv(ACCEPT + "per_region_coverage.csv", keep_default_na=False, dtype={"idp": str, "chromosome": str}).set_index("region_id")
    TT = T.set_index("region_id")
    rec("complete receipts on disk for coloc and reverse", all(os.path.exists(COLOC + r + "/COMPLETE.json") and os.path.exists(REVERSE + r + "/COMPLETE.json") for r in plan))
    rec("model_linked count from plan", int(T.model_linked.sum()) == sum(bool(v["model_linked"]) for v in plan.values()), int(T.model_linked.sum()))
    rec("chrX count from plan", int(T.chrX.sum()) == sum(str(v["chromosome"]) == "23" for v in plan.values()), int(T.chrX.sum()))
    for c in ("protein_state", "imaging_state", "coloc_state", "reverse_state"):
        rec(f"{c} matches acceptance table for all rows", (TT[c].astype(str) == cov[c].reindex(TT.index).astype(str)).all())
    rec("no region labelled with a verdict column", not any(k in c.lower() for c in T.columns for k in ("supported", "significant", "verdict")))
    ids = sorted(plan)
    sample = list(rng.choice(ids, min(a.sample, len(ids)), replace=False))
    bad_c, bad_r, n_c, n_r = [], [], 0, 0
    for rid in sample:
        q = pd.read_csv(COLOC + rid + "/coloc_prior_sensitivity.csv")
        row = TT.loc[rid]
        for r in q.to_dict("records"):
            for h in range(5):
                col = f"coloc_w{int(r['window_bp']) // 1000}k_p12_{r['p12_requested']:.3e}_PPH{h}"
                n_c += 1
                if not near(row[col], r[f"PP.H{h}.abf"], 1e-12):
                    bad_c.append((rid, col))
        e = pd.read_csv(REVERSE + rid + "/reverse_direction_estimates.csv")
        for r in e.to_dict("records"):
            t = f"rev_{r['palindrome_policy']}_{r['imaging_effect_set']}"
            n_r += 1
            if str(row[t + "_status"]) != str(r["status"]) or not near(row[t + "_beta"], r["beta"], 1e-12) or not near(row[t + "_p"], r["p"], 1e-15):
                bad_r.append((rid, t))
        if int(row.rev_n_estimated) != int((e.status == "estimated").sum()):
            bad_r.append((rid, "count"))
    rec(f"coloc posteriors, {len(sample)} random regions", not bad_c, f"{n_c} values, {len(bad_c)} differ {bad_c[:3]}")
    rec(f"reverse status and estimates, {len(sample)} random regions", not bad_r, f"{n_r} rows, {len(bad_r)} differ {bad_r[:3]}")
    agg = {"not_estimated_no_usable_instrument": 0, "estimated": 0, "not_estimated_no_published_peak": 0}
    for rid in ids:
        for k, v in json.loads(cov.loc[rid, "reverse_status_counts"]).items():
            agg[k] = agg.get(k, 0) + int(v)
    acc = json.load(open(ACCEPT + "ACCEPTANCE.json"))["reverse_scientific_status_totals"]
    rec("reverse status totals in table equal ACCEPTANCE.json", int(T.rev_n_estimated.sum()) == acc["estimated"] and int(T.rev_n_no_usable_instrument.sum()) == acc["not_estimated_no_usable_instrument"] and int(T.rev_n_no_published_peak.sum()) == acc["not_estimated_no_published_peak"], {k: int(v) for k, v in acc.items()})
    # routes and audits from unit receipts and raw audit files
    bad = []
    for rid in sample:
        r = plan[rid]
        row = TT.loc[rid]
        pc = json.load(open(REGIONAL + "sources/" + r["accession"] + "/COMPLETE.json"))
        pr = pc.get("acquisition_route", "full_file")
        if row.protein_route != pr:
            bad.append((rid, "protein", row.protein_route, pr))
        if str(row.family) == "BIG40":
            ic = json.load(open(REGIONAL + f"imaging/BIG40/{row.idp}/COMPLETE.json"))
            ir = ic.get("acquisition_route", "full_file")
            if row.imaging_route != ir:
                bad.append((rid, "imaging", row.imaging_route, ir))
    rec(f"routes match unit COMPLETE.json, {len(sample)} random regions", not bad, bad[:3])
    res = {}
    n_lines = json.load(open(HERE + "/SOURCE_IDENTITIES.json"))["audit_result_lines"]
    for kind in ("protein", "imaging"):
        s = {}
        for f in (FAST + f"audit/results_{kind}.jsonl", FAST + f"audit_targeted_H4_20261119/results_{kind}.jsonl"):
            if os.path.exists(f):
                # the audit may still be appending: use the lines the table was built from
                for line in open(f).read().split("\n")[:n_lines.get(f, 0)]:
                    if line.strip():
                        d = json.loads(line)
                        s.setdefault(str(d["unit"]), []).append(d)
        res[kind] = s
    wrong = []
    for rid, row in TT.iterrows():
        for kind, unit, route, st in (("protein", row.accession, row.protein_route, row.protein_audit_status), ("imaging", row.idp, row.imaging_route, row.imaging_audit_status)):
            if route != "fast_range_v1_20261118":
                if str(st).startswith("range_"):
                    wrong.append((rid, kind, st))
                continue
            got = res[kind].get(str(unit), [])
            if not got and str(st) not in ("range_not_audited", "range_not_audited_in_targeted_queue"):
                wrong.append((rid, kind, st))
            if got and all(x.get("outputs_match") is True for x in got) and st != "range_audited_match":
                wrong.append((rid, kind, st))
            if got and any(x.get("outputs_match") is False for x in got) and st != "range_audited_mismatch":
                wrong.append((rid, kind, st))
    rec("audit status agrees with raw audit files for every range-route unit", not wrong, wrong[:3])
    # assay level values
    m = pd.read_csv(RESEARCH + "integration/atlas/assay_maxT.csv").set_index("exposure_id")
    q = pd.read_csv(RESEARCH + "mr_imaging/assay_level_calibration.csv.gz").set_index("exposure_id")
    bad = [r for r in sample if not near(TT.loc[r, "p_maxT"], m.loc[plan[r]["exposure_ids"][0], "p_maxT"], 1e-12)]
    bad2 = [r for r in sample if plan[r]["exposure_ids"][0] in q.index and not near(TT.loc[r, "assay_q"], q.loc[plan[r]["exposure_ids"][0], "assay_q"], 1e-12)]
    rec("p_maxT and assay_q match their source tables in the sample", not bad and not bad2, f"{len(bad)} {len(bad2)}")
    rec("assay_q empty exactly when exposure is absent from the calibration", all(pd.isna(TT.loc[r, "assay_q"]) == (plan[r]["exposure_ids"][0] not in q.index) for r in ids))
    # selection null, independent re-simulation
    keys = np.load(HERE + "/cache/panel_keys.npy", mmap_mode="r")
    man = json.load(open(HERE + "/cache/panel_manifest.json"))
    panel_ids = np.array(man["ids"])
    Z = np.memmap(HERE + "/cache/panel_z.f32", dtype=np.float32, mode="r", shape=tuple(man["shape"]))
    maxz = np.load(RESEARCH + "integration/atlas/maxnull_draws.npy").astype(float)
    done = T[T.selection_null_computable & T.selnull_w1000k_instrument_in_coloc_snps.astype(str).eq("True")]
    rec("selection p present exactly for computable regions", (T.selection_p_H4.notna() == T.selection_null_computable).all(), int(T.selection_null_computable.sum()))
    rec("not-computable regions all carry a reason", T.loc[~T.selection_null_computable, "selnull_reason"].astype(str).str.len().gt(3).all(), int((~T.selection_null_computable).sum()))
    pick = list(rng.choice(done.region_id.values, min(a.resim, len(done)), replace=False))
    high = done.sort_values("H4_w1000k_1e5", ascending=False).region_id.head(2).tolist()
    bad, lines = [], []
    for rid in dict.fromkeys(pick + high):
        row = TT.loc[rid]
        v = pd.read_csv(COLOC + rid + "/window_1000000_snp_posteriors.csv.gz", keep_default_na=False, usecols=["variant_key", "beta_x", "standard_error_x", "beta_y", "standard_error_y", "log_abf_protein"])
        k = v.variant_key.values.astype(keys.dtype)
        ix = np.searchsorted(keys, k)
        ix[ix >= len(keys)] = len(keys) - 1
        hit = np.asarray(keys[ix]) == k
        v = v[hit].reset_index(drop=True)
        ix = ix[hit]
        n_all = int(hit.size)
        p1 = 1e-4 if n_all * 1e-4 < 1 else 1 / (n_all + 1)
        use = panel_ids != str(row.idp) if str(row.family) == "BIG40" else np.ones(len(panel_ids), bool)
        M = np.nan_to_num(np.asarray(Z[ix])[:, use].astype(float)).T
        inst = f"{row.chromosome_label}:{int(row.center_grch37)}:" + "_".join(sorted([row.instrument_effect_allele, row.instrument_other_allele]))
        pos = int(np.where(v.variant_key.values == inst)[0][0])
        M = M[np.abs(M).max(axis=1) < 5.45]
        slope = M.T.dot(M[:, pos]) / M[:, pos].dot(M[:, pos])
        l1 = v.log_abf_protein.values.astype(float)
        se2 = v.standard_error_y.values.astype(float) ** 2
        shrink = 0.0225 / (0.0225 + se2)
        lab = lambda zz: 0.5 * (np.log1p(-shrink) + shrink * zz ** 2)
        zo = v.beta_y.values.astype(float) / v.standard_error_y.values.astype(float)
        H_obs = float(coloc_h4(l1, lab(zo)[None, :], p1, 1e-5)[0])
        R = np.random.default_rng(int(hashlib.md5(rid.encode()).hexdigest()[:8], 16))
        out = []
        for _ in range(12):
            jj = R.integers(0, len(M), 250)
            zs = R.choice(maxz, 250) * 1.025
            out.append(coloc_h4(l1, lab(M[jj] + np.outer(zs - M[jj, pos], slope)), p1, 1e-5))
        sim = np.concatenate(out)
        p = (1 + (sim >= H_obs).sum()) / (len(sim) + 1)
        tp, tm = float(row.selection_p_H4), float(row.expected_H4_under_selection)
        sep = np.sqrt(max(tp * (1 - tp), 1e-4) * (1 / len(sim) + 1 / 2000))
        sem = sim.std() * np.sqrt(1 / len(sim) + 1 / 2000)
        ok = abs(p - tp) <= 4 * sep + 2e-3 and abs(sim.mean() - tm) <= 4 * sem + 2e-3 and abs(H_obs - float(row.selnull_w1000k_1e5_H4_obs_used)) < 1e-8
        lines.append(f"{row.gene}: p {tp:.4f} vs {p:.4f}, mean {tm:.4f} vs {sim.mean():.4f}")
        if not ok:
            bad.append(rid)
    rec("independent re-simulation of selection p and expected H4", not bad, "; ".join(lines))
    # panel equals production z for the region's own phenotype
    bad = []
    pool = T[(T.family == "BIG40") & T.idp.astype(str).isin(panel_ids)].region_id.values
    for rid in rng.choice(pool, min(10, len(pool)), replace=False):
        v = pd.read_csv(COLOC + rid + "/shared_snps.csv.gz", usecols=["variant_key", "beta_y", "standard_error_y"], keep_default_na=False)
        k = v.variant_key.values.astype(keys.dtype)
        ix = np.searchsorted(keys, k)
        j = int(np.where(panel_ids == str(TT.loc[rid, "idp"]))[0][0])
        d = np.abs(np.abs(np.asarray(Z[ix])[:, j].astype(float)) - np.abs(v.beta_y.values.astype(float) / v.standard_error_y.values.astype(float))).max()
        if d > 1e-4:
            bad.append((rid, float(d)))
    rec("panel |z| equals production |z| for a region's own phenotype, 10 regions", not bad, bad[:3])
    # group table recount
    G = pd.read_csv(HERE + "/h4_vs_selection_by_model_route.csv")
    g = G[(G.group == "all") & (G.unit == "region") & (G.window == "w1000k") & (G.p12 == "p12_1e-5") & (G.null_set == "s1p025")].iloc[0]
    ok = T[T.selnull_w1000k_1e5_H4_obs_used.notna()]
    rec("group table count of H4>=0.8 recomputed from main table", int(g.observed_H4ge080) == int((ok.selnull_w1000k_1e5_H4_obs_used >= .8).sum()), int(g.observed_H4ge080))
    rec("group table expected count equals sum of per-region probabilities (within 1.0)", abs(float(g.expected_H4ge080_under_selection) - float(ok.selection_P_H4ge080.sum())) < 1.0, f"{g.expected_H4ge080_under_selection:.2f} vs {ok.selection_P_H4ge080.sum():.2f}")
    rp = SCAN + "results/big40_repro_20261119/region_repro.csv.gz"
    if os.path.exists(rp):
        r = pd.read_csv(rp, low_memory=False)
        rec("repro_available matches region_repro.csv.gz", set(T.region_id[T.repro_available]) == set(r.region_id) & set(T.region_id), len(r))
    else:
        rec("repro_available is false everywhere while the repro file is absent", not T.repro_available.any())
    # ---- counting units, from the plan and raw coloc files for every region
    pe = {r: v["exposure_ids"][0] for r, v in plan.items()}
    T_exp = T.region_id.map(pe)
    rec("exposure_id equals the plan for all rows", (T_exp == T.exposure_id).all())
    sz = T_exp.map(T_exp.value_counts())
    rec("n_regions_for_exposure recomputed from the plan", (sz == T.n_regions_for_exposure).all(), f"{T_exp.nunique()} exposures")
    pts = sorted({(23 if str(v["chromosome"]).lstrip("0") in ("23", "X") else int(v["chromosome"]), int(v["center_grch37"])) for v in plan.values()})
    lab, cur, prev = {}, 0, None
    for c_, pos_ in pts:
        if prev is not None and (c_ != prev[0] or pos_ - prev[1] > 1_000_000):
            cur += 1
        lab[(c_, pos_)] = cur
        prev = (c_, pos_)
    mine = T.region_id.map(lambda r: lab[(23 if str(plan[r]["chromosome"]).lstrip("0") in ("23", "X") else int(plan[r]["chromosome"]), int(plan[r]["center_grch37"]))])
    same = T.groupby("locus_id").apply(lambda x: mine[x.index].nunique()).eq(1).all() and mine.groupby(mine).apply(lambda x: T.locus_id[x.index].nunique()).eq(1).all()
    rec("locus_id is the same partition as my own 1 Mb chaining of plan positions", bool(same), f"{cur + 1} loci")
    raw = {}
    for rid in T.region_id:
        q = pd.read_csv(COLOC + rid + "/coloc_prior_sensitivity.csv", usecols=["window_bp", "p12_requested", "PP.H4.abf"])
        q = q[q.window_bp == 1000000]
        raw[rid] = float(q.iloc[(np.log(q.p12_requested) - np.log(1e-5)).abs().values.argmin()]["PP.H4.abf"])
    h4 = T.region_id.map(raw).values
    ml = T.model_linked.values
    bad = []
    for gname, m in (("all", np.ones(len(T), bool)), ("model_linked=True", ml), ("model_linked=False", ~ml)):
        for unit, codes in (("region", np.arange(len(T))), ("exposure", T.exposure_id.values), ("locus", T.locus_id.values)):
            o = int(any_by(h4[m] >= .8, codes[m]).sum())
            n = len(np.unique(codes[m]))
            r_ = G[(G.group == gname) & (G.unit == unit) & (G.window == "w1000k") & (G.p12 == "p12_1e-5") & (G.null_set == "s1p025")].iloc[0]
            if int(r_.observed_H4ge080) != o or int(r_.n_units) != n:
                bad.append((gname, unit, o, int(r_.observed_H4ge080), n, int(r_.n_units)))
    a_ = {u: int(any_by(h4 >= .8, c).sum()) for u, c in (("region", np.arange(len(T))), ("exposure", T.exposure_id.values), ("locus", T.locus_id.values))}
    rec("observed H4>=0.8 per unit from raw coloc files equals the group table", not bad, f"{a_} {bad[:3]}")
    # ---- expected counts per unit, recounted from the cached simulations with separate code
    bad, shown = [], []
    for nm in ("s1p025", "lambda_median_each", "realvar_signflip"):
        k = NULL_ORDER.index(nm)
        Hm = np.stack([np.load(HERE + f"/cache/null_v2/{r}.npy", mmap_mode="r")[0, 0, k] for r in T.region_id]) >= .8
        for gname, m in (("all", np.ones(len(T), bool)), ("model_linked=True", ml), ("model_linked=False", ~ml)):
            for unit, codes in (("region", np.arange(len(T))), ("exposure", T.exposure_id.values), ("locus", T.locus_id.values)):
                c_ = any_by(Hm[m], codes[m]).sum(0)
                o = int(any_by(h4[m] >= .8, codes[m]).sum())
                pv = ((c_ >= o).sum() + 1) / (len(c_) + 1)
                r_ = G[(G.group == gname) & (G.unit == unit) & (G.window == "w1000k") & (G.p12 == "p12_1e-5") & (G.null_set == nm)].iloc[0]
                if abs(c_.mean() - r_.expected_H4ge080_under_selection) > 1e-9 or abs(pv - r_.P_null_count_ge_observed_H4ge080) > 1e-12:
                    bad.append((nm, gname, unit))
                if gname == "all":
                    shown.append(f"{nm}/{unit} {c_.mean():.1f} P={pv:.4f}")
    rec("expected counts and P per unit recounted from cached simulations", not bad, "; ".join(shown) + f" | differ {bad[:3]}")
    # ---- range columns and regression against the first delivery
    bad = []
    for lab_, suf in (("selection_p_H4", "p"), ("expected_H4_under_selection", "expected_H4"), ("selection_P_H4ge080", "P_H4ge080")):
        m = np.column_stack([T[f"selnull_w1000k_1e5_{n}_{suf}"].values for n in REF])
        if np.abs(m.min(1) - T[lab_ + "_ref_min"].values).max() > 1e-12 or np.abs(m.max(1) - T[lab_ + "_ref_max"].values).max() > 1e-12:
            bad.append(lab_)
    rec("reference-range columns equal min and max over the seven reference nulls", not bad, bad)
    prev_f = HERE + "/cache/previous_v1/all_target_region_table.csv.gz"
    if os.path.exists(prev_f):
        cols = [f"selnull_{w}_{t}_{n}_{s}" for w in ("w1000k", "w500k") for t in ("1e5", "1e6") for n in ("s1p000", "s1p025", "s1p050", "s1p100") for s in ("p", "expected_H4", "P_H4ge080")]
        P0 = pd.read_csv(prev_f, usecols=["region_id"] + cols).set_index("region_id").reindex(T.region_id)
        d = max(float(np.nanmax(np.abs(P0[c].values - T[c].values))) for c in cols)
        rec("the four original scales are unchanged from the first delivery", d == 0, f"largest difference {d}")
    # ---- new nulls, re-simulated with separate code, including regions whose instrument is outside the coloc SNPs
    outside = T[T.selnull_w1000k_instrument_in_coloc_snps.astype(str).eq("False")].region_id.values
    inside = done.region_id.values
    pick = list(rng.choice(inside, 3, replace=False)) + list(rng.choice(outside, 3, replace=False))
    bad, lines = [], []
    pools = {n: draws_for(n) for n in ("lambda_median_family", "lambda_meanchi2_each", "realvar_signflip", "s1p030")}
    for rid in pick:
        row = TT.loc[rid]
        v = pd.read_csv(COLOC + rid + "/window_1000000_snp_posteriors.csv.gz", keep_default_na=False, usecols=["variant_key", "beta_y", "standard_error_y", "log_abf_protein"])
        k = v.variant_key.values.astype(keys.dtype)
        ix = np.searchsorted(keys, k)
        use = panel_ids != str(row.idp) if str(row.family) == "BIG40" else np.ones(len(panel_ids), bool)
        M = np.nan_to_num(np.asarray(Z[ix])[:, use].astype(float)).T
        inst = f"{row.chromosome_label}:{int(row.center_grch37)}:" + "_".join(sorted([row.instrument_effect_allele, row.instrument_other_allele]))
        w_ = np.where(v.variant_key.values == inst)[0]
        if len(w_):
            cs = M[:, w_[0]]
        else:
            j = int(np.searchsorted(keys, np.array([inst]).astype(keys.dtype))[0])
            cs = np.nan_to_num(np.asarray(Z[j])[use].astype(float))
        keep = (np.abs(M).max(1) < 5.45) & (np.abs(cs) < 5.45)
        M, cs = M[keep], cs[keep]
        slope = M.T.dot(cs) / cs.dot(cs)
        if len(w_):
            slope[w_[0]] = 1.0
        l1 = v.log_abf_protein.values.astype(float)
        se2 = v.standard_error_y.values.astype(float) ** 2
        shrink = 0.0225 / (0.0225 + se2)
        lab2 = lambda zz: 0.5 * (np.log1p(-shrink) + shrink * zz ** 2)
        p1 = 1e-4 if len(v) * 1e-4 < 1 else 1 / (len(v) + 1)
        R = np.random.default_rng(int(hashlib.md5(("chk" + rid).encode()).hexdigest()[:8], 16))
        for nm, pool in pools.items():
            sim = []
            for _ in range(8):
                jj = R.integers(0, len(M), 250)
                sim.append(coloc_h4(l1, lab2(M[jj] + np.outer(R.choice(pool, 250) - cs[jj], slope)), p1, 1e-5))
            sim = np.concatenate(sim)
            tm = float(row[f"selnull_w1000k_1e5_{nm}_expected_H4"])
            tf = float(row[f"selnull_w1000k_1e5_{nm}_P_H4ge080"])
            sem = max(sim.std(), 1e-3) * np.sqrt(1 / len(sim) + 1 / 2000)
            sef = np.sqrt(max(tf * (1 - tf), 1e-3) * (1 / len(sim) + 1 / 2000))
            if abs(sim.mean() - tm) > 4.5 * sem + 2e-3 or abs((sim >= .8).mean() - tf) > 4.5 * sef + 2e-3:
                bad.append((rid, nm, round(tm, 4), round(float(sim.mean()), 4)))
        lines.append(f"{row.gene}{'(outside)' if not len(w_) else ''}")
    rec("independent re-simulation of the added nulls, 3 regions inside and 3 outside the coloc SNPs", not bad, f"{lines} differ {bad[:4]}")
    # ---- decoy calibration
    Dp = pd.read_csv(HERE + "/decoy_pairs_this_simulator.csv.gz", dtype={"decoy": str})
    Bn = pd.read_csv(HERE + "/decoy_calibration_bins.csv")
    edges = [2.5, 3, 3.25, 3.5, 3.75, 4, 4.5, 5.45, np.inf]
    labs = ["2.5-3", "3-3.25", "3.25-3.5", "3.5-3.75", "3.75-4", "4-4.5", "4.5-5.45", "5.45-99"]
    b_ = pd.cut(Dp.z_inst, edges, right=False, labels=labs)
    bad = []
    for tag, pl_ in (("1e5", "p12_1e-5"), ("1e6", "p12_1e-6")):
        for th in (.5, .8, .9):
            x = Bn[(Bn.calibration_source == "this_simulator") & (Bn.p12 == pl_) & np.isclose(Bn.threshold, th)].set_index("z_bin")
            hits = (Dp[f"H4_{tag}"] >= th).groupby(b_, observed=False).sum()
            pred = Dp[f"P{int(th * 100)}_{tag}"].groupby(b_, observed=False).sum()
            if np.abs(hits.reindex(labs).values - x.real_hits.reindex(labs).values).max() > 0 or np.abs(pred.reindex(labs).values - x.model_expected.reindex(labs).values).max() > 1e-8:
                bad.append((tag, th))
    rec("calibration bin totals recomputed from the pair file", not bad, f"{len(Dp)} pairs {bad[:3]}")
    bad = []
    for r_ in Dp.sample(12, random_state=a.seed).itertuples(index=False):
        row = TT.loc[r_.region_id]
        v = pd.read_csv(COLOC + r_.region_id + "/window_1000000_snp_posteriors.csv.gz", keep_default_na=False, usecols=["variant_key", "standard_error_y", "log_abf_protein"])
        ix = np.searchsorted(keys, v.variant_key.values.astype(keys.dtype))
        zz = np.asarray(Z[ix])[:, int(np.where(panel_ids == r_.decoy)[0][0])].astype(float)
        se2 = v.standard_error_y.values.astype(float) ** 2
        sh = 0.0225 / (0.0225 + se2)
        p1 = 1e-4 if len(v) * 1e-4 < 1 else 1 / (len(v) + 1)
        got = float(coloc_h4(v.log_abf_protein.values.astype(float), (0.5 * (np.log1p(-sh) + sh * zz ** 2))[None, :], p1, 1e-5)[0])
        if abs(got - r_.H4_1e5) > 1e-6:
            bad.append((r_.region_id, r_.decoy, got, r_.H4_1e5))
    rec("real decoy H4 recomputed from the panel, 12 pairs", not bad, bad[:3])
    Cc = pd.read_csv(HERE + "/h4_vs_selection_decoy_calibrated.csv")
    x = Cc[(Cc.group == "all") & (Cc.unit == "region") & (Cc.p12 == "p12_1e-5") & np.isclose(Cc.threshold, .8) & (Cc.calibration_source == "this_simulator")]
    W = Bn[(Bn.calibration_source == "this_simulator") & (Bn.p12 == "p12_1e-5") & np.isclose(Bn.threshold, .8)].set_index("z_bin").weight_used.reindex(labs).values
    bad, shown = [], []
    ex_list = list(pd.unique(T.exposure_id))
    U = np.stack([np.random.default_rng(int.from_bytes(hashlib.sha256(("selnull|" + e).encode()).digest()[:8], "little")).random(2000) for e in ex_list])
    ei = pd.Series(range(len(ex_list)), index=ex_list).reindex(T.exposure_id).values
    for nm in ("s1p025", "realvar_signflip"):
        k = NULL_ORDER.index(nm)
        S_ = np.sort(draws_for(nm))
        zs = S_[(U * len(S_)).astype(int)]
        bn = np.clip(np.searchsorted(np.array(edges), zs, side="right") - 1, 0, 7)
        Hm = np.stack([np.load(HERE + f"/cache/null_v2/{r}.npy", mmap_mode="r")[0, 0, k] for r in T.region_id]) >= .8
        e_w = float((Hm * W[bn[ei]]).sum() / 2000)
        t_ = float(x[x.null_set == nm].expected_calibrated.iloc[0])
        shown.append(f"{nm}: weighted {e_w:.1f}, table {t_:.1f}")
        if abs(e_w - t_) > 4:
            bad.append(nm)
    rec("calibrated expected count agrees with a direct weighted sum (within 4 hits; the table uses random thinning)", not bad, "; ".join(shown))
    # ---- reverse MR without cis instruments, every region, own clump and own IVW
    n_est = n_cis_est = 0
    cis_regions = set()
    bad_n, bad_v, bad_t, bad_prod = [], [], [], []
    for rid, row in TT.iterrows():
        d = REVERSE + rid + "/"
        ra = pd.read_csv(d + "all_peak_reverse_ratios_and_coverage.csv", low_memory=False)
        au = pd.read_csv(d + "all_peaks_exposure_only_LD_audit.csv", low_memory=False)
        e = pd.read_csv(d + "reverse_direction_estimates.csv")
        chrom = 23 if str(row.chromosome_label) == "X" else int(row.chromosome_label)
        cen = int(row.center_grch37)
        def is_cis(x):
            return (pd.to_numeric(x.chromosome.astype(str).str.replace("X", "23"), errors="coerce") == chrom) & ((x.position_grch37 - cen).abs() <= 1_000_000)
        cis_k = set(au.variant_key[is_cis(au)]) if len(au) else set()
        kept_prod = set(au.variant_key[au.ld_selection_status == "retained_by_exposure_only_LD_pruning"]) if len(au) else set()
        kept_nc = kept_prod
        if cis_k & set(au.variant_key[au.ld_selection_status.isin(["retained_by_exposure_only_LD_pruning", "LD_pruned", "LD_not_estimable"])]):
            pr = pd.read_csv(d + "all_peak_pair_reference_LD.csv")
            if my_clump(au, pr, set()) != kept_prod:
                bad_prod.append(rid)
            kept_nc = my_clump(au, pr, cis_k)
        for r_ in e.itertuples(index=False):
            t = f"rev_{r_.palindrome_policy}_{r_.imaging_effect_set}"
            if r_.status != "estimated":
                if row[t + "_nocis_status"] != r_.status or int(row[t + "_n_cis_instruments"]) != 0:
                    bad_n.append((rid, t))
                continue
            n_est += 1
            g_ = ra[(ra.palindrome_policy == r_.palindrome_policy) & (ra.imaging_effect_set == r_.imaging_effect_set) & (ra.status == "estimated")]
            sel = g_[g_.selected_after_LD.astype(str) == "True"]
            nc = int(sel.variant_key.isin(cis_k).sum())
            n_cis_est += nc > 0
            if nc:
                cis_regions.add(rid)
            if nc != int(row[t + "_n_cis_instruments"]):
                bad_n.append((rid, t, nc))
            if str(row[t + "_evidence_tier"]) != ("B" if nc else "A"):
                bad_t.append((rid, t))
            b0, s0, p0 = my_mr(sel)
            if abs(b0 - r_.beta) > 1e-9 * max(1, abs(r_.beta)) or abs(p0 - r_.p) > 1e-9:
                bad_prod.append((rid, t))
            w_ = g_[g_.variant_key.isin(kept_nc)]
            if not len(w_):
                if row[t + "_nocis_status"] != "not_estimated_no_instrument_outside_cis_window":
                    bad_v.append((rid, t, "status"))
                continue
            b1, s1, p1_ = my_mr(w_)
            if row[t + "_nocis_status"] != "estimated" or int(row[t + "_nocis_n_instruments"]) != len(w_) or abs(b1 - float(row[t + "_nocis_beta"])) > 1e-9 * max(1, abs(b1)) or abs(s1 - float(row[t + "_nocis_se"])) > 1e-9 or abs(p1_ - float(row[t + "_nocis_p"])) > 1e-9:
                bad_v.append((rid, t))
    rec("reverse estimates recomputed with own IVW and Wald equal production, own clump reproduces production", not bad_prod, f"{n_est} estimates {bad_prod[:3]}")
    rec("cis instrument counts per estimate from raw peak positions", not bad_n, f"{n_cis_est} estimates with a cis instrument in {len(cis_regions)} regions {bad_n[:3]}")
    rec("reverse estimates without cis peaks recomputed with own clump and own IVW", not bad_v, bad_v[:3])
    rec("reverse evidence tier is B exactly when a cis instrument is used", not bad_t, bad_t[:3])
    bad = [c["check"] for c in checks if not c["passed"]]
    with open(HERE + "/check_aggregate_report.json", "w") as f:
        json.dump(dict(n_checks=len(checks), n_failed=len(bad), failed=bad, checks=checks, sample_size=a.sample, resimulated=a.resim + 2, seed=a.seed), f, indent=1)
    print(f"{len(checks) - len(bad)} of {len(checks)} checks passed")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
