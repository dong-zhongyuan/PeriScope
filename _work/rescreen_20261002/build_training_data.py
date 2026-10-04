"""Build repaired CCWM tensors from the existing purified cell-state maps.

Brain GSE178265, blood GSE223138 and MB2025, and the rebuilt raw-count Hao CITE
bridge share log1p(CP10k) RNA features. Centers use training donors only. Brain
model-selection, calibration and locked-test donor identities stay separate;
blood and CITE split masks are explicit. Cell-state vocabularies come from the
purification report, not old literal state codes.

Outputs: brain_train, brain_calibration, brain_locked_test, blood and citeseq
NPZ files, together with gene order, center vectors, donor metadata and report.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from input_utils import metadata, symbols, load_expression, BRAIN_VAL, BLOOD_VAL, BLOOD_TEST, donor_mean

from pdproduct.provenance import finalize_run_manifest, new_run_manifest, write_run_manifest

# 状态词表唯一真源 = purification_report.json（用户 2026-09-09 指令：全部动态读取，
# 禁止旧版硬编码残留）。main() 里 load_state_vocab() 填充 VOCAB；此前此处的
# BLOOD_STATES/BRAIN_STATES 字面量与纯化产物已不一致（多 inflammatory_mono / dam）。
VOCAB = None   # dict: blood_states / brain_states / state_vocab / n_states / codes
PURIFY_DIR = ""


def load_state_vocab(purify_dir):
    """Derive the combined state vocab from the purification report.

    codes: blood 0..len(blood)-1, brain follows (astro first, then microglia),
    then neural_other, then none. Every consumer reads the vocab from the npz
    this script writes — no downstream literals.
    """
    rep = json.loads((Path(purify_dir) / "purification_report.json").read_text())
    blood = sorted({s for k, r in rep.items() if k.endswith("blood_myeloid") for s in r["state_counts"]})
    astro = sorted({s for k, r in rep.items() if k.endswith("brain_astro") for s in r["state_counts"]})
    micro = sorted({s for k, r in rep.items() if k.endswith("brain_microglia") for s in r["state_counts"]
                    if s not in astro})
    brain = astro + micro
    return {
        "blood_states": blood,
        "brain_states": brain,
        "state_vocab": blood + brain + ["neural_other", "none"],
        "n_states": len(blood) + len(brain) + 2,
        "neural_other": len(blood) + len(brain),
        "state_none": len(blood) + len(brain) + 1,
        "brain_offset": len(blood),
    }


def state_from_map(obs_names, cohort, comp_family):
    """state codes for all cells of one h5ad from the purification mapping."""
    df = pd.read_csv(Path(PURIFY_DIR) / f"{cohort}_{comp_family}.csv", index_col=0)
    st = pd.Series(VOCAB["state_none"], index=obs_names, dtype=np.int64)
    vocab = VOCAB["blood_states"] if comp_family == "blood_myeloid" else VOCAB["brain_states"]
    kept = df[df["kept"]]
    offset = 0 if comp_family == "blood_myeloid" else VOCAB["brain_offset"]
    mapped = kept["state_pure"].astype(str).map({s: i + offset for i, s in enumerate(vocab)})
    common = st.index.intersection(mapped.index)  # map may cover cells filtered out upstream
    st.loc[common] = mapped.loc[common].astype(np.int64)
    return st.values


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def densify(x) -> np.ndarray:
    """h5ad layers may be scipy sparse; np.asarray alone yields a 0-d object array."""
    import scipy.sparse as sp

    return x.toarray() if sp.issparse(x) else np.asarray(x)


def gene_symbols(adata) -> np.ndarray:
    if "gene_symbol" in adata.var.columns:
        return adata.var["gene_symbol"].astype(str).values
    return adata.var_names.astype(str).values


def block_npz(X: np.ndarray, cond2: np.ndarray, comp: np.ndarray, ct: np.ndarray, donor: np.ndarray,
              mu: np.ndarray, sd: np.ndarray, names: dict) -> dict:
    return dict(X=X.astype(np.float32), cond2=cond2.astype(np.int64), comp=comp.astype(np.int64),
                ct=ct.astype(np.int64), donor=donor.astype(np.int64), gene_mu=mu, gene_sd=sd, **names)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--brain-h5ad", required=True)
    ap.add_argument("--split-json", required=True)
    ap.add_argument("--blood-h5ad", required=True)
    ap.add_argument("--blood2-h5ad", required=True, help="MB2025 (pooled per prereg v1.2)")
    ap.add_argument("--purify-dir", required=True)
    ap.add_argument("--cite-npz", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--run-dir", required=True)
    args = ap.parse_args()
    global PURIFY_DIR
    PURIFY_DIR = args.purify_dir
    global VOCAB
    VOCAB = load_state_vocab(args.purify_dir)

    manifest = new_run_manifest(
        "ccwm-build-training-data-v12",
        approved_scope="CCWM prereg v1.2 fine-state strata (a9e73a1)",
        command="python scripts/step1_world_model/ccwm/build_training_data.py",
    )
    write_run_manifest(manifest, args.run_dir)

    inputs = {k: {"path": v, "sha256": sha256_file(v)} for k, v in
              {"brain": args.brain_h5ad, "split": args.split_json, "blood223": args.blood_h5ad,
               "bloodmb": args.blood2_h5ad, "cite": args.cite_npz}.items()}

    split = json.loads(Path(args.split_json).read_text())
    donor_split = split["donors"]
    brain_donors_train = sorted(d for d, r in donor_split.items() if r in ("train","model_selection"))
    brain_donors_test = sorted(d for d, r in donor_split.items() if r == "locked_test")

    cite = np.load(args.cite_npz, allow_pickle=True)
    g_brain = symbols(metadata(args.brain_h5ad)[1])
    g_blood = symbols(metadata(args.blood_h5ad)[1])
    g_blood2 = symbols(metadata(args.blood2_h5ad)[1])
    g_cite = cite['genes'].astype(str)
    shared = sorted(set(g_brain) & set(g_blood) & set(g_blood2) & set(g_cite))
    print('shared genes',len(shared),'CR1 present', 'CR1' in shared,flush=True)
    brain = load_expression(args.brain_h5ad, shared)
    blood = load_expression(args.blood_h5ad, shared)
    blood2 = load_expression(args.blood2_h5ad, shared)
    idx_brain = idx_blood = idx_blood2 = np.arange(len(shared))
    idx_cite = np.array([np.where(g_cite == g)[0][0] for g in shared])

    def center_only(X: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """v2.2: center only, NO std division — preserves sparse expression structure
        that Zero-Masked MSE relies on (blueprint par.3.2 '中心化'; ad engine gene_mean)."""
        mu = X.mean(axis=0)
        return (X - mu).astype(np.float32), mu.astype(np.float32), np.zeros(X.shape[1], dtype=np.float32)

    out = {}
    # ---- brain (lognorm layer, split donors only; read column-subset once) ----
    brain_obs = brain.obs
    sub_brain = brain[:, idx_brain].copy()
    Xb_all = densify(sub_brain.layers["lognorm"])
    ct_all = brain_obs["cell_type"].values
    cond_all = brain_obs["condition"].values
    donor_all = brain_obs["donor"].values
    st_astro = state_from_map(brain_obs.index.astype(str), "c178", "brain_astro")
    st_mg = state_from_map(brain_obs.index.astype(str), "c178", "brain_microglia")
    st_brain = np.where(st_astro != VOCAB["state_none"], st_astro, st_mg)
    st_brain = np.where(st_brain == VOCAB["state_none"], VOCAB["neural_other"], st_brain)

    keep_train = np.isin(donor_all, brain_donors_train) & (ct_all != "unassigned")
    keep_test = np.isin(donor_all, brain_donors_test) & (ct_all != "unassigned")
    fit_brain = keep_train & ~np.isin(donor_all, BRAIN_VAL)
    mu = donor_mean(Xb_all, donor_all, fit_brain)
    sd = np.ones_like(mu)
    Xb_train = (Xb_all[keep_train] - mu).astype(np.float32)
    Xb_test = (Xb_all[keep_test] - mu).astype(np.float32)
    del Xb_all
    cond_map = {"Disease": 1, "PD": 1, "PDD": 1, "Ctrl": 0, "Control": 0, "Normal": 0, "Early PD": 1, "Late PD": 1, "CTRL": 0}
    donor_names_brain = sorted(set(donor_all))
    dcode = {d: i for i, d in enumerate(donor_names_brain)}
    ct_names_brain = sorted(set(ct_all[keep_train]) | set(ct_all[keep_test]))
    ccode_b = {c: i for i, c in enumerate(ct_names_brain)}
    out["brain_train.npz"] = dict(
        X=Xb_train,
        is_val=np.isin(donor_all[keep_train],BRAIN_VAL),
        barcodes=brain_obs.index.astype(str).to_numpy()[keep_train],
        donor_names=np.array(donor_names_brain,dtype=str),
        cond2=np.array([cond_map[c] for c in cond_all[keep_train]]),
        state=st_brain[keep_train].astype(np.int64),
        ct=np.array([ccode_b[c] for c in ct_all[keep_train]]),
        donor=np.array([dcode[d] for d in donor_all[keep_train]]),
        gene_mu=mu, gene_sd=sd,
    )
    out["brain_locked_test.npz"] = dict(
        X=Xb_test,
        barcodes=brain_obs.index.astype(str).to_numpy()[keep_test],
        donor_names=np.array(donor_names_brain,dtype=str),
        cond2=np.array([cond_map[c] for c in cond_all[keep_test]]),
        state=st_brain[keep_test].astype(np.int64),
        ct=np.array([ccode_b[c] for c in ct_all[keep_test]]),
        donor=np.array([dcode[d] for d in donor_all[keep_test]]),
        gene_mu=mu, gene_sd=sd,
    )

    calibration = np.isin(donor_all,[d for d,r in donor_split.items() if r=='calibration']) & (ct_all!='unassigned')
    # Read calibration RNA separately; neither features nor training nor validation selection uses it.
    Xcal = densify(sub_brain.layers['lognorm'])[calibration]-mu
    out['brain_calibration.npz']=dict(X=Xcal.astype(np.float32),
        cond2=np.array([cond_map[c] for c in cond_all[calibration]]),
        state=st_brain[calibration].astype(np.int64),ct=np.array([ccode_b[c] for c in ct_all[calibration]]),
        donor=np.array([dcode[d] for d in donor_all[calibration]]),
        barcodes=brain_obs.index.astype(str).to_numpy()[calibration],
        donor_names=np.array(donor_names_brain,dtype=str),gene_mu=mu,gene_sd=sd)

    # ---- blood: two cohorts pooled, single z-score scale (prereg v1.2) ----
    Xs, conds, cts, donors, sts, cohorts, barcodes_bl = [], [], [], [], [], [], []
    for a, gidx, cname in ((blood, idx_blood, "b223"), (blood2, idx_blood2, "bmb")):
        sub = a[:, gidx].copy()
        # binary cond2 axis is PD-vs-Ctrl: MSA/PSP donors excluded (prereg v1.2)
        keepc = np.isin(sub.obs["condition"].astype(str).values,
                        ["Disease", "PD", "PDD", "Early PD", "Late PD", "Ctrl", "Control", "Normal", "CTRL"])
        sub = sub[keepc].copy()
        Xs.append(densify(sub.layers["lognorm"]))
        conds.append(sub.obs["condition"].astype(str).values)
        cts.append(sub.obs["cell_type"].astype(str).values)
        donors.append(np.array([cname+":"+d for d in sub.obs["donor"].astype(str)]))
        barcodes_bl.extend([cname+":"+b for b in sub.obs.index.astype(str)])
        sts.append(state_from_map(sub.obs.index.astype(str), cname, "blood_myeloid"))
        cohorts.append(np.full(sub.n_obs, 0 if cname == "b223" else 1, dtype=np.int64))

    Xbl = np.concatenate(Xs)
    del Xs
    cond_bl = np.concatenate(conds)
    ct_bl = np.concatenate(cts)
    donor_bl = np.concatenate(donors)
    st_bl = np.concatenate(sts)
    cohort_bl = np.concatenate(cohorts)
    # One blood-side center is used by both disease RNA and CITE-seq.
    # Every training donor contributes equal weight within each source.
    Xc = cite["X"][:, idx_cite].astype(np.float32)
    donors_c = cite["donor"].astype(str)
    fit_bl = ~np.isin(donor_bl, BLOOD_VAL + BLOOD_TEST)
    fit_c = np.isin(donors_c,["P1","P2","P3","P4"])
    means = [donor_mean(Xbl,donor_bl,fit_bl & (cohort_bl==i)) for i in (0,1)]
    means.append(donor_mean(Xc,donors_c,fit_c))
    mu_bl = np.stack(means).mean(0).astype(np.float32)
    sd_bl = np.ones_like(mu_bl)
    Xbl = (Xbl - mu_bl).astype(np.float32)
    ct_names_blood = sorted(set(ct_bl))
    ccode_bl = {c: i for i, c in enumerate(ct_names_blood)}
    donor_names_blood = sorted(set(donor_bl))
    dcode_bl = {d: i for i, d in enumerate(donor_names_blood)}
    out["blood.npz"] = dict(
        X=Xbl,
        barcodes=np.array(barcodes_bl,dtype=str),
        is_val=np.isin(donor_bl,BLOOD_VAL),
        is_test=np.isin(donor_bl,BLOOD_TEST),
        donor_names=np.array(donor_names_blood,dtype=str),
        cond2=np.array([cond_map[c] for c in cond_bl]),
        state=st_bl.astype(np.int64),
        ct=np.array([ccode_bl[c] for c in ct_bl]),
        donor=np.array([dcode_bl[d] for d in donor_bl]),
        cohort=cohort_bl,
        gene_mu=mu_bl, gene_sd=sd_bl,
    )

    # ---- CITE-seq (Hao): X/Y already scaled by upstream; subset genes to shared ----
    Xc = cite["X"][:, idx_cite].astype(np.float32)
    Yc = cite["Y"].astype(np.float32)
    Xc = Xc - mu_bl
    mu_c, sd_c = mu_bl, sd_bl
    ct_labels = [str(x) for x in cite["ct_labels"]]

    donors_c = cite["donor"].astype(str)
    donor_names_cite = sorted(set(donors_c))
    held_out = donor_names_cite[-3:]
    dcode_c = {d: i for i, d in enumerate(donor_names_cite)}
    is_test = np.isin(donors_c, held_out)
    out["citeseq.npz"] = dict(
        X=Xc, Y=Yc,
        barcodes=cite["barcodes"], subtype=cite["subtype"],
        is_val=(donors_c=="P5"),
        donor_names=np.array(donor_names_cite,dtype=str),
        proteins=cite["proteins"],
        adt_counts=cite["adt_counts"],
        cond2=np.zeros(len(Xc), dtype=np.int64),  # healthy donors -> Ctrl-like
        state=np.array([VOCAB["state_vocab"].index({"CD14 Mono":"classical_mono", "CD16 Mono":"nonclassical_mono", "cDC1":"cDC", "cDC2":"cDC", "ASDC":"cDC", "pDC":"pDC"}.get(t,"none")) for t in cite["subtype"]],dtype=np.int64),
        ct=cite["ct"].astype(np.int64),
        donor=np.array([dcode_c[d] for d in donors_c]),
        is_test=is_test,
        gene_mu=mu_c, gene_sd=sd_c,
    )

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    for name, arrays in out.items():
        # 每个npz自描述状态词表，下游一律从此读取（无字面量）
        arrays["genes"] = np.array(shared,dtype=str)
        arrays.setdefault("state_vocab", np.array(VOCAB["state_vocab"], dtype=str))
        arrays.setdefault("n_states", np.int64(VOCAB["n_states"]))
        arrays = {k:(v.astype(str) if isinstance(v,np.ndarray) and v.dtype.kind=="O" else v) for k,v in arrays.items()}
        np.savez(outdir / name, **arrays)

    report = {
        "run": "rescreen_20261002; amended input and training recipe, not original preregistration",
        "normalization": "raw counts summed per duplicate symbol; log1p(CP10k); shared blood/CITE donor-balanced training center; separate brain training center",
        "brain_validation_donors": BRAIN_VAL, "blood_validation_donors": BLOOD_VAL, "blood_test_donors": BLOOD_TEST,
        "cite_validation_donors": ["P5"],
        "inputs": inputs,
        "purify_dir": PURIFY_DIR,
        "shared_genes": len(shared),
        "gene_sources": {"brain": len(g_brain), "blood223": len(g_blood), "bloodmb": len(g_blood2), "cite": len(g_cite)},
        "brain_train": {"cells": int(Xb_train.shape[0]), "donors": brain_donors_train,
                        "state_counts": {int(k): int(v) for k, v in zip(*np.unique(st_brain[keep_train], return_counts=True))}},
        "brain_locked_test": {"cells": int(Xb_test.shape[0]), "donors": brain_donors_test,
                              "state_counts": {int(k): int(v) for k, v in zip(*np.unique(st_brain[keep_test], return_counts=True))}},
        "blood": {"cells": int(Xbl.shape[0]), "donors": donor_names_blood,
                  "state_counts": {int(k): int(v) for k, v in zip(*np.unique(st_bl, return_counts=True))},
                  "cohorts": {"b223": int((cohort_bl == 0).sum()), "bmb": int((cohort_bl == 1).sum())}},
        "cite": {"cells": len(Xc), "donors": donor_names_cite, "ct_labels": ct_labels, "held_out_donors": held_out,
                 "train_cells": int((~is_test & (donors_c!="P5")).sum()), "validation_cells": int((donors_c=="P5").sum()), "test_cells": int(is_test.sum())},
        "state_vocab": VOCAB["state_vocab"],
        "state_vocab_source": "purification_report.json (dynamic, 2026-09-09)",
    }
    (outdir / "build_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False, default=str))
    finalize_run_manifest(manifest, status="completed", exit_code=0,
                          output_paths=[str(outdir / n) for n in out] + [str(outdir / "build_report.json")],
                          run_dir=args.run_dir)
    print(json.dumps({k: report[k] for k in ("shared_genes", "brain_train", "blood", "cite")}, ensure_ascii=False, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
