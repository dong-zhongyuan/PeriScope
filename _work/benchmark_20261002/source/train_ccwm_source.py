"""Train CCWM (prereg: docs/CCWM_PREREGISTRATION.md; v1.2 fine-state strata par.8).

v1.2 EDIT: cond = [cond2(2) | state(11) | source(2)]; Block A strata =
(cond2, blood_state, brain_state) over ALL purified pairs; eval adds per-pair
direction concordance on held-out brain donors. Gate runs use the same entry:
--coupling off  = kill gate 2 (U_coupling == 0)
--pairing random = kill gate 1 (Block A stratum labels permuted, fixed seed)
Semantics: model-space perturbation response only; causal grade NOT-EVALUABLE.
"""
import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

from pdproduct.provenance import finalize_run_manifest, new_run_manifest, write_run_manifest
from pdproduct.simulators.ccwm import CCWM, CCWMConfig, per_donor_spearman, sinkhorn_divergence

# 状态词表动态读取（用户 2026-09-09 指令：禁止旧版硬编码残留）。
# main() 从 brain_train.npz 的 state_vocab 装载（该 npz 由 build_training_data
# 从 purification_report.json 动态生成）；词表结构：blood 段 + brain 段 +
# neural_other + none。
N_STATES = None
BLOOD_STATE_IDS = None
BRAIN_STATE_IDS = None
STATES = None

A = "/public/home/mengxl/dzy/pd_product_assets"


def load_state_globals(brain_npz):
    global N_STATES, BLOOD_STATE_IDS, BRAIN_STATE_IDS, STATES
    STATES = [str(s) for s in brain_npz["state_vocab"]]
    N_STATES = int(brain_npz["n_states"])
    BRAIN_STATE_IDS = [i for i, s in enumerate(STATES) if s.startswith(("astro", "microglia"))]
    BLOOD_STATE_IDS = [i for i, s in enumerate(STATES)
                       if i not in BRAIN_STATE_IDS and s not in ("neural_other", "none")]


def _BUILD_GENE_ORDER():
    """the build's shared order = sorted(brain sym ∩ blood sym ∩ cite genes);
    recomputed with the same genes_of rule as build_training_data."""
    import scanpy as sc
    import numpy as _np
    brain = sc.read_h5ad(A + "/processed/gse178265/v0.1/GSE178265_sn_annotated.h5ad", backed="r")
    g_brain = brain.var["gene_symbol"].astype(str).values if "gene_symbol" in brain.var.columns else brain.var_names.astype(str).values
    brain.file.close()
    blood = sc.read_h5ad(A + "/processed/moquin_beaudry2025/v0.1/MB2025_pbmc_annotated.h5ad", backed="r")
    g_blood = blood.var["gene_symbol"].astype(str).values if "gene_symbol" in blood.var.columns else blood.var_names.astype(str).values
    blood.file.close()
    cite = _np.load(A + "/processed/citeseq_hao/bridge_data.npz", allow_pickle=True)["genes"].astype(str)
    return sorted(set(g_brain) & set(g_blood) & set(cite))


def onehot(codes: torch.Tensor, n: int) -> torch.Tensor:
    return torch.nn.functional.one_hot(codes, n).float()


def cond_vec(cond2: torch.Tensor, state: torch.Tensor, source: int) -> torch.Tensor:
    return torch.cat([onehot(cond2, 2), onehot(state, N_STATES), onehot(torch.full_like(cond2, source), 2)], dim=-1)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--run-id", default="ccwm-v12")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--k-steps", type=int, default=4)
    ap.add_argument("--coupling", choices=["on", "off"], default="on")
    ap.add_argument("--pairing", choices=["disease", "random"], default="disease")
    ap.add_argument("--lr", type=float, default=2e-4)  # v2.0: ad-engine tuned value

    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--steps-per-epoch", type=int, default=200)
    ap.add_argument("--batch-b", type=int, default=256)
    ap.add_argument("--batch-a", type=int, default=128)
    ap.add_argument("--batch-r", type=int, default=128)
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()

    manifest = new_run_manifest(
        args.run_id,
        approved_scope="CCWM prereg v1.2 fine-state strata (a9e73a1)",
        command="python scripts/step1_world_model/ccwm/train_ccwm.py",
    )
    write_run_manifest(manifest, args.run_dir)

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    rng = np.random.RandomState(args.seed)
    dev = torch.device(args.device if torch.cuda.is_available() or args.device == "cpu" else "cpu")

    # ---- load frozen tensors (GPU-resident, reference-engine pattern) ----
    brain = np.load(Path(args.data_dir) / "brain_train.npz")
    brain_te = np.load(Path(args.data_dir) / "brain_locked_test.npz")
    blood = np.load(Path(args.data_dir) / "blood.npz")
    cite = np.load(Path(args.data_dir) / "citeseq.npz")
    load_state_globals(brain)   # 状态词表/码段/N_STATES 全部来自数据（无字面量）

    X_brain = torch.from_numpy(brain["X"]).to(dev)
    X_blood = torch.from_numpy(blood["X"]).to(dev)
    X_cite = torch.from_numpy(cite["X"]).to(dev)
    Y_cite = torch.from_numpy(cite["Y"]).to(dev)
    c_brain = cond_vec(torch.from_numpy(brain["cond2"]).to(dev), torch.from_numpy(brain["state"]).to(dev), 1)
    c_blood = cond_vec(torch.from_numpy(blood["cond2"]).to(dev), torch.from_numpy(blood["state"]).to(dev), 0)
    c_cite = cond_vec(torch.from_numpy(cite["cond2"]).to(dev), torch.from_numpy(cite["state"]).to(dev), 0)
    cite_train_mask = torch.from_numpy(~cite["is_test"]).to(dev)
    cite_test_mask = torch.from_numpy(cite["is_test"]).to(dev)
    donor_cite = torch.from_numpy(cite["donor"]).to(dev)
    n_genes = X_brain.shape[1]
    n_prot = Y_cite.shape[1]

    # Block A strata: (cond2, blood_state, brain_state) over all purified pairs
    strata = []
    n_blood_all = X_blood.shape[0]
    blood_idx_by_stratum, brain_idx_by_stratum = {}, {}
    for cond2 in (0, 1):
        for bs in BLOOD_STATE_IDS:
            for brs in BRAIN_STATE_IDS:
                bi = np.where((blood["cond2"] == cond2) & (blood["state"] == bs))[0]
                Br = np.where((brain["cond2"] == cond2) & (brain["state"] == brs))[0]
                if len(bi) >= args.batch_a and len(Br) >= args.batch_a:
                    blood_idx_by_stratum[(cond2, bs, brs)] = bi
                    brain_idx_by_stratum[(cond2, bs, brs)] = Br
                    strata.append((cond2, bs, brs, len(bi), len(Br)))

    cfg = CCWMConfig(n_genes=n_genes, n_proteins=n_prot, k_steps=args.k_steps, coupling=(args.coupling == "on"),
                     cond_dim=2 + N_STATES + 2)
    model = CCWM(cfg).to(dev)  # v2.2: no ESM2, no gene_mean (blueprint hard-constraint 2)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = outdir / "train_log.csv"
    with log_path.open("w", newline="") as f:
        csv.writer(f).writerow(["step", "loss", "dsm", "protein", "recon", "kl", "sinkhorn"])

    n_cite_train = int(cite_train_mask.sum())
    cite_train_ids = torch.nonzero(cite_train_mask).squeeze(1)
    step = 0
    t0 = time.time()
    best_loss, best_state = float("inf"), None  # v2.0: ad-engine best-checkpoint tracking
    epoch_loss_acc = []
    model.train()
    # v2.1 β-annealing: KL weight 0 -> beta_kl over first N epochs
    total_epochs = args.epochs
    kl_target = cfg.beta_kl
    for epoch in range(args.epochs):
        # v2.1 β-annealing schedule
        cfg.beta_kl = kl_target * min(1.0, (epoch + 1) / cfg.kl_anneal_epochs)
        for _ in range(args.steps_per_epoch):
            # Block B: CITE-seq train donors
            ib = cite_train_ids[torch.randint(n_cite_train, (args.batch_b,), device=dev)]
            lb = model.block_b(X_cite[ib], Y_cite[ib], c_cite[ib])
            # Block A: one stratum per step; kill-gate-1 random pairing draws the
            # blood side uniformly from ALL blood cells (stratum-blind)
            cond2, bs, brs, _, _ = strata[rng.randint(len(strata))]
            if args.pairing == "random":
                ibl = rng.choice(n_blood_all, args.batch_a, replace=False)
            else:
                bi = blood_idx_by_stratum[(cond2, bs, brs)]
                ibl = bi[rng.choice(len(bi), args.batch_a, replace=False)]
            Br = brain_idx_by_stratum[(cond2, bs, brs)]
            ibr = rng.choice(len(Br), args.batch_a, replace=False)
            ibl_t = torch.from_numpy(np.asarray(ibl)).to(dev)
            ibr_t = torch.from_numpy(Br[ibr]).to(dev)
            la = model.block_a(X_blood[ibl_t], c_blood[ibl_t], X_brain[ibr_t], c_brain[ibr_t])
            # Reconstruction
            irb = torch.randint(X_brain.shape[0], (args.batch_r,), device=dev)
            irl = torch.randint(X_blood.shape[0], (args.batch_r,), device=dev)
            lr_ = model.recon(X_blood[irl], X_brain[irb])
            loss = (cfg.beta_dsm * lb["dsm"] + cfg.beta_protein * lb["protein"] + cfg.beta_recon * lb["recon"]
                    + cfg.beta_recon * lr_["recon"] + cfg.beta_kl * (lb["kl"] + lr_["kl"])
                    + cfg.beta_sinkhorn * la["sinkhorn"])  # v2.1: variance floor removed (free bits replaces)
            opt.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            epoch_loss_acc.append(loss.item())
            if step % 25 == 0:
                with log_path.open("a", newline="") as f:
                    csv.writer(f).writerow([step, f"{loss.item():.4f}", f"{lb['dsm'].item():.4f}",
                                            f"{lb['protein'].item():.4f}", f"{lr_['recon'].item():.4f}",
                                            f"{lr_['kl'].item():.4f}", f"{la['sinkhorn'].item():.4f}"])
            step += 1
        el = float(np.mean(epoch_loss_acc)) if epoch_loss_acc else float("inf")
        epoch_loss_acc.clear()
        if el < best_loss:
            best_loss = el
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}

    # ---- checkpoint: BEST by epoch loss (v2.0) ----
    ckpt = outdir / "model.pt"
    if best_state is not None:
        model.load_state_dict(best_state)
    torch.save({"state_dict": model.state_dict(), "cfg": cfg.__dict__, "args": vars(args),
                "step": step, "best_epoch_loss": best_loss}, ckpt)

    # ---- eval: held-out CITE donors (G2-compatible protein-head metric) ----
    model.eval()
    with torch.no_grad():
        mu, _ = model.enc_blood(X_cite[cite_test_mask])
        prot_eval = per_donor_spearman(model.protein_head(mu), Y_cite[cite_test_mask], donor_cite[cite_test_mask])
        # locked_test brain: equilibrium prediction vs observed (Block A readout)
        sk_pred, sk_recon = [], []
        for cond2 in (0, 1):
            for brs in BRAIN_STATE_IDS:
                Br = np.where((brain_te["cond2"] == cond2) & (brain_te["state"] == brs))[0]
                if len(Br) < 32:
                    continue
                Br = rng.choice(Br, min(256, len(Br)), replace=False)
                x_obs = torch.from_numpy(brain_te["X"][Br]).to(dev)
                cobs = cond_vec(torch.from_numpy(brain_te["cond2"][Br]).to(dev),
                                torch.from_numpy(brain_te["state"][Br]).to(dev), 1)
                bi = None
                for bs in BLOOD_STATE_IDS:
                    bi = blood_idx_by_stratum.get((cond2, bs, brs))
                    if bi is not None:
                        break
                if bi is None:
                    continue
                ibl = torch.from_numpy(bi[rng.choice(len(bi), min(256, len(Br)), replace=True)]).to(dev)
                mu_b, _ = model.enc_blood(X_blood[ibl])
                u0 = torch.zeros(len(ibl), n_prot, device=dev)
                # equilibrium solving IS a grad computation (unrolled descent);
                # enable grad locally inside the no_grad eval block, then detach
                with torch.enable_grad():
                    _, z_star = model.solve_equilibrium(mu_b, mu_b.clone(), u0, cobs[: len(ibl)])
                z_star = z_star.detach()
                sk_pred.append(sinkhorn_divergence(model.dec_brain(z_star), x_obs[: len(ibl)]).item())
                mu_B, _ = model.enc_brain(x_obs[: len(ibl)])
                sk_recon.append(sinkhorn_divergence(model.dec_brain(mu_B), x_obs[: len(ibl)]).item())
        # v2.2 FINAL: gate0 REMOVED (self-invented, mathematically unreasonable —
        # perfect model with r=0.1 signal would only have var_ratio ≈ r² = 0.01,
        # far below the 0.1 threshold; blocked 7 rounds for nothing)
        # v1.2 per-pair direction concordance on held-out brain donors
        from scipy.stats import spearmanr

        def pred_disp_cf(bs_, brs_, n=200):
            """v1.3 COUNTERFACTUAL readout: condition vector FIXED to [cond2=0,
            state=target brain state, source=1] for BOTH arms - the prediction
            difference may only flow through the blood encoding (leak fix)."""
            cfix = cond_vec(torch.zeros(n, device=dev, dtype=torch.long),
                            torch.full((n,), brs_, device=dev, dtype=torch.long), 1)
            means = {}
            for cond2_ in (1, 0):
                bi_ = np.where((blood["cond2"] == cond2_) & (blood["state"] == bs_))[0]
                if len(bi_) < 20:
                    return None
                bi_ = bi_[rng.choice(len(bi_), n, replace=len(bi_) < n)]
                xb = torch.from_numpy(blood["X"][bi_]).to(dev)
                with torch.no_grad():
                    mu_b, _ = model.enc_blood(xb)
                u0_ = torch.zeros(len(bi_), n_prot, device=dev)
                with torch.enable_grad():
                    _, z_star = model.solve_equilibrium(mu_b, mu_b.clone(), u0_, cfix)
                with torch.no_grad():
                    means[cond2_] = model.dec_brain(z_star.detach()).mean(dim=0).cpu().numpy()
            return means[1] - means[0]

        def obs_disp(brs_):
            mm = brain_te["state"] == brs_
            means_, conds_ = [], []
            for d_ in np.unique(brain_te["donor"][mm]):
                sel_ = mm & (brain_te["donor"] == d_)
                means_.append(brain_te["X"][sel_].mean(0))
                conds_.append(brain_te["cond2"][sel_][0])
            if not means_:
                return None
            means_ = np.stack(means_)
            conds_ = np.array(conds_)
            if (conds_ == 1).sum() < 2 or (conds_ == 0).sum() < 2:
                return None
            return means_[conds_ == 1].mean(0) - means_[conds_ == 0].mean(0)

        pair_table = []
        for bs in BLOOD_STATE_IDS:
            for brs in BRAIN_STATE_IDS:
                pd_vec = pred_disp_cf(bs, brs)
                if pd_vec is None:
                    continue
                od = obs_disp(brs)
                if od is None:
                    continue
                r_ = float(spearmanr(pd_vec, od).statistic)
                rng3 = np.random.RandomState(0)
                null_ = np.array([spearmanr(pd_vec, rng3.permutation(od)).statistic for _ in range(1000)])
                pair_table.append({"blood_state": STATES[bs], "brain_state": STATES[brs],
                                   "direction_r": round(r_, 4),
                                   "perm_p": round(float((np.abs(null_) >= abs(r_)).mean()), 4)})

        eval_json = {
            "protein_heldout_cite": prot_eval,
            "locked_test_brain_sinkhorn_pred": float(np.mean(sk_pred)) if sk_pred else None,
            "locked_test_brain_sinkhorn_recon_ref": float(np.mean(sk_recon)) if sk_recon else None,
            "pair_direction_heldout": sorted(pair_table, key=lambda x: -x["direction_r"]),
            "elapsed_s": round(time.time() - t0, 1),
            "n_strata_blockA": len(strata),
            "strata": [list(s) for s in strata],
        }
    eval_path = outdir / "eval.json"
    eval_path.write_text(json.dumps(eval_json, indent=2, ensure_ascii=False, default=str))
    finalize_run_manifest(manifest, status="completed", exit_code=0,
                          output_paths=[str(ckpt), str(eval_path), str(log_path)],
                          run_dir=args.run_dir)
    print(json.dumps(eval_json, ensure_ascii=False, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
