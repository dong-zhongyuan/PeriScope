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
from ccwm import CCWM, CCWMConfig, per_donor_spearman, sinkhorn_divergence

from training_utils import BalancedSampler, validation

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
    ap.add_argument("--sample-latent-a",action="store_true")
    ap.add_argument("--sinkhorn-iters",type=int,default=50)
    ap.add_argument("--potential-beta",type=float,default=5.0)
    ap.add_argument("--dsm-mode",choices=['joint_legacy','partial_blood','blood_marginal'],default='partial_blood')
    ap.add_argument("--patience",type=int,default=12)
    ap.add_argument("--min-epochs",type=int,default=50)
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
        approved_scope="User-authorized full candidate rescreen with audited input/training repairs (2026-10-02)",
        command="python scripts/step1_world_model/ccwm/train_ccwm.py",
    )
    write_run_manifest(manifest, args.run_dir)

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    rng = np.random.RandomState(args.seed)
    torch.set_num_threads(2)
    dev = torch.device(args.device if torch.cuda.is_available() or args.device == "cpu" else "cpu")

    # ---- load frozen tensors (GPU-resident, reference-engine pattern) ----
    brain = np.load(Path(args.data_dir) / "brain_train.npz")
    # Locked test is evaluated only after recipe and checkpoints are frozen.
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
    cite_train_mask = torch.from_numpy(~(cite["is_test"] | cite["is_val"])).to(dev)
    cite_test_mask = torch.from_numpy(cite["is_test"]).to(dev)
    donor_cite = torch.from_numpy(cite["donor"]).to(dev)
    n_genes = X_brain.shape[1]
    n_prot = Y_cite.shape[1]

    brain_train=~brain['is_val']
    blood_train=~(blood['is_val']|blood['is_test'])
    brain_groups=np.array([str(d)+'|'+str(st) for d,st in zip(brain['donor'],brain['state'])])
    blood_groups=np.array([str(d)+'|'+str(st) for d,st in zip(blood['donor'],blood['state'])])
    cite_groups=np.array([str(d)+'|'+str(st) for d,st in zip(cite['donor'],cite['subtype'])])
    sample_brain=BalancedSampler(brain_groups,np.where(brain_train)[0])
    sample_blood=BalancedSampler(blood_groups,np.where(blood_train)[0])
    sample_cite=BalancedSampler(cite_groups,np.where(~(cite['is_test']|cite['is_val']))[0])
    # Block A strata: (cond2, blood_state, brain_state) over all purified pairs
    strata = []
    n_blood_all = X_blood.shape[0]
    blood_idx_by_stratum, brain_idx_by_stratum = {}, {}
    for cond2 in (0, 1):
        for bs in BLOOD_STATE_IDS:
            for brs in BRAIN_STATE_IDS:
                bi = np.where(blood_train & (blood["cond2"] == cond2) & (blood["state"] == bs))[0]
                Br = np.where(brain_train & (brain["cond2"] == cond2) & (brain["state"] == brs))[0]
                if sum(np.unique(blood['donor'][bi],return_counts=True)[1]>=25)>=2 and sum(np.unique(brain['donor'][Br],return_counts=True)[1]>=25)>=2:
                    blood_idx_by_stratum[(cond2, bs, brs)] = BalancedSampler(blood["donor"],bi)
                    brain_idx_by_stratum[(cond2, bs, brs)] = BalancedSampler(brain["donor"],Br)
                    strata.append((cond2, bs, brs, len(bi), len(Br)))

    if not strata: raise ValueError('No adequately sampled training strata')
    # Cost scale fitted once on donor/state-balanced training brain observations.
    si=sample_brain.draw(4096,np.random.RandomState(918))
    ot_scale=float(2*X_brain[si].var(0).sum())
    prot_scale=float(Y_cite[cite_train_mask].var(0).mean())
    cfg = CCWMConfig(n_genes=n_genes, n_proteins=n_prot, k_steps=args.k_steps, coupling=(args.coupling == "on"),
                     cond_dim=2 + N_STATES + 2, potential_beta=args.potential_beta, ot_scale=ot_scale, sample_latent_a=args.sample_latent_a, sinkhorn_iters=args.sinkhorn_iters, dsm_mode=args.dsm_mode)
    model = CCWM(cfg).to(dev)  # v2.2: no ESM2, no gene_mean (blueprint hard-constraint 2)
    model.blood_center.copy_(torch.as_tensor(blood['gene_mu'],device=dev))
    model.brain_center.copy_(torch.as_tensor(brain['gene_mu'],device=dev))
    initial_p=model.u_brain.net[0].weight[:,cfg.z_dim:cfg.z_dim+cfg.n_proteins].detach().clone()
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir/'training_inputs.json').write_text(json.dumps(dict(args=vars(args),strata=strata,
        n_genes=n_genes,n_proteins=n_prot,ot_scale=ot_scale,protein_variance=prot_scale,
        brain_train_donors=np.unique(brain['donor'][brain_train]).tolist(),
        blood_train_donors=np.unique(blood['donor'][blood_train]).tolist(),
        cite_train_donors=np.unique(cite['donor'][~(cite['is_val']|cite['is_test'])]).tolist()),indent=2))
    log_path = outdir / "train_log.csv"
    with log_path.open("w", newline="") as f:
        csv.writer(f).writerow(["step", "loss", "dsm", "protein", "recon", "kl", "sinkhorn"])

    n_cite_train = int(cite_train_mask.sum())
    cite_train_ids = torch.nonzero(cite_train_mask).squeeze(1)
    step = 0
    t0 = time.time()
    best_loss, best_state = float("inf"), None  # v2.0: ad-engine best-checkpoint tracking
    epoch_loss_acc = []
    best_epoch=-1
    validation_log=[]
    model.train()
    # v2.1 β-annealing: KL weight 0 -> beta_kl over first N epochs
    total_epochs = args.epochs
    kl_target = cfg.beta_kl
    for epoch in range(args.epochs):
        # v2.1 β-annealing schedule
        cfg.beta_kl = kl_target * min(1.0, (epoch + 1) / cfg.kl_anneal_epochs)
        for _ in range(args.steps_per_epoch):
            # Block B: CITE-seq train donors
            ib = sample_cite.draw(args.batch_b,rng)
            lb = model.block_b(X_cite[ib], Y_cite[ib], c_cite[ib])
            # Block A: one stratum per step; kill-gate-1 random pairing draws the
            # blood side uniformly from ALL blood cells (stratum-blind)
            cond2, bs, brs, _, _ = strata[rng.randint(len(strata))]
            if args.pairing == "random":
                ibl = sample_blood.draw(args.batch_a,rng)
            else:
                bi = blood_idx_by_stratum[(cond2, bs, brs)]
                ibl = bi.draw(args.batch_a,rng)
            Br = brain_idx_by_stratum[(cond2, bs, brs)]
            ibr = Br.draw(args.batch_a,rng)
            ibl_t = torch.from_numpy(np.asarray(ibl)).to(dev)
            ibr_t = torch.from_numpy(ibr).to(dev)
            la = model.block_a(X_blood[ibl_t], c_blood[ibl_t], X_brain[ibr_t], c_brain[ibr_t])
            # Reconstruction
            irb = sample_brain.draw(args.batch_r,rng)
            irl = sample_blood.draw(args.batch_r,rng)
            lr_ = model.recon(X_blood[irl], X_brain[irb])
            loss = (cfg.beta_dsm * lb["dsm"] + cfg.beta_protein * lb["protein"] + cfg.beta_recon * lb["recon"]
                    + cfg.beta_recon * lr_["recon"] + cfg.beta_kl * (lb["kl"] + lr_["kl"])
                    + cfg.beta_sinkhorn * la["sinkhorn"])  # v2.1: variance floor removed (free bits replaces)
            opt.zero_grad(set_to_none=True)
            if not torch.isfinite(loss): raise FloatingPointError((step,float(loss)))
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
        # Validate every 5 epochs. Checkpoints never use training loss or locked tests for selection.
        if (epoch+1)%5==0 or epoch==0:
            vr=validation(model,brain,blood,cite,X_brain,X_blood,X_cite,Y_cite,c_brain,BRAIN_STATE_IDS,BLOOD_STATE_IDS,prot_scale)
            vr.update(epoch=epoch+1,step=step,train_loss=el,seconds=time.time()-t0,
                      protein_weight_update=float((model.u_brain.net[0].weight[:,cfg.z_dim:cfg.z_dim+cfg.n_proteins]-initial_p).abs().max()))
            validation_log.append(vr)
            (outdir/'validation.json').write_text(json.dumps(validation_log,indent=2))
            print(json.dumps({k:v for k,v in vr.items() if k!='strata'}),flush=True)
            if vr['selection_score'] < best_loss-1e-5:
                best_loss=vr['selection_score'];best_epoch=epoch+1
                best_state={k:v.detach().cpu().clone() for k,v in model.state_dict().items()}
                torch.save({'state_dict':best_state,'cfg':cfg.__dict__,'args':vars(args),'step':step,
                            'best_epoch':best_epoch,'best_validation_score':best_loss},outdir/'model.pending.pt')
                (outdir/'model.pending.pt').replace(outdir/'model.pt')
            if epoch+1>=args.min_epochs and (epoch+1-best_epoch)>=args.patience*5:
                break
    if best_state is None: raise RuntimeError('no validated checkpoint')
    summary=dict(seed=args.seed,steps=step,epochs=epoch+1,best_epoch=best_epoch,
                 best_validation_score=best_loss,seconds=time.time()-t0,
                 checkpoint='model.pt',status='training_complete')
    (outdir/'completed.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary),flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
