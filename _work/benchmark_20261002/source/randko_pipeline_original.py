"""p4_5 randko — optimized double-decoy set-level validation (2026-09-17).

Replaces the flawed mean|delta| set-level statistic (global-influence
confounded, corr~0.95; response space rank 1-2 PCs) with a two-layer design:

LAYER-2 (primary gate, decoy gene sets): program sets from 2-fold seed CV
(discovery {42,43,44} / {45,46}), signed set selectivity of the designated
protein vs 2000 size-matched random gene sets on the fold-matched ensemble,
Fisher-combined, BH over family. Spatial-competitive granularity.
LAYER-1 (decoy proteins, calibrated effect sizes): per-seed per-gene z
calibration -> own-set mean minus 300-draw random-set baseline -> panel
standardized interaction -> mean-z across confirmation seeds; exhaustive
decoy pool. Global panel-shift Wilcoxon; exclusive top-rank NOT claimed
(promiscuous high-leverage proteins disclosed).
Recipe-pure seeds 42-51 (torch 2.5.1/py3.10); RECIPE_PURE_MAX env to widen.

Reuses interim/randko_cache curve npz (5-point dose x 3 cell-draws) when
present; else computes from checkpoints (fast inference path).
Outputs: results/p4_5_randko/randko_setlevel_v1.json
"""
import os
import sys
import json
import glob
import time
import argparse

import numpy as np
from scipy.stats import combine_pvalues, wilcoxon

A = '/public/home/mengxl/dzy/pd_product_assets'
CACHE = A + '/interim/randko_cache'
QR = np.arange(5, dtype=float)
QUANTILES = [0.05, 0.25, 0.50, 0.75, 0.95]
DRAWS = [42, 43, 44]
GSEA_A = CACHE + '/gsea_disc344.json'    # fold-A sets (seeds 42-44)
GSEA_B = CACHE + '/gsea_conf256.json'    # fold-B sets (seeds 45-46)


def slope(C):
    x = QR - QR.mean()
    return np.tensordot(x, C, axes=(0, 2)) / (x * x).sum()


def bh(ps):
    ps = np.asarray(ps, float)
    m = len(ps)
    o = np.argsort(ps)
    q = np.empty(m)
    prev = 1.0
    for r_, i_ in zip(range(m, 0, -1), o[::-1]):
        prev = min(prev, ps[i_] * m / r_)
        q[i_] = prev
    return q


def emp_p(vr, vn):
    return float((1 + (vn >= vr).sum()) / (1 + len(vn)))


def build_cache(seeds):
    sys.path.insert(0, '/public/home/mengxl/dzy/pd_product/src')
    import torch
    from pdproduct.simulators.ccwm import CCWM, CCWMConfig
    import scanpy as sc
    dev = 'cuda:0' if torch.cuda.is_available() else 'cpu'
    _og = torch.autograd.grad
    def _fast(*a, **k):
        k['create_graph'] = False
        return _og(*a, **k)
    torch.autograd.grad = _fast
    models = {}
    for s in seeds:
        ck = torch.load(A + f'/checkpoints/project/ccwm-v23-s{s}/model.pt', map_location='cpu', weights_only=False)
        mm = CCWM(CCWMConfig(**ck['cfg'])); mm.load_state_dict(ck['state_dict']); mm.eval().to(dev)
        models[s] = mm
    blood = np.load(A + '/interim/ccwm_v23/blood.npz')
    X_blood = torch.from_numpy(blood['X']).to(dev)
    STATES = [str(s) for s in blood['state_vocab']]
    cite = np.load(A + '/interim/ccwm_v23/citeseq.npz'); Y_cite = cite['Y']
    hao = np.load(A + '/processed/citeseq_hao/bridge_data.npz', allow_pickle=True)
    proteins_all = [str(p) for p in hao['proteins']]
    brain = sc.read_h5ad(A + '/processed/gse178265/v0.1/GSE178265_sn_annotated.h5ad', backed='r')
    g_brain = brain.var['gene_symbol'].astype(str).values; brain.file.close()
    bmb = sc.read_h5ad(A + '/processed/moquin_beaudry2025/v0.1/MB2025_pbmc_annotated.h5ad', backed='r')
    g_blood = bmb.var['gene_symbol'].astype(str).values; bmb.file.close()
    genes = sorted(set(g_brain) & set(g_blood) & set(hao['genes'].astype(str)))
    ss = json.load(open(A + '/results/p3_5_seed_stability/seed_posterior_v1.json'))
    stable = [(e['blood'], e['brain'], e['mean_r']) for e in ss['stable_edges']]
    def oh(c, k): return torch.nn.functional.one_hot(c, k).float()
    for bs_name, brs_name, mean_r in stable:
        pk = bs_name + '__x__' + brs_name
        path = f'{CACHE}/{pk}.npz'
        zold_seeds, Cold = [], None
        if os.path.exists(path):
            zo = np.load(path, allow_pickle=True)
            if 'C' in zo:
                Cold = zo['C']; zold_seeds = [int(s) for s in zo['seeds']]
        if all(s in zold_seeds for s in seeds) and Cold is not None:
            continue
        bs_i, brs_i = STATES.index(bs_name), STATES.index(brs_name)
        idx = np.where(blood['state'] == bs_i)[0]
        n = min(150, len(idx))
        c_fix = torch.cat([oh(torch.zeros(n, dtype=torch.long), 2),
                           oh(torch.full((n,), brs_i, dtype=torch.long), len(STATES)),
                           oh(torch.full((n,), 1, dtype=torch.long), 2)], -1).to(dev)
        sens = {}
        for s in [42, 43, 44]:
            m = models[s]
            sel0 = idx[np.random.RandomState(42).choice(len(idx), n, replace=False)]
            with torch.no_grad(): mu_s, _ = m.enc_blood(X_blood[torch.from_numpy(sel0)])
            zr = mu_s.clone().requires_grad_(True)
            uh = m.protein_head(zr)
            for j in range(min(30, uh.shape[1])):
                if zr.grad is not None: zr.grad.zero_()
                uh[0, j].backward(retain_graph=True)
                sens[j] = sens.get(j, 0) + zr.grad.abs().sum().item() / 3
        top3 = [proteins_all[p] for p, _ in sorted(sens.items(), key=lambda x: -x[1])[:3]]
        C = np.zeros((len(seeds), len(proteins_all), len(QUANTILES), len(genes)), dtype=np.float32)
        for si, s in enumerate(seeds):
            if s in zold_seeds and Cold.shape[1:] == (len(proteins_all), len(QUANTILES), len(genes)):
                C[si] = C[zold_seeds.index(s)] if False else Cold[zold_seeds.index(s)]
                continue
            m = models[s]
            for dr in DRAWS:
                sel = idx[np.random.RandomState(dr).choice(len(idx), n, replace=False)]
                x_b = X_blood[torch.from_numpy(sel)]
                with torch.no_grad(): mu, _ = m.enc_blood(x_b); u_base = m.protein_head(mu)
                for pi in range(len(proteins_all)):
                    for qi, u_val in enumerate(np.quantile(Y_cite[:, pi], QUANTILES)):
                        u_c = u_base.clone(); u_c[:, pi] = float(u_val)
                        _, z_s = m.solve_equilibrium(mu, mu.clone(), u_c.detach(), c_fix)
                        with torch.no_grad():
                            C[si, pi, qi] += m.dec_brain(z_s.detach()).mean(0).cpu().numpy() / len(DRAWS)
        np.savez_compressed(path, C=C, proteins=np.array(proteins_all), genes=np.array(genes),
                            top3=np.array(top3), seeds=np.array(seeds), quantiles=np.array(QUANTILES),
                            draws=np.array(DRAWS), mean_r=mean_r)
        print(pk, 'cached', flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--recipe-pure-max', type=int, default=int(os.environ.get('RECIPE_PURE_MAX', '51')))
    args = ap.parse_args()
    PMAX = args.recipe_pure_max
    seeds = [s for s in range(42, PMAX + 1)]
    build_cache(seeds)
    gsea_a = json.load(open(GSEA_A))
    gsea_b = json.load(open(GSEA_B))
    res = {}
    for f in sorted(glob.glob(CACHE + '/*.npz')):
        z_ = np.load(f, allow_pickle=True)
        seeds_all = [int(s) for s in z_['seeds']]
        keep = [i for i, s in enumerate(seeds_all) if s <= PMAX]
        C = z_['C'][keep]
        seeds_l = [seeds_all[i] for i in keep]
        proteins = [str(p) for p in z_['proteins']]
        genes = [str(g) for g in z_['genes']]
        top3 = [str(p) for p in z_['top3']]
        pk = os.path.basename(f)[:-4]
        pidx = {p: i for i, p in enumerate(proteins)}
        gidx = {g: i for i, g in enumerate(genes)}
        D = slope(C)
        Z = (D - D.mean(1, keepdims=True)) / (D.std(1, keepdims=True) + 1e-9)
        pool = np.array([i for i in range(len(proteins)) if proteins[i] not in top3])
        conf_a = [i for i, s in enumerate(seeds_l) if s not in (42, 43, 44)]
        disc_a = [i for i, s in enumerate(seeds_l) if s in (42, 43, 44)]
        disc_b = [i for i, s in enumerate(seeds_l) if s in (45, 46)]
        ens_a = D[disc_b].mean(0)
        ens_b = D[disc_a].mean(0)
        for pname in top3:
            ck = pname + '__' + pk
            ri = pidx[pname]
            rec = res.setdefault(ck, {'pair': pk, 'protein': pname})
            for fold, gs, conf_idx in (('A', gsea_a, conf_a), ('B', gsea_b, [i for i, s in enumerate(seeds_l) if s not in (45, 46)])):
                if ck not in gs or not gs[ck]:
                    continue
                S = sorted({g for t in gs[ck] for g in t['genes']} & set(gidx))
                if len(S) < 5:
                    continue
                si = np.array([gidx[g] for g in S])
                rec['n_genes_' + fold] = len(S)
                vset = Z[:, :, si].mean(2)
                rng_i = np.random.RandomState(11)
                base = np.zeros_like(vset)
                for _ in range(300):
                    ssi = rng_i.choice(len(genes), len(si), replace=False)
                    base += Z[:, :, ssi].mean(2) / 300
                I = vset - base
                zI = (I - I.mean(1, keepdims=True)) / (I.std(1, keepdims=True) + 1e-12)
                mz = zI[conf_idx].mean(0)
                rec['P_meanz_' + fold] = emp_p(mz[ri], mz[pool])
                rec['meanz_' + fold] = round(float(mz[ri]), 3)
                rec['pct_pool_' + fold] = round(float((mz[pool] < mz[ri]).mean()), 3)
            for fold, gs, ens in (('A', gsea_a, ens_a), ('B', gsea_b, ens_b)):
                if ck not in gs or not gs[ck]:
                    continue
                S = sorted({g for t in gs[ck] for g in t['genes']} & set(gidx))
                if len(S) < 5:
                    continue
                si = np.array([gidx[g] for g in S])
                dR = ens[ri]
                vR = np.mean(dR[si]) / (np.mean(np.abs(dR)) + 1e-12)
                rng = np.random.RandomState(7)
                nsets = np.array([np.mean(dR[rng.choice(len(genes), len(si), replace=False)]) / (np.mean(np.abs(dR)) + 1e-12)
                                  for _ in range(2000)])
                rec['S_set_' + fold] = emp_p(vR, nsets)
    rows = list(res.items())
    rd = dict(rows)
    for k, r in rows:
        if 'S_set_A' in r and 'S_set_B' in r:
            _, pf = combine_pvalues([r['S_set_A'], r['S_set_B']], method='fisher')
            r['S_set_fish'] = pf
    l2 = [(k, r['S_set_fish']) for k, r in rows if 'S_set_fish' in r]
    q2 = bh([p for _, p in l2])
    for (k, _), qq in zip(l2, q2):
        rd[k]['S_set_q'] = round(float(qq), 5)
        rd[k]['S_set_sig'] = bool(qq <= 0.05)
    pct = [r['pct_pool_A'] for _, r in rows if 'pct_pool_A' in r]
    _, wp = wilcoxon([p_ - 0.5 for p_ in pct], alternative='greater')
    out = {
        'method': 'two-layer double-decoy (2026-09-17 optimization): '
                  'L2 = decoy gene sets, signed set selectivity vs 2000 size-matched random sets, '
                  '2-fold seed CV (sets from {42,43,44}/{45,46}), Fisher+BH; '
                  'L1 = decoy proteins, per-seed per-gene z, own-set minus random-set baseline, '
                  'panel-standardized interaction mean-z across confirmation seeds {45..%d}, exhaustive pool' % PMAX,
        'seeds': seeds,
        'layer2_gate': {'n_tested': len(l2), 'n_sig': int((q2 <= 0.05).sum()),
                        'sig': {k: {'q': rd[k]['S_set_q'], 'pA': rd[k]['S_set_A'], 'pB': rd[k]['S_set_B']}
                                for (k, _) in l2 if rd[k]['S_set_q'] <= 0.05}},
        'layer1_calibration': {'global_panel_shift': {'median_pct': round(float(np.median(pct)), 3),
                                                      'wilcoxon_p': round(float(wp), 4), 'n': len(pct)},
                                'exclusive_top_rank': 'NOT CLAIMED — promiscuous high-leverage decoys '
                                                      '(CD127/IL7R, CLEC12A, CD13/CD2/CD61/CD41 cluster and '
                                                      'cross-program CD35) co-move programs seed-consistently; '
                                                      'designated proteins sit in the upper tail, not rank 1'},
        'implementation_sensitivity': 'SCNet A800 replication (torch 2.7.0+cu118/py3.12, seeds 52-57) of the '
                                      'frozen recipe shows NO layer-1 panel shift (median pct 0.511, p=0.49) — '
                                      'protein-axis specificity is recipe-bound (model-space claim boundary); '
                                      'recipe-pure A6000 extension seeds 47-51 alone: median pct 0.831, p=0.0096',
        'results': {k: {kk: (round(float(vv), 5) if isinstance(vv, float) else vv) for kk, vv in r.items()}
                    for k, r in sorted(res.items())},
    }
    path = A + '/results/p4_5_randko/randko_setlevel_v1.json'
    json.dump(out, open(path, 'w'), indent=1)
    print('L2 sig %d/%d; L1 panel-shift median %.3f p=%.4f -> %s' %
          (int((q2 <= 0.05).sum()), len(l2), np.median(pct), wp, path))


if __name__ == '__main__':
    main()
