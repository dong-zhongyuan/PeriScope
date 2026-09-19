"""p4 dose-response + gene-level randko + GSEA on the curve cache (v2, 2026-09-17).

Upgrade: seed ensemble 42-51 (recipe-pure, frozen torch 2.5.1/py3.10) x 3
cell-draws x 5 quantile points, read from interim/randko_cache curve npz
(built by scripts/p4_5_randko/pipeline.py build_cache; values bitwise identical
to inline solves). Gene-level randko now uses the EXHAUSTIVE decoy pool
(228 minus the pair's top3) instead of the std-matched n300 sample.
GSEA family definition: top-200 up genes of the 5-point SLOPE ensemble delta (readout upgrade 2026-09-18: slope proven lower-noise than hi/lo in the randko CV work), enrichr GO_BP+KEGG, q<0.05, top6 terms.
Outputs (overwrite in place, no version variants):
  p4_dose_response/dose_response_stable_v1.json
  p4_5_randko/randko_stable_v1.json (gene level)
  p5_gsea/gsea_stable_v1.json
Env: CCWM_SEEDS (default 42,43,44,45,46,47,48,49,50,51), RECIPE_PURE_MAX=51.
"""
import os
import sys
import json

import numpy as np

A = '/public/home/mengxl/dzy/pd_product_assets'
CACHE = A + '/interim/randko_cache'
sys.path.insert(0, '/public/home/mengxl/dzy/pd_product/src')
PMAX = int(os.environ.get('RECIPE_PURE_MAX', '51'))
SEEDS = [int(x) for x in os.environ.get('CCWM_SEEDS', '42,43,44,45,46,47,48,49,50,51').split(',')]
SEEDS = [s for s in SEEDS if s <= PMAX]
QUANTILES = [0.05, 0.25, 0.50, 0.75, 0.95]

import gseapy as gp

ss = json.load(open(A + '/results/p3_5_seed_stability/seed_posterior_v1.json'))
stable_pairs = [(e['blood'], e['brain'], e['mean_r']) for e in ss['stable_edges']]
print('Seed-stable pairs:', len(stable_pairs))
print('Ensemble seeds (recipe-pure):', SEEDS)

all_dr = {}
all_randko = {}
all_gsea = {}

for bs_name, brs_name, mean_r in stable_pairs:
    pk = bs_name + '__x__' + brs_name
    path = f'{CACHE}/{pk}.npz'
    if not os.path.exists(path):
        raise SystemExit(f'curve cache missing for {pk}; run scripts/p4_5_randko/pipeline.py first (build_cache)')
    z = np.load(path, allow_pickle=True)
    seeds_avail = [int(s) for s in z['seeds']]
    if not set(SEEDS) <= set(seeds_avail):
        raise SystemExit(f'{pk}: cache lacks seeds {sorted(set(SEEDS) - set(seeds_avail))}')
    C = z['C'][[seeds_avail.index(s) for s in SEEDS]]      # (nseed, prot, q, gene)
    proteins = [str(p) for p in z['proteins']]
    genes = [str(g) for g in z['genes']]
    top3_names = [str(p) for p in z['top3']]
    pidx = {p: i for i, p in enumerate(proteins)}
    gidx = {g: i for i, g in enumerate(genes)}

    ens = C.mean(0)                                        # (prot, q, gene) ensemble mean
    QRx = np.arange(len(QUANTILES), dtype=float); XS = QRx - QRx.mean()
    D = np.tensordot(XS, ens, axes=(0, 1)) / (XS * XS).sum()  # slope delta (5-point), 10-seed x 3-draw
    pool = [i for i in range(len(proteins)) if proteins[i] not in top3_names]

    for pname in top3_names:
        ck = pname + '__' + pk
        ri = pidx[pname]
        rd = D[ri]
        # gene-level randko: top-10 |delta| genes vs exhaustive decoy pool
        top_gi = np.argsort(-np.abs(rd))[:10]
        rec = {}
        for gi in top_gi:
            nv = np.abs(D[pool, gi])
            q = float((1 + (nv >= abs(rd[gi])).sum()) / (1 + len(nv)))
            rec[genes[gi]] = {'delta': round(float(rd[gi]), 4), 'q': round(q, 4), 'sig': bool(q <= 0.05)}
        all_randko[ck] = rec
        # GSEA family: top-200 up genes of ensemble hi/lo delta
        ranked = sorted(zip(genes, rd), key=lambda x: -x[1])
        try:
            enr = gp.enrichr(gene_list=[g for g, _ in ranked[:200]],
                             gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
                             organism='human', outdir=None, no_plot=True, cutoff=0.05)
            sig_df = enr.results[enr.results['Adjusted P-value'] < 0.05]
            all_gsea[ck] = [
                {'term': row['Term'], 'p_adj': float(row['Adjusted P-value']),
                 'genes': row['Genes'].split(';')[:6]}
                for _, row in sig_df.head(6).iterrows()]
        except Exception:
            all_gsea[ck] = []
        nr_sig = sum(1 for d in rec.values() if d['sig'])
        print(f'  {pk} | {pname}: randko {nr_sig}/10 sig, gsea {len(all_gsea[ck])} terms', flush=True)

    # dose-response curves for the pair's top proteins (ensemble mean)
    all_dr[pk] = {
        'blood_state': bs_name, 'brain_state': brs_name, 'seed_mean_r': mean_r,
        'proteins': top3_names,
        'curves': {}, 'top_up': [], 'top_down': []}
    rd0 = D[pidx[top3_names[0]]]
    all_dr[pk]['top_up'] = [(genes[i], round(float(rd0[i]), 4)) for i in np.argsort(-rd0)[:5]]
    all_dr[pk]['top_down'] = [(genes[i], round(float(rd0[i]), 4)) for i in np.argsort(rd0)[:5]]
    Y_cite = np.load(A + '/interim/ccwm_v23/citeseq.npz')['Y']
    hao = np.load(A + '/processed/citeseq_hao/bridge_data.npz', allow_pickle=True)
    prot_all = [str(p) for p in hao['proteins']]
    for pname in top3_names:
        pi = prot_all.index(pname)
        all_dr[pk]['curves'][pname] = {
            'u': np.quantile(Y_cite[:, pi], QUANTILES).tolist(),
            'x_mean': ens[pidx[pname], :, :].tolist()}

json.dump(all_dr, open(A + '/results/p4_dose_response/dose_response_stable_v1.json', 'w'), indent=1)
json.dump(all_randko, open(A + '/results/p4_5_randko/randko_stable_v1.json', 'w'), indent=1)
json.dump(all_gsea, open(A + '/results/p5_gsea/gsea_stable_v1.json', 'w'), indent=1)

print(f'\n=== Summary (seeds {SEEDS[0]}-{SEEDS[-1]}, {len(SEEDS)} seeds, 3 draws) ===')
print(f'Stable pairs: {len(all_dr)}')
print(f'Gene-level randko sig (q<=0.05): {sum(1 for v in all_randko.values() for d in v.values() if d["sig"])}'
      f'/{sum(len(v) for v in all_randko.values())}')
print(f'GSEA sig combos: {sum(1 for v in all_gsea.values() if v)}')
