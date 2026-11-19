"""Benchmark on the fixed targets: locked 3v2 and external cohorts. All rho recomputed from per-seed predictions."""
from bm_common import *
import time
t0 = time.time(); NPERM = 5000
EV = SRC / 'eval'
Q = pickle.load(open(RES / 'verify-baselines-B4-a/preds.pkl', 'rb')); AXP, BP = Q['AXP'], Q['P']
REG = []
def reg(method, group, training, tier, readouts, n_seeds, note='', earlier='', source=''):
    REG.append(dict(method=method, group=group, training_donors=training, tier=tier, tier_earlier=earlier, readouts=readouts, n_seeds=n_seeds, note=note, source=source))
def axis_model(method, group, training, tier, D, seeds, note='', earlier='', source=''):
    """D[(seed, axis)] -> predicted PD-minus-control vector"""
    R = {}
    for rd in ('single_seed_mean', 'ensemble_z', 'ensemble_rawdelta_per_axis'):
        R[rd] = {}
    for st in ST:
        ax = [a for a in AX if a.endswith('__' + st)]
        U = np.stack([D[(s, a)] for s in seeds for a in ax]).astype(float)
        R['single_seed_mean'][st] = U
        R['ensemble_z'][st] = np.mean([Z(u) for u in U], 0)[None]
        R['ensemble_rawdelta_per_axis'][st] = np.stack([np.mean([D[(s, a)] for s in seeds], 0) for a in ax]).astype(float)
    reg(method, group, training, tier, R, len(seeds), note, earlier, source)
def brain_model(method, group, training, tier, D, seeds, note='', earlier='', source=''):
    R = {k: {} for k in ('single_seed_mean', 'ensemble_z', 'ensemble_rawdelta_per_axis')}
    for st in ST:
        U = np.stack([D[(s, st)] for s in seeds]).astype(float)
        R['single_seed_mean'][st] = U; R['ensemble_z'][st] = np.mean([Z(u) for u in U], 0)[None]; R['ensemble_rawdelta_per_axis'][st] = U.mean(0)[None]
    reg(method, group, training, tier, R, len(seeds), note, earlier, source)
S3 = (42, 43, 44); S10 = tuple(range(42, 52)); S7 = (42, 43, 44, 45, 46, 50, 51); TR9 = '9 fit donors (5 PD, 4 control)'; TR13 = '13 non-test donors (7 PD, 6 control)'
PS = {(s, a): AXP[('periscope', s, a)] for s in S10 for a in AX}
axis_model('PeriScope, seeds 42-51', 'PeriScope', TR9, 'A*', PS, S10, 'published model, all ten stored seeds. Seeds 47-49 were retrained on 2026-10-03 under the post-warmup rule after their first checkpoints (epochs 20/15/15) had been scored on the locked test, so the ten-seed rows are A*; epoch choice is validation only', 'A', source='benchmark_20261002/results/seed_*/periscope__*.npz via preds.pkl')
axis_model('PeriScope, 7 seeds never repaired (42-46, 50, 51)', 'PeriScope', TR9, 'A', PS, S7, 'the ten published seeds without the three repaired ones; added in the corrected version as the tier A level', source='same')
axis_model('PeriScope, seeds 42-44', 'PeriScope', TR9, 'A', PS, S3, 'benchmark seeds; the three seeds run high relative to the ten', source='same')
NOTE_SHARED = 'A*: the shared-representation comparison was revised after the initial results were inspected (benchmark_20261002/benchmark_scope_index.json post_result_revision=true; 评估与运行说明.txt; pd_product/_work/benchmark_shared_20261003/protocol.txt). Earlier comparator outputs were replaced without an archive (superseded_files.json, cleanup_record.json), so the size of the change cannot be measured'
NOTE_CELLOT = 'A*: CellOT was in the non-shared set removed by cleanup_obsolete.py and its current result files were rebuilt on 2026-10-04 on its own 64-component PCA projection; cellot_restored/execution_modes.json shows three numerical execution modes (seed 42 cDC__astro locked rho 0.1609, 0.2145, 0.1910) and the mode was picked by matching archived values, so small orderings against CellOT are not stable. CellOT pDC axes diverge (see REPORT 4.4)'
for m in ['conditional_mlp', 'conditional_icnn', 'cellot', 'shared_label_marginal', 'shared_label_mean', 'ot_ridge']:
    axis_model(m, 'published comparator', TR9, 'A*', {(s, a): AXP[(m, s, a)] for s in S3 for a in AX}, S3, NOTE_CELLOT if m == 'cellot' else NOTE_SHARED, 'A', source='benchmark_20261002/results via preds.pkl')
brain_model('CPA', 'added model', TR9, 'A', {(s, st): BP[('cpa', s)][st] for s in S3 for st in ST}, S3, 'trained on 9 fit donors, stopped on val pair', source='research_20261118/baselines/runs')
brain_model('scGen', 'added model', TR9, 'A', {(s, st): BP[('scgen', s)][st] for s in S3 for st in ST}, S3, 'trained on 9 fit donors, stopped on val pair', source='research_20261118/baselines/runs')
brain_model('scGen (decode_means readout)', 'added model', TR9, 'A', {(s, st): BP[('scgen_decode_means', s)][st] for s in S3 for st in ST}, S3, source='research_20261118/baselines/runs')
for tag in ['gene_cmlp', 'gene_cmlp_contrast', 'gene_cmlp_full_lam10', 'gene_cmlp_full_nolam', 'pca_cmlp', 'cfm']:
    Zs = {s: np.load(RES / f'baselines/runs/{tag}_s{s}.npz') for s in S3}
    axis_model(tag, 'added model', TR9, 'A', {(s, a): Zs[s][f'{a}|ca'] for s in S3 for a in AX}, S3, 'gene-space or PCA baseline from baselines lane', source='research_20261118/baselines/runs')
    if f'{AX[0]}|ca_last' in Zs[42].files:
        axis_model(tag + ' (last step)', 'added model', TR9, 'A', {(s, a): Zs[s][f'{a}|ca_last'] for s in S3 for a in AX}, S3, source='research_20261118/baselines/runs')
# ablation arms under the corrected checkpoint rule, and the two published ineligible variants for reference
def arm(run_seed): return {(s, a): v for (r, s) in run_seed for a, v in np.load(EV / 'deltas' / f'{r}__s{s}.npz').items()}
ARMD = {}
NE = 'not eligible (before warmup end)'
def arms(name, rs, tier, note='', eligible=True):
    D = arm(rs); ARMD[name] = (D, tuple(s for _, s in rs)); axis_model(name, 'ablation arm (post-warmup rule)' if eligible else 'ablation arm (published checkpoint, NOT eligible)', TR9, tier, D, tuple(s for _, s in rs), note, '' if eligible else 'A', source='benchmark_20261119/src/eval/deltas')
arms('random_pairing', [('random_pairing__retrained', s) for s in S3], 'A*', 'all seeds retrained after the published test collapse was seen; epoch choice by validation only')
arms('coupling_off', [('coupling_off__published_ckpt', 42), ('coupling_off__published_ckpt', 43), ('coupling_off__retrained', 44)], 'A*', 'mean includes seed 44, retrained after its published test collapse was seen; epoch choice by validation only')
arms('no_blood', [('no_blood', 42), ('no_blood', 43)], 'A', 'two seeds; seed 44 dropped from queue')
arms('no_cond2', [('no_cond2', s) for s in S3], 'A'); arms('no_cond2_random', [('no_cond2_random', s) for s in S3], 'A'); arms('label_shuffle', [('label_shuffle', 42)], 'A', 'one seed')
arms('random_pairing (published checkpoints, epoch 10-15)', [('random_pairing__published_ckpt', s) for s in S3], NE, 'published checkpoint before KL warmup end; kept for comparison only', eligible=False)
arms('coupling_off (published checkpoints, seed 44 at epoch 15)', [('coupling_off__published_ckpt', s) for s in S3], NE, 'published checkpoint; seed 44 before KL warmup end; kept for comparison only', eligible=False)
# PeriScope with added donors or contrast loss (train_upper_bound runs; fixed epoch, not validation selected)
def tub(run, seeds, ep):
    D = {}
    for s in seeds:
        z = np.load(RES / f'train_upper_bound/runs/{run}_s{s}/deltas/ep{ep:04d}.npz')['ca']; D.update({(s, a): z[i] for i, a in enumerate(AX)})
    return D
for name, run, seeds, ep, tr, tier, earlier, note in [
    ('PeriScope fine-tune: 13 donors + contrast(modt)', 'F6_addDonors_modt', S10, 25, TR13, 'A*', 'A (analyst) / B (verifier)', 'fine-tune 25 epochs from published checkpoint; donor set chosen after test results were seen'),
    ('PeriScope fine-tune: 13 donors, no contrast loss', 'F5_addDonors', S10, 25, TR13, 'A*', 'A', 'fine-tune 25 epochs from published checkpoint; donor set chosen after test results were seen'),
    ('PeriScope fine-tune: 9 donors + contrast(modt)', 'F2_modt', S3, 25, TR9, 'A*', 'A', 'fine-tune 25 epochs from the published checkpoints (train_upper_bound/runs/F2_modt_s*/args.json init_ckpt), adding a contrast loss on 2026-11-18, after those checkpoints already had locked scores; the contrast target was chosen by training-donor CV'),
    ('PeriScope fine-tune: 9 donors + contrast(same_all3z)', 'F3_sameall3z', (42,), 25, TR9, 'B', 'B', 'contrast target picked after seeing the test'),
    ('PeriScope fine-tune: 13 donors + contrast(same_all3z)', 'F7_allnt_sameall3z', (42,), 25, TR13, 'B', 'B', 'contrast target picked after seeing the test'),
    ('PeriScope fine-tune: 13 donors + contrast(all3z)', 'F8_allnt_all3z', (42,), 25, TR13, 'B', 'B', 'contrast target picked after seeing the test'),
    ('PeriScope from scratch: 13 donors + contrast(modt)', 'S2_allnt_modt', (42,), 100, TR13, 'A*', 'A', '100 epochs from scratch, one seed')]:
    axis_model(name, 'PeriScope variant', tr, tier, tub(run, seeds, ep), seeds, note, earlier, 'research_20261118/train_upper_bound/runs')
# 13-donor CPA / scGen finals
for m, n, lab, earlier in [('cpa', 44, 'A', 'A'), ('cpa', 58, 'A*', 'A*'), ('scgen', 60, 'A', 'A'), ('scgen', 59, 'A*', 'A*')]:
    zs = {s: np.load(RES / f'donor_cv/runs_final/{m}_final13ep{n}_s{s}.npz') for s in S3}
    brain_model(f'{"CPA" if m == "cpa" else "scGen"} 13 donors, fixed epoch {n}', 'added model', TR13, 'A*', {(s, st): zs[s][f'{st}|ca'] for s in S3 for st in ST}, S3,
                'epoch from inner splits of the 13 donors' if earlier == 'A' else 'epoch from CV medians; the CV folds include the test donors', earlier, 'research_20261118/donor_cv/runs_final')
print('registered methods', len(REG), int(time.time() - t0), 's', flush=True)

# ------------------------------------------------------------------- targets
TG = {}
# the locked target is the harness target (float32 donor profiles in the benchmark cache); pb18 was summed in float64 and differs by up to 1.7e-5 per gene, which moves rho in the 4th decimal
PROF = {}
for s in ST:
    b = np.load(BENCH / f'cache/cDC__{s}.npz', allow_pickle=True); X = np.vstack([b['profiles1'], b['profiles0']]).astype(float); assert len(b['profiles1']) == 3 and len(b['profiles0']) == 2
    for d in TEST:
        i = np.abs(X - MU[KEY[(d, s)]]).max(1).argmin(); assert np.abs(X[i] - MU[KEY[(d, s)]]).max() < 5e-4, (d, s); PROF[(d, s)] = X[i]
def AT(ds, s): return np.stack([PROF[(d, s)] for d in ds])
TG['locked_3v2'] = dict(block='locked_3v2', Y={s: AT(TP, s).mean(0) - AT(TC, s).mean(0) for s in ST}, mask={s: np.ones(NG, bool) for s in ST}, desc={s: '3 PD v 2 control' for s in ST})
assert max(np.abs(TG['locked_3v2']['Y'][s] - Q['OBS'][s]).max() for s in ST) < 1e-6
pd.DataFrame([dict(state=s, max_abs_gene_diff_pb18_vs_harness_target=float(np.abs((A(TP, s).mean(0) - A(TC, s).mean(0)) - TG['locked_3v2']['Y'][s]).max()), spearman_between_the_two_targets=float(rz(A(TP, s).mean(0) - A(TC, s).mean(0)) @ rz(TG['locked_3v2']['Y'][s]))) for s in ST]).to_csv(SRC / 'target_pb18_vs_harness.csv', index=False)
Ti = pd.read_csv(RES / 'baselines/external_pb_index.csv'); E = np.load(RES / 'baselines/external_pb.npz'); EXT = {}
for ds in ['smajic_midbrain', 'gse184950']:
    Yy, Mk, De, Don = {}, {}, {}, {}
    for s, ex in (('astro', 'astro'), ('microglia_mhc2', 'microglia')):
        q = Ti[(Ti.dataset == ds) & (Ti.state == ex) & (Ti.n >= 25) & (Ti.cond >= 0)]; X = np.stack([E[f'{ds}|{d}|{ex}'] for d in q.donor]).astype(float); y = q.cond.values.astype(int); Mk[s] = X.sum(0) > 0; Yy[s] = X[y == 1].mean(0) - X[y == 0].mean(0); De[s] = f'{(y == 1).sum()} PD v {(y == 0).sum()} control, {int(Mk[s].sum())} genes (astro-to-astro, microglia_mhc2-to-all-microglia)'; Don[s] = (q.donor.tolist(), X, y)
    TG[ds] = dict(block='external', Y=Yy, mask=Mk, desc=De); EXT[ds] = Don
I = pd.read_csv(RES / 'verify-baselines-B2-a/ma_pb_index.csv'); Pma = np.fromfile(RES / 'verify-baselines-B2-a/ma_pb.bin').reshape(len(I), NG); I['y'] = (I.cond == 'PD').astype(int)
for tn, mp in (('ma_substantia_nigra', {'astro': 'astro', 'microglia_mhc2': 'microglia_mhc2'}), ('ma_sensitivity_all_myeloid', {'astro': 'astro', 'microglia_mhc2': 'microglia'})):
    Yy, Mk, De, Don = {}, {}, {}, {}
    for s, ex in mp.items():
        q = I[(I.state == ex) & (I.n >= 25)]; X = Pma[q.index.values]; y = q.y.values; Mk[s] = np.ones(NG, bool); Yy[s] = X[y == 1].mean(0) - X[y == 0].mean(0); De[s] = f'{(y == 1).sum()} PD v {(y == 0).sum()} control, 4741 panel genes (4630 found in Ma)'; Don[s] = (q.donor.tolist(), X, y)
    TG[tn] = dict(block='external', Y=Yy, mask=Mk, desc=De); EXT[tn] = Don
Yy, Mk, De, Don = {}, {}, {}, {}   # Ma with the same gene rule as Smajic and GSE184950 (genes that are zero in all donors used are masked)
for s in ST:
    q = I[(I.state == s) & (I.n >= 25)]; X = Pma[q.index.values]; y = q.y.values; Mk[s] = X.sum(0) > 0; Yy[s] = X[y == 1].mean(0) - X[y == 0].mean(0); De[s] = f'{(y == 1).sum()} PD v {(y == 0).sum()} control, {int(Mk[s].sum())} genes (zero-sum genes masked, as for Smajic and GSE184950)'; Don[s] = (q.donor.tolist(), X, y)
TG['ma_sensitivity_zero_genes_masked'] = dict(block='external', Y=Yy, mask=Mk, desc=De); EXT['ma_sensitivity_zero_genes_masked'] = Don

# ------------------------------------------------------------------- scoring
rows_idx = [(i, rd) for i, m in enumerate(REG) for rd in m['readouts']]
def stack(s, mask):
    V, id_ = [], []
    for j, (i, rd) in enumerate(rows_idx):
        v = REG[i]['readouts'][rd][s]; V.append(v[:, mask]); id_ += [j] * len(v)
    return rz(np.vstack(V)), np.array(id_)
def scores(Rm, id_, y):
    return np.bincount(id_, weights=Rm @ y, minlength=len(rows_idx)) / np.bincount(id_, minlength=len(rows_idx))
CACHE = {}
def RR(tn, s):
    if (tn, s) not in CACHE:
        CACHE[(tn, s)] = stack(s, TG[tn]['mask'][s])
    return CACHE[(tn, s)]
RES_ = {}
for tn, T in TG.items():
    a = {s: scores(*RR(tn, s), rz(T['Y'][s][T['mask'][s]])) for s in ST}
    RES_[tn] = dict(astro=a['astro'], microglia_mhc2=a['microglia_mhc2'], mean=(a['astro'] + a['microglia_mhc2']) / 2)
print('observed scores done', int(time.time() - t0), flush=True)
# permutation of disease labels among the target donors: locked = all 10 relabellings (p floor 0.1); external = random, same labels for both states
PV = {}
rng = np.random.default_rng(20261119)
ob = RES_['locked_3v2']['mean']; n = 0; ge = np.zeros(len(ob))
for trip in itertools.combinations(TEST, 3):
    o = [d for d in TEST if d not in trip]; v = np.zeros(len(ob))
    for s in ST:
        v += scores(*RR('locked_3v2', s), rz(AT(list(trip), s).mean(0) - AT(o, s).mean(0))) / 2
    n += 1; ge += v >= ob - 1e-12
PV['locked_3v2'] = (ge / n, n)
for tn in ['smajic_midbrain', 'gse184950', 'ma_substantia_nigra', 'ma_sensitivity_all_myeloid', 'ma_sensitivity_zero_genes_masked']:
    D = EXT[tn]; don = sorted(set(D['astro'][0]) | set(D['microglia_mhc2'][0])); lab = {}
    for s in ST:
        lab.update(dict(zip(D[s][0], D[s][2])))
    y0 = np.array([lab[d] for d in don]); ob = RES_[tn]['mean']; ge = np.zeros(len(ob)); N = 0
    for b in range(NPERM):
        pl = dict(zip(don, rng.permutation(y0))); v = np.zeros(len(ob)); ok = True
        for s in ST:
            nm, X, _ = D[s]; yy = np.array([pl[d] for d in nm])
            if yy.sum() == 0 or yy.sum() == len(yy):
                ok = False; break
            v += scores(*RR(tn, s), rz((X[yy == 1].mean(0) - X[yy == 0].mean(0))[TG[tn]['mask'][s]])) / 2
        if ok:
            N += 1; ge += v >= ob - 1e-12
    PV[tn] = ((1 + ge) / (N + 1), N); print('perm', tn, N, int(time.time() - t0), flush=True)
# ------------------------------------------------------------------- long table
out = []
for tn, T in TG.items():
    r = RES_[tn]; pm, np_ = PV[tn]
    for j, (i, rd) in enumerate(rows_idx):
        m = REG[i]
        out.append(dict(block='locked_3v2 (fixed test, 5 donors)' if tn == 'locked_3v2' else f'external cohort', target=tn, method=m['method'], method_group=m['group'], added_2026_11_18_or_later='no' if m['group'] in ('PeriScope', 'published comparator', 'ablation arm (published checkpoint, NOT eligible)') else 'yes', readout=rd, training_donors=m['training_donors'], n_seeds=m['n_seeds'], tier=m['tier'], tier_earlier=m['tier_earlier'],
                        rho_astro=r['astro'][j], rho_microglia_mhc2=r['microglia_mhc2'][j], rho_mean_of_states=r['mean'][j],
                        perm_p_one_sided_rho=pm[j], n_label_permutations=np_, target_design_astro=T['desc']['astro'], target_design_microglia=T['desc']['microglia_mhc2'], note=m['note'], source=m['source']))
L = pd.DataFrame(out); L.to_csv(SRC / 'bm_fixed_targets_long.csv', index=False)
# leave-one-test-donor-out on the locked target (all methods, mean of states), to expose the five-donor limit
rows = []
for dr in TEST:
    keep = [d for d in TEST if d != dr]; v = {s: scores(*RR('locked_3v2', s), rz(AT([d for d in keep if COND[d] == 1], s).mean(0) - AT([d for d in keep if COND[d] == 0], s).mean(0))) for s in ST}
    for j, (i, rd) in enumerate(rows_idx):
        rows.append(dict(dropped_test_donor=dr, dropped_condition='PD' if COND[dr] == 1 else 'control', method=REG[i]['method'], readout=rd, tier=REG[i]['tier'], rho_astro=v['astro'][j], rho_microglia_mhc2=v['microglia_mhc2'][j], rho_mean=(v['astro'][j] + v['microglia_mhc2'][j]) / 2))
pd.DataFrame(rows).to_csv(OUT / 'TABLE_locked_leave_one_test_donor.csv', index=False)
pickle.dump(dict(rows_idx=[(REG[i]['method'], rd) for i, rd in rows_idx]), open(SRC / 'bm_index.pkl', 'wb'))
# every ablation arm under every readout (locked and external targets), with seed agreement and prediction size: the single-seed and pooled readouts can disagree
def seedstat(D, seeds):
    o = {}
    for st in ST:
        ax = [a for a in AX if a.endswith('__' + st)]; rr = []; ms = []
        for a in ax:
            U = np.stack([D[(s_, a)] for s_ in seeds]).astype(float); ms.append(float(np.sqrt((U ** 2).mean(1)).mean()))
            if len(seeds) > 1:
                V = rz(U); q = V @ V.T; rr.append(float(q[np.triu_indices(len(seeds), 1)].mean()))
        o[st] = (float(np.mean(rr)) if rr else np.nan, float(np.mean(ms)))
    return o
AB = [(n, ARMD[n][0], ARMD[n][1]) for n in ['random_pairing', 'coupling_off', 'no_blood', 'no_cond2', 'no_cond2_random', 'label_shuffle']] + [('PeriScope, seeds 42-44', PS, S3), ('PeriScope, 7 seeds never repaired (42-46, 50, 51)', PS, S7), ('PeriScope, seeds 42-51', PS, S10)]
rows = []; F3 = L[L.method == 'PeriScope, seeds 42-44'].set_index(['target', 'readout'])
for n, D, sd in AB:
    q = seedstat(D, sd)
    for tn in ['locked_3v2', 'smajic_midbrain', 'gse184950', 'ma_substantia_nigra']:
        for rd in ['single_seed_mean', 'ensemble_rawdelta_per_axis', 'ensemble_z']:
            x = L[(L.method == n) & (L.target == tn) & (L.readout == rd)].iloc[0]; f = F3.loc[(tn, rd)]
            rows.append(dict(arm=n, tier=x.tier, n_seeds=len(sd), seeds=' '.join(map(str, sd)), target=tn, readout=rd, rho_astro=x.rho_astro, rho_microglia_mhc2=x.rho_microglia_mhc2, rho_mean_of_states=x.rho_mean_of_states, full_seeds_42_44_same_readout_astro=f.rho_astro, full_seeds_42_44_same_readout_microglia_mhc2=f.rho_microglia_mhc2, full_seeds_42_44_same_readout_mean=f.rho_mean_of_states, diff_mean_vs_full_seeds_42_44=x.rho_mean_of_states - f.rho_mean_of_states, perm_p_one_sided_rho=x.perm_p_one_sided_rho, n_label_permutations=x.n_label_permutations, mean_interseed_spearman_astro=q['astro'][0], mean_interseed_spearman_microglia_mhc2=q['microglia_mhc2'][0], predicted_delta_rms_astro=q['astro'][1], predicted_delta_rms_microglia_mhc2=q['microglia_mhc2'][1]))
pd.DataFrame(rows).to_csv(OUT / 'TABLE_ablation_by_readout.csv', index=False); print('ablation by readout rows', len(rows), flush=True)
print('done', int(time.time() - t0), flush=True)
