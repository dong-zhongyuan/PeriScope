"""Collect per-pair SuSiE-coloc results, attach the existing single-signal H4, and write the delivery tables.
Usage: s04_collect.py [output_dir]   (default: the lane directory)"""
from common import *
import sys, glob
from multiprocessing import Pool

DEST = (sys.argv[1].rstrip('/') + '/') if len(sys.argv) > 1 else OUT
os.makedirs(DEST, exist_ok=True)
RESD = WORK + 'results/'
LD_POOR = 0.2
H4_CUT = 0.8
COLOC = RES + 'all_target_coloc_20261118/'
STR = {'pair_id': str, 'gid': str, 'idp': str}


def read_result(f):
    return pd.read_csv(f, dtype={'pair_id': str}, low_memory=False)


def read_prep(f):
    return pd.read_csv(f, dtype={'pair_id': str, 'accession': str, 'region_id': str}, low_memory=False)


def production_h3(rid):
    """H3 at p12=1e-5 and 1e-6 for both production windows of one region."""
    d = pd.read_csv(COLOC + rid + '/coloc_prior_sensitivity.csv', usecols=['window_bp', 'p12_requested', 'effect_prior_sd_fraction', 'PP.H3.abf', 'PP.H4.abf'])
    d = d[d.effect_prior_sd_fraction.eq(.15)]
    o = {'region_id': rid}
    for w, lab in [(500000, '500k'), (1000000, '1000k')]:
        for p, pl in [(1e-5, '1e5'), (1e-6, '1e6')]:
            x = d[d.window_bp.eq(w) & np.isclose(d.p12_requested, p)]
            o[f'H3_{pl}_{lab}'] = float(x['PP.H3.abf'].iloc[0]) if len(x) else np.nan
            o[f'H4file_{pl}_{lab}'] = float(x['PP.H4.abf'].iloc[0]) if len(x) else np.nan
    return o


def pick_window(a500, a1000, b500, b1000):
    """Larger H4 of the two windows; returns H4, the paired H3 and the window label."""
    a500, a1000 = np.asarray(a500, float), np.asarray(a1000, float)
    use1000 = np.where(np.isnan(a500), True, np.where(np.isnan(a1000), False, a1000 >= a500))
    h4 = np.where(use1000, a1000, a500)
    h3 = np.where(use1000, np.asarray(b1000, float), np.asarray(b500, float))
    lab = np.where(np.isnan(h4), '', np.where(use1000, '1000k', '500k'))
    return h4, h3, lab


if __name__ == '__main__':
    P = pd.read_csv(WORK + 'plan_pairs.csv.gz', keep_default_na=False, low_memory=False, dtype=STR)
    for c in ['mr_hit', 'in_region_table', 'in_ukb_results', 'MHC']:
        P[c] = P[c].astype(str).eq('True')
    P['q_bh_global'] = pd.to_numeric(P.q_bh_global, errors='coerce')
    with Pool(6) as X:
        files = [f for f in glob.glob(RESD + '*/*.csv') if not f.endswith('.sig.csv')]
        R = pd.concat(X.map(read_result, files, chunksize=50), ignore_index=True)
        sig = glob.glob(RESD + '*/*.sig.csv')
        L = pd.concat(X.map(read_result, sig, chunksize=50), ignore_index=True) if sig else pd.DataFrame()
        D = pd.concat(X.map(read_prep, glob.glob(WORK + 'groups/*/pairs.csv'), chunksize=50), ignore_index=True)
        rids = sorted(set(P.region_id) - {''})
        H3 = pd.DataFrame(X.map(production_h3, rids, chunksize=20))
    assert not R.pair_id.duplicated().any() and not D.pair_id.duplicated().any()
    res_cols = [c for c in R.columns if c not in ('pair_id', 'status', 'reason', 'seconds')]
    A = P.merge(R.rename(columns={'status': 'run_status', 'reason': 'run_reason'}), on='pair_id', how='left')
    A['status'] = np.where(A.pre_status.ne(''), A.pre_status, A.run_status.fillna('not_evaluable'))
    A['status_reason'] = np.where(A.pre_status.ne(''), A.reason, np.where(A.run_status.isna(), 'no result file from the SuSiE step', A.run_reason.fillna('')))
    n_missing = int((A.pre_status.eq('') & A.run_status.isna()).sum())
    # pairs whose SuSiE run was stopped before finishing: status and reason come from an explicit record, never from a result file
    if os.path.exists(WORK + 'unfinished_pairs.csv'):
        u = pd.read_csv(WORK + 'unfinished_pairs.csv', dtype=str).set_index('pair_id')
        m = A.pre_status.eq('') & A.run_status.isna() & A.pair_id.isin(u.index)
        A.loc[m, 'status'] = A.pair_id[m].map(u.status)
        A.loc[m, 'status_reason'] = A.pair_id[m].map(u.reason)

    # existing single-signal H4
    prep = ['pair_id'] + [c for c in D.columns if c.startswith(('H4_1e', 'H3_1e'))] + ['H4_primary', 'H3_primary', 'H4_primary_1e6', 'H3_primary_1e6', 'H4_any_signal', 'cs_mass_primary', 'n_signals', 'n_imaging', 'n_imaging_source', 'n_protein']
    A = A.merge(D[prep], on='pair_id', how='left')
    T = pd.read_csv(RES + 'all_target_delivery_20261119/all_target_region_table.csv.gz', usecols=['region_id', 'H4_w1000k_1e5', 'H4_w500k_1e5', 'H4_w1000k_1e6', 'H4_w500k_1e6'], low_memory=False)
    T = T.merge(H3, on='region_id', how='left')
    h4_diff = float(np.nanmax(np.abs(np.r_[T.H4_w500k_1e5 - T.H4file_1e5_500k, T.H4_w1000k_1e5 - T.H4file_1e5_1000k])))
    A = A.merge(T.drop(columns=[c for c in T.columns if c.startswith('H4file')]).rename(columns={c: 'prod_' + c for c in T.columns if c.startswith('H3_')}), on='region_id', how='left')
    U = pd.read_csv(RM + 'ukb_cs_coloc_results.csv.gz', usecols=['exposure_id', 'phenotype', 'H4_primary', 'H3_primary', 'max_H4_any_signal']).rename(columns={'H4_primary': 'st_H4', 'H3_primary': 'st_H3', 'max_H4_any_signal': 'st_any'})
    A = A.merge(U, on=['exposure_id', 'phenotype'], how='left')
    for c in ['H4_1e5_500k', 'H4_1e5_1000k', 'H3_1e5_500k', 'H3_1e5_1000k', 'H4_1e6_500k', 'H4_1e6_1000k', 'H3_1e6_500k', 'H3_1e6_1000k']:
        if c not in A:
            A[c] = np.nan
    isA, isU, inT = A.platform.eq('AGES'), A.platform.eq('UKB'), A.region_id.ne('')
    a5, a5h, a5w = pick_window(A.H4_w500k_1e5, A.H4_w1000k_1e5, A.prod_H3_1e5_500k, A.prod_H3_1e5_1000k)
    b5, b5h, b5w = pick_window(A.H4_1e5_500k, A.H4_1e5_1000k, A.H3_1e5_500k, A.H3_1e5_1000k)
    a6, a6h, _ = pick_window(A.H4_w500k_1e6, A.H4_w1000k_1e6, A.prod_H3_1e6_500k, A.prod_H3_1e6_1000k)
    b6, b6h, _ = pick_window(A.H4_1e6_500k, A.H4_1e6_1000k, A.H3_1e6_500k, A.H3_1e6_1000k)
    useSt = isU & A.st_H4.notna()
    A['H4_single'] = np.select([isA & inT, isA, useSt, isU], [a5, b5, A.st_H4, A.H4_primary], np.nan)
    A['H3_single'] = np.select([isA & inT, isA, useSt, isU], [a5h, b5h, A.st_H3, A.H3_primary], np.nan)
    A['H4_single_p12_1e6'] = np.select([isA & inT, isA, isU], [a6, b6, A.H4_primary_1e6], np.nan)
    A['H3_single_p12_1e6'] = np.select([isA & inT, isA, isU], [a6h, b6h, A.H3_primary_1e6], np.nan)
    A['single_window'] = np.select([isA & inT, isA], [a5w, b5w], '')
    A['H4_single_any_signal'] = np.where(isU, A.st_any.where(A.st_any.notna(), A.H4_any_signal), np.nan)
    src = np.select([isA & inT, isA, useSt, isU], ['production_all_target_coloc', 'same_recipe_computed_here', 'stored_ukb_cs_coloc', 'same_recipe_computed_here'], '')
    A['single_source'] = np.where(A.H4_single.isna(), 'none', src)

    # unified value
    s, m = A.H4_single.values.astype(float), A.H4_susie_mr_signal.values.astype(float)
    A['unified_H4'] = np.where(np.isnan(s) & np.isnan(m), np.nan, np.fmax(s, m))
    A['unified_source'] = np.select([np.isnan(s) & np.isnan(m), np.isnan(m), np.isnan(s), m > s], ['none', 'single_signal', 'susie_mr_signal', 'susie_mr_signal'], 'single_signal')

    # LD flags
    def poor(x): return np.where(x.isna(), None, x > LD_POOR)
    A['ld_poor_protein'] = poor(A.lambda_protein)
    A['ld_poor_imaging'] = poor(A.lambda_imaging)
    anyp = A.lambda_protein.gt(LD_POOR) | A.lambda_imaging.gt(LD_POOR)
    A['ld_qc_flag'] = np.select([anyp, A.lambda_protein.notna() | A.lambda_imaging.notna()], ['poor', 'ok'], 'not_assessed')

    # several AGES accessions can map to one assay: rows stay per accession, one is marked per exposure-phenotype-platform
    K = ['exposure_id', 'phenotype', 'platform']
    A['n_accessions_for_pair'] = A.groupby(K).pair_id.transform('size')
    o = A.assign(_h=A.unified_H4.fillna(-1), _r=A.region_id.eq('')).sort_values(['_h', '_r', 'accession'], ascending=[False, True, True], kind='stable')
    A['is_best_accession_for_pair'] = A.pair_id.isin(o.drop_duplicates(K).pair_id)
    rg = A[isA & inT].groupby(['exposure_id', 'phenotype']).region_id.agg(lambda x: ';'.join(sorted(set(x)))).to_dict()
    A.loc[isU, 'region_id'] = [rg.get(k, '') for k in zip(A.exposure_id[isU], A.phenotype[isU])]
    A = A.rename(columns={'instrument': 'mr_instrument', 'chrom': 'instrument_chrom', 'pos37': 'instrument_pos37', 'ea': 'instrument_ea', 'oa': 'instrument_oa', 'MHC': 'instrument_in_MHC', 'img_src': 'imaging_source', 'n_snps': 'n_snps_susie', 'q_bh_global': 'mr_q_bh_global'})
    A['H4_single_ge_080'] = A.H4_single.ge(H4_CUT)
    A['unified_ge_080'] = A.unified_H4.ge(H4_CUT)
    A['instrument_in_signal'] = A.instrument_in_signal.astype('boolean').fillna(False)
    cols = ['pair_id', 'region_id', 'exposure_id', 'accession', 'gene', 'phenotype', 'family', 'platform', 'status', 'status_reason',
            'n_protein_cs', 'n_imaging_cs', 'n_protein_cs_qualifying', 'n_protein_signals_low_overlap', 'instrument_in_signal', 'instrument_basis', 'instrument_max_r2_to_cs', 'instrument_in_panel', 'instrument_in_snpset',
            'H4_susie_mr_signal', 'H3_susie_mr_signal', 'H4_susie_mr_signal_p12_1e6', 'H3_susie_mr_signal_p12_1e6', 'mr_signal_protein_cs', 'mr_signal_imaging_cs', 'mr_signal_hit_protein', 'mr_signal_hit_imaging',
            'H4_susie_any', 'H3_susie_any', 'H4_susie_any_p12_1e6', 'H3_susie_any_p12_1e6', 'any_protein_cs', 'any_imaging_cs', 'any_hit_protein', 'any_hit_imaging',
            'lambda_protein', 'lambda_imaging', 'ld_poor_protein', 'ld_poor_imaging', 'ld_qc_flag', 'n_af_discordant_protein', 'n_af_discordant_imaging', 'conv_protein', 'conv_imaging', 'min_p_protein', 'min_p_imaging', 'n_snps_susie',
            'H4_single', 'H3_single', 'H4_single_p12_1e6', 'H3_single_p12_1e6', 'H4_single_any_signal', 'single_window', 'single_source', 'unified_H4', 'unified_source',
            'mr_hit', 'mr_q_bh_global', 'in_region_table', 'in_ukb_results', 'n_accessions_for_pair', 'is_best_accession_for_pair', 'mr_instrument', 'instrument_chrom', 'instrument_pos37', 'instrument_ea', 'instrument_oa', 'instrument_in_MHC', 'imaging_source', 'n_imaging', 'n_imaging_source', 'n_protein']
    atomic_csv(A[cols].sort_values(['platform', 'gene', 'phenotype', 'accession']), DEST + 'susie_coloc_pairs.csv.gz', compression='gzip')
    if len(L):
        L = L.merge(A[['pair_id', 'platform', 'exposure_id', 'accession', 'gene', 'phenotype']], on='pair_id', how='left')
        atomic_csv(L, DEST + 'susie_coloc_signal_pairs.csv.gz', compression='gzip')

    # status summary
    rows = []
    scopes = {'all_pairs': A.pair_id.notna(), 'region_table': A.in_region_table, 'ukb_cs_results': A.in_ukb_results, 'mr_hit': A.mr_hit}
    for unit, um in [('accession_row', A.pair_id.notna()), ('best_accession_per_pair', A.is_best_accession_for_pair)]:
        for sc, sm in scopes.items():
            for (pl, st), g in A[um & sm].groupby(['platform', 'status']):
                rows.append(dict(unit=unit, scope=sc, platform=pl, status=st, n=len(g), n_instrument_in_signal=int(g.instrument_in_signal.sum()), n_H4_single_ge_080=int(g.H4_single_ge_080.sum()), n_H4_susie_mr_ge_080=int(g.H4_susie_mr_signal.ge(H4_CUT).sum()), n_H4_susie_any_ge_080=int(g.H4_susie_any.ge(H4_CUT).sum()), n_unified_ge_080=int(g.unified_ge_080.sum()), n_unified_ge_080_only_by_susie=int((g.unified_ge_080 & ~g.H4_single_ge_080).sum()), n_ld_poor=int(g.ld_qc_flag.eq('poor').sum())))
    atomic_csv(pd.DataFrame(rows), DEST + 'status_summary.csv')
    gain = A[A.unified_ge_080 & ~A.H4_single_ge_080]
    atomic_csv(gain[['platform', 'gene', 'accession', 'exposure_id', 'phenotype', 'region_id', 'mr_hit', 'H4_single', 'H4_susie_mr_signal', 'H3_susie_mr_signal', 'instrument_basis', 'lambda_protein', 'lambda_imaging', 'ld_qc_flag', 'instrument_in_MHC']].sort_values(['gene', 'phenotype']), DEST + 'pairs_passing_only_by_susie.csv')
    loss = A[A.H4_single_ge_080 & A.status.eq('ran') & A.H4_susie_mr_signal.notna() & A.H4_susie_mr_signal.lt(H4_CUT)]
    print('rows', len(A), 'results', len(R), 'signal rows', len(L), 'pre-run pairs without result', n_missing, 'group errors', len(glob.glob(RESD + '*.error')))
    print('max |table H4 - production file H4|', h4_diff)
    print(A.groupby(['platform', 'status']).size().to_string())
    print('unified source', A.groupby(['platform', 'unified_source']).size().to_dict())
    b = A[A.is_best_accession_for_pair]
    print('per pair (best accession) single>=.8', b.groupby('platform').H4_single_ge_080.sum().to_dict(), 'susie_mr>=.8', b[b.H4_susie_mr_signal.ge(H4_CUT)].groupby('platform').size().to_dict(), 'unified>=.8', b.groupby('platform').unified_ge_080.sum().to_dict())
    print('accession rows passing only by susie', len(gain), 'single>=.8 but ran with susie_mr<.8', len(loss), 'ld flag', A.ld_qc_flag.value_counts().to_dict())
