#!/usr/bin/env python
"""pd-v3 final target funnel (gates decided by the user on 2026-11-19).

One command rebuilds every output:
    /public/home/mengxl/dzy/envs/pdenv/bin/python build_target_funnel_final.py
Alternative program-direction source (e.g. 15-seed direction), written elsewhere:
    ... build_target_funnel_final.py --direction-csv FILE --direction-key combo \
        --direction-col restorative_direction_15seed --out .../with_15_seed_direction

Gates per antigen, in order. Gates 1-5 hold on the same program.
 1 C1 identity and source-cell context      2 C2 measured dose
 3 C3 specific model response             4 C4 matched-cell-state PD brain support
 5 C5 surface/extracellular accessibility
 6 MR: pair BH q<0.05 (atlas global BH, drop policy), or the pair is the assay's
   best trait and the assay has maxT p<0.05 or assay q<0.05
 7 Coloc on the same pair: unified_H4>=0.8 on AGES or UKB (best accession per pair)
 8 Direction-matched drug. Model direction over the antigen's gate 1-5 programs,
   ignoring unresolved; all increase or all decrease, else fail (mixed / none).
   Primary (expanded): usable record (drug lane), route direct, ligand_route,
   complex_partner or gene_silencing, approved or clinical, not withdrawn,
   up<->increase or down<->decrease. Direct-only subset kept as a labelled column.
No heatmap.
"""
import argparse, json, os, sys, time
import numpy as np
import pandas as pd

ROOT = '/public/home/mengxl/dzy/pd-v3'
P1 = ROOT + '/papers/01_mr_imaging_scan'
P3 = ROOT + '/papers/03_world_model_integration'
RES = P3 + '/results'
SIX = '/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003/six_criteria_20261003'
IN = dict(
    programs=SIX + '/all_6849_programs_six_criteria.csv',
    links=RES + '/cross_domain_inputs_20261117/all_model_target_assay_links.csv',
    mr_summary=RES + '/new_mr_integration_20261117/all_model_target_MR_summary.csv',
    mr_rows=RES + '/new_mr_integration_20261117/all_linked_assay_imaging_rows_drop.csv.gz',
    phenotypes=RES + '/new_mr_integration_20261117/all_3963_imaging_phenotypes.csv',
    hits=RES + '/research_20261118/mr_imaging/hits_q05_drop.csv.gz',
    maxT=RES + '/research_20261118/integration/atlas/assay_maxT.csv',
    assay_cal=RES + '/research_20261118/mr_imaging/assay_level_calibration.csv.gz',
    region=P1 + '/results/all_target_delivery_20261119/all_target_region_table.csv.gz',
    coloc=P1 + '/results/all_target_susie_coloc_20261119/susie_coloc_pairs.csv.gz',
    drugs=RES + '/drug_evidence_expanded_20261119/drug_records.csv.gz',
    drug_summary=RES + '/drug_evidence_expanded_20261119/antigen_drug_summary.csv',
)
DEFAULT_OUT = RES + '/target_funnel_final_20261119'
EXPANDED_ROUTES = ['direct', 'ligand_route', 'complex_partner', 'gene_silencing']
H4_MIN = 0.8
MATCH = {'increase': 'up', 'decrease': 'down'}
OPP = {'increase': 'down', 'decrease': 'up'}

ap = argparse.ArgumentParser()
ap.add_argument('--out', default=DEFAULT_OUT)
ap.add_argument('--direction-csv', default=None, help='program-level direction table replacing predicted_restorative_protein_change')
ap.add_argument('--direction-key', default='combo')
ap.add_argument('--direction-col', default='predicted_restorative_protein_change')
ap.add_argument('--direction-label', default=None)
a = ap.parse_args()
OUT = a.out
os.makedirs(OUT, exist_ok=True)
DIR_LABEL = a.direction_label or ('frozen six-criteria predicted_restorative_protein_change (confirmation seeds 47-51)'
                                  if a.direction_csv is None else f'{a.direction_csv}:{a.direction_col}')
truth = lambda s: s.astype(str).str.strip().str.lower().isin(['true', '1', 'yes'])
log = []
def note(*x):
    s = ' '.join(str(v) for v in x); print(s, flush=True); log.append(s)

# ---------------------------------------------------------------- programs, gates 1-5
E = pd.read_csv(IN['programs'], low_memory=False)
assert E.combo.is_unique, 'combo is not a unique program key'
for k in range(1, 6):
    E[f'p_C{k}'] = E[f'C{k}'].astype(str).str.strip().str.lower().eq('pass')
E['cum1'] = E.p_C1
for k in range(2, 6):
    E[f'cum{k}'] = E[f'cum{k-1}'] & E[f'p_C{k}']
E['pass_1to5'] = E.cum5
if a.direction_csv:
    Dd = pd.read_csv(a.direction_csv, low_memory=False)[[a.direction_key, a.direction_col]].drop_duplicates(a.direction_key)
    E = E.merge(Dd.rename(columns={a.direction_key: 'combo', a.direction_col: '_dir'}), on='combo', how='left')
    note('direction override:', a.direction_csv, 'programs matched', int(E._dir.notna().sum()), 'of', len(E),
         '; gate1-5 programs without a value:', int((E.pass_1to5 & E._dir.isna()).sum()))
else:
    E['_dir'] = E.predicted_restorative_protein_change
E['model_dir'] = E._dir.astype(str).str.strip().str.lower().where(E._dir.notna(), 'unresolved')
E.loc[~E.model_dir.isin(['increase', 'decrease']), 'model_dir'] = 'unresolved'

antigens = sorted(E.protein.unique())
A = pd.DataFrame({'antigen': antigens})
A = A.merge(E.groupby('protein').agg(genes=('genes', 'first'), identity_mapping=('identity_mapping', 'first'),
                                     n_programs=('combo', 'size')).reset_index().rename(columns={'protein': 'antigen'}), on='antigen')
for k in range(1, 6):
    g = E.groupby('protein')[f'cum{k}'].sum()
    A[f'n_prog_gate1to{k}'] = A.antigen.map(g).fillna(0).astype(int)
    A[f'gate{k}'] = A[f'n_prog_gate1to{k}'] > 0
# individual (non-cumulative) criteria, for context
for k in range(1, 6):
    A[f'any_program_C{k}'] = A.antigen.map(E.groupby('protein')[f'p_C{k}'].any()).astype(bool)

P5 = E[E.pass_1to5]
dcount = P5.groupby(['protein', 'model_dir']).size().unstack(fill_value=0)
for c in ['increase', 'decrease', 'unresolved']:
    A[f'n_prog_1to5_{c}'] = A.antigen.map(dcount[c] if c in dcount else pd.Series(dtype=int)).fillna(0).astype(int)
def resolve(r):
    if not r.gate5: return 'no_gate1to5_program'
    i, d = r.n_prog_1to5_increase, r.n_prog_1to5_decrease
    if i and d: return 'mixed'
    if i: return 'increase'
    if d: return 'decrease'
    return 'no_resolved_program'
A['model_restorative_direction'] = A.apply(resolve, axis=1)
prog_txt = P5.assign(t=P5.combo + '[' + P5.model_dir.str[:3] + ']').groupby('protein').t.apply(lambda s: '; '.join(sorted(s)))
A['programs_gate1to5'] = A.antigen.map(prog_txt).fillna('')
A['axes_gate1to5'] = A.antigen.map(P5.groupby('protein')['axis'].apply(lambda s: '; '.join(sorted(set(s))))).fillna('')

# ---------------------------------------------------------------- MR pairs, gate 6
L = pd.read_csv(IN['links'], keep_default_na=False)
L = L[(~truth(L.is_control)) & (L.exposure_id != '')][['model_target', 'exposure_id', 'UKB_gene']].drop_duplicates()
L = L.rename(columns={'model_target': 'antigen'})
A['linked_exposures'] = A.antigen.map(L.groupby('antigen').exposure_id.apply(lambda s: ';'.join(sorted(s)))).fillna('')
A['has_linked_assay'] = A.linked_exposures != ''
ex = set(L.exposure_id)

H = pd.read_csv(IN['hits'], low_memory=False)
H = H[H.exposure_id.isin(ex)]
H = H[['exposure_id', 'phenotype', 'trait', 'family', 'variant', 'beta', 'se', 'p', 'q_bh_global', 'F']].assign(mr_pair_q_pass=True)
MX = pd.read_csv(IN['maxT'])[['exposure_id', 'best_phenotype', 'p_maxT', 'q_maxT_all_assays']]
CAL = pd.read_csv(IN['assay_cal'])[['exposure_id', 'assay_q', 'assay_p']]
AS = MX.merge(CAL, on='exposure_id', how='left')
AS = AS[AS.exposure_id.isin(ex)]
B = AS[(AS.p_maxT < 0.05) | (AS.assay_q < 0.05)].rename(columns={'best_phenotype': 'phenotype'})[['exposure_id', 'phenotype']].assign(mr_best_trait_pass=True)
M = H.merge(B, on=['exposure_id', 'phenotype'], how='outer')
M['mr_pair_q_pass'] = M.mr_pair_q_pass.fillna(False).astype(bool)
M['mr_best_trait_pass'] = M.mr_best_trait_pass.fillna(False).astype(bool)
M = M.merge(AS.rename(columns={'best_phenotype': 'assay_best_phenotype'}), on='exposure_id', how='left')
M['is_assay_best_trait'] = M.phenotype == M.assay_best_phenotype
# MR estimates for best-trait pairs that are not atlas hits
need = M[M.beta.isna()][['exposure_id', 'phenotype']]
if len(need):
    R = pd.read_csv(IN['mr_rows'], low_memory=False, usecols=['exposure_id', 'phenotype', 'trait', 'family', 'variant', 'beta', 'se', 'p', 'F'])
    R = R.merge(need, on=['exposure_id', 'phenotype']).drop_duplicates(['exposure_id', 'phenotype'])
    M = M.set_index(['exposure_id', 'phenotype'])
    R = R.set_index(['exposure_id', 'phenotype'])
    for c in ['trait', 'family', 'variant', 'beta', 'se', 'p', 'F']:
        M[c] = M[c].fillna(R[c]) if c in R else M[c]
    M = M.reset_index()
PH = pd.read_csv(IN['phenotypes'], usecols=['phenotype', 'description']).drop_duplicates('phenotype')
M = M.merge(PH.rename(columns={'description': 'phenotype_description'}), on='phenotype', how='left')
M['mr_route'] = np.select([M.mr_pair_q_pass & M.mr_best_trait_pass, M.mr_pair_q_pass], ['pair_q_and_best_trait', 'pair_q'], 'best_trait_maxT_or_assay_q')

# ---------------------------------------------------------------- coloc, gate 7
S = pd.read_csv(IN['coloc'], low_memory=False)
S = S[S.is_best_accession_for_pair == True]
keep = ['unified_H4', 'unified_source', 'H4_single', 'H4_susie_mr_signal', 'status', 'status_reason', 'ld_qc_flag', 'accession']
W = None
for plat in ['AGES', 'UKB']:
    s = S[S.platform == plat][['exposure_id', 'phenotype'] + keep].drop_duplicates(['exposure_id', 'phenotype'])
    s = s.rename(columns={c: f'{c}_{plat}' for c in keep})
    W = s if W is None else W.merge(s, on=['exposure_id', 'phenotype'], how='outer')
W['in_coloc_table'] = True
M = M.merge(W, on=['exposure_id', 'phenotype'], how='left')
M['in_coloc_table'] = M.in_coloc_table.fillna(False).astype(bool)
M['H4_max'] = M[['unified_H4_AGES', 'unified_H4_UKB']].max(axis=1)
M['coloc_pass'] = M.H4_max >= H4_MIN
M['coloc_pass_both_platforms'] = (M.unified_H4_AGES >= H4_MIN) & (M.unified_H4_UKB >= H4_MIN)
M['coloc_status'] = np.select([M.coloc_pass, M.H4_max.notna()], ['pass', 'fail_evaluated'], 'undefined')
# annotations from the region table (not gates)
RT = pd.read_csv(IN['region'], low_memory=False, usecols=['exposure_id', 'phenotype', 'phenotype_description', 'H4_w1000k_1e5', 'selection_p_H4', 'reverse_state', 'repro_available'])
RT = RT.sort_values('H4_w1000k_1e5', ascending=False).drop_duplicates(['exposure_id', 'phenotype'])
M['phenotype_description'] = M.phenotype_description.fillna(M.phenotype.map(RT.drop_duplicates('phenotype').set_index('phenotype').phenotype_description)).fillna(M.trait)
M = M.merge(RT[['exposure_id', 'phenotype', 'selection_p_H4', 'reverse_state', 'repro_available']].rename(
    columns={'selection_p_H4': 'annot_selection_p_H4', 'reverse_state': 'annot_reverse_state', 'repro_available': 'annot_repro_available'}),
    on=['exposure_id', 'phenotype'], how='left')
M = L.merge(M, on='exposure_id', how='inner')
M = M.sort_values(['antigen', 'coloc_pass', 'H4_max', 'p'], ascending=[True, False, False, True])

g6 = M.groupby('antigen').size()
A['n_mr_pairs'] = A.antigen.map(g6).fillna(0).astype(int)
A['any_mr_pair'] = A.n_mr_pairs > 0
A['n_mr_pairs_coloc_pass'] = A.antigen.map(M[M.coloc_pass].groupby('antigen').size()).fillna(0).astype(int)
A['n_mr_pairs_coloc_undefined'] = A.antigen.map(M[M.coloc_status == 'undefined'].groupby('antigen').size()).fillna(0).astype(int)
A['any_mr_coloc_pair'] = A.n_mr_pairs_coloc_pass > 0
A['gate6'] = A.gate5 & A.any_mr_pair
A['gate7'] = A.gate6 & A.any_mr_coloc_pair
A['gate7_status'] = np.select([~A.gate6, A.gate7, A.n_mr_pairs_coloc_undefined == A.n_mr_pairs, A.n_mr_pairs_coloc_undefined > 0],
                              ['not_reached', 'pass', 'undefined_all_pairs_not_evaluable', 'fail_some_pairs_undefined'], 'fail_all_pairs_evaluated')

# ---------------------------------------------------------------- drugs, gate 8
D = pd.read_csv(IN['drugs'], low_memory=False)
D['withdrawn_b'] = truth(D.withdrawn)
D['usable_b'] = truth(D.usable)
base = D.route_class.isin(EXPANDED_ROUTES) & D.direction.isin(['up', 'down'])
D['qual_expanded'] = base & D.usable_b & ~D.withdrawn_b & D.level.isin(['approved', 'clinical'])
D['qual_direct'] = D.qual_expanded & (D.route_class == 'direct')
D['qual_basic_not_usable'] = base & ~D.usable_b & ~D.withdrawn_b & D.level.isin(['approved', 'clinical'])
D['qual_preclinical'] = base & (D.level == 'preclinical_tool') & ~D.withdrawn_b
dmap = A.set_index('antigen').model_restorative_direction
D['model_dir'] = D.antigen.map(dmap)
D['matched'] = D.direction == D.model_dir.map(MATCH)
D['opposite'] = D.direction == D.model_dir.map(OPP)
lvl_rank = {'approved': 2, 'clinical': 1}

def drug_txt(d):
    d = d.sort_values(['level', 'max_phase', 'drug'], ascending=[True, False, True])
    return '; '.join(f'{r.drug} ({r.route_class}, {r.direction}, {r.level}' + (f' ph{int(r.max_phase)}' if pd.notna(r.max_phase) and r.level == 'clinical' else '') + ')'
                     for r in d.drop_duplicates(['drug', 'route_class']).itertuples())
def summarize(mask, prefix):
    d = D[mask]
    if d.empty:
        A[f'{prefix}_n'] = 0
        for c in ['drugs', 'routes', 'best_level']: A[f'{prefix}_{c}'] = ''
        return
    A[f'{prefix}_n'] = A.antigen.map(d.groupby('antigen').drug.nunique()).fillna(0).astype(int)
    A[f'{prefix}_drugs'] = A.antigen.map(d.groupby('antigen').apply(drug_txt)).fillna('')
    A[f'{prefix}_routes'] = A.antigen.map(d.groupby('antigen').route_class.apply(lambda s: ';'.join(sorted(set(s))))).fillna('')
    A[f'{prefix}_best_level'] = A.antigen.map(d.assign(r=d.level.map(lvl_rank)).groupby('antigen').r.max().map({2: 'approved', 1: 'clinical'})).fillna('')
summarize(D.qual_expanded & D.matched, 'matched_expanded')
summarize(D.qual_direct & D.matched, 'matched_direct')
summarize(D.qual_expanded & D.opposite, 'opposite_expanded')
summarize(D.qual_preclinical & D.matched, 'matched_preclinical_tool')
summarize(D.qual_basic_not_usable & D.matched, 'matched_but_not_usable')
A['matched_but_not_usable_reasons'] = A.antigen.map(D[D.qual_basic_not_usable & D.matched].groupby('antigen').not_usable_reason.apply(
    lambda s: ';'.join(f'{k}:{v}' for k, v in s.value_counts().items()))).fillna('')
A['usable_any_direction_n'] = A.antigen.map(D[D.qual_expanded].groupby('antigen').drug.nunique()).fillna(0).astype(int)
A['usable_directions'] = A.antigen.map(D[D.qual_expanded].groupby('antigen').direction.apply(lambda s: ';'.join(sorted(set(s))))).fillna('')
A['clinical_depletion_n'] = A.antigen.map(D[(D.direction == 'depletion') & D.level.isin(['approved', 'clinical'])].groupby('antigen').drug.nunique()).fillna(0).astype(int)
A['clinical_unknown_direction_n'] = A.antigen.map(D[(D.direction == 'unknown') & D.level.isin(['approved', 'clinical'])].groupby('antigen').drug.nunique()).fillna(0).astype(int)
A['direction_resolved'] = A.model_restorative_direction.isin(['increase', 'decrease'])
A['gate8a_direction'] = A.gate7 & A.direction_resolved
A['gate8_expanded'] = A.gate8a_direction & (A.matched_expanded_n > 0)
A['gate8_direct_only'] = A.gate8a_direction & (A.matched_direct_n > 0)
def g8_reason(r):
    if not r.gate7: return ''
    if r.model_restorative_direction == 'mixed': return 'mixed_model_direction'
    if r.model_restorative_direction == 'no_resolved_program': return 'no_resolved_program_direction'
    if r.matched_expanded_n > 0: return 'pass'
    if r.opposite_expanded_n > 0: return 'only_opposite_direction_drugs'
    if r.clinical_depletion_n > 0 or r.clinical_unknown_direction_n > 0: return 'only_depletion_or_unknown_direction_drugs'
    return 'no_clinical_drug_record'
A['gate8_reason'] = A.apply(g8_reason, axis=1)

# ---------------------------------------------------------------- outputs
def first_fail(r):
    for k, lab in [(1, 'C1'), (2, 'C2'), (3, 'C3'), (4, 'C4'), (5, 'C5'), (6, 'MR'), (7, 'coloc')]:
        if not r[f'gate{k}']: return lab
    if not r.gate8a_direction: return 'direction'
    if not r.gate8_expanded: return 'drug'
    return 'pass'
A['first_failed_gate'] = A.apply(first_fail, axis=1)
A['final_target'] = A.gate8_expanded
A['final_target_direct_only'] = A.gate8_direct_only
A.to_csv(f'{OUT}/antigen_gate_matrix.csv', index=False)

pair_cols = ['antigen', 'UKB_gene', 'exposure_id', 'phenotype', 'phenotype_description', 'trait', 'family', 'mr_route', 'beta', 'se', 'p', 'q_bh_global', 'F',
             'is_assay_best_trait', 'p_maxT', 'assay_q', 'variant', 'in_coloc_table', 'H4_max', 'coloc_status', 'coloc_pass_both_platforms',
             'unified_H4_AGES', 'unified_source_AGES', 'H4_single_AGES', 'H4_susie_mr_signal_AGES', 'status_AGES', 'accession_AGES', 'ld_qc_flag_AGES',
             'unified_H4_UKB', 'unified_source_UKB', 'H4_single_UKB', 'H4_susie_mr_signal_UKB', 'status_UKB', 'ld_qc_flag_UKB',
             'status_reason_AGES', 'status_reason_UKB', 'annot_selection_p_H4', 'annot_reverse_state', 'annot_repro_available']
M = M.merge(A[['antigen', 'gate5', 'model_restorative_direction']], on='antigen', how='left')
M[pair_cols[:1] + ['gate5', 'model_restorative_direction'] + pair_cols[1:]].to_csv(f'{OUT}/mr_coloc_pairs_all_antigens.csv', index=False)

F = A[A.final_target].copy()
best = M[M.coloc_pass].drop_duplicates('antigen').set_index('antigen')
pairs_txt = M[M.coloc_pass].assign(t=lambda d: d.phenotype + ' ' + d.trait.astype(str) + ' (H4 A=' + d.unified_H4_AGES.round(3).astype(str) + ', U=' + d.unified_H4_UKB.round(3).astype(str) + ')').groupby('antigen').t.apply(lambda s: '; '.join(s))
F['n_mr_coloc_pairs'] = F.n_mr_pairs_coloc_pass
F['mr_coloc_pairs'] = F.antigen.map(pairs_txt)
for c in ['exposure_id', 'phenotype', 'phenotype_description', 'trait', 'mr_route', 'beta', 'se', 'p', 'q_bh_global', 'p_maxT', 'assay_q',
          'unified_H4_AGES', 'unified_H4_UKB', 'H4_max', 'coloc_pass_both_platforms', 'unified_source_AGES', 'unified_source_UKB', 'annot_selection_p_H4', 'annot_reverse_state']:
    F[f'top_pair_{c}'] = F.antigen.map(best[c])
F['any_pair_coloc_both_platforms'] = F.antigen.map(M[M.coloc_pass].groupby('antigen').coloc_pass_both_platforms.any())
F['passes_direct_only'] = F.gate8_direct_only
fcols = ['antigen', 'genes', 'model_restorative_direction', 'n_prog_gate1to5', 'n_prog_1to5_increase', 'n_prog_1to5_decrease', 'n_prog_1to5_unresolved',
         'axes_gate1to5', 'programs_gate1to5', 'matched_expanded_n', 'matched_expanded_best_level', 'matched_expanded_routes', 'matched_expanded_drugs',
         'passes_direct_only', 'matched_direct_n', 'matched_direct_drugs', 'opposite_expanded_n', 'opposite_expanded_drugs',
         'matched_preclinical_tool_drugs', 'matched_but_not_usable_drugs', 'n_mr_coloc_pairs', 'any_pair_coloc_both_platforms'] + \
        [c for c in F.columns if c.startswith('top_pair_')] + ['mr_coloc_pairs', 'linked_exposures']
F = F.sort_values(['passes_direct_only', 'top_pair_H4_max'], ascending=False)
F[fcols].to_csv(f'{OUT}/final_targets.csv', index=False)
DM = D[D.qual_expanded & D.matched & D.antigen.isin(F.antigen)].copy()
DM['counts_in_direct_only'] = DM.route_class == 'direct'
DM[['antigen', 'model_dir', 'drug', 'chembl_id', 'drug_type', 'route_class', 'acted_gene', 'link', 'link_quality', 'direction', 'direction_basis',
    'mechanism_text', 'level', 'max_phase', 'withdrawn', 'usable_tier', 'counts_in_direct_only', 'sources', 'references', 'evidence_links']
   ].sort_values(['antigen', 'route_class', 'level', 'drug']).to_csv(f'{OUT}/final_target_drug_matches.csv', index=False)
M[M.antigen.isin(F.antigen) & M.coloc_pass][pair_cols].to_csv(f'{OUT}/final_target_mr_coloc_pairs.csv', index=False)

NM = A[A.gate7 & ~A.gate8_expanded].copy()
NM[['antigen', 'genes', 'gate8_reason', 'model_restorative_direction', 'n_prog_1to5_increase', 'n_prog_1to5_decrease', 'n_prog_1to5_unresolved',
    'opposite_expanded_drugs', 'matched_but_not_usable_drugs', 'matched_but_not_usable_reasons', 'matched_preclinical_tool_drugs',
    'clinical_depletion_n', 'clinical_unknown_direction_n', 'usable_directions', 'n_mr_pairs_coloc_pass', 'axes_gate1to5', 'programs_gate1to5']
   ].to_csv(f'{OUT}/near_misses.csv', index=False)

steps = [('0', 'all PeriScope antigens', pd.Series(True, index=A.index), len(E)),
         ('1', 'C1 identity and source-cell context', A.gate1, int(E.cum1.sum())),
         ('2', 'C2 measured dose (same program)', A.gate2, int(E.cum2.sum())),
         ('3', 'C3 specific model response (same program)', A.gate3, int(E.cum3.sum())),
         ('4', 'C4 matched-cell-state PD brain support (same program)', A.gate4, int(E.cum4.sum())),
         ('5', 'C5 surface/extracellular accessibility (same program)', A.gate5, int(E.cum5.sum())),
         ('6', 'MR: pair BH q<0.05 or best trait with maxT p<0.05 / assay q<0.05', A.gate6, None),
         ('7', 'coloc on the same pair: unified_H4>=0.8 on AGES or UKB', A.gate7, None),
         ('8a', 'single resolved restorative direction (mixed fails)', A.gate8a_direction, None),
         ('8', 'PRIMARY: direction-matched drug, expanded routes', A.gate8_expanded, None),
         ('8_direct', 'labelled subset: direction-matched drug, direct route only', A.gate8_direct_only, None)]
rows, prev = [], len(A)
for sid, lab, m, nprog in steps:
    n = int(m.sum())
    rows.append(dict(step=sid, gate=lab, antigens=n, removed_from_previous=(prev - n) if sid not in ('0', '8_direct') else '',
                     programs=nprog if nprog is not None else '', antigens_list=';'.join(sorted(A.antigen[m])) if n <= 40 else ''))
    if sid != '8_direct': prev = n
FC = pd.DataFrame(rows)
FC.loc[FC.step == '8_direct', 'removed_from_previous'] = int(A.gate8_expanded.sum() - A.gate8_direct_only.sum())
FC['level'] = np.where(FC.step == '8_direct', 'strict (direct only)', np.where(FC.step == '8', 'expanded (primary)', 'shared'))
FC.to_csv(f'{OUT}/funnel_counts.csv', index=False)

# context counts for notes
ctx = {}
ctx['gate7_status_among_gate6'] = A[A.gate6].gate7_status.value_counts().to_dict()
ctx['dir_counts_gate5'] = A[A.gate5].model_restorative_direction.value_counts().to_dict()
ctx['dir_counts_gate7'] = A[A.gate7].model_restorative_direction.value_counts().to_dict()
ctx['g8_reasons'] = A[A.gate7].gate8_reason.value_counts().to_dict()
ctx['mr_pairs'] = dict(total=int(len(M)), pairs_gate5_antigens=int(M.gate5.sum()), coloc_status=M.coloc_status.value_counts().to_dict(),
                       coloc_status_gate5=M[M.gate5].coloc_status.value_counts().to_dict(), route=M.mr_route.value_counts().to_dict())
ctx['finals'] = F[['antigen', 'model_restorative_direction', 'matched_expanded_n', 'matched_expanded_best_level', 'matched_expanded_routes',
                   'passes_direct_only', 'top_pair_H4_max', 'any_pair_coloc_both_platforms']].to_dict('records')

# ---------------------------------------------------------------- headline metric candidates (descriptive, with denominators)
hm = []
def add(key, num, den, what, caution, ref=''):
    hm.append(dict(metric=key, numerator=int(num), denominator=int(den), percent=round(100 * num / den, 1) if den else np.nan,
                   what_it_measures=what, comparison_or_null=ref, how_to_read=caution))
# background: all atlas MR-hit pairs (any exposure), coloc pass among pairs with a defined H4
Sb = W.merge(pd.read_csv(IN['hits'], usecols=['exposure_id', 'phenotype']), on=['exposure_id', 'phenotype'])
Sb['H4_max'] = Sb[['unified_H4_AGES', 'unified_H4_UKB']].max(axis=1)
bg_ev = Sb[Sb.H4_max.notna()]
bg_txt = f'all atlas MR-hit pairs: {int((bg_ev.H4_max >= H4_MIN).sum())}/{len(bg_ev)} = {100 * (bg_ev.H4_max >= H4_MIN).mean():.1f}% among pairs with a defined H4'
bg_ex = bg_ev.groupby('exposure_id').H4_max.max()
bg_ex_txt = f'all atlas assays with an MR hit and a defined H4: {int((bg_ex >= H4_MIN).sum())}/{len(bg_ex)} = {100 * (bg_ex >= H4_MIN).mean():.1f}% colocalise (assay level, best-trait pairs not included)'
# background: MR support among linked antigens outside the model gates
la = A[A.has_linked_assay]
mr_bg = f'linked antigens failing C1-C5: {int(la[~la.gate5].any_mr_pair.sum())}/{int((~la.gate5).sum())}'
# null for drug-direction agreement: permute model directions among the same antigens
rng = np.random.default_rng(20261119)
def agree_null(u):
    if not len(u): return ''
    has_up = (u.usable_directions.str.contains('up')).values; has_dn = (u.usable_directions.str.contains('down')).values
    dirs = u.model_restorative_direction.values.copy(); obs = int(((dirs == 'increase') & has_up | (dirs == 'decrease') & has_dn).sum())
    sims = np.array([int(((p == 'increase') & has_up | (p == 'decrease') & has_dn).sum()) for p in (rng.permutation(dirs) for _ in range(20000))])
    exp_half = float(np.where(has_up & has_dn, 1.0, 0.5).sum())
    return f'permuted directions: mean {sims.mean():.1f}/{len(u)} ({100 * sims.mean() / len(u):.1f}%), P(>=obs {obs})={((sims >= obs).sum() + 1) / (len(sims) + 1):.3f}; coin-flip expectation {exp_half:.1f}'
add('panel_with_model_program', A.gate5.sum(), len(A), 'antigens with >=1 program passing C1-C5 on the same program', 'model-side gates only; C4 not yet re-verified')
g5 = A[A.gate5]; g5r = g5[g5.direction_resolved]
add('gate5_single_direction', len(g5r), len(g5), 'antigens passing C1-C5 whose resolved programs all point one way', 'unresolved programs ignored')
u = g5r[g5r.usable_any_direction_n > 0]
add('drug_direction_agreement_gate5', (u.matched_expanded_n > 0).sum(), len(u), 'gate1-5 antigens with a usable directional drug that have >=1 drug in the model restorative direction', 'existing drugs were developed for other diseases; agreement is availability, not PD efficacy', agree_null(u))
g6a = A[A.gate6]
add('gate6_antigens_colocalised', g6a.gate7.sum(), len(g6a), 'model-nominated antigens with MR support that colocalise (H4>=0.8) on the same pair', 'antigens with all pairs undefined count in the denominator')
g6e = g6a[g6a.gate7_status != 'undefined_all_pairs_not_evaluable']
add('gate6_antigens_colocalised_evaluable', g6e.gate7.sum(), len(g6e), 'same, among antigens with at least one pair whose H4 is defined', 'undefined antigens removed from the denominator', bg_ex_txt)
add('mr_support_among_nominated_linked', la[la.gate5].any_mr_pair.sum(), int(la.gate5.sum()), 'C1-C5 antigens with a linked UKB-PPP assay that have an MR pair', 'low: most assays have no MR pair', mr_bg)
g7 = A[A.gate7]
add('gate7_with_matched_drug', g7.gate8_expanded.sum(), len(g7), 'antigens passing gates 1-7 that also pass gate 8 (expanded)', 'gate 8 primary definition')
g7r = g7[g7.direction_resolved]
add('gate7_resolved_with_matched_drug', g7r.gate8_expanded.sum(), len(g7r), 'gate1-7 antigens with a single resolved direction that have a matched drug', 'excludes mixed / unresolved')
if len(F):
    add('final_with_approved_drug', (F.matched_expanded_best_level == 'approved').sum(), len(F), 'final targets whose matched drugs include an approved drug', 'approved for some indication, not PD')
    add('final_direct_only', F.passes_direct_only.sum(), len(F), 'final targets that also pass with direct-route drugs only', 'strict subset')
    add('final_coloc_both_platforms', F.any_pair_coloc_both_platforms.fillna(False).sum(), len(F), 'final targets with a pair at H4>=0.8 on both AGES and UKB', 'platforms are separate cohorts and assays')
    add('final_top_pair_H4_ge_0.9', (F.top_pair_H4_max >= 0.9).sum(), len(F), 'final targets whose best pair has unified H4>=0.9', 'unified_H4 = max(single, SuSiE MR-signal)')
mp = M[M.gate5]
ev = mp[mp.coloc_status != 'undefined']
add('mr_pairs_coloc_pass_among_evaluable_gate5', (ev.coloc_status == 'pass').sum(), len(ev), 'MR-passing pairs of gate1-5 antigens with H4>=0.8, among pairs with a defined H4', 'pairs without H4 are undefined, not negative', bg_txt)
HM = pd.DataFrame(hm)
HM.to_csv(f'{OUT}/headline_metric_candidates.csv', index=False)

meta = dict(script=os.path.abspath(__file__), run_time=time.strftime('%Y-%m-%d %H:%M:%S %z'), out=OUT, direction_source=DIR_LABEL,
            direction_csv=a.direction_csv, direction_key=a.direction_key, direction_col=a.direction_col,
            inputs={k: dict(path=v, size=os.path.getsize(v), mtime=time.strftime('%Y-%m-%d %H:%M', time.localtime(os.path.getmtime(v)))) for k, v in IN.items()},
            thresholds=dict(H4=H4_MIN, mr_pair_q=0.05, maxT_p=0.05, assay_q=0.05, expanded_routes=EXPANDED_ROUTES, levels=['approved', 'clinical']),
            funnel=FC[['step', 'antigens']].to_dict('records'), context=ctx, log=log)
json.dump(meta, open(f'{OUT}/run_metadata.json', 'w'), ensure_ascii=False, indent=1, default=str)
print(FC[['step', 'antigens', 'programs', 'antigens_list']].to_string(index=False))
print(json.dumps({k: v for k, v in ctx.items() if k != 'finals'}, default=str))
print(HM[['metric', 'numerator', 'denominator', 'percent', 'comparison_or_null']].to_string(index=False))

# ---------------------------------------------------------------- GATES.md (plain Chinese, regenerated on every run)
def md_table(df):
    h = '| ' + ' | '.join(df.columns) + ' |\n|' + '---|' * len(df.columns) + '\n'
    return h + ''.join('| ' + ' | '.join('' if pd.isna(v) else str(v) for v in r) + ' |\n' for r in df.itertuples(index=False))
zh = {'0': '全部抗原', '1': 'C1 身份与来源细胞情境', '2': 'C2 实测剂量', '3': 'C3 模型特异响应', '4': 'C4 匹配细胞状态的 PD 脑部支持',
      '5': 'C5 细胞表面或胞外可及', '6': 'MR', '7': 'coloc（同一配对）', '8a': '模型修复方向唯一', '8': '方向匹配的药物（扩大定义，主结果）',
      '8_direct': '方向匹配的药物（只算直接作用，标注子集）'}
ft = FC[['step', 'antigens', 'programs', 'antigens_list']].copy()
ft.insert(1, '门槛', ft.step.map(zh))
ft.columns = ['步骤', '门槛', '抗原数', '程序数（1-5 步）', '抗原']
g7s = A[A.gate6].gate7_status.value_counts().to_dict()
txt = f"""# 靶点漏斗门槛

本文件由 `build_target_funnel_final.py` 每次运行时重写。运行时间：{meta['run_time']}。输出目录：`{OUT}`。

**这些门槛是在看过初步结果之后才定下的**（包括 NRP1 这类方向不一致的靶点判为不通过、药物门槛采用扩大定义）。**目前所有结果都还没有通过统一核验。**

模型修复方向的来源：{DIR_LABEL}。

## 门槛（按顺序，逐个抗原判断）

第 1 到 5 步必须在同一个程序上同时成立。

1. C1：身份明确，且在来源细胞情境中直接测到 ADT（`C1 == pass`）。
2. C2：实测剂量（`C2 == pass`）。
3. C3：模型特异响应，冻结程序的基因诱饵 q≤0.05（`C3 == pass`）。
4. C4：匹配细胞状态的 PD 脑部支持（`C4 == pass`，沿用 six_criteria 原判定）。
5. C5：细胞表面或胞外可及（`C5 == pass`）。
6. MR：配对层面全局 BH q<0.05（图谱 drop 策略，`hits_q05_drop.csv.gz`）；或者该配对是该检测的最佳性状，且检测的 maxT p<0.05 或检测层面 q<0.05。
7. coloc：与第 6 步同一个蛋白-性状配对，`unified_H4`≥0.8，AGES 或 UKB 任一平台即可（每个配对取最佳 accession）。`unified_H4` 取单信号 H4 和携带 MR 工具变量的 SuSiE 信号 H4 中较大者。
8. 方向匹配的药物。模型修复方向取该抗原通过 1-5 步的程序的 `predicted_restorative_protein_change`，忽略 unresolved：全部为 increase 或全部为 decrease 才算有方向；两种都有（mixed）不通过，没有任何已判定程序也不通过。药物方向 up 配 increase、down 配 decrease 算匹配；清除类（depletion）和方向未知永远不算。
   - 主结果（扩大定义）：药物证据目录里 `usable=True` 的记录，路线为 direct、ligand_route、complex_partner 或 gene_silencing，已批准或进入临床，未撤市。
   - 标注子集（只算直接作用）：同上，但路线只限 direct。
   - 临床前工具药只作注释，不计入。

选择校正 p、BIG40 复现、反向 MR、平台一致性只作注释，不是门槛。

## 各步剩余

{md_table(ft)}

第 6 步通过的抗原在第 7 步的状态：{', '.join(f'{k} {v}' for k, v in g7s.items())}。`undefined_all_pairs_not_evaluable` 表示所有 MR 配对都没有可用的 H4（例如 MHC 区域或影像窗口不在磁盘上），这是“未定义”，不是阴性。

## 需要注意

- C4 沿用 six_criteria 的原判定（Kamath/GSE157783 匹配细胞状态的竞争性检验与合并 q）。`brain_gate_rebuild_20261119` 用置换检验重建另一套脑部门控时发现竞争性基因集检验在打乱标签后仍有约一半程序“显著”。这个问题是否同样影响 C4，留给统一核验。
- 药物方向是药物对靶点活性的作用方向，模型方向是细胞表面蛋白量的变化方向，两者不完全等价。已有药物大多为其他疾病开发，匹配只说明有方向一致的可用药物，不说明对 PD 有效。
- 第 6 步“最佳性状”路线的配对没有配对层面的显著性，只靠检测层面的 maxT 或 q。
"""
open(f'{OUT}/GATES.md', 'w').write(txt)
print('wrote', f'{OUT}/GATES.md')
