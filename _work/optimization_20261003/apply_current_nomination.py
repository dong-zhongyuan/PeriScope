"""Remove model-evaluation and antigen-panel gates from frozen evidence.

The original joint table remains the historical evidence/rule record.
current_nomination is the active nomination output. No model or test is refit.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd

POLICY_ID = '20261003_remove_model_evaluation_gates'
GATES = ['mapping_pass', 'measured_dose_pass', 'specificity_pass',
         'brain_pass', 'independent_direction_pass']
SORT = ['tier', 'q_gene_decoy', 'brain_meta_q', 'combo']
ASC = [True] * len(SORT)
PANEL_COLUMNS = ['protein_percentile', 'q_protein_decoy', 'p_protein_decoy', 'protein_selectivity']
MODEL_COLUMNS = ['protein_rho', 'positive_protein_donors', 'protein_prediction_pass',
                 'confirmation_positive_fraction', 'confirmation_monotonic_fraction',
                 'reproducible_pass', 'technical_pass']


def apply_policy(source):
    e = source.copy()
    e['tier_before_policy_revision'] = source.tier
    e['specificity_pass_before_policy_revision'] = source.specificity_pass
    e['legacy_technical_pass'] = source.technical_pass
    e['legacy_reproducible_pass'] = source.reproducible_pass
    # Compatibility field: under the current policy this tests gene-program
    # response only. Antigen-panel measurements have no decision role.
    e['specificity_pass'] = e.q_gene_decoy.le(.05)
    e['gene_program_response_pass'] = e.specificity_pass
    e['technical_pass'] = e.mapping_pass & e.measured_dose_pass
    core = e.technical_pass & e.specificity_pass & e.brain_pass
    eligible = core & e.independent_direction_pass
    e['tier'] = np.select([eligible, core], ['A', 'C'], default='D')
    if 'tier_primary_meta' in e:
        e['tier_primary_meta'] = e.tier
    flags = [('gene_or_complex_identity', 'mapping_pass'),
             ('measured_dose_range', 'measured_dose_pass'),
             ('gene_program_response', 'specificity_pass'),
             ('brain_disease_support', 'brain_pass'),
             ('independent_direction', 'independent_direction_pass')]
    e['unmet_nomination_criteria'] = [';'.join(label for label, key in flags if not row[key])
                                     for row in e.to_dict('records')]
    e['unmet_A_criteria'] = e.unmet_nomination_criteria
    e['all_nonbrain_gates_pass'] = e.technical_pass & e.specificity_pass
    incomplete = e.unassessed_evidence.fillna('').str.len().gt(0)
    e['availability_status'] = np.where(eligible, 'nominated',
                                      np.where(incomplete, 'incomplete_coverage', 'tested_below_requirements'))
    e['nomination_policy'] = POLICY_ID
    return e.sort_values(SORT, ascending=ASC, kind='stable').reset_index(drop=True)


def main(result_dir):
    root = Path(result_dir)
    out = root / 'current_nomination'
    out.mkdir(exist_ok=True)
    source_path = root / 'joint_analysis/all_program_evidence.csv'
    before_hash = hashlib.sha256(source_path.read_bytes()).hexdigest()
    source = pd.read_csv(source_path)
    e = apply_policy(source)
    registry = json.loads((root / 'candidate_registry.json').read_text())
    targets = {r['target'] for r in registry if not r['is_control']}
    assert source.combo.is_unique and set(e.combo) == set(source.combo)
    assert e[GATES].all(axis=1).equals(e.tier.isin(['A', 'B']))
    assert (e.tier.eq('A') == e[GATES].all(axis=1)).all()
    # Removed diagnostic metrics and legacy composite flags cannot affect
    # membership or ordering, even when replaced by arbitrary values.
    altered = source.copy()
    model_columns = sorted(set(MODEL_COLUMNS + [c for c in source if 'monotonic' in c]))
    for i, col in enumerate(PANEL_COLUMNS + model_columns):
        if col in altered:
            altered[col] = np.random.default_rng(100 + i).random(len(altered))
    check = apply_policy(altered)
    assert e[['combo', 'tier']].equals(check[['combo', 'tier']])
    derived = {'tier', 'tier_primary_meta', 'technical_pass', 'specificity_pass', 'unmet_nomination_criteria',
               'unmet_A_criteria', 'all_nonbrain_gates_pass', 'availability_status'}
    untouched = [c for c in source.columns if c not in derived]
    pd.testing.assert_frame_equal(source.set_index('combo').sort_index()[[c for c in untouched if c != 'combo']],
        e.set_index('combo').sort_index()[[c for c in untouched if c != 'combo']], check_dtype=False)
    assert hashlib.sha256(source_path.read_bytes()).hexdigest() == before_hash
    e.to_csv(out / 'all_program_evidence.csv', index=False)
    diagnostic = ['combo', 'protein', 'blood_state', 'brain_state'] + [
        c for c in model_columns + PANEL_COLUMNS if c in source]
    source[diagnostic].to_csv(out / 'model_evaluation_only.csv', index=False)
    best = e.drop_duplicates('protein')
    assert set(best.protein) == targets
    best.to_csv(out / 'all_213_targets.csv', index=False)
    priority = e[e.tier.isin(['A', 'B'])].drop_duplicates('protein')
    priority.to_csv(out / 'priority_targets.csv', index=False)
    priority.drop_duplicates('genes').to_csv(out / 'priority_gene_groups.csv', index=False)
    e[e.tier.eq('C')].drop_duplicates('protein').to_csv(out / 'followup_targets_C.csv', index=False)
    e[e.protein.isin(['CD22', 'CD35'])].to_csv(out / 'CD22_CD35_all_programs.csv', index=False)
    brief = {'protein': '抗原', 'genes': '对应基因', 'tier': '当前等级',
        'tier_before_policy_revision': '原程序等级', 'blood_state': '血细胞状态', 'brain_state': '脑细胞状态',
        'protein_rho': '蛋白预测rho', 'confirmation_positive_fraction': '确认种子同向比例',
        'q_gene_decoy': '程序随机基因集合q', 'brain_meta_q': '脑疾病程序q',
        'above_isotype_controls': '超过同型对照', 'independent_p': '独立队列p', 'independent_q': '独立队列q'}
    priority[list(brief)].rename(columns=brief).to_csv(out / '当前候选汇总.csv', index=False, encoding='utf-8-sig')
    changes = e[e.tier.ne(e.tier_before_policy_revision)]
    changes[['combo', 'protein', 'tier_before_policy_revision', 'tier', 'unmet_A_criteria']].to_csv(
        out / 'classification_changes.csv', index=False)
    keep = pd.Series(True, index=e.index)
    flow = [dict(step='all_programs', n_programs=len(e), n_antigens=e.protein.nunique())]
    for gate in GATES:
        keep &= e[gate]
        flow.append(dict(step='gene_program_response' if gate == 'specificity_pass' else gate,
                         n_programs=int(keep.sum()), n_antigens=int(e.loc[keep, 'protein'].nunique())))
    pd.DataFrame(flow).to_csv(out / 'nomination_flow.csv', index=False)
    original_rules = json.loads((root / 'ranking_rules.json').read_text())
    rules = dict(original_rules)
    rules.update(policy_id=POLICY_ID,
        technical='unique gene or defined complex identity and nonzero measured dose IQR; no prediction-performance threshold',
        reproducible='Model evaluation only; no seed-direction-fraction or monotonicity threshold in eligibility or order',
        specificity='gene-decoy q<=0.05; no antigen-panel percentile or P/q requirement',
        ranking='tier, gene-decoy q, brain-meta q, combo as deterministic tie-break; no removed model or panel metrics',
        tier_C='identity and measured range eligible, gene-decoy q<=0.05 and historical brain criterion; independent direction incomplete',
        change_scope='Remove protein rho>=0.20, positive protein donor count>=2, confirmation fraction>=0.8, monotonicity>=0.75 and antigen-panel rank from eligibility and ordering.',
        biological_rule_status='Brain/external rules retained for an isolated mechanical amendment; reviewed separately in nomination_policy_review.txt. These labels are not a final biological priority ranking.',
        panel_values='Retained solely as historical diagnostic measurements; not used by this policy.',
        provenance='User-directed rule amendment after reviewing existing results; not a new independent validation.')
    (out / 'ranking_rules.json').write_text(json.dumps(rules, indent=2))
    summary = dict(status='complete', policy_id=POLICY_ID, source_sha256=before_hash,
        n_programs=len(e), n_antigens=len(best), program_tiers=e.tier.value_counts().to_dict(),
        best_antigen_tiers=best.tier.value_counts().to_dict(),
        priority_targets=priority.protein.tolist(), n_changed_program_tiers=len(changes),
        nomination_flow=flow,
        checks=dict(source_unchanged=True, all_other_measurements_unchanged=True,
                    removed_model_and_panel_values_cannot_affect_tier_or_order=True, complete_antigen_coverage=True),
        source_table='joint_analysis/all_program_evidence.csv',
        changes='Model-evaluation thresholds and antigen-panel rank removed from admission and sorting. No retraining or recomputed P/q values.',
        interpretation='Mechanical rule-amendment counts with historical biological gates; not a final first-choice target list.')
    (out / 'completed.json').write_text(json.dumps(summary, indent=2))
    (root / 'active_nomination.json').write_text(json.dumps(dict(active_directory='current_nomination',
        policy_id=POLICY_ID, historical_rule_directory='joint_analysis',
        preferred_term='peripherally druggable targets candidates'), indent=2))
    manifest = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir())
                if p.is_file() and p.name != 'manifest.json'}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--result-dir', required=True)
    main(ap.parse_args().result_dir)
