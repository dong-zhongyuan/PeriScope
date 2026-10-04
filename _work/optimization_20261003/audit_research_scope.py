"""Separate the manuscript's candidate-discovery aim from later hard gates.

This retrospective scope review uses unchanged, full-family corrected results.
It creates an evidence view, not new A/B assignments or independent tests.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path('/Users/dawnmeng/Desktop/董仲元')
D = ROOT / 'delivery/optimization_20261003'
S = D / 'scope_review'
S.mkdir(exist_ok=True)
source = D / 'joint_analysis/all_program_evidence.csv'
e = pd.read_csv(source)
gates = ['technical_pass', 'reproducible_pass', 'specificity_pass',
         'brain_pass', 'independent_direction_pass']
e['existing_strict_tier'] = e.tier
e['measured_reproducible_program_response'] = (
    e.technical_pass & e.reproducible_pass & e.q_gene_decoy.le(.05))
e['same_cell_state_disease_evidence'] = e.cell_state_sensitivity_matched_state_support
e['response_and_cell_state_support'] = (
    e.measured_reproducible_program_response & e.same_cell_state_disease_evidence)
e['same_cell_state_meta_q'] = e[['cell_state_sensitivity_meta_q_up',
                                'cell_state_sensitivity_meta_q_down']].min(axis=1)
e.to_csv(S / 'all_program_evidence_roles.csv', index=False)
supported = e[e.response_and_cell_state_support]
supported.to_csv(S / 'cell_state_supported_programs.csv', index=False)
cd35 = supported[genetic.protein.eq('CD35')]
brief = {
    'protein': '抗原', 'genes': '对应基因', 'blood_state': '血细胞状态',
    'brain_state': '脑细胞状态', 'n_genes': '本轮程序基因数',
    'response_direction': '蛋白增加时的程序响应方向',
    'protein_rho': '留出供体蛋白预测rho',
    'confirmation_positive_fraction': '确认种子同向比例',
    'confirmation_monotonic_fraction': '确认剂量单调比例',
    'q_gene_decoy': '随机基因对照q',
    'same_cell_state_meta_q': '同类细胞跨队列合并q',
    'above_isotype_controls': '超过同型对照',
    'protein_percentile': '其他抗原参照百分位',
    'ma_donor_score_effect': '混合空间组织PD减对照',
    'cell_state_sensitivity_kamath_effect': 'Kamath同类细胞PD减对照',
    'cell_state_sensitivity_gse157783_effect': 'GSE157783同类细胞PD减对照',
    'existing_strict_tier': '原严格筛法等级',
}
cd35[list(brief)].rename(columns=brief).to_csv(S / 'CD35同类细胞证据.csv', index=False, encoding='utf-8-sig')
leave_out = []
for omitted in gates:
    mask = e[[g for g in gates if g != omitted]].all(axis=1)
    leave_out.append(dict(omitted_gate=omitted, n_programs=int(mask.sum()),
                          antigens=sorted(e.loc[mask, 'protein'].unique().tolist())))
counts = {}
for name, mask in [
    ('response_and_cell_state_support', e.response_and_cell_state_support)]:
    counts[name] = dict(n_programs=int(mask.sum()), n_antigens=int(e.loc[mask, 'protein'].nunique()))
summary = dict(
    status='complete', scope='Retrospective audit of evidence roles against the manuscript candidate-discovery aim.',
    unchanged='No retraining, new P values, altered full-family q values, or reassignment of existing A/B tiers.',
    technical_and_repeatability_requirements='Retained here to isolate the effect of evidence-role changes.',
    brain_comparison='Matched fine brain cell states in locked Kamath and independent GSE157783, using already computed pooled-family tests.',
    spatial_role='Localization and tissue-context evidence; its bulk disease sign is recorded separately.',
    protein_panel_role='Relative preference among measured antigens, retained as a comparative dimension rather than a universal candidate veto.',
    program_generation='Current retrained programs; not restoration of the old manuscript gene sets or old statistics.',
    source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    source_manuscript_sha256=hashlib.sha256((ROOT / 'delivery/manuscript_gm_v2.tex').read_bytes()).hexdigest(),
    counts=counts, one_gate_leave_out=leave_out,
    CD35_programs=cd35.combo.tolist(),
    interpretation='An evidence view for refining the research scope, not a newly preregistered validation or proof of intervention efficacy.')
(S / 'scope_audit.json').write_text(json.dumps(summary, indent=2))
assert (e.tier == pd.read_csv(source).tier).all()
print(json.dumps(dict(counts=counts, CD35_programs=summary['CD35_programs']), indent=2))
