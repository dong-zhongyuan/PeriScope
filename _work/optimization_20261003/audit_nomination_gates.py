"""Audit existing nomination gates without changing scores, thresholds or tiers."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
D = ROOT / 'delivery/optimization_20261003'
SOURCE = D / 'joint_analysis/all_program_evidence.csv'
OUT = D / 'gate_audit'
OUT.mkdir(exist_ok=True)
e = pd.read_csv(SOURCE)
before_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
sign = np.where(e.brain_meta_q_up < e.brain_meta_q_down, 1., -1.)

# Ordered only to decompose the current composite gates. The code evaluates the
# component flags jointly; first-failure counts depend on this stated order.
specs = [
 ('identity', '抗原唯一基因映射或已定义复合物', e.mapping_pass,
  '保证干预身份、基因注释及证据追踪一致', '保留身份核实；多基因表位不应因不能唯一映射而自动否定其抗原候选资格'),
 ('protein_rho', '留出供体蛋白预测平均rho≥0.20', e.protein_rho.ge(.2),
  '评估RNA到目标蛋白输入的可预测性', '输入质量是核心；0.20不是经校准的生物学界线，应补供体不确定性与简单预测基线'),
 ('positive_donors', '至少2位留出供体rho>0', e.positive_protein_donors.ge(2),
  '避免平均相关由单个供体支撑', '保留供体层面质量检查；本次在rho门槛后无额外淘汰'),
 ('dose_range', '来源细胞实测剂量IQR>1e-6', e.measured_dose_pass,
  '保证存在非退化的来源细胞测量范围', '保留；非零IQR本身不等于特异表达或药理有效剂量'),
 ('seed_direction', '确认种子至少4/5与发现方向一致', e.confirmation_positive_fraction.ge(.8),
  '独立训练重复中保持冻结程序的响应方向', '保留重复性证据；4/5是工作阈值，不能替代效应量、区间或供体外验证'),
 ('monotonic', '确认种子相邻剂量步同向比例≥0.75', e.confirmation_monotonic_fraction.ge(.75),
  '选择在采样范围内近似单调的剂量响应', '适合作曲线类型描述及稳健性检查，不宜否决所有非单调但可复现响应'),
 ('gene_decoy', '随机基因集合q≤0.05', e.q_gene_decoy.le(.05),
  '冻结程序的确认响应超过同规模随机基因集合', '保留程序响应集中性证据与全家族校正；这不是实验KO，也不等于干预有效'),
 ('protein_percentile', '相对其余抗原排名≥95百分位', e.protein_percentile.ge(.95),
  '衡量该程序对目标抗原的相对偏好', '改作比较和排序维度；其他生物抗原不是全为无效干预，前5%不宜作为统一否决条件'),
 ('brain_meta', 'Ma/Kamath竞争性脑程序meta q≤0.05', e.brain_meta_q.le(.05),
  '评估程序与PD脑疾病变化的联系', '疾病支持重要；优先采用同细胞状态的比较，缺失队列须标为未充分评估'),
 ('ma_relative', 'Ma相对背景方向与合并方向一致', pd.Series(sign * e.ma_relative_effect > 0, index=e.index),
  '确认空间组织竞争性检验与合并方向一致', '保留效应及异质性描述；混合组织不是细胞内效应的同义测量'),
 ('kamath_relative', 'Kamath相对背景方向与合并方向一致', pd.Series(sign * e.kamath_relative_effect > 0, index=e.index),
  '确认脑细胞竞争性检验方向', '同细胞状态方向有意义；本次该步骤淘汰的5个程序均为供体不足，不能判为相反生物学'),
 ('mixed_actual_direction', 'Ma混合组织与Kamath细胞实际疾病方向一致', e.ma_donor_score_effect.mul(e.kamath_donor_score_effect).gt(0),
  '试图要求跨测量层级的疾病改变同向', '撤出普遍否决条件；先作细胞组成和组织语境分析，细胞内疾病支持由匹配细胞队列评价'),
 ('external_direction', 'GSE157783与Ma/Kamath共同方向一致', e.independent_direction_pass,
  '检验额外脑队列的程序方向', '独立支持重要；当前仅同号而不要求显著，且被上一方向条件连带置零，应改为独立的匹配细胞证据'),
]
checks = pd.DataFrame({k: m.astype(bool) for k, _, m, _, _ in specs}, index=e.index)
assert len(e) == 6849 and e.protein.nunique() == 213
assert (checks.iloc[:, :4].all(axis=1) == e.technical_pass).all()
assert (checks.iloc[:, 4:6].all(axis=1) == e.reproducible_pass).all()
assert (checks.iloc[:, 6:8].all(axis=1) == e.specificity_pass).all()
assert (checks.iloc[:, 8:12].all(axis=1) == e.brain_pass).all()
assert (checks.all(axis=1) == e.tier.eq('A')).all()

keep = pd.Series(True, index=e.index)
furthest = pd.Series(0, index=e.index)
first_failure = pd.Series('', index=e.index)
flow = []
for i, (key, label, mask, meaning, recommendation) in enumerate(specs, start=1):
    old = keep.copy()
    previous_antigens = set(e.loc[old, 'protein'])
    keep &= mask
    remaining_antigens = set(e.loc[keep, 'protein'])
    lost = previous_antigens - remaining_antigens
    newly_failed = old & ~mask
    first_failure.loc[newly_failed] = key
    furthest.loc[keep] = i
    unavailable = newly_failed & (e.ma_test_status.ne('ok') | e.kamath_test_status.ne('ok')) if 9 <= i <= 12 else pd.Series(False, index=e.index)
    flow.append(dict(步骤=i, 条件=label, 条件标识=key,
        本步进入程序数=int(old.sum()), 本步淘汰程序数=int(newly_failed.sum()), 累积保留程序数=int(keep.sum()),
        本步进入抗原数=len(previous_antigens), 本步淘汰抗原数=len(lost), 累积保留抗原数=len(remaining_antigens),
        本步失去全部合格程序的抗原=';'.join(sorted(lost)),
        单独满足本条件程序数=int(mask.sum()), 单独至少一程序满足的抗原数=int(e.loc[mask, 'protein'].nunique()),
        本步淘汰且脑数据不完整程序数=int(unavailable.sum()),
        生物学或统计作用=meaning, 审查建议=recommendation))
flow = pd.DataFrame(flow)
flow.to_csv(OUT / '逐步通过与淘汰.csv', index=False, encoding='utf-8-sig')

per_program = e[['combo', 'protein', 'blood_state', 'brain_state', 'program_definition', 'tier',
    'ma_test_status', 'kamath_test_status', 'unmet_A_criteria']].copy()
per_program['first_failure_in_stated_order'] = first_failure
per_program['furthest_completed_step'] = furthest
per_program = pd.concat([per_program, checks.add_prefix('pass_')], axis=1)
per_program.to_csv(OUT / '全部程序逐项核查.csv', index=False)

antigen = []
for target, rows in per_program.groupby('protein', sort=True):
    stage = int(rows.furthest_completed_step.max())
    best = rows[rows.furthest_completed_step.eq(stage)]
    antigen.append(dict(抗原=target, 程序数=len(rows), 最远通过步骤=stage,
        最后合格程序首次失败条件=specs[stage][1] if stage < len(specs) else '全部通过',
        对应程序=';'.join(best.combo), 现严格最好等级=sorted(rows.tier)[0]))
pd.DataFrame(antigen).to_csv(OUT / '213个抗原去向.csv', index=False, encoding='utf-8-sig')

leave = []
for key, label, _, _, _ in specs:
    mask = checks.drop(columns=key).all(axis=1)
    leave.append(dict(仅去掉条件=label, 其余条件全通过程序数=int(mask.sum()),
                      抗原数=int(e.loc[mask, 'protein'].nunique()), 抗原=';'.join(sorted(e.loc[mask, 'protein'].unique()))))
pd.DataFrame(leave).to_csv(OUT / '单项门槛影响.csv', index=False, encoding='utf-8-sig')

prioritized = e[e.tier.eq('B') | e.combo.isin([
    'raw::CD35__nonclassical_mono__x__astro__down', 'raw::CD35__classical_mono__x__astro__down'])].copy()
cols = ['combo', 'protein', 'tier', 'protein_rho', 'positive_protein_donors', 'dose_iqr',
    'confirmation_positive_fraction', 'confirmation_monotonic_fraction', 'q_gene_decoy',
    'protein_percentile', 'q_protein_decoy', 'above_isotype_controls', 'brain_meta_q',
    'ma_relative_effect', 'kamath_relative_effect', 'ma_donor_score_effect', 'kamath_donor_score_effect',
    'independent_effect', 'independent_p', 'independent_q',
    'cell_state_sensitivity_meta_q_up', 'cell_state_sensitivity_meta_q_down',
    'cell_state_sensitivity_matched_state_support', 'unmet_A_criteria']
prioritized[cols].to_csv(OUT / 'B类与CD35证据明细.csv', index=False, encoding='utf-8-sig')

reg = json.loads((D / 'candidate_registry.json').read_text())
counts = dict(original_ADT_channels=sum(len(x['channels']) for x in reg),
    grouped_antigens=len(reg), isotype_groups=sum(x['is_control'] for x in reg),
    biological_antigens=sum(not x['is_control'] for x in reg), blood_brain_axes=int(e.axis.nunique()))
core = e.technical_pass & e.reproducible_pass & e.q_gene_decoy.le(.05)
matched = core & e.cell_state_sensitivity_matched_state_support
summary = dict(status='complete', preferred_term='peripherally druggable targets candidates',
    source=str(SOURCE.relative_to(ROOT)), source_sha256=before_hash, registry=counts,
    stages_are_an_ordered_decomposition=True,
    antigen_removal='An antigen is removed at a step only if none of its complete program paths remains eligible.',
    program_tiers=e.tier.value_counts().to_dict(),
    best_antigen_tiers=e.sort_values('tier').drop_duplicates('protein').tier.value_counts().to_dict(),
    n_strict_A_programs=int(e.tier.eq('A').sum()),
    pre_program_slots_per_definition=counts['biological_antigens']*counts['blood_brain_axes']*2,
    program_definitions=e.program_definition.value_counts().to_dict(),
    response_and_matched_state=dict(programs=int(matched.sum()), antigens=int(e.loc[matched, 'protein'].nunique())),
    prebrain_insufficient_programs=int((e.technical_pass & e.reproducible_pass & e.specificity_pass & e.kamath_test_status.ne('ok')).sum()),
    entire_homeostatic_state_unassessed=bool(e.loc[e.brain_state.eq('microglia_homeostatic'), 'kamath_test_status'].eq('insufficient_donors').all()),
    classification_changed=False, statistics_changed=False,
    branch_note='final_results retains the 1840-program relative branch; joint_analysis is the 6849-program pooled evidence table used here.')
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == before_hash
(OUT / 'audit_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2))
(OUT / 'audit_manifest.json').write_text(json.dumps({p.name: hashlib.sha256(p.read_bytes()).hexdigest()
    for p in sorted(OUT.iterdir()) if p.is_file() and p.name != 'audit_manifest.json'}, ensure_ascii=False, indent=2))
print(json.dumps(summary, ensure_ascii=False, indent=2))
