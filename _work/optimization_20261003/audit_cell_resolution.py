"""Audit resolution and prior exclusions without changing frozen statistics.

Run after audit_cell_resolution_inputs.py on the server and copy its outputs
into result_dir/cell_resolution_audit_20261003. Uses only pandas/numpy.
"""
from pathlib import Path
import argparse, hashlib, json
import numpy as np
import pandas as pd

def run(root):
    root=Path(root);out=root/'cell_resolution_audit_20261003';out.mkdir(exist_ok=True)
    source=root/'joint_analysis/all_program_evidence.csv'
    current=root/'six_criteria_20261003/all_6849_programs_six_criteria.csv'
    hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,current]}
    e=pd.read_csv(current); assert len(e)==6849 and e.protein.nunique()==213
    p='cell_state_sensitivity_'
    e['audit_spatial_opposes_matched_brain']=e.ma_donor_score_effect.mul(e[p+'kamath_effect']).lt(0)
    e['audit_matched_brain_support']=e[p+'matched_state_support']
    e['audit_old_spatial_sign_rule_pass']=e.ma_donor_score_effect.mul(e.kamath_donor_score_effect).gt(0)
    e['audit_only_C1_blocks_six']=e.C1.ne('pass')&e[['C2','C3','C4','C5','C6']].eq('pass').all(axis=1)
    e['audit_nomination_changed']=False
    e['audit_source_resolution']=e.blood_state.map({'cDC':'cDC1+cDC2+ASDC mapped together; disease atlases marker-labelled cDC',
        'pDC':'pDC state','classical_mono':'classical monocyte state','nonclassical_mono':'nonclassical monocyte state'})
    e['audit_recipient_resolution']=e.brain_state.map({'astro':'astrocytes, multiple clusters combined',
        'microglia_mhc2':'marker-assigned MHC-II microglia state','microglia_homeostatic':'marker-assigned homeostatic microglia state'})
    e['audit_followup']=np.where(e.C4.eq('unresolved'),'coverage or detectable-gene limitation; not biological failure',
        np.where(e.audit_only_C1_blocks_six,'review source-context decision at measured subtype and activation-state level',
        np.where(e.audit_matched_brain_support&e.audit_spatial_opposes_matched_brain,'matched-cell support retained; tissue sign cannot veto','retain measured evidence and current status')))
    cols=['combo','protein','blood_state','brain_state','historical_strict_tier','C1','C1_note','C2','C3','C4','C5','C6',
          'six_criteria_pass','q_gene_decoy','matched_brain_q','ma_donor_score_effect',p+'kamath_effect',
          p+'gse157783_effect',p+'kamath_status',p+'gse157783_status']+[c for c in e if c.startswith('audit_')]
    e[cols].to_csv(out/'全部6849程序分辨率与排除复查.csv',index=False,encoding='utf-8-sig')
    subsets={
        '已通过六项但空间异向':e.six_criteria_pass&e.audit_spatial_opposes_matched_brain,
        '同脑状态支持但空间异向':e.audit_matched_brain_support&e.audit_spatial_opposes_matched_brain,
        '仅被来源语境C1暂缓':e.audit_only_C1_blocks_six,
        '验证覆盖不足或基因不足':e.C4.eq('unresolved'),
        'CD22_CD35逐程序复查':e.protein.isin(['CD22','CD35'])}
    for name,mask in subsets.items():e.loc[mask,cols].to_csv(out/(name+'.csv'),index=False,encoding='utf-8-sig')
    curation=json.loads((root/'six_criteria_20261003/six_criteria_curation.json').read_text())
    categories=[]
    for i,item in enumerate(curation['source_context_review'],1):
        mask=e.protein.isin(item['antigens'])
        if 'states' in item:mask &= e.blood_state.isin(item['states'])
        hit=e[mask];only=hit[hit.audit_only_C1_blocks_six]
        categories.append(dict(rule=i,antigens=';'.join(item['antigens']),states=';'.join(item.get('states',['all four model states'])),
            programs=len(hit),antigen_count=hit.protein.nunique(),only_C1_programs=len(only),only_C1_antigens=only.protein.nunique(),
            original_reason=item['reason'],source=item['source'],audit='Curated source-context concern, not a computed proof of absence. Revisit subtype, activation, clone, and cell association separately.'))
    pd.DataFrame(categories).to_csv(out/'全部来源语境规则复查.csv',index=False,encoding='utf-8-sig')
    detail=pd.read_csv(root/'assay_identity_by_heldout_subtype.csv')
    statepairs=e.loc[e.audit_only_C1_blocks_six,['protein','blood_state']].drop_duplicates().rename(columns={'protein':'target'})
    detail.merge(statepairs,on=['target','blood_state'],validate='many_to_one').to_csv(out/'C1暂缓候选逐供者亚型实测.csv',index=False,encoding='utf-8-sig')
    detail[detail.target.isin(['CD22','CD35','CD152'])].to_csv(out/'CD22_CD35_CD152逐供者亚型.csv',index=False,encoding='utf-8-sig')
    short=detail[detail.target.eq('CD22') & detail.citeseq_subtype.eq('cDC2')]
    assert set(short.donor)=={'P6','P7','P8'}
    coverage=pd.read_csv(out/'cell_state_donor_coverage.csv')
    c=coverage[coverage.meets_25_cells].copy()
    c.groupby(['cohort','state','split','condition']).size().rename('n_evaluable_donors').reset_index().to_csv(out/'可评估供者按状态汇总.csv',index=False,encoding='utf-8-sig')
    cite=pd.read_csv(out/'cite_subtype_to_model_state.csv')
    cite['sampled_by_current_B_sampler']=cite.split.eq('train')&cite.n_cells.ge(25)
    cite['included_in_dose_pool']=cite.split.eq('train')&cite.model_state.ne('none')
    cite['included_in_heldout_context_table']=cite.split.eq('test')&cite.n_cells.ge(25)&cite.model_state.ne('none')
    cite.to_csv(out/'CITE亚型抽样与剂量覆盖.csv',index=False,encoding='utf-8-sig')
    units=[
      ('01 输入配对','Hao逐细胞RNA/ADT/metadata条形码交集','paired single cells','配对审计通过；细胞身份继承原l2注释','preprocess_citeseq_hao.py'),
      ('02 疾病图谱注释','各队列谱系内聚类后统一state_pure标签','cluster-labelled single cells','共同标签不等于每个队列组成完全相同；现存映射和聚类评分可核查，生成聚类脚本本次未定位','purification/*.csv; purification_report.json'),
      ('03 特征选择与中心化','训练供者内选择特征，共用基因词表','shared features, not a cell-type result','共用词表不混合下游细胞状态的PD对照检验','input_utils.py; build_training_data.py'),
      ('04 蛋白桥接训练','全图谱共享编码器/蛋白头；供者×l2亚型平衡抽样','shared model with subtype-balanced batches','不是每亚型独立训练；每组25细胞使ASDC全未入桥接抽样，cDC1仅P4进入训练','train_ccwm.py; training_utils.py'),
      ('05 耦合训练','疾病×血状态×脑状态分层，组内分布匹配','cell-state-conditioned distributions','4血状态×3有训练支持脑状态；并非一一真实配对的血细胞与脑细胞','train_ccwm.py'),
      ('06 剂量定义','训练CITE供者内按4个模型血状态取分位数','model-state dose pools','cDC合并cDC1/cDC2/ASDC；文档写subtype-specific过细；稀有亚型剂量池与训练抽样覆盖不同','candidate_registry.py'),
      ('07 扰动与程序','每条血脑轴内配对源细胞和噪声；每个剂量输出细胞均值','axis-specific response programs','raw与relative两种定义都保留轴；CLR补偿不改变细胞分层；不是对空间混合组织定义程序','screen_candidates.py; dose_programs.py; target_specific/'),
      ('08 随机基因/种子评价','固定轴与固定程序，确认种子响应','program-level statistics within axis','没有将脑状态混为一个测试；同一抗原跨轴不能拼证据','randko.py; run_joint_analysis.py'),
      ('09 血端疾病关联','供者×血状态内RNA与推断ADT','within-state donor summaries','细胞内均值与细胞数量比例是不同结果；ADT推断不同于PD队列实测ADT','disease_association.py; blood_covariate_sensitivity.py'),
      ('10 留出脑疾病评价','Kamath黑质、供者×指定脑状态','within-state donor summaries','astro/MHC-II各3PD+2Control；稳态仅1Control达到25细胞；供者聚合不丢失状态分辨率','brain_competitive.py; program_effect_sizes.py'),
      ('11 外部脑评价','GSE157783供者×同名state_pure','within-state donor summaries','astro/MHC-II各5PD+6Control；无homeostatic标签；有DAM但训练Kamath无DAM，不能静默替换','cell_state_validation.py; independent_validation.py'),
      ('12 空间验证','11859混合点位、10供者','mixed spatial spots and donor tissue means','无obsm/uns解卷积结果；marker高分点位仍混合；定位证据保留，组织方向不能否决同细胞支持','spatial_analysis.py; brain_competitive.py'),
      ('13 组织疾病轴','血脑各状态供者；Ma整组织；PDD宽astro','explicitly mixed resolution comparisons','Ma与PDD敏感性不能称所有均最细状态独立验证；该步不是候选门槛','tissue_axis.py'),
      ('15 可接近性和药物','蛋白/复合物/表位及已知药物机制','molecular annotation','不是细胞类型测量；蛋白在细胞膜不自动证明该来源细胞的表面表达','screen_six_criteria.py'),
      ('16 来源语境C1','实际ADT+人工典型谱系排查','measurement plus manual biological curation','典型谱系非排他；41抗原因人工语境被全来源暂缓，另2个身份未解决；CD152已有反例；不把待核实说成阴性','six_criteria_curation.json'),
      ('17 跨层合并及论文','每条完整轴评价；最终按抗原汇总','program first, antigen summary second','最新版六项C4已不用Ma；论文仍为旧24程序/1469特征版，未自动变成本轮6849程序/4741特征版','screen_six_criteria.py; delivery/manuscript_gm_v2.tex'),
    ]
    pd.DataFrame(units,columns=['步骤','真实分析单位','分辨率','审查结论','代码与输入']).to_csv(out/'全流程分辨率审查.csv',index=False,encoding='utf-8-sig')
    count=lambda m:dict(programs=int(m.sum()),antigens=int(e.loc[m,'protein'].nunique()))
    s=dict(source_hashes=hashes,total=count(pd.Series(True,index=e.index)),
        subsets={name:count(mask) for name,mask in subsets.items()},
        C4_unresolved_causes=e[e.C4.eq('unresolved')].groupby(['brain_state',p+'kamath_status',p+'gse157783_status']).size().rename('programs').reset_index().to_dict('records'),
        existing_A6=count(e.six_criteria_pass),
        C1_only_additional_antigens=sorted(set(e.loc[e.audit_only_C1_blocks_six,'protein'])-set(e.loc[e.six_criteria_pass,'protein'])),
        counts_are_not_new_nominations=True,all_statistics_unchanged=True,
        original_purification_generator='not located in current project; retained mapping files and cluster-score report inspected',
        CD22_cDC2_measured=short[['donor','n_cells','RNA_fraction_any_detected','ADT_mapped_RNA_partial_rho_isotype']].to_dict('records'))
    (out/'summary.json').write_text(json.dumps(s,ensure_ascii=False,indent=2))
    for rel,h in hashes.items():assert hashlib.sha256((root/rel).read_bytes()).hexdigest()==h
    (out/'manifest.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir()) if p.is_file() and p.name!='manifest.json'},indent=2))
    print(json.dumps(s,ensure_ascii=False,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--result-dir',required=True);args=ap.parse_args();run(args.result_dir)
