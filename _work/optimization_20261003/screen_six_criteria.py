"""Six biological/drug criteria on frozen full-family program evidence.

No model-performance metrics enter eligibility or ordering. Drug-target
mechanism evidence is kept at record level, with accession and source links.
"""
from pathlib import Path
import argparse, hashlib, json, shutil, subprocess, sys
import numpy as np
import pandas as pd

POLICY='20261003_six_biological_and_drug_criteria'

def main(result_dir):
    root=Path(result_dir); out=root/'six_criteria_20261003'; out.mkdir(exist_ok=True)
    load=lambda n:json.loads((out/n).read_text())
    registry=json.loads((root/'candidate_registry.json').read_text())
    registry=[x for x in registry if not x['is_control']]
    uni=load('uniprot_targets.json')
    proteins={}
    for u in uni:
        for g in u.get('genes',[]):
            if 'geneName' in g: proteins.setdefault(g['geneName']['value'],[]).append(u)
    targets={x['target_chembl_id']:x for x in load('chembl_targets.json')}
    molecules={x['molecule_chembl_id']:x for x in load('chembl_molecules.json')}
    mechanisms=load('chembl_mechanisms.json')
    curation=json.loads(Path(__file__).with_name('six_criteria_curation.json').read_text())
    rows=[]
    for a in registry:
        gs=set(a['genes'])
        for r in mechanisms:
            t=targets[r['target_chembl_id']]
            components={s['component_synonym'] for c in t['target_components'] for s in c['target_component_synonyms'] if s['syn_type']=='GENE_SYMBOL'}
            if not gs.intersection(components):continue
            mol=molecules[r['molecule_chembl_id']]
            phase=float(r.get('max_phase') or 0)
            phase=max(phase,float(mol.get('max_phase') or 0))
            exact=bool(components) and components.issubset(gs) and t['target_type'] in ['SINGLE PROTEIN','PROTEIN COMPLEX']
            resolved=curation['binding_site_gene_resolution'].get(str(r['mec_id']))
            if resolved:
                exact=bool(gs.intersection(resolved['binding_genes']))
            status='verified_direct_clinical'
            note='Direct clinical molecular mechanism in ChEMBL; not a drug response prediction.'
            if r.get('direct_interaction')!=1 or r.get('molecular_mechanism')!=1:
                status='not_direct_molecular_mechanism'
            elif not exact:status='family_or_complex_member_requires_resolution'
            elif phase<.5:status='clinical_stage_not_established'
            if a['target'] in ['CD45RA','CD45RB','CD45RO']:
                status='pan_PTPRC_drug_isoform_specificity_unresolved'
            if a['target']=='CD35':
                status='CR1_functional_analogue_not_direct_CR1_targeting'
                note='Mirococept is CR1-derived and inhibits complement convertases; not an inhibitor binding endogenous CR1.'
            if a['target']=='CD16':
                status='Fc_receptor_engagement_not_independent_CD16_intervention'
                note='Fc recruitment by antibodies does not isolate the modelled CD16 abundance intervention.'
            if a['target']=='CD44' and 'BIVATUZUMAB' in (mol.get('pref_name') or ''):
                status='CD44v6_epitope_specificity_unresolved'
                note='Drug targets CD44v6, whereas panel CD44 is not established as v6-specific.'
            if (mol.get('pref_name') or '')=='HUMAX-IL15':
                status='incorrect_source_target_annotation'
                note='HuMax-IL15/AMG714 targets IL15; source CD38 annotation is not accepted.'
            if (mol.get('pref_name') or '')=='NANGIBOTIDE':
                status='ligand_trap_not_direct_TREM1_binding'
                note='Nangibotide binds a TREM1 agonist ligand; retain as a pathway-modulation route, not direct receptor binding.'
            if a['target']=='HLA-DR' and (mol.get('pref_name') or '')=='LYM-1':
                status='HLA_DR_allele_epitope_requires_resolution'
                note='LYM-1 recognizes an HLA-DR10-associated epitope; generic HLA-DR panel signal does not resolve it.'
            refs=';'.join(x['ref_url'] for x in r.get('mechanism_refs',[]) if x.get('ref_url'))
            rows.append(dict(protein=a['target'],genes=';'.join(a['genes']),drug=mol.get('pref_name') or r['molecule_chembl_id'],
                molecule_chembl_id=r['molecule_chembl_id'],target_chembl_id=r['target_chembl_id'],
                molecular_target=t['pref_name'],target_genes=';'.join(sorted(components)),target_type=t['target_type'],
                action=r['action_type'],mechanism=r['mechanism_of_action'],clinical_max_phase=phase,
                binding_site_comment=r.get('binding_site_comment'),site_id=r.get('site_id'),
                binding_genes_resolved=';'.join(resolved['binding_genes']) if resolved else '',
                ever_approved=phase==4,withdrawn_flag=bool(mol.get('withdrawn_flag')),
                molecule_type=mol.get('molecule_type'),direct_interaction=r.get('direct_interaction'),
                match_status=status,interpretation=note,mechanism_refs=refs,
                source_url='https://www.ebi.ac.uk/chembl/explore/compound/'+r['molecule_chembl_id']))
    drugs=pd.DataFrame(rows).drop_duplicates(['protein','molecule_chembl_id','target_chembl_id','action'])
    drugs.to_csv(out/'all_drug_target_mechanisms.csv',index=False)
    drugs[drugs.match_status.isin(['CR1_functional_analogue_not_direct_CR1_targeting',
          'ligand_trap_not_direct_TREM1_binding'])].to_csv(out/'related_pharmacology_routes.csv',index=False)
    verified=drugs[drugs.match_status.eq('verified_direct_clinical')].copy()
    # Ever-approved and historical clinical-stage strata, not a claim that every
    # agent remains marketed or is in an ongoing trial at retrieval time.
    verified=verified.sort_values(['protein','withdrawn_flag','clinical_max_phase','drug'],ascending=[True,True,False,True])
    verified.to_csv(out/'verified_drug_target_mechanisms.csv',index=False)
    annotation=[]
    for a in registry:
        records=[u for g in a['genes'] for u in proteins.get(g,[])]
        loc=[]; tissue=[]; extra=[]
        for u in records:
            for c in u.get('comments',[]):
                if c['commentType']=='SUBCELLULAR LOCATION':
                    loc += [s['location']['value'] for s in c.get('subcellularLocations',[]) if 'location' in s]
                if c['commentType']=='TISSUE SPECIFICITY':
                    tissue += [t['value'] for t in c.get('texts',[])]
            extra += [f.get('description','') for f in u.get('features',[]) if f.get('type')=='Topological domain']
        accessible=any(x in ['Cell membrane','Secreted','Cell surface'] for x in loc) or 'Extracellular' in extra or a['target'] in ['CD15','CD57']
        hits=verified[verified.protein.eq(a['target'])]
        approved=hits[hits.ever_approved & ~hits.withdrawn_flag]
        annotation.append(dict(protein=a['target'],identity_mapping=a['mapping'],
            identity_interpretable=bool(a['genes']) or a['target'] in ['CD15','CD57'],
            surface_accessible=accessible,subcellular_location=';'.join(sorted(set(loc))),
            tissue_context=' | '.join(dict.fromkeys(tissue)),
            uniprot_accessions=';'.join(sorted({u['primaryAccession'] for u in records})),
            direct_clinical_drug_count=int(hits.molecule_chembl_id.nunique()),
            approved_not_flagged_withdrawn_count=int(approved.molecule_chembl_id.nunique()),
            drug_examples=';'.join(hits.drug.drop_duplicates().head(6)),
            approved_examples=';'.join(approved.drug.drop_duplicates().head(6))))
    ann=pd.DataFrame(annotation);ann.to_csv(out/'all_antigen_annotations.csv',index=False)
    source=root/'joint_analysis/all_program_evidence.csv'
    h=hashlib.sha256(source.read_bytes()).hexdigest()
    e=pd.read_csv(source)
    e=e.rename(columns={'tier':'historical_strict_tier'})
    e=e.merge(ann,on='protein',validate='many_to_one')
    original=e.copy()
    # All numeric support remains visible; no rho, donor-count, lineage RNA
    # fraction or monotonicity cutoff is introduced.
    e['C1']=np.where(e.identity_interpretable & e.heldout_identity_n_donors.gt(0)
                    & e.heldout_observed_fraction_ADT_counts_positive.gt(0),'pass','unresolved')
    e['C1_note']='Resolved assay object and directly measured ADT in matched source-cell context; RNA/ADT context retained.'
    for item in curation['source_context_review']:
        mask=e.protein.isin(item['antigens'])
        if item.get('states'):mask &= e.blood_state.isin(item['states'])
        e.loc[mask,'C1']='unresolved'
        e.loc[mask,'C1_note']=item['reason']
    e['C2']=np.where(e.dose_iqr.gt(1e-6),'pass','not_supported')
    e['C3']=np.where(e.q_gene_decoy.le(.05),'pass','not_supported')
    prefix='cell_state_sensitivity_'
    e['matched_brain_q']=e[[prefix+'meta_q_up',prefix+'meta_q_down']].min(axis=1)
    e['C4']=np.where(e[prefix+'matched_state_support'],'pass',np.where(e[prefix+'both_available'],'not_supported','unresolved'))
    e['C5']=np.where(e.surface_accessible,'pass','unresolved')
    e['C5_note']='Extracellular/cell-surface route supported by UniProt; specific clinical route is recorded under C6.'
    e['C6']=np.where(e.direct_clinical_drug_count.gt(0),'pass','not_established')
    e['C6_note']='Direct molecular target with documented clinical development; stages and withdrawal flags kept separately.'
    mask=e.protein.eq('CD35')
    e.loc[mask,'C6_note']='CR1-derived Mirococept supports a functional-analogue route, not direct targeting of endogenous CD35.'
    for antigen in ['CD45RA','CD45RB','CD45RO','CD44','CD16']:
        e.loc[e.protein.eq(antigen),'C6_note']='Potential drug exists but antigen/subunit/isoform or intervention correspondence requires resolution; see mechanism table.'
    e.loc[e.protein.eq('CD354'),'C6_note']='Nangibotide is a clinical TREM1-pathway ligand-trap route; direct TREM1 receptor binding is not established.'
    gates=['C1','C2','C3','C4','C5','C6']
    e['six_criteria_pass']=e[gates].eq('pass').all(axis=1)
    e['criteria_pass_prefix_length']=e[gates].eq('pass').astype(int).cumprod(axis=1).sum(axis=1)
    e['five_criteria_pass']=e[gates[:5]].eq('pass').all(axis=1)
    e['six_criteria_approved_subset']=e.six_criteria_pass & e.approved_not_flagged_withdrawn_count.gt(0)
    e['status']=np.where(e.six_criteria_pass,'A6_drug_supported',np.where(e.five_criteria_pass,'five_criteria_without_verified_direct_clinical_drug',
                   np.where(e[gates].eq('unresolved').any(axis=1),'unresolved_evidence','tested_not_supported')))
    e['unmet_six_criteria']=[';'.join(c+':'+r[c] for c in gates if r[c]!='pass') for r in e.to_dict('records')]
    # This scalar program-mean hypothesis uses a matched-cell-state effect.
    # The existing predicted_restorative_protein_change also uses matched
    # cells, but is based on a multigene projection and seed agreement.
    # It is a different estimand, not a mixed-tissue predecessor.
    e['matched_brain_disease_direction']=np.where(e[prefix+'kamath_effect']>0,'up','down')
    mean_effect=e.confirmation_high_minus_low_program_mean
    e['matched_program_abundance_hypothesis']=np.where(mean_effect*e[prefix+'kamath_effect']>0,'decrease',
        np.where(mean_effect*e[prefix+'kamath_effect']<0,'increase','unresolved'))
    e['drug_response_predicted']=False
    e=e.sort_values(['six_criteria_pass','six_criteria_approved_subset','criteria_pass_prefix_length','q_gene_decoy','matched_brain_q','protein','combo'],ascending=[False,False,False,True,True,True,True])
    e.to_csv(out/'all_6849_programs_six_criteria.csv',index=False)
    best=e.drop_duplicates('protein');assert len(best)==213
    best.to_csv(out/'all_213_antigen_dispositions.csv',index=False)
    nominated=e[e.six_criteria_pass]
    nominated.to_csv(out/'A6_all_programs.csv',index=False)
    nominated.drop_duplicates('protein').to_csv(out/'A6_antigens.csv',index=False)
    nominated[nominated.six_criteria_approved_subset].drop_duplicates('protein').to_csv(out/'A6_approved_drug_subset.csv',index=False)
    e[e.protein.isin(['CD22','CD35'])].to_csv(out/'CD22_CD35_six_criteria.csv',index=False)
    e[e[gates].eq('unresolved').any(axis=1)].to_csv(out/'unresolved_programs.csv',index=False)
    count=lambda df:dict(programs=len(df),antigens=int(df.protein.nunique()))
    flow=[dict(step='all',**count(e))];keep=pd.Series(True,index=e.index)
    for c in gates:
        previous=e[keep]
        keep &= e[c].eq('pass')
        flow.append(dict(step=c,**count(e[keep]),removed_programs=len(previous)-int(keep.sum()),
          removed_antigens=previous.protein.nunique()-e.loc[keep,'protein'].nunique()))
    pd.DataFrame(flow).to_csv(out/'six_criteria_funnel.csv',index=False)
    brief={'protein':'抗原','genes':'基因','blood_state':'代表血细胞状态','brain_state':'代表脑细胞状态',
      'drug_examples':'临床药物示例','approved_examples':'曾获批且库中未标撤回的药物示例',
      'q_gene_decoy':'随机程序q','matched_brain_q':'同细胞脑疾病q','gwas_q':'程序遗传q_非门槛',
      'matched_program_abundance_hypothesis':'丰度调节假说_不等于药物效应',
      'protein_rho':'蛋白预测rho_非门槛','confirmation_positive_fraction':'种子同向比例_非门槛'}
    compact=nominated.drop_duplicates('protein')[list(brief)].rename(columns=brief)
    aggregate=nominated.groupby('protein').agg(支持程序数=('combo','size'),
      血细胞状态=('blood_state',lambda x:';'.join(sorted(set(x)))),
      脑细胞状态=('brain_state',lambda x:';'.join(sorted(set(x)))))
    compact=compact.merge(aggregate,left_on='抗原',right_index=True,validate='one_to_one')
    compact.to_csv(out/'六项通过候选汇总.csv',index=False,encoding='utf-8-sig')
    summary=dict(policy=POLICY,status='complete',source_sha256=h,total=count(e),flow=flow,
      five_criteria=count(e[e.five_criteria_pass]),six_criteria=count(nominated),
      approved_subset=count(e[e.six_criteria_approved_subset]),
      nominated_antigens=sorted(nominated.protein.unique().tolist()),
      approved_subset_antigens=sorted(e.loc[e.six_criteria_approved_subset,'protein'].unique().tolist()),
      note='Retrospective six-criterion evidence screen. A6 differs from historical A. Clinical phases are historical maximum phases, not current active/marketed status. No drug-specific brain response is claimed.',
      cd22=count(nominated[nominated.protein.eq('CD22')]),cd35=count(nominated[nominated.protein.eq('CD35')]))
    assert hashlib.sha256(source.read_bytes()).hexdigest()==h
    # Reconfirm all original measured quantities and P/q by key.
    old=pd.read_csv(source).set_index('combo').sort_index()
    now=e.set_index('combo').sort_index()
    for c in old.columns:
        if c=='tier':continue
        pd.testing.assert_series_equal(old[c],now[c],check_names=False)
    assert not set(['protein_rho','confirmation_positive_fraction','protein_percentile']).intersection(gates)
    (out/'completed.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
    approved_names=summary['approved_subset_antigens']
    clinical_only=sorted(set(summary['nominated_antigens'])-set(approved_names))
    cd35_five=e[e.protein.eq('CD35') & e.five_criteria_pass]
    text=[
      '六项标准全候选重筛：结果与复现说明', '2026-10-03', '',
      f"范围：{len(e)}个冻结程序、{e.protein.nunique()}个生物学抗原。前五项通过{summary['five_criteria']['antigens']}个抗原；加入已有临床药物证据后通过{summary['six_criteria']['antigens']}个抗原、{summary['six_criteria']['programs']}个程序，记为A6。",
      f"其中{len(approved_names)}个抗原有药物获批记录且ChEMBL未标记撤回；其余{len(clinical_only)}个只有临床开发记录。药物阶段是历史最高阶段，不等于当前仍在研或仍上市；这里不把研究试剂自动计作临床药物。",
      'A6与旧A的定义不同：GWAS、蛋白预测相关性、种子同向比例、单调性、抗原相对排名均不作为本次资格或排序条件。', '',
      '逐步筛选（程序数 / 抗原数）：']
    labels=['全部','1 身份与来源语境','2 实测剂量范围','3 程序响应','4 同细胞状态PD脑支持','5 外周可接近性','6 已有临床药物']
    text += [f"{label}：{r['programs']} / {r['antigens']}" for label,r in zip(labels,flow)]
    text += ['',
      '第1项：已知抗原/复合物/表位身份、匹配来源细胞中的直接ADT测量。明显的旁系细胞、血小板关联或来源语境冲突列为待核实，详见six_criteria_curation.json与C1_note；不设RNA检出率、RNA–ADT相关性或模型预测相关性的数值切点。测量可追溯不等同于已经排除所有混合事件或背景。',
      '第2项：沿用实测IQR>1e-6，排除退化剂量输入。',
      '第3项：冻结程序随机基因对照全家族q<=0.05。',
      '第4项：采用此前已计算的Kamath留出供体与GSE157783匹配细胞状态检验；两者可评估、实际方向一致、相对背景方向一致，合并全家族q<=0.05。Ma只作空间语境，不再参与方向否决。合并显著不等于两队列分别显著。',
      '第5项：UniProt细胞膜、分泌或胞外拓扑证据；明确细胞表面糖基表位单独保留。已有药物的具体机制列于第6项，不能把调节丰度直接等同于抑制/激动。',
      '第6项：从全部候选查询ChEMBL临床机制和阶段，依据UniProt accession匹配人类分子靶点，并与已有DGIdb记录核对。保留直接分子作用且有临床阶段记录的条目。复合物要核对结合亚基；泛家族命中、表位不对应、单纯Fc招募和蛋白作为药物支架的记录分别标记。', '',
      '具有获批药物记录的抗原：', '、'.join(approved_names), '',
      '仅有临床开发记录的抗原：', '、'.join(clinical_only), '',
      'CD22与CD35：',
      'CD22通过六项，保留cDC→MHC-II小胶质和cDC→星形胶质两条程序。直接药物包括Inotuzumab ozogamicin及有临床研究记录的Epratuzumab等；抗体偶联毒素、放射偶联和受体调节分别保留机制。两个程序的丰度方向假说不同，不能合并成一个不分脑程序的药物方向。',
      f"CD35有{len(cd35_five)}条程序通过前五项，其中{int(cd35_five.gwas_q.le(.05).sum())}条有程序遗传支持。本轮没有确认直接靶向内源CD35的临床药物；Mirococept是CR1片段构成的补体抑制剂，单独保留为功能模拟路线，不误写成CD35抑制剂。这个差异来自新增药物标准。", '',
      '结果的用途：',
      'A6是按六项计算和公开药理注释筛出的候选集，不是已经验证能治疗PD的药物靶点。各抗原与对应程序保持关联，不把不同程序上的遗传、疾病或响应优点拼成一条证据。已知细胞语境冲突标为未解决，不作生物学阴性。',
      f"{len(nominated)}个通过程序涉及的脑状态为：{';'.join(sorted(nominated.brain_state.unique()))}。稳态小胶质在既有同细胞状态验证中覆盖不足，不能据此解释为其无生物学作用。",
      '药物例子按获批记录、阶段及名称展示；候选按通过标准、随机程序q、同细胞脑q排序。该顺序不是药物疗效或安全性排名。模型性能和种子诊断仍完整保留。', '',
      '文件：',
      f"六项通过候选汇总.csv：{nominated.protein.nunique()}个抗原的中文汇总，含全部支持轴与代表程序。",
      'A6_approved_drug_subset.csv：有获批药物记录的子集。',
      'all_213_antigen_dispositions.csv：每个抗原最远通过记录与未满足条件。',
      'all_6849_programs_six_criteria.csv：全部程序的六项状态及原始证据。',
      'all_drug_target_mechanisms.csv：包含未接受记录与具体原因。',
      'verified_drug_target_mechanisms.csv：直接靶向且有临床开发记录的药物。',
      'related_pharmacology_routes.csv：CR1功能模拟与TREM1配体捕获等相关药理路线，单独保留。',
      'unresolved_programs.csv：身份/细胞语境或数据覆盖尚未解决的记录。',
      'CD22_CD35_six_criteria.csv：两个原候选的全部程序。', '',
      '复现：',
      '源码在上一级reproduction_source.tar.gz中；本地工作源码位于pd_product_figure_design/manuscript_review/optimization_20261003。',
      '离线重筛：python screen_six_criteria.py --result-dir <optimization_20261003结果目录>。',
      '获取参考数据：python fetch_six_criteria_evidence.py --result-dir <同目录>。缓存命中时复用原始返回；要刷新应使用新的版本目录，不覆盖此次药物证据。',
      '输入为服务器冻结统计结果、本项目实测身份表，以及逐请求记录URL和哈希的UniProt/ChEMBL数据。原表及全部P/q经逐列比对保持不变；没有重训、扩大剂量范围或重调显著性阈值。',
      '完整来源在reference_manifest.json、机制表的mechanism_refs、six_criteria_curation.json和manifest.json。', '',
      '关键资料：',
      'ChEMBL：https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services',
      'CD22抗体与外周DC：https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.b.20469',
      'Mirococept构成与作用：https://pmc.ncbi.nlm.nih.gov/articles/PMC1576234/',
      'Inotuzumab说明书：https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=cc7014b1-c775-411d-b374-8113248b4077']
    (out/'筛选结果与复现说明.txt').write_text('\n'.join(text)+'\n')
    (root/'active_nomination.json').write_text(json.dumps(dict(active_directory='six_criteria_20261003',
        policy_id=POLICY, preferred_term='peripherally druggable targets candidates',
        historical_rule_directories=['joint_analysis','current_nomination'],
        candidate_table='A6_antigens.csv', approved_subset='A6_approved_drug_subset.csv'),indent=2))
    for name in ['six_criteria_protocol.json','six_criteria_curation.json']:
        shutil.copyfile(Path(__file__).with_name(name),out/name)
    manifest={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob('*')) if p.is_file() and p.name!='manifest.json'}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps(summary,ensure_ascii=False,indent=2))
    # Restore only an independently verified spatial assessment of identical inputs.
    spatial_hook=Path(__file__).with_name('reapply_completed_spatial.py')
    if spatial_hook.exists() and (root/'seven_criteria_20261003/completed.json').exists():
        subprocess.run([sys.executable,str(spatial_hook),'--result-dir',str(root)],check=True)
    mr_hook=root/'mr_parallel_gate_20261003/apply_gate.py'
    if mr_hook.exists() and (root.parent/'blood_brain_MR_20261003/all_brain_MR_estimates.csv').exists():
        subprocess.run([sys.executable,str(mr_hook),'--results-root',str(root.parent)],check=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--result-dir',required=True)
    main(ap.parse_args().result_dir)
