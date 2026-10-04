"""Frozen-program spatial tests on subtype expression and seven-criterion join."""
from pathlib import Path
import json,hashlib,sys,os
from datetime import datetime,timezone
from itertools import combinations
import numpy as np,pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control

A=Path('/public/home/mengxl/dzy/pd_product_assets');P=Path(os.environ.get('PERISCOPE_RESULT_ROOT',str(A/'results/optimization_20261003')))
S=Path(os.environ.get('SPATIAL_OUTPUT_DIR',str(P/'spatial_subtype_20261003')));O=P/'seven_criteria_20261003';O.mkdir(exist_ok=True)
source=P/'six_criteria_20261003/all_6849_programs_six_criteria.csv'
source_hash=hashlib.sha256(source.read_bytes()).hexdigest();e=pd.read_csv(source)
protocol=json.loads((S/'protocol.json').read_text())
qc=pd.read_csv(S/'state_identifiability_QC.csv').set_index('state')
audit=S/'candidate_evaluation_started.json'
if not audit.exists():
    audit.write_text(json.dumps(dict(started_utc=datetime.now(timezone.utc).isoformat(),protocol_sha256=hashlib.sha256((S/'protocol.json').read_bytes()).hexdigest(),six_criteria_sha256=source_hash,reference_and_QC_fixed_before_candidate_evaluation=True),indent=2))
programs={}
for tag,rel in [('raw','programs.json'),('relative','target_specific/programs.json')]:
    programs.update({tag+'::'+k:v for k,v in json.loads((P/rel).read_text()).items()})
assert set(e.combo)==set(programs)
deconv=set((S/'inputs/deconvolution_genes.txt').read_text().splitlines())
evaluation=set((S/'inputs/evaluation_genes.txt').read_text().splitlines())
assert not evaluation.intersection(deconv)
md=pd.read_csv(S/'inputs/spatial_metadata.csv');samples=md[['sample_gsm','donor','condition']].drop_duplicates().set_index('sample_gsm')
results=[];gene_results=[];support=[];donor_scores=[];permutation_cache={}
mapping=protocol.get('recipient_mapping',{s:[s] for s in protocol['states']})
parent_of={leaf:parent for parent,leaves in mapping.items() for leaf in leaves}
for state in protocol['states']:
    parent=parent_of[state]
    columns=[];donors=[];conditions=[]
    for sample,row in samples.iterrows():
        d=S/'samples'/sample;path=d/'cell_state_log_expression.csv'
        if not path.exists():continue
        x=pd.read_csv(path,index_col=0);ok=pd.read_csv(d/'cell_state_gene_converged.csv',index_col=0)
        if state not in x:continue
        y=x[state].where(ok[state].fillna(False).astype(bool)&np.isfinite(x[state]))
        columns.append(y.rename(row.donor));donors.append(row.donor);conditions.append(row.condition=='PD')
    M=pd.concat(columns,axis=1).T if columns else pd.DataFrame()
    cond=np.array(conditions,dtype=bool)
    quality=state in qc.index and bool(qc.loc[state,'quality_pass'])
    enough=min(int(cond.sum()),int((~cond).sum()))>=3
    base=dict(brain_state=parent,spatial_subtype=state,state_identifiability_pass=quality,n_PD=int(cond.sum()),n_control=int((~cond).sum()))
    support.append(dict(**base,n_genes=int(M.shape[1]),donors=';'.join(donors)))
    if enough:
        M.to_csv(S/(state+'_donor_log_expression.csv'))
        x=M.to_numpy(dtype=float);n1=np.isfinite(x[cond]).sum(0);n0=np.isfinite(x[~cond]).sum(0)
        with np.errstate(invalid='ignore',divide='ignore'):
            mu1=np.nanmean(x[cond],axis=0);mu0=np.nanmean(x[~cond],axis=0)
            se=np.sqrt(np.nanvar(x[cond],axis=0,ddof=1)/n1+np.nanvar(x[~cond],axis=0,ddof=1)/n0)
            t=(mu1-mu0)/se
        usable=(n1>=3)&(n0>=3)&np.isfinite(t)&(se>1e-10)
        gs=M.columns.to_numpy()[usable];x=x[:,usable];t=t[usable]
        mean=np.nanmean(x,axis=0);sd=np.nanstd(x,axis=0,ddof=1)
        Z=(x-mean)/np.maximum(sd,1e-10);gi={g:i for i,g in enumerate(gs)}
        quant=np.unique(np.quantile(mean,np.linspace(0,1,21)));bins=np.digitize(mean,quant[1:-1]);nullcache={}
        for j,g in enumerate(gs):gene_results.append(dict(state=state,gene=g,mean_log_expression=mean[j],t=t[j],effect=mu1[usable][j]-mu0[usable][j],n_PD=int(n1[usable][j]),n_control=int(n0[usable][j])))
    for r in e[e.brain_state.eq(parent)].itertuples():
        members=programs[r.combo];rec=dict(combo=r.combo,protein=r.protein,**base,n_program_genes=len(members),n_used=0,coverage=0.,p_up=1.,p_down=1.,effect=np.nan)
        if not enough:
            rec['spatial_status']='insufficient_donors';results.append(rec);continue
        ii=np.array([gi[g] for g in members if g in gi],dtype=int);rec.update(n_used=len(ii),coverage=len(ii)/len(members))
        if len(ii)<4 or len(ii)/len(members)<.5:
            rec['spatial_status']='insufficient_program_gene_coverage';results.append(rec);continue
        # Every donor score requires half the used genes; no imputation of fits.
        zm=Z[:,ii];score=np.nanmean(zm,axis=1);score[np.isfinite(zm).sum(1)<max(4,int(np.ceil(len(ii)/2)))]=np.nan
        n1s=np.isfinite(score[cond]).sum();n0s=np.isfinite(score[~cond]).sum()
        if min(n1s,n0s)<3:
            rec['spatial_status']='insufficient_donor_program_scores';results.append(rec);continue
        # Descriptive exact donor-label test; it is not an extra candidate gate.
        valid=np.isfinite(score);v=score[valid];c=cond[valid];key=(len(v),int(c.sum()))
        if key not in permutation_cache:
            masks=np.zeros((len(list(combinations(range(key[0]),key[1]))),key[0]),dtype=bool)
            for j,ix in enumerate(combinations(range(key[0]),key[1])):masks[j,list(ix)]=True
            permutation_cache[key]=masks
        masks=permutation_cache[key]
        observed=float(v[c].mean()-v[~c].mean())
        perm=masks@v/key[1]-(~masks)@v/(key[0]-key[1])
        rec['spatial_donor_permutation_p']=float(np.mean(np.abs(perm)>=abs(observed)-1e-12))
        obs=float(t[ii].mean());null=np.zeros(5000)
        for b,k in zip(*np.unique(bins[ii],return_counts=True)):
            pool=np.where(bins==b)[0];key=(int(b),int(k))
            if key not in nullcache:
                rng=np.random.default_rng(42+1009*int(b)+9176*int(k))
                nullcache[key]=np.array([t[rng.choice(pool,int(k),replace=False)].sum() for _ in range(5000)])
            null+=nullcache[key]
        null/=len(ii)
        rec.update(spatial_status='ok' if quality else 'state_identifiability_not_met',
            effect=float(np.nanmean(score[cond])-np.nanmean(score[~cond])),competitive_stat=obs,null_mean=float(null.mean()),
            relative_effect=obs-float(null.mean()),p_up=float((1+(null>=obs).sum())/5001),p_down=float((1+(null<=obs).sum())/5001))
        results.append(rec)
        for donor,condition,value in zip(donors,cond,score):donor_scores.append(dict(combo=r.combo,state=state,donor=donor,condition='PD' if condition else 'Control',score=value))
s=pd.DataFrame(results)
assert len(s)==sum(len(mapping[st]) for st in e.brain_state)
s[['q_up','q_down']]=false_discovery_control(s[['p_up','p_down']].to_numpy().ravel()).reshape(-1,2)
direction=e.set_index('combo').cell_state_sensitivity_kamath_effect
s['matched_brain_effect']=s.combo.map(direction)
s['spatial_q']=np.where(s.matched_brain_effect>0,s.q_up,s.q_down)
s['leaf_direction_agrees']=s.effect*s.matched_brain_effect>0
s['leaf_C7']=np.where(s.spatial_status.ne('ok'),'unresolved',np.where(s.spatial_q.le(.05)&s.leaf_direction_agrees,'pass','not_supported'))
s.to_csv(S/'all_leaf_program_spatial_tests.csv',index=False)
# Preserve every leaf test, then aggregate using the prospectively fixed rule.
aggregated=[]
for combo,g in s.groupby('combo',sort=False):
    status='pass' if g.leaf_C7.eq('pass').any() else ('unresolved' if g.leaf_C7.eq('unresolved').any() else 'not_supported')
    choices=g[g.leaf_C7.eq(status)]
    r=choices.sort_values('spatial_q').iloc[0].to_dict()
    r['C7_aggregated']=status
    r['spatial_supported_subtypes']=';'.join(g.loc[g.leaf_C7.eq('pass'),'spatial_subtype'])
    r['spatial_unresolved_subtypes']=';'.join(g.loc[g.leaf_C7.eq('unresolved'),'spatial_subtype'])
    aggregated.append(r)
agg=pd.DataFrame(aggregated);assert len(agg)==len(e)
e=e.merge(agg,on=['combo','protein','brain_state'],validate='one_to_one')
same=e.effect*e.cell_state_sensitivity_kamath_effect>0
e['C7_direction_agrees_matched_brain']=same
e['old_mixed_spatial_direction_opposed_matched_brain']=e.ma_donor_score_effect*e.cell_state_sensitivity_kamath_effect<0
e['C7']=e.C7_aggregated
e['seven_criteria_pass']=e.six_criteria_pass&e.C7.eq('pass')
e['spatial_nomination_status']=np.where(e.seven_criteria_pass,'A7_spatial_supported',np.where(e.six_criteria_pass,np.where(e.C7.eq('unresolved'),'A6_spatial_unresolved','A6_spatial_not_supported'),'other_six_criterion_disposition'))
e['seven_criteria_approved_subset']=e.seven_criteria_pass&e.six_criteria_approved_subset
e['unmet_seven_criteria']=[(';'.join([r.unmet_six_criteria if isinstance(r.unmet_six_criteria,str) else '', 'C7:'+r.C7 if r.C7!='pass' else ''])).strip(';') for r in e.itertuples()]
e['criteria_pass_prefix_length_seven']=e[[f'C{i}' for i in range(1,8)]].eq('pass').astype(int).cumprod(axis=1).sum(axis=1)
e['antigen_C7_supported_programs']=e.groupby('protein').C7.transform(lambda x:int(x.eq('pass').sum()))
e['antigen_C7_unresolved_programs']=e.groupby('protein').C7.transform(lambda x:int(x.eq('unresolved').sum()))
e=e.sort_values(['seven_criteria_pass','seven_criteria_approved_subset','criteria_pass_prefix_length_seven','q_gene_decoy','matched_brain_q','protein','combo'],ascending=[False,False,False,True,True,True,True])
e.to_csv(O/'all_6849_programs_seven_criteria.csv',index=False)
e.drop_duplicates('protein').to_csv(O/'all_213_antigen_dispositions.csv',index=False)
a=e[e.seven_criteria_pass];a.to_csv(O/'A7_all_programs.csv',index=False);a.drop_duplicates('protein').to_csv(O/'A7_antigens.csv',index=False)
e[e.seven_criteria_approved_subset].drop_duplicates('protein').to_csv(O/'A7_approved_drug_subset.csv',index=False)
e[e.six_criteria_pass].to_csv(O/'six_to_seven_changes.csv',index=False)
e[e.protein.isin(['CD22','CD35'])].to_csv(O/'CD22_CD35_seven_criteria.csv',index=False)
e[e.five_criteria_pass&e.C7.eq('pass')].to_csv(O/'five_biological_plus_spatial_without_drug_requirement.csv',index=False)
e[e.C7.eq('unresolved')].to_csv(O/'C7_unresolved_programs.csv',index=False)
e[e.six_criteria_pass&e.C7.eq('unresolved')].to_csv(O/'A6_pending_spatial_resolution.csv',index=False)
e[e.six_criteria_pass&e.C7.eq('not_supported')].to_csv(O/'A6_tested_without_spatial_support.csv',index=False)
e[e.C7.eq('pass')&~e.six_criteria_pass].to_csv(O/'spatial_supported_previously_excluded_programs.csv',index=False)
e[['combo','protein','blood_state','brain_state','six_criteria_pass',
   'ma_donor_score_effect','cell_state_sensitivity_kamath_effect','effect',
   'old_mixed_spatial_direction_opposed_matched_brain','C7_direction_agrees_matched_brain',
   'spatial_status','spatial_q','C7','seven_criteria_pass']].to_csv(O/'mixed_tissue_vs_matched_state.csv',index=False)
pd.DataFrame(donor_scores).to_csv(S/'spatial_cell_state_donor_program_scores.csv',index=False)
pd.DataFrame(gene_results).to_csv(S/'spatial_cell_state_gene_effects.csv',index=False)
pd.DataFrame(support).to_csv(S/'state_expression_support.csv',index=False)
s.to_csv(S/'all_program_spatial_tests.csv',index=False)
count=lambda x:dict(programs=len(x),antigens=int(x.protein.nunique()))
flow=[dict(step='all',**count(e))];keep=pd.Series(True,index=e.index)
for i in range(1,8):
    before=e[keep];keep &= e[f'C{i}'].eq('pass')
    flow.append(dict(step=f'C{i}',**count(e[keep]),removed_programs=len(before)-int(keep.sum()),removed_antigens=before.protein.nunique()-e.loc[keep,'protein'].nunique()))
pd.DataFrame(flow).to_csv(O/'seven_criteria_funnel.csv',index=False)
# A global seventh gate must be usable for every currently nominated state.
# A technically unresolved state must not disappear through global adoption.
required_states=sorted(e.loc[e.six_criteria_pass,'brain_state'].unique())
usable_states=sorted(e.loc[e.spatial_status.eq('ok'),'brain_state'].unique())
required_subtypes=sorted({leaf for parent in required_states for leaf in mapping[parent]})
usable_subtypes=sorted(s.loc[s.spatial_status.eq('ok'),'spatial_subtype'].unique())
# A partially observed parent cannot support global negative exclusion.
# Positive evidence in one well-resolved leaf still supports the A7 tier.
adopt=set(required_subtypes)<=set(usable_subtypes)
summary=dict(status='complete',method='RCTD full + C-SIDE donor intercepts',spatial_layer_adoptable=adopt,
    scope='existing six criteria unchanged; C7 applied to the same program and recipient state',
    required_states_for_global_adoption=required_states,usable_states=usable_states,
    required_subtypes_for_global_adoption=required_subtypes,usable_subtypes=usable_subtypes,
    matched_parent_states_analyzable=set(required_states)<=set(usable_states),
    adoption_reason='All mapped subtypes covered' if adopt else 'Some mapped subtypes unresolved; retain six-criterion screen and report strict seven-criterion supported tier',
    source_six_criteria_sha256=source_hash,protocol_sha256=hashlib.sha256((S/'protocol.json').read_bytes()).hexdigest(),
    total=count(e),A6=count(e[e.six_criteria_pass]),A7=count(a),approved_subset=count(e[e.seven_criteria_approved_subset]),flow=flow,
    C7_status=e.C7.value_counts().to_dict(),A6_C7_status=e[e.six_criteria_pass].C7.value_counts().to_dict(),
    A7_antigens=sorted(a.protein.unique()),quality=qc.reset_index().to_dict('records'),state_coverage=support,
    A7_previously_opposed_in_mixed_tissue=count(a[a.old_mixed_spatial_direction_opposed_matched_brain]),
    CD22=count(a[a.protein.eq('CD22')]),CD35=count(a[a.protein.eq('CD35')]))
(O/'completed.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2));(S/'completed.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
old=pd.read_csv(source).set_index('combo').sort_index();new=e.set_index('combo').sort_index()
for c in old.columns:pd.testing.assert_series_equal(old[c],new[c],check_names=False)
assert hashlib.sha256(source.read_bytes()).hexdigest()==source_hash
if adopt:
    pointer=json.loads((P/'active_nomination.json').read_text());pointer.update(active_directory='seven_criteria_20261003',
        policy_id='20261003_seven_criteria_spatial_subtype',candidate_table='A7_antigens.csv',approved_subset='A7_approved_drug_subset.csv',
        spatial_evidence_directory=S.name,pending_spatial_table='A6_pending_spatial_resolution.csv',previous_six_criteria_directory='six_criteria_20261003')
    (P/'active_nomination.json').write_text(json.dumps(pointer,indent=2))
else:
    pointer=json.loads((P/'active_nomination.json').read_text())
    for obsolete in ['pending_spatial_table','previous_six_criteria_directory']:pointer.pop(obsolete,None)
    pointer.update(active_directory='six_criteria_20261003',policy_id='20261003_six_biological_and_drug_criteria',
        candidate_table='A6_antigens.csv',approved_subset='A6_approved_drug_subset.csv',
        spatial_evidence_directory=S.name,spatial_supported_candidate_table='A7_antigens.csv',
        spatial_pending_program_table='A6_pending_spatial_resolution.csv',
        spatial_assessment_directory='seven_criteria_20261003',
        spatial_layer_status='Partial validation; global six-criterion screen retained')
    (P/'active_nomination.json').write_text(json.dumps(pointer,indent=2))
brief=a.drop_duplicates('protein')[['protein','genes','blood_state','brain_state','q_gene_decoy','matched_brain_q','spatial_q','effect','drug_examples','approved_examples']].copy()
brief=brief.rename(columns={'protein':'抗原','genes':'基因','blood_state':'血来源状态','brain_state':'脑受体状态',
    'q_gene_decoy':'响应q','matched_brain_q':'单细胞脑q','spatial_q':'亚型空间q','effect':'亚型空间PD减对照',
    'drug_examples':'临床药物示例','approved_examples':'获批记录示例'})
brief.to_csv(O/'七项通过候选汇总.csv',index=False,encoding='utf-8-sig')
report=['亚型分辨空间验证与七项筛选','',
    '输入：服务器原始空间计数11859点位、10供者；RCTD完整混合模式；C-SIDE按供者估计状态内表达。',
    '空间参考使用Ma原文作者标注：Fibrous/Protoplasmic星形胶质、Baseline/Activated小胶质；Monocyte_Derived BAM独立建模，并保留其他脑谱系。',
    '参考排除与空间供者重叠者及身份未明确的OJ5，四位独立供者用于质量检查；亚型映射有原文注释和标记表达记录。',
    '组成标记与实际用于程序检验的基因完全分开。具体标记选择见本轮protocol.json；保留原程序成员，按剩余可评估基因检查覆盖。',
    '400个已知组成的参考混合样本仅用于可辨识性检查，不作为项目生物学证据。候选检验只使用真实空间供者。','',
    '质量检查：']
for st in protocol['states']:
    if st in qc.index:
        q=qc.loc[st];report.append(f"{st}: rho={q.spearman:.3f}, RMSE={q.RMSE:.3f}, 可辨识性工作标准通过={bool(q.quality_pass)}")
report+=['','C7条件：对应亚型质量合格、每组至少3个可评估供者、至少4个且覆盖原程序至少一半的基因、表达匹配随机程序全家族q<=0.05、实际疾病方向与此前同状态Kamath方向一致。',
    '可辨识性rho>=0.5、RMSE<=0.15及覆盖门槛是此次提前写入的分析工作标准，不是通用生物学真值界线。',
    '原C1-C6及所有旧p/q保持不变；不按本轮通过数量修改标准。C7缺覆盖/不可分辨记unresolved，已完成检验但未支持记not_supported。','',
    f"六项：{summary['A6']['antigens']}抗原/{summary['A6']['programs']}程序。",
    f"七项：{summary['A7']['antigens']}抗原/{summary['A7']['programs']}程序。",
    f"其中有获批药物记录：{summary['approved_subset']['antigens']}抗原。",
    f"已有六项候选涉及的脑状态：{required_states}；空间可评估状态：{usable_states}。",
    f"作为全局硬门槛需要覆盖其映射亚型：{required_subtypes}；实际可评估亚型：{usable_subtypes}。",
    f"能否作为全局新增筛选层：{adopt}。未覆盖全部映射亚型时，七项表仅表示空间支持子集，原六项仍是全局入口。",
    '原六项程序的C7状态：'+json.dumps(summary['A6_C7_status'],ensure_ascii=False),
    '七项抗原：'+'、'.join(summary['A7_antigens']),'',
    'CD22、CD35详见逐程序表。CD35若未满足直接药物条件，即使空间支持也不计入完整A7；前五项加空间的结果另表保留。',
    '前一轮人工来源语境C1的复核问题仍在，不能把本轮空间分析当成已经解决该问题。','',
    '复现顺序：prepare_spatial_inputs.py导出原空间计数 → prepare_author_reference.R → freeze_author_protocol.py → launch_spatial.py → evaluate_seventh_criterion.py → plot_spatial_review.R；最终分支SPATIAL_OUTPUT_DIR=spatial_subtype_author_20261003。',
    '完整环境、输入角色、固定参数、spacexr版本及数据哈希随交付保存。',
    '方法来源：https://github.com/dmcable/spacexr；https://pmc.ncbi.nlm.nih.gov/articles/PMC10463137/']
(O/'七项筛选结果与说明.txt').write_text('\n'.join(report)+'\n')
print(json.dumps(summary,ensure_ascii=False,indent=2))
