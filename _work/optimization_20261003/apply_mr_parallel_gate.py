"""Join unchanged MR and perturbation FDR results before downstream prioritization."""
from pathlib import Path
import json,time,hashlib,argparse
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parents[3]
ap=argparse.ArgumentParser();ap.add_argument('--results-root',type=Path,default=ROOT/'delivery');args=ap.parse_args();A=args.results_root/'optimization_20261003';M=args.results_root/'blood_brain_MR_20261003';O=A/'mr_parallel_gate_20261003';O.mkdir(exist_ok=True)
t0=time.perf_counter();d=pd.read_csv(A/'six_criteria_20261003/all_6849_programs_six_criteria.csv');input_rows=len(d)
spatial_record=json.loads((A/'seven_criteria_20261003/completed.json').read_text())
if hashlib.sha256((A/'six_criteria_20261003/all_6849_programs_six_criteria.csv').read_bytes()).hexdigest()==spatial_record['source_six_criteria_sha256']:
 spatial=pd.read_csv(A/'seven_criteria_20261003/all_6849_programs_seven_criteria.csv')[['combo','C7']];d=d.merge(spatial,on='combo',validate='one_to_one')
else:d['C7']='unresolved_changed_program_inputs'
mr=pd.read_csv(M/'all_brain_MR_estimates.csv');reg=json.loads((A/'candidate_registry.json').read_text());mapping={a['target']:a['genes'] for a in reg}
# Single gene objects only. Multi-gene complexes/isoform assays do not inherit whole-gene MR automatically.
def status(target):
 genes=mapping.get(target,[])
 if len(genes)!=1:return dict(MR_gate='unresolved_assay_mapping',MR_gene=';'.join(genes),MR_q=np.nan,MR_outcomes='')
 gene=genes[0];z=mr[mr.gene.eq(gene)]
 if z.empty:return dict(MR_gate='unassessed_no_usable_MR',MR_gene=gene,MR_q=np.nan,MR_outcomes='')
 if target in ['CD45RA','CD45RB','CD45RO']:return dict(MR_gate='unresolved_isoform_mapping',MR_gene=gene,MR_q=float(z.q_BH.min()),MR_outcomes='')
 hit=z[z.q_BH.lt(.05)]
 return dict(MR_gate='pass' if len(hit) else 'tested_not_significant',MR_gene=gene,MR_q=float(z.q_BH.min()),MR_outcomes=';'.join(hit.outcome_trait),MR_source=';'.join(z.source.unique()),MR_n_outcomes=len(z))
lookup=pd.DataFrame([dict(protein=p,**status(p)) for p in sorted(d.protein.unique())]);d=d.merge(lookup,on='protein',validate='many_to_one')
d['perturbation_gate']=np.where(d.q_gene_decoy.le(.05),'pass',np.where(d.q_gene_decoy.notna(),'tested_not_significant','unassessed'))
d['MR_and_perturbation_pass']=d.MR_gate.eq('pass')&d.perturbation_gate.eq('pass')
d['priority_pass']=d.MR_and_perturbation_pass&d[['C1','C2','C4','C5','C6']].eq('pass').all(axis=1)
d['spatial_supported_priority']=d.priority_pass&d.C7.eq('pass')
d['priority_status']=np.where(d.priority_pass,'MR_and_perturbation_supported_priority',np.where(d.MR_gate.str.startswith(('unassessed','unresolved')),'MR_pending',np.where(d.MR_and_perturbation_pass,'joint_signal_pending_downstream_criteria','tested_not_jointly_supported')))
d.to_csv(O/'all_programs_parallel_gate.csv',index=False);lookup.to_csv(O/'all_antigen_MR_dispositions.csv',index=False)
intersection=d[d.MR_and_perturbation_pass];priority=d[d.priority_pass].sort_values(['MR_q','q_gene_decoy','protein','combo'])
intersection.to_csv(O/'MR_perturbation_intersection.csv',index=False);priority.to_csv(O/'priority_all_programs.csv',index=False);priority.drop_duplicates('protein').to_csv(O/'priority_antigens.csv',index=False)
d[d.MR_gate.str.startswith(('unassessed','unresolved'))].to_csv(O/'MR_pending_programs.csv',index=False)
priority[['combo','protein','MR_gene','axis','response_direction','q_gene_decoy','MR_q','MR_outcomes','C1','C2','C4','C5','C6','C7','drug_examples','predicted_restorative_protein_change']].to_csv(O/'priority_evidence.csv',index=False)
count=lambda z:dict(programs=len(z),antigens=int(z.protein.nunique()))
steps=[('all_frozen_programs',d),('perturbation_pass',d[d.perturbation_gate.eq('pass')]),('MR_pass',d[d.MR_gate.eq('pass')]),('MR_AND_perturbation',intersection)]
f=intersection
for gate in ['C1','C2','C4','C5','C6']:
 f=f[f[gate].eq('pass')];steps.append((gate,f))
pd.DataFrame([dict(step=s,**count(z)) for s,z in steps]).to_csv(O/'parallel_gate_flow.csv',index=False)
summary=dict(status='complete',all=count(d),MR_and_perturbation=count(intersection),priority=count(priority),priority_targets=sorted(priority.protein.unique()),spatial_supported_priority=count(d[d.spatial_supported_priority]),MR_available_antigens=int(lookup.MR_gate.isin(['pass','tested_not_significant']).sum()),MR_pending_antigens=int(lookup.MR_gate.str.startswith(('unassessed','unresolved')).sum()),elapsed_seconds=time.perf_counter()-t0,MR_threshold='q_BH<0.05 across original576 selected exposure-outcome tests',perturbation_threshold='q_gene_decoy<=0.05 from unchanged full frozen program family',MRI_direction_not_forced_against_cell_state=True,downstream_evidence_reused=True,models_retrained=False,spatial_rule='retain existing partial evaluation; not newly forced as a global gate',MR_scope='current59-target exposure search, not all213-antigen proteome-wide search')
(O/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
assert (d.C3.eq('pass')==d.perturbation_gate.eq('pass')).all();assert priority.MR_and_perturbation_pass.all();assert len(d)==input_rows;assert d.combo.is_unique
old=json.loads((A/'active_nomination.json').read_text());old.update(active_directory=O.name,policy_id='MR_parallel_perturbation_20261003',candidate_table='priority_antigens.csv',all_program_table='priority_all_programs.csv',pending_table='MR_pending_programs.csv',previous_active_directory='six_criteria_20261003',parallel_gates=['MR','perturbation'],spatial_layer_status='Existing partial cell-subtype evaluation retained; no new global C7 gate')
for k in ['approved_subset','spatial_supported_candidate_table','spatial_pending_program_table']:old.pop(k,None)
(A/'active_nomination.json').write_text(json.dumps(old,indent=2,ensure_ascii=False))
lines=['MR与扰动并列筛选（2026-10-03）','', '进入重点后续分析的逻辑：MR显著 AND 现有计算扰动显著；之后核对身份、可测剂量、同亚型疾病证据、外周可及性和药物证据。空间沿用原有部分可评估规则。', 'MR回答血端蛋白与整体脑影像表型的遗传联系；计算扰动回答特定血细胞来源到脑细胞亚型程序的模型响应。两者按靶点合流，保留具体程序和影像结局，不把两者当成同一表型的重复实验。', 'MR与扰动P值、FDR、发现程序均不因交集筛选重新计算。MR未覆盖不记阴性。MRI体积/T2*方向不强制与亚型程序表达同向。', '当前MR只检索了此前59个候选的工具；本次是对这些已有证据的重优先排序，不声称完成全部213个抗原的MR独立全筛。', '','计数：']
for s,z in steps:lines.append(f'{s}: {len(z)}个程序，{z.protein.nunique()}个抗原')
lines+=['','优先靶点：'+', '.join(summary['priority_targets']),f'目前空间支持的优先程序：{summary["spatial_supported_priority"]["programs"]}；逐项C7状态保留在priority_evidence.csv，不新增强制同向条件。',f'本轮用时{summary["elapsed_seconds"]:.2f}秒，只复用真实结果表，无模型训练、随机扰动或空间重算。', '建议：后续昂贵的靶点细化优先在交集上做；候选发现、程序定义、全家族FDR维持原全量范围。共定位可优先对MR显著的三个基因开展，现阶段未完成，不把它写成已通过。']
(O/'筛选结果与运行说明.txt').write_text('\n'.join(lines)+'\n')
(O/'manifest.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.iterdir() if p.is_file() and p.name!='manifest.json'},indent=2))
print(json.dumps(summary,ensure_ascii=False))
