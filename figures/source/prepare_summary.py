"""Summary inputs: each evidence column is one unchanged frozen program."""
from pathlib import Path
import pandas as pd, json, hashlib
from collections import Counter
W=Path(__file__).resolve().parents[1]; ROOT=W.parents[2]
O=ROOT/'delivery/optimization_20261003'
e=pd.read_csv(O/'mr_parallel_gate_20261003/all_programs_parallel_gate.csv')
targets=sorted(e.loc[e.priority_pass,'protein'].unique())
layers=['C1','C2','C3','C4','C5','C6','MR_gate']
w=e[e.protein.isin(targets)].sort_values(['protein','blood_state','brain_state','program_definition','response_direction']).copy()
w['layers_passed']=w[layers].eq('pass').sum(axis=1)
assert (w.layers_passed.eq(7)==w.priority_pass).all()
w[['combo','protein','blood_state','brain_state','program_definition','response_direction',*layers,'layers_passed','priority_pass']].to_csv(W/'data/evidence_wheel.csv',index=False)
programs=json.loads((O/'joint_analysis/programs.json').read_text())
sets={t:set(g for c in e.loc[e.priority_pass&e.protein.eq(t),'combo'] for g in programs[c]) for t in targets}
membership=pd.DataFrame([dict(gene=g,**{t:g in sets[t] for t in targets}) for g in sorted(set.union(*sets.values()))])
membership.to_csv(W/'data/priority_gene_membership.csv',index=False)
for name in ['evidence_wheel','priority_gene_membership']:
    pd.read_csv(W/'data'/f'{name}.csv').to_csv(ROOT/'delivery/tables'/f'Figure7_{name}.csv',index=False)
report={'targets':targets,'programs':w.groupby('protein').size().to_dict(),'joint_pass':w.groupby('protein').priority_pass.sum().to_dict(),'genes':{t:len(s) for t,s in sets.items()},'gene_union':len(membership),'all_three_genes':int(membership[targets].all(axis=1).sum()),'criteria':layers,'data_unit':'One frozen program per angular column; no pooling of positive criteria across programs.'}
(W/'audit_density_20261004').mkdir(parents=True,exist_ok=True)
(W/'audit_density_20261004/summary_data_checks.json').write_text(json.dumps(report,indent=2))
print(report)
