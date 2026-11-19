import sys; sys.path.insert(0,'/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/drug_evidence_expanded_20261119/code')
from common import *
I=O+'/interim/'
EV=pd.read_csv(I+'evidence_long.csv.gz',low_memory=False).fillna({'target_genes':'','action':'','target_type':'','mechanism':''}); DM=pd.read_csv(I+'drug_master.csv.gz',low_memory=False).set_index('drug_key')
EV=EV[EV.source.isin(['OpenTargets_26.09','ChEMBL_REST_20261003','GtoPdb_2026.3.1','DrugCentral_2026-09','DGIdb_live'])]
EV['dirn']=[action_direction(a)[0] for a in EV.action]; EV['sil']=[action_direction(a)[1] for a in EV.action]
EV['stage']=EV.drug_key.map(DM[['ot_stage','chembl_phase']].max(axis=1)); EV['appr']=EV.drug_key.map(DM.gtop_approved|DM.dgidb_approved); EV['name']=EV.drug_key.map(DM.name); EV['dtype']=EV.drug_key.map(DM.drug_type)
EV['stage2']=np.where(EV.appr.fillna(False)&~(EV.stage>=0.5),4,EV.stage)
X=EV.assign(g=EV.target_genes.str.split(';')).explode('g')
X=X[X.target_type.isin(['single protein','ligand_target','protein complex','protein complex group','protein-protein interaction','chimeric protein'])&(X.stage2>=0.5)&(X.primary|X.source.isin(['OpenTargets_26.09','ChEMBL_REST_20261003'])|X.source.eq('DGIdb_live'))]
def summ(d):
    o=[]
    for k,y in d.groupby('drug_key'):
        dirs=set(y.dirn)-{'unknown'}; o.append((-y.stage2.max(),f"{str(y.name.iloc[0])[:18]}[{'/'.join(sorted(dirs)) or '?'};{y.stage2.max():g}]"))
    o.sort(); return len(o),' '.join(x[1] for x in o[:7])
LL=pd.read_csv(I+'ligand_link_candidates.csv').fillna(''); PP=pd.read_csv(I+'complex_link_candidates.csv').fillna('')
g={k:v for k,v in X.groupby('g')}
r=[]
for x in LL.itertuples():
    if x.ligand_gene in g:
        n,s=summ(g[x.ligand_gene]); r.append(dict(ant=x.antigens,R=x.receptor_gene,L=x.ligand_gene,src=x.link_sources.replace('GtoPdb_endogenous','G').replace('CellPhoneDB','C').replace('OmniPath','O'),ty=x.ligand_type[:30],n=n,drugs=s,refs=x.refs[:60]))
r=pd.DataFrame(r); r.to_csv(I+'review_ligand_pairs.tsv',sep='\t',index=False); print(len(r),r.R.nunique())
r2=[]
for x in PP.itertuples():
    if x.partner_gene in g:
        n,s=summ(g[x.partner_gene]); r2.append(dict(ant=x.antigens,G=x.antigen_gene,P=x.partner_gene,src=x.link_sources[:40],cx=x.complexes[:60],n=n,drugs=s))
r2=pd.DataFrame(r2); r2.to_csv(I+'review_complex_pairs.tsv',sep='\t',index=False); print(len(r2),r2.G.nunique())
