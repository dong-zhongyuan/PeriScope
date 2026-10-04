from evidence_cache import hit,finish,acquire
lock=acquire(__file__)
if hit(__file__):
    print('IDENTICAL NEW-RUN EVIDENCE ALREADY COMPLETE',flush=True)
    raise SystemExit(0)
import json,sys,io,zipfile,itertools
import numpy as np,pandas as pd,anndata as ad
from pathlib import Path
from input_utils import canonical_symbols
from external_gene_symbols import external_symbols
O=Path("/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003/annotation_control")
from scipy import sparse
programs=json.load(open("/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003/annotation_control/programs.json"))
a=ad.read_h5ad('/public/home/mengxl/dzy/pd_product_assets/processed/gse157783/v0.1/gse157783_qc.h5ad')
var=external_symbols(a.var)
unique=sorted(set(var));idx={g:i for i,g in enumerate(unique)}
collapse=sparse.csr_matrix((np.ones(len(var)),(np.arange(len(var)),[idx[g] for g in var])),shape=(len(var),len(unique)))
X=a.X@collapse
ct=a.obs.cell_type.astype(str).values;don=a.obs.donor.astype(str).values;cond=a.obs.condition.astype(str).values
# Subtype assignments are the existing independently built purification maps.
subtype=pd.Series('other',index=a.obs_names)
for family in ['brain_astro','brain_microglia']:
 mapping=pd.read_csv('/public/home/mengxl/dzy/pd_product_assets/interim/v0.1/purification/c157_'+family+'.csv',index_col=0)
 mapping=mapping[mapping.kept];common=subtype.index.intersection(mapping.index)
 subtype.loc[common]=mapping.loc[common,'state_pure'].astype(str)
ct=subtype.to_numpy()
rows=[];stats=[];cache={};perm_cache={};size_rows=[];size_stats=[]
rankings=json.loads((O/'program_gene_rankings.json').read_text())
for combo,genes in programs.items():
 celltype=combo.split('__x__')[1].split('__')[0]
 if celltype not in cache:
  donors=sorted(set(don[ct==celltype]));mat=[];conditions=[];ns=[]
  for who in donors:
   m=(ct==celltype)&(don==who)
   if m.sum()<25:continue
   mat.append(np.asarray(X[m].sum(0)).ravel());conditions.append(cond[m][0]);ns.append(m.sum())
  donors=[who for who in donors if ((ct==celltype)&(don==who)).sum()>=25]
  if not mat or min(sum(c=='PD' for c in conditions),sum(c!='PD' for c in conditions))<2:cache[celltype]=None
  else:
   mat=np.stack(mat);log=np.log1p(mat/np.maximum(mat.sum(1)[:,None],1)*1e4);cache[celltype]=(donors,conditions,ns,log)
 if cache[celltype] is None:continue
 donors,conditions,ns,log=cache[celltype]
 gi=[idx[g] for g in genes if g in idx];raw=log[:,gi];sd=raw.std(0);valid=sd>1e-8
 if valid.sum()<4:continue
 score=((raw[:,valid]-raw[:,valid].mean(0))/sd[valid]).mean(1);pdmask=np.array(conditions)=='PD';effect=score[pdmask].mean()-score[~pdmask].mean();null=[]
 key=(len(score),int(pdmask.sum()))
 if key not in perm_cache:
  weights=[]
  for sel in itertools.combinations(range(len(score)),int(pdmask.sum())):
   mm=np.zeros(len(score),bool);mm[list(sel)]=True;weights.append(mm/mm.sum()-(~mm)/(~mm).sum())
  perm_cache[key]=np.stack(weights)
 null=perm_cache[key]@score
 pval=np.mean(np.abs(null)>=abs(effect)-1e-12)
 ranked=[g for g in rankings[combo] if g in idx and log[:,idx[g]].std()>1e-8]
 for k in sorted(set([min(v,len(ranked)) for v in [5,10,20,30]]+[len(ranked)])):
  rr=log[:,[idx[g] for g in ranked[:k]]];ss=((rr-rr.mean(0))/rr.std(0)).mean(1);ee=float(ss[pdmask].mean()-ss[~pdmask].mean())
  pp=float(np.mean(np.abs(perm_cache[key]@ss)>=abs(ee)-1e-12))
  size_stats.append(dict(combo=combo,k=k,n_available=len(ranked),effect_PD_minus_control=ee,p_exact_two_sided=pp))
  for j,who in enumerate(donors):size_rows.append(dict(combo=combo,k=k,donor=who,condition=conditions[j],score=float(ss[j])))
 stats.append(dict(combo=combo,cell_type=celltype,n_genes=len(gi),n_nonconstant=int(valid.sum()),effect_PD_minus_control=effect,p_exact_two_sided=pval,n_permutations=len(null)))
 for j,who in enumerate(donors):rows.append(dict(combo=combo,donor=who,condition=conditions[j],score=score[j],n_cells=int(ns[j]),n_genes=len(gi),cohort='GSE157783'))
s=pd.DataFrame(stats,columns=['combo','cell_type','n_genes','n_nonconstant','effect_PD_minus_control','p_exact_two_sided','n_permutations']);pv=s.p_exact_two_sided.values;o=np.argsort(pv);q=np.minimum.accumulate((pv[o]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1];s['q_BH_all_programs']=0.;s.loc[o,'q_BH_all_programs']=np.minimum(q,1)
pd.DataFrame(rows).to_csv(O/'independent_donor_scores.csv',index=False)
s.to_csv(O/'independent_validation_stats.csv',index=False)
from randko import bh
ss=pd.DataFrame(size_stats,columns=['combo','k','n_available','effect_PD_minus_control','p_exact_two_sided']);ss['q_all_program_sizes']=bh(ss.p_exact_two_sided)
ss.to_csv(O/'signature_size_sensitivity.csv',index=False)
pd.DataFrame(size_rows,columns=['combo','k','donor','condition','score']).to_csv(O/'signature_size_donor_scores.csv',index=False)
(O/'signature_size_definition.json').write_text(json.dumps(dict(ranking='signed mean actual-dose slope across discovery seeds42-46, within frozen program genes',sizes='5,10,20,30 and full detectable program, capped at available nonconstant genes',normalization='fixed gene scaling across all donors; no gene-set-size optimization',family='all program-size combinations'),indent=2))
print('independent programs',len(s),'size sensitivity',len(ss),flush=True)

finish(__file__,['independent_validation_stats.csv', 'independent_donor_scores.csv','signature_size_sensitivity.csv','signature_size_donor_scores.csv','signature_size_definition.json'])
