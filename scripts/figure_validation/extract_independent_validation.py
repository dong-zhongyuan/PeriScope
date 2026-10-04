import json,sys,io,zipfile,itertools
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
programs=json.load(open("/public/home/mengxl/dzy/pd_product_assets/results/figure_program_registry/programs.json"))
a=ad.read_h5ad('/public/home/mengxl/dzy/pd_product_assets/processed/gse157783/v0.1/gse157783_qc.h5ad')
var=a.var['gene_symbol'].astype(str).values;idx={g:i for i,g in enumerate(var)}
ct=a.obs.cell_type.astype(str).values;don=a.obs.donor.astype(str).values;cond=a.obs.condition.astype(str).values
rows=[];stats=[]
for combo,genes in programs.items():
 celltype='Astrocytes' if combo.endswith('__astro') else 'Microglia';donors=sorted(set(don[ct==celltype]));mat=[];conditions=[];ns=[]
 for who in donors:
  m=(ct==celltype)&(don==who);mat.append(np.asarray(a.X[m].sum(0)).ravel());conditions.append(cond[m][0]);ns.append(m.sum())
 mat=np.stack(mat);log=np.log1p(mat/np.maximum(mat.sum(1)[:,None],1)*1e4);gi=[idx[g] for g in genes if g in idx];raw=log[:,gi];sd=raw.std(0);valid=sd>1e-8;score=((raw[:,valid]-raw[:,valid].mean(0))/sd[valid]).mean(1);pdmask=np.array(conditions)=='PD';effect=score[pdmask].mean()-score[~pdmask].mean();null=[]
 for sel in itertools.combinations(range(len(score)),int(pdmask.sum())):
  mm=np.zeros(len(score),bool);mm[list(sel)]=True;null.append(score[mm].mean()-score[~mm].mean())
 pval=np.mean(np.abs(null)>=abs(effect)-1e-12)
 stats.append(dict(combo=combo,cell_type=celltype,n_genes=len(gi),n_nonconstant=int(valid.sum()),effect_PD_minus_control=effect,p_exact_two_sided=pval,n_permutations=len(null)))
 for j,who in enumerate(donors):rows.append(dict(combo=combo,donor=who,condition=conditions[j],score=score[j],n_cells=int(ns[j]),n_genes=len(gi),cohort='GSE157783'))
s=pd.DataFrame(stats);pv=s.p_exact_two_sided.values;o=np.argsort(pv);q=np.minimum.accumulate((pv[o]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1];s['q_BH_all_programs']=0.;s.loc[o,'q_BH_all_programs']=np.minimum(q,1)
b=io.BytesIO()
with zipfile.ZipFile(b,'w',zipfile.ZIP_DEFLATED) as z:
 z.writestr('independent_donor_scores.csv',pd.DataFrame(rows).to_csv(index=False));z.writestr('independent_validation_stats.csv',s.to_csv(index=False))
sys.stdout.buffer.write(b.getvalue())
