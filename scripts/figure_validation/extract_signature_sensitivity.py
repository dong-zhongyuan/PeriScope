import sys,io,zipfile,json
import anndata as ad,numpy as np,pandas as pd
ranked=json.load(open('/public/home/mengxl/dzy/pd_product_assets/results/figure_program_registry/discovery_rankings.json'))['CD22__cDC__x__astro']
a=ad.read_h5ad('/public/home/mengxl/dzy/pd_product_assets/processed/gse157783/v0.1/gse157783_qc.h5ad');var=a.var.gene_symbol.astype(str).values;idx={g:i for i,g in enumerate(var)};genes=[g for g in ranked if g in idx];obs=a.obs;donors=sorted(obs.loc[obs.cell_type=='Astrocytes','donor'].astype(str).unique());mat=[];conditions=[]
for who in donors:
 m=(obs.cell_type=='Astrocytes')&(obs.donor==who);v=np.asarray(a.X[m.values].sum(0)).ravel();mat.append(np.log1p(v/max(v.sum(),1)*1e4));conditions.append(str(obs.loc[m,'condition'].iloc[0]))
mat=np.stack(mat)[:,[idx[g] for g in genes]];sd=mat.std(0);Z=(mat-mat.mean(0))/np.maximum(sd,1e-8);rows=[]
for k in sorted(set([5,10,20,min(30,len(genes)),len(genes)])):
 for j,who in enumerate(donors):rows.append(dict(donor=who,condition=conditions[j],k=k,score=Z[j,:k].mean(),n_genes_total=len(genes)))
print(pd.DataFrame(rows).to_csv(index=False),end='')
