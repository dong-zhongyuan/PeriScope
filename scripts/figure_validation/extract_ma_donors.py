import anndata as ad,numpy as np,json,base64,io,zipfile
from pathlib import Path
A=Path('/public/home/mengxl/dzy/pd_product_assets');d=ad.read_h5ad(A/'processed/gse253975/v0.1/GSE253975_geomx.h5ad')
programs=json.load(open("/public/home/mengxl/dzy/pd_product_assets/results/figure_program_registry/programs.json"))
names=list(map(str,d.var_names));wanted=sorted(set(sum(programs.values(),[]))|{'GFAP','AQP4','SLC1A3','ALDH1L1'});present=[g for g in wanted if g in names];indices=[names.index(g)for g in present];idx={g:i for i,g in enumerate(present)}
X=d.X;lib=np.asarray(X.sum(axis=1)).ravel();V=X[:,indices];V=V.toarray()if hasattr(V,'toarray')else np.asarray(V);L=np.log1p(V/np.maximum(lib[:,None],1)*1e6)
astro=[idx[g]for g in ['GFAP','AQP4','SLC1A3','ALDH1L1']if g in idx];a=L[:,astro];score=((a-a.mean(0))/(a.std(0)+1e-9)).mean(1)
donor=d.obs['donor'].astype(str).values;condition=d.obs['condition'].astype(str).values;ds=sorted(set(donor));Mfull=np.stack([L[donor==x].mean(0)for x in ds]);mu=Mfull.mean(0);sd=Mfull.std(0,ddof=1);rows=[]
for frac in [.3,.5,.7,1.]:
 for who in ds:
  mask=(donor==who);keep=mask&(score>=np.quantile(score[mask],1-frac));vec=L[keep].mean(0);z=(vec-mu)/(sd+1e-9)
  for ck,gs in programs.items():
   gi=[idx[g]for g in gs if g in idx];rows.append(dict(combo=ck,donor=who,condition=str(condition[mask][0]),roi_fraction=frac,n_roi=int(keep.sum()),n_genes=len(gi),score=round(float(z[gi].mean()),6)))
out=dict(method='per donor mean log1p CPM; fixed per-gene center/SD from 10 all-ROI donor pseudobulks; average over discovery_union_v1 genes; descriptive score, not original competitive-test statistic',unit='one donor; 5 PD and 5 Control',astro_markers=['GFAP','AQP4','SLC1A3','ALDH1L1'],records=rows)
b=io.BytesIO()
with zipfile.ZipFile(b,'w',zipfile.ZIP_DEFLATED)as z:z.writestr('ma_donor_strength.json',json.dumps(out,separators=(',',':')))
print('DONOR_START'+base64.b64encode(b.getvalue()).decode()+'DONOR_END',flush=True)
