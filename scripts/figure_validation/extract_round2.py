import os,sys,io,json,zipfile,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import spearmanr
root='/public/home/mengxl/dzy/pd_product_assets';out={}
def put(n,x):out[n]=x.to_csv(index=False).encode()
a=ad.read_h5ad(root+'/processed/gse253975/v0.1/GSE253975_geomx.h5ad')
programs=json.load(open("/public/home/mengxl/dzy/pd_product_assets/results/figure_program_registry/programs.json"))
cases=['CD22__cDC__x__astro','CD35__nonclassical_mono__x__astro','CD62P__classical_mono__x__astro']
markers=['GFAP','AQP4','SLC1A3','ALDH1L1','P2RY12','TYROBP','MBP','TH']
genes=sorted(set(sum([programs[k] for k in cases],[])+markers));genes=[g for g in genes if g in a.var_names];idx=[a.var_names.get_loc(g) for g in genes]
X=a.X.tocsr();V=X[:,idx].toarray();total=np.asarray(X.sum(1)).ravel();L=np.log1p(V/np.maximum(total[:,None],1)*1e6)
don=a.obs.donor.astype(str).to_numpy();ds=sorted(set(don));M=np.stack([L[don==d].mean(0) for d in ds]);mu=M.mean(0);sd=M.std(0,ddof=1);z=(L-mu)/(sd+1e-9)
d=a.obs.copy();d.insert(0,'spot_id',a.obs_names)
for g in markers:d[g]=L[:,genes.index(g)]
for k in cases:d[k]=z[:,[genes.index(g) for g in programs[k] if g in genes]].mean(1)
put('spatial_expression_scores.csv',d)
out['spatial_score_definition.json']=json.dumps({'source':root+'/processed/gse253975/v0.1/GSE253975_geomx.h5ad','source_sha256':__import__('subprocess').check_output(['sha256sum',root+'/processed/gse253975/v0.1/GSE253975_geomx.h5ad']).decode().split()[0],'registry':'/public/home/mengxl/dzy/pd_product_assets/results/figure_program_registry/programs.json','registry_sha256':hashlib.sha256(open('/public/home/mengxl/dzy/pd_product_assets/results/figure_program_registry/programs.json','rb').read()).hexdigest(),'normalization':'log1p CPM per spot','programs':{k:programs[k] for k in cases},'scale':'each gene centered and scaled by mean and sample SD of the10 all-spot donor mean logCPM values; then mean frozen members per spot','n_spots':a.n_obs,'markers':markers},indent=2).encode()
# Independent single-cell accuracy for each protein, donor and cell type.
import torch
sys.path.insert(0,'/public/home/mengxl/dzy/pd_product/src')
from pdproduct.simulators.ccwm import CCWM,CCWMConfig
torch.set_num_threads(4);s=torch.load(root+'/checkpoints/project/ccwm-v23-s42/model.pt',map_location='cpu',weights_only=False);model=CCWM(CCWMConfig(**s['cfg']));model.load_state_dict(s['state_dict']);model.eval()
f0=np.load(root+'/interim/ccwm_v23/citeseq.npz',allow_pickle=True);f={k:f0[k] for k in ['X','Y','donor','ct','is_test']};rng=np.random.default_rng(6674);rows=[];sample=[]
for donor in sorted(set(f['donor'][f['is_test']])):
 for ct in sorted(set(f['ct'][f['is_test']])):
  ix=np.flatnonzero(f['is_test']&(f['donor']==donor)&(f['ct']==ct));ix=np.sort(rng.choice(ix,min(150,len(ix)),replace=False));sample.extend(ix.tolist())
  with torch.no_grad():mu,_=model.enc_blood(torch.tensor(f['X'][ix],dtype=torch.float32));pred=model.protein_head(mu).numpy()
  for j in range(f['Y'].shape[1]):
   real=f['Y'][ix,j];pr=pred[:,j];valid=np.ptp(real)>1e-10 and np.ptp(pr)>1e-10
   rows.append({'donor':'P'+str(int(donor)+1),'cell_type_index':int(ct),'protein_index':j,'n_cells':len(ix),'rho':float(spearmanr(real,pr).statistic) if valid else np.nan,'estimable':valid,'observed_sd':float(np.std(real))})
put('protein_single_cell_accuracy.csv',pd.DataFrame(rows));put('protein_single_cell_sample.csv',pd.DataFrame({'row_index':sample}))
b=io.BytesIO()
with zipfile.ZipFile(b,'w',zipfile.ZIP_DEFLATED) as zz:
 for k,v in out.items():zz.writestr(k,v)
sys.stdout.buffer.write(b.getvalue())
