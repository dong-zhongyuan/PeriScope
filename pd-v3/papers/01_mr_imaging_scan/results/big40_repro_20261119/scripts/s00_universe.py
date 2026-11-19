"""Step 0: read one policy of the frozen full-atlas family table (read-only) and keep the estimated rows with the columns needed here."""
import pandas as pd, numpy as np, sys, time
OUT='/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan/results/big40_repro_20261119/'
A='/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan/results/full_atlas_20261116/full/'
pol=sys.argv[1]; t=time.time()
Q=['q_bh_global','q_by_global','q_storey_lambda_0p25','q_storey_lambda_0p5','q_storey_lambda_0p75','q_bh_cis_family','q_by_cis_family','q_bh_family']
U=['exposure_id','gene','assay','phenotype','family','status','variant','chrom','pos37','ea_x','oa_x','beta_x','se_x','eaf_x','ea_y','oa_y','beta_y','se_y','eaf_y','n_y','beta_y_aligned','eaf_y_aligned','harmonization_status','allele_flipped','hm_code','orientation_recovery','beta','se','log_p','F']+Q
ty={c:'category' for c in ['exposure_id','gene','assay','phenotype','family','status','variant','ea_x','oa_x','ea_y','oa_y','harmonization_status','hm_code','orientation_recovery']}
parts=[]; cnt={}; n=0
for ch in pd.read_csv(A+f'full_family_palindrome_{pol}.csv.gz',usecols=U,dtype=ty,chunksize=1500000,low_memory=False):
    n+=len(ch); vc=ch.status.astype(str).value_counts()
    for k,v in vc.items(): cnt[k]=cnt.get(k,0)+int(v)
    e=ch[ch.status.astype(str)=='estimated'].copy()
    for c in ty: e[c]=e[c].astype(str)
    parts.append(e); print(pol,n,round(time.time()-t),flush=True)
d=pd.concat(parts,ignore_index=True); del parts
for c in ty:
    if c!='status': d[c]=d[c].astype('category')
d=d.drop(columns=['status']); d.to_pickle(OUT+f'work/est_{pol}.pkl',protocol=4)
print('DONE',pol,'rows',len(d),'status',cnt,'sec',round(time.time()-t),flush=True)
