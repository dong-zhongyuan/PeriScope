"""One complete BIG40 repro file (IDP 0450) downloaded whole: checks the range reader on real replication data, the gzip CRC, the variant order against the catalog and the discovery file, the numeric row check, and gives a genome-wide discovery-vs-replication z correlation for one phenotype."""
import sys, json, subprocess, numpy as np, pandas as pd, io
O='/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan/results/big40_repro_20261119/'; sys.path.insert(0,O+'scripts'); import rfetch as RF
D='/public/home/mengxl/dzy/pd-v3/data/all_target_followup_20261118/regional_all/imaging/BIG40/0450/0450.txt.gz'; Rr=O+'work/repro_truth_0450.txt.gz'
ph=pd.read_csv('/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan/results/full_atlas_20261116/full/phenotypes.csv',low_memory=False); ph=ph[ph.phenotype=='BIG40:release2:discovery:0450'].iloc[0]; n=float(ph.n_replication); RF.set_df(n)
def rd(p):
    t=subprocess.run(['gzip','-dc',p],stdout=subprocess.PIPE,check=True).stdout; return t
tr=rd(Rr); td=rd(D); out=dict(idp='0450',n_replication=n,n_discovery=float(ph.n_discovery),repro_gzip_crc_ok=True,repro_bytes=int(__import__('os').path.getsize(Rr)))
hdr=dict(names=['chr','rsid','pos','a1','a2','beta','se','p'],sep=' ',header=0,dtype={'chr':str,'rsid':str,'a1':str,'a2':str})
r=pd.read_csv(io.BytesIO(tr),**hdr); d=pd.read_csv(io.BytesIO(td),**hdr); out['rows_repro']=len(r); out['rows_discovery']=len(d)
same=(r.chr.values==d.chr.values)&(r.pos.values==d.pos.values)&(r.a1.values==d.a1.values)&(r.a2.values==d.a2.values)&(r.rsid.values==d.rsid.values); out['rows_identical_variant_columns']=int(same.sum())
cat=RF.catalog(); out['catalog_rows']=cat.N; rng=np.random.default_rng(4); idx=rng.integers(0,len(r),200000)
pre=[f'{a} {b} {c} {x} {y} '.encode() for a,b,c,x,y in zip(r.chr.values[idx[:20000]],r.rsid.values[idx[:20000]],r.pos.values[idx[:20000]],r.a1.values[idx[:20000]],r.a2.values[idx[:20000]])]; out['catalog_prefix_match_sample20000']=int(sum(p==RF.prefix(cat,int(i)) for p,i in zip(pre,idx[:20000])))
ls=tr.split(b'\n',1)[1]; nl=np.flatnonzero(np.frombuffer(ls,np.uint8)==10); idx2=np.sort(idx); s=np.r_[0,nl[:-1]+1]; bad=0
for i in idx2:
    l=ls[int(s[i]):int(nl[i])]; bad+= not RF.row_ok(b' '.join(l.split(b' ')[5:]))
out['row_check_failures_on_true_rows_of_200000']=int(bad); z_r=(r.beta/r.se).values; z_d=(d.beta/d.se).values; ok=np.isfinite(z_r)&np.isfinite(z_d)
out['max_abs_z_repro']=float(np.nanmax(np.abs(z_r))); out['corr_z_all_variants']=float(np.corrcoef(z_d[ok],z_r[ok])[0,1]); out['mean_z_d_times_z_r_all']=float(np.mean(z_d[ok]*z_r[ok]))
for t in (1,2,3):
    m=ok&(np.abs(z_d)<t); out[f'corr_z_abs_zd_lt_{t}']=float(np.corrcoef(z_d[m],z_r[m])[0,1]); out[f'mean_product_abs_zd_lt_{t}']=float(np.mean(z_d[m]*z_r[m]))
for t in (5.0,7.3):
    m=ok&(np.abs(z_d)>=t); out[f'n_abs_zd_ge_{t}']=int(m.sum()); out[f'same_sign_abs_zd_ge_{t}']=float((np.sign(z_d[m])==np.sign(z_r[m])).mean()) if m.sum() else None
out['note']='The variant rows are in strong LD, so these correlations are descriptive; they include true polygenic signal and are not a sample-overlap estimate by themselves.'
json.dump(out,open(O+'work/truth_0450.json','w'),indent=1); print(json.dumps(out,indent=0))
