"""Pilot list: rows already read by the earlier independent decoder (F2-r: symgz + streamed values), a few chrX rows, some near file start/end."""
import pandas as pd, numpy as np, sys
O='/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan/results/big40_repro_20261119/'; sys.path.insert(0,O+'scripts'); import rfetch as RF
R='/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/research_20261118/verify-mr_imaging-F2-r/'
t=pd.read_csv(R+'replication_targets.csv.gz',dtype={'idp':str},low_memory=False); r=pd.read_csv(R+'replication_raw.csv',dtype={'idp':str},low_memory=False); m=t.merge(r[r.status=='ok'][['exposure_id','phs','a1','a2','beta','se','mlogp','z']],on=['exposure_id','phs'])
cat=RF.catalog(); s=m.sample(60,random_state=7).copy(); L=[RF.find_row(cat,int(c),int(p),e,o) for c,p,e,o in zip(s.chrom,s.pos37,s.ea_x,s.oa_x)]; s['catalog_row']=[l[0] for l in L]; s['idp']=s.idp.str.zfill(4)
q=pd.read_pickle(O+'work/query_table.pkl'); x=q[(q.chrom==23)&q.testable].sample(8,random_state=1)
A=pd.DataFrame(dict(idp=s.idp,catalog_row=s.catalog_row,prio=0,n_rep=s.n_replication,chrom=s.chrom,pos37=s.pos37,a1=s.a1,a2=s.a2,release='repro')); B=pd.DataFrame(dict(idp=x.idp,catalog_row=x.catalog_row.astype(int),prio=1,n_rep=x.n_rep,chrom=x.chrom,pos37=x.pos37,a1=x.a1,a2=x.a2,release='repro'))
pd.concat([A,B]).to_csv(O+'work/pilot_targets.csv.gz',index=False); s.rename(columns={'z_y':'z_f2'})[['idp','catalog_row','beta','se','mlogp','z_f2','a1','a2']].to_csv(O+'work/pilot_truth.csv',index=False); print(len(A),len(B))
