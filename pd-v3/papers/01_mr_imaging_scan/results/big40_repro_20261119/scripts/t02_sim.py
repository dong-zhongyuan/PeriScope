"""Offline truth test of rfetch.fetch_targets: complete local file (discovery or the whole repro file) versus decoded rows."""
import sys, time, subprocess, numpy as np, pandas as pd
O='/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan/results/big40_repro_20261119/'; sys.path.insert(0,O+'scripts')
import rfetch as RF, pfz; F=RF.F
from multiprocessing import Pool
path,n,NR,L,E=sys.argv[1],float(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])*1024,int(sys.argv[5])*1024; TAG=sys.argv[6]; seed=int(sys.argv[7]) if len(sys.argv)>7 else 11
cat=RF.catalog(); RF.set_df(n)
T=subprocess.run(['gzip','-dc',path],stdout=subprocess.PIPE,check=True).stdout; A=np.frombuffer(T,np.uint8); nl=np.flatnonzero(A==10); print('lines',len(nl)-1,flush=True); assert len(nl)-1==cat.N
# prefix identity of every row with the catalog is checked on a sample
def job(rows):
    out=[]; RF.set_df(n)
    for r in rows:
        r=int(r); sim=pfz.RangeSim(path); t=time.time()
        tr=T[int(nl[r])+1:int(nl[r+1])].split(b' '); tb,ts=float(tr[5]),float(tr[6]); tz=tb/ts
        try:
            o,info=RF.fetch_targets(sim,cat,[r],L,E); s=o[r]
            if s.get('failed'): rec=dict(st='failed')
            else: rec=dict(st='ok',beta=s['beta'],se=s['se'],mlogp=s['mlogp'],z=s['z'],z_lo=s['z_lo'],z_hi=s['z_hi'],w=s['z_width'],nu=s['n_unknown'],nf=s['n_feasible_values'],ds=s['digit_status'],tries=s['tries'],beta_exact=s['beta']==tb,se_exact=s['se']==ts,covered=s['z_lo']-1e-9<=tz<=s['z_hi']+1e-9)
        except Exception as e: rec=dict(st='err:'+repr(e)[:60])
        rec.update(row=r,tz=tz,bytes=sim.bytes,req=sim.requests,sec=time.time()-t,prefix_ok=b' '.join(tr[:5])+b' '==RF.prefix(cat,r)); out.append(rec); sim.close()
    return out
if __name__=='__main__':
    rng=np.random.default_rng(seed); rows=np.sort(rng.integers(3000,cat.N-3000,NR)); parts=np.array_split(rows,6)
    with Pool(6) as P: res=P.map(job,parts,1)
    d=pd.DataFrame([x for r in res for x in r]); d.to_csv(O+f'work/sim_{TAG}.csv',index=False)
    print(d.st.value_counts().to_dict()); k=d[d.st=='ok']
    print('n',len(d),'ok',len(k),'prefix ok',d.prefix_ok.mean(),'mean KiB',round(d.bytes.mean()/1024,1),'p50',round(d.bytes.median()/1024,1),'p90',round(d.bytes.quantile(.9)/1024,1),'max',round(d.bytes.max()/1024,1),'req',round(d.req.mean(),2),'sec',round(d.sec.mean(),2))
    print('digit status',k.ds.value_counts().to_dict()); print('n_unknown mean',round(k.nu.mean(),2),'tries',k.tries.value_counts().to_dict())
    print('z covered',k.covered.mean(),'beta exact',k.beta_exact.mean(),'se exact',k.se_exact.mean(),'|z-tz| max',float((k.z-k.tz).abs().max()),'p99',float((k.z-k.tz).abs().quantile(.99)),'width p50/p90/p99/max',k.w.quantile([.5,.9,.99,1]).round(4).tolist())
    print(k[~k.covered].head(10).to_string())
