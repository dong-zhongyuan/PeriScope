"""Offline tuning/verification: range-read random rows out of complete local BIG40 files through the same code path and compare with the true decompressed rows."""
import sys, time, subprocess, numpy as np, pandas as pd
O='/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan/results/big40_repro_20261119/'; sys.path.insert(0,O+'scripts')
import rfetch as RF, pfz; F=RF.F
from multiprocessing import Pool
B='/public/home/mengxl/dzy/pd-v3/data/all_target_followup_20261118/regional_all/imaging/BIG40/'
FILES=sys.argv[1].split(','); NR=int(sys.argv[2]); SETS=[tuple(int(x) for x in s.split(':')) for s in sys.argv[3].split(',')]; TAG=sys.argv[4] if len(sys.argv)>4 else 'x'
cat=RF.catalog()
def run(a):
    idp,L,E,ms,X,rows=a; L*=1024; E*=1024; RF.MS[0]=ms; p=B+idp+'/'+idp+'.txt.gz'; RF.set_df(21000)
    T=subprocess.run(['gzip','-dc',p],stdout=subprocess.PIPE,check=True).stdout; A=np.frombuffer(T,np.uint8); nl=np.flatnonzero(A==10)
    out=[]
    for r in rows:
        r=int(r); sim=pfz.RangeSim(p); t=time.time(); r0=max(0,r-X//2); r1=min(cat.N-1,r+X//2)
        try:
            seg,info=F.fetch_rows(sim,cat,r0,r1,L,E,step=1<<17); st='ok'; ls=seg.split(b'\n'); got=ls[r-r0]+b'\n'
            truth=T[int(nl[r])+1:int(nl[r+1])+1]; same=got==truth; allsame=seg==T[int(nl[r0])+1:int(nl[r1+1])+1]
        except Exception as e:
            st='fail:'+repr(e)[:50]; same=False; allsame=False; info={}
        out.append(dict(idp=idp,L=L,E=E,ms=ms,X=X,row=r,st=st,same=same,allsame=allsame,bytes=sim.bytes,req=sim.requests,sec=time.time()-t,unk=info.get('unknown_symbols',-1),uidx=info.get('unknown_idx',-1),wid=info.get('widen',-1))); sim.close()
    return out
if __name__=='__main__':
    rng=np.random.default_rng(1); tasks=[]
    rows=np.sort(rng.integers(5000,cat.N-5000,NR))
    for idp in FILES:
        for L,E,ms,X in SETS: tasks.append((idp,L,E,ms,X,rows))
    with Pool(min(6,len(tasks))) as P: res=P.map(run,tasks,1)
    d=pd.DataFrame([x for r in res for x in r]); d.to_csv(O+f'work/sim_tune_{TAG}.csv',index=False)
    g=d.groupby(['L','E','ms','X']).agg(n=('row','size'),ok=('st',lambda s:(s=='ok').mean()),same=('same','mean'),allsame=('allsame','mean'),mean_KiB=('bytes',lambda s:s.mean()/1024),p90_KiB=('bytes',lambda s:s.quantile(.9)/1024),max_KiB=('bytes',lambda s:s.max()/1024),req=('req','mean'),sec=('sec','mean'),unk=('unk','mean'),uidx=('uidx','mean')).round(3); print(g.to_string()); print(d[d.st!='ok'].st.value_counts().head(10).to_string())
