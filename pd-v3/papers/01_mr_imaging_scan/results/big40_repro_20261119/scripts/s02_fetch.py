"""Step 2: byte-range reads of the frozen target rows from BIG40 release2 repro (and, for palindromic hit pairs only, the native discovery file stats/).
One HTTPS connection per worker process. Resumable: finished rows are skipped. Hard cap on bytes."""
import sys, os, json, time, glob, random
import numpy as np, pandas as pd
from multiprocessing import Pool, Value
O='/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan/results/big40_repro_20261119/'; sys.path.insert(0,O+'scripts')
import rfetch as RF; F=RF.F
TARGETS=sys.argv[1]; SUB=sys.argv[2]; NW=int(sys.argv[3]); CAP=int(float(sys.argv[4])*(1<<30)); L0=128<<10; E0=32<<10; GAP=96<<10
OUT=O+f'work/{SUB}/'; os.makedirs(OUT,exist_ok=True)
cat=RF.catalog(); total=Value('q',0); conn=None; fh=None
def init(tv):
    global total,conn,fh
    total=tv; conn=RF.Conn(); fh=open(OUT+f'part_{os.getpid()}.jsonl','a',buffering=1)
def clusters(rows,size):
    rows=sorted(rows); out=[[rows[0]]]
    for r in rows[1:]:
        if (cat.comp(r,size)-cat.comp(out[-1][-1],size))<GAP and r-out[-1][0]<6000: out[-1].append(r)
        else: out.append([r])
    return out
def emit(rec): fh.write(json.dumps(rec)+'\n')
def work(task):
    rel,idp,n,rows=task; t0=time.time(); base=dict(release=rel,idp=idp); got=0
    RF.set_df(n if n==n else 11000.0)
    with total.get_lock(): over=total.value>=CAP
    if over:
        for r in rows: emit({**base,'row':int(r),'status':'cap_reached'})
        return dict(idp=idp,n=len(rows),ok=0,bytes=0)
    try: rem=RF.Remote(conn,rel,idp)
    except FileNotFoundError:
        for r in rows: emit({**base,'row':int(r),'status':'no_remote_file'})
        return dict(idp=idp,n=len(rows),ok=0,bytes=0)
    except Exception as e:
        for r in rows: emit({**base,'row':int(r),'status':'error_head','error':repr(e)[:120]})
        return dict(idp=idp,n=len(rows),ok=0,bytes=0)
    okc=0
    for cl in clusters(rows,rem.size):
        b0=rem.bytes; res={}; note=None
        for (Lx,Ex,tries) in [(L0,E0,3),(3*L0,E0,3)]:
            todo=[r for r in cl if r not in res or res[r].get('failed')]
            if not todo: break
            try:
                o,info=RF.fetch_targets(rem,cat,todo,L=Lx,E=Ex,max_tries=tries)
                for r,v in o.items(): res[r]=v
            except F.Failure as e: note='failure:'+str(e)[:80]
            except Exception as e: note='error:'+repr(e)[:100]
        for r in [r for r in cl if r not in res or res[r].get('failed')]:   # last resort: frozen production reader with its own widening, at most 4 MiB
            try:
                seg,info=F.fetch_rows(rem,cat,r,r,L0,E0,step=1<<17,max_total=rem.bytes+(4<<20)); f=seg.split(b'\n')[0].split(b' '); b_,s_,m_=float(f[5]),float(f[6]),float(f[7]); res[r]=dict(beta=b_,se=s_,mlogp=m_,z=b_/s_,z_lo=b_/s_,z_hi=b_/s_,n_feasible_values=1,beta_se_unique=True,n_unknown=int(info.get('unknown_symbols',0)),z_width=0.0,digit_status='production_reader_exact',tries=int(info.get('reads',0)))
            except Exception as e: note='fallback:'+repr(e)[:80]
        used=rem.bytes-b0; got+=used
        with total.get_lock(): total.value+=used
        for r in cl:
            v=res.get(r,{})
            if v and not v.get('failed'):
                okc+=1; emit({**base,'row':int(r),'status':'ok',**{k:(bool(x) if isinstance(x,(np.bool_,)) else x) for k,x in v.items() if k not in('row_text_ok',)},'cluster_size':len(cl),'cluster_bytes':used,'file_size':rem.size,'etag':rem.etag})
            else: emit({**base,'row':int(r),'status':'failed','note':note,'cluster_size':len(cl),'cluster_bytes':used,'file_size':rem.size})
    return dict(idp=idp,n=len(rows),ok=okc,bytes=got,sec=time.time()-t0,req=rem.requests,retries=conn.retries)
def done_set():
    D={}
    for f in glob.glob(OUT+'part_*.jsonl'):
        for l in open(f):
            try: r=json.loads(l)
            except Exception: continue
            k=(r['release'],r['idp'],r['row']); s=r['status']
            if s in('ok','no_remote_file') or k not in D: D[k]=s
            if s=='ok': D[k]='ok'
    return D
if __name__=='__main__':
    T=pd.read_csv(TARGETS,dtype={'idp':str}); T['idp']=T.idp.str.zfill(4)
    D=done_set(); T['k']=list(zip(T.release,T.idp,T.catalog_row)); T['st']=T.k.map(D); ntot=len(T); T=T[~T.st.isin(['ok','no_remote_file'])]; print('targets',ntot,'to do',len(T),flush=True)
    random.seed(20261119); tasks=[]
    for (rel,idp,pr),g in T.groupby(['release','idp','prio']): tasks.append((pr,random.random(),rel,idp,float(g.n_rep.iloc[0]),sorted(g.catalog_row.astype(int).tolist())))
    tasks.sort(); print('tasks',len(tasks),flush=True)
    t0=time.time(); ok=0; nr=0; b=0
    with Pool(NW,initializer=init,initargs=(total,)) as P:
        for i,r in enumerate(P.imap_unordered(work,[t[2:] for t in tasks],1)):
            ok+=r['ok']; nr+=r['n']; b+=r['bytes']
            if (i+1)%50==0 or i==len(tasks)-1: print(f'{i+1}/{len(tasks)} tasks rows {nr} ok {ok} fail {nr-ok} GiB {b/2**30:.3f} sec {time.time()-t0:.0f} MiB/s {b/2**20/(time.time()-t0+1e-9):.2f}',flush=True)
    print('DONE rows',nr,'ok',ok,'bytes',b,flush=True)
