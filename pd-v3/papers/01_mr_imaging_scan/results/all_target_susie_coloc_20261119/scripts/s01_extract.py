"""Standardised imaging windows (+/-1 Mb) from whole BIG40 files on disk and from the QSM tarball on disk. Nothing is downloaded."""
from common import *
import sys, gzip, tarfile, io, time
from multiprocessing import Pool
P=pd.read_csv(WORK+'plan_pairs.csv.gz',keep_default_na=False,low_memory=False,dtype={'idp':str}); P=P[P.pre_status.eq('')]
OI=WORK+'imaging/'; os.makedirs(OI,exist_ok=True)
CH={**{str(i):i for i in range(1,23)},**{f'{i:02d}':i for i in range(1,23)}}
def win_path(fam,idp,ch,c): return f'{OI}{fam}_{idp}/chr{ch}_{c}.csv.gz'
def emit(d,fam,idp,wins,src):
    # d: chrom(int),pos,a1,a2,beta,se,lp ; a2 is the effect allele in BIG40/QSM release files
    os.makedirs(f'{OI}{fam}_{idp}',exist_ok=True); n=0
    for ch,c in wins:
        w=d[(d.chrom==ch)&(abs(d.pos-c)<=FULL)]
        o=pd.DataFrame(dict(chromosome=ch,base_pair_location=w.pos.values,effect_allele=w.a2.astype(str).values,other_allele=w.a1.astype(str).values,beta=w.beta.values,standard_error=w.se.values,p_value=10.0**(-w.lp.values.astype(float))))
        o=o[(o.effect_allele.str.len()==1)&(o.other_allele.str.len()==1)]
        atomic_csv(o,win_path(fam,idp,ch,c),compression='gzip'); n+=len(o)
    return n
def whole(a):
    fam,idp,path,wins=a; wins=[w for w in wins if not os.path.exists(win_path(fam,idp,*w))]
    if not wins: return idp,'cached',0
    t=time.time()
    try:
        d=pd.read_csv(path,sep=' ',header=0,names=['chr','rsid','pos','a1','a2','beta','se','lp'],usecols=['chr','pos','a1','a2','beta','se','lp'],dtype={'chr':'category','pos':'int32','a1':'category','a2':'category','beta':'float32','se':'float32','lp':'float32'})
        d['chrom']=d['chr'].astype(str).map(CH).fillna(0).astype('int8').values
        n=emit(d,fam,idp,wins,path); return idp,'ok',n,round(time.time()-t)
    except Exception as e:
        os.makedirs(WORK+'imaging_errors',exist_ok=True); open(f'{WORK}imaging_errors/{fam}_{idp}.txt','w').write(repr(e)[:500]); return idp,'error',repr(e)[:120]
def tar_pass(req):
    with tarfile.open(TAR,'r|gz') as tf:
        seen=0
        for m in tf:
            b=os.path.basename(m.name); idp=b.split('.')[0]; seen+=1
            if seen<=5 or seen%50==0: print('tar member',m.name,m.size,flush=True)
            if idp not in req or not m.isfile() or not b.endswith('.txt') or not re.fullmatch(r'\d{4}',idp) or 'QSM' not in m.name.upper() and False: continue
            data=tf.fileobj.read(m.size-1); nl=data.index(b'\n'); head=data[:nl].decode().strip(); print('member',m.name,'header',head[:100],flush=True)
            cols=head.split(); d=pd.read_csv(io.BytesIO(data[nl+1:]),sep=' ',header=None,names=cols,low_memory=False); del data
            d=d.rename(columns={cols[0]:'chr'}); ren={c:'lp' for c in d.columns if c.lower().startswith('pval')}; d=d.rename(columns=ren)
            d['chrom']=d['chr'].astype(str).map(CH).fillna(0).astype(int); print('rows',len(d),'n',emit(d,'OxfordQSM2022',idp,req[idp],m.name),flush=True)
TAR=ROOT+'data/oxford_qsm/download/QSMandT2s_IDPs.tar.gz'
if __name__=='__main__':
    mode=sys.argv[1]
    w=P[P.img_src.eq('whole')].groupby(['family','idp','img_path']).apply(lambda g:sorted({(int(a),int(b)) for a,b in zip(g.chrom,g.pos37)})).reset_index(name='w')
    if mode=='whole':
        args=[(r.family,r.idp,r.img_path,r.w) for r in w.itertuples()]; print('whole idps',len(args),flush=True)
        with Pool(int(sys.argv[2])) as X:
            for i,r in enumerate(X.imap_unordered(whole,args)): print(i,r,flush=True)
    else:
        q=P[P.img_src.eq('qsm_tar')]; req={k:sorted({(int(a),int(b)) for a,b in zip(g.chrom,g.pos37)}) for k,g in q.groupby('idp')}; print('qsm idps',len(req),flush=True); tar_pass(req)
    print('DONE',flush=True)
