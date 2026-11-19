"""Byte-range access to BIG40 per-IDP gzip files, single rows or tight row clusters.
Re-uses (imports, does not modify) the frozen production decoder papers/01_mr_imaging_scan/fast_range/{pfz.py,libpfz.so,imgfetch.py} and its shared variant catalog.
Only the numeric row check is replaced, because the production check assumes the discovery sample size; the replication files have about half that sample."""
import sys, os, re, time, ssl, random, http.client
import numpy as np
from scipy import special
MR='/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan/'
sys.path.insert(0,MR); sys.path.insert(0,MR+'fast_range')
import pfz, imgfetch as F, fast_range_common as C
HOST='open.oxcin.ox.ac.uk'; BASE='/ukbiobank/big40/release2/'; UA='research-client'
NUM=F.NUM; NUMPOS=F.NUMPOS; hulp=F.hulp
DF=[8000.0,13000.0]        # replaced per file by set_df
def set_df(n):
    DF[0]=0.70*float(n); DF[1]=1.02*float(n)
def _nlp_t(z,df):
    p=special.stdtr(df,-abs(z))
    return float(-np.log10(2*p)) if p>0 else float('inf')
def row_ok(seg):
    p=seg.split(b' ')
    if len(p)!=3 or not NUM.match(p[0]) or not NUMPOS.match(p[1]) or not NUMPOS.match(p[2]): return False
    try:
        b,s,l=float(p[0]),float(p[1]),float(p[2])
        if s<=0 or not np.isfinite([b,s,l]).all(): return False
        z=b/s; dz=hulp(p[0])/s+abs(z)*hulp(p[1])/s
        lo=_nlp_t(max(abs(z)-dz,0),DF[0]); hi=_nlp_t(abs(z)+dz,DF[1])
        tol=hulp(p[2])*2+0.01+0.001*l
        return bool(lo-tol<=l<=hi+tol)
    except (OverflowError,ZeroDivisionError,ValueError): return False
F.row_ok=row_ok
def catalog(): return C.imaging_catalog(F)
class Conn:
    """one keep-alive HTTPS connection per worker process"""
    def __init__(self,timeout=60): self.c=None; self.timeout=timeout; self.bytes=0; self.requests=0; self.retries=0
    def close(self):
        try:
            if self.c is not None: self.c.close()
        except Exception: pass
        self.c=None
    def req(self,method,path,headers):
        if self.c is None: self.c=http.client.HTTPSConnection(HOST,timeout=self.timeout,context=ssl.create_default_context())
        self.c.request(method,path,headers={'User-Agent':UA,'Accept-Encoding':'identity',**headers}); r=self.c.getresponse(); b=r.read(); self.requests+=1
        if r.will_close: self.close()
        return r.status,r.headers,b
    def retry(self,method,path,headers,check=None,attempts=12):
        last=None
        for n in range(attempts):
            try:
                st,h,b=self.req(method,path,headers)
                if st in (408,425,429) or st>=500: raise ConnectionError('HTTP %d'%st)
                if check: check(st,h,b)
                return st,h,b
            except (OSError,http.client.HTTPException,EOFError,ssl.SSLError) as e:
                last=e; self.close(); self.retries+=1; time.sleep(min(2**n,60)+random.uniform(0,2))
        raise RuntimeError('gave up: '+repr(last)[:120])
class Remote:
    """pinned range view of one remote file; interface matches what imgfetch.fetch_rows needs"""
    def __init__(self,conn,release,idp):
        self.conn=conn; self.path=f'{BASE}{release}/{idp}.txt.gz'
        st,h,_=conn.retry('HEAD',self.path,{})
        if st!=200: raise FileNotFoundError(f'HEAD {st} {self.path}')
        self.size=int(h['Content-Length']); self.etag=h.get('ETag'); self.lastmod=h.get('Last-Modified'); self.bytes=0; self.requests=0
    def read(self,a,b):
        a=max(0,int(a)); b=min(self.size,int(b))
        if b<=a: return b''
        want='bytes %d-%d/%d'%(a,b-1,self.size)
        def check(st,h,body):
            if st!=206 or h.get('Content-Range')!=want: raise RuntimeError(f'bad range response {st} {h.get("Content-Range")}')
            if self.etag and h.get('ETag')!=self.etag: raise RuntimeError('etag changed')
            if len(body)!=b-a: raise EOFError('short body')
        _,_,body=self.conn.retry('GET',self.path,{'Range':'bytes=%d-%d'%(a,b-1),'If-Range':self.etag},check)
        self.bytes+=len(body); self.requests+=1; self.conn.bytes+=len(body)
        return body
def prefix(cat,r):
    return bytes(cat.T[int(cat.start[r]):int(cat.pend[r])])
def parse_rows(seg):
    out=[]
    for l in seg.split(b'\n'):
        if not l: continue
        f=l.split(b' ')
        out.append((f[0].decode(),f[1].decode(),int(f[2]),f[3].decode(),f[4].decode(),float(f[5]),float(f[6]),float(f[7])))
    return out
def find_row(cat,chrom,pos,ea,oa):
    """catalog rows at chrom:pos with exactly these two alleles. returns (row or None, n_rows_at_position, a1, a2)"""
    r0,r1=cat.rows(int(chrom),int(pos),int(pos)); n=max(0,r1-r0+1)
    for r in range(r0,r1+1):
        f=prefix(cat,r).split(b' '); a1=f[3].decode(); a2=f[4].decode()
        if {a1,a2}=={ea,oa}: return r,n,a1,a2
    return None,n,None,None
# For a single target row there are few later rows sharing a placeholder. A placeholder is accepted only when exactly one candidate character passes the numeric check on every evidence row
# (the true character always passes); the production extra requirement of >=30 occurrences is relaxed to MS[0] so that one row does not force a 2 MiB lead.
MS=[3]; _orig_solve=F.solve_unknown
def _solve(S2,rs,ns,a,b,last,**k):
    k['minsupport']=MS[0]; return _orig_solve(S2,rs,ns,a,b,last,**k)
F.solve_unknown=_solve

# ---------------------------------------------------------------------------------------------------------------
# Single-row reader. A deflate entry leaves some characters as placeholders (references into the unseen 32 KiB before the entry).
# Placeholders in the variant columns are resolved from the catalog by pfz.align. Those left in the numeric columns are treated as unknowns with a feasible set:
# a candidate character survives only if it passes the numeric check on every row where that placeholder occurs. The result is an interval for z over ALL surviving
# assignments, so nothing is guessed: if a digit cannot be determined the reported interval says how much it could matter.
CANDS=list(F.CANDS)
def _tmpl(arr):
    b=bytearray(np.minimum(arr,255).astype(np.uint8).tobytes()); idx=np.flatnonzero(arr>=256)
    return b,[(int(j),int(arr[j])-256) for j in idx]
def _ok_with(b,ph,assign):
    for j,v in ph:
        c=assign.get(v)
        if c is None: return False
        b[j]=c
    return row_ok(bytes(b))
def _any_ok(b,ph,dom,fix):
    """is there an assignment of the row's other placeholders (within domains) that passes, given fix?"""
    ids=sorted({v for _,v in ph}); free=[v for v in ids if v not in fix]
    if not free: return _ok_with(b,ph,fix)
    n=1
    for v in free: n*=len(dom[v])
    if n>400: return True    # too many combinations to use this row as evidence
    a=dict(fix)
    def rec(i):
        if i==len(free): return _ok_with(b,ph,a)
        for c in dom[free[i]]:
            a[free[i]]=c
            if rec(i+1): return True
        return False
    return rec(0)
def solve_target(S2,rs,ns,t,first,last,max_rows=40,max_combo=20000):
    arr=S2[int(ns[t]):int(rs[t+1])-1]; b0,ph0=_tmpl(arr); U=sorted({v for _,v in ph0})
    if not U:
        s=bytes(b0); return dict(ok=row_ok(s),n_unknown=0,feasible=[s] if row_ok(s) else [])
    lo=int(rs[first]); hi=int(rs[last]); P=np.flatnonzero(S2[lo:hi]>=256)+lo; I=S2[P].astype(np.int64)-256
    starts=rs[first:last+1]; R=np.searchsorted(starts,P,'right')-1+first
    inU=np.isin(I,U); bad=set(R[~inU].tolist()); ev={v:[] for v in U}
    order=np.argsort(np.abs(R[inU]-t),kind='stable'); seen=set()
    for r_,v in zip(R[inU][order].tolist(),I[inU][order].tolist()):
        if r_==t or r_ in bad or len(ev[v])>=max_rows or (r_,v) in seen: continue
        seen.add((r_,v)); ev[v].append(r_)
    T={}
    def tm(r_):
        if r_ not in T: T[r_]=_tmpl(S2[int(ns[r_]):int(rs[r_+1])-1])
        return T[r_]
    dom={v:list(CANDS) for v in U}
    for it in range(3):
        ch=False
        for v in U:
            keep=[]
            for c in dom[v]:
                good=True
                for r_ in ev[v]:
                    b,ph=tm(r_)
                    if not _any_ok(b,ph,dom,{v:c}): good=False; break
                if good and not _any_ok(b0,ph0,dom,{v:c}): good=False
                if good: keep.append(c)
            if len(keep)!=len(dom[v]): dom[v]=keep; ch=True
            if not keep: return dict(ok=False,n_unknown=len(ph0),feasible=[],why='empty_domain')
        if not ch: break
    n=1
    for v in U: n*=len(dom[v])
    if n>max_combo: return dict(ok=False,n_unknown=len(ph0),feasible=[],why='too_many_combinations',n_combo=n)
    F_=[]; a={}
    def rec(i):
        if i==len(U):
            if _ok_with(b0,ph0,a): F_.append(bytes(b0))
            return
        for c in dom[U[i]]:
            a[U[i]]=c; rec(i+1)
    rec(0)
    return dict(ok=len(F_)>0,n_unknown=len(ph0),feasible=F_,why=None if F_ else 'no_feasible')
def summarise(fe):
    vals=sorted({tuple(float(x) for x in s.split(b' ')) for s in fe},key=lambda v:v[0]/v[1])
    zs=[v[0]/v[1] for v in vals]; m=vals[len(vals)//2]
    return dict(beta=m[0],se=m[1],mlogp=m[2],z=m[0]/m[1],z_lo=min(zs),z_hi=max(zs),n_feasible_values=len(vals),beta_se_unique=len({(v[0],v[1]) for v in vals})==1)
def fetch_targets(sim,cat,rows,L=192<<10,E=48<<10,step=1<<17,max_tries=4,max_width=0.10,extra_rows=1500):
    """rows: sorted catalog rows close together in the file. Returns {row:result}, info."""
    size=sim.size; rows=sorted(rows); a=max(0,int(cat.comp(rows[0],size))-L-E); b=min(size,int(cat.comp(min(rows[-1]+1,cat.N-1),size))+E)
    D=sim.read(a,b); out={}; pend=list(rows); info=dict(reads=1,tries=0)
    for attempt in range(max_tries):
        info['tries']=attempt+1
        Rr=pfz.decode_range(D)
        if Rr['stop'] in (0,-1,2): act='both'; S=None
        else:
            S=Rr['out']; entry=a+Rr['entry_bit']/8; ra=max(0,cat.row_at(entry,size)-3000); rb=min(cat.N,pend[-1]+extra_rows+2)
            st=cat.start[ra:rb]; en=cat.pend[ra:rb]; A=pfz.align(S,cat.T,st,en)
            if A['conflicts']: raise F.Failure('alignment conflicts')
            rs=A['rowstart']; W=A['W']; S2=S.copy(); ph=np.flatnonzero(S2>=256); wv=W[S2[ph].astype(np.int64)-256]; S2[ph[wv>=0]]=wv[wv>=0].astype(np.uint16)
            ns=rs+(en-st).astype(np.int64); still=[]; need_back=False; need_fwd=False
            good=np.flatnonzero(rs>=0)
            for r in pend:
                t=r-ra
                if len(good)==0 or t<0 or t+1>=len(rs): need_back=True; still.append(r); continue
                if rs[t]<0 or rs[t+1]<0:
                    if good[0]>t: need_back=True
                    elif good[-1]<t+1: need_fwd=True
                    else: need_back=True
                    still.append(r); continue
                f=t
                while f>0 and rs[f-1]>=0: f-=1
                l=t+1
                while l+1<len(rs) and rs[l+1]>=0: l+=1
                sol=solve_target(S2,rs,ns,t,f,l)
                if sol['ok']:
                    s=summarise(sol['feasible']); w=s['z_hi']-s['z_lo']
                    if w<=max_width or attempt==max_tries-1:
                        s.update(n_unknown=sol['n_unknown'],z_width=w,digit_status='exact' if sol['n_unknown']==0 else ('unknown_chars_no_effect_on_beta_se' if s['beta_se_unique'] else ('interval_narrow' if w<=max_width else 'interval_wide')),tries=attempt+1,row_text_ok=True); out[r]=s; continue
                need_back=True; still.append(r)
            pend=still; act='done' if not pend else ('fwd' if need_fwd and not need_back else 'back')
        if act=='done': break
        if attempt==max_tries-1: break
        if act in('fwd',):
            n=sim.read(b,b+step); b+=len(n); D+=n; info['reads']+=1; continue
        na=max(0,a-step*(attempt+1))
        if na==a and act=='back': break
        n=sim.read(na,a); D=n+D; a=na; info['reads']+=1
        if act=='both':
            n=sim.read(b,b+step); b+=len(n); D+=n; info['reads']+=1
    for r in pend: out.setdefault(r,dict(failed=True))
    return out,info
