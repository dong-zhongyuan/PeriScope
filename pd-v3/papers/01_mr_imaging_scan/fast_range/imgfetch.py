"""Random-access extraction of target rows from a BIG40 gzip file through byte ranges only."""
import re,time,numpy as np
from scipy import special
import pfz
NUM=re.compile(rb'^-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:e[-+]?[0-9]+)?$')
NUMPOS=re.compile(rb'^(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:e[-+]?[0-9]+)?$')
CANDS=[ord(c) for c in '0123456789.-e+ ']
def nlp(z):return -(special.log_ndtr(-abs(z))+np.log(2))/np.log(10)
def hulp(x):
    m=re.match(rb'^-?[0-9]+(?:\.([0-9]+))?(?:e([-+]?[0-9]+))?$',x);return 0.5*10.0**(int(m.group(2) or 0)-len(m.group(1) or b''))
def row_ok(seg):
    p=seg.split(b' ')
    if len(p)!=3 or not NUM.match(p[0]) or not NUMPOS.match(p[1]) or not NUMPOS.match(p[2]):return False
    try:
        b,s,l=float(p[0]),float(p[1]),float(p[2])
        if s<=0 or not np.isfinite([b,s,l]).all():return False
        z=b/s;pred=nlp(z);dz=hulp(p[0])/s+abs(z)*hulp(p[1])/s
        bound=abs(nlp(abs(z)+dz)-nlp(max(abs(z)-dz,0)))+hulp(p[2])
        # BIG40 P values come from a t statistic; the normal formula drifts by about z^4/(4*df*ln10) at strong signals
        return bool(abs(pred-l)<=3*bound+0.01+0.001*l+1.5*z**4/(4*20000*2.302585))
    except (OverflowError,ZeroDivisionError,ValueError):return False
class Catalog:
    def __init__(self,text,keys,ref_idx):
        self.T=np.frombuffer(text,np.uint8);self.start=keys['start'].astype(np.uint32);self.N=keys['rows']
        self.key=(keys['chrom'].astype(np.uint64)<<32)|keys['pos']
        sp=np.flatnonzero(self.T==32).reshape(-1,7)[1:];self.pend=(sp[:,4]+1).astype(np.uint32)
        self.nls=ref_idx['nls'].astype(np.float64);self.cbyte=ref_idx['bits']/8.0;self.size=float(ref_idx['size'])
    def rows(self,chrom,a,b):
        return int(np.searchsorted(self.key,(chrom<<32)|a,'left')),int(np.searchsorted(self.key,(chrom<<32)|b,'right'))-1
    def comp(self,r,size):  # predicted compressed byte of data row r, scaled to the target file size
        return float(np.interp(r+1,self.nls,self.cbyte))*size/self.size
    def row_at(self,c,size):
        return int(np.interp(c*self.size/size,self.cbyte,self.nls))-1
class Failure(Exception):pass
def solve_unknown(S2,rowstart,ns,tfirst,tlast,last,minsupport=30,rows_sample=120,max_idx=60):
    """Resolve placeholders that occur inside target rows [tfirst,tlast). Evidence: target rows and aligned rows after them."""
    lo=int(ns[tfirst]);hi=int(rowstart[tlast])
    tu=S2[lo:hi];tu=tu[tu>=256]
    if len(tu)==0:return {},dict(unknown_idx=0,unknown_symbols=0,ambiguous=0,low_support=0,min_count=0)
    ids=np.unique(tu).astype(np.int64)-256
    if len(ids)>max_idx:return {},dict(unknown_idx=len(ids),unknown_symbols=len(tu),ambiguous=len(ids),low_support=len(ids),min_count=0,too_many=1)
    hi2=int(rowstart[last]);unk=np.flatnonzero(S2[lo:hi2]>=256)+lo
    u=S2[unk].astype(np.int64)-256;keep=np.isin(u,ids)
    rid=np.searchsorted(rowstart[tfirst:last+1],unk,'right')-1+tfirst
    # rows that contain an unknown outside the solved set cannot serve as evidence
    badrows=set(np.unique(rid[~keep]).tolist())
    unk=unk[keep];u=u[keep];rid=rid[keep]
    counts=np.array([(u==v).sum() for v in ids.tolist()]);assign={}
    def text(i,a):
        p,q=int(ns[i]),int(rowstart[i+1])-1;s=S2[p:q].copy();m=np.flatnonzero(s>=256)
        for j in m.tolist():
            c=a.get(int(s[j])-256)
            if c is None:return None
            s[j]=c
        return s.astype(np.uint8).tobytes()
    isd=lambda x:(x>=48)&(x<=57)
    occ={v:np.flatnonzero(u==v) for v in ids.tolist()}
    for v in ids.tolist():
        pos=unk[occ[v]][:2000];pv=S2[pos-1];nx=S2[np.minimum(pos+1,len(S2)-1)]
        st=(pv==32)|(pv==10)|(pv>=256)
        if (isd(pv)&(isd(nx)|(nx>=256))).mean()>0.5:g=46
        elif (st&(nx==46)).mean()>0.3:g=48
        elif (st&isd(nx)).mean()>0.3:g=45
        elif (nx==46).mean()>0.3:g=48
        else:g=46
        assign[v]=g
    rng=np.random.default_rng(0);order=ids[np.argsort(-counts)].tolist();samples={}
    for v in order:
        r=np.unique(rid[occ[v]]);r=np.array([x for x in r.tolist() if x not in badrows and x<last],dtype=np.int64)
        samples[v]=r if len(r)<=rows_sample else rng.choice(r,rows_sample,replace=False)
    amb=set()
    for it in range(4):
        changed=False;amb=set()
        for v in order:
            res={}
            for c in CANDS:
                a=dict(assign);a[v]=c;bad=0
                for i in samples[v].tolist():
                    t=text(i,a)
                    if t is None or not row_ok(t):bad+=1
                res[c]=bad
            best=min(res.values());win=[c for c,b in res.items() if b==best]
            if assign[v] not in win:assign[v]=win[0];changed=True
            if len(win)>1 or best>0 or len(samples[v])==0:amb.add(v)
        if not changed:break
    return assign,dict(unknown_idx=len(ids),unknown_symbols=len(tu),ambiguous=len(amb),min_count=int(counts.min()),low_support=int((counts<minsupport).sum()))
def fetch_rows(sim,cat,r0,r1,L,E,extra_rows=1500,step=1<<18,max_total=64<<20,verbose=False):
    """Return (row bytes r0..r1 inclusive, info). Row text ends with \\n."""
    size=sim.size;r_end=min(cat.N-1,r1+extra_rows)
    a=int(cat.comp(r0,size))-L-E;b=int(cat.comp(min(r1+1,cat.N-1),size))+E
    a=max(0,a);b=min(size,b);D=sim.read(a,b);info=dict(reads=1,widen=0);lead=L
    t0=time.time()
    for attempt in range(12):
        Rr=pfz.decode_range(D)
        if Rr['stop'] in(0,-1,2):raise Failure('decode stop=%d entry=%d'%(Rr['stop'],Rr['entry_bit']))
        S=Rr['out'];entry=a+Rr['entry_bit']/8
        ra=max(0,cat.row_at(entry,size)-3000)
        rb=min(cat.N,r_end+2)
        st=cat.start[ra:rb];en=cat.pend[ra:rb]
        A=pfz.align(S,cat.T,st,en)
        info.update(entry_byte=entry,sym=len(S),aligned=A['aligned'],conflicts=A['conflicts'],attempt=attempt)
        if A['conflicts']:raise Failure('alignment conflicts %d'%A['conflicts'])
        rs=A['rowstart'];tfirst=r0-ra;tlast=r1+1-ra;good=np.flatnonzero(rs>=0)
        if len(good)==0:act='both'
        else:
            seed=int(good[0]);lastgood=int(good[-1]);span=rs[tfirst:tlast+1]
            if (span>=0).all():
                act='ok';last=tlast
                while last+1<len(rs) and rs[last+1]>=0:last+=1
            elif seed>tfirst:act='back'
            elif lastgood<tlast:act='fwd'
            else:act='moreback'
        info['segments']=A['broke']
        if act=='ok':
            W=A['W'];S2=S.copy();ph=np.flatnonzero(S2>=256);wv=W[S2[ph].astype(np.int64)-256];S2[ph[wv>=0]]=wv[wv>=0].astype(np.uint16)
            plen=(en-st).astype(np.int64);ns=rs+plen  # numeric section starts for aligned rows
            asg,diag=solve_unknown(S2,rs,ns,tfirst,tlast,last)
            info.update(diag)
            if diag['ambiguous'] or diag.get('low_support',0)>0:
                info['widen']+=1;act='moreback'
            else:
                for v,c in asg.items():S2[S2==256+v]=c
                info['solved_chars']=''.join(sorted(chr(c) for c in asg.values()))
                lo=int(rs[tfirst]);hi=int(rs[tlast])
                seg=S2[lo:hi]
                if (seg>=256).any():raise Failure('unresolved symbols in target')
                bad=0;p0=ns[tfirst:tlast];q0=rs[tfirst+1:tlast+1]-1
                for p,q in zip(p0.tolist(),q0.tolist()):
                    if not row_ok(S2[p:q].astype(np.uint8).tobytes()):bad+=1
                info.update(identity_bad=bad,lead_used=a,seconds=round(time.time()-t0,3),last_step=attempt)
                if bad:raise Failure('identity failures %d'%bad)
                return seg.astype(np.uint8).tobytes(),info
        total=sim.bytes
        if total>max_total:raise Failure('budget exceeded')
        if act=='fwd':
            n=sim.read(b,b+step);b+=len(n);D+=n;info['reads']+=1;continue
        na=max(0,a-max(step,lead));n=sim.read(na,a)
        if not n and act!='both':raise Failure('reached file start')
        lead=max(step,lead)*2 if act=='moreback' else lead
        D=n+D;a=na;info['reads']+=1
        if act=='both':
            n=sim.read(b,b+step);b+=len(n);D+=n;info['reads']+=1
    raise Failure('attempts exhausted '+str(info))
