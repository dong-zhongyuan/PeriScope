"""ctypes wrapper and simulated range reader for libpfz."""
import ctypes as C, os, subprocess
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
L=C.CDLL(str(HERE/'libpfz.so'))
u64p=np.ctypeslib.ndpointer(np.uint64,flags='C');u8p=np.ctypeslib.ndpointer(np.uint8,flags='C');u16p=np.ctypeslib.ndpointer(np.uint16,flags='C');u32p=np.ctypeslib.ndpointer(np.uint32,flags='C')
L.pf_decode_range.restype=C.c_long
L.pf_index_file.restype=C.c_long
L.pf_keyscan.restype=C.c_long
L.pf_tsvize.restype=C.c_uint64
TEXT=np.zeros(256,np.uint8);TEXT[[9,10]]=1;TEXT[32:127]=1
WSIZE=32768
def index_file(path,maxb=200000):
    bits=np.zeros(maxb,np.uint64);outs=np.zeros(maxb,np.uint64);nls=np.zeros(maxb,np.uint64);ty=np.zeros(maxb,np.uint8)
    tot=C.c_uint64();end=C.c_uint64();fs=C.c_uint64();st=C.c_int()
    n=L.pf_index_file(str(path).encode(),bits.ctypes.data_as(C.c_void_p),outs.ctypes.data_as(C.c_void_p),nls.ctypes.data_as(C.c_void_p),ty.ctypes.data_as(C.c_void_p),C.c_long(maxb),C.byref(tot),C.byref(end),C.byref(fs),C.byref(st))
    assert n<=maxb
    return dict(bits=bits[:n],outs=outs[:n],nls=nls[:n],types=ty[:n],total_out=tot.value,end_bit=end.value,size=fs.value,status=st.value)
def decode_range(data,from_bit=0,to_bit=None,cap=None,need=3,allow=TEXT):
    """data: bytes for one contiguous compressed range. returns dict with out (uint16, no window prefix)."""
    n=len(data);buf=np.zeros(n+32,np.uint8);buf[:n]=np.frombuffer(data,np.uint8)
    cap=cap or int(n*12+1<<16) if False else (cap or (n*12+(1<<17)))
    out=np.empty(cap,np.uint16);blk=np.zeros(2*20000,np.uint64)
    eb=C.c_uint64();ee=C.c_uint64();st=C.c_int();nb=C.c_long();tested=C.c_uint64()
    if to_bit is None:to_bit=n*8
    k=L.pf_decode_range(buf.ctypes.data_as(C.c_void_p),C.c_uint64(n),C.c_uint64(from_bit),C.c_uint64(to_bit),out.ctypes.data_as(C.c_void_p),C.c_uint64(cap),allow.ctypes.data_as(C.c_void_p),C.c_int(need),C.byref(eb),C.byref(ee),C.byref(st),blk.ctypes.data_as(C.c_void_p),C.c_long(20000),C.byref(nb),C.byref(tested))
    return dict(out=out[WSIZE:WSIZE+k],entry_bit=eb.value,end_bit=ee.value,stop=st.value,blocks=blk[:2*nb.value].reshape(-1,2),tested=tested.value)
def keyscan(text,family,has_header,maxr):
    """text: bytes-like (numpy uint8 array or bytes)."""
    arr=np.frombuffer(text,np.uint8) if not isinstance(text,np.ndarray) else text
    ch=np.zeros(maxr,np.uint8);pos=np.zeros(maxr,np.uint32);h=np.zeros(maxr,np.uint64);st=np.zeros(maxr,np.uint64)
    d=C.c_uint64();d2=C.c_uint64();cc=np.zeros(32,np.uint64);un=C.c_long();bad=C.c_long()
    r=L.pf_keyscan(arr.ctypes.data_as(C.c_void_p),C.c_uint64(len(arr)),C.c_int(family),C.c_int(int(has_header)),ch.ctypes.data_as(C.c_void_p),pos.ctypes.data_as(C.c_void_p),h.ctypes.data_as(C.c_void_p),st.ctypes.data_as(C.c_void_p),C.c_long(maxr),C.byref(d),C.byref(d2),cc.ctypes.data_as(C.c_void_p),C.byref(un),C.byref(bad))
    return dict(rows=r,chrom=ch[:min(r,maxr)],pos=pos[:min(r,maxr)],hash=h[:min(r,maxr)],start=st[:min(r,maxr)],digest=d.value,digest2=d2.value,chrcnt=cc.copy(),unsorted=un.value,bad=bad.value)
def gunzip(path):
    return subprocess.run(['gzip','-dc',str(path)],stdout=subprocess.PIPE,check=True).stdout
def tsvize(text):
    arr=np.frombuffer(text,np.uint8);out=np.empty(len(arr)+16,np.uint8);c=C.c_uint64()
    n=L.pf_tsvize(arr.ctypes.data_as(C.c_void_p),C.c_uint64(len(arr)),out.ctypes.data_as(C.c_void_p),C.byref(c))
    return out[:n],c.value
class RangeSim:
    """Simulated HTTP byte-range server backed by a local file; every read is accounted."""
    def __init__(self,path):
        self.fd=os.open(path,os.O_RDONLY);self.size=os.fstat(self.fd).st_size;self.requests=0;self.bytes=0;self.log=[]
    def read(self,a,b):
        a=max(0,int(a));b=min(self.size,int(b));d=os.pread(self.fd,b-a,a);self.requests+=1;self.bytes+=len(d);self.log.append((a,b));return d
    def close(self):os.close(self.fd)
class MemRange(RangeSim):
    def __init__(self,arr):self.arr=arr;self.size=len(arr);self.requests=0;self.bytes=0;self.log=[]
    def read(self,a,b):
        a=max(0,int(a));b=min(self.size,int(b));d=self.arr[a:b].tobytes();self.requests+=1;self.bytes+=len(d);self.log.append((a,b));return d
    def close(self):pass

L.pf_align.restype=C.c_long
def align(S,blob,st,en,step=500,minnum=6,maxnum=64,W=None):
    """S uint16 symbols; blob uint8; st,en uint32 row prefix bounds. Returns dict."""
    m=len(st);W=(-np.ones(32768,np.int16)) if W is None else W;rs=np.empty(m,np.int64);seed=C.c_long();conf=C.c_long();brk=C.c_long()
    n=L.pf_align(S.ctypes.data_as(C.c_void_p),C.c_long(len(S)),blob.ctypes.data_as(C.c_void_p),st.ctypes.data_as(C.c_void_p),en.ctypes.data_as(C.c_void_p),C.c_long(m),C.c_long(step),C.c_int(minnum),C.c_int(maxnum),W.ctypes.data_as(C.c_void_p),rs.ctypes.data_as(C.c_void_p),C.byref(seed),C.byref(conf),C.byref(brk))
    return dict(aligned=n,seed=seed.value,conflicts=conf.value,broke=brk.value,rowstart=rs,W=W)
