"""Shared paths and helpers for the all-target SuSiE-coloc run."""
import os, re, json, hashlib
import numpy as np, pandas as pd
ROOT='/public/home/mengxl/dzy/pd-v3/'
OUT=ROOT+'papers/01_mr_imaging_scan/results/all_target_susie_coloc_20261119/'
WORK=OUT+'work/'
REG=ROOT+'data/all_target_followup_20261118/regional_all/'
IM=REG+'imaging/'
RES=ROOT+'papers/01_mr_imaging_scan/results/'
RM=ROOT+'papers/03_world_model_integration/results/research_20261118/mr_imaging/'
BED='/public/home/mengxl/dzy/tools/magma/ref/EUR'
NIND=503; NB=(NIND+3)//4
HALF=500000; FULL=1000000
def h12(*x): return hashlib.sha1('|'.join(map(str,x)).encode()).hexdigest()[:12]
def snv_key(ch,pos,a,b):
    if len(a)!=1 or len(b)!=1 or a not in 'ACGT' or b not in 'ACGT' or a==b: return None
    return f'{ch}:{int(pos)}:'+'_'.join(sorted((a,b)))
def atomic_csv(df,path,**kw):
    tmp=path+f'.tmp{os.getpid()}'; df.to_csv(tmp,index=False,**kw); os.replace(tmp,path)
