from common import *
P=pd.read_csv(WORK+'plan_pairs.csv.gz',keep_default_na=False,low_memory=False,dtype={'idp':str}); Q=P[P.pre_status.eq('')]
qg=set(Q.gid[Q.img_src.eq('qsm_tar')]); a=sorted(set(Q.gid)-qg); b=sorted(qg); open(WORK+'groups_pass1.txt','w').write('\n'.join(a)+'\n'); open(WORK+'groups_pass2.txt','w').write('\n'.join(b)+'\n'); print(len(a),len(b))
