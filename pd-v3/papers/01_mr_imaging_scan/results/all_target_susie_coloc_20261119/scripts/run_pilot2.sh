O=/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan/results/all_target_susie_coloc_20261119; cd $O/scripts
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
R=/public/home/mengxl/dzy/envs/rcoloc/bin/Rscript
/public/home/mengxl/dzy/pd-v3_runtime/bin/python - <<PY
import pandas as pd
s=pd.read_csv('$O/work/pilot/sample.csv'); r=pd.read_csv('/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/research_20261118/mr_imaging/susie_results.csv')
b=r[r.err.notna()].region_id; x=s[s.region_id.isin(b)]; open('$O/work/pilot/pairs_nc.txt','w').write('\n'.join(x.pair_id)+'\n'); open('$O/work/pilot/groups_nc.txt','w').write('\n'.join(x.gid)+'\n')
PY
SUSIE_MAXIT=2000 nice -n 5 $R s03_susie.R $O/work $O/work/pilot/res_maxit2000 3 $O/work/pilot/groups_nc.txt $O/work/pilot/pairs_nc.txt; echo FIN_MAXIT
SUSIE_L=10 nice -n 5 $R s03_susie.R $O/work $O/work/pilot/res_L10 8 $O/work/pilot/groups.txt $O/work/pilot/pairs.txt; echo FIN_L10
