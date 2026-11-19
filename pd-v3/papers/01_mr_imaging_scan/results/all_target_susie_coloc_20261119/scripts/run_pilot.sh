O=/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan/results/all_target_susie_coloc_20261119; cd $O/scripts
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
rm -f $O/work/groups/*.error
/public/home/mengxl/dzy/pd-v3_runtime/bin/python -u s02_prep.py 4 $O/work/pilot/groups.txt 2>&1 | tail -2 | cut -c1-200
R=/public/home/mengxl/dzy/envs/rcoloc/bin/Rscript
SUSIE_PILOT_N=1 nice -n 5 $R s03_susie.R $O/work $O/work/pilot/res_pilotN 8 $O/work/pilot/groups.txt $O/work/pilot/pairs.txt
echo FINISHED_PILOTN
SUSIE_PILOT_N=0 nice -n 5 $R s03_susie.R $O/work $O/work/pilot/res_trueN 8 $O/work/pilot/groups.txt $O/work/pilot/pairs.txt
echo FINISHED_TRUEN
