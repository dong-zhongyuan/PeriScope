O=/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan/results/all_target_susie_coloc_20261119; cd $O/scripts
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
SUSIE_L=5 SUSIE_MAXIT=2000 SUSIE_PILOT_N=0 nice -n 5 /public/home/mengxl/dzy/envs/rcoloc/bin/Rscript s03_susie.R $O/work $O/work/results 14
echo FIN_FULL
